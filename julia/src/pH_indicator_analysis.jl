#=
pH Indicator Analysis Module
============================

Analysis of NMR chemical shift data for pH indicator calibration.

This module provides functions to:
1. Load and organise chemical shift data from CSV
2. Fit Henderson-Hasselbalch models to determine pKa and limiting shifts
3. Bootstrap a refined pH scale using phosphate as a reference
4. Produce diagnostic plots and summary tables

Usage:
    include("pH_indicator_analysis.jl")

    # Load data
    data = load_indicator_data("chemical_shifts.csv")

    # Select one (T, I) condition
    subset = select_condition(data, temperature=298.0, ionic_strength="low")

    # Initial fits using electrode pH
    initial_results = fit_all_indicators(subset)
    plot_initial_fits(subset, initial_results)

    # Bootstrap using literature phosphate pKa
    bootstrap_results = bootstrap_pH(subset, initial_results, pKa_phosphate_lit=7.20)
    plot_bootstrap_results(subset, bootstrap_results)
=#

using CSV
using DataFrames
using LsqFit
using Measurements
using Plots
using Printf
using Statistics

# ============================================================================
# HENDERSON-HASSELBALCH MODEL
# ============================================================================

"""
Henderson-Hasselbalch model for chemical shift vs pH.

    δ_obs = δ_HA + (δ_A - δ_HA) / (1 + 10^(pKa - pH))

Parameters:
    p[1] = pKa
    p[2] = δ_HA (limiting shift of acid/protonated form)
    p[3] = δ_A (limiting shift of base/deprotonated form)
"""
function henderson_hasselbalch(pH, p)
    pKa, δ_HA, δ_A = p
    return @. δ_HA + (δ_A - δ_HA) / (1 + 10^(pKa - pH))
end

"""
Henderson-Hasselbalch model with fixed limiting shifts (for bootstrap fitting).
Only pKa is fitted.
"""
function hh_fixed_limits(pH, p, δ_HA, δ_A)
    pKa = p[1]
    return @. δ_HA + (δ_A - δ_HA) / (1 + 10^(pKa - pH))
end

"""
Calculate pH from observed chemical shift by inverting Henderson-Hasselbalch.

Returns pH as a Measurement (value ± uncertainty) or missing if outside valid range.

Arguments:
    δ_obs: observed chemical shift (can be Measurement or Number)
    pKa: pKa value (can be Measurement or Number)
    δ_HA: limiting shift of acid form
    δ_A: limiting shift of base form
    σ_δ: measurement uncertainty in δ_obs (used if δ_obs is not a Measurement)
"""
function pH_from_shift(δ_obs, pKa, δ_HA, δ_A; σ_δ=0.01)
    # Extract values for range checking
    δ_obs_val = Measurements.value(δ_obs)
    δ_HA_val = Measurements.value(δ_HA)
    δ_A_val = Measurements.value(δ_A)

    # Check if δ_obs is between the limits
    δ_min, δ_max = minmax(δ_HA_val, δ_A_val)
    if δ_obs_val <= δ_min || δ_obs_val >= δ_max
        return missing
    end

    # Convert to Measurements if needed
    δ_obs_m = δ_obs isa Measurement ? δ_obs : δ_obs ± σ_δ
    pKa_m = pKa isa Measurement ? pKa : pKa ± 0.0
    δ_HA_m = δ_HA isa Measurement ? δ_HA : δ_HA ± 0.0
    δ_A_m = δ_A isa Measurement ? δ_A : δ_A ± 0.0

    # Calculate pH with uncertainty propagation
    ratio = (δ_obs_m - δ_HA_m) / (δ_A_m - δ_obs_m)
    return pKa_m + log10(ratio)
end

"""
Check if a pH estimate is reliable (uncertainty below threshold).
"""
function is_reliable(pH_estimate; threshold=0.15)
    ismissing(pH_estimate) && return false
    return Measurements.uncertainty(pH_estimate) < threshold
end

# ============================================================================
# DATA LOADING AND SELECTION
# ============================================================================

"""
    load_indicator_data(filename) -> DataFrame

Load chemical shift data from CSV file.

Expected columns:
    - expt_number: experiment identifier (ignored in analysis)
    - ionic_strength: ionic strength condition (e.g., "low", "150", "300")
    - temperature: temperature in K
    - nominal_pH: pH measured by electrode at 298 K
    - P_xxx, F_xxx, H_xxx: chemical shifts for each indicator
      (prefix indicates nucleus: P=³¹P, F=¹⁹F, H=¹H)

Returns a DataFrame with the data reorganised into long format.
"""
function load_indicator_data(filename)
    raw = CSV.read(filename, DataFrame)

    # Identify indicator columns (those starting with P_, F_, H_)
    indicator_cols = filter(names(raw)) do col
        startswith(string(col), "P_") ||
            startswith(string(col), "F_") ||
            startswith(string(col), "H_")
    end

    # Build long-format dataframe
    rows = []
    for row in eachrow(raw)
        for col in indicator_cols
            col_str = string(col)
            # Extract nucleus and indicator name
            nucleus = col_str[1:1]  # P, F, or H
            indicator = col_str[3:end]  # everything after "X_"

            # Get the chemical shift value
            δ = row[col]

            # Skip missing values
            ismissing(δ) && continue

            push!(rows, (
                expt_number=row.expt_number,
                ionic_strength=Float64(row.ionic_strength),
                temperature=Float64(row.temperature),
                nominal_pH=Float64(row.nominal_pH),
                nucleus=nucleus,
                indicator=indicator,
                delta_obs=Float64(δ)
            ))
        end
    end

    return DataFrame(rows)
end

"""
    select_condition(data, temperature, ionic_strength; tol_T=1.0, tol_I=0.01) -> DataFrame

Select data for a single (T, I) condition.

Arguments:
    data: DataFrame from load_indicator_data
    temperature: target temperature in K (will select closest)
    ionic_strength: target ionic strength in M (will select closest)
    tol_T: tolerance for temperature matching (K), warns if closest is further
    tol_I: tolerance for ionic strength matching (M), warns if closest is further
"""
function select_condition(data, temperature, ionic_strength; tol_T=1.0, tol_I=0.01)
    # Find closest ionic strength
    I_values = unique(data.ionic_strength)
    closest_I = I_values[argmin(abs.(I_values .- ionic_strength))]

    if abs(closest_I - ionic_strength) > tol_I
        @warn "Requested I=$ionic_strength M, using closest available I=$closest_I M"
    end

    subset = filter(row -> row.ionic_strength == closest_I, data)

    # Find closest temperature
    temps = unique(subset.temperature)
    closest_T = temps[argmin(abs.(temps .- temperature))]

    if abs(closest_T - temperature) > tol_T
        @warn "Requested T=$temperature K, using closest available T=$closest_T K"
    end

    subset = filter(row -> row.temperature == closest_T, subset)

    println("Selected $(nrow(subset)) measurements at T=$(closest_T) K, I=$(closest_I) M")
    println("Indicators: ", join(sort(unique(subset.indicator)), ", "))

    return subset
end

"""
    list_conditions(data)

Print available (temperature, ionic_strength) combinations in the data.
"""
function list_conditions(data)
    conditions = unique(select(data, [:temperature, :ionic_strength]))
    sort!(conditions, [:temperature, :ionic_strength])

    println("\nAvailable conditions:")
    println("─"^30)
    for row in eachrow(conditions)
        n = count(r -> r.temperature == row.temperature &&
                r.ionic_strength == row.ionic_strength, eachrow(data))
        println("  T = $(row.temperature) K, I = $(row.ionic_strength)  ($n measurements)")
    end
    println()
end

# ============================================================================
# INITIAL FITTING
# ============================================================================

"""
Result structure for a single indicator fit.
"""
struct IndicatorFit
    indicator::String
    nucleus::String
    pKa::Measurement{Float64}
    δ_HA::Measurement{Float64}
    δ_A::Measurement{Float64}
    n_points::Int
    rmsd::Float64
    pH_range::Tuple{Float64,Float64}
end

"""
    fit_indicator(pH_values, delta_values; indicator="", nucleus="") -> IndicatorFit

Fit Henderson-Hasselbalch model to chemical shift vs pH data.
"""
function fit_indicator(pH_values, delta_values; indicator="unknown", nucleus="")
    n = length(pH_values)

    if n < 4
        error("Need at least 4 data points to fit, got $n")
    end

    # Initial guesses
    δ_min, δ_max = extrema(delta_values)
    mid_δ = (δ_min + δ_max) / 2

    # Find pH closest to midpoint of transition
    mid_idx = argmin(abs.(delta_values .- mid_δ))
    pKa_guess = pH_values[mid_idx]

    # Determine which limit corresponds to acid vs base
    # At low pH, we have the acid form
    low_pH_idx = argmin(pH_values)
    high_pH_idx = argmax(pH_values)

    δ_HA_guess = delta_values[low_pH_idx]
    δ_A_guess = delta_values[high_pH_idx]

    p0 = [pKa_guess, δ_HA_guess, δ_A_guess]

    # Fit
    try
        fit = curve_fit(henderson_hasselbalch, pH_values, delta_values, p0)

        params = coef(fit)
        errors = stderror(fit)

        # Calculate RMSD
        predicted = henderson_hasselbalch(pH_values, params)
        rmsd = sqrt(mean((delta_values .- predicted) .^ 2))

        return IndicatorFit(
            indicator,
            nucleus,
            params[1] ± errors[1],
            params[2] ± errors[2],
            params[3] ± errors[3],
            n,
            rmsd,
            extrema(pH_values)
        )
    catch e
        error("Fitting failed for $indicator: $e")
    end
end

"""
    fit_all_indicators(data) -> Dict{String, IndicatorFit}

Fit Henderson-Hasselbalch model to all indicators in the dataset.

Returns a dictionary mapping indicator names to IndicatorFit results.
"""
function fit_all_indicators(data)
    results = Dict{String,IndicatorFit}()

    indicators = unique(data.indicator)

    println("\n" * "="^70)
    println("INITIAL HENDERSON-HASSELBALCH FITS (using electrode pH)")
    println("="^70)

    for ind in sort(indicators)
        subset = filter(row -> row.indicator == ind, data)
        nucleus = subset.nucleus[1]

        pH_vals = Float64.(subset.nominal_pH)
        delta_vals = Float64.(subset.delta_obs)

        try
            fit = fit_indicator(pH_vals, delta_vals; indicator=ind, nucleus=nucleus)
            results[ind] = fit
        catch e
            @warn "Could not fit $ind: $e"
        end
    end

    # Print summary table
    print_fit_summary(results)

    return results
end

"""
Print a formatted summary table of fit results.
"""
function print_fit_summary(results::Dict{String,IndicatorFit})
    println("\n" * "─"^70)
    @printf("%-15s %3s %8s %12s %12s %5s %8s\n",
        "Indicator", "Nuc", "pKa", "δ_HA (ppm)", "δ_A (ppm)", "n", "RMSD")
    println("─"^70)

    for ind in sort(collect(keys(results)))
        r = results[ind]
        @printf("%-15s %3s %8.3f %12.4f %12.4f %5d %8.4f\n",
            r.indicator, r.nucleus,
            Measurements.value(r.pKa),
            Measurements.value(r.δ_HA),
            Measurements.value(r.δ_A),
            r.n_points, r.rmsd)
        @printf("%-15s     ±%-6.3f ±%-10.4f ±%-10.4f\n",
            "",
            Measurements.uncertainty(r.pKa),
            Measurements.uncertainty(r.δ_HA),
            Measurements.uncertainty(r.δ_A))
    end
    println("─"^70)
end

# ============================================================================
# PLOTTING - INITIAL FITS
# ============================================================================

"""
    plot_initial_fits(data, results; save_path=nothing) -> Plot

Create a multi-panel plot showing all indicator fits.
"""
function plot_initial_fits(data, results::Dict{String,IndicatorFit}; save_path=nothing)
    indicators = sort(collect(keys(results)))
    n_ind = length(indicators)

    # Determine grid layout
    n_cols = min(3, n_ind)
    n_rows = ceil(Int, n_ind / n_cols)

    plots = []

    for ind in indicators
        r = results[ind]
        subset = filter(row -> row.indicator == ind, data)

        pH_data = Float64.(subset.nominal_pH)
        δ_data = Float64.(subset.delta_obs)

        # Generate smooth curve
        pH_curve = range(minimum(pH_data) - 0.3, maximum(pH_data) + 0.3, length=100)
        params = [Measurements.value(r.pKa),
            Measurements.value(r.δ_HA),
            Measurements.value(r.δ_A)]
        δ_curve = henderson_hasselbalch(collect(pH_curve), params)

        p = plot(pH_curve, δ_curve,
            label="Fit (pKa=$(round(Measurements.value(r.pKa), digits=2)))",
            linewidth=2, color=:blue)
        scatter!(p, pH_data, δ_data,
            label="Data", markersize=5, color=:red)

        xlabel!(p, "pH (electrode)")
        ylabel!(p, "δ (ppm)")
        title!(p, "$(r.indicator) ($(r.nucleus))")

        push!(plots, p)
    end

    # Combine into single figure
    fig = plot(plots..., layout=(n_rows, n_cols), size=(400 * n_cols, 350 * n_rows),
        margin=5Plots.mm, legend=:best)

    if !isnothing(save_path)
        savefig(fig, save_path)
        println("Saved plot to $save_path")
    end

    return fig
end

# ============================================================================
# BOOTSTRAP PROCEDURE
# ============================================================================

"""
Result structure for bootstrap analysis.
"""
struct BootstrapResult
    indicator::String
    nucleus::String
    pKa_initial::Measurement{Float64}      # from electrode pH
    pKa_bootstrap::Measurement{Float64}    # from bootstrap
    δ_HA::Measurement{Float64}
    δ_A::Measurement{Float64}
    n_points_used::Int                      # points with reliable pH_true
    is_reference::Bool                      # true for phosphate (fixed pKa)
end

"""
    bootstrap_pH(data, initial_fits, pKa_phosphate_lit; σ_δ_phosphate=0.01) -> Dict

Bootstrap pH scale using phosphate as reference.

Arguments:
    data: DataFrame with chemical shift data
    initial_fits: Dict from fit_all_indicators
    pKa_phosphate_lit: literature pKa for phosphate at this (T, I) condition
    σ_δ_phosphate: measurement uncertainty for ³¹P chemical shifts (ppm)

Returns a NamedTuple with:
    - fits: Dict of BootstrapResult for each indicator
    - pH_data: DataFrame with nominal_pH, pH_true, and per-indicator pH estimates
"""
function bootstrap_pH(data, initial_fits::Dict{String,IndicatorFit};
    pKa_phosphate_lit, σ_δ_phosphate=0.01, σ_δ_other=0.005)

    # Check phosphate is present
    if !haskey(initial_fits, "phosphate")
        error("Phosphate not found in initial fits. Available: $(keys(initial_fits))")
    end

    phos_fit = initial_fits["phosphate"]
    pKa_lit = pKa_phosphate_lit ± 0.01  # assume small uncertainty on literature value

    println("\n" * "="^70)
    println("BOOTSTRAP ANALYSIS")
    println("="^70)
    println("Reference: phosphate with literature pKa = $pKa_lit")

    # -------------------------------------------------------------------------
    # Step 1: Calculate pH_true from phosphate for all samples
    # -------------------------------------------------------------------------

    # Get unique pH values from the data
    pH_nominal_unique = sort(unique(data.nominal_pH))

    # Get phosphate chemical shifts at each pH
    phos_data = filter(row -> row.indicator == "phosphate", data)

    # Create pH_data DataFrame
    pH_data = DataFrame(
        nominal_pH=Float64[],
        pH_phosphate=Union{Missing,Measurement{Float64}}[],
        σ_pH_phosphate=Union{Missing,Float64}[]
    )

    for pH_nom in pH_nominal_unique
        phos_row = filter(row -> row.nominal_pH == pH_nom, phos_data)

        if nrow(phos_row) == 0
            push!(pH_data, (nominal_pH=pH_nom, pH_phosphate=missing, σ_pH_phosphate=missing))
            continue
        end

        δ_obs = phos_row.delta_obs[1]
        pH_true = pH_from_shift(δ_obs, pKa_lit, phos_fit.δ_HA, phos_fit.δ_A; σ_δ=σ_δ_phosphate)

        if ismissing(pH_true)
            push!(pH_data, (nominal_pH=pH_nom, pH_phosphate=missing, σ_pH_phosphate=missing))
        else
            push!(pH_data, (
                nominal_pH=pH_nom,
                pH_phosphate=pH_true,
                σ_pH_phosphate=Measurements.uncertainty(pH_true)
            ))
        end
    end

    # Count reliable phosphate pH estimates
    n_reliable_phos = count(row -> !ismissing(row.pH_phosphate) &&
            is_reliable(row.pH_phosphate),
        eachrow(pH_data))
    println("\nPhosphate provides reliable pH for $n_reliable_phos / $(nrow(pH_data)) samples")

    # -------------------------------------------------------------------------
    # Step 2: Refit other indicators using phosphate-derived pH
    # -------------------------------------------------------------------------

    bootstrap_fits = Dict{String,BootstrapResult}()

    # Add phosphate as reference (not refitted)
    bootstrap_fits["phosphate"] = BootstrapResult(
        "phosphate", "P",
        phos_fit.pKa,      # initial pKa from electrode pH
        pKa_lit,           # literature value used
        phos_fit.δ_HA,
        phos_fit.δ_A,
        n_reliable_phos,
        true               # is reference
    )

    # Process other indicators
    other_indicators = filter(ind -> ind != "phosphate", collect(keys(initial_fits)))

    for ind in sort(other_indicators)
        init_fit = initial_fits[ind]

        # Get data for this indicator
        ind_data = filter(row -> row.indicator == ind, data)

        # Match with pH_true values
        pH_true_vals = Measurement{Float64}[]
        δ_vals = Float64[]

        for row in eachrow(ind_data)
            pH_row = filter(r -> r.nominal_pH == row.nominal_pH, pH_data)

            if nrow(pH_row) == 0 || ismissing(pH_row.pH_phosphate[1])
                continue
            end

            pH_true = pH_row.pH_phosphate[1]

            if !is_reliable(pH_true)
                continue
            end

            push!(pH_true_vals, pH_true)
            push!(δ_vals, row.delta_obs)
        end

        n_usable = length(pH_true_vals)

        if n_usable < 3
            @warn "Insufficient reliable pH values for $ind ($n_usable points), skipping bootstrap"
            bootstrap_fits[ind] = BootstrapResult(
                ind, init_fit.nucleus,
                init_fit.pKa, init_fit.pKa,  # no change
                init_fit.δ_HA, init_fit.δ_A,
                n_usable, false
            )
            continue
        end

        # Fit pKa using pH_true, with fixed limiting shifts
        pH_true_values = [Measurements.value(pH) for pH in pH_true_vals]
        δ_HA_val = Measurements.value(init_fit.δ_HA)
        δ_A_val = Measurements.value(init_fit.δ_A)

        # Initial guess
        pKa_guess = Measurements.value(init_fit.pKa)

        try
            model(pH, p) = hh_fixed_limits(pH, p, δ_HA_val, δ_A_val)
            fit = curve_fit(model, pH_true_values, δ_vals, [pKa_guess])

            pKa_new = coef(fit)[1] ± stderror(fit)[1]

            bootstrap_fits[ind] = BootstrapResult(
                ind, init_fit.nucleus,
                init_fit.pKa,   # initial
                pKa_new,        # bootstrap
                init_fit.δ_HA,
                init_fit.δ_A,
                n_usable,
                false
            )
        catch e
            @warn "Bootstrap fitting failed for $ind: $e"
            bootstrap_fits[ind] = BootstrapResult(
                ind, init_fit.nucleus,
                init_fit.pKa, init_fit.pKa,
                init_fit.δ_HA, init_fit.δ_A,
                n_usable, false
            )
        end
    end

    # -------------------------------------------------------------------------
    # Step 3: Calculate pH from each indicator and combine
    # -------------------------------------------------------------------------

    # Add columns for each indicator's pH estimate
    for ind in keys(bootstrap_fits)
        pH_data[!, Symbol("pH_", ind)] = Union{Missing,Measurement{Float64}}[missing for _ in 1:nrow(pH_data)]
    end

    # Calculate pH from each indicator
    for ind in keys(bootstrap_fits)
        br = bootstrap_fits[ind]
        ind_data = filter(row -> row.indicator == ind, data)

        pKa_use = br.pKa_bootstrap
        σ_δ = br.nucleus == "P" ? σ_δ_phosphate : σ_δ_other

        for i in 1:nrow(pH_data)
            pH_nom = pH_data.nominal_pH[i]
            ind_row = filter(row -> row.nominal_pH == pH_nom, ind_data)

            if nrow(ind_row) == 0
                continue
            end

            δ_obs = ind_row.delta_obs[1]
            pH_est = pH_from_shift(δ_obs, pKa_use, br.δ_HA, br.δ_A; σ_δ=σ_δ)

            pH_data[i, Symbol("pH_", ind)] = pH_est
        end
    end

    # Combine pH estimates using inverse-variance weighting
    pH_data[!, :pH_combined] = Union{Missing,Measurement{Float64}}[missing for _ in 1:nrow(pH_data)]

    for i in 1:nrow(pH_data)
        estimates = Measurement{Float64}[]

        for ind in keys(bootstrap_fits)
            pH_est = pH_data[i, Symbol("pH_", ind)]
            if !ismissing(pH_est) && is_reliable(pH_est)
                push!(estimates, pH_est)
            end
        end

        if length(estimates) > 0
            # Inverse-variance weighted combination
            values = [Measurements.value(e) for e in estimates]
            uncertainties = [Measurements.uncertainty(e) for e in estimates]
            weights = 1 ./ uncertainties .^ 2

            pH_combined = sum(values .* weights) / sum(weights)
            σ_combined = 1 / sqrt(sum(weights))

            pH_data[i, :pH_combined] = pH_combined ± σ_combined
        end
    end

    # Print results
    print_bootstrap_summary(bootstrap_fits)
    print_pH_comparison(pH_data)

    return (fits=bootstrap_fits, pH_data=pH_data)
end

"""
Print summary of bootstrap fit results.
"""
function print_bootstrap_summary(results::Dict{String,BootstrapResult})
    println("\n" * "─"^70)
    @printf("%-15s %3s %10s %10s %10s %5s\n",
        "Indicator", "Nuc", "pKa(elec)", "pKa(boot)", "Δ(pKa)", "n")
    println("─"^70)

    for ind in sort(collect(keys(results)))
        r = results[ind]
        pKa_init = Measurements.value(r.pKa_initial)
        pKa_boot = Measurements.value(r.pKa_bootstrap)
        Δ = pKa_boot - pKa_init

        marker = r.is_reference ? " (REF)" : ""

        @printf("%-15s %3s %10.3f %10.3f %+10.3f %5d%s\n",
            r.indicator, r.nucleus,
            pKa_init, pKa_boot, Δ, r.n_points_used, marker)
        @printf("%-15s     ±%-8.3f ±%-8.3f\n",
            "",
            Measurements.uncertainty(r.pKa_initial),
            Measurements.uncertainty(r.pKa_bootstrap))
    end
    println("─"^70)
end

"""
Print comparison of nominal vs fitted pH values.
"""
function print_pH_comparison(pH_data)
    println("\n" * "─"^70)
    println("pH COMPARISON: Nominal (electrode) vs Bootstrap (combined)")
    println("─"^70)
    @printf("%8s %12s %12s %10s\n", "Nominal", "Bootstrap", "σ", "Δ(pH)")
    println("─"^70)

    for row in eachrow(pH_data)
        if ismissing(row.pH_combined)
            @printf("%8.2f %12s %12s %10s\n", row.nominal_pH, "—", "—", "—")
        else
            pH_boot = Measurements.value(row.pH_combined)
            σ = Measurements.uncertainty(row.pH_combined)
            Δ = pH_boot - row.nominal_pH
            @printf("%8.2f %12.3f %12.3f %+10.3f\n", row.nominal_pH, pH_boot, σ, Δ)
        end
    end
    println("─"^70)
end

# ============================================================================
# PLOTTING - BOOTSTRAP RESULTS
# ============================================================================

"""
    plot_bootstrap_results(data, bootstrap_result; save_path=nothing) -> Plot

Create diagnostic plots for bootstrap analysis.

Returns a combined figure with:
1. pH comparison: nominal vs bootstrap
2. Multi-panel fits using bootstrap pH
3. Residuals
"""
function plot_bootstrap_results(data, bootstrap_result; save_path=nothing)
    fits = bootstrap_result.fits
    pH_data = bootstrap_result.pH_data

    # -------------------------------------------------------------------------
    # Plot 1: pH comparison (nominal vs bootstrap)
    # -------------------------------------------------------------------------

    # Extract data for plotting
    pH_nom = Float64[]
    pH_boot = Float64[]
    pH_err = Float64[]

    for row in eachrow(pH_data)
        if !ismissing(row.pH_combined)
            push!(pH_nom, row.nominal_pH)
            push!(pH_boot, Measurements.value(row.pH_combined))
            push!(pH_err, Measurements.uncertainty(row.pH_combined))
        end
    end

    p_compare = scatter(pH_nom, pH_boot, yerror=pH_err,
        xlabel="pH (electrode)", ylabel="pH (bootstrap)",
        label="Data", markersize=6,
        title="pH Comparison")

    # Add 1:1 line
    pH_range = [minimum(pH_nom) - 0.5, maximum(pH_nom) + 0.5]
    plot!(p_compare, pH_range, pH_range, label="1:1", linestyle=:dash, color=:gray)

    # -------------------------------------------------------------------------
    # Plot 2: Residuals (bootstrap pH - nominal pH)
    # -------------------------------------------------------------------------

    residuals = pH_boot .- pH_nom

    p_resid = scatter(pH_nom, residuals,
        xlabel="pH (nominal)", ylabel="pH(boot) - pH(nom)",
        label=nothing, markersize=6,
        title="pH Residuals")
    hline!(p_resid, [0], linestyle=:dash, color=:gray, label=nothing)

    # Add ±0.1 reference lines
    hline!(p_resid, [-0.1, 0.1], linestyle=:dot, color=:lightgray, label=nothing)

    # -------------------------------------------------------------------------
    # Plot 3: Individual indicator fits using bootstrap pH
    # -------------------------------------------------------------------------

    indicators = sort(collect(keys(fits)))
    n_ind = length(indicators)
    n_cols = min(3, n_ind)
    n_rows = ceil(Int, n_ind / n_cols)

    ind_plots = []

    for ind in indicators
        br = fits[ind]
        ind_data_raw = filter(row -> row.indicator == ind, data)

        # Get bootstrap pH values for this indicator's data points
        pH_boot_vals = Float64[]
        δ_vals = Float64[]

        for row in eachrow(ind_data_raw)
            pH_row = filter(r -> r.nominal_pH == row.nominal_pH, pH_data)
            if nrow(pH_row) > 0 && !ismissing(pH_row.pH_combined[1])
                push!(pH_boot_vals, Measurements.value(pH_row.pH_combined[1]))
                push!(δ_vals, row.delta_obs)
            end
        end

        if length(pH_boot_vals) < 2
            continue
        end

        # Generate smooth curve
        pH_curve = range(minimum(pH_boot_vals) - 0.3, maximum(pH_boot_vals) + 0.3, length=100)
        params = [Measurements.value(br.pKa_bootstrap),
            Measurements.value(br.δ_HA),
            Measurements.value(br.δ_A)]
        δ_curve = henderson_hasselbalch(collect(pH_curve), params)

        p = plot(pH_curve, δ_curve,
            label="pKa=$(round(Measurements.value(br.pKa_bootstrap), digits=2))",
            linewidth=2, color=:blue)
        scatter!(p, pH_boot_vals, δ_vals,
            label="Data", markersize=5, color=:red)

        xlabel!(p, "pH (bootstrap)")
        ylabel!(p, "δ (ppm)")

        ref_marker = br.is_reference ? " [REF]" : ""
        title!(p, "$(br.indicator)$ref_marker")

        push!(ind_plots, p)
    end

    # -------------------------------------------------------------------------
    # Combine all plots
    # -------------------------------------------------------------------------

    # Top row: comparison and residuals
    p_top = plot(p_compare, p_resid, layout=(1, 2), size=(800, 350))

    # Bottom: indicator fits
    p_bottom = plot(ind_plots..., layout=(n_rows, n_cols),
        size=(400 * n_cols, 300 * n_rows))

    # Final combined figure
    fig = plot(p_top, p_bottom, layout=@layout([a{0.35h}; b]),
        size=(max(800, 400 * n_cols), 350 + 300 * n_rows))

    if !isnothing(save_path)
        savefig(fig, save_path)
        println("\nSaved plot to $save_path")
    end

    return fig
end

# ============================================================================
# CONVENIENCE FUNCTION: RUN COMPLETE ANALYSIS
# ============================================================================

"""
    run_analysis(filename; temperature, ionic_strength, pKa_phosphate_lit, save_plots=true)

Run complete analysis pipeline on a CSV file.

Arguments:
    filename: path to CSV file with chemical shift data
    temperature: target temperature in K
    ionic_strength: ionic strength condition (e.g., "low", "150")
    pKa_phosphate_lit: literature pKa for phosphate at this (T, I) condition
    save_plots: if true, save plots to PNG files

Returns a NamedTuple with all results.
"""
function run_analysis(filename; temperature, ionic_strength, pKa_phosphate_lit, save_plots=true)
    println("\n" * "="^70)
    println("pH INDICATOR ANALYSIS")
    println("="^70)
    println("File: $filename")
    println("Condition: T = $temperature K, I = $ionic_strength")
    println("Reference pKa (phosphate): $pKa_phosphate_lit")
    println("="^70)

    # Load data
    println("\nLoading data...")
    data = load_indicator_data(filename)
    list_conditions(data)

    # Select condition
    subset = select_condition(data, temperature, ionic_strength)

    # Initial fits
    initial_fits = fit_all_indicators(subset)

    if save_plots
        plot_initial_fits(subset, initial_fits; save_path="initial_fits.png")
    end

    # Bootstrap
    bootstrap_result = bootstrap_pH(subset, initial_fits; pKa_phosphate_lit=pKa_phosphate_lit)

    if save_plots
        plot_bootstrap_results(subset, bootstrap_result; save_path="bootstrap_results.png")
    end

    println("\n" * "="^70)
    println("ANALYSIS COMPLETE")
    println("="^70)

    return (
        data=subset,
        initial_fits=initial_fits,
        bootstrap=bootstrap_result
    )
end

# # Export public functions
# export load_indicator_data, select_condition, list_conditions
# export fit_all_indicators, plot_initial_fits
# export bootstrap_pH, plot_bootstrap_results
# export run_analysis
