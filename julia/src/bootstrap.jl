# Iterated bootstrap algorithm for pH calibration

# Default measurement uncertainties by nucleus (ppm)
const DEFAULT_σ_δ = Dict(
    "H"  => 0.003,   # 1H - high resolution, good shimming
    "C"  => 0.05,    # 13C - broader lines
    "N"  => 0.1,     # 15N - broader lines
    "F"  => 0.03,    # 19F - intermediate
    "P"  => 0.01     # 31P - intermediate
)

"""
Get the measurement uncertainty for a given nucleus.
Falls back to 0.01 ppm if nucleus not found in dictionary.
"""
function get_σ_δ(σ_δ_dict::Dict{String,Float64}, nucleus::String)
    return get(σ_δ_dict, nucleus, 0.01)
end

"""
    iterated_bootstrap_pH(data, initial_fits; reference_indicator, pKa_reference, ...)

Perform iterated bootstrap calibration of pH indicators.

The algorithm:
1. Start with a reference indicator (e.g., phosphate) with known pKa
2. Calculate pH values using only the reference
3. Find the next uncalibrated indicator with pKa closest to the reference
4. Re-fit that indicator using the current calibrated pH values (using only data within pKa_window of the pKa)
5. Update pH estimates using all calibrated indicators (inverse-variance weighted)
6. Repeat steps 3-5 until all indicators are calibrated

Arguments:
    data: DataFrame with chemical shift data (from select_condition)
    initial_fits: Dict from fit_all_indicators
    reference_indicator: name of reference indicator (default "phosphate")
    pKa_reference: literature pKa for reference at this (T, I) condition
    σ_δ: Dict mapping nucleus names to measurement uncertainties (ppm).
         Defaults to DEFAULT_σ_δ. Keys should be "H", "C", "N", "F", "P" etc.
    pKa_window: only use data within ±pKa_window pH units of the pKa for bootstrap
                corrections (default 1.5). Set to Inf to use all data.

Returns a NamedTuple with:
    - fits: Dict of BootstrapResult for each indicator
    - pH_data: DataFrame with nominal_pH, per-indicator pH estimates, and pH_combined
    - calibration_order: Vector of indicator names in calibration order
"""
function iterated_bootstrap_pH(data, initial_fits::Dict{String,IndicatorFit};
    reference_indicator="phosphate",
    pKa_reference,
    σ_δ::Dict{String,Float64}=DEFAULT_σ_δ,
    pKa_window::Float64=1.5)

    # Check reference indicator is present
    if !haskey(initial_fits, reference_indicator)
        error("Reference indicator '$reference_indicator' not found. Available: $(keys(initial_fits))")
    end

    ref_fit = initial_fits[reference_indicator]
    pKa_lit = pKa_reference ± 0.01  # assume small uncertainty on literature value
    σ_δ_ref = get_σ_δ(σ_δ, ref_fit.nucleus)

    println("\n" * "="^70)
    println("ITERATED BOOTSTRAP ANALYSIS")
    println("="^70)
    println("Reference: $reference_indicator with literature pKa = $pKa_lit")
    println("pKa window: ±$pKa_window pH units")
    println("Shift uncertainties (σ_δ): ", σ_δ)

    # -------------------------------------------------------------------------
    # Initialize pH_data DataFrame with unique pH values
    # -------------------------------------------------------------------------
    pH_nominal_unique = sort(unique(data.nominal_pH))
    n_samples = length(pH_nominal_unique)

    pH_data = DataFrame(
        nominal_pH = pH_nominal_unique
    )

    # Add columns for each indicator's pH estimate
    for ind in keys(initial_fits)
        pH_data[!, Symbol("pH_", ind)] = Vector{Union{Missing,Measurement{Float64}}}(missing, n_samples)
    end
    pH_data[!, :pH_combined] = Vector{Union{Missing,Measurement{Float64}}}(missing, n_samples)

    # -------------------------------------------------------------------------
    # Initialize results storage
    # -------------------------------------------------------------------------
    bootstrap_fits = Dict{String,BootstrapResult}()
    calibrated = Set{String}()
    calibration_order = String[]

    # -------------------------------------------------------------------------
    # Step 1: Calculate pH from reference indicator
    # -------------------------------------------------------------------------
    ref_data = filter(row -> row.indicator == reference_indicator, data)

    for i in 1:n_samples
        pH_nom = pH_data.nominal_pH[i]
        ref_row = filter(row -> row.nominal_pH == pH_nom, ref_data)

        if nrow(ref_row) == 0
            continue
        end

        δ_obs = ref_row.delta_obs[1]

        if ref_fit.n_pKa == 1
            pH_est = pH_from_shift(δ_obs, pKa_lit, get_δ_HA(ref_fit), get_δ_A(ref_fit);
                                   σ_δ=σ_δ_ref)
        else
            # For 2-pKa reference indicator, use numerical inversion
            # For now, only support 1-pKa references
            error("2-pKa reference indicators not yet supported")
        end

        pH_data[i, Symbol("pH_", reference_indicator)] = pH_est
    end

    # Add reference to calibrated set
    push!(calibrated, reference_indicator)
    push!(calibration_order, reference_indicator)

    bootstrap_fits[reference_indicator] = BootstrapResult(
        reference_indicator, ref_fit.nucleus,
        ref_fit.pKa[1],     # initial pKa from electrode pH
        pKa_lit,            # literature value used
        get_δ_HA(ref_fit),
        get_δ_A(ref_fit),
        count(!ismissing, pH_data[!, Symbol("pH_", reference_indicator)]),
        true;               # is reference
        calibration_tier=0
    )

    # Update combined pH using reference only
    update_combined_pH!(pH_data, bootstrap_fits, calibrated)

    n_reliable = count(row -> !ismissing(row.pH_combined) && is_reliable(row.pH_combined),
                      eachrow(pH_data))
    println("\nTier 0: $reference_indicator provides reliable pH for $n_reliable / $n_samples samples")

    # -------------------------------------------------------------------------
    # Step 2: Iteratively calibrate other indicators
    # -------------------------------------------------------------------------
    uncalibrated = setdiff(Set(keys(initial_fits)), calibrated)
    tier = 1

    while !isempty(uncalibrated)
        # Find uncalibrated indicator with pKa closest to any calibrated indicator
        best_ind = nothing
        best_distance = Inf

        for ind in uncalibrated
            init_fit = initial_fits[ind]
            # Use the first pKa for distance calculation
            ind_pKa = Measurements.value(init_fit.pKa[1])

            for cal_ind in calibrated
                cal_fit = bootstrap_fits[cal_ind]
                cal_pKa = Measurements.value(cal_fit.pKa_bootstrap[1])

                distance = abs(ind_pKa - cal_pKa)
                if distance < best_distance
                    best_distance = distance
                    best_ind = ind
                end
            end
        end

        if best_ind === nothing
            break
        end

        # Calibrate this indicator
        init_fit = initial_fits[best_ind]
        ind_data = filter(row -> row.indicator == best_ind, data)

        # Get the initial pKa estimate for filtering by proximity
        # For 2-pKa indicators, use both pKa values for the window
        pKa_centers = [Measurements.value(pKa) for pKa in init_fit.pKa]

        # Collect data points with reliable combined pH that are within pKa_window of any pKa
        pH_true_vals = Measurement{Float64}[]
        δ_vals = Float64[]
        n_outside_window = 0

        for row in eachrow(ind_data)
            pH_row_idx = findfirst(r -> r.nominal_pH == row.nominal_pH, eachrow(pH_data))

            if pH_row_idx === nothing
                continue
            end

            pH_combined = pH_data.pH_combined[pH_row_idx]

            if ismissing(pH_combined) || !is_reliable(pH_combined)
                continue
            end

            # Check if pH is within pKa_window of any pKa
            pH_val = Measurements.value(pH_combined)
            within_window = any(abs(pH_val - pKa_c) <= pKa_window for pKa_c in pKa_centers)

            if !within_window
                n_outside_window += 1
                continue
            end

            push!(pH_true_vals, pH_combined)
            push!(δ_vals, row.delta_obs)
        end

        n_usable = length(pH_true_vals)

        if n_usable < 3
            @warn "Insufficient reliable pH values for $best_ind ($n_usable points within pKa window, $n_outside_window excluded), skipping"
            delete!(uncalibrated, best_ind)
            continue
        end

        # Fit pKa using calibrated pH values
        pH_true_values = [Measurements.value(pH) for pH in pH_true_vals]

        if init_fit.n_pKa == 1
            result = refit_indicator_1pKa(best_ind, init_fit, pH_true_values, δ_vals,
                                          n_usable, tier)
        else
            result = refit_indicator_2pKa(best_ind, init_fit, pH_true_values, δ_vals,
                                          n_usable, tier)
        end

        bootstrap_fits[best_ind] = result
        push!(calibrated, best_ind)
        push!(calibration_order, best_ind)
        delete!(uncalibrated, best_ind)

        # Update pH estimates from this indicator
        σ_δ_ind = get_σ_δ(σ_δ, init_fit.nucleus)
        update_indicator_pH!(pH_data, best_ind, ind_data, result, σ_δ_ind)

        # Update combined pH
        update_combined_pH!(pH_data, bootstrap_fits, calibrated)

        n_reliable = count(row -> !ismissing(row.pH_combined) && is_reliable(row.pH_combined),
                          eachrow(pH_data))
        pKa_str = result.n_pKa == 1 ?
            @sprintf("%.3f", Measurements.value(result.pKa_bootstrap[1])) :
            @sprintf("%.3f, %.3f", Measurements.value(result.pKa_bootstrap[1]),
                                   Measurements.value(result.pKa_bootstrap[2]))

        window_note = n_outside_window > 0 ? ", $n_outside_window pts outside pKa window" : ""
        println("Tier $tier: $best_ind (pKa = $pKa_str, $n_usable pts$window_note) → $n_reliable reliable pH values")

        tier += 1
    end

    # Print summary
    print_bootstrap_summary(bootstrap_fits, calibration_order)
    print_pH_comparison(pH_data)

    return (fits=bootstrap_fits, pH_data=pH_data, calibration_order=calibration_order)
end

"""
Refit a single-pKa indicator using calibrated pH values.
"""
function refit_indicator_1pKa(indicator, init_fit, pH_values, δ_values, n_usable, tier)
    δ_HA_val = Measurements.value(get_δ_HA(init_fit))
    δ_A_val = Measurements.value(get_δ_A(init_fit))
    pKa_guess = Measurements.value(init_fit.pKa[1])

    try
        model(pH, p) = hh_fixed_limits(pH, p, δ_HA_val, δ_A_val)
        fit = curve_fit(model, pH_values, δ_values, [pKa_guess])

        pKa_new = coef(fit)[1] ± stderror(fit)[1]

        return BootstrapResult(
            indicator, init_fit.nucleus,
            init_fit.pKa[1],   # initial
            pKa_new,           # bootstrap
            get_δ_HA(init_fit),
            get_δ_A(init_fit),
            n_usable,
            false;
            calibration_tier=tier
        )
    catch e
        @warn "Bootstrap fitting failed for $indicator: $e"
        return BootstrapResult(
            indicator, init_fit.nucleus,
            init_fit.pKa[1], init_fit.pKa[1],
            get_δ_HA(init_fit), get_δ_A(init_fit),
            n_usable, false;
            calibration_tier=tier
        )
    end
end

"""
Refit a two-pKa indicator using calibrated pH values.
"""
function refit_indicator_2pKa(indicator, init_fit, pH_values, δ_values, n_usable, tier)
    δ_0_val = Measurements.value(init_fit.δ_limits[1])
    δ_1_val = Measurements.value(init_fit.δ_limits[2])
    δ_2_val = Measurements.value(init_fit.δ_limits[3])
    pKa1_guess = Measurements.value(init_fit.pKa[1])
    pKa2_guess = Measurements.value(init_fit.pKa[2])

    try
        model(pH, p) = hh_2pKa_fixed_limits(pH, p, δ_0_val, δ_1_val, δ_2_val)
        fit = curve_fit(model, pH_values, δ_values, [pKa1_guess, pKa2_guess])

        pKa1_new = coef(fit)[1] ± stderror(fit)[1]
        pKa2_new = coef(fit)[2] ± stderror(fit)[2]

        return BootstrapResult(
            indicator, init_fit.nucleus,
            2,  # n_pKa
            init_fit.pKa,                    # initial
            [pKa1_new, pKa2_new],            # bootstrap
            init_fit.δ_limits,
            n_usable,
            false,
            tier
        )
    catch e
        @warn "Bootstrap fitting failed for $indicator: $e"
        return BootstrapResult(
            indicator, init_fit.nucleus,
            2,
            init_fit.pKa, init_fit.pKa,
            init_fit.δ_limits,
            n_usable, false, tier
        )
    end
end

"""
Update pH estimates from a single indicator in pH_data.
"""
function update_indicator_pH!(pH_data, indicator, ind_data, result, σ_δ)
    col = Symbol("pH_", indicator)

    for i in 1:nrow(pH_data)
        pH_nom = pH_data.nominal_pH[i]
        ind_row = filter(row -> row.nominal_pH == pH_nom, ind_data)

        if nrow(ind_row) == 0
            continue
        end

        δ_obs = ind_row.delta_obs[1]

        if result.n_pKa == 1
            pH_est = pH_from_shift(δ_obs,
                                   result.pKa_bootstrap[1],
                                   result.δ_limits[1],
                                   result.δ_limits[2];
                                   σ_δ=σ_δ)
        else
            # For 2-pKa, use numerical inversion
            pH_est = pH_from_shift_2pKa(δ_obs,
                                        result.pKa_bootstrap[1],
                                        result.pKa_bootstrap[2],
                                        result.δ_limits[1],
                                        result.δ_limits[2],
                                        result.δ_limits[3];
                                        σ_δ=σ_δ)
        end

        pH_data[i, col] = pH_est
    end
end

"""
Update combined pH using inverse-variance weighting of all calibrated indicators.
"""
function update_combined_pH!(pH_data, bootstrap_fits, calibrated)
    for i in 1:nrow(pH_data)
        estimates = Measurement{Float64}[]

        for ind in calibrated
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
        else
            pH_data[i, :pH_combined] = missing
        end
    end
end

"""
Print summary of bootstrap fit results.
"""
function print_bootstrap_summary(results::Dict{String,BootstrapResult}, calibration_order)
    println("\n" * "─"^80)
    println("BOOTSTRAP RESULTS (in calibration order)")
    println("─"^80)
    @printf("%-15s %3s %4s %10s %10s %10s %5s\n",
        "Indicator", "Nuc", "Tier", "pKa(elec)", "pKa(boot)", "Δ(pKa)", "n")
    println("─"^80)

    for ind in calibration_order
        r = results[ind]

        if r.n_pKa == 1
            pKa_init = Measurements.value(r.pKa_initial[1])
            pKa_boot = Measurements.value(r.pKa_bootstrap[1])
            Δ = pKa_boot - pKa_init

            marker = r.is_reference ? " (REF)" : ""
            tier_str = r.is_reference ? "REF" : string(r.calibration_tier)

            @printf("%-15s %3s %4s %10.3f %10.3f %+10.3f %5d%s\n",
                r.indicator, r.nucleus, tier_str,
                pKa_init, pKa_boot, Δ, r.n_points_used, marker)
            @printf("%-15s         ±%-8.3f ±%-8.3f\n",
                "",
                Measurements.uncertainty(r.pKa_initial[1]),
                Measurements.uncertainty(r.pKa_bootstrap[1]))
        else
            # Two-pKa case
            for j in 1:2
                pKa_init = Measurements.value(r.pKa_initial[j])
                pKa_boot = Measurements.value(r.pKa_bootstrap[j])
                Δ = pKa_boot - pKa_init
                tier_str = r.is_reference ? "REF" : string(r.calibration_tier)
                suffix = j == 1 ? " (pKa1)" : " (pKa2)"

                @printf("%-15s %3s %4s %10.3f %10.3f %+10.3f %5d%s\n",
                    r.indicator * suffix, r.nucleus, tier_str,
                    pKa_init, pKa_boot, Δ, r.n_points_used, "")
            end
        end
    end
    println("─"^80)
end

"""
Print comparison of nominal vs fitted pH values.
"""
function print_pH_comparison(pH_data)
    println("\n" * "─"^70)
    println("pH COMPARISON: Nominal (electrode) vs Bootstrap (combined)")
    println("─"^70)
    @printf("%8s %12s %12s %10s %6s\n", "Nominal", "Bootstrap", "σ", "Δ(pH)", "n_ind")
    println("─"^70)

    for row in eachrow(pH_data)
        if ismissing(row.pH_combined)
            @printf("%8.2f %12s %12s %10s %6s\n", row.nominal_pH, "—", "—", "—", "0")
        else
            pH_boot = Measurements.value(row.pH_combined)
            σ = Measurements.uncertainty(row.pH_combined)
            Δ = pH_boot - row.nominal_pH

            # Count number of indicators contributing
            n_ind = 0
            for col in names(row)
                if startswith(string(col), "pH_") && col != :pH_combined
                    val = row[col]
                    if !ismissing(val) && is_reliable(val)
                        n_ind += 1
                    end
                end
            end

            @printf("%8.2f %12.3f %12.3f %+10.3f %6d\n", row.nominal_pH, pH_boot, σ, Δ, n_ind)
        end
    end
    println("─"^70)
end

# Keep the old bootstrap_pH function for backwards compatibility
"""
    bootstrap_pH(data, initial_fits; pKa_phosphate_lit, ...)

Single-pass bootstrap (legacy function). Use `iterated_bootstrap_pH` for improved calibration.
"""
function bootstrap_pH(data, initial_fits::Dict{String,IndicatorFit};
    pKa_phosphate_lit,
    σ_δ::Dict{String,Float64}=DEFAULT_σ_δ,
    pKa_window::Float64=1.5)

    return iterated_bootstrap_pH(data, initial_fits;
        reference_indicator="phosphate",
        pKa_reference=pKa_phosphate_lit,
        σ_δ=σ_δ,
        pKa_window=pKa_window)
end
