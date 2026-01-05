# Fitting functions for pH indicator analysis

"""
    fit_indicator(pH_values, delta_values; indicator="", nucleus="", n_pKa=1, pKa_guess=nothing) -> IndicatorFit

Fit Henderson-Hasselbalch model to chemical shift vs pH data.

Arguments:
    pH_values: array of pH values
    delta_values: array of chemical shift values
    indicator: indicator name for labeling
    nucleus: nucleus type (P, F, H)
    n_pKa: number of pKa values (1 or 2)
    pKa_guess: initial guess for pKa (or vector of 2 for 2-pKa model)
"""
function fit_indicator(pH_values, delta_values;
                      indicator="unknown", nucleus="",
                      n_pKa=1, pKa_guess=nothing)
    n = length(pH_values)

    if n_pKa == 1
        return fit_indicator_1pKa(pH_values, delta_values;
                                  indicator=indicator, nucleus=nucleus,
                                  pKa_guess=pKa_guess)
    elseif n_pKa == 2
        return fit_indicator_2pKa(pH_values, delta_values;
                                  indicator=indicator, nucleus=nucleus,
                                  pKa_guess=pKa_guess)
    else
        error("n_pKa must be 1 or 2, got $n_pKa")
    end
end

"""
Fit single-pKa Henderson-Hasselbalch model.
"""
function fit_indicator_1pKa(pH_values, delta_values;
                           indicator="unknown", nucleus="",
                           pKa_guess=nothing)
    n = length(pH_values)

    if n < 4
        error("Need at least 4 data points to fit, got $n")
    end

    # Initial guesses
    δ_min, δ_max = extrema(delta_values)
    mid_δ = (δ_min + δ_max) / 2

    # Find pH closest to midpoint of transition
    if pKa_guess === nothing
        mid_idx = argmin(abs.(delta_values .- mid_δ))
        pKa_guess = pH_values[mid_idx]
    end

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
Fit two-pKa Henderson-Hasselbalch model.
"""
function fit_indicator_2pKa(pH_values, delta_values;
                           indicator="unknown", nucleus="",
                           pKa_guess=nothing)
    n = length(pH_values)

    if n < 6
        error("Need at least 6 data points for 2-pKa fit, got $n")
    end

    # Sort by pH for better initial guesses
    order = sortperm(pH_values)
    pH_sorted = pH_values[order]
    δ_sorted = delta_values[order]

    # Initial guesses for limiting shifts
    δ_0_guess = δ_sorted[1]      # low pH (fully protonated)
    δ_2_guess = δ_sorted[end]    # high pH (fully deprotonated)
    δ_1_guess = mean(delta_values)  # middle state

    # Initial guesses for pKa values
    if pKa_guess === nothing
        # Find inflection points in the curve
        pH_range = extrema(pH_values)
        pKa1_guess = pH_range[1] + (pH_range[2] - pH_range[1]) / 3
        pKa2_guess = pH_range[1] + 2 * (pH_range[2] - pH_range[1]) / 3
    else
        pKa1_guess, pKa2_guess = pKa_guess
    end

    # Ensure pKa1 < pKa2
    if pKa1_guess > pKa2_guess
        pKa1_guess, pKa2_guess = pKa2_guess, pKa1_guess
    end

    p0 = [pKa1_guess, pKa2_guess, δ_0_guess, δ_1_guess, δ_2_guess]

    # Fit with bounds to ensure pKa1 < pKa2
    try
        # Use lower/upper bounds
        lower = [0.0, 0.0, -Inf, -Inf, -Inf]
        upper = [14.0, 14.0, Inf, Inf, Inf]

        fit = curve_fit(henderson_hasselbalch_2pKa, pH_values, delta_values, p0;
                       lower=lower, upper=upper)

        params = coef(fit)
        errors = stderror(fit)

        # Ensure pKa1 < pKa2 in result
        if params[1] > params[2]
            # Swap pKa values and corresponding limiting shifts
            params[[1,2]] = params[[2,1]]
            errors[[1,2]] = errors[[2,1]]
            params[[3,5]] = params[[5,3]]
            errors[[3,5]] = errors[[5,3]]
        end

        # Calculate RMSD
        predicted = henderson_hasselbalch_2pKa(pH_values, params)
        rmsd = sqrt(mean((delta_values .- predicted) .^ 2))

        return IndicatorFit(
            indicator,
            nucleus,
            2,  # n_pKa
            [params[1] ± errors[1], params[2] ± errors[2]],  # pKa values
            [params[3] ± errors[3], params[4] ± errors[4], params[5] ± errors[5]],  # δ limits
            n,
            rmsd,
            extrema(pH_values)
        )
    catch e
        error("2-pKa fitting failed for $indicator: $e")
    end
end

"""
    fit_all_indicators(data; indicator_props=nothing) -> Dict{String, IndicatorFit}

Fit Henderson-Hasselbalch model to all indicators in the dataset.

Arguments:
    data: DataFrame from select_condition
    indicator_props: optional Dict{String, IndicatorProperties} specifying n_pKa for each indicator

Returns a dictionary mapping indicator names to IndicatorFit results.
"""
function fit_all_indicators(data; indicator_props=nothing)
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

        # Determine n_pKa from properties if provided
        n_pKa = 1
        pKa_guess = nothing
        if indicator_props !== nothing && haskey(indicator_props, ind)
            props = indicator_props[ind]
            n_pKa = props.n_pKa
            if !isempty(props.pKa_literature)
                pKa_guess = n_pKa == 1 ? props.pKa_literature[1] : props.pKa_literature
            end
        end

        try
            fit = fit_indicator(pH_vals, delta_vals;
                               indicator=ind, nucleus=nucleus,
                               n_pKa=n_pKa, pKa_guess=pKa_guess)
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
    println("\n" * "─"^80)

    # Separate 1-pKa and 2-pKa indicators
    single_pKa = filter(kv -> kv.second.n_pKa == 1, results)
    double_pKa = filter(kv -> kv.second.n_pKa == 2, results)

    if !isempty(single_pKa)
        println("Single-pKa Indicators:")
        println("─"^80)
        @printf("%-15s %3s %8s %12s %12s %5s %8s\n",
            "Indicator", "Nuc", "pKa", "δ_HA (ppm)", "δ_A (ppm)", "n", "RMSD")
        println("─"^80)

        for ind in sort(collect(keys(single_pKa)))
            r = single_pKa[ind]
            @printf("%-15s %3s %8.3f %12.4f %12.4f %5d %8.4f\n",
                r.indicator, r.nucleus,
                Measurements.value(r.pKa[1]),
                Measurements.value(r.δ_limits[1]),
                Measurements.value(r.δ_limits[2]),
                r.n_points, r.rmsd)
            @printf("%-15s     ±%-6.3f ±%-10.4f ±%-10.4f\n",
                "",
                Measurements.uncertainty(r.pKa[1]),
                Measurements.uncertainty(r.δ_limits[1]),
                Measurements.uncertainty(r.δ_limits[2]))
        end
    end

    if !isempty(double_pKa)
        println("\nTwo-pKa Indicators:")
        println("─"^80)
        @printf("%-12s %3s %7s %7s %10s %10s %10s %4s %7s\n",
            "Indicator", "Nuc", "pKa1", "pKa2", "δ_0", "δ_1", "δ_2", "n", "RMSD")
        println("─"^80)

        for ind in sort(collect(keys(double_pKa)))
            r = double_pKa[ind]
            @printf("%-12s %3s %7.3f %7.3f %10.4f %10.4f %10.4f %4d %7.4f\n",
                r.indicator, r.nucleus,
                Measurements.value(r.pKa[1]),
                Measurements.value(r.pKa[2]),
                Measurements.value(r.δ_limits[1]),
                Measurements.value(r.δ_limits[2]),
                Measurements.value(r.δ_limits[3]),
                r.n_points, r.rmsd)
        end
    end

    println("─"^80)
end
