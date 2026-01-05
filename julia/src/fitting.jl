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

    # Fit with bounds to constrain pKa
    try
        # Constrain pKa to ±1 unit of initial guess
        lower = [pKa_guess - 1.0, -Inf, -Inf]
        upper = [pKa_guess + 1.0, Inf, Inf]

        # Also constrain chemical shifts for 1H to ±1 ppm of initial guess
        if nucleus == "H"
            lower[2] = δ_HA_guess - 1.0
            upper[2] = δ_HA_guess + 1.0
            lower[3] = δ_A_guess - 1.0
            upper[3] = δ_A_guess + 1.0
        end

        fit = curve_fit(henderson_hasselbalch, pH_values, delta_values, p0;
                       lower=lower, upper=upper)

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

    # Fit with bounds to ensure pKa1 < pKa2 and constrain to initial guesses
    try
        # Constrain pKa values to ±1 unit of initial guesses
        lower = [max(0.0, pKa1_guess - 1.0), max(0.0, pKa2_guess - 1.0), -Inf, -Inf, -Inf]
        upper = [min(14.0, pKa1_guess + 1.0), min(14.0, pKa2_guess + 1.0), Inf, Inf, Inf]

        # Also constrain chemical shifts for 1H to ±1 ppm of initial guesses
        if nucleus == "H"
            lower[3] = δ_0_guess - 1.0
            upper[3] = δ_0_guess + 1.0
            lower[4] = δ_1_guess - 1.0
            upper[4] = δ_1_guess + 1.0
            lower[5] = δ_2_guess - 1.0
            upper[5] = δ_2_guess + 1.0
        end

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
    fit_grouped_indicators_2pKa(pH_values_list, delta_values_list, indicator_names, nucleus_list;
                                 pKa_guess=nothing) -> (pKa1, pKa2, δ_limits_list)

Fit multiple signals from the same indicator molecule with shared pKa values.

Arguments:
    pH_values_list: vector of pH arrays (one per signal)
    delta_values_list: vector of chemical shift arrays (one per signal)
    indicator_names: vector of indicator names
    nucleus_list: vector of nucleus types for each signal
    pKa_guess: initial guess for [pKa1, pKa2]

Returns:
    pKa1, pKa2: shared pKa values (as Measurements)
    δ_limits_list: vector of limiting shift arrays (one per signal)
"""
function fit_grouped_indicators_2pKa(pH_values_list, delta_values_list, indicator_names, nucleus_list;
                                     pKa_guess=nothing)
    n_signals = length(indicator_names)

    # Initial guesses for pKa values
    if pKa_guess === nothing
        # Use first signal to estimate pKa values
        pH_range = extrema(pH_values_list[1])
        pKa1_guess = pH_range[1] + (pH_range[2] - pH_range[1]) / 3
        pKa2_guess = pH_range[1] + 2 * (pH_range[2] - pH_range[1]) / 3
    else
        pKa1_guess, pKa2_guess = pKa_guess
    end

    # Ensure pKa1 < pKa2
    if pKa1_guess > pKa2_guess
        pKa1_guess, pKa2_guess = pKa2_guess, pKa1_guess
    end

    # Initial guesses for limiting shifts (3 per signal)
    δ_guesses = []
    for i in 1:n_signals
        order = sortperm(pH_values_list[i])
        δ_sorted = delta_values_list[i][order]
        push!(δ_guesses, [δ_sorted[1], mean(delta_values_list[i]), δ_sorted[end]])
    end

    # Build parameter vector: [pKa1, pKa2, δ_0_1, δ_1_1, δ_2_1, δ_0_2, δ_1_2, δ_2_2, ...]
    p0 = vcat([pKa1_guess, pKa2_guess], vcat(δ_guesses...)...)

    # Concatenate all data
    pH_all = vcat(pH_values_list...)
    δ_all = vcat(delta_values_list...)

    # Create signal index vector to track which data point belongs to which signal
    signal_idx = vcat([fill(i, length(pH_values_list[i])) for i in 1:n_signals]...)

    # Model function that uses shared pKa but signal-specific limiting shifts
    function grouped_model(pH, p)
        pKa1, pKa2 = p[1], p[2]
        result = similar(pH)

        for i in 1:n_signals
            # Extract limiting shifts for this signal
            δ_0 = p[2 + (i-1)*3 + 1]
            δ_1 = p[2 + (i-1)*3 + 2]
            δ_2 = p[2 + (i-1)*3 + 3]

            # Find data points for this signal
            mask = signal_idx .== i
            pH_i = pH[mask]

            # Calculate predictions
            α = @. 10^(pH_i - pKa1)
            β = @. 10^(2 * pH_i - pKa1 - pKa2)
            D = @. 1 + α + β
            result[mask] = @. (δ_0 + δ_1 * α + δ_2 * β) / D
        end

        return result
    end

    # Set up bounds
    lower = [max(0.0, pKa1_guess - 1.0), max(0.0, pKa2_guess - 1.0)]
    upper = [min(14.0, pKa1_guess + 1.0), min(14.0, pKa2_guess + 1.0)]

    # Add bounds for limiting shifts
    for i in 1:n_signals
        if nucleus_list[i] == "H"
            # Constrain 1H shifts to ±1 ppm
            for j in 1:3
                push!(lower, δ_guesses[i][j] - 1.0)
                push!(upper, δ_guesses[i][j] + 1.0)
            end
        else
            # No constraints for other nuclei
            append!(lower, [-Inf, -Inf, -Inf])
            append!(upper, [Inf, Inf, Inf])
        end
    end

    try
        fit = curve_fit(grouped_model, pH_all, δ_all, p0; lower=lower, upper=upper)
        params = coef(fit)
        errors = stderror(fit)

        # Extract results
        pKa1 = params[1] ± errors[1]
        pKa2 = params[2] ± errors[2]

        δ_limits_list = []
        for i in 1:n_signals
            δ_0 = params[2 + (i-1)*3 + 1] ± errors[2 + (i-1)*3 + 1]
            δ_1 = params[2 + (i-1)*3 + 2] ± errors[2 + (i-1)*3 + 2]
            δ_2 = params[2 + (i-1)*3 + 3] ± errors[2 + (i-1)*3 + 3]
            push!(δ_limits_list, [δ_0, δ_1, δ_2])
        end

        return pKa1, pKa2, δ_limits_list
    catch e
        error("Grouped fitting failed: $e")
    end
end

"""
    fit_all_indicators(data; indicator_props=nothing, indicator_groups=nothing) -> Dict{String, IndicatorFit}

Fit Henderson-Hasselbalch model to all indicators in the dataset.

Arguments:
    data: DataFrame from select_condition
    indicator_props: optional Dict{String, IndicatorProperties} specifying n_pKa for each indicator
    indicator_groups: optional Dict{String, Vector{String}} mapping group name to list of indicator names
                      that should share the same pKa values (e.g., "HEPES" => ["HEPES_1", "HEPES_2", "HEPES_3", "HEPES_4"])

Returns a dictionary mapping indicator names to IndicatorFit results.
"""
function fit_all_indicators(data; indicator_props=nothing, indicator_groups=nothing)
    results = Dict{String,IndicatorFit}()

    indicators = unique(data.indicator)

    # Track which indicators are part of groups
    grouped_indicators = Set{String}()
    if indicator_groups !== nothing
        for (group_name, group_members) in indicator_groups
            union!(grouped_indicators, group_members)
        end
    end

    println("\n" * "="^70)
    println("INITIAL HENDERSON-HASSELBALCH FITS (using electrode pH)")
    println("="^70)

    # First, fit grouped indicators
    if indicator_groups !== nothing
        for (group_name, group_members) in indicator_groups
            # Filter to only members that exist in the data
            available_members = filter(m -> m in indicators, group_members)

            if length(available_members) < 2
                @warn "Group $group_name has fewer than 2 members in data, fitting individually"
                continue
            end

            println("\nFitting grouped indicators: $(join(available_members, ", "))")

            # Collect data for all group members
            pH_values_list = []
            delta_values_list = []
            nucleus_list = []

            for ind in available_members
                subset = filter(row -> row.indicator == ind, data)
                push!(pH_values_list, Float64.(subset.nominal_pH))
                push!(delta_values_list, Float64.(subset.delta_obs))
                push!(nucleus_list, subset.nucleus[1])
            end

            # Get pKa guess from first member's properties
            pKa_guess = nothing
            n_pKa = 2  # Assume groups are for 2-pKa indicators
            if indicator_props !== nothing && haskey(indicator_props, available_members[1])
                props = indicator_props[available_members[1]]
                n_pKa = props.n_pKa
                if !isempty(props.pKa_literature)
                    pKa_guess = props.pKa_literature
                end
            end

            try
                if n_pKa == 2
                    pKa1, pKa2, δ_limits_list = fit_grouped_indicators_2pKa(
                        pH_values_list, delta_values_list, available_members, nucleus_list;
                        pKa_guess=pKa_guess
                    )

                    # Create individual IndicatorFit objects for each member
                    for (i, ind) in enumerate(available_members)
                        n_points = length(pH_values_list[i])
                        pH_range = extrema(pH_values_list[i])

                        # Calculate RMSD for this signal
                        params = [Measurements.value(pKa1), Measurements.value(pKa2),
                                 Measurements.value(δ_limits_list[i][1]),
                                 Measurements.value(δ_limits_list[i][2]),
                                 Measurements.value(δ_limits_list[i][3])]
                        predicted = henderson_hasselbalch_2pKa(pH_values_list[i], params)
                        rmsd = sqrt(mean((delta_values_list[i] .- predicted) .^ 2))

                        results[ind] = IndicatorFit(
                            ind,
                            nucleus_list[i],
                            2,
                            [pKa1, pKa2],
                            δ_limits_list[i],
                            n_points,
                            rmsd,
                            pH_range
                        )
                    end
                else
                    @warn "Grouped fitting only implemented for 2-pKa indicators, fitting individually"
                end
            catch e
                @warn "Grouped fitting failed for $group_name: $e"
                @warn "Fitting members individually instead"
                # Fall back to individual fitting
                setdiff!(grouped_indicators, available_members)
            end
        end
    end

    # Then fit ungrouped indicators individually
    for ind in sort(indicators)
        # Skip if already fitted as part of a group
        if ind in grouped_indicators
            continue
        end

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
