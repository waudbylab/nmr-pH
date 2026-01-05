# Plotting functions for pH indicator analysis

"""
    print_bootstrap_table(bootstrap_result)

Print a nicely formatted summary table of indicator pKas and chemical shifts after bootstrapping.
"""
function print_bootstrap_table(bootstrap_result)
    fits = bootstrap_result.fits
    calibration_order = haskey(bootstrap_result, :calibration_order) ?
                        bootstrap_result.calibration_order :
                        sort(collect(keys(fits)))

    println("\n" * "="^100)
    println("BOOTSTRAP PARAMETER SUMMARY")
    println("="^100)

    # Header for single-pKa indicators
    println("\n─── Single-pKa Indicators ───")
    println()
    @printf("%-15s %4s %4s  %12s  %12s  %12s  %12s\n",
        "Indicator", "Nuc", "Tier", "pKa", "δ(HA) ppm", "δ(A⁻) ppm", "Δδ ppm")
    println("─"^85)

    for ind in calibration_order
        br = fits[ind]
        br.n_pKa == 1 || continue

        tier_str = br.is_reference ? "REF" : string(br.calibration_tier)

        pKa_val = Measurements.value(br.pKa_bootstrap[1])
        pKa_err = Measurements.uncertainty(br.pKa_bootstrap[1])

        δ_HA_val = Measurements.value(br.δ_limits[1])
        δ_HA_err = Measurements.uncertainty(br.δ_limits[1])

        δ_A_val = Measurements.value(br.δ_limits[2])
        δ_A_err = Measurements.uncertainty(br.δ_limits[2])

        Δδ = abs(δ_A_val - δ_HA_val)

        @printf("%-15s %4s %4s  %5.3f ± %4.3f  %6.3f ± %4.3f  %6.3f ± %4.3f  %6.3f\n",
            br.indicator, br.nucleus, tier_str,
            pKa_val, pKa_err,
            δ_HA_val, δ_HA_err,
            δ_A_val, δ_A_err,
            Δδ)
    end

    # Check if there are any 2-pKa indicators
    has_2pKa = any(fits[ind].n_pKa == 2 for ind in calibration_order)

    if has_2pKa
        println("\n─── Two-pKa Indicators ───")
        println()
        @printf("%-15s %4s %4s  %12s  %12s  %12s  %12s  %12s\n",
            "Indicator", "Nuc", "Tier", "pKa₁", "pKa₂", "δ(H₂A) ppm", "δ(HA⁻) ppm", "δ(A²⁻) ppm")
        println("─"^105)

        for ind in calibration_order
            br = fits[ind]
            br.n_pKa == 2 || continue

            tier_str = br.is_reference ? "REF" : string(br.calibration_tier)

            pKa1_val = Measurements.value(br.pKa_bootstrap[1])
            pKa1_err = Measurements.uncertainty(br.pKa_bootstrap[1])
            pKa2_val = Measurements.value(br.pKa_bootstrap[2])
            pKa2_err = Measurements.uncertainty(br.pKa_bootstrap[2])

            δ_0_val = Measurements.value(br.δ_limits[1])
            δ_0_err = Measurements.uncertainty(br.δ_limits[1])
            δ_1_val = Measurements.value(br.δ_limits[2])
            δ_1_err = Measurements.uncertainty(br.δ_limits[2])
            δ_2_val = Measurements.value(br.δ_limits[3])
            δ_2_err = Measurements.uncertainty(br.δ_limits[3])

            @printf("%-15s %4s %4s  %5.3f ± %4.3f  %5.3f ± %4.3f  %6.3f ± %4.3f  %6.3f ± %4.3f  %6.3f ± %4.3f\n",
                br.indicator, br.nucleus, tier_str,
                pKa1_val, pKa1_err,
                pKa2_val, pKa2_err,
                δ_0_val, δ_0_err,
                δ_1_val, δ_1_err,
                δ_2_val, δ_2_err)
        end
    end

    println("="^100)
    println()
end

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

        if r.n_pKa == 1
            params = [Measurements.value(r.pKa[1]),
                Measurements.value(r.δ_limits[1]),
                Measurements.value(r.δ_limits[2])]
            δ_curve = henderson_hasselbalch(collect(pH_curve), params)
            pKa_label = "pKa=$(round(Measurements.value(r.pKa[1]), digits=2))"
        else
            params = [Measurements.value(r.pKa[1]),
                Measurements.value(r.pKa[2]),
                Measurements.value(r.δ_limits[1]),
                Measurements.value(r.δ_limits[2]),
                Measurements.value(r.δ_limits[3])]
            δ_curve = henderson_hasselbalch_2pKa(collect(pH_curve), params)
            pKa_label = "pKa=$(round(Measurements.value(r.pKa[1]), digits=2)), $(round(Measurements.value(r.pKa[2]), digits=2))"
        end

        p = plot(pH_curve, δ_curve,
            label="Fit ($pKa_label)",
            linewidth=2, color=:blue)
        scatter!(p, pH_data, δ_data,
            label="Data", markersize=5, markerstrokewidth=0, color=:red)

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

    p_compare = scatter(pH_nom, pH_boot .± pH_err,
        xlabel="pH (electrode)", ylabel="pH (bootstrap)",
        label="Data", markersize=3, markerstrokewidth=0.5,
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
        label=nothing, markersize=6, markerstrokewidth=0,
        title="pH Residuals")
    hline!(p_resid, [0], linestyle=:dash, color=:gray, label=nothing)

    # Add ±0.1 reference lines
    hline!(p_resid, [-0.1, 0.1], linestyle=:dot, color=:lightgray, label=nothing)

    # -------------------------------------------------------------------------
    # Plot 3: Individual indicator fits using bootstrap pH
    # -------------------------------------------------------------------------

    # Use calibration order if available
    indicators = haskey(bootstrap_result, :calibration_order) ?
                 bootstrap_result.calibration_order :
                 sort(collect(keys(fits)))

    n_ind = length(indicators)
    n_cols = min(3, n_ind)
    n_rows = ceil(Int, n_ind / n_cols)

    ind_plots = []

    for ind in indicators
        br = fits[ind]
        ind_data_raw = filter(row -> row.indicator == ind, data)

        # Get bootstrap pH values and uncertainties for this indicator's data points
        pH_boot_vals = Float64[]
        pH_boot_errs = Float64[]
        δ_vals = Float64[]

        for row in eachrow(ind_data_raw)
            pH_row = filter(r -> r.nominal_pH == row.nominal_pH, pH_data)
            if nrow(pH_row) > 0 && !ismissing(pH_row.pH_combined[1])
                push!(pH_boot_vals, Measurements.value(pH_row.pH_combined[1]))
                push!(pH_boot_errs, Measurements.uncertainty(pH_row.pH_combined[1]))
                push!(δ_vals, row.delta_obs)
            end
        end

        if length(pH_boot_vals) < 2
            continue
        end

        # Generate smooth curve
        pH_curve = range(minimum(pH_boot_vals) - 0.3, maximum(pH_boot_vals) + 0.3, length=100)

        if br.n_pKa == 1
            params = [Measurements.value(br.pKa_bootstrap[1]),
                Measurements.value(br.δ_limits[1]),
                Measurements.value(br.δ_limits[2])]
            δ_curve = henderson_hasselbalch(collect(pH_curve), params)
            pKa_label = "pKa=$(round(Measurements.value(br.pKa_bootstrap[1]), digits=2))"
        else
            params = [Measurements.value(br.pKa_bootstrap[1]),
                Measurements.value(br.pKa_bootstrap[2]),
                Measurements.value(br.δ_limits[1]),
                Measurements.value(br.δ_limits[2]),
                Measurements.value(br.δ_limits[3])]
            δ_curve = henderson_hasselbalch_2pKa(collect(pH_curve), params)
            pKa_label = "pKa=$(round(Measurements.value(br.pKa_bootstrap[1]), digits=2)), $(round(Measurements.value(br.pKa_bootstrap[2]), digits=2))"
        end

        p = plot(pH_curve, δ_curve,
            label=pKa_label,
            linewidth=2, color=:blue)
        scatter!(p, pH_boot_vals .± pH_boot_errs, δ_vals,
            label="Data", markersize=3, markerstrokewidth=0.5, color=:red)

        xlabel!(p, "pH (bootstrap)")
        ylabel!(p, "δ (ppm)")

        tier_str = br.is_reference ? "REF" : "T$(br.calibration_tier)"
        title!(p, "$(br.indicator) [$tier_str]")

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
