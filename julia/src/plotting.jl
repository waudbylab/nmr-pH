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
    plot_initial_fits(data, results; save_path=nothing)

Create multi-panel plots showing all indicator fits, organized by nucleus.
Saves separate files for each nucleus.
"""
function plot_initial_fits(data, results::Dict{String,IndicatorFit}; save_path=nothing)
    # Group indicators by nucleus
    by_nucleus = Dict{String,Vector{String}}()
    for (ind, fit) in results
        nuc = fit.nucleus
        if !haskey(by_nucleus, nuc)
            by_nucleus[nuc] = String[]
        end
        push!(by_nucleus[nuc], ind)
    end

    # Sort indicators within each nucleus group
    for nuc in keys(by_nucleus)
        sort!(by_nucleus[nuc])
    end

    # Create plots for each nucleus
    for (nuc, indicators) in sort(by_nucleus)
        n_ind = length(indicators)
        n_cols = min(4, n_ind)
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
                pKa_label = @sprintf("%.2f", Measurements.value(r.pKa[1]))
            else
                params = [Measurements.value(r.pKa[1]),
                    Measurements.value(r.pKa[2]),
                    Measurements.value(r.δ_limits[1]),
                    Measurements.value(r.δ_limits[2]),
                    Measurements.value(r.δ_limits[3])]
                δ_curve = henderson_hasselbalch_2pKa(collect(pH_curve), params)
                pKa_label = @sprintf("%.2f, %.2f",
                    Measurements.value(r.pKa[1]),
                    Measurements.value(r.pKa[2]))
            end

            p = plot(pH_curve, δ_curve,
                label="",
                linewidth=2, color=:blue,
                legend=false)
            scatter!(p, pH_data, δ_data,
                label="", markersize=4, markerstrokewidth=0, color=:red,
                markeralpha=0.6)

            xlabel!(p, "pH")
            ylabel!(p, "δ (ppm)")
            title!(p, "$(ind)\npKa = $pKa_label", titlefontsize=10)

            push!(plots, p)
        end

        # Combine into single figure
        fig = plot(plots...,
            layout=(n_rows, n_cols),
            size=(300 * n_cols, 250 * n_rows),
            margin=4Plots.mm,
            plot_title="Initial Fits: $(nuc) nucleus (electrode pH)",
            plot_titlefontsize=14)

        if !isnothing(save_path)
            # Create separate file for each nucleus
            base_path = splitext(save_path)[1]
            nuc_path = "$(base_path)_$(nuc).png"
            savefig(fig, nuc_path)
            println("Saved $(nuc) nucleus plots to $nuc_path")
        end
    end
end

"""
    plot_bootstrap_results(data, bootstrap_result; save_path=nothing)

Create diagnostic plots for bootstrap analysis.
Saves separate files for:
1. pH corrections (comparison + residuals)
2. Bootstrap fits by nucleus
"""
function plot_bootstrap_results(data, bootstrap_result; save_path=nothing)
    fits = bootstrap_result.fits
    pH_data = bootstrap_result.pH_data

    # -------------------------------------------------------------------------
    # Plot 1: pH comparison and corrections
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

    # Comparison plot
    p_compare = scatter(pH_nom, pH_boot,
        xlabel="Electrode pH", ylabel="Corrected pH",
        label="", markersize=5, markerstrokewidth=0,
        color=:blue, markeralpha=0.6,
        title="pH Calibration")

    # Add error bars
    for i in 1:length(pH_nom)
        plot!(p_compare, [pH_nom[i], pH_nom[i]],
              [pH_boot[i] - pH_err[i], pH_boot[i] + pH_err[i]],
              color=:blue, alpha=0.3, label="")
    end

    # Add 1:1 line
    pH_range = [minimum(pH_nom) - 0.5, maximum(pH_nom) + 0.5]
    plot!(p_compare, pH_range, pH_range, label="1:1",
          linestyle=:dash, color=:gray, linewidth=2)

    # Residuals plot
    residuals = pH_boot .- pH_nom

    p_resid = scatter(pH_nom, residuals,
        xlabel="Electrode pH", ylabel="ΔpH (corrected - electrode)",
        label="", markersize=5, markerstrokewidth=0,
        color=:red, markeralpha=0.6,
        title="pH Correction")
    hline!(p_resid, [0], linestyle=:dash, color=:gray, linewidth=2, label="")

    # Add ±0.1 reference lines
    hline!(p_resid, [-0.1, 0.1], linestyle=:dot, color=:lightgray, linewidth=1.5, label="")

    # Combine pH correction plots
    fig_pH = plot(p_compare, p_resid,
        layout=(1, 2),
        size=(900, 400),
        margin=5Plots.mm,
        plot_title="Bootstrap pH Correction",
        plot_titlefontsize=14)

    if !isnothing(save_path)
        base_path = splitext(save_path)[1]
        pH_path = "$(base_path)_pH_correction.png"
        savefig(fig_pH, pH_path)
        println("Saved pH correction plot to $pH_path")
    end

    # -------------------------------------------------------------------------
    # Plot 2: Individual indicator fits by nucleus
    # -------------------------------------------------------------------------

    # Use calibration order if available
    indicators = haskey(bootstrap_result, :calibration_order) ?
                 bootstrap_result.calibration_order :
                 sort(collect(keys(fits)))

    # Group by nucleus
    by_nucleus = Dict{String,Vector{String}}()
    for ind in indicators
        nuc = fits[ind].nucleus
        if !haskey(by_nucleus, nuc)
            by_nucleus[nuc] = String[]
        end
        push!(by_nucleus[nuc], ind)
    end

    # Create plots for each nucleus
    for (nuc, nuc_indicators) in sort(by_nucleus)
        n_ind = length(nuc_indicators)
        n_cols = min(4, n_ind)
        n_rows = ceil(Int, n_ind / n_cols)

        ind_plots = []

        for ind in nuc_indicators
            br = fits[ind]
            ind_data_raw = filter(row -> row.indicator == ind, data)

            # Get bootstrap pH values and uncertainties
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
                pKa_label = @sprintf("%.2f", Measurements.value(br.pKa_bootstrap[1]))
            else
                params = [Measurements.value(br.pKa_bootstrap[1]),
                    Measurements.value(br.pKa_bootstrap[2]),
                    Measurements.value(br.δ_limits[1]),
                    Measurements.value(br.δ_limits[2]),
                    Measurements.value(br.δ_limits[3])]
                δ_curve = henderson_hasselbalch_2pKa(collect(pH_curve), params)
                pKa_label = @sprintf("%.2f, %.2f",
                    Measurements.value(br.pKa_bootstrap[1]),
                    Measurements.value(br.pKa_bootstrap[2]))
            end

            p = plot(pH_curve, δ_curve,
                label="",
                linewidth=2, color=:blue,
                legend=false)
            scatter!(p, pH_boot_vals, δ_vals,
                label="", markersize=4, markerstrokewidth=0,
                color=:red, markeralpha=0.6)

            xlabel!(p, "pH")
            ylabel!(p, "δ (ppm)")

            tier_str = br.is_reference ? "REF" : "T$(br.calibration_tier)"
            title!(p, "$(ind) [$tier_str]\npKa = $pKa_label", titlefontsize=10)

            push!(ind_plots, p)
        end

        # Combine into single figure for this nucleus
        fig_nuc = plot(ind_plots...,
            layout=(n_rows, n_cols),
            size=(300 * n_cols, 250 * n_rows),
            margin=4Plots.mm,
            plot_title="Bootstrap Fits: $(nuc) nucleus (corrected pH)",
            plot_titlefontsize=14)

        if !isnothing(save_path)
            base_path = splitext(save_path)[1]
            nuc_path = "$(base_path)_fits_$(nuc).png"
            savefig(fig_nuc, nuc_path)
            println("Saved $(nuc) nucleus bootstrap fits to $nuc_path")
        end
    end
end
