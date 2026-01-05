"""
PHIndicators - pH Indicator Calibration from NMR Chemical Shifts

This module provides functions to:
1. Load and organise chemical shift data from CSV
2. Fit Henderson-Hasselbalch models (1-pKa and 2-pKa) to determine pKa and limiting shifts
3. Bootstrap a refined pH scale using an iterated algorithm
4. Produce diagnostic plots and summary tables

Usage:
    using PHIndicators

    # Load data
    data = load_indicator_data("chemical_shifts.csv")

    # Select one (T, I) condition
    subset = select_condition(data, 298.0, 0.15)

    # Initial fits using electrode pH
    initial_fits = fit_all_indicators(subset)
    plot_initial_fits(subset, initial_fits)

    # Iterated bootstrap using literature phosphate pKa
    result = iterated_bootstrap_pH(subset, initial_fits;
        reference_indicator="phosphate",
        pKa_reference=6.7)
    plot_bootstrap_results(subset, result)
"""
module PHIndicators

using CSV
using DataFrames
using LsqFit
using Measurements
using Plots
using Printf
using Statistics

# Include submodules
include("types.jl")
include("models.jl")
include("data.jl")
include("fitting.jl")
include("bootstrap.jl")
include("plotting.jl")

# Export types
export IndicatorProperties, IndicatorFit, BootstrapResult

# Export model functions
export henderson_hasselbalch, henderson_hasselbalch_2pKa
export hh_fixed_limits, hh_2pKa_fixed_limits
export pH_from_shift, pH_from_shift_2pKa, is_reliable
export population_fractions_2pKa

# Export data functions
export load_indicator_data, select_condition, list_conditions
export get_indicator_data

# Export fitting functions
export fit_indicator, fit_all_indicators
export get_pKa, get_δ_HA, get_δ_A, get_δ_mid

# Export bootstrap functions
export iterated_bootstrap_pH, bootstrap_pH

# Export plotting functions
export plot_initial_fits, plot_bootstrap_results, print_bootstrap_table

# ============================================================================
# CONVENIENCE FUNCTION: RUN COMPLETE ANALYSIS
# ============================================================================

"""
    run_analysis(filename; temperature, ionic_strength, pKa_reference, ...)

Run complete analysis pipeline on a CSV file.

Arguments:
    filename: path to CSV file with chemical shift data
    temperature: target temperature in K
    ionic_strength: ionic strength in M
    pKa_reference: literature pKa for reference indicator at this (T, I) condition
    reference_indicator: name of reference indicator (default "phosphate")
    indicator_props: optional Dict{String,IndicatorProperties} for multi-pKa indicators
    save_plots: if true, save plots to PNG files

Returns a NamedTuple with all results.
"""
function run_analysis(filename;
    temperature,
    ionic_strength,
    pKa_reference,
    reference_indicator="phosphate",
    indicator_props=nothing,
    save_plots=true)

    println("\n" * "="^70)
    println("pH INDICATOR ANALYSIS")
    println("="^70)
    println("File: $filename")
    println("Condition: T = $temperature K, I = $ionic_strength M")
    println("Reference: $reference_indicator with pKa = $pKa_reference")
    println("="^70)

    # Load data
    println("\nLoading data...")
    data = load_indicator_data(filename)
    list_conditions(data)

    # Select condition
    subset = select_condition(data, temperature, ionic_strength)

    # Initial fits
    initial_fits = fit_all_indicators(subset; indicator_props=indicator_props)

    if save_plots
        plot_initial_fits(subset, initial_fits; save_path="initial_fits.png")
    end

    # Iterated bootstrap
    bootstrap_result = iterated_bootstrap_pH(subset, initial_fits;
        reference_indicator=reference_indicator,
        pKa_reference=pKa_reference)

    # Print parameter summary table
    print_bootstrap_table(bootstrap_result)

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

export run_analysis

end # module
