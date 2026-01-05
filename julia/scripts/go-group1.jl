#=
Example: pH Indicator Analysis
==============================

This script demonstrates how to use the PHIndicators module.

Before running, ensure you have the required packages:
    using Pkg
    Pkg.add(["CSV", "DataFrames", "LsqFit", "Measurements", "Plots", "Printf"])
=#

# Add the src directory to the load path
push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))

using PHIndicators

# ============================================================================
# OPTION 1: Run complete analysis in one step (using iterated bootstrap)
# ============================================================================

# Define indicator properties for multi-pKa indicators
# Maleate has two pKa values in the measured range
indicator_props = Dict(
    # Example of a 2-pKa indicator (uncomment to use)
    # "maleate" => IndicatorProperties("maleate", "H", 2, 0, [1.9, 6.2], false)
)

# Run analysis with default single-pKa fitting
results = run_analysis(
    joinpath(@__DIR__, "..", "data", "group1-rough-h2o.csv");
    temperature=298.0,
    ionic_strength=0.3,
    pKa_reference=6.7,  # Literature value for phosphate at this T, I
    reference_indicator="phosphate",
    # indicator_props=indicator_props,  # Uncomment to use multi-pKa fitting
    save_plots=true
)

# Access individual results
println("\n\nAccessing results programmatically:")
println("─"^40)

# Initial fit for acetate
acetate_init = results.initial_fits["acetate"]
println("Acetate pKa (electrode pH): $(get_pKa(acetate_init))")

# Bootstrap fit for acetate
acetate_boot = results.bootstrap.fits["acetate"]
println("Acetate pKa (bootstrap):    $(acetate_boot.pKa_bootstrap[1])")

# Show calibration order
println("\nCalibration order: ", join(results.bootstrap.calibration_order, " → "))

# pH data
println("\nFirst few rows of pH data:")
println(first(results.bootstrap.pH_data, 5))


# ============================================================================
# OPTION 2: Run analysis step by step
# ============================================================================

#=
# Step 1: Load data
data = load_indicator_data(joinpath(@__DIR__, "..", "data", "group1-rough-h2o.csv"))

# See what conditions are available
list_conditions(data)

# Step 2: Select one (T, I) condition
subset = select_condition(data, 298.0, 0.3)

# Step 3: Initial fits using electrode pH
# For 2-pKa indicators, provide indicator_props
initial_fits = fit_all_indicators(subset)

# Plot initial fits
plot_initial_fits(subset, initial_fits; save_path="my_initial_fits.png")

# Step 4: Iterated bootstrap using literature phosphate pKa
# The algorithm starts with phosphate, then calibrates indicators
# in order of proximity to already-calibrated pKa values
bootstrap_result = iterated_bootstrap_pH(
    subset,
    initial_fits;
    reference_indicator = "phosphate",
    pKa_reference = 6.7,      # Literature value at this T, I
    σ_δ_reference = 0.01,     # ³¹P measurement uncertainty (ppm)
    σ_δ_other = 0.005         # ¹H, ¹⁹F measurement uncertainty (ppm)
)

# Plot bootstrap results
plot_bootstrap_results(subset, bootstrap_result; save_path="my_bootstrap.png")

# Access the refined pH values
pH_data = bootstrap_result.pH_data
println(pH_data)

# Access the refined pKa values
for ind in bootstrap_result.calibration_order
    fit = bootstrap_result.fits[ind]
    tier = fit.is_reference ? "REF" : "T$(fit.calibration_tier)"
    println("[$tier] $ind: pKa = $(fit.pKa_bootstrap[1])")
end
=#


# ============================================================================
# OPTION 3: Using 2-pKa indicators (e.g., maleate)
# ============================================================================

#=
# Define indicator properties with n_pKa=2 for diprotic indicators
indicator_props = Dict(
    "maleate" => IndicatorProperties(
        "maleate",      # name (must match column name)
        "H",            # nucleus
        2,              # n_pKa
        0,              # z_acid (charge of most protonated form)
        [1.9, 6.2],     # literature pKa values (initial guesses)
        false           # is_reference
    )
)

# Run analysis with 2-pKa fitting enabled
results_2pKa = run_analysis(
    joinpath(@__DIR__, "..", "data", "group1-rough-h2o.csv");
    temperature=298.0,
    ionic_strength=0.3,
    pKa_reference=6.7,
    indicator_props=indicator_props,
    save_plots=true
)
=#


# ============================================================================
# NOTES ON THE ITERATED BOOTSTRAP ALGORITHM
# ============================================================================

#=
The iterated bootstrap algorithm improves pH calibration by:

1. Starting with a reference indicator (phosphate) with known literature pKa
2. Calculating pH values using only the reference indicator
3. Finding the next uncalibrated indicator with estimated pKa closest to
   any already-calibrated indicator
4. Re-fitting that indicator's pKa using the current calibrated pH values
5. Updating combined pH estimates using inverse-variance weighting
6. Repeating steps 3-5 until all indicators are calibrated

This approach ensures that:
- Indicators are calibrated using pH values from the most reliable sources
- pH estimates improve iteratively as more indicators are calibrated
- The calibration "spreads out" from the reference in pH space

The output shows:
- Calibration order: which indicators were calibrated in which order
- Tier numbers: 0=reference, 1=first calibrated, etc.
- pKa changes: difference between electrode-based and bootstrap pKa
=#
