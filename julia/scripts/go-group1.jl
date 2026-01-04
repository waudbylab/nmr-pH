#=
Example: pH Indicator Analysis
==============================

This script demonstrates how to use the pH indicator analysis module.

Before running, ensure you have the required packages:
    using Pkg
    Pkg.add(["CSV", "DataFrames", "LsqFit", "Measurements", "Plots", "Printf"])
=#

# Include the analysis module
include("../src/pH_indicator_analysis.jl")

# ============================================================================
# OPTION 1: Run complete analysis in one step
# ============================================================================

# This is the simplest way to run the analysis
results = run_analysis(
    "../data/group1-rough-h2o.csv";
    temperature=298.0,
    ionic_strength=0.3,
    pKa_phosphate_lit=6.7,  # Literature value at this T, I
    save_plots=true
)

# Access individual results
println("\n\nAccessing results programmatically:")
println("─"^40)

# Initial fit for acetate
acetate_init = results.initial_fits["acetate"]
println("Acetate pKa (electrode pH): $(acetate_init.pKa)")

# Bootstrap fit for acetate
acetate_boot = results.bootstrap.fits["acetate"]
println("Acetate pKa (bootstrap):    $(acetate_boot.pKa_bootstrap)")

# pH data
println("\nFirst few rows of pH data:")
println(first(results.bootstrap.pH_data, 5))


# ============================================================================
# OPTION 2: Run analysis step by step
# ============================================================================

#=
# Step 1: Load data
data = load_indicator_data("your_data.csv")

# See what conditions are available
list_conditions(data)

# Step 2: Select one (T, I) condition
subset = select_condition(data, temperature=298.0, ionic_strength="low")

# Step 3: Initial fits using electrode pH
initial_fits = fit_all_indicators(subset)

# Plot initial fits
plot_initial_fits(subset, initial_fits; save_path="my_initial_fits.png")

# Step 4: Bootstrap using literature phosphate pKa
# The literature pKa should be corrected for your specific T and I
# See the ionic strength explanation document for how to calculate this
pKa_phosphate_lit = 7.20  # Example value - adjust for your conditions!

bootstrap_result = bootstrap_pH(
    subset, 
    initial_fits; 
    pKa_phosphate_lit = pKa_phosphate_lit,
    σ_δ_phosphate = 0.01,  # ³¹P measurement uncertainty (ppm)
    σ_δ_other = 0.005      # ¹H, ¹⁹F measurement uncertainty (ppm)
)

# Plot bootstrap results
plot_bootstrap_results(subset, bootstrap_result; save_path="my_bootstrap.png")

# Access the refined pH values
pH_data = bootstrap_result.pH_data
println(pH_data)

# Access the refined pKa values
for (ind, fit) in bootstrap_result.fits
    println("$ind: pKa = $(fit.pKa_bootstrap)")
end
=#


# ============================================================================
# NOTES ON INTERPRETING OUTPUT
# ============================================================================

#=
The analysis produces several outputs:

1. INITIAL FITS (using electrode pH)
   - pKa, δ_HA, δ_A for each indicator
   - These use the nominal pH from the electrode as the independent variable
   - Good for seeing the overall shape of the data

2. BOOTSTRAP FITS (using phosphate-derived pH)
   - Phosphate pKa is FIXED to the literature value
   - Other indicators are REFITTED using phosphate-derived pH
   - The pKa values may shift slightly from initial fits

3. pH COMPARISON TABLE
   - Shows nominal (electrode) pH vs bootstrap (combined) pH
   - The difference (Δ) shows systematic offset between electrode and NMR-derived pH
   - At 298 K these should agree well; at other temperatures expect systematic shifts

4. PLOTS
   - initial_fits.png: Each indicator plotted vs electrode pH
   - bootstrap_results.png: 
     - Top left: electrode pH vs bootstrap pH (should be ~1:1)
     - Top right: residuals (should be small and random)
     - Bottom: Each indicator plotted vs bootstrap pH

WHAT TO CHECK:
- Initial fits should capture the sigmoidal shape well
- Bootstrap pKa values should be similar to initial (within ~0.1 units)
- pH comparison should show good correlation (R² > 0.99)
- Residuals should be small (<0.1 pH units) and not show systematic patterns
=#