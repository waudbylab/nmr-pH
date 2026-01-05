using DataFrames, CSV, Printf

# Load the data
data = CSV.read("julia/data/group1-rough-h2o.csv", DataFrame)

# Look at the high pH data points
high_ph_data = filter(row -> row.nominal_pH >= 9.0, data)

println("="^80)
println("HIGH pH DATA ANALYSIS")
println("="^80)

# From the bootstrap results (user provided)
indicators = Dict(
    "phosphate" => (pKa_boot=6.800, δ_HA=0.1028, δ_A=2.5618),
    "Tris" => (pKa_boot=8.731, δ_HA=3.7274, δ_A=3.5052),
    "imidazole1" => (pKa_boot=7.634, δ_HA=8.6865, δ_A=7.7579),
    "imidazole2" => (pKa_boot=7.633, δ_HA=7.4704, δ_A=7.1188),
)

column_map = Dict(
    "phosphate" => "P_phosphate",
    "Tris" => "H_Tris",
    "imidazole1" => "H_imidazole1",
    "imidazole2" => "H_imidazole2",
)

for row in eachrow(high_ph_data)
    println("\nNominal pH = $(row.nominal_pH), T = $(row.temperature) K")
    println("-"^80)
    @printf("%-15s %8s %8s %10s %12s %10s\n",
            "Indicator", "δ_obs", "δ_A", "Diff", "% to limit", "Valid?")
    println("-"^80)

    for (name, params) in indicators
        col = column_map[name]
        δ_obs = row[col]
        δ_HA = params.δ_HA
        δ_A = params.δ_A

        # Calculate expected shift at this pH using HH equation
        pKa = params.pKa_boot
        pH = row.nominal_pH
        frac_deprot = 1 / (1 + 10^(pKa - pH))
        δ_expected = δ_HA + (δ_A - δ_HA) * frac_deprot

        # Check if within valid range
        δ_min, δ_max = minmax(δ_HA, δ_A)
        valid = (δ_obs > δ_min && δ_obs < δ_max)

        # Calculate how close to limit (as percentage)
        if δ_A > δ_HA  # normal case
            pct_to_limit = (δ_obs - δ_HA) / (δ_A - δ_HA) * 100
        else  # inverted
            pct_to_limit = (δ_HA - δ_obs) / (δ_HA - δ_A) * 100
        end

        diff = δ_obs - δ_A

        @printf("%-15s %8.3f %8.3f %+10.4f %11.1f%% %10s\n",
                name, δ_obs, δ_A, diff, pct_to_limit, valid ? "YES" : "NO")
    end
end

println("\n" * "="^80)
println("ANALYSIS")
println("="^80)
println("At high pH, several indicators have observed shifts BEYOND their limiting")
println("shifts (δ_A), which causes them to be excluded from pH calculation.")
println("\nThis happens because:")
println("  1. Small measurement uncertainties near limits")
println("  2. Temperature/ionic strength effects not fully modeled")
println("  3. Additional ionization equilibria at extreme pH")
println("\nWhen indicators are excluded, fewer contribute to combined pH,")
println("leading to less reliable and potentially biased pH estimates.")
