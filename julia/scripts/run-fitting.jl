include("../src/PHIndicators.jl")
using .PHIndicators
##
# For 2-pKa indicators like maleate:
indicator_props = Dict(
    "HEPES1" => IndicatorProperties("HEPES1", "H", 2, 0, [1.9, 7.5], false)
)

results = run_analysis("../data/group1-rough-h2o.csv";
    temperature=283.0,
    ionic_strength=0.3,
    pKa_reference=6.8,
    indicator_props=indicator_props
)

# Check calibration order
println(results.bootstrap.calibration_order)
