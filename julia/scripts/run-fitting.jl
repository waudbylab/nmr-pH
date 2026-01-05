include("../src/PHIndicators.jl")
using .PHIndicators
##
# For 2-pKa indicators like maleate:
indicator_props = Dict(
    "HEPES_1" => IndicatorProperties("HEPES_1", "H", 2, 0, [3.5, 7.5], false),
    "HEPES_2" => IndicatorProperties("HEPES_2", "H", 2, 0, [3.5, 7.5], false),
    "HEPES_3" => IndicatorProperties("HEPES_3", "H", 2, 0, [3.5, 7.5], false),
    "HEPES_4" => IndicatorProperties("HEPES_4", "H", 2, 0, [3.5, 7.5], false),
    "piperazine" => IndicatorProperties("piperazine", "H", 2, 0, [5.6, 9.5], false),
)

results = run_analysis("data/group1-rough-h2o.csv";
    temperature=283.0,
    ionic_strength=0.0,
    pKa_reference=7.2,
    indicator_props=indicator_props
)

# Check calibration order
println(results.bootstrap.calibration_order)
