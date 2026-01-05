# Type definitions for pH indicator analysis

"""
Properties of an indicator molecule.

Fields:
    name: indicator name (must match column name suffix in data)
    nucleus: NMR nucleus (P, F, H)
    n_pKa: number of pKa values (1 or 2)
    z_acid: charge of most protonated form
    pKa_literature: literature pKa values (starting guesses for fitting)
    is_reference: if true, pKa is fixed to literature value
"""
struct IndicatorProperties
    name::String
    nucleus::String
    n_pKa::Int
    z_acid::Int
    pKa_literature::Vector{Float64}
    is_reference::Bool
end

# Convenience constructor for single-pKa indicators
function IndicatorProperties(name::String, nucleus::String;
                            z_acid::Int=0,
                            pKa_literature::Float64=7.0,
                            is_reference::Bool=false)
    IndicatorProperties(name, nucleus, 1, z_acid, [pKa_literature], is_reference)
end

"""
Result structure for a single indicator fit (supports 1 or 2 pKa values).
"""
struct IndicatorFit
    indicator::String
    nucleus::String
    n_pKa::Int
    pKa::Vector{Measurement{Float64}}           # pKa values (1 or 2 elements)
    δ_limits::Vector{Measurement{Float64}}      # limiting shifts (2 or 3 elements)
    n_points::Int
    rmsd::Float64
    pH_range::Tuple{Float64,Float64}
end

# Convenience constructor for single-pKa fits
function IndicatorFit(indicator::String, nucleus::String,
                     pKa::Measurement{Float64},
                     δ_HA::Measurement{Float64},
                     δ_A::Measurement{Float64},
                     n_points::Int, rmsd::Float64,
                     pH_range::Tuple{Float64,Float64})
    IndicatorFit(indicator, nucleus, 1, [pKa], [δ_HA, δ_A], n_points, rmsd, pH_range)
end

# Access helpers for IndicatorFit
function get_pKa(fit::IndicatorFit, idx::Int=1)
    return fit.pKa[idx]
end

function get_δ_HA(fit::IndicatorFit)
    return fit.δ_limits[1]
end

function get_δ_A(fit::IndicatorFit)
    return fit.δ_limits[end]
end

function get_δ_mid(fit::IndicatorFit)
    fit.n_pKa == 2 || error("get_δ_mid only valid for 2-pKa indicators")
    return fit.δ_limits[2]
end

"""
Result structure for bootstrap analysis (supports 1 or 2 pKa values).
"""
struct BootstrapResult
    indicator::String
    nucleus::String
    n_pKa::Int
    pKa_initial::Vector{Measurement{Float64}}      # from electrode pH
    pKa_bootstrap::Vector{Measurement{Float64}}    # from bootstrap
    δ_limits::Vector{Measurement{Float64}}         # limiting shifts
    n_points_used::Int                             # points with reliable pH_true
    is_reference::Bool                             # true if pKa fixed to literature
    calibration_tier::Int                          # order in iterated bootstrap
end

# Convenience constructor for single-pKa bootstrap results
function BootstrapResult(indicator::String, nucleus::String,
                        pKa_initial::Measurement{Float64},
                        pKa_bootstrap::Measurement{Float64},
                        δ_HA::Measurement{Float64},
                        δ_A::Measurement{Float64},
                        n_points_used::Int, is_reference::Bool;
                        calibration_tier::Int=0)
    BootstrapResult(indicator, nucleus, 1,
                   [pKa_initial], [pKa_bootstrap], [δ_HA, δ_A],
                   n_points_used, is_reference, calibration_tier)
end
