# Henderson-Hasselbalch models for pH indicator analysis

"""
Henderson-Hasselbalch model for chemical shift vs pH (single pKa).

    δ_obs = δ_HA + (δ_A - δ_HA) / (1 + 10^(pKa - pH))

Parameters:
    p[1] = pKa
    p[2] = δ_HA (limiting shift of acid/protonated form)
    p[3] = δ_A (limiting shift of base/deprotonated form)
"""
function henderson_hasselbalch(pH, p)
    pKa, δ_HA, δ_A = p
    return @. δ_HA + (δ_A - δ_HA) / (1 + 10^(pKa - pH))
end

"""
Henderson-Hasselbalch model with fixed limiting shifts (for bootstrap fitting).
Only pKa is fitted.
"""
function hh_fixed_limits(pH, p, δ_HA, δ_A)
    pKa = p[1]
    return @. δ_HA + (δ_A - δ_HA) / (1 + 10^(pKa - pH))
end

"""
Two-pKa Henderson-Hasselbalch model for diprotic indicators.

For a diprotic species H2A ⇌ HA⁻ ⇌ A²⁻:

    D = 1 + 10^(pH - pKa1) + 10^(2·pH - pKa1 - pKa2)
    δ_obs = (δ_0 + δ_1·10^(pH-pKa1) + δ_2·10^(2·pH-pKa1-pKa2)) / D

Parameters:
    p[1] = pKa1 (first ionisation)
    p[2] = pKa2 (second ionisation)
    p[3] = δ_0 (limiting shift of fully protonated form H2A)
    p[4] = δ_1 (limiting shift of singly deprotonated form HA⁻)
    p[5] = δ_2 (limiting shift of fully deprotonated form A²⁻)
"""
function henderson_hasselbalch_2pKa(pH, p)
    pKa1, pKa2, δ_0, δ_1, δ_2 = p
    α = @. 10^(pH - pKa1)
    β = @. 10^(2 * pH - pKa1 - pKa2)
    D = @. 1 + α + β
    return @. (δ_0 + δ_1 * α + δ_2 * β) / D
end

"""
Two-pKa model with fixed limiting shifts (for bootstrap fitting).
Only pKa1 and pKa2 are fitted.
"""
function hh_2pKa_fixed_limits(pH, p, δ_0, δ_1, δ_2)
    pKa1, pKa2 = p
    α = @. 10^(pH - pKa1)
    β = @. 10^(2 * pH - pKa1 - pKa2)
    D = @. 1 + α + β
    return @. (δ_0 + δ_1 * α + δ_2 * β) / D
end

"""
Calculate population fractions for a diprotic species.

Returns (f_0, f_1, f_2) where:
    f_0 = fraction in fully protonated form H2A
    f_1 = fraction in singly deprotonated form HA⁻
    f_2 = fraction in fully deprotonated form A²⁻
"""
function population_fractions_2pKa(pH, pKa1, pKa2)
    α = 10^(pH - pKa1)
    β = 10^(2 * pH - pKa1 - pKa2)
    D = 1 + α + β
    return (1 / D, α / D, β / D)
end

"""
Calculate pH from observed chemical shift for a single-pKa indicator.

Returns pH as a Measurement (value ± uncertainty) or missing if outside valid range.

Arguments:
    δ_obs: observed chemical shift (can be Measurement or Number)
    pKa: pKa value (can be Measurement or Number)
    δ_HA: limiting shift of acid form
    δ_A: limiting shift of base form
    σ_δ: measurement uncertainty in δ_obs (used if δ_obs is not a Measurement)
"""
function pH_from_shift(δ_obs, pKa, δ_HA, δ_A; σ_δ=0.01)
    # Extract values for range checking
    δ_obs_val = Measurements.value(δ_obs)
    δ_HA_val = Measurements.value(δ_HA)
    δ_A_val = Measurements.value(δ_A)

    # Check if δ_obs is between the limits
    δ_min, δ_max = minmax(δ_HA_val, δ_A_val)
    if δ_obs_val <= δ_min || δ_obs_val >= δ_max
        return missing
    end

    # Convert to Measurements if needed
    δ_obs_m = δ_obs isa Measurement ? δ_obs : δ_obs ± σ_δ
    pKa_m = pKa isa Measurement ? pKa : pKa ± 0.0
    δ_HA_m = δ_HA isa Measurement ? δ_HA : δ_HA ± 0.0
    δ_A_m = δ_A isa Measurement ? δ_A : δ_A ± 0.0

    # Calculate pH with uncertainty propagation
    ratio = (δ_obs_m - δ_HA_m) / (δ_A_m - δ_obs_m)
    return pKa_m + log10(ratio)
end

"""
Calculate pH from observed chemical shift for a two-pKa indicator.

Uses bisection root-finding since the two-pKa equation cannot be inverted analytically.

Returns pH as a Measurement (value ± uncertainty) or missing if no valid solution.
"""
function pH_from_shift_2pKa(δ_obs, pKa1, pKa2, δ_0, δ_1, δ_2; σ_δ=0.01, pH_range=(2.0, 12.0), tol=1e-6)
    δ_obs_val = Measurements.value(δ_obs)
    pKa1_val = Measurements.value(pKa1)
    pKa2_val = Measurements.value(pKa2)
    δ_0_val = Measurements.value(δ_0)
    δ_1_val = Measurements.value(δ_1)
    δ_2_val = Measurements.value(δ_2)

    # Check if δ_obs is within the limiting shifts
    δ_limits = [δ_0_val, δ_1_val, δ_2_val]
    δ_min, δ_max = extrema(δ_limits)
    if δ_obs_val <= δ_min || δ_obs_val >= δ_max
        return missing
    end

    # Define objective function
    function residual(pH)
        pred = henderson_hasselbalch_2pKa([pH], [pKa1_val, pKa2_val, δ_0_val, δ_1_val, δ_2_val])[1]
        return pred - δ_obs_val
    end

    # Simple bisection root-finding
    try
        a, b = pH_range
        fa, fb = residual(a), residual(b)

        # Check that root is bracketed
        if fa * fb > 0
            return missing
        end

        # Bisection iteration
        for _ in 1:100
            mid = (a + b) / 2
            fmid = residual(mid)

            if abs(fmid) < tol || (b - a) / 2 < tol
                pH_sol = mid

                # Estimate uncertainty via numerical gradient
                ε = 0.001
                dδ_dpH = (henderson_hasselbalch_2pKa([pH_sol + ε], [pKa1_val, pKa2_val, δ_0_val, δ_1_val, δ_2_val])[1] -
                          henderson_hasselbalch_2pKa([pH_sol - ε], [pKa1_val, pKa2_val, δ_0_val, δ_1_val, δ_2_val])[1]) / (2ε)

                σ_δ_val = δ_obs isa Measurement ? Measurements.uncertainty(δ_obs) : σ_δ
                σ_pH = abs(σ_δ_val / dδ_dpH)

                return pH_sol ± σ_pH
            end

            if fa * fmid < 0
                b, fb = mid, fmid
            else
                a, fa = mid, fmid
            end
        end

        return missing  # Failed to converge
    catch
        return missing
    end
end

"""
Check if a pH estimate is reliable (uncertainty below threshold).
"""
function is_reliable(pH_estimate; threshold=0.15)
    ismissing(pH_estimate) && return false
    return Measurements.uncertainty(pH_estimate) < threshold
end
