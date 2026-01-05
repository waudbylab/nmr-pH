# Bootstrap pH Correction Issue at High pH

## Problem Summary

The bootstrap pH correction algorithm shows increasingly poor performance at high pH (>9), with the pH correction shifting from +0.19 at pH 8 to -0.54 at pH 10.2. The number of contributing indicators drops from 9 at pH 7-8 to only 2 at pH 10.2.

## Root Cause Analysis

### 1. **Strict Range Checking in `pH_from_shift`**

The function `julia/src/models.jl:89-110` calculates pH from chemical shift using:

```julia
if δ_obs_val <= δ_min || δ_obs_val >= δ_max
    return missing
end
```

This **strict inequality** excludes any observation outside the limiting chemical shifts `[δ_HA, δ_A]`.

### 2. **Observed Shifts Beyond Limits at Extreme pH**

Looking at the data at pH 10.2:

| Indicator | δ_obs (ppm) | δ_A (ppm) | Status |
|-----------|-------------|-----------|---------|
| Tris | 3.509-3.513 | 3.5052 | **OUTSIDE** (δ_obs < δ_A) |
| phosphate | 2.659-2.769 | 2.5618 | **OUTSIDE** (δ_obs > δ_A) |
| imidazole1 | 7.763-7.772 | 7.7579 | **OUTSIDE** (δ_obs > δ_A) |
| imidazole2 | 7.118-7.134 | 7.1188 | **OUTSIDE** (δ_obs > δ_A) |

At pH 10.2 (far above most pKa values), indicators are nearly fully deprotonated and their observed shifts are very close to or slightly beyond `δ_A` due to:

- **Measurement noise** (~0.005 ppm typical uncertainty)
- **Temperature/ionic strength effects** not fully captured in the model
- **Additional ionization equilibria** at extreme pH
- **Systematic calibration errors**

When δ_obs exceeds the limiting shift by even 0.001 ppm, the indicator is excluded (`returns missing`).

### 3. **Cascade Effect**

As pH increases:

| pH | n_indicators | pH_correction | Notes |
|----|--------------|---------------|-------|
| 7-8 | 9-10 | +0.17 | Most indicators in transition region |
| 8.5-9 | 8-9 | +0.16-0.19 | Some indicators saturating |
| 9.42 | 4 | +0.01 | Tris saturating, others excluded |
| 9.75 | 3 | -0.14 | Mostly phosphate + 2 others |
| 10.2 | 2 | **-0.54** | Only 2 indicators remain |

With fewer indicators:
1. **Reduced averaging** → higher uncertainty
2. **Systematic bias** from remaining indicators
3. **Poor pKa coverage** at extreme pH

### 4. **Tris-Specific Issues**

Tris shows:
- **Higher pKa uncertainty**: ±0.022 vs. ±0.001-0.010 for others
- **Only 23/27 samples** usable (others have 21-27)
- **Calibration tier 8** (later in bootstrap, depends on earlier fits)

At pH 10.2:
- Tris is 1.5 pH units above its pKa (97% deprotonated)
- Observed shift is at the noise floor relative to δ_A
- Small errors push it outside the valid range → excluded
- Without Tris, no indicator covers the pH 9-10 range well

## Why This Matters

The Henderson-Hasselbalch equation for calculating pH from shift:

```
pH = pKa + log₁₀[(δ_obs - δ_HA) / (δ_A - δ_obs)]
```

Near the limits:
- As δ_obs → δ_A, the denominator → 0
- Small measurement errors → large pH errors
- Uncertainty σ_pH → ∞

The `is_reliable(pH_estimate; threshold=0.15)` function correctly identifies these as unreliable, but the algorithm doesn't gracefully handle the case where *all* indicators become unreliable.

## Solutions

### Solution 1: **Relax Range Checking** (Quick Fix)

Modify `julia/src/models.jl:96-99`:

```julia
# Allow small extrapolation (5% beyond limits)
δ_range = abs(δ_A_val - δ_HA_val)
δ_min = min(δ_HA_val, δ_A_val) - 0.05 * δ_range
δ_max = max(δ_HA_val, δ_A_val) + 0.05 * δ_range

if δ_obs_val <= δ_min || δ_obs_val >= δ_max
    return missing
end
```

**Pros:** Simple, allows indicators near limits to contribute
**Cons:** Extrapolation can have large uncertainties

### Solution 2: **Distance-from-pKa Weighting** (Better)

Modify `julia/src/bootstrap.jl:342-367` to weight indicators by their proximity to pKa:

```julia
function update_combined_pH!(pH_data, bootstrap_fits, calibrated)
    for i in 1:nrow(pH_data)
        estimates = Measurement{Float64}[]
        weights_quality = Float64[]

        for ind in calibrated
            pH_est = pH_data[i, Symbol("pH_", ind)]
            if !ismissing(pH_est) && is_reliable(pH_est)
                push!(estimates, pH_est)

                # Weight by distance from pKa (optimal at pKa ± 1)
                pKa = Measurements.value(bootstrap_fits[ind].pKa_bootstrap[1])
                pH_val = Measurements.value(pH_est)
                dist = abs(pH_val - pKa)

                # Gaussian weight: max at pKa, drops off at extremes
                quality_weight = exp(-(dist / 1.5)^2)
                push!(weights_quality, quality_weight)
            end
        end

        if length(estimates) > 0
            values = [Measurements.value(e) for e in estimates]
            uncertainties = [Measurements.uncertainty(e) for e in estimates]

            # Combined weighting: inverse variance × quality
            weights_var = 1 ./ uncertainties .^ 2
            weights = weights_var .* weights_quality

            pH_combined = sum(values .* weights) / sum(weights)
            # ... rest of calculation
        end
    end
end
```

**Pros:** Statistically sound, reduces influence of extrapolated indicators
**Cons:** More complex, requires tuning of distance threshold

### Solution 3: **Bootstrap in pH-Ordered Tiers** (Best)

Modify the calibration algorithm to bootstrap in pH order rather than pKa proximity:

1. Start with reference (phosphate, pKa=6.8)
2. Calibrate indicators with pKa close to pH values with good coverage
3. For high pH, use **only** indicators with pKa > 7.5
4. For low pH, use **only** indicators with pKa < 6

This prevents low-pKa indicators (saturated at high pH) from biasing the estimate.

### Solution 4: **Flag and Report** (Minimum)

At minimum, add warnings for extreme pH:

```julia
if n_ind < 4
    @warn "pH $pH_nom: Only $n_ind indicators contribute. Result may be unreliable."
end

if any(dist > 2 for dist in distances_from_pKa)
    @warn "pH $pH_nom: Indicators are >2 pH units from their pKa. Consider extrapolation."
end
```

## Recommended Fix

**Implement Solutions 1 + 2**:

1. Allow 5-10% extrapolation beyond limits (Solution 1)
2. Add distance-from-pKa weighting (Solution 2)
3. Add warnings for extreme conditions (Solution 4)

This provides a robust solution without major algorithm changes.

## Testing Strategy

After implementing fixes:

1. Verify pH corrections are smooth across full pH range
2. Check that Tris contributes to pH 9-10 estimates
3. Ensure n_indicators doesn't drop below 4-5 at any pH
4. Validate against independent pH measurements at extreme pH
5. Test with synthetic data at pH 2-12 to check edge cases

## Additional Recommendations

Consider expanding the buffer database to include:
- **High-pKa buffers** (pKa 9-11) like CAPS, CHES, or carbonate
- **Low-pKa buffers** (pKa 2-4) for acidic samples
- Multi-nucleus indicators (19F, 31P) which may be more reliable at extremes
