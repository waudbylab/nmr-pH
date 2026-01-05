# Fix for pH Correction Issues at High pH

## Changes Made

### 1. Relaxed Range Checking in Chemical Shift Calculations

**File:** `julia/src/models.jl`

**Functions Modified:**
- `pH_from_shift()` (line 89)
- `pH_from_shift_2pKa()` (line 119)

**Change:** Allow 10% extrapolation beyond limiting chemical shifts

**Before:**
```julia
if δ_obs_val <= δ_min || δ_obs_val >= δ_max
    return missing
end
```

**After:**
```julia
# Allow 10% extrapolation tolerance
δ_range = abs(δ_A_val - δ_HA_val)
δ_min = min(δ_HA_val, δ_A_val) - 0.10 * δ_range
δ_max = max(δ_HA_val, δ_A_val) + 0.10 * δ_range

if δ_obs_val <= δ_min || δ_obs_val >= δ_max
    return missing
end
```

**Rationale:** At extreme pH, observed chemical shifts can slightly exceed limiting shifts due to measurement noise, temperature/ionic strength effects, or additional ionization equilibria. The strict range check was excluding valid data points, especially for Tris at pH >9.

### 2. Added Distance-from-pKa Weighting

**File:** `julia/src/bootstrap.jl`

**Function Modified:** `update_combined_pH!()` (line 349)

**Change:** Combined inverse-variance weighting with quality weighting based on distance from pKa

**New Weighting Scheme:**
```julia
# Inverse variance weight (precision)
weights_variance = 1 ./ uncertainties .^ 2

# Quality weight (distance from pKa)
quality = exp(-(dist / 1.5)^2)

# Combined weight
weights = weights_variance .* weights_quality
```

**Rationale:** Indicators are most reliable within ±1.5 pH units of their pKa. At greater distances, they approach limiting shifts where small measurement errors cause large pH uncertainties. The Gaussian weighting scheme (σ = 1.5 pH units) gives:
- Full weight at pKa
- 60% weight at pKa ± 1
- 10% weight at pKa ± 2.5
- Near-zero weight beyond ±3

### 3. Added Warnings for Extreme pH

**File:** `julia/src/bootstrap.jl`

**Function Modified:** `print_pH_comparison()` (line 445)

**Change:** Added warning messages for:
- Low indicator count (n < 4)
- High uncertainty (σ > 0.10)
- Large pH corrections (|Δ| > 0.30)

**Example Output:**
```
WARNINGS:
  ⚠ pH 9.42: Only 3 indicators contributing (low confidence)
  ⚠ pH 9.75: Large correction (Δ = -0.140)
  ⚠ pH 10.20: Only 2 indicators contributing (low confidence)
  ⚠ pH 10.20: Large correction (Δ = -0.535)

Note: Results at extreme pH may be less reliable due to:
  - Fewer indicators in their sensitive range
  - Extrapolation beyond measured chemical shift limits
  - Consider adding buffers with pKa values at extremes
```

## Expected Impact

### Before Fix:

| pH | n_indicators | Δ(pH) | Issue |
|----|--------------|-------|-------|
| 8.0 | 9 | +0.197 | ✓ Good |
| 9.0 | 8 | +0.191 | ✓ Good |
| 9.42 | 4 | +0.009 | ⚠ Unreliable |
| 9.75 | 3 | -0.140 | ❌ Poor |
| 10.20 | 2 | -0.535 | ❌ Very poor |

### After Fix (Expected):

1. **More indicators contribute at high pH** - Tris and other indicators near limiting shifts can now contribute with appropriate weighting

2. **Smoother pH corrections** - Distance weighting prevents saturated indicators from dominating

3. **Better uncertainty estimates** - Quality weighting properly inflates uncertainties when extrapolating

4. **Clear warnings** - Users are alerted to potential issues at extreme pH

## Testing Recommendations

After deploying these changes:

1. **Run the fitting script** and verify:
   - Tris contributes to pH 9-10 estimates
   - n_indicators stays ≥ 4-5 across full range
   - pH corrections are smooth (no jumps)
   - Warnings appear for pH <4 and pH >9.5

2. **Validate results** against:
   - Literature pKa values
   - Independent pH measurements
   - Cross-validation with different buffer subsets

3. **Check edge cases**:
   - Very high pH (>10)
   - Very low pH (<4)
   - Single-indicator scenarios
   - Missing data patterns

## Future Improvements

1. **Add high-pKa buffers** - Include buffers with pKa 9-11 (CAPS, CHES, carbonate) for better coverage at high pH

2. **Temperature/ionic strength corrections** - Improve modeling of these effects on limiting shifts

3. **Adaptive reliability threshold** - Make the 0.15 threshold pH-dependent

4. **Bootstrap validation** - Add cross-validation to assess systematic biases

## Files Modified

- `julia/src/models.jl` - Relaxed range checking
- `julia/src/bootstrap.jl` - Added distance weighting and warnings

## Files Added

- `julia/HIGH_PH_ISSUE.md` - Detailed analysis of the problem
- `julia/scripts/diagnose_high_ph.jl` - Diagnostic script (requires Julia)
- `julia/CHANGES.md` - This file
