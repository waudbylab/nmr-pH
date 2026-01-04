# pH Indicator Calibration from NMR Chemical Shifts

## Project Overview

This project analyses NMR chemical shift data for pH indicator molecules to build models relating chemical shift to pH, temperature, and ionic strength. The goal is to enable inference of sample conditions (pH, T, I) from measured chemical shifts (implemented in the parent folder!).

## Repository Structure

```
src/
  pH_indicator_analysis.jl    # Core analysis functions
scripts/
  go-group1.jl                # Analysis script for group 1 data
data/
  group1.csv                  # Experimental data
```

## Dependencies

```julia
using Pkg
Pkg.add(["CSV", "DataFrames", "LsqFit", "Measurements", "Plots", "Printf", "Statistics"])
```

## Data Format

CSV files with columns:
- `expt_number`: experiment identifier (for reference only)
- `ionic_strength`: ionic strength in M (calculated from charge balance)
- `temperature`: temperature in K
- `nominal_pH`: pH measured by electrode at 298 K
- `P_xxx`, `F_xxx`, `H_xxx`: chemical shifts (ppm, DSS-referenced) where prefix indicates nucleus (³¹P, ¹⁹F, ¹H)

---

## Current Implementation

### The Henderson-Hasselbalch Model

For an indicator in fast exchange between protonated (HA) and deprotonated (A) forms:

```
δ_obs = δ_HA + (δ_A - δ_HA) / (1 + 10^(pKa - pH))
```

Inverting to obtain pH from chemical shift:

```
pH = pKa + log₁₀[(δ_obs - δ_HA) / (δ_A - δ_obs)]
```

Uncertainty propagation:

```
σ_pH = σ_δ · (δ_A - δ_HA) / [ln(10) · (δ_obs - δ_HA) · (δ_A - δ_obs)]
```

The uncertainty diverges as δ_obs approaches either limiting shift. An indicator provides reliable pH only within ~1.5 units of its pKa.

### Bootstrap Algorithm (Current)

The current implementation performs a **single-pass bootstrap**:

1. **Initial fits**: Fit Henderson-Hasselbalch to all indicators using electrode pH → obtain pKa, δ_HA, δ_A for each
2. **Anchor to phosphate**: Fix phosphate pKa to literature value, calculate pH_true from ³¹P shifts
3. **Refit other indicators**: Using pH_true (from phosphate only), refit pKa for other indicators with limiting shifts held fixed
4. **Combine pH estimates**: Inverse-variance weighted combination of pH from all calibrated indicators

### Key Functions

```julia
# Load and reshape data
data = load_indicator_data("data/group1.csv")

# Select one (T, I) condition
subset = select_condition(data, temperature=298.0, ionic_strength=0.15)

# Initial fits using electrode pH
initial_fits = fit_all_indicators(subset)
plot_initial_fits(subset, initial_fits)

# Bootstrap using literature phosphate pKa
bootstrap_result = bootstrap_pH(subset, initial_fits; pKa_phosphate_lit=7.20)
plot_bootstrap_results(subset, bootstrap_result)

# Or run everything in one call
results = run_analysis("data/group1.csv"; 
    temperature=298.0, ionic_strength=0.025, pKa_phosphate_lit=7.20)
```

### Outputs

- Formatted tables of fitted parameters with uncertainties
- Comparison of electrode pH vs bootstrap pH
- Multi-panel diagnostic plots

---

## Next Steps for Implementation

### 1. Iterated Bootstrap Algorithm

**Problem**: The current implementation only uses phosphate to establish pH_true. Indicators with pKa outside the phosphate-sensitive range (roughly pH 5.5–8.5) are calibrated using few reliable points.

**Solution**: Iterate outwards in pH from the fixed phosphate pKa based on pKa.

**Proposed algorithm**:

Take steps of X pH units, e.g. 0.5 units, away from the phosphate reference.

```
Tier 0 (reference):
  - Phosphate (pKa FIXED to literature)

Tier 1 (initial pKa estimate within X units of tier 0):
  - Calibrate using phosphate-derived pH
  - Update pH_combined using tiers 0+1

Tier 2 (initial pKa estimate within X units of tier 1):
  - Calibrate using pH_combined from tiers 0+1
  - Update pH_combined using tiers 0+1+2

Tier 3: iterate as needed
```

Actually, a better approach may simply be to select the next closest estimated pKa, re-fit that using calibrated pH values, update pH values, then iterate. Do this!

**Implementation tasks**:
- [ ] Define indicator ordering based on ΔpKa from phosphate reference
- [ ] Process indicators sequentially, updating pH_combined after each
- [ ] Track calibration status per indicator


### 2. Indicators with Multiple pKa Values

**Problem**: Some indicators (maleate, phosphate itself for pKa1/pKa3, piperazine, etc...) have multiple ionisations within the measured pH range. The single-pKa Henderson-Hasselbalch model doesn't fit these.

**Solution**: Implement extended models for multi-pKa indicators.

**Two-pKa model**:
```
D = 1 + 10^(pH - pKa1) + 10^(2·pH - pKa1 - pKa2)
δ_obs = (δ_0 + δ_1·10^(pH-pKa1) + δ_2·10^(2·pH-pKa1-pKa2)) / D
```

Where δ_0, δ_1, δ_2 are the limiting shifts of the fully protonated, singly deprotonated, and doubly deprotonated forms.

**Implementation tasks**:
- [ ] Add `henderson_hasselbalch_2pKa(pH, p)` model function
- [ ] Add indicator metadata specifying number of pKa values
- [ ] Modify `fit_indicator` to select appropriate model
- [ ] Handle pH inversion for 2-pKa case (may need numerical solution)
- [ ] Extend `IndicatorFit` struct to accommodate multiple pKa values

**Data structure**:
```julia
struct IndicatorProperties
    name::String
    n_pKa::Int                    # 1 or 2
    z_acid::Int                   # charge of most protonated form
    pKa_literature::Vector{Float64}  # literature values if known (starting guesses for fitting)
end
```

Some input will need to be worked out for this, as well as mapping buffers with multiple resonances or multiple nuclei. Something simple like data column names = X_name_Y, where X=nucleus and Y=resonance number perhaps.

### 3. Cross-Condition Analysis

**Problem**: After analysing each (T, I) condition independently, we have pKa and limiting shift values at multiple conditions (e.g. 3 or 4 temperatures × 3 ionic strengths). These need to be combined into unified models.

**Solution**: Fit temperature and ionic strength dependence to build predictive models.

#### 3.1 Temperature Dependence of pKa

**Van 't Hoff equation**:
```
pKa(T) = pKa(T_ref) + (ΔH° / (R·ln10)) · (1/T - 1/T_ref)
```

Where R = 8.314 J/(mol·K), T_ref = 298.15 K.

**Implementation tasks**:
- [ ] Collect pKa values across temperatures at fixed I
- [ ] Fit van 't Hoff to extract ΔH° for each indicator
- [ ] Compare ΔH° to literature values where available
- [ ] Expected ranges: carboxylic acids 0–5 kJ/mol, amines 20–50 kJ/mol, phosphate 4–8 kJ/mol

#### 3.2 Ionic Strength Dependence of pKa

**Davies equation**:
```
pKa(I) = pKa° - A · Δz² · (√I / (1 + √I) - 0.3·I)
```

Where A ≈ 0.51 at 298 K (Debye-Hückel constant), Δz² = z_base² - z_acid².

**Implementation tasks**:
- [ ] Collect pKa values across ionic strengths at fixed T
- [ ] Fit Davies equation to extract pKa° and effective A
- [ ] Verify Δz² matches expected from indicator chemistry
- [ ] Consider whether A should be fitted globally or per-indicator

#### 3.3 Temperature Dependence of Limiting Shifts

Chemical shifts have weak but measurable temperature dependence:
```
δ_HA(T) = δ_HA° + β_HA · (T - T_ref)
δ_A(T) = δ_A° + β_A · (T - T_ref)
```

**Implementation tasks**:
- [ ] Fit linear T-dependence of limiting shifts
- [ ] Verify limiting shifts are independent of I (they should be)
- [ ] Typical β values: ~0.001–0.01 ppm/K for ¹H

#### 3.4 Combined Model

The final model for each indicator:
```
pKa(T, I) = pKa° + (ΔH° / R·ln10) · (1/T - 1/T_ref) - A·Δz² · f(I)
δ_HA(T) = δ_HA° + β_HA · (T - T_ref)
δ_A(T) = δ_A° + β_A · (T - T_ref)
δ(pH, T, I) = δ_HA(T) + [δ_A(T) - δ_HA(T)] / [1 + 10^(pKa(T,I) - pH)]
```

**Implementation tasks**:
- [ ] Define `IndicatorModel` struct with all parameters
- [ ] Implement `predict_shift(model, pH, T, I)`
- [ ] Implement `predict_pH(model, δ_obs, T, I)` with uncertainty
- [ ] Export models to JSON for downstream use

**Proposed struct**:
```julia
struct IndicatorModel
    name::String
    nucleus::String
    
    # Charge information
    z_acid::Int
    z_base::Int
    Δz²::Int
    
    # pKa parameters
    pKa°::Measurement{Float64}      # at T_ref, I=0
    ΔH°::Measurement{Float64}       # J/mol
    A_eff::Measurement{Float64}     # Debye-Hückel parameter
    
    # Limiting shift parameters
    δ_HA°::Measurement{Float64}     # at T_ref
    β_HA::Measurement{Float64}      # ppm/K
    δ_A°::Measurement{Float64}      # at T_ref  
    β_A::Measurement{Float64}       # ppm/K
end
```

### 4. Validation and Quality Control

**Implementation tasks**:
- [ ] Compare all fitted pKa° to literature compilation
- [ ] Check internal consistency (multi-signal indicators should give same pKa)
- [ ] Residual analysis: no systematic trends with pH, T, or I
- [ ] Cross-validation: hold out one condition, predict from others
- [ ] Uncertainty coverage: do 68% of measurements fall within 1σ of predictions?

---

## Physical Constants

```julia
const R_gas = 8.314          # J/(mol·K)
const ln10 = log(10)         # 2.303
const T_ref = 298.15         # K
const A_debye_298 = 0.510    # Debye-Hückel constant at 298 K
```

## Literature pKa Values (298 K, I → 0)

| Indicator | pKa | ΔH° (kJ/mol) | Δz² | Source |
|-----------|-----|--------------|-----|--------|
| Phosphate (pKa2) | 7.198 | 3.6 | 3 | NIST |
| Acetate | 4.756 | 0.4 | 1 | NIST |
| Formate | 3.75 | 0 | 1 | NIST |
| Imidazole | 6.99 | 36.6 | -1 | NIST |
| Tris | 8.07 | 47.4 | -1 | NIST |
| HEPES | 7.55 | 20.4 | 1 | Good et al. |
| PIPES | 6.76 | 11.2 | 1 | Good et al. |

---

## Notes

- All chemical shifts are referenced to DSS
- Ionic strength is calculated from charge balance (see ionic_strength_calculator.xlsx)
- The electrode pH was measured at 298 K; this is NOT the true pH at other temperatures
- Phosphate pKa values at different (T, I) should be calculated from literature parameters, not fitted
