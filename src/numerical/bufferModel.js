/**
 * Buffer Model Module
 *
 * Handles pKa corrections for temperature and ionic strength,
 * ionisation state populations, and chemical shift predictions.
 */

// Physical constants
const R = 8.314462618; // Gas constant (J/mol/K)
const LN10 = Math.log(10);

/**
 * Davies equation coefficient A at 25°C
 * A = 0.5085 for water at 25°C
 */
const DAVIES_A = 0.5085;

/**
 * Extract value from a value that may include uncertainty.
 * Values can be either a number or [value, uncertainty] array.
 *
 * @param {number|Array<number>} valueWithUncertainty - Value or [value, uncertainty]
 * @returns {number} The value
 */
export function getValue(valueWithUncertainty) {
  if (Array.isArray(valueWithUncertainty)) {
    return valueWithUncertainty[0];
  }
  return valueWithUncertainty;
}

/**
 * Extract uncertainty from a value that may include uncertainty.
 *
 * @param {number|Array<number>} valueWithUncertainty - Value or [value, uncertainty]
 * @returns {number} The uncertainty (0 if not provided)
 */
export function getUncertainty(valueWithUncertainty) {
  if (Array.isArray(valueWithUncertainty)) {
    return valueWithUncertainty[1] ?? 0;
  }
  return 0;
}

/**
 * Calculate pKa at a given temperature using van't Hoff equation with heat capacity correction.
 *
 * pKa(T) = pKa(T_ref) + (ΔH/R/ln10)(1/T - 1/T_ref) + (ΔCp/R/ln10)(T_ref/T - 1 + ln(T/T_ref))
 *
 * @param {Object} pKaParams - pKa parameters object from database
 * @param {number} pKaParams.pKa - pKa at reference temperature
 * @param {number} [pKaParams.dH_kJ_mol] - Enthalpy of ionisation (kJ/mol)
 * @param {number} [pKaParams.dCp_J_mol_K] - Heat capacity change (J/mol/K)
 * @param {number} temperature - Temperature (K)
 * @param {number} referenceTemperature - Reference temperature (K), typically 298.15
 * @returns {number} pKa at the specified temperature
 */
export function calculatePKaTemperature(pKaParams, temperature, referenceTemperature = 298.15) {
  const pKaRef = getValue(pKaParams.pKa);
  const dH = getValue(pKaParams.dH_kJ_mol ?? 0) * 1000; // Convert kJ to J
  const dCp = getValue(pKaParams.dCp_J_mol_K ?? 0);

  const T = temperature;
  const Tref = referenceTemperature;

  // van't Hoff with heat capacity correction
  let pKa = pKaRef;

  // Enthalpy term
  if (dH !== 0) {
    pKa += (dH / (R * LN10)) * (1/T - 1/Tref);
  }

  // Heat capacity term
  if (dCp !== 0) {
    pKa += (dCp / (R * LN10)) * (Tref/T - 1 + Math.log(T/Tref));
  }

  return pKa;
}

/**
 * Calculate pKa correction for ionic strength using Davies equation.
 *
 * ΔpKa = A × Δz² × (√I/(1+√I) - 0.3I)
 *
 * where Δz² = z_acid² - z_base² (change in squared charge upon deprotonation)
 *
 * @param {Object} pKaParams - pKa parameters object
 * @param {number} ionicStrength - Ionic strength (M)
 * @returns {number} pKa correction to add
 */
export function calculatePKaIonicStrengthDavies(pKaParams, ionicStrength) {
  if (ionicStrength <= 0) return 0;

  const protonatedCharge = pKaParams.protonated_charge ?? 0;
  const deprotonatedCharge = protonatedCharge - 1;

  // Δz² = z_products² - z_reactants² for: HA ⇌ H⁺ + A⁻
  // Products: H⁺ (charge +1) and A⁻ (deprotonatedCharge)
  // Reactant: HA (protonatedCharge)
  const deltaZSquared = 1 + deprotonatedCharge * deprotonatedCharge - protonatedCharge * protonatedCharge;

  const sqrtI = Math.sqrt(ionicStrength);
  const daviesTerm = sqrtI / (1 + sqrtI) - 0.3 * ionicStrength;

  return DAVIES_A * deltaZSquared * daviesTerm;
}

/**
 * Calculate pKa correction for ionic strength using extended Debye-Hückel equation.
 *
 * @param {Object} pKaParams - pKa parameters object
 * @param {number} ionicStrength - Ionic strength (M)
 * @returns {number} pKa correction to add
 */
export function calculatePKaIonicStrengthExtendedDH(pKaParams, ionicStrength) {
  if (ionicStrength <= 0) return 0;

  const protonatedCharge = pKaParams.protonated_charge ?? 0;
  const deprotonatedCharge = protonatedCharge - 1;
  const ionSize = pKaParams.ion_size_angstrom ?? 4.5; // Default ion size

  const deltaZSquared = 1 + deprotonatedCharge * deprotonatedCharge - protonatedCharge * protonatedCharge;

  const sqrtI = Math.sqrt(ionicStrength);
  // B ≈ 0.328 Å⁻¹ at 25°C for water
  const B = 0.328;
  const dhTerm = sqrtI / (1 + B * ionSize * sqrtI);

  return DAVIES_A * deltaZSquared * dhTerm;
}

/**
 * Calculate pKa correction for ionic strength using empirical linear model.
 *
 * @param {Object} pKaParams - pKa parameters object
 * @param {number} ionicStrength - Ionic strength (M)
 * @returns {number} pKa correction to add
 */
export function calculatePKaIonicStrengthEmpirical(pKaParams, ionicStrength) {
  const coefficient = getValue(pKaParams.ionic_strength_coefficient_per_M ?? 0);
  return coefficient * ionicStrength;
}

/**
 * Calculate pKa at given temperature and ionic strength.
 *
 * @param {Object} pKaParams - pKa parameters object from database
 * @param {number} temperature - Temperature (K)
 * @param {number} ionicStrength - Ionic strength (M)
 * @param {number} referenceTemperature - Reference temperature (K)
 * @returns {number} Corrected pKa value
 */
export function calculatePKa(pKaParams, temperature, ionicStrength, referenceTemperature = 298.15) {
  // Temperature correction
  let pKa = calculatePKaTemperature(pKaParams, temperature, referenceTemperature);

  // Ionic strength correction based on model
  const model = pKaParams.ionic_strength_model ?? 'davies';

  switch (model) {
    case 'davies':
      pKa += calculatePKaIonicStrengthDavies(pKaParams, ionicStrength);
      break;
    case 'extended_debye_huckel':
      pKa += calculatePKaIonicStrengthExtendedDH(pKaParams, ionicStrength);
      break;
    case 'empirical':
      pKa += calculatePKaIonicStrengthEmpirical(pKaParams, ionicStrength);
      break;
    case 'none':
    default:
      // No ionic strength correction
      break;
  }

  return pKa;
}

/**
 * Calculate ionisation state fractions at a given pH.
 *
 * For a molecule with N ionisation states (0 to N-1) and N-1 pKa values:
 * The fraction in each state is determined by the Henderson-Hasselbalch relationships.
 *
 * @param {number} pH - The pH value
 * @param {Array<number>} pKaValues - Array of pKa values (sorted, ascending)
 * @returns {Array<number>} Array of fractions for each ionisation state (0 = most protonated)
 */
export function ionisationFractions(pH, pKaValues) {
  const n = pKaValues.length + 1; // Number of ionisation states

  if (n === 1) {
    return [1.0]; // No ionisation
  }

  // Calculate relative populations using cumulative products
  // For state idx: relative population = product of 10^(pKa_k - pH) for k > idx
  const relatives = new Array(n).fill(1);

  for (let stateIdx = 0; stateIdx < n; stateIdx++) {
    for (let pKaIdx = stateIdx; pKaIdx < pKaValues.length; pKaIdx++) {
      relatives[stateIdx] *= Math.pow(10, pKaValues[pKaIdx] - pH);
    }
  }

  // Normalize to get fractions
  const total = relatives.reduce((sum, r) => sum + r, 0);
  return relatives.map(r => r / total);
}

/**
 * Calculate the limiting chemical shift for a resonance at given conditions.
 * Applies temperature and ionic strength corrections to the base shift.
 *
 * @param {Object} limitingShift - Limiting shift object from database
 * @param {number} temperature - Temperature (K)
 * @param {number} ionicStrength - Ionic strength (M)
 * @param {number} referenceTemperature - Reference temperature (K)
 * @param {number} referenceIonicStrength - Reference ionic strength (M)
 * @returns {number} Corrected chemical shift (ppm)
 */
export function calculateLimitingShift(
  limitingShift,
  temperature,
  ionicStrength,
  referenceTemperature = 298.15,
  referenceIonicStrength = 0
) {
  let shift = getValue(limitingShift.shift_ppm);

  // Temperature correction
  const tempCoeff = getValue(limitingShift.temperature_coefficient_ppm_per_K ?? 0);
  shift += tempCoeff * (temperature - referenceTemperature);

  // Ionic strength correction
  const ionicCoeff = getValue(limitingShift.ionic_strength_coefficient_ppm_per_M ?? 0);
  shift += ionicCoeff * (ionicStrength - referenceIonicStrength);

  return shift;
}

/**
 * Predict the observed chemical shift for a resonance at given conditions.
 * The observed shift is a population-weighted average of limiting shifts.
 *
 * δ_obs = Σ(f_i × δ_i)
 *
 * @param {Object} resonance - Resonance object from database
 * @param {Array<number>} pKaValues - Corrected pKa values at current conditions
 * @param {number} pH - The pH value
 * @param {number} temperature - Temperature (K)
 * @param {number} ionicStrength - Ionic strength (M)
 * @param {number} referenceTemperature - Reference temperature (K)
 * @param {number} referenceIonicStrength - Reference ionic strength (M)
 * @returns {number} Predicted chemical shift (ppm)
 */
export function predictShift(
  resonance,
  pKaValues,
  pH,
  temperature,
  ionicStrength,
  referenceTemperature = 298.15,
  referenceIonicStrength = 0
) {
  // Get ionisation fractions
  const fractions = ionisationFractions(pH, pKaValues);

  // Calculate weighted average shift
  let observedShift = 0;

  for (const limitingShift of resonance.limiting_shifts) {
    const stateIndex = limitingShift.ionisation_state;
    const fraction = fractions[stateIndex] ?? 0;

    const shift = calculateLimitingShift(
      limitingShift,
      temperature,
      ionicStrength,
      referenceTemperature,
      referenceIonicStrength
    );

    observedShift += fraction * shift;
  }

  return observedShift;
}

/**
 * Get all corrected pKa values for a buffer at given conditions.
 *
 * @param {Object} buffer - Buffer object from database
 * @param {number} temperature - Temperature (K)
 * @param {number} ionicStrength - Ionic strength (M)
 * @param {number} referenceTemperature - Reference temperature (K)
 * @returns {Array<number>} Array of corrected pKa values (sorted)
 */
export function getBufferPKaValues(buffer, temperature, ionicStrength, referenceTemperature = 298.15) {
  return buffer.pKa_parameters
    .map(p => calculatePKa(p, temperature, ionicStrength, referenceTemperature))
    .sort((a, b) => a - b);
}

/**
 * Predict all chemical shifts for a buffer at given conditions.
 *
 * @param {Object} buffer - Buffer object from database
 * @param {number} pH - The pH value
 * @param {number} temperature - Temperature (K)
 * @param {number} ionicStrength - Ionic strength (M)
 * @param {Object} sample - Sample object for reference conditions
 * @returns {Object} Object mapping nucleus -> array of {resonance_id, shift}
 */
export function predictBufferShifts(buffer, pH, temperature, ionicStrength, sample) {
  const refTemp = sample?.reference_temperature_K ?? 298.15;
  const refIonic = sample?.reference_ionic_strength_M ?? 0;

  const pKaValues = getBufferPKaValues(buffer, temperature, ionicStrength, refTemp);
  const predictions = {};

  for (const [nucleus, resonances] of Object.entries(buffer.chemical_shifts)) {
    predictions[nucleus] = resonances.map(resonance => ({
      resonance_id: resonance.resonance_id,
      description: resonance.description,
      shift: predictShift(resonance, pKaValues, pH, temperature, ionicStrength, refTemp, refIonic)
    }));
  }

  return predictions;
}

/**
 * Calculate derivative of shift with respect to pH (for uncertainty propagation).
 * Uses numerical differentiation.
 *
 * @param {Object} resonance - Resonance object from database
 * @param {Array<number>} pKaValues - Corrected pKa values
 * @param {number} pH - The pH value
 * @param {number} temperature - Temperature (K)
 * @param {number} ionicStrength - Ionic strength (M)
 * @param {number} referenceTemperature - Reference temperature (K)
 * @param {number} referenceIonicStrength - Reference ionic strength (M)
 * @returns {number} dδ/dpH
 */
function dShiftDpH(resonance, pKaValues, pH, temperature, ionicStrength, refTemp, refIonic) {
  const delta = 0.001;
  const shiftPlus = predictShift(resonance, pKaValues, pH + delta, temperature, ionicStrength, refTemp, refIonic);
  const shiftMinus = predictShift(resonance, pKaValues, pH - delta, temperature, ionicStrength, refTemp, refIonic);
  return (shiftPlus - shiftMinus) / (2 * delta);
}

/**
 * Predict chemical shift with uncertainty propagation from database parameters.
 *
 * Propagates uncertainties from:
 * - pKa values
 * - Limiting shifts
 * - Temperature coefficients
 * - Ionic strength coefficients
 *
 * Uses linear error propagation: σ² = Σ(∂f/∂x_i)² × σ²_xi
 *
 * @param {Object} resonance - Resonance object from database
 * @param {Object} buffer - Buffer object (for pKa uncertainties)
 * @param {Array<number>} pKaValues - Corrected pKa values at current conditions
 * @param {number} pH - The pH value
 * @param {number} temperature - Temperature (K)
 * @param {number} ionicStrength - Ionic strength (M)
 * @param {number} referenceTemperature - Reference temperature (K)
 * @param {number} referenceIonicStrength - Reference ionic strength (M)
 * @returns {Object} { shift, uncertainty }
 */
export function predictShiftWithUncertainty(
  resonance,
  buffer,
  pKaValues,
  pH,
  temperature,
  ionicStrength,
  referenceTemperature = 298.15,
  referenceIonicStrength = 0
) {
  const shift = predictShift(
    resonance, pKaValues, pH, temperature, ionicStrength,
    referenceTemperature, referenceIonicStrength
  );

  const fractions = ionisationFractions(pH, pKaValues);
  let varianceTotal = 0;

  // 1. Uncertainty from limiting shifts
  for (const limitingShift of resonance.limiting_shifts) {
    const stateIndex = limitingShift.ionisation_state;
    const fraction = fractions[stateIndex] ?? 0;

    // Uncertainty in the shift value itself
    const shiftUncertainty = getUncertainty(limitingShift.shift_ppm);
    varianceTotal += (fraction * shiftUncertainty) ** 2;

    // Uncertainty from temperature coefficient (if temperature differs from reference)
    const tempCoeffUncertainty = getUncertainty(limitingShift.temperature_coefficient_ppm_per_K ?? 0);
    const tempDiff = temperature - referenceTemperature;
    varianceTotal += (fraction * tempCoeffUncertainty * tempDiff) ** 2;

    // Uncertainty from ionic strength coefficient
    const ionicCoeffUncertainty = getUncertainty(limitingShift.ionic_strength_coefficient_ppm_per_M ?? 0);
    const ionicDiff = ionicStrength - referenceIonicStrength;
    varianceTotal += (fraction * ionicCoeffUncertainty * ionicDiff) ** 2;
  }

  // 2. Uncertainty from pKa values (affects fractions)
  // Use numerical differentiation to get sensitivity to each pKa
  for (let pKaIdx = 0; pKaIdx < buffer.pKa_parameters.length; pKaIdx++) {
    const pKaUncertainty = getUncertainty(buffer.pKa_parameters[pKaIdx].pKa);
    if (pKaUncertainty > 0) {
      // Approximate dShift/dpKa by shifting the pKa
      const delta = 0.01;
      const pKaValuesPlus = [...pKaValues];
      const pKaValuesMinus = [...pKaValues];
      pKaValuesPlus[pKaIdx] += delta;
      pKaValuesMinus[pKaIdx] -= delta;

      const shiftPlus = predictShift(
        resonance, pKaValuesPlus, pH, temperature, ionicStrength,
        referenceTemperature, referenceIonicStrength
      );
      const shiftMinus = predictShift(
        resonance, pKaValuesMinus, pH, temperature, ionicStrength,
        referenceTemperature, referenceIonicStrength
      );

      const dShiftDpKa = (shiftPlus - shiftMinus) / (2 * delta);
      varianceTotal += (dShiftDpKa * pKaUncertainty) ** 2;
    }
  }

  return {
    shift,
    uncertainty: Math.sqrt(varianceTotal)
  };
}

/**
 * Generate chemical shift vs pH curve data for plotting.
 *
 * @param {Object} buffer - Buffer object from database
 * @param {string} nucleus - Nucleus type (e.g., '1H', '19F')
 * @param {number} temperature - Temperature (K)
 * @param {number} ionicStrength - Ionic strength (M)
 * @param {Object} sample - Sample object for reference conditions
 * @param {number} pHMin - Minimum pH for curve
 * @param {number} pHMax - Maximum pH for curve
 * @param {number} pHStep - pH step size
 * @returns {Array<Object>} Array of curve data objects
 */
export function generateShiftCurves(
  buffer,
  nucleus,
  temperature,
  ionicStrength,
  sample,
  pHMin = 2,
  pHMax = 12,
  pHStep = 0.05
) {
  const refTemp = sample?.reference_temperature_K ?? 298.15;
  const refIonic = sample?.reference_ionic_strength_M ?? 0;

  const pKaValues = getBufferPKaValues(buffer, temperature, ionicStrength, refTemp);
  const resonances = buffer.chemical_shifts[nucleus] ?? [];

  const curves = [];

  for (const resonance of resonances) {
    const pHValues = [];
    const shifts = [];

    for (let pH = pHMin; pH <= pHMax; pH += pHStep) {
      pHValues.push(pH);
      shifts.push(
        predictShift(resonance, pKaValues, pH, temperature, ionicStrength, refTemp, refIonic)
      );
    }

    curves.push({
      buffer_id: buffer.buffer_id,
      buffer_name: buffer.buffer_name,
      resonance_id: resonance.resonance_id,
      description: resonance.description,
      nucleus,
      pHValues,
      shifts
    });
  }

  return curves;
}
