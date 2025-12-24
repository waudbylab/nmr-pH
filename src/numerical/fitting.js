/**
 * Fitting Module
 *
 * Nonlinear least-squares fitting using Nelder-Mead simplex optimization.
 * Estimates pH, temperature, ionic strength, and reference offsets from chemical shifts.
 * Uses proper χ² formulation with combined observed and predicted uncertainties.
 */

import { nelderMead } from './nelderMead.js';
import { predictShift, getBufferPKaValues, predictShiftWithUncertainty } from './bufferModel.js';
import { assignPeaks, getAssignedPeaksForFitting } from './peakAssignment.js';

/**
 * Default observed shift uncertainties per nucleus (ppm).
 * These represent typical experimental precision.
 */
export const DEFAULT_SHIFT_UNCERTAINTIES = {
  '1H': 0.01,
  '13C': 0.1,
  '15N': 0.1,
  '19F': 0.05,
  '31P': 0.1
};

/**
 * Updated IUPAC Xi ratios from BMRB/Iowa State.
 * These are fractional frequencies relative to 1H = 1.0
 */
const XI_RATIOS = {
  '1H': 1.0,
  '19F': 0.94094008,
  '31P': 0.404808636,
  '13C': 0.251449530,
  '15N': 0.101329118
};

/**
 * Calculate temperature-dependent water reference shift for 1H.
 * When DSS is not present, we use water as a secondary reference.
 *
 * @param {number} temperature - Temperature in Kelvin
 * @returns {number} Initial guess for 1H reference offset (ppm)
 */
export function calculateWaterReference(temperature) {
  return temperature / 96.9 - 3.13;
}

/**
 * Calculate reference offset for nucleus X from 1H reference offset.
 * Uses spectrometer frequency relationship via Xi ratios.
 */
function calculateLinkedOffset(nucleus, protonFrequencyMHz, protonOffsetPpm) {
  if (nucleus === '1H') return protonOffsetPpm;

  const xiX = XI_RATIOS[nucleus];
  const xiH = XI_RATIOS['1H'];

  if (!xiX) return 0;

  const nuDSS_H = protonFrequencyMHz * (1 - protonOffsetPpm / 1e6);
  const nu0_X = nuDSS_H * (xiX / xiH);
  const bfX = protonFrequencyMHz * (xiX / xiH);
  const deltaRef_X = ((bfX - nu0_X) / bfX) * 1e6;

  return deltaRef_X;
}

/**
 * Default fitting options.
 */
const DEFAULT_OPTIONS = {
  refineTemperature: false,
  refineIonicStrength: false,
  refineReferences: {},
  referenceBounds: {},
  linkedToProton: [],
  protonFrequency: null,
  maxIterations: 500,
  tolerance: 1e-8,
  initialPH: 7.0,
  useGridSearch: true,
  shiftUncertainties: DEFAULT_SHIFT_UNCERTAINTIES
};

/**
 * Parameter bounds for optimization.
 */
const PARAM_BOUNDS = {
  pH: { min: 0, max: 14 },
  temperature: { min: 273, max: 373 },
  ionicStrength: { min: 0, max: 1 },
  refDefault1H: { min: -1, max: 1 },
  refDefaultHetero: { min: -10, max: 10 }
};

/**
 * Build parameter vector from conditions.
 */
export function buildParameterVector(conditions, options) {
  const params = [conditions.pH];
  const parameterMap = {
    pH: { index: 0, name: 'pH' }
  };

  let index = 1;

  if (options.refineTemperature) {
    params.push(conditions.temperature);
    parameterMap.temperature = { index, name: 'Temperature (K)' };
    index++;
  }

  if (options.refineIonicStrength) {
    params.push(conditions.ionicStrength);
    parameterMap.ionicStrength = { index, name: 'Ionic strength (M)' };
    index++;
  }

  for (const [nucleus, refine] of Object.entries(options.refineReferences)) {
    if (refine) {
      params.push(conditions.referenceOffsets?.[nucleus] ?? 0);
      parameterMap[`ref_${nucleus}`] = { index, name: `${nucleus} reference offset (ppm)` };
      index++;
    }
  }

  return { params, parameterMap };
}

/**
 * Extract conditions from parameter vector.
 */
export function extractConditions(params, parameterMap, baseConditions, options = {}) {
  const conditions = {
    pH: params[parameterMap.pH.index],
    temperature: parameterMap.temperature
      ? params[parameterMap.temperature.index]
      : baseConditions.temperature,
    ionicStrength: parameterMap.ionicStrength
      ? params[parameterMap.ionicStrength.index]
      : baseConditions.ionicStrength,
    referenceOffsets: { ...baseConditions.referenceOffsets }
  };

  for (const [key, mapping] of Object.entries(parameterMap)) {
    if (key.startsWith('ref_')) {
      const nucleus = key.slice(4);
      conditions.referenceOffsets[nucleus] = params[mapping.index];
    }
  }

  if (options.linkedToProton && options.linkedToProton.length > 0 && options.protonFrequency) {
    const protonOffset = conditions.referenceOffsets['1H'] ?? 0;
    for (const nucleus of options.linkedToProton) {
      conditions.referenceOffsets[nucleus] = calculateLinkedOffset(
        nucleus,
        options.protonFrequency,
        protonOffset
      );
    }
  }

  return conditions;
}

/**
 * Get bounds for a parameter.
 */
function getParameterBounds(paramKey, options, initialValue) {
  if (paramKey === 'pH') return PARAM_BOUNDS.pH;
  if (paramKey === 'temperature') return PARAM_BOUNDS.temperature;
  if (paramKey === 'ionicStrength') return PARAM_BOUNDS.ionicStrength;

  if (paramKey.startsWith('ref_')) {
    const nucleus = paramKey.slice(4);
    const customBounds = options.referenceBounds?.[nucleus];
    if (customBounds) return customBounds;
    return nucleus === '1H' ? PARAM_BOUNDS.refDefault1H : PARAM_BOUNDS.refDefaultHetero;
  }

  return { min: -100, max: 100 };
}

/**
 * Transform parameters from unbounded to bounded space.
 * Uses logit transform: x_bounded = min + (max-min) * sigmoid(x_unbounded)
 */
function transformToBounded(unboundedParams, parameterMap, options, initialParams) {
  const result = [];
  for (let idx = 0; idx < unboundedParams.length; idx++) {
    const val = unboundedParams[idx];
    const key = Object.keys(parameterMap).find(k => parameterMap[k].index === idx);
    const bounds = getParameterBounds(key, options, initialParams[idx]);
    const sigmoid = 1 / (1 + Math.exp(-val));
    result.push(bounds.min + (bounds.max - bounds.min) * sigmoid);
  }
  return result;
}

/**
 * Transform parameters from bounded to unbounded space.
 * Uses inverse logit (logit): x_unbounded = log((x - min) / (max - x))
 */
function transformToUnbounded(boundedParams, parameterMap, options, initialParams) {
  const result = [];
  for (let idx = 0; idx < boundedParams.length; idx++) {
    const val = boundedParams[idx];
    const key = Object.keys(parameterMap).find(k => parameterMap[k].index === idx);
    const bounds = getParameterBounds(key, options, initialParams[idx]);
    // Clamp to avoid infinities
    const clamped = Math.max(bounds.min + 1e-10, Math.min(bounds.max - 1e-10, val));
    const normalized = (clamped - bounds.min) / (bounds.max - bounds.min);
    result.push(Math.log(normalized / (1 - normalized)));
  }
  return result;
}

/**
 * Wrapper to call predictShiftWithUncertainty from bufferModel.
 * Passes the buffer object needed for pKa uncertainties.
 */
function getPredictedShiftWithUncertainty(resonance, buffer, pKaValues, pH, temperature, ionicStrength, refTemp, refIonic) {
  return predictShiftWithUncertainty(
    resonance, buffer, pKaValues, pH, temperature, ionicStrength, refTemp, refIonic
  );
}

/**
 * Calculate χ² for given parameters.
 * χ² = Σ[(δ_obs - δ_pred)² / (σ_obs² + σ_pred²)]
 */
function calculateChiSquared(
  params,
  parameterMap,
  assignedPeaks,
  buffersMap,
  samplesMap,
  baseConditions,
  options
) {
  const conditions = extractConditions(params, parameterMap, baseConditions, options);
  let chiSq = 0;

  for (const peak of assignedPeaks) {
    const buffer = buffersMap.get(peak.buffer_id);
    const sample = samplesMap.get(buffer.sample_id);
    const resonances = buffer.chemical_shifts[peak.nucleus] ?? [];
    const resonance = resonances.find(r => r.resonance_id === peak.resonance_id);

    if (!resonance) continue;

    const refTemp = sample?.reference_temperature_K ?? 298.15;
    const refIonic = sample?.reference_ionic_strength_M ?? 0;
    const pKaValues = getBufferPKaValues(buffer, conditions.temperature, conditions.ionicStrength, refTemp);

    const { shift: predictedShift, uncertainty: sigmaPred } = getPredictedShiftWithUncertainty(
      resonance,
      buffer,
      pKaValues,
      conditions.pH,
      conditions.temperature,
      conditions.ionicStrength,
      refTemp,
      refIonic
    );

    // Apply reference offset
    const refOffset = conditions.referenceOffsets[peak.nucleus] ?? 0;
    const totalPredicted = predictedShift + refOffset;

    // Get observed uncertainty
    const sigmaObs = options.shiftUncertainties?.[peak.nucleus] ?? DEFAULT_SHIFT_UNCERTAINTIES[peak.nucleus] ?? 0.01;

    // Combined variance
    const sigmaTotal = Math.sqrt(sigmaObs * sigmaObs + sigmaPred * sigmaPred);

    // Residual as z-score
    const residual = peak.observed_shift - totalPredicted;
    const zScore = residual / sigmaTotal;

    chiSq += zScore * zScore;
  }

  return chiSq;
}

/**
 * Perform grid search to find initial parameter estimates.
 */
export function gridSearchInitialParameters(observedShifts, buffers, samplesMap, baseConditions, options) {
  const gridOptions = {
    ...options,
    refineTemperature: false,
    refineIonicStrength: false
  };

  const { parameterMap } = buildParameterVector({ ...baseConditions, pH: 7 }, gridOptions);

  const refNuclei = Object.entries(options.refineReferences)
    .filter(([_, refine]) => refine)
    .map(([nucleus]) => nucleus);

  // Collect all resonances for SSR calculation
  const buffersMap = new Map(buffers.map(b => [b.buffer_id, b]));

  const pHStep = 0.2;
  const pHMin = 0;
  const pHMax = 14;

  let bestParams = null;
  let bestChiSq = Infinity;

  // Simple SSR-based grid search (not full χ² yet, for speed)
  const calculateSSR = (pH, refOffsets) => {
    const testConditions = { ...baseConditions, pH, referenceOffsets: { ...baseConditions.referenceOffsets, ...refOffsets } };
    let ssr = 0;

    for (const [nucleus, shifts] of Object.entries(observedShifts)) {
      if (!shifts || shifts.length === 0) continue;

      const refOffset = testConditions.referenceOffsets[nucleus] ?? 0;
      const predictedShifts = [];

      for (const buffer of buffers) {
        const resonances = buffer.chemical_shifts[nucleus];
        if (!resonances) continue;

        const sample = samplesMap.get(buffer.sample_id);
        const refTemp = sample?.reference_temperature_K ?? 298.15;
        const refIonic = sample?.reference_ionic_strength_M ?? 0;
        const pKaValues = getBufferPKaValues(buffer, testConditions.temperature, testConditions.ionicStrength, refTemp);

        for (const resonance of resonances) {
          const predicted = predictShift(resonance, pKaValues, pH, testConditions.temperature, testConditions.ionicStrength, refTemp, refIonic);
          predictedShifts.push(predicted + refOffset);
        }
      }

      for (const observed of shifts) {
        if (predictedShifts.length === 0) continue;
        let minDist = Infinity;
        for (const predicted of predictedShifts) {
          const dist = Math.abs(observed - predicted);
          if (dist < minDist) minDist = dist;
        }
        ssr += minDist * minDist;
      }
    }

    return ssr;
  };

  if (refNuclei.length === 0) {
    for (let pH = pHMin; pH <= pHMax; pH += pHStep) {
      const ssr = calculateSSR(pH, {});
      if (ssr < bestChiSq) {
        bestChiSq = ssr;
        bestParams = [pH];
      }
    }
  } else if (refNuclei.length === 1 && refNuclei[0] === '1H') {
    const h1Bounds = options.referenceBounds?.['1H'] || { min: -0.5, max: 0.5 };
    const h1Step = 0.05;

    for (let pH = pHMin; pH <= pHMax; pH += pHStep) {
      for (let h1Ref = h1Bounds.min; h1Ref <= h1Bounds.max; h1Ref += h1Step) {
        const ssr = calculateSSR(pH, { '1H': h1Ref });
        if (ssr < bestChiSq) {
          bestChiSq = ssr;
          bestParams = [pH, h1Ref];
        }
      }
    }
  } else {
    // Multi-parameter grid search with coarser steps
    const nParams = 1 + refNuclei.length;
    const paramRanges = [{ min: pHMin, max: pHMax, step: pHStep }];

    for (const nucleus of refNuclei) {
      const bounds = options.referenceBounds?.[nucleus] || (nucleus === '1H' ? { min: -0.5, max: 0.5 } : { min: -5, max: 5 });
      const step = nucleus === '1H' ? 0.1 : 1.0;
      paramRanges.push({ min: bounds.min, max: bounds.max, step });
    }

    function searchGrid(paramIndex, currentParams) {
      if (paramIndex === nParams) {
        const refOffsets = {};
        for (let nucIdx = 0; nucIdx < refNuclei.length; nucIdx++) {
          refOffsets[refNuclei[nucIdx]] = currentParams[nucIdx + 1];
        }
        const ssr = calculateSSR(currentParams[0], refOffsets);
        if (ssr < bestChiSq) {
          bestChiSq = ssr;
          bestParams = [...currentParams];
        }
        return;
      }

      const range = paramRanges[paramIndex];
      for (let value = range.min; value <= range.max; value += range.step) {
        currentParams[paramIndex] = value;
        searchGrid(paramIndex + 1, currentParams);
      }
    }

    searchGrid(0, new Array(nParams));
  }

  return { params: bestParams, chiSq: bestChiSq };
}

/**
 * Compute numerical Hessian of χ² at the minimum.
 * Uses central finite differences.
 *
 * @returns {Array<Array<number>>} Hessian matrix
 */
function computeHessian(chiSqFn, params, delta = 1e-5) {
  const n = params.length;
  const H = Array(n).fill(null).map(() => Array(n).fill(0));

  for (let row = 0; row < n; row++) {
    for (let col = row; col < n; col++) {
      const p_pp = [...params];
      const p_pm = [...params];
      const p_mp = [...params];
      const p_mm = [...params];

      p_pp[row] += delta; p_pp[col] += delta;
      p_pm[row] += delta; p_pm[col] -= delta;
      p_mp[row] -= delta; p_mp[col] += delta;
      p_mm[row] -= delta; p_mm[col] -= delta;

      const f_pp = chiSqFn(p_pp);
      const f_pm = chiSqFn(p_pm);
      const f_mp = chiSqFn(p_mp);
      const f_mm = chiSqFn(p_mm);

      H[row][col] = (f_pp - f_pm - f_mp + f_mm) / (4 * delta * delta);
      H[col][row] = H[row][col]; // Symmetric
    }
  }

  return H;
}

/**
 * Invert a matrix using Gaussian elimination with partial pivoting.
 */
function invertMatrix(matrix) {
  const n = matrix.length;
  // Create augmented matrix [A | I]
  const aug = [];
  for (let rowIdx = 0; rowIdx < n; rowIdx++) {
    const newRow = [...matrix[rowIdx]];
    for (let colIdx = 0; colIdx < n; colIdx++) {
      newRow.push(rowIdx === colIdx ? 1 : 0);
    }
    aug.push(newRow);
  }

  for (let col = 0; col < n; col++) {
    // Find pivot
    let maxRow = col;
    for (let row = col + 1; row < n; row++) {
      if (Math.abs(aug[row][col]) > Math.abs(aug[maxRow][col])) {
        maxRow = row;
      }
    }
    [aug[col], aug[maxRow]] = [aug[maxRow], aug[col]];

    if (Math.abs(aug[col][col]) < 1e-12) {
      // Matrix is singular or nearly singular
      return null;
    }

    // Scale pivot row
    const scale = aug[col][col];
    for (let j = 0; j < 2 * n; j++) {
      aug[col][j] /= scale;
    }

    // Eliminate column
    for (let row = 0; row < n; row++) {
      if (row !== col) {
        const factor = aug[row][col];
        for (let j = 0; j < 2 * n; j++) {
          aug[row][j] -= factor * aug[col][j];
        }
      }
    }
  }

  return aug.map(row => row.slice(n));
}

/**
 * Calculate parameter uncertainties from Hessian.
 * Covariance matrix C = 2 * H^(-1) (factor of 2 from χ² definition)
 * Standard errors are sqrt of diagonal elements.
 */
function calculateUncertaintiesFromHessian(H) {
  const Hinv = invertMatrix(H);
  if (!Hinv) {
    return null; // Matrix was singular
  }

  // Covariance = 2 * H^(-1) for χ² minimization
  // Standard errors = sqrt(diag(Cov))
  const uncertainties = [];
  for (let idx = 0; idx < Hinv.length; idx++) {
    const variance = 2 * Hinv[idx][idx];
    uncertainties.push(variance > 0 ? Math.sqrt(variance) : null);
  }
  return uncertainties;
}

/**
 * Create residual data for each peak (for reporting).
 */
function calculateResiduals(params, parameterMap, assignedPeaks, buffersMap, samplesMap, baseConditions, options) {
  const conditions = extractConditions(params, parameterMap, baseConditions, options);
  const residuals = [];

  for (const peak of assignedPeaks) {
    const buffer = buffersMap.get(peak.buffer_id);
    const sample = samplesMap.get(buffer.sample_id);
    const resonances = buffer.chemical_shifts[peak.nucleus] ?? [];
    const resonance = resonances.find(r => r.resonance_id === peak.resonance_id);

    if (!resonance) continue;

    const refTemp = sample?.reference_temperature_K ?? 298.15;
    const refIonic = sample?.reference_ionic_strength_M ?? 0;
    const pKaValues = getBufferPKaValues(buffer, conditions.temperature, conditions.ionicStrength, refTemp);

    const { shift: predictedShift, uncertainty: sigmaPred } = getPredictedShiftWithUncertainty(
      resonance, buffer, pKaValues, conditions.pH, conditions.temperature, conditions.ionicStrength, refTemp, refIonic
    );

    const refOffset = conditions.referenceOffsets[peak.nucleus] ?? 0;
    const totalPredicted = predictedShift + refOffset;

    const sigmaObs = options.shiftUncertainties?.[peak.nucleus] ?? DEFAULT_SHIFT_UNCERTAINTIES[peak.nucleus] ?? 0.01;
    const sigmaTotal = Math.sqrt(sigmaObs * sigmaObs + sigmaPred * sigmaPred);

    const residual = peak.observed_shift - totalPredicted;

    residuals.push({
      nucleus: peak.nucleus,
      observed: peak.observed_shift,
      predicted: totalPredicted,
      residual,
      sigmaObs,
      sigmaPred,
      sigmaTotal,
      zScore: residual / sigmaTotal
    });
  }

  return residuals;
}

/**
 * Perform nonlinear least-squares fitting using Nelder-Mead.
 */
export function fitParameters(observedShifts, buffers, samplesMap, initialConditions, options = {}) {
  const opts = { ...DEFAULT_OPTIONS, ...options };

  // Initial assignment
  let assignments = assignPeaks(
    observedShifts,
    buffers,
    samplesMap,
    opts.initialPH ?? initialConditions.pH ?? 7.0,
    initialConditions.temperature,
    initialConditions.ionicStrength
  );

  let assignedPeaks = getAssignedPeaksForFitting(assignments);

  if (assignedPeaks.length === 0) {
    return {
      success: false,
      error: 'No peaks could be assigned to buffer resonances',
      assignments
    };
  }

  const buffersMap = new Map(buffers.map(b => [b.buffer_id, b]));

  const baseConditions = {
    temperature: initialConditions.temperature,
    ionicStrength: initialConditions.ionicStrength,
    referenceOffsets: initialConditions.referenceOffsets ?? {},
    pH: opts.initialPH ?? initialConditions.pH ?? 7.0
  };

  // Grid search for initial parameters
  let initialParams;
  let parameterMap;

  if (opts.useGridSearch) {
    const gridResult = gridSearchInitialParameters(observedShifts, buffers, samplesMap, baseConditions, opts);
    if (gridResult.params) {
      const gridConditions = { ...baseConditions, pH: gridResult.params[0] };

      const refNuclei = Object.entries(opts.refineReferences)
        .filter(([_, refine]) => refine)
        .map(([nucleus]) => nucleus);

      for (let nucIdx = 0; nucIdx < refNuclei.length; nucIdx++) {
        gridConditions.referenceOffsets[refNuclei[nucIdx]] = gridResult.params[nucIdx + 1];
      }

      const result = buildParameterVector(gridConditions, opts);
      initialParams = result.params;
      parameterMap = result.parameterMap;

      // Re-assign with grid search pH
      assignments = assignPeaks(observedShifts, buffers, samplesMap, gridConditions.pH, initialConditions.temperature, initialConditions.ionicStrength);
      assignedPeaks = getAssignedPeaksForFitting(assignments);
    } else {
      const result = buildParameterVector(baseConditions, opts);
      initialParams = result.params;
      parameterMap = result.parameterMap;
    }
  } else {
    const result = buildParameterVector(baseConditions, opts);
    initialParams = result.params;
    parameterMap = result.parameterMap;
  }

  const nParams = initialParams.length;
  const nObs = assignedPeaks.length;
  const dof = nObs - nParams;

  if (dof < 0) {
    return {
      success: false,
      error: `Underdetermined system: ${nObs} observations, ${nParams} parameters (DoF = ${dof})`,
      assignments,
      nObservations: nObs,
      nParameters: nParams,
      degreesOfFreedom: dof
    };
  }

  // Create χ² function for optimization
  const chiSqFn = (params) => calculateChiSquared(
    params, parameterMap, assignedPeaks, buffersMap, samplesMap, baseConditions, opts
  );

  // Transform to unbounded space for Nelder-Mead
  const unboundedInitial = transformToUnbounded(initialParams, parameterMap, opts, initialParams);

  // Objective function in unbounded space
  const objective = (unboundedParams) => {
    const boundedParams = transformToBounded(unboundedParams, parameterMap, opts, initialParams);
    return chiSqFn(boundedParams);
  };

  try {
    // Run Nelder-Mead optimization
    const result = nelderMead(objective, unboundedInitial);

    // Transform back to bounded space
    const fittedParams = transformToBounded(result.x, parameterMap, opts, initialParams);
    const finalChiSq = chiSqFn(fittedParams);

    // Compute Hessian for uncertainties
    const H = computeHessian(chiSqFn, fittedParams);
    const uncertainties = calculateUncertaintiesFromHessian(H);

    const fittedConditions = extractConditions(fittedParams, parameterMap, baseConditions, opts);

    // Re-assign with fitted conditions
    const finalAssignments = assignPeaks(
      observedShifts, buffers, samplesMap,
      fittedConditions.pH, fittedConditions.temperature, fittedConditions.ionicStrength
    );

    // Calculate detailed residuals
    const residualsData = calculateResiduals(fittedParams, parameterMap, assignedPeaks, buffersMap, samplesMap, baseConditions, opts);
    const residuals = residualsData.map(r => r.residual);
    const sumSquares = residuals.reduce((sum, r) => sum + r * r, 0);
    const rmsd = Math.sqrt(sumSquares / residuals.length);
    const reducedChiSq = dof > 0 ? finalChiSq / dof : finalChiSq;

    // Build parameter results
    const parameterResults = {};
    for (const [key, mapping] of Object.entries(parameterMap)) {
      parameterResults[key] = {
        value: fittedParams[mapping.index],
        uncertainty: uncertainties ? uncertainties[mapping.index] : null,
        name: mapping.name
      };
    }

    // Add linked reference offsets
    if (opts.linkedToProton && opts.linkedToProton.length > 0) {
      const protonOffset = fittedConditions.referenceOffsets['1H'] ?? 0;
      for (const nucleus of opts.linkedToProton) {
        const linkedOffset = calculateLinkedOffset(nucleus, opts.protonFrequency, protonOffset);
        parameterResults[`ref_${nucleus}_linked`] = {
          value: linkedOffset,
          uncertainty: null,
          name: `${nucleus} reference offset (linked to ¹H)`,
          isLinked: true
        };
      }
    }

    return {
      success: true,
      parameters: parameterResults,
      conditions: fittedConditions,
      assignments: finalAssignments,
      residuals,
      residualsDetailed: residualsData,
      statistics: {
        nObservations: nObs,
        nParameters: nParams,
        degreesOfFreedom: dof,
        sumSquares,
        rmsd,
        chiSquared: finalChiSq,
        reducedChiSquared: reducedChiSq
      },
      convergence: {
        converged: true,
        finalValue: result.fx
      }
    };
  } catch (error) {
    return {
      success: false,
      error: `Fitting failed: ${error.message}`,
      assignments
    };
  }
}

/**
 * Perform iterative fitting with reassignment.
 */
export function fitWithReassignment(
  observedShifts,
  buffers,
  samplesMap,
  initialConditions,
  options = {},
  maxRounds = 3
) {
  let conditions = { ...initialConditions };
  let result = null;

  for (let round = 0; round < maxRounds; round++) {
    result = fitParameters(observedShifts, buffers, samplesMap, conditions, {
      ...options,
      initialPH: conditions.pH ?? options.initialPH ?? 7.0
    });

    if (!result.success) {
      return result;
    }

    const pHChange = Math.abs(result.conditions.pH - conditions.pH);

    if (pHChange < 0.1) {
      break;
    }

    conditions = result.conditions;
  }

  return result;
}
