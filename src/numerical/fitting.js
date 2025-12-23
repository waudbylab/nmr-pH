/**
 * Fitting Module
 *
 * Nonlinear least-squares fitting using Levenberg-Marquardt algorithm.
 * Estimates pH, temperature, ionic strength, and reference offsets from chemical shifts.
 * Includes grid search for initial parameter estimation.
 */

import { levenbergMarquardt } from 'ml-levenberg-marquardt';
import { predictBufferShifts, getBufferPKaValues, predictShift } from './bufferModel.js';
import { assignPeaks, getAssignedPeaksForFitting } from './peakAssignment.js';
import { calculateParameterUncertainties } from './uncertainty.js';

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
  // Temperature-dependent water shift formula: T/96.9 - 3.13
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

  // Calculate DSS frequency for 1H
  const nuDSS_H = protonFrequencyMHz * (1 - protonOffsetPpm / 1e6);

  // Calculate true zero frequency for X
  const nu0_X = nuDSS_H * (xiX / xiH);

  // Calculate bf(X)
  const bfX = protonFrequencyMHz * (xiX / xiH);

  // Calculate X reference offset
  const deltaRef_X = ((bfX - nu0_X) / bfX) * 1e6;

  return deltaRef_X;
}

/**
 * Default fitting options.
 */
const DEFAULT_OPTIONS = {
  refineTemperature: false,
  refineIonicStrength: false,
  refineReferences: {}, // { nucleus: boolean }
  referenceBounds: {}, // { nucleus: { min, max } }
  linkedToProton: [], // nuclei whose reference is linked to 1H via spectrometer frequency
  protonFrequency: null, // spectrometer 1H frequency in MHz
  maxIterations: 100,
  tolerance: 1e-8,
  initialPH: 7.0,
  useGridSearch: true // Enable grid search for initial parameters
};

/**
 * Build parameter vector from conditions.
 *
 * @param {Object} conditions - Current conditions
 * @param {Object} options - Fitting options
 * @returns {Object} { params: Array, parameterMap: Object }
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

  // Reference offsets per nucleus
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
 *
 * @param {Array<number>} params - Parameter vector
 * @param {Object} parameterMap - Map of parameter names to indices
 * @param {Object} baseConditions - Base conditions for fixed parameters
 * @param {Object} options - Fitting options (for linked references)
 * @returns {Object} Conditions object
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

  // Extract fitted reference offsets
  for (const [key, mapping] of Object.entries(parameterMap)) {
    if (key.startsWith('ref_')) {
      const nucleus = key.slice(4);
      conditions.referenceOffsets[nucleus] = params[mapping.index];
    }
  }

  // Calculate linked reference offsets (scenario 3)
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
 * Calculate sum of squared residuals for given parameters.
 * Used for grid search optimization.
 *
 * @param {Array<number>} params - Parameter vector
 * @param {Object} parameterMap - Map of parameter names to indices
 * @param {Object} observedShifts - Observed chemical shifts per nucleus
 * @param {Array<Object>} buffers - Selected buffers
 * @param {Map<string, Object>} samplesMap - Sample objects
 * @param {Object} baseConditions - Base conditions
 * @param {Object} options - Fitting options
 * @returns {number} Sum of squared residuals
 */
function calculateSSR(params, parameterMap, observedShifts, buffers, samplesMap, baseConditions, options) {
  const conditions = extractConditions(params, parameterMap, baseConditions, options);
  let ssr = 0;

  // For each nucleus with observed shifts
  for (const [nucleus, shifts] of Object.entries(observedShifts)) {
    if (!shifts || shifts.length === 0) continue;

    const refOffset = conditions.referenceOffsets[nucleus] ?? 0;

    // Collect all predicted shifts for this nucleus from all buffers
    const predictedShifts = [];
    for (const buffer of buffers) {
      const resonances = buffer.chemical_shifts[nucleus];
      if (!resonances) continue;

      const sample = samplesMap.get(buffer.sample_id);
      const refTemp = sample?.reference_temperature_K ?? 298.15;
      const refIonic = sample?.reference_ionic_strength_M ?? 0;
      const pKaValues = getBufferPKaValues(buffer, conditions.temperature, conditions.ionicStrength, refTemp);

      for (const resonance of resonances) {
        const predicted = predictShift(
          resonance,
          pKaValues,
          conditions.pH,
          conditions.temperature,
          conditions.ionicStrength,
          refTemp,
          refIonic
        );
        // Apply reference offset to predicted shift
        predictedShifts.push(predicted + refOffset);
      }
    }

    // For each observed shift, find closest predicted and accumulate residual
    for (const observed of shifts) {
      if (predictedShifts.length === 0) continue;

      // Find closest predicted shift
      let minDist = Infinity;
      for (const predicted of predictedShifts) {
        const dist = Math.abs(observed - predicted);
        if (dist < minDist) {
          minDist = dist;
        }
      }
      ssr += minDist * minDist;
    }
  }

  return ssr;
}

/**
 * Perform grid search to find initial parameter estimates.
 * Fixes temperature and ionic strength at nominal values.
 *
 * @param {Object} observedShifts - Observed chemical shifts per nucleus
 * @param {Array<Object>} buffers - Selected buffers
 * @param {Map<string, Object>} samplesMap - Sample objects
 * @param {Object} baseConditions - Base conditions
 * @param {Object} options - Fitting options
 * @returns {Object} Best parameters found { params, ssr }
 */
export function gridSearchInitialParameters(observedShifts, buffers, samplesMap, baseConditions, options) {
  // Build parameter map for the grid search (pH + reference offsets only)
  // Temperature and ionic strength are fixed during grid search
  const gridOptions = {
    ...options,
    refineTemperature: false,
    refineIonicStrength: false
  };

  const { parameterMap } = buildParameterVector({ ...baseConditions, pH: 7 }, gridOptions);

  // Determine which reference offsets need fitting
  const refNuclei = Object.entries(options.refineReferences)
    .filter(([_, refine]) => refine)
    .map(([nucleus]) => nucleus);

  // Grid parameters
  const pHStep = 0.2;
  const pHMin = 0;
  const pHMax = 14;

  let bestParams = null;
  let bestSSR = Infinity;

  // Determine grid strategy based on what needs fitting
  if (refNuclei.length === 0) {
    // Only pH to search - simple 1D grid
    for (let pH = pHMin; pH <= pHMax; pH += pHStep) {
      const params = [pH];
      const ssr = calculateSSR(params, parameterMap, observedShifts, buffers, samplesMap, baseConditions, gridOptions);
      if (ssr < bestSSR) {
        bestSSR = ssr;
        bestParams = [...params];
      }
    }
  } else if (refNuclei.length === 1 && refNuclei[0] === '1H') {
    // pH + 1H reference - 2D grid
    const h1Bounds = options.referenceBounds?.['1H'] || { min: -0.5, max: 0.5 };
    const h1Step = 0.05;

    for (let pH = pHMin; pH <= pHMax; pH += pHStep) {
      for (let h1Ref = h1Bounds.min; h1Ref <= h1Bounds.max; h1Ref += h1Step) {
        const params = [pH, h1Ref];
        const testConditions = { ...baseConditions, referenceOffsets: { ...baseConditions.referenceOffsets, '1H': h1Ref } };
        const ssr = calculateSSR(params, parameterMap, observedShifts, buffers, samplesMap, testConditions, gridOptions);
        if (ssr < bestSSR) {
          bestSSR = ssr;
          bestParams = [...params];
        }
      }
    }
  } else {
    // Multiple reference offsets - do coarse grid for each
    // This is more complex; we'll do a nested grid but with coarser steps

    // Build initial params array structure
    const nParams = 1 + refNuclei.length; // pH + ref offsets
    const paramRanges = [{ min: pHMin, max: pHMax, step: pHStep }]; // pH

    for (const nucleus of refNuclei) {
      const bounds = options.referenceBounds?.[nucleus] || (nucleus === '1H' ? { min: -0.5, max: 0.5 } : { min: -5, max: 5 });
      const step = nucleus === '1H' ? 0.1 : 0.5; // Coarser grid for multi-parameter search
      paramRanges.push({ min: bounds.min, max: bounds.max, step });
    }

    // Recursive grid search
    function searchGrid(paramIndex, currentParams) {
      if (paramIndex === nParams) {
        // Evaluate at this point
        const testConditions = { ...baseConditions };
        testConditions.referenceOffsets = { ...baseConditions.referenceOffsets };
        for (let i = 0; i < refNuclei.length; i++) {
          testConditions.referenceOffsets[refNuclei[i]] = currentParams[i + 1];
        }
        const ssr = calculateSSR(currentParams, parameterMap, observedShifts, buffers, samplesMap, testConditions, gridOptions);
        if (ssr < bestSSR) {
          bestSSR = ssr;
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

  return { params: bestParams, ssr: bestSSR };
}

/**
 * Create residual function for fitting.
 *
 * @param {Array<Object>} assignedPeaks - Assigned peaks for fitting
 * @param {Map<string, Object>} buffersMap - Map of buffer_id to buffer object
 * @param {Map<string, Object>} samplesMap - Map of sample_id to sample object
 * @param {Object} parameterMap - Map of parameter names to indices
 * @param {Object} baseConditions - Base conditions for fixed parameters
 * @param {Object} options - Fitting options
 * @returns {Function} Residual function for optimizer
 */
export function createResidualFunction(assignedPeaks, buffersMap, samplesMap, parameterMap, baseConditions, options) {
  return function(params) {
    const conditions = extractConditions(params, parameterMap, baseConditions, options);
    const residuals = [];

    for (const peak of assignedPeaks) {
      const buffer = buffersMap.get(peak.buffer_id);
      const sample = samplesMap.get(buffer.sample_id);

      if (!buffer) {
        throw new Error(`Buffer not found: ${peak.buffer_id}`);
      }

      // Find the resonance
      const resonances = buffer.chemical_shifts[peak.nucleus] ?? [];
      const resonance = resonances.find(r => r.resonance_id === peak.resonance_id);

      if (!resonance) {
        throw new Error(`Resonance not found: ${peak.resonance_id} in ${peak.buffer_id}`);
      }

      // Calculate predicted shift
      const refTemp = sample?.reference_temperature_K ?? 298.15;
      const refIonic = sample?.reference_ionic_strength_M ?? 0;
      const pKaValues = getBufferPKaValues(buffer, conditions.temperature, conditions.ionicStrength, refTemp);

      let predictedShift = predictShift(
        resonance,
        pKaValues,
        conditions.pH,
        conditions.temperature,
        conditions.ionicStrength,
        refTemp,
        refIonic
      );

      // Apply reference offset if applicable
      const refOffset = conditions.referenceOffsets[peak.nucleus] ?? 0;
      predictedShift += refOffset;

      // Residual = observed - predicted
      residuals.push(peak.observed_shift - predictedShift);
    }

    return residuals;
  };
}

/**
 * Create model function for Levenberg-Marquardt.
 * Returns predicted shifts given parameters.
 *
 * @param {Array<Object>} assignedPeaks - Assigned peaks for fitting
 * @param {Map<string, Object>} buffersMap - Map of buffer_id to buffer object
 * @param {Map<string, Object>} samplesMap - Map of sample_id to sample object
 * @param {Object} parameterMap - Map of parameter names to indices
 * @param {Object} baseConditions - Base conditions for fixed parameters
 * @param {Object} options - Fitting options
 * @returns {Function} Model function
 */
export function createModelFunction(assignedPeaks, buffersMap, samplesMap, parameterMap, baseConditions, options) {
  // ml-levenberg-marquardt expects a function that takes (params) and returns array of predicted values
  // It minimizes the sum of squares of (data - predicted)
  return function(params) {
    const conditions = extractConditions(params, parameterMap, baseConditions, options);
    const predicted = [];

    for (const peak of assignedPeaks) {
      const buffer = buffersMap.get(peak.buffer_id);
      const sample = samplesMap.get(buffer.sample_id);

      const resonances = buffer.chemical_shifts[peak.nucleus] ?? [];
      const resonance = resonances.find(r => r.resonance_id === peak.resonance_id);

      const refTemp = sample?.reference_temperature_K ?? 298.15;
      const refIonic = sample?.reference_ionic_strength_M ?? 0;
      const pKaValues = getBufferPKaValues(buffer, conditions.temperature, conditions.ionicStrength, refTemp);

      let predictedShift = predictShift(
        resonance,
        pKaValues,
        conditions.pH,
        conditions.temperature,
        conditions.ionicStrength,
        refTemp,
        refIonic
      );

      const refOffset = conditions.referenceOffsets[peak.nucleus] ?? 0;
      predictedShift += refOffset;

      predicted.push(predictedShift);
    }

    return predicted;
  };
}

/**
 * Build parameter bounds based on options.
 *
 * @param {Array<number>} initialParams - Initial parameter values
 * @param {Object} parameterMap - Parameter mapping
 * @param {Object} options - Fitting options with referenceBounds
 * @returns {Object} { minValues, maxValues }
 */
function buildParameterBounds(initialParams, parameterMap, options) {
  const minValues = [];
  const maxValues = [];

  for (let i = 0; i < initialParams.length; i++) {
    // pH bounds
    if (i === parameterMap.pH.index) {
      minValues.push(0);
      maxValues.push(14);
      continue;
    }

    // Temperature bounds
    if (parameterMap.temperature && i === parameterMap.temperature.index) {
      minValues.push(273);
      maxValues.push(373);
      continue;
    }

    // Ionic strength bounds
    if (parameterMap.ionicStrength && i === parameterMap.ionicStrength.index) {
      minValues.push(0);
      maxValues.push(1);
      continue;
    }

    // Reference offset bounds - check for custom bounds
    let foundBounds = false;
    for (const [key, mapping] of Object.entries(parameterMap)) {
      if (key.startsWith('ref_') && mapping.index === i) {
        const nucleus = key.slice(4);
        const bounds = options.referenceBounds?.[nucleus];
        if (bounds) {
          minValues.push(bounds.min);
          maxValues.push(bounds.max);
        } else if (nucleus === '1H') {
          // Default 1H bounds: initial guess ± 0.5
          const guess = initialParams[i];
          minValues.push(guess - 0.5);
          maxValues.push(guess + 0.5);
        } else {
          // Default heteronuclear bounds: ±5 ppm
          minValues.push(-5);
          maxValues.push(5);
        }
        foundBounds = true;
        break;
      }
    }

    if (!foundBounds) {
      // Fallback
      minValues.push(-10);
      maxValues.push(10);
    }
  }

  return { minValues, maxValues };
}

/**
 * Perform nonlinear least-squares fitting to estimate parameters.
 *
 * @param {Object} observedShifts - Object mapping nucleus -> array of observed shifts
 * @param {Array<Object>} buffers - Array of selected buffer objects
 * @param {Map<string, Object>} samplesMap - Map of sample_id to sample object
 * @param {Object} initialConditions - Initial conditions (pH, temperature, ionicStrength)
 * @param {Object} [options] - Fitting options
 * @returns {Object} Fitting results
 */
export function fitParameters(observedShifts, buffers, samplesMap, initialConditions, options = {}) {
  const opts = { ...DEFAULT_OPTIONS, ...options };

  // Initial assignment at initial conditions
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

  // Build maps for quick lookup
  const buffersMap = new Map(buffers.map(b => [b.buffer_id, b]));

  // Build base conditions
  const baseConditions = {
    temperature: initialConditions.temperature,
    ionicStrength: initialConditions.ionicStrength,
    referenceOffsets: initialConditions.referenceOffsets ?? {},
    pH: opts.initialPH ?? initialConditions.pH ?? 7.0
  };

  // Optionally run grid search for initial parameters
  let initialParams;
  let parameterMap;

  if (opts.useGridSearch) {
    const gridResult = gridSearchInitialParameters(observedShifts, buffers, samplesMap, baseConditions, opts);
    if (gridResult.params) {
      // Update base conditions with grid search results
      const gridConditions = { ...baseConditions, pH: gridResult.params[0] };

      // Extract reference offsets from grid search
      const refNuclei = Object.entries(opts.refineReferences)
        .filter(([_, refine]) => refine)
        .map(([nucleus]) => nucleus);

      for (let i = 0; i < refNuclei.length; i++) {
        gridConditions.referenceOffsets[refNuclei[i]] = gridResult.params[i + 1];
      }

      // Build parameter vector from grid search result
      const result = buildParameterVector(gridConditions, opts);
      initialParams = result.params;
      parameterMap = result.parameterMap;

      // Re-assign peaks with grid search pH for better starting assignments
      assignments = assignPeaks(
        observedShifts,
        buffers,
        samplesMap,
        gridConditions.pH,
        initialConditions.temperature,
        initialConditions.ionicStrength
      );
      assignedPeaks = getAssignedPeaksForFitting(assignments);
    } else {
      // Grid search failed, use default initial params
      const result = buildParameterVector(baseConditions, opts);
      initialParams = result.params;
      parameterMap = result.parameterMap;
    }
  } else {
    const result = buildParameterVector(baseConditions, opts);
    initialParams = result.params;
    parameterMap = result.parameterMap;
  }

  // Check degrees of freedom
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

  // Create model function
  const modelFn = createModelFunction(assignedPeaks, buffersMap, samplesMap, parameterMap, baseConditions, opts);

  // Build parameter bounds using custom bounds
  const { minValues, maxValues } = buildParameterBounds(initialParams, parameterMap, opts);

  let fittedParams;
  let iterations = 0;

  // Levenberg-Marquardt requires at least 2 data points
  // For 1 observation, use 1D root-finding to find exact pH intercept
  if (nObs < 2) {
    // For a single observation, find the exact pH where predicted shift equals observed shift
    // Use bisection method for robust 1D root-finding
    fittedParams = [...initialParams];

    if (nObs === 1) {
      // Get the single assigned peak
      const peak = assignedPeaks[0];
      const buffer = buffersMap.get(peak.buffer_id);
      const sample = samplesMap.get(buffer.sample_id);
      const resonances = buffer.chemical_shifts[peak.nucleus] ?? [];
      const resonance = resonances.find(r => r.resonance_id === peak.resonance_id);

      if (resonance) {
        const refTemp = sample?.reference_temperature_K ?? 298.15;
        const refIonic = sample?.reference_ionic_strength_M ?? 0;

        // Function to compute residual (observed - predicted) at given pH
        const computeResidual = (testPH) => {
          // Build conditions with test pH
          const testParams = [...fittedParams];
          testParams[0] = testPH;
          const testConditions = extractConditions(testParams, parameterMap, baseConditions, opts);

          const pKaValues = getBufferPKaValues(buffer, testConditions.temperature, testConditions.ionicStrength, refTemp);
          let predictedShift = predictShift(
            resonance,
            pKaValues,
            testPH,
            testConditions.temperature,
            testConditions.ionicStrength,
            refTemp,
            refIonic
          );

          // Apply reference offset
          const refOffset = testConditions.referenceOffsets[peak.nucleus] ?? 0;
          predictedShift += refOffset;

          return peak.observed_shift - predictedShift;
        };

        // Bisection method to find pH where residual is zero
        let pHLow = 0;
        let pHHigh = 14;
        const tolerance = 1e-6;
        const maxIter = 50;

        let residualLow = computeResidual(pHLow);
        let residualHigh = computeResidual(pHHigh);

        // Check if solution exists in the range
        if (residualLow * residualHigh > 0) {
          // No sign change - solution may not exist or is at an extreme
          // Fall back to finding minimum residual via golden section search
          const phi = (1 + Math.sqrt(5)) / 2;
          let a = pHLow, b = pHHigh;
          let c = b - (b - a) / phi;
          let d = a + (b - a) / phi;

          for (let i = 0; i < maxIter; i++) {
            const fc = Math.abs(computeResidual(c));
            const fd = Math.abs(computeResidual(d));

            if (fc < fd) {
              b = d;
              d = c;
              c = b - (b - a) / phi;
            } else {
              a = c;
              c = d;
              d = a + (b - a) / phi;
            }

            if (Math.abs(b - a) < tolerance) break;
          }

          fittedParams[0] = (a + b) / 2;
        } else {
          // Bisection when there's a sign change
          for (let i = 0; i < maxIter; i++) {
            const pHMid = (pHLow + pHHigh) / 2;
            const residualMid = computeResidual(pHMid);

            if (Math.abs(residualMid) < tolerance || (pHHigh - pHLow) / 2 < tolerance) {
              fittedParams[0] = pHMid;
              break;
            }

            if (residualMid * residualLow < 0) {
              pHHigh = pHMid;
              residualHigh = residualMid;
            } else {
              pHLow = pHMid;
              residualLow = residualMid;
            }

            fittedParams[0] = pHMid;
          }
        }

        iterations = maxIter;
      }
    }
  } else {
    // Prepare data for ml-levenberg-marquardt
    // x values are just indices (we need them but they're not really used)
    const xData = assignedPeaks.map((_, i) => i);
    const yData = assignedPeaks.map(p => p.observed_shift);

    // Wrap model function to match expected signature
    const wrappedModel = params => (x => modelFn(params)[x]);

    try {
      // Run Levenberg-Marquardt
      const result = levenbergMarquardt(
        { x: xData, y: yData },
        wrappedModel,
        {
          initialValues: initialParams,
          minValues,
          maxValues,
          maxIterations: opts.maxIterations,
          errorTolerance: opts.tolerance,
          damping: 1.5,
          dampingStepUp: 11,
          dampingStepDown: 9
        }
      );

      fittedParams = result.parameterValues;
      iterations = result.iterations;
    } catch (error) {
      return {
        success: false,
        error: `Fitting failed: ${error.message}`,
        assignments
      };
    }
  }

  try {
    const fittedConditions = extractConditions(fittedParams, parameterMap, baseConditions, opts);

    // Re-assign peaks with fitted conditions
    const finalAssignments = assignPeaks(
      observedShifts,
      buffers,
      samplesMap,
      fittedConditions.pH,
      fittedConditions.temperature,
      fittedConditions.ionicStrength
    );

    // Calculate residuals at final parameters
    const residualFn = createResidualFunction(assignedPeaks, buffersMap, samplesMap, parameterMap, baseConditions, opts);
    const residuals = residualFn(fittedParams);
    const sumSquares = residuals.reduce((sum, r) => sum + r * r, 0);
    const rmsd = Math.sqrt(sumSquares / residuals.length);

    // Calculate chi-squared (assuming unit variance for now)
    const chiSquared = sumSquares;
    const reducedChiSquared = dof > 0 ? chiSquared / dof : chiSquared;

    // Calculate parameter uncertainties
    const uncertainties = calculateParameterUncertainties(
      fittedParams,
      residualFn,
      reducedChiSquared
    );

    // Build result object
    const parameterResults = {};
    for (const [key, mapping] of Object.entries(parameterMap)) {
      parameterResults[key] = {
        value: fittedParams[mapping.index],
        uncertainty: uncertainties[mapping.index],
        name: mapping.name
      };
    }

    // Add linked reference offsets to results (calculated, not fitted)
    if (opts.linkedToProton && opts.linkedToProton.length > 0) {
      const protonOffset = fittedConditions.referenceOffsets['1H'] ?? 0;
      for (const nucleus of opts.linkedToProton) {
        const linkedOffset = calculateLinkedOffset(nucleus, opts.protonFrequency, protonOffset);
        parameterResults[`ref_${nucleus}_linked`] = {
          value: linkedOffset,
          uncertainty: null, // Derived value, uncertainty would need error propagation
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
      statistics: {
        nObservations: nObs,
        nParameters: nParams,
        degreesOfFreedom: dof,
        sumSquares,
        rmsd,
        chiSquared,
        reducedChiSquared,
        iterations
      },
      convergence: {
        converged: iterations < opts.maxIterations,
        iterations
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
 * Fits, reassigns peaks with new conditions, and refits until stable.
 *
 * @param {Object} observedShifts - Object mapping nucleus -> array of observed shifts
 * @param {Array<Object>} buffers - Array of selected buffer objects
 * @param {Map<string, Object>} samplesMap - Map of sample_id to sample object
 * @param {Object} initialConditions - Initial conditions
 * @param {Object} [options] - Fitting options
 * @param {number} [maxRounds=3] - Maximum number of fit-reassign rounds
 * @returns {Object} Final fitting results
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

    // Check if pH changed significantly (would affect assignments)
    const pHChange = Math.abs(result.conditions.pH - conditions.pH);

    if (pHChange < 0.1) {
      // Assignments unlikely to change, stop iterating
      break;
    }

    // Update conditions for next round
    conditions = result.conditions;
  }

  return result;
}
