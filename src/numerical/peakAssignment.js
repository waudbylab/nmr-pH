/**
 * Peak Assignment Module
 *
 * Handles automatic matching of observed chemical shifts to buffer resonances.
 */

import { predictBufferShifts } from './bufferModel.js';

/**
 * Default tolerance for considering a peak match (ppm).
 * Varies by nucleus - 1H is more precise than 19F or 31P.
 */
const DEFAULT_TOLERANCES = {
  '1H': 0.5,
  '13C': 2.0,
  '15N': 2.0,
  '19F': 3.0,
  '31P': 2.0
};

/**
 * Calculate assignment confidence based on distance and uniqueness.
 *
 * @param {number} distance - Distance between observed and predicted shift (ppm)
 * @param {number} tolerance - Tolerance for this nucleus (ppm)
 * @param {number} nextBestDistance - Distance to next best match (ppm), Infinity if none
 * @returns {string} Confidence level: 'high', 'medium', 'low', or 'none'
 */
export function calculateConfidence(distance, tolerance, nextBestDistance = Infinity) {
  const absDistance = Math.abs(distance);

  // No match if beyond tolerance
  if (absDistance > tolerance) {
    return 'none';
  }

  // High confidence: close match and well-separated from alternatives
  if (absDistance < tolerance * 0.3 && nextBestDistance > tolerance * 0.6) {
    return 'high';
  }

  // Medium confidence: reasonable match but may have nearby alternatives
  if (absDistance < tolerance * 0.6) {
    return 'medium';
  }

  // Low confidence: within tolerance but not great
  return 'low';
}

/**
 * Generate all possible predictions for selected buffers at given conditions.
 *
 * @param {Array<Object>} buffers - Array of buffer objects
 * @param {Map<string, Object>} samplesMap - Map of sample_id to sample object
 * @param {number} pH - Initial pH estimate
 * @param {number} temperature - Temperature (K)
 * @param {number} ionicStrength - Ionic strength (M)
 * @param {Object} referenceOffsets - Object mapping nucleus -> reference offset (ppm) to subtract
 * @returns {Object} Object mapping nucleus -> array of predictions
 */
export function generatePredictions(buffers, samplesMap, pH, temperature, ionicStrength, referenceOffsets = {}) {
  const allPredictions = {};

  for (const buffer of buffers) {
    const sample = samplesMap.get(buffer.sample_id);
    const predictions = predictBufferShifts(buffer, pH, temperature, ionicStrength, sample, referenceOffsets);

    for (const [nucleus, resonancePredictions] of Object.entries(predictions)) {
      if (!allPredictions[nucleus]) {
        allPredictions[nucleus] = [];
      }

      for (const pred of resonancePredictions) {
        allPredictions[nucleus].push({
          buffer_id: buffer.buffer_id,
          buffer_name: buffer.buffer_name,
          resonance_id: pred.resonance_id,
          description: pred.description,
          predicted_shift: pred.shift
        });
      }
    }
  }

  return allPredictions;
}

/**
 * Assign a single observed shift to the best matching prediction.
 *
 * @param {number} observedShift - Observed chemical shift (ppm)
 * @param {Array<Object>} predictions - Array of prediction objects for this nucleus
 * @param {number} tolerance - Assignment tolerance (ppm)
 * @returns {Object} Assignment object with match details
 */
export function assignSingleShift(observedShift, predictions, tolerance) {
  if (!predictions || predictions.length === 0) {
    return {
      observed_shift: observedShift,
      assigned: false,
      confidence: 'none',
      message: 'No predictions available'
    };
  }

  // Sort predictions by distance to observed shift
  const withDistances = predictions.map(pred => ({
    ...pred,
    distance: observedShift - pred.predicted_shift,
    absDistance: Math.abs(observedShift - pred.predicted_shift)
  }));

  withDistances.sort((a, b) => a.absDistance - b.absDistance);

  const best = withDistances[0];
  const nextBest = withDistances[1];
  const nextBestDistance = nextBest ? nextBest.absDistance : Infinity;

  const confidence = calculateConfidence(best.absDistance, tolerance, nextBestDistance);

  if (confidence === 'none') {
    return {
      observed_shift: observedShift,
      assigned: false,
      confidence: 'none',
      nearest: {
        buffer_id: best.buffer_id,
        buffer_name: best.buffer_name,
        resonance_id: best.resonance_id,
        predicted_shift: best.predicted_shift,
        distance: best.distance
      },
      message: `No prediction within tolerance (nearest: ${best.buffer_name} ${best.resonance_id} at ${best.predicted_shift.toFixed(3)} ppm)`
    };
  }

  return {
    observed_shift: observedShift,
    assigned: true,
    confidence,
    buffer_id: best.buffer_id,
    buffer_name: best.buffer_name,
    resonance_id: best.resonance_id,
    description: best.description,
    predicted_shift: best.predicted_shift,
    residual: best.distance,
    alternatives: confidence !== 'high' && nextBest && nextBest.absDistance < tolerance ? [{
      buffer_id: nextBest.buffer_id,
      buffer_name: nextBest.buffer_name,
      resonance_id: nextBest.resonance_id,
      predicted_shift: nextBest.predicted_shift,
      distance: nextBest.distance
    }] : []
  };
}

/**
 * Assign all observed shifts to predictions using greedy matching.
 * Each prediction can only be assigned once.
 *
 * @param {Object} observedShifts - Object mapping nucleus -> array of observed shifts
 * @param {Array<Object>} buffers - Array of selected buffer objects
 * @param {Map<string, Object>} samplesMap - Map of sample_id to sample object
 * @param {number} pH - pH value for predictions
 * @param {number} temperature - Temperature (K)
 * @param {number} ionicStrength - Ionic strength (M)
 * @param {Object} [tolerances] - Optional custom tolerances by nucleus
 * @param {Object} [referenceOffsets] - Optional reference offsets by nucleus (ppm) to subtract from predictions
 * @returns {Object} Assignment results by nucleus
 */
export function assignPeaks(
  observedShifts,
  buffers,
  samplesMap,
  pH,
  temperature,
  ionicStrength,
  tolerances = {},
  referenceOffsets = {}
) {
  // Generate all predictions (with reference offsets applied)
  const predictions = generatePredictions(buffers, samplesMap, pH, temperature, ionicStrength, referenceOffsets);

  const assignments = {};
  const usedPredictions = new Set();

  for (const [nucleus, shifts] of Object.entries(observedShifts)) {
    const nucleusPredictions = predictions[nucleus] ?? [];
    const tolerance = tolerances[nucleus] ?? DEFAULT_TOLERANCES[nucleus] ?? 1.0;

    assignments[nucleus] = [];

    // Sort observed shifts for consistent ordering
    const sortedShifts = [...shifts].sort((a, b) => a - b);

    for (const observedShift of sortedShifts) {
      // Filter out already-used predictions
      const availablePredictions = nucleusPredictions.filter(
        p => !usedPredictions.has(`${p.buffer_id}:${p.resonance_id}`)
      );

      const assignment = assignSingleShift(observedShift, availablePredictions, tolerance);

      if (assignment.assigned) {
        usedPredictions.add(`${assignment.buffer_id}:${assignment.resonance_id}`);
      }

      assignments[nucleus].push({
        ...assignment,
        nucleus
      });
    }
  }

  return assignments;
}

/**
 * Calculate assignment quality metrics.
 *
 * @param {Object} assignments - Assignment results from assignPeaks
 * @returns {Object} Quality metrics
 */
export function calculateAssignmentQuality(assignments) {
  let totalAssigned = 0;
  let totalUnassigned = 0;
  let highConfidence = 0;
  let mediumConfidence = 0;
  let lowConfidence = 0;
  let totalResidualSquared = 0;
  const residuals = [];

  for (const nucleusAssignments of Object.values(assignments)) {
    for (const assignment of nucleusAssignments) {
      if (assignment.assigned) {
        totalAssigned++;
        residuals.push(assignment.residual);
        totalResidualSquared += assignment.residual ** 2;

        switch (assignment.confidence) {
          case 'high':
            highConfidence++;
            break;
          case 'medium':
            mediumConfidence++;
            break;
          case 'low':
            lowConfidence++;
            break;
        }
      } else {
        totalUnassigned++;
      }
    }
  }

  const n = residuals.length;
  const rmsd = n > 0 ? Math.sqrt(totalResidualSquared / n) : 0;

  return {
    totalObserved: totalAssigned + totalUnassigned,
    totalAssigned,
    totalUnassigned,
    highConfidence,
    mediumConfidence,
    lowConfidence,
    rmsd,
    residuals
  };
}

/**
 * Get flat array of all assignments for fitting.
 *
 * @param {Object} assignments - Assignment results from assignPeaks
 * @returns {Array<Object>} Flat array of assigned peaks
 */
export function getAssignedPeaksForFitting(assignments) {
  const peaks = [];

  for (const [nucleus, nucleusAssignments] of Object.entries(assignments)) {
    for (const assignment of nucleusAssignments) {
      if (assignment.assigned) {
        peaks.push({
          nucleus,
          observed_shift: assignment.observed_shift,
          buffer_id: assignment.buffer_id,
          resonance_id: assignment.resonance_id,
          predicted_shift: assignment.predicted_shift
        });
      }
    }
  }

  return peaks;
}
