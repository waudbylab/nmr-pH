/**
 * Validation Module
 *
 * Checks degrees of freedom, extrapolation warnings, and fit quality.
 */

/**
 * Validate degrees of freedom for fitting based on the new referencing model.
 *
 * DOF Calculation:
 * For each nucleus:
 *   effective_shifts[nucleus] = n_shifts[nucleus] - (1 if needs_ref_fitting else 0)
 *
 * total_dof = sum(effective_shifts) - n_fitting_params
 *
 * where n_fitting_params includes:
 *   - 1 (pH, always)
 *   + (1 if refining temperature)
 *   + (1 if refining ionic strength)
 *   (reference offsets are already accounted for in effective_shifts)
 *
 * If DoF is insufficient, the function will suggest reduced options by
 * disabling ionic strength first, then temperature, to achieve DoF >= 0.
 *
 * @param {Object} shiftCounts - Object mapping nucleus -> number of observed shifts
 * @param {Object} refineReferences - Object mapping nucleus -> boolean (whether to fit reference)
 * @param {Array<string>} linkedToProton - Nuclei linked to 1H (don't need fitting)
 * @param {boolean} refineTemperature - Whether temperature is being refined
 * @param {boolean} refineIonicStrength - Whether ionic strength is being refined
 * @returns {Object} Validation result with suggested reduced options if needed
 */
export function validateDegreesOfFreedom(
  shiftCounts,
  refineReferences,
  linkedToProton = [],
  refineTemperature = false,
  refineIonicStrength = false
) {
  const details = [];
  let totalEffectiveShifts = 0;

  // Calculate effective shifts per nucleus
  for (const [nucleus, nShifts] of Object.entries(shiftCounts)) {
    if (nShifts === 0) continue;

    const needsRefFitting = refineReferences[nucleus] && !linkedToProton.includes(nucleus);
    const effectiveShifts = nShifts - (needsRefFitting ? 1 : 0);

    if (needsRefFitting) {
      details.push(`${nucleus}: ${nShifts} shift(s) - 1 for reference = ${effectiveShifts} effective`);
    } else {
      details.push(`${nucleus}: ${nShifts} shift(s) = ${effectiveShifts} effective`);
    }

    totalEffectiveShifts += effectiveShifts;
  }

  // Count fitting parameters (excluding reference offsets, which are in effective_shifts)
  let nFittingParams = 1; // pH always
  const paramDetails = ['pH'];

  if (refineTemperature) {
    nFittingParams++;
    paramDetails.push('Temperature');
  }

  if (refineIonicStrength) {
    nFittingParams++;
    paramDetails.push('Ionic strength');
  }

  details.push(`Fitting parameters: ${paramDetails.join(', ')}`);

  // Calculate total DOF
  const totalDof = totalEffectiveShifts - nFittingParams;
  details.push(`Degrees of freedom: ${totalEffectiveShifts} - ${nFittingParams} = ${totalDof}`);

  // Determine required DOF
  // 0 DOF is acceptable (exactly determined system)
  const requiredDof = 0;

  // Check if valid with current settings
  let valid = totalDof >= requiredDof;
  let message = '';
  const warnings = [];

  // If not valid, try reducing parameters
  let effectiveRefineTemperature = refineTemperature;
  let effectiveRefineIonicStrength = refineIonicStrength;
  let effectiveNFittingParams = nFittingParams;
  let effectiveTotalDof = totalDof;

  if (!valid) {
    // First try: disable ionic strength refinement
    if (refineIonicStrength) {
      effectiveRefineIonicStrength = false;
      effectiveNFittingParams--;
      effectiveTotalDof = totalEffectiveShifts - effectiveNFittingParams;
      warnings.push('Ionic strength refinement disabled due to insufficient degrees of freedom');

      if (effectiveTotalDof >= requiredDof) {
        valid = true;
      }
    }

    // Second try: also disable temperature refinement
    if (!valid && refineTemperature) {
      effectiveRefineTemperature = false;
      effectiveNFittingParams--;
      effectiveTotalDof = totalEffectiveShifts - effectiveNFittingParams;
      warnings.push('Temperature refinement disabled due to insufficient degrees of freedom');

      if (effectiveTotalDof >= requiredDof) {
        valid = true;
      }
    }

    // Still not valid - truly underdetermined
    if (!valid) {
      message = `Underdetermined system: need more chemical shifts. Current DOF: ${effectiveTotalDof}`;
    }
  }

  return {
    valid,
    totalEffectiveShifts,
    nFittingParams: effectiveNFittingParams,
    totalDof: effectiveTotalDof,
    requiredDof,
    message,
    details,
    warnings,
    // Return effective settings for use by caller
    effectiveRefineTemperature,
    effectiveRefineIonicStrength
  };
}

/**
 * Check if sufficient chemical shifts are provided for each nucleus based on referencing mode.
 *
 * Scenarios:
 * 1. DSS-referenced: 1 shift minimum (no reference offset to fit)
 * 2. Not referenced, independent offset: 2 shifts minimum (1 for pH, 1 for reference)
 * 3. Linked to 1H (via spectrometer frequency): treated as DSS-referenced (1 shift minimum)
 *
 * @param {Object} observedShifts - Object mapping nucleus -> array of shifts
 * @param {Object} referenceConfigs - Configuration for each nucleus
 * @param {Array<string>} linkedToProton - Nuclei whose reference is linked to 1H
 * @returns {Object} Validation result
 */
export function validateMinimumShifts(observedShifts, referenceConfigs, linkedToProton = []) {
  const issues = [];
  let valid = true;

  for (const [nucleus, shifts] of Object.entries(observedShifts)) {
    const config = referenceConfigs[nucleus] || { mode: 'not_referenced', refineOffset: true };
    const nShifts = shifts.length;

    let minRequired;
    let reason;

    if (config.mode === 'referenced') {
      // Scenario 2: Referenced to DSS
      minRequired = 1;
      reason = 'referenced to DSS';
    } else if (linkedToProton.includes(nucleus)) {
      // Scenario 3: Linked to 1H via spectrometer frequency
      minRequired = 1;
      reason = 'reference linked to ¹H';
    } else {
      // Scenario 1 or 4: Independent reference offset needs fitting
      minRequired = 2;
      reason = 'independent reference offset';
    }

    if (nShifts < minRequired) {
      valid = false;
      issues.push({
        nucleus,
        nShifts,
        minRequired,
        message: `${nucleus}: ${nShifts} shift(s) provided, but ${minRequired} required (${reason})`
      });
    }
  }

  return {
    valid,
    issues,
    summary: valid
      ? 'Sufficient shifts for all nuclei'
      : `Insufficient shifts: ${issues.map(i => i.nucleus).join(', ')}`
  };
}

/**
 * Check if system has sufficient degrees of freedom.
 *
 * @param {number} nObservations - Number of observations
 * @param {number} nParameters - Number of parameters to fit
 * @returns {Object} Validation result
 */
export function checkDegreesOfFreedom(nObservations, nParameters) {
  const dof = nObservations - nParameters;

  return {
    valid: dof > 0,
    warning: dof === 1,
    nObservations,
    nParameters,
    degreesOfFreedom: dof,
    message: dof <= 0
      ? `Underdetermined system: need at least ${nParameters + 1} observations for ${nParameters} parameters`
      : dof === 1
        ? 'Marginal degrees of freedom (DoF = 1). Consider adding more data or fixing parameters.'
        : `Degrees of freedom: ${dof}`
  };
}

/**
 * Count parameters based on fitting options.
 *
 * @param {Object} options - Fitting options
 * @returns {Object} Parameter counts
 */
export function countParameters(options) {
  let count = 1; // pH is always fitted
  const details = ['pH'];

  if (options.refineTemperature) {
    count++;
    details.push('Temperature');
  }

  if (options.refineIonicStrength) {
    count++;
    details.push('Ionic strength');
  }

  for (const [nucleus, refine] of Object.entries(options.refineReferences ?? {})) {
    if (refine) {
      count++;
      details.push(`${nucleus} reference`);
    }
  }

  return { count, details };
}

/**
 * Check if conditions are within measurement ranges.
 *
 * @param {Object} conditions - Current conditions (pH, temperature, ionicStrength)
 * @param {Array<Object>} samples - Sample objects with measurement_ranges
 * @returns {Object} Extrapolation warnings
 */
export function checkExtrapolation(conditions, samples) {
  const warnings = [];

  for (const sample of samples) {
    const ranges = sample.measurement_ranges;
    if (!ranges) continue;

    const sampleWarnings = [];

    // Check pH
    if (ranges.pH) {
      if (conditions.pH < ranges.pH.min) {
        sampleWarnings.push(`pH ${conditions.pH.toFixed(2)} is below measured range (min: ${ranges.pH.min})`);
      } else if (conditions.pH > ranges.pH.max) {
        sampleWarnings.push(`pH ${conditions.pH.toFixed(2)} is above measured range (max: ${ranges.pH.max})`);
      }
    }

    // Check temperature
    if (ranges.temperature_K) {
      if (conditions.temperature < ranges.temperature_K.min) {
        sampleWarnings.push(`Temperature ${conditions.temperature.toFixed(1)} K is below measured range (min: ${ranges.temperature_K.min} K)`);
      } else if (conditions.temperature > ranges.temperature_K.max) {
        sampleWarnings.push(`Temperature ${conditions.temperature.toFixed(1)} K is above measured range (max: ${ranges.temperature_K.max} K)`);
      }
    }

    // Check ionic strength
    if (ranges.ionic_strength_M) {
      if (conditions.ionicStrength < ranges.ionic_strength_M.min) {
        sampleWarnings.push(`Ionic strength ${conditions.ionicStrength.toFixed(3)} M is below measured range (min: ${ranges.ionic_strength_M.min} M)`);
      } else if (conditions.ionicStrength > ranges.ionic_strength_M.max) {
        sampleWarnings.push(`Ionic strength ${conditions.ionicStrength.toFixed(3)} M is above measured range (max: ${ranges.ionic_strength_M.max} M)`);
      }
    }

    if (sampleWarnings.length > 0) {
      warnings.push({
        sample_id: sample.sample_id,
        warnings: sampleWarnings
      });
    }
  }

  return {
    hasWarnings: warnings.length > 0,
    warnings
  };
}

/**
 * Validate assignment quality.
 *
 * @param {Object} assignments - Assignment results from peakAssignment
 * @returns {Object} Quality assessment
 */
export function validateAssignments(assignments) {
  const issues = [];
  let totalAssigned = 0;
  let totalUnassigned = 0;
  let lowConfidenceCount = 0;
  let ambiguousCount = 0;

  for (const [nucleus, nucleusAssignments] of Object.entries(assignments)) {
    for (const assignment of nucleusAssignments) {
      if (assignment.assigned) {
        totalAssigned++;

        if (assignment.confidence === 'low') {
          lowConfidenceCount++;
          issues.push({
            type: 'low_confidence',
            nucleus,
            observed_shift: assignment.observed_shift,
            message: `Low confidence assignment: ${assignment.observed_shift.toFixed(3)} ppm → ${assignment.buffer_name} ${assignment.resonance_id}`
          });
        }

        if (assignment.alternatives && assignment.alternatives.length > 0) {
          ambiguousCount++;
          issues.push({
            type: 'ambiguous',
            nucleus,
            observed_shift: assignment.observed_shift,
            message: `Ambiguous assignment: ${assignment.observed_shift.toFixed(3)} ppm could match multiple resonances`
          });
        }
      } else {
        totalUnassigned++;
        issues.push({
          type: 'unassigned',
          nucleus,
          observed_shift: assignment.observed_shift,
          message: assignment.message || `Unassigned peak at ${assignment.observed_shift.toFixed(3)} ppm`
        });
      }
    }
  }

  return {
    valid: totalAssigned > 0,
    totalAssigned,
    totalUnassigned,
    lowConfidenceCount,
    ambiguousCount,
    issues,
    summary: totalUnassigned > 0
      ? `${totalAssigned} peaks assigned, ${totalUnassigned} unassigned`
      : `All ${totalAssigned} peaks assigned`
  };
}

/**
 * Validate residuals and identify outliers.
 *
 * @param {Array<Object>} assignments - Flat array of assignments with residuals
 * @param {number} rmsd - Root mean square deviation
 * @param {number} [threshold=2] - Z-score threshold for outliers
 * @returns {Object} Residual validation results
 */
export function validateResiduals(assignments, rmsd, threshold = 2) {
  const outliers = [];
  const residualStats = {
    min: Infinity,
    max: -Infinity,
    mean: 0,
    count: 0
  };

  const flatAssignments = [];

  for (const [nucleus, nucleusAssignments] of Object.entries(assignments)) {
    for (const assignment of nucleusAssignments) {
      if (assignment.assigned && assignment.residual !== undefined) {
        flatAssignments.push({ ...assignment, nucleus });
        residualStats.count++;
        residualStats.mean += assignment.residual;
        residualStats.min = Math.min(residualStats.min, assignment.residual);
        residualStats.max = Math.max(residualStats.max, assignment.residual);
      }
    }
  }

  if (residualStats.count > 0) {
    residualStats.mean /= residualStats.count;
  }

  // Identify outliers based on z-score
  if (rmsd > 0) {
    for (const assignment of flatAssignments) {
      const zScore = Math.abs(assignment.residual) / rmsd;
      if (zScore > threshold) {
        outliers.push({
          nucleus: assignment.nucleus,
          observed_shift: assignment.observed_shift,
          buffer_name: assignment.buffer_name,
          resonance_id: assignment.resonance_id,
          residual: assignment.residual,
          zScore
        });
      }
    }
  }

  return {
    hasOutliers: outliers.length > 0,
    outliers,
    statistics: residualStats,
    rmsd
  };
}

/**
 * Check for physically unrealistic parameter values.
 *
 * @param {Object} parameters - Fitted parameters
 * @returns {Object} Validation result
 */
export function validateParameters(parameters) {
  const issues = [];

  if (parameters.pH) {
    const pH = parameters.pH.value;
    if (pH < 0 || pH > 14) {
      issues.push({
        parameter: 'pH',
        value: pH,
        message: `pH ${pH.toFixed(2)} is outside physical range (0-14)`
      });
    }
  }

  if (parameters.temperature) {
    const T = parameters.temperature.value;
    if (T < 273) {
      issues.push({
        parameter: 'temperature',
        value: T,
        message: `Temperature ${T.toFixed(1)} K is below freezing point of water`
      });
    } else if (T > 373) {
      issues.push({
        parameter: 'temperature',
        value: T,
        message: `Temperature ${T.toFixed(1)} K is above boiling point of water`
      });
    }
  }

  if (parameters.ionicStrength) {
    const I = parameters.ionicStrength.value;
    if (I < 0) {
      issues.push({
        parameter: 'ionicStrength',
        value: I,
        message: `Ionic strength ${I.toFixed(3)} M is negative (physically impossible)`
      });
    } else if (I > 1) {
      issues.push({
        parameter: 'ionicStrength',
        value: I,
        message: `Ionic strength ${I.toFixed(3)} M is very high (may be unreliable)`
      });
    }
  }

  return {
    valid: issues.length === 0,
    issues
  };
}

/**
 * Compare fitted values to nominal values and flag large deviations.
 *
 * @param {Object} fitted - Fitted conditions
 * @param {Object} nominal - Nominal conditions
 * @returns {Object} Deviation warnings
 */
export function checkDeviations(fitted, nominal) {
  const deviations = [];

  if (nominal.temperature && fitted.temperature) {
    const diff = Math.abs(fitted.temperature - nominal.temperature);
    if (diff > 2) {
      deviations.push({
        parameter: 'temperature',
        nominal: nominal.temperature,
        fitted: fitted.temperature,
        difference: diff,
        message: `Refined temperature differs from nominal by ${diff.toFixed(1)} K`
      });
    }
  }

  if (nominal.ionicStrength !== undefined && fitted.ionicStrength !== undefined) {
    const diff = Math.abs(fitted.ionicStrength - nominal.ionicStrength);
    if (diff > 0.05) {
      deviations.push({
        parameter: 'ionicStrength',
        nominal: nominal.ionicStrength,
        fitted: fitted.ionicStrength,
        difference: diff,
        message: `Refined ionic strength differs from nominal by ${diff.toFixed(3)} M`
      });
    }
  }

  return {
    hasDeviations: deviations.length > 0,
    deviations
  };
}

/**
 * Perform comprehensive validation of fitting results.
 *
 * @param {Object} result - Fitting result from fitting.js
 * @param {Object} nominalConditions - Nominal conditions for comparison
 * @param {Array<Object>} samples - Sample objects for extrapolation check
 * @returns {Object} Comprehensive validation report
 */
export function validateFitResult(result, nominalConditions, samples) {
  if (!result.success) {
    return {
      valid: false,
      error: result.error,
      warnings: [],
      issues: []
    };
  }

  const warnings = [];
  const issues = [];

  // Check parameters
  const paramValidation = validateParameters(result.parameters);
  if (!paramValidation.valid) {
    issues.push(...paramValidation.issues);
  }

  // Check deviations from nominal
  const deviationCheck = checkDeviations(result.conditions, nominalConditions);
  if (deviationCheck.hasDeviations) {
    warnings.push(...deviationCheck.deviations.map(d => d.message));
  }

  // Check extrapolation
  const extrapolationCheck = checkExtrapolation(result.conditions, samples);
  if (extrapolationCheck.hasWarnings) {
    for (const warning of extrapolationCheck.warnings) {
      warnings.push(...warning.warnings);
    }
  }

  // Validate residuals
  const residualCheck = validateResiduals(result.assignments, result.statistics.rmsd);
  if (residualCheck.hasOutliers) {
    for (const outlier of residualCheck.outliers) {
      warnings.push(`Outlier: ${outlier.buffer_name} ${outlier.resonance_id} has z-score ${outlier.zScore.toFixed(1)}`);
    }
  }

  return {
    valid: issues.length === 0,
    warnings,
    issues,
    statistics: result.statistics,
    residualCheck,
    parameterCheck: paramValidation,
    extrapolationCheck,
    deviationCheck
  };
}
