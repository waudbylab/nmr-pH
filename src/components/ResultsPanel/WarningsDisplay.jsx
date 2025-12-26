/**
 * WarningsDisplay component.
 * Shows validation warnings and errors.
 */
export function WarningsDisplay({ result, validation }) {
  const warnings = [];
  const errors = [];

  // Collect errors
  if (!result?.success && result?.error) {
    errors.push(result.error);
  }

  if (validation) {
    // Parameter validation issues
    if (validation.issues) {
      for (const issue of validation.issues) {
        errors.push(issue.message);
      }
    }

    // Warnings
    if (validation.warnings) {
      warnings.push(...validation.warnings);
    }
  }

  // Check for outliers in residuals
  if (validation?.residualCheck?.hasOutliers) {
    for (const outlier of validation.residualCheck.outliers) {
      warnings.push(
        `Large residual for ${outlier.buffer_name} ${outlier.resonance_id}: ` +
        `${outlier.residual.toFixed(3)} ppm (z = ${outlier.zScore.toFixed(1)})`
      );
    }
  }

  // Check for high reduced chi-square
  if (result?.statistics?.reducedChiSquared > 2) {
    warnings.push(
      `High reduced χ² (${result.statistics.reducedChiSquared.toFixed(2)}) suggests poor fit quality or underestimated uncertainties`
    );
  }

  // Include fitting warnings (e.g., parameters disabled due to insufficient DoF)
  if (result?.warnings) {
    warnings.push(...result.warnings);
  }

  if (errors.length === 0 && warnings.length === 0) {
    return null;
  }

  return (
    <div className="warnings-display">
      {errors.length > 0 && (
        <div className="error-section">
          <h4>Errors</h4>
          <ul>
            {errors.map((error, i) => (
              <li key={i} className="error-item">{error}</li>
            ))}
          </ul>
        </div>
      )}

      {warnings.length > 0 && (
        <div className="warning-section">
          <h4>Warnings</h4>
          <ul>
            {warnings.map((warning, i) => (
              <li key={i} className="warning-item">{warning}</li>
            ))}
          </ul>
        </div>
      )}
    </div>
  );
}

export default WarningsDisplay;
