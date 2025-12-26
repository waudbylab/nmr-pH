/**
 * FittedParameters component.
 * Displays fitted parameter values with uncertainties.
 */
export function FittedParameters({ result, nominalConditions }) {
  if (!result || !result.success) {
    return null;
  }

  const { parameters, conditions, statistics } = result;

  // Check for significant deviations from nominal
  const tempDeviation = parameters.temperature
    ? Math.abs(parameters.temperature.value - nominalConditions.temperature)
    : 0;
  const ionicDeviation = parameters.ionicStrength
    ? Math.abs(parameters.ionicStrength.value - nominalConditions.ionicStrength)
    : 0;

  return (
    <div className="fitted-parameters">
      <h3>Fitted Parameters</h3>

      <div className="parameter-list">
        {/* pH - always shown */}
        <div className="parameter-row primary">
          <span className="parameter-name">pH</span>
          <span className="parameter-value">
            {parameters.pH.value.toFixed(2)}
            {parameters.pH.uncertainty > 0 && (
              <span className="uncertainty"> ± {parameters.pH.uncertainty.toFixed(2)}</span>
            )}
          </span>
        </div>

        {/* Temperature */}
        {parameters.temperature ? (
          <div className={`parameter-row ${tempDeviation > 2 ? 'warning' : ''}`}>
            <span className="parameter-name">Temperature</span>
            <span className="parameter-value">
              {parameters.temperature.value.toFixed(1)} K
              {parameters.temperature.uncertainty > 0 && (
                <span className="uncertainty"> ± {parameters.temperature.uncertainty.toFixed(1)} K</span>
              )}
            </span>
            <span className="nominal">
              (nominal: {nominalConditions.temperature.toFixed(1)} K)
            </span>
            {tempDeviation > 2 && (
              <span className="warning-text">
                Differs by {tempDeviation.toFixed(1)} K
              </span>
            )}
          </div>
        ) : (
          <div className="parameter-row fixed">
            <span className="parameter-name">Temperature</span>
            <span className="parameter-value">
              {nominalConditions.temperature.toFixed(1)} K (fixed)
            </span>
          </div>
        )}

        {/* Ionic strength */}
        {parameters.ionicStrength ? (
          <div className={`parameter-row ${ionicDeviation > 0.05 ? 'warning' : ''}`}>
            <span className="parameter-name">Ionic Strength</span>
            <span className="parameter-value">
              {parameters.ionicStrength.value.toFixed(3)} M
              {parameters.ionicStrength.uncertainty > 0 && (
                <span className="uncertainty"> ± {parameters.ionicStrength.uncertainty.toFixed(3)} M</span>
              )}
            </span>
            <span className="nominal">
              (nominal: {nominalConditions.ionicStrength.toFixed(3)} M)
            </span>
            {ionicDeviation > 0.05 && (
              <span className="warning-text">
                Differs by {ionicDeviation.toFixed(3)} M
              </span>
            )}
          </div>
        ) : (
          <div className="parameter-row fixed">
            <span className="parameter-name">Ionic Strength</span>
            <span className="parameter-value">
              {nominalConditions.ionicStrength.toFixed(3)} M (fixed)
            </span>
          </div>
        )}

        {/* Reference offsets */}
        {Object.entries(parameters)
          .filter(([key]) => key.startsWith('ref_'))
          .map(([key, param]) => {
            const nucleus = key.slice(4);
            return (
              <div key={key} className="parameter-row">
                <span className="parameter-name">
                  <sup>{nucleus.match(/^\d+/)?.[0]}</sup>
                  {nucleus.replace(/^\d+/, '')} Reference Offset
                </span>
                <span className="parameter-value">
                  {param.value >= 0 ? '+' : ''}{param.value.toFixed(3)} ppm
                  {param.uncertainty > 0 && (
                    <span className="uncertainty"> ± {param.uncertainty.toFixed(3)} ppm</span>
                  )}
                </span>
              </div>
            );
          })}
      </div>

      {/* Fit statistics */}
      <div className="fit-statistics">
        <h4>Fit Quality</h4>
        <div className="stats-grid">
          <div className="stat">
            <span className="stat-label">RMSD</span>
            <span className="stat-value">{statistics.rmsd.toFixed(4)} ppm</span>
          </div>
          <div className={`stat ${statistics.reducedChiSquared > 2 ? 'warning' : ''}`}>
            <span className="stat-label">χ² / DoF</span>
            <span className="stat-value">
              {statistics.reducedChiSquared.toFixed(2)}
              {statistics.reducedChiSquared > 2 && (
                <span className="chi-warning" title="High reduced chi-square suggests poor fit or underestimated uncertainties">⚠</span>
              )}
            </span>
          </div>
          <div className="stat">
            <span className="stat-label">Observations</span>
            <span className="stat-value">{statistics.nObservations}</span>
          </div>
          <div className="stat">
            <span className="stat-label">Parameters</span>
            <span className="stat-value">{statistics.nParameters}</span>
          </div>
          <div className="stat">
            <span className="stat-label">DoF</span>
            <span className="stat-value">{statistics.degreesOfFreedom}</span>
          </div>
        </div>
        {statistics.reducedChiSquared > 2 && (
          <div className="chi-square-warning">
            High reduced χ² ({statistics.reducedChiSquared.toFixed(2)}) suggests poor fit quality.
            Check for outliers or consider if uncertainties are underestimated.
          </div>
        )}
      </div>
    </div>
  );
}

export default FittedParameters;
