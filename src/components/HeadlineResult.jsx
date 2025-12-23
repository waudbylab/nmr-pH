import { calculateWaterReference } from './ReferencingPanel';

/**
 * HeadlineResult component.
 * Displays the main pH result prominently with option to view full results.
 */
export function HeadlineResult({
  result,
  calculating,
  onScrollToResults,
  usingAssumedWaterRef = false,
  temperature = 298.15
}) {
  if (calculating) {
    return (
      <div className="headline-result calculating">
        <div className="headline-spinner" />
        <span className="headline-text">Calculating...</span>
      </div>
    );
  }

  if (!result) {
    return null;
  }

  if (!result.success) {
    return (
      <div className="headline-result error">
        <span className="headline-label">Error</span>
        <span className="headline-error-text">{result.error}</span>
      </div>
    );
  }

  const pH = result.parameters.pH.value;
  const uncertainty = result.parameters.pH.uncertainty;
  const waterRef = calculateWaterReference(temperature);

  return (
    <div className={`headline-result ${usingAssumedWaterRef ? 'warning' : 'success'}`}>
      <div className="headline-ph">
        <span className="headline-label">Fitted pH</span>
        <span className="headline-value">
          {pH.toFixed(2)}
          {uncertainty != null && (
            <span className="headline-uncertainty">
              {' '}± {uncertainty.toFixed(2)}
            </span>
          )}
        </span>
        {usingAssumedWaterRef && (
          <span className="headline-warning">
            Using assumed water reference ({waterRef.toFixed(3)} ppm at {temperature.toFixed(1)} K)
          </span>
        )}
      </div>
      <button
        type="button"
        className="view-details-button"
        onClick={onScrollToResults}
      >
        View full results ↓
      </button>
    </div>
  );
}

export default HeadlineResult;
