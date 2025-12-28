/**
 * ConditionsPanel component.
 * Input fields for temperature and ionic strength with refinement toggles.
 */
export function ConditionsPanel({
  temperature,
  ionicStrength,
  refineTemperature,
  refineIonicStrength,
  onTemperatureChange,
  onIonicStrengthChange,
  onRefineTemperatureChange,
  onRefineIonicStrengthChange
}) {
  const handleTemperatureChange = (e) => {
    const value = parseFloat(e.target.value);
    if (!isNaN(value)) {
      onTemperatureChange(value);
    }
  };

  const handleIonicStrengthChange = (e) => {
    const value = parseFloat(e.target.value);
    if (!isNaN(value) && value >= 0) {
      onIonicStrengthChange(value);
    }
  };

  return (
    <div className="conditions-panel">
      <h3>Experimental Conditions</h3>

      <div className="condition-row">
        <div className="condition-input">
          <label htmlFor="temperature">Temperature (K)</label>
          <input
            type="number"
            id="temperature"
            value={temperature}
            onChange={handleTemperatureChange}
            min="273"
            max="373"
            step="0.1"
          />
        </div>
        <div className="condition-checkbox">
          <label>
            <input
              type="checkbox"
              checked={refineTemperature}
              onChange={(e) => onRefineTemperatureChange(e.target.checked)}
            />
            Refine during fitting
          </label>
        </div>
      </div>

      <div className="condition-row">
        <div className="condition-input">
          <label htmlFor="ionic-strength">
            Ionic Strength (M)
            <span className="help-icon" title="Ionic strength is calculated as I = 0.5 × Σ(ci × zi²) where ci is the molar concentration and zi is the charge of each ion.&#10;&#10;Examples:&#10;• 50 mM HEPES, 100 mM NaCl: I ≈ 0.10 M&#10;• 100 mM phosphate buffer (pH 7): I ≈ 0.15 M&#10;• 20 mM Tris, 150 mM NaCl: I ≈ 0.15 M">ⓘ</span>
          </label>
          <input
            type="number"
            id="ionic-strength"
            value={ionicStrength}
            onChange={handleIonicStrengthChange}
            min="0"
            max="1"
            step="0.01"
          />
        </div>
        <div className="condition-checkbox">
          <label>
            <input
              type="checkbox"
              checked={refineIonicStrength}
              onChange={(e) => onRefineIonicStrengthChange(e.target.checked)}
            />
            Refine during fitting
          </label>
        </div>
      </div>
    </div>
  );
}

export default ConditionsPanel;
