/**
 * ReferencingPanel component.
 * Configure chemical shift referencing with DSS detection workflow.
 */

import { useState, useEffect } from 'react';

/**
 * Updated IUPAC Xi ratios relative to DSS (from BMRB/Iowa State).
 * These are the fractional frequencies relative to 1H.
 */
export const XI_RATIOS = {
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
 * Calculate the spectrometer frequency for nucleus X given 1H frequency.
 * bf(X) = bf(1H) * Xi(X)
 */
export function calculateSpectrometerFrequency(nucleus, protonFrequencyMHz) {
  const xiRatio = XI_RATIOS[nucleus];
  if (!xiRatio) return null;
  return protonFrequencyMHz * xiRatio;
}

/**
 * Calculate the reference offset for nucleus X from 1H reference offset.
 * Given:
 *   - bf(1H) and bf(X): spectrometer frequencies
 *   - delta_ref(1H): 1H reference offset in ppm
 *
 * Algorithm:
 *   1. Calculate frequency of DSS for 1H: nu_DSS(1H) = bf(1H) * (1 - delta_ref(1H)/1e6)
 *   2. Calculate true zero frequency for X: nu_0(X) = nu_DSS(1H) * Xi(X) / Xi(1H)
 *   3. Calculate X reference offset: delta_ref(X) = (bf(X) - nu_0(X)) / bf(X) * 1e6
 */
export function calculateLinkedReferenceOffset(nucleus, protonFrequencyMHz, protonReferenceOffsetPpm) {
  if (nucleus === '1H') return protonReferenceOffsetPpm;

  const xiX = XI_RATIOS[nucleus];
  const xiH = XI_RATIOS['1H'];

  if (!xiX) return 0;

  // Step 1: Calculate DSS frequency for 1H
  const nuDSS_H = protonFrequencyMHz * (1 - protonReferenceOffsetPpm / 1e6);

  // Step 2: Calculate true zero frequency for X
  const nu0_X = nuDSS_H * (xiX / xiH);

  // Step 3: Calculate bf(X)
  const bfX = protonFrequencyMHz * (xiX / xiH);

  // Step 4: Calculate X reference offset
  const deltaRef_X = ((bfX - nu0_X) / bfX) * 1e6;

  return deltaRef_X;
}

/**
 * Step 1: DSS Detection question - compact inline layout
 */
function DSSDetectionStep({
  hasDSS,
  dssShift,
  onHasDSSChange,
  onDSSShiftChange
}) {
  // Use local state for text input to allow typing negative numbers freely
  const [dssText, setDssText] = useState(dssShift?.toString() ?? '0');
  // Track if user is actively editing to prevent external sync during typing
  const [isEditing, setIsEditing] = useState(false);

  // Sync local state when prop changes externally (but not while editing)
  useEffect(() => {
    if (!isEditing) {
      setDssText(dssShift?.toString() ?? '0');
    }
  }, [dssShift, isEditing]);

  const handleDssTextChange = (e) => {
    const text = e.target.value;
    setDssText(text);
    setIsEditing(true);

    // Parse and update if valid number (including partial input like "-" or "-0.")
    const parsed = parseFloat(text);
    if (!isNaN(parsed)) {
      onDSSShiftChange(parsed);
    }
  };

  const handleDssBlur = () => {
    setIsEditing(false);
    // On blur, normalize the display and ensure valid value
    const parsed = parseFloat(dssText);
    if (isNaN(parsed)) {
      setDssText('0');
      onDSSShiftChange(0);
    } else {
      setDssText(parsed.toString());
      onDSSShiftChange(parsed);
    }
  };

  return (
    <div className="referencing-step compact">
      <div className="step-row">
        <span className="step-question">Does your sample contain DSS?</span>
        <div className="radio-group inline">
          <label>
            <input
              type="radio"
              name="hasDSS"
              checked={hasDSS === true}
              onChange={() => onHasDSSChange(true)}
            />
            Yes
          </label>
          <label>
            <input
              type="radio"
              name="hasDSS"
              checked={hasDSS === false}
              onChange={() => onHasDSSChange(false)}
            />
            No
          </label>
        </div>
        {hasDSS === true && (
          <div className="dss-shift-input inline">
            <label>
              DSS shift:
              <input
                type="text"
                value={dssText}
                onChange={handleDssTextChange}
                onBlur={handleDssBlur}
                placeholder="0.00"
              />
              ppm
            </label>
          </div>
        )}
      </div>
    </div>
  );
}

/**
 * Step 2: Heteronuclear DSS Referencing question - compact inline layout
 */
function HeteronuclearReferencingStep({
  heteroReferencedToDSS,
  onHeteroReferencedChange
}) {
  return (
    <div className="referencing-step compact">
      <div className="step-row">
        <span className="step-question">Are heteronuclear spectra referenced to DSS?</span>
        <div className="radio-group inline">
          <label>
            <input
              type="radio"
              name="heteroRef"
              checked={heteroReferencedToDSS === true}
              onChange={() => onHeteroReferencedChange(true)}
            />
            Yes
          </label>
          <label>
            <input
              type="radio"
              name="heteroRef"
              checked={heteroReferencedToDSS === false}
              onChange={() => onHeteroReferencedChange(false)}
            />
            No
          </label>
        </div>
      </div>
    </div>
  );
}

/**
 * Brief message about spectrometer frequencies (replaces Step 3)
 */
function SpectrometerFrequencyMessage({ hasDSS }) {
  if (hasDSS) {
    return (
      <div className="referencing-message">
        <span className="hint">
          Enter spectrometer frequencies below each nucleus to link heteronuclear references to <sup>1</sup>H.
          Without frequencies, heteronuclear references will be fitted independently (±5 ppm).
        </span>
      </div>
    );
  } else {
    return (
      <div className="referencing-message">
        <span className="hint">
          The <sup>1</sup>H reference will be fitted using the water signal. Enter spectrometer frequencies
          to link heteronuclear references to <sup>1</sup>H, or leave blank to fit independently.
        </span>
      </div>
    );
  }
}

/**
 * Determine referencing state for each nucleus based on user configuration.
 * This is used by App.jsx to set up fitting parameters.
 *
 * @param {Array<string>} nuclei - List of nuclei present
 * @param {boolean} hasDSS - Whether sample contains DSS
 * @param {number} dssShift - DSS chemical shift (if present)
 * @param {boolean} heteroReferencedToDSS - Whether heteronuclei are referenced to DSS
 * @param {Object} spectrometerFreqs - Object mapping nucleus -> frequency in MHz
 * @param {number} temperature - Temperature in K (for water reference guess)
 * @returns {Object} Referencing configuration for fitting
 */
export function buildReferencingConfig(nuclei, hasDSS, dssShift, heteroReferencedToDSS, spectrometerFreqs, temperature) {
  const config = {
    // Per-nucleus referencing info
    nucleusConfigs: {},
    // Which nuclei have references linked to 1H
    linkedToProton: [],
    // Initial reference offsets
    referenceOffsets: {},
    // Which reference offsets to refine
    refineReferences: {},
    // Parameter bounds for reference offsets
    referenceBounds: {},
    // Proton frequency for linked calculations
    protonFrequency: spectrometerFreqs?.['1H'] || null
  };

  for (const nucleus of nuclei) {
    if (nucleus === '1H') {
      if (hasDSS) {
        // 1H is referenced to DSS
        config.nucleusConfigs['1H'] = { mode: 'referenced', dssShift: dssShift ?? 0 };
        config.referenceOffsets['1H'] = dssShift ?? 0;
        config.refineReferences['1H'] = false;
      } else {
        // 1H not referenced - use water as initial guess
        const waterRef = calculateWaterReference(temperature);
        config.nucleusConfigs['1H'] = { mode: 'not_referenced', initialGuess: waterRef };
        config.referenceOffsets['1H'] = waterRef;
        config.refineReferences['1H'] = true;
        config.referenceBounds['1H'] = { min: waterRef - 0.5, max: waterRef + 0.5 };
      }
    } else {
      // Heteronucleus
      if (hasDSS && heteroReferencedToDSS) {
        // X is directly referenced to DSS
        config.nucleusConfigs[nucleus] = { mode: 'referenced', dssShift: 0 };
        config.referenceOffsets[nucleus] = 0;
        config.refineReferences[nucleus] = false;
      } else if (spectrometerFreqs?.['1H'] && spectrometerFreqs?.[nucleus]) {
        // X is linked to 1H via spectrometer frequency
        config.nucleusConfigs[nucleus] = { mode: 'linked', linkedTo: '1H' };
        config.linkedToProton.push(nucleus);
        // Reference offset will be calculated dynamically from 1H
        config.referenceOffsets[nucleus] = 0; // Placeholder
        config.refineReferences[nucleus] = false;
      } else {
        // X has independent reference that needs fitting
        config.nucleusConfigs[nucleus] = { mode: 'not_referenced', initialGuess: 0 };
        config.referenceOffsets[nucleus] = 0;
        config.refineReferences[nucleus] = true;
        config.referenceBounds[nucleus] = { min: -5, max: 5 };
      }
    }
  }

  return config;
}

/**
 * ReferencingPanel component - simplified compact version.
 */
export function ReferencingPanel({
  nuclei,
  hasDSS,
  dssShift,
  heteroReferencedToDSS,
  temperature,
  onHasDSSChange,
  onDSSShiftChange,
  onHeteroReferencedChange
}) {
  if (nuclei.length === 0) {
    return null;
  }

  const hasHeteronuclei = nuclei.some(n => n !== '1H');
  const showStep2 = hasDSS === true && hasHeteronuclei;
  const showFreqMessage = hasHeteronuclei && hasDSS !== null && !(hasDSS && heteroReferencedToDSS);

  return (
    <div className="referencing-panel compact">
      <h3>Chemical Shift Referencing</h3>

      <DSSDetectionStep
        hasDSS={hasDSS}
        dssShift={dssShift}
        onHasDSSChange={onHasDSSChange}
        onDSSShiftChange={onDSSShiftChange}
      />

      {showStep2 && (
        <HeteronuclearReferencingStep
          heteroReferencedToDSS={heteroReferencedToDSS}
          onHeteroReferencedChange={onHeteroReferencedChange}
        />
      )}

      {showFreqMessage && (
        <SpectrometerFrequencyMessage hasDSS={hasDSS} />
      )}
    </div>
  );
}

/**
 * Reference configuration summary for results panel.
 */
export function ReferenceConfigSummary({
  nuclei,
  hasDSS,
  dssShift,
  heteroReferencedToDSS,
  spectrometerFreqs,
  temperature,
  fittedReferenceOffsets
}) {
  const getStatusForNucleus = (nucleus) => {
    if (nucleus === '1H') {
      if (hasDSS) {
        return { status: 'Referenced', detail: `DSS at ${dssShift ?? 0} ppm` };
      } else {
        const fittedOffset = fittedReferenceOffsets?.['1H'];
        if (fittedOffset !== undefined) {
          return { status: 'Fitted', detail: `${fittedOffset.toFixed(3)} ppm` };
        }
        return { status: 'Fitting', detail: 'from water signal' };
      }
    } else {
      if (hasDSS && heteroReferencedToDSS) {
        return { status: 'Referenced', detail: 'DSS at 0 ppm' };
      } else if (spectrometerFreqs?.['1H'] && spectrometerFreqs?.[nucleus]) {
        const fittedOffset = fittedReferenceOffsets?.[nucleus];
        if (fittedOffset !== undefined) {
          return { status: 'Linked', detail: `${fittedOffset.toFixed(3)} ppm (from ¹H)` };
        }
        return { status: 'Linked', detail: 'to ¹H' };
      } else {
        const fittedOffset = fittedReferenceOffsets?.[nucleus];
        if (fittedOffset !== undefined) {
          return { status: 'Fitted', detail: `${fittedOffset.toFixed(3)} ppm` };
        }
        return { status: 'Fitting', detail: 'independently' };
      }
    }
  };

  return (
    <div className="reference-config-summary">
      <h4>Reference Configuration</h4>
      <table className="summary-table compact">
        <thead>
          <tr>
            <th>Nucleus</th>
            <th>Status</th>
            <th>Offset</th>
          </tr>
        </thead>
        <tbody>
          {nuclei.map(nucleus => {
            const info = getStatusForNucleus(nucleus);
            return (
              <tr key={nucleus}>
                <td>
                  <sup>{nucleus.match(/^\d+/)?.[0]}</sup>
                  {nucleus.replace(/^\d+/, '')}
                </td>
                <td>{info.status}</td>
                <td>{info.detail}</td>
              </tr>
            );
          })}
        </tbody>
      </table>
    </div>
  );
}

export default ReferencingPanel;
