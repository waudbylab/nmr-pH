import { useState, useMemo, useEffect, useRef } from 'react';
import { ChemicalShiftPlot } from './ChemicalShiftPlot';
import { ShiftInputArea } from './ShiftInputArea';

/**
 * Mass order for nuclei (lower = shown first).
 */
const NUCLEUS_MASS_ORDER = {
  '1H': 1,
  '13C': 13,
  '15N': 15,
  '19F': 19,
  '31P': 31
};

/**
 * Sort nuclei by mass number.
 */
function sortNucleiByMass(nuclei) {
  return [...nuclei].sort((a, b) => {
    const massA = NUCLEUS_MASS_ORDER[a] ?? parseInt(a.match(/^\d+/)?.[0] || '999');
    const massB = NUCLEUS_MASS_ORDER[b] ?? parseInt(b.match(/^\d+/)?.[0] || '999');
    return massA - massB;
  });
}

/**
 * NucleusTabPanel component.
 * Tabbed interface for each nucleus type with plot and input area.
 */
export function NucleusTabPanel({
  nuclei,
  buffers,
  samplesMap,
  temperature,
  ionicStrength,
  observedShifts,
  onShiftsChange,
  shiftUncertainties = {},
  onShiftUncertaintyChange,
  spectrometerFreqs = {},
  onSpectrometerFreqChange,
  showFrequencyInputs = false,
  fittedPH = null,
  phUncertainty = null,
  assignments = null,
  referenceOffsets = {}, // Reference offsets to subtract from buffer data
  fittedReferenceOffsets = {} // Fitted reference offsets from result (for display only)
}) {
  // Sort nuclei by mass
  const sortedNuclei = useMemo(() => sortNucleiByMass(nuclei), [nuclei]);

  const [activeTab, setActiveTab] = useState(sortedNuclei[0] || null);
  const containerRef = useRef(null);

  // Ensure active tab is valid
  if (activeTab && !sortedNuclei.includes(activeTab)) {
    setActiveTab(sortedNuclei[0] || null);
  }

  // Trigger Plotly resize when tab changes or window resizes
  useEffect(() => {
    // Dispatch a resize event after tab change to update Plotly plots
    const handleResize = () => {
      window.dispatchEvent(new Event('resize'));
    };

    // Small delay to ensure DOM has updated
    const timeoutId = setTimeout(handleResize, 50);

    return () => clearTimeout(timeoutId);
  }, [activeTab]);

  if (sortedNuclei.length === 0) {
    return (
      <div className="nucleus-tab-panel empty">
        <p>Select buffers to view chemical shift data</p>
      </div>
    );
  }

  return (
    <div className="nucleus-tab-panel">
      <div className="tab-header">
        {sortedNuclei.map(nucleus => (
          <button
            key={nucleus}
            className={`tab-button ${activeTab === nucleus ? 'active' : ''}`}
            onClick={() => setActiveTab(nucleus)}
            type="button"
          >
            <sup>{nucleus.match(/^\d+/)?.[0]}</sup>
            {nucleus.replace(/^\d+/, '')}
          </button>
        ))}
      </div>

      <div className="tab-content" ref={containerRef}>
        {sortedNuclei.map(nucleus => (
          <div
            key={nucleus}
            className={`tab-pane ${activeTab === nucleus ? 'active' : 'hidden'}`}
          >
            <div className="nucleus-content">
              <div className="plot-section">
                <ChemicalShiftPlot
                  nucleus={nucleus}
                  buffers={buffers.filter(b =>
                    Object.keys(b.chemical_shifts).includes(nucleus)
                  )}
                  samplesMap={samplesMap}
                  temperature={temperature}
                  ionicStrength={ionicStrength}
                  observedShifts={observedShifts[nucleus] || []}
                  fittedPH={fittedPH}
                  phUncertainty={phUncertainty}
                  assignments={assignments?.[nucleus]}
                  referenceOffset={referenceOffsets[nucleus] ?? 0}
                />
              </div>

              <div className="input-section">
                <ShiftInputArea
                  nucleus={nucleus}
                  value={observedShifts[nucleus] || []}
                  onChange={(shifts) => onShiftsChange(nucleus, shifts)}
                  shiftUncertainty={shiftUncertainties[nucleus]}
                  onShiftUncertaintyChange={onShiftUncertaintyChange}
                  spectrometerFreq={spectrometerFreqs[nucleus]}
                  protonFreq={spectrometerFreqs['1H']}
                  protonReferenceOffset={referenceOffsets['1H'] ?? 0}
                  onSpectrometerFreqChange={onSpectrometerFreqChange}
                  showFrequencyInput={showFrequencyInputs}
                  fittedReferenceOffset={fittedReferenceOffsets[nucleus]}
                />
              </div>
            </div>
          </div>
        ))}
      </div>
    </div>
  );
}

export default NucleusTabPanel;
