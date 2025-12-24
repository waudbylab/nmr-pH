import { useState, useEffect, useCallback } from 'react';
import { calculateSpectrometerFrequency } from './ReferencingPanel';
import { DEFAULT_SHIFT_UNCERTAINTIES } from '../numerical/fitting';

/**
 * Debounce hook.
 */
function useDebounce(value, delay) {
  const [debouncedValue, setDebouncedValue] = useState(value);

  useEffect(() => {
    const timer = setTimeout(() => setDebouncedValue(value), delay);
    return () => clearTimeout(timer);
  }, [value, delay]);

  return debouncedValue;
}

/**
 * Parse chemical shifts from text input.
 * One shift per line, numbers only.
 */
function parseShifts(text) {
  return text
    .split('\n')
    .map(line => line.trim())
    .filter(line => line.length > 0)
    .map(line => parseFloat(line))
    .filter(value => !isNaN(value));
}

/**
 * ShiftInputArea component.
 * Text area for entering observed chemical shifts for a single nucleus,
 * with optional spectrometer frequency input for heteronuclear linking.
 */
export function ShiftInputArea({
  nucleus,
  value,
  onChange,
  shiftUncertainty,
  onShiftUncertaintyChange,
  spectrometerFreq,
  protonFreq,
  onSpectrometerFreqChange,
  showFrequencyInput = false,
  debounceMs = 500
}) {
  const [text, setText] = useState('');
  const defaultUncertainty = DEFAULT_SHIFT_UNCERTAINTIES[nucleus] ?? 0.01;

  // Local state for uncertainty input to allow free typing
  const [uncertaintyText, setUncertaintyText] = useState(
    shiftUncertainty != null ? shiftUncertainty.toString() : ''
  );

  // Initialize text from value
  useEffect(() => {
    if (value && value.length > 0) {
      setText(value.join('\n'));
    }
  }, []); // Only on mount

  const debouncedText = useDebounce(text, debounceMs);

  // Parse and propagate changes after debounce
  useEffect(() => {
    const shifts = parseShifts(debouncedText);
    onChange(shifts);
  }, [debouncedText]); // eslint-disable-line react-hooks/exhaustive-deps

  const handleChange = useCallback((e) => {
    setText(e.target.value);
  }, []);

  // Calculate expected frequency for heteronuclei
  const expectedFreq = nucleus !== '1H' && protonFreq
    ? calculateSpectrometerFrequency(nucleus, protonFreq)
    : null;

  return (
    <div className="shift-input-area">
      <label
        htmlFor={`shifts-${nucleus}`}
        dangerouslySetInnerHTML={{
          __html: `<sup>${nucleus.match(/^\d+/)?.[0] || ''}</sup>${nucleus.replace(/^\d+/, '')} Chemical Shifts (ppm)`
        }}
      />
      <textarea
        id={`shifts-${nucleus}`}
        value={text}
        onChange={handleChange}
        placeholder={`Enter ${nucleus} shifts, one per line:\n2.45\n3.82\n3.15`}
        rows={6}
        spellCheck={false}
      />
      <div className="shift-count">
        {parseShifts(text).length} peaks entered
      </div>

      <div className="shift-uncertainty-input">
        <label>
          <span className="uncertainty-label">
            Shift uncertainty (ppm):
          </span>
          <input
            type="text"
            value={uncertaintyText}
            onChange={(e) => {
              const inputText = e.target.value;
              setUncertaintyText(inputText);
              // Parse and propagate valid non-negative numbers
              const val = parseFloat(inputText);
              if (!isNaN(val) && val >= 0) {
                onShiftUncertaintyChange?.(nucleus, val);
              } else if (inputText === '') {
                onShiftUncertaintyChange?.(nucleus, null);
              }
            }}
            onBlur={() => {
              // On blur, normalize display
              const val = parseFloat(uncertaintyText);
              if (isNaN(val) || val < 0) {
                setUncertaintyText('');
                onShiftUncertaintyChange?.(nucleus, null);
              } else {
                setUncertaintyText(val.toString());
              }
            }}
            placeholder={defaultUncertainty.toString()}
          />
        </label>
        <span className="hint">
          Default: {defaultUncertainty} ppm
        </span>
      </div>

      {showFrequencyInput && onSpectrometerFreqChange && (
        <div className="spectrometer-freq-input">
          <label>
            <span className="freq-label">
              <sup>{nucleus.match(/^\d+/)?.[0]}</sup>
              {nucleus.replace(/^\d+/, '')} spectrometer frequency (MHz):
            </span>
            <input
              type="number"
              value={spectrometerFreq ?? ''}
              onChange={(e) => {
                const val = parseFloat(e.target.value);
                onSpectrometerFreqChange(nucleus, isNaN(val) ? null : val);
              }}
              placeholder={expectedFreq ? expectedFreq.toFixed(3) : 'e.g., 600.13'}
              step="0.001"
            />
          </label>
          {expectedFreq && !spectrometerFreq && (
            <span className="expected-freq hint">
              Expected from ¹H: {expectedFreq.toFixed(3)} MHz
            </span>
          )}
          {nucleus !== '1H' && protonFreq && spectrometerFreq && (
            <span className="linked-status hint">
              Reference linked to ¹H
            </span>
          )}
        </div>
      )}
    </div>
  );
}

export default ShiftInputArea;
