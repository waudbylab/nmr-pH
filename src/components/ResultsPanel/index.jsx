import { FittedParameters } from './FittedParameters';
import { AssignmentsTable } from './AssignmentsTable';
import { WarningsDisplay } from './WarningsDisplay';
import { CitationsSection } from './CitationsSection';
import { DownloadButtons } from './DownloadButtons';
import { ReferenceConfigSummary } from '../ReferencingPanel';

/**
 * ResultsPanel component.
 * Combines all result display components.
 */
export function ResultsPanel({
  result,
  validation,
  conditions,
  buffers,
  samplesMap,
  observedShifts,
  nuclei,
  hasDSS,
  dssShift,
  heteroReferencedToDSS,
  spectrometerFreqs,
  fittedReferenceOffsets
}) {
  if (!result) {
    return null;
  }

  return (
    <div className="results-panel">
      <h2>Results</h2>

      <WarningsDisplay result={result} validation={validation} />

      {result.success && (
        <>
          <FittedParameters result={result} nominalConditions={conditions} />

          {nuclei && nuclei.length > 0 && (
            <ReferenceConfigSummary
              nuclei={nuclei}
              hasDSS={hasDSS}
              dssShift={dssShift}
              heteroReferencedToDSS={heteroReferencedToDSS}
              spectrometerFreqs={spectrometerFreqs}
              temperature={conditions.temperature}
              fittedReferenceOffsets={fittedReferenceOffsets}
            />
          )}

          <AssignmentsTable
            assignments={result.assignments}
            residualsDetailed={result.residualsDetailed}
          />
          <CitationsSection buffers={buffers} samplesMap={samplesMap} />
          <DownloadButtons
            result={result}
            conditions={conditions}
            buffers={buffers}
            samplesMap={samplesMap}
            observedShifts={observedShifts}
            nuclei={nuclei}
            hasDSS={hasDSS}
            dssShift={dssShift}
            heteroReferencedToDSS={heteroReferencedToDSS}
            spectrometerFreqs={spectrometerFreqs}
            fittedReferenceOffsets={fittedReferenceOffsets}
          />
        </>
      )}
    </div>
  );
}

export { FittedParameters } from './FittedParameters';
export { AssignmentsTable } from './AssignmentsTable';
export { WarningsDisplay } from './WarningsDisplay';
export { CitationsSection } from './CitationsSection';
export { DownloadButtons } from './DownloadButtons';

export default ResultsPanel;
