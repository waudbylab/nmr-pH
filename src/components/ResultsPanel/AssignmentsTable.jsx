/**
 * AssignmentsTable component.
 * Displays peak assignments in a table format with z-scores.
 */
export function AssignmentsTable({ assignments, residualsDetailed }) {
  if (!assignments) {
    return null;
  }

  // Build lookup map from residualsDetailed for z-scores
  const zScoreMap = new Map();
  if (residualsDetailed) {
    for (const r of residualsDetailed) {
      // Key by nucleus + observed shift (rounded to match)
      const key = `${r.nucleus}:${r.observed.toFixed(4)}`;
      zScoreMap.set(key, r);
    }
  }

  // Flatten assignments into a single array
  const rows = [];
  for (const [nucleus, nucleusAssignments] of Object.entries(assignments)) {
    for (const assignment of nucleusAssignments) {
      const key = `${nucleus}:${assignment.observed_shift.toFixed(4)}`;
      const detailedResidual = zScoreMap.get(key);
      rows.push({
        nucleus,
        ...assignment,
        zScore: detailedResidual?.zScore ?? null,
        sigmaTotal: detailedResidual?.sigmaTotal ?? null
      });
    }
  }

  if (rows.length === 0) {
    return null;
  }

  // Count outliers (|z| > 2)
  const outlierCount = rows.filter(r => r.zScore !== null && Math.abs(r.zScore) > 2).length;

  return (
    <div className="assignments-table">
      <h3>Peak Assignments</h3>

      <table>
        <thead>
          <tr>
            <th>Nucleus</th>
            <th>Observed (ppm)</th>
            <th>Assigned To</th>
            <th>Predicted (ppm)</th>
            <th>Residual (ppm)</th>
            <th>z-score</th>
          </tr>
        </thead>
        <tbody>
          {rows.map((row, index) => {
            const isOutlier = row.zScore !== null && Math.abs(row.zScore) > 2;
            return (
              <tr
                key={index}
                className={`
                  ${!row.assigned ? 'unassigned' : ''}
                  ${isOutlier ? 'outlier' : ''}
                `}
              >
                <td>
                  <sup>{row.nucleus.match(/^\d+/)?.[0]}</sup>
                  {row.nucleus.replace(/^\d+/, '')}
                </td>
                <td>{row.observed_shift.toFixed(3)}</td>
                <td>
                  {row.assigned ? (
                    <>
                      {row.buffer_name}
                      <br />
                      <small>{row.resonance_id}</small>
                    </>
                  ) : (
                    <span className="unassigned-text">Unassigned</span>
                  )}
                </td>
                <td>
                  {row.assigned ? row.predicted_shift.toFixed(3) : '-'}
                </td>
                <td>
                  {row.assigned ? (
                    <>
                      {row.residual >= 0 ? '+' : ''}{row.residual.toFixed(3)}
                    </>
                  ) : '-'}
                </td>
                <td className={isOutlier ? 'outlier-zscore' : ''}>
                  {row.zScore !== null ? (
                    <>
                      {row.zScore >= 0 ? '+' : ''}{row.zScore.toFixed(2)}
                      {isOutlier && <span className="outlier-warning" title="Large residual">⚠</span>}
                    </>
                  ) : (
                    row.assigned ? '-' : '-'
                  )}
                </td>
              </tr>
            );
          })}
        </tbody>
      </table>

      {/* Summary */}
      <div className="assignments-summary">
        <span className="assigned-count">
          {rows.filter(r => r.assigned).length} assigned
        </span>
        {rows.some(r => !r.assigned) && (
          <span className="unassigned-count">
            {rows.filter(r => !r.assigned).length} unassigned
          </span>
        )}
        {outlierCount > 0 && (
          <span className="outlier-count">
            {outlierCount} outlier{outlierCount > 1 ? 's' : ''} (|z| &gt; 2)
          </span>
        )}
      </div>
    </div>
  );
}

export default AssignmentsTable;
