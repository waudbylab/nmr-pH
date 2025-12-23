import { jsPDF } from 'jspdf';
import html2canvas from 'html2canvas';

/**
 * Generate JSON data for download.
 */
function generateJSON(result, conditions, buffers, observedShifts, referencingInfo) {
  return JSON.stringify({
    timestamp: new Date().toISOString(),
    input: {
      conditions,
      observedShifts,
      referencing: referencingInfo,
      buffers: buffers.map(b => ({
        buffer_id: b.buffer_id,
        buffer_name: b.buffer_name
      }))
    },
    output: result.success ? {
      parameters: result.parameters,
      conditions: result.conditions,
      statistics: result.statistics,
      assignments: result.assignments
    } : {
      error: result.error
    }
  }, null, 2);
}

/**
 * Format nucleus for display.
 */
function formatNucleus(nucleus) {
  const mass = nucleus.match(/^\d+/)?.[0] || '';
  const element = nucleus.replace(/^\d+/, '');
  return `${mass}${element}`;
}

/**
 * Generate PDF report with plots.
 */
async function generatePDF(result, conditions, buffers, samplesMap, nuclei, referencingInfo) {
  const doc = new jsPDF();
  let y = 20;

  // Title
  doc.setFontSize(18);
  doc.text('NMR Buffer pH Estimation Report', 20, y);
  y += 15;

  // Timestamp
  doc.setFontSize(10);
  doc.text(`Generated: ${new Date().toLocaleString()}`, 20, y);
  y += 15;

  // Results
  if (result.success) {
    doc.setFontSize(14);
    doc.text('Fitted Parameters', 20, y);
    y += 10;

    doc.setFontSize(11);
    const params = result.parameters;

    // pH
    const phUnc = params.pH.uncertainty;
    doc.text(`pH: ${params.pH.value.toFixed(3)} ${phUnc ? `± ${phUnc.toFixed(3)}` : ''}`, 25, y);
    y += 7;

    // Temperature
    if (params.temperature) {
      const tUnc = params.temperature.uncertainty;
      doc.text(`Temperature: ${params.temperature.value.toFixed(2)} ${tUnc ? `± ${tUnc.toFixed(2)}` : ''} K`, 25, y);
    } else {
      doc.text(`Temperature: ${conditions.temperature.toFixed(2)} K (fixed)`, 25, y);
    }
    y += 7;

    // Ionic strength
    if (params.ionicStrength) {
      const iUnc = params.ionicStrength.uncertainty;
      doc.text(`Ionic Strength: ${params.ionicStrength.value.toFixed(4)} ${iUnc ? `± ${iUnc.toFixed(4)}` : ''} M`, 25, y);
    } else {
      doc.text(`Ionic Strength: ${conditions.ionicStrength.toFixed(4)} M (fixed)`, 25, y);
    }
    y += 12;

    // Reference Configuration
    if (nuclei && nuclei.length > 0 && referencingInfo) {
      doc.setFontSize(14);
      doc.text('Reference Configuration', 20, y);
      y += 10;

      doc.setFontSize(11);
      for (const nucleus of nuclei) {
        let status = '';
        if (nucleus === '1H') {
          if (referencingInfo.hasDSS) {
            status = `Referenced to DSS at ${referencingInfo.dssShift ?? 0} ppm`;
          } else {
            const offset = referencingInfo.fittedReferenceOffsets?.['1H'];
            status = offset !== undefined
              ? `Fitted: ${offset.toFixed(3)} ppm`
              : 'Fitting from water signal';
          }
        } else {
          if (referencingInfo.hasDSS && referencingInfo.heteroReferencedToDSS) {
            status = 'Referenced to DSS at 0 ppm';
          } else if (referencingInfo.spectrometerFreqs?.['1H'] && referencingInfo.spectrometerFreqs?.[nucleus]) {
            const offset = referencingInfo.fittedReferenceOffsets?.[nucleus];
            status = offset !== undefined
              ? `Linked to ¹H: ${offset.toFixed(3)} ppm`
              : 'Linked to ¹H';
          } else {
            const offset = referencingInfo.fittedReferenceOffsets?.[nucleus];
            status = offset !== undefined
              ? `Fitted: ${offset.toFixed(3)} ppm`
              : 'Fitting independently';
          }
        }
        doc.text(`${formatNucleus(nucleus)}: ${status}`, 25, y);
        y += 7;
      }
      y += 5;
    }

    // Fit statistics
    doc.setFontSize(14);
    doc.text('Fit Statistics', 20, y);
    y += 10;

    doc.setFontSize(11);
    const stats = result.statistics;
    doc.text(`RMSD: ${stats.rmsd.toFixed(4)} ppm`, 25, y);
    y += 7;
    doc.text(`Reduced chi-squared: ${stats.reducedChiSquared.toFixed(3)}`, 25, y);
    y += 7;
    doc.text(`Observations: ${stats.nObservations}`, 25, y);
    y += 7;
    doc.text(`Parameters: ${stats.nParameters}`, 25, y);
    y += 7;
    doc.text(`Degrees of freedom: ${stats.degreesOfFreedom}`, 25, y);
    y += 12;

    // Peak assignments
    doc.setFontSize(14);
    doc.text('Peak Assignments', 20, y);
    y += 10;

    doc.setFontSize(10);
    for (const [nucleus, assignments] of Object.entries(result.assignments)) {
      if (!assignments || assignments.length === 0) continue;

      doc.setFontSize(11);
      doc.text(`${formatNucleus(nucleus)}:`, 25, y);
      y += 6;

      doc.setFontSize(9);
      for (const assignment of assignments) {
        if (assignment.assigned) {
          const residual = assignment.residual !== undefined
            ? ` (residual: ${assignment.residual.toFixed(4)} ppm)`
            : '';
          doc.text(
            `  ${assignment.observed_shift.toFixed(3)} ppm → ${assignment.buffer_name} ${assignment.resonance_id}${residual}`,
            25,
            y
          );
        } else {
          doc.text(`  ${assignment.observed_shift.toFixed(3)} ppm → unassigned`, 25, y);
        }
        y += 5;

        // Check for page break
        if (y > 270) {
          doc.addPage();
          y = 20;
        }
      }
      y += 3;
    }

    // Buffers used
    if (y > 240) {
      doc.addPage();
      y = 20;
    }

    doc.setFontSize(14);
    doc.text('Buffers Used', 20, y);
    y += 10;

    doc.setFontSize(11);
    for (const buffer of buffers) {
      doc.text(`- ${buffer.buffer_name} (${buffer.buffer_id})`, 25, y);
      y += 7;
    }

    // Try to capture plots
    try {
      const plotElements = document.querySelectorAll('.chemical-shift-plot .js-plotly-plot');
      if (plotElements.length > 0) {
        doc.addPage();
        y = 20;

        doc.setFontSize(14);
        doc.text('Chemical Shift Plots', 20, y);
        y += 15;

        for (const plotEl of plotElements) {
          if (y > 200) {
            doc.addPage();
            y = 20;
          }

          try {
            const canvas = await html2canvas(plotEl, {
              scale: 2,
              useCORS: true,
              logging: false
            });
            const imgData = canvas.toDataURL('image/png');

            // Calculate dimensions to fit on page
            const imgWidth = 170;
            const imgHeight = (canvas.height / canvas.width) * imgWidth;

            doc.addImage(imgData, 'PNG', 20, y, imgWidth, Math.min(imgHeight, 120));
            y += Math.min(imgHeight, 120) + 15;
          } catch (err) {
            console.warn('Could not capture plot:', err);
          }
        }
      }
    } catch (err) {
      console.warn('Could not capture plots:', err);
    }

  } else {
    doc.setFontSize(12);
    doc.text('Fitting failed:', 20, y);
    y += 10;
    doc.text(result.error, 25, y);
  }

  return doc;
}

/**
 * DownloadButtons component.
 * Buttons to download results in various formats.
 */
export function DownloadButtons({
  result,
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

  const referencingInfo = {
    hasDSS,
    dssShift,
    heteroReferencedToDSS,
    spectrometerFreqs,
    fittedReferenceOffsets
  };

  const downloadJSON = () => {
    const json = generateJSON(result, conditions, buffers, observedShifts, referencingInfo);
    const blob = new Blob([json], { type: 'application/json' });
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = 'nmr-ph-results.json';
    a.click();
    URL.revokeObjectURL(url);
  };

  const downloadPDF = async () => {
    try {
      const doc = await generatePDF(result, conditions, buffers, samplesMap, nuclei, referencingInfo);
      doc.save('nmr-ph-report.pdf');
    } catch (err) {
      console.error('PDF generation failed:', err);
      alert('Failed to generate PDF. Please try again.');
    }
  };

  const downloadCitations = () => {
    const citations = buffers.map(buffer => {
      const sample = samplesMap.get(buffer.sample_id);
      const authors = sample?.authors?.map(a => a.name).join(', ') || 'Unknown';
      const date = sample?.date_measured ? new Date(sample.date_measured).getFullYear() : 'n.d.';
      return `${authors} (${date}). ${buffer.buffer_name} buffer parameters. ${buffer.buffer_id}`;
    }).join('\n\n');

    const blob = new Blob([citations], { type: 'text/plain' });
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = 'nmr-buffer-citations.txt';
    a.click();
    URL.revokeObjectURL(url);
  };

  return (
    <div className="download-buttons">
      <h3>Download Results</h3>
      <div className="button-group">
        <button onClick={downloadPDF} type="button" className="download-button">
          PDF Report
        </button>
        <button onClick={downloadJSON} type="button" className="download-button">
          JSON Data
        </button>
        <button onClick={downloadCitations} type="button" className="download-button">
          Citations
        </button>
      </div>
    </div>
  );
}

export default DownloadButtons;
