import { useMemo } from 'react';
import Plot from 'react-plotly.js';
import { generateShiftCurves } from '../numerical/bufferModel';

/**
 * Color palette for buffers.
 */
const BUFFER_COLORS = [
  '#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd',
  '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf'
];

/**
 * Get color for a buffer based on index.
 */
function getBufferColor(index) {
  return BUFFER_COLORS[index % BUFFER_COLORS.length];
}

/**
 * Format nucleus label for axis.
 */
function formatNucleusLabel(nucleus) {
  const mass = nucleus.match(/^\d+/)?.[0] || '';
  const element = nucleus.replace(/^\d+/, '');
  return `<sup>${mass}</sup>${element}`;
}

/**
 * ChemicalShiftPlot component.
 * Displays chemical shift vs pH curves for selected buffers.
 */
export function ChemicalShiftPlot({
  nucleus,
  buffers,
  samplesMap,
  temperature,
  ionicStrength,
  observedShifts = [],
  fittedPH = null,
  phUncertainty = null,
  assignments = null,
  pHRange = [2, 12],
  height = 600 // 50% taller than original 400
}) {
  // Generate curve data for all selected buffers
  const curveData = useMemo(() => {
    const allCurves = [];

    buffers.forEach((buffer, bufferIndex) => {
      const sample = samplesMap.get(buffer.sample_id);
      const curves = generateShiftCurves(
        buffer,
        nucleus,
        temperature,
        ionicStrength,
        sample,
        pHRange[0],
        pHRange[1],
        0.05
      );

      curves.forEach(curve => {
        allCurves.push({
          ...curve,
          color: getBufferColor(bufferIndex)
        });
      });
    });

    return allCurves;
  }, [buffers, samplesMap, nucleus, temperature, ionicStrength, pHRange]);

  // Build Plotly traces
  const traces = useMemo(() => {
    const plotTraces = [];

    // Add buffer curves
    curveData.forEach((curve, i) => {
      plotTraces.push({
        x: curve.shifts,
        y: curve.pHValues,
        type: 'scatter',
        mode: 'lines',
        name: `${curve.buffer_name} ${curve.resonance_id}`,
        line: {
          color: curve.color,
          width: 2
        },
        hovertemplate: `${curve.buffer_name}<br>${curve.resonance_id}<br>δ: %{x:.3f} ppm<br>pH: %{y:.2f}<extra></extra>`
      });
    });

    // Add observed shift vertical lines
    if (observedShifts.length > 0) {
      observedShifts.forEach((shift, i) => {
        plotTraces.push({
          x: [shift, shift],
          y: pHRange,
          type: 'scatter',
          mode: 'lines',
          name: i === 0 ? 'Observed' : undefined,
          showlegend: i === 0,
          line: {
            color: 'rgba(0, 0, 0, 0.5)',
            width: 1,
            dash: 'dash'
          },
          hovertemplate: `Observed: ${shift.toFixed(3)} ppm<extra></extra>`
        });
      });
    }

    // Add fitted pH horizontal line
    if (fittedPH !== null) {
      // Uncertainty band
      if (phUncertainty !== null && phUncertainty > 0) {
        const shiftMin = Math.min(...curveData.flatMap(c => c.shifts));
        const shiftMax = Math.max(...curveData.flatMap(c => c.shifts));

        plotTraces.push({
          x: [shiftMin, shiftMax, shiftMax, shiftMin],
          y: [
            fittedPH - phUncertainty,
            fittedPH - phUncertainty,
            fittedPH + phUncertainty,
            fittedPH + phUncertainty
          ],
          type: 'scatter',
          fill: 'toself',
          fillcolor: 'rgba(255, 0, 0, 0.1)',
          line: { color: 'transparent' },
          name: 'pH uncertainty',
          showlegend: false,
          hoverinfo: 'skip'
        });
      }

      // Fitted pH line
      const shiftMin = curveData.length > 0
        ? Math.min(...curveData.flatMap(c => c.shifts))
        : 0;
      const shiftMax = curveData.length > 0
        ? Math.max(...curveData.flatMap(c => c.shifts))
        : 10;

      plotTraces.push({
        x: [shiftMin, shiftMax],
        y: [fittedPH, fittedPH],
        type: 'scatter',
        mode: 'lines',
        name: `Fitted pH: ${fittedPH.toFixed(2)}`,
        line: {
          color: 'red',
          width: 2
        },
        hovertemplate: `Fitted pH: ${fittedPH.toFixed(2)}<extra></extra>`
      });
    }

    // Add assignment markers
    if (assignments && assignments.length > 0) {
      const assignedShifts = [];
      const assignedPHs = [];
      const markerTexts = [];

      for (const assignment of assignments) {
        if (assignment.assigned) {
          assignedShifts.push(assignment.observed_shift);
          assignedPHs.push(fittedPH ?? 7);
          markerTexts.push(`${assignment.buffer_name}<br>${assignment.resonance_id}`);
        }
      }

      if (assignedShifts.length > 0) {
        plotTraces.push({
          x: assignedShifts,
          y: assignedPHs,
          type: 'scatter',
          mode: 'markers',
          name: 'Assigned peaks',
          marker: {
            color: 'red',
            size: 10,
            symbol: 'circle'
          },
          text: markerTexts,
          hovertemplate: '%{text}<br>δ: %{x:.3f} ppm<extra></extra>'
        });
      }
    }

    return plotTraces;
  }, [curveData, observedShifts, fittedPH, phUncertainty, assignments, pHRange]);

  // Layout configuration with legend below and axis labels
  const layout = useMemo(() => ({
    xaxis: {
      title: {
        text: `${formatNucleusLabel(nucleus)} chemical shift (ppm)`,
        font: { size: 14 }
      },
      autorange: 'reversed', // NMR convention: high field to low field
      showgrid: true,
      gridcolor: '#eee'
    },
    yaxis: {
      title: {
        text: 'pH',
        font: { size: 14 }
      },
      range: pHRange,
      showgrid: true,
      gridcolor: '#eee'
    },
    legend: {
      orientation: 'h',
      yanchor: 'top',
      y: -0.15,
      xanchor: 'center',
      x: 0.5
    },
    margin: { t: 20, r: 20, b: 120, l: 60 },
    hovermode: 'closest',
    showlegend: true
  }), [nucleus, pHRange]);

  const config = {
    responsive: true,
    displayModeBar: true,
    modeBarButtonsToRemove: ['lasso2d', 'select2d'],
    toImageButtonOptions: {
      format: 'svg',
      filename: `${nucleus}_shift_vs_pH`
    }
  };

  if (buffers.length === 0) {
    return (
      <div className="chemical-shift-plot empty" style={{ height }}>
        <p>Select buffers to see chemical shift curves</p>
      </div>
    );
  }

  return (
    <div className="chemical-shift-plot">
      <Plot
        data={traces}
        layout={layout}
        config={config}
        style={{ width: '100%', height }}
        useResizeHandler
      />
    </div>
  );
}

export default ChemicalShiftPlot;
