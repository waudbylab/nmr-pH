#!/usr/bin/env node
/**
 * Generate buffer documentation HTML from the database JSON.
 *
 * Usage: node scripts/generate-buffer-docs.js
 *
 * Reads: public/database/current/database.json
 * Writes: public/docs/buffers.html
 */

import { readFileSync, writeFileSync, mkdirSync } from 'fs';
import { dirname, join } from 'path';
import { fileURLToPath } from 'url';

const __dirname = dirname(fileURLToPath(import.meta.url));
const rootDir = join(__dirname, '..');

// Read database
const dbPath = join(rootDir, 'public/database/current/database.json');
const db = JSON.parse(readFileSync(dbPath, 'utf-8'));

// Helper to format value with uncertainty
function formatValue(val, decimals = 2) {
  if (Array.isArray(val)) {
    return `${val[0].toFixed(decimals)} ± ${val[1].toFixed(decimals)}`;
  }
  return typeof val === 'number' ? val.toFixed(decimals) : val;
}

// Helper to format solvent name
function formatSolvent(solvent) {
  const names = {
    '10pct_D2O': '10% D₂O / 90% H₂O',
    '100pct_D2O': '100% D₂O',
    'H2O': 'H₂O',
    'other': 'Other'
  };
  return names[solvent] || solvent;
}

// Helper to format nucleus with superscript
function formatNucleus(nucleus) {
  const match = nucleus.match(/^(\d+)(\w+)$/);
  if (match) {
    return `<sup>${match[1]}</sup>${match[2]}`;
  }
  return nucleus;
}

// Build samples map
const samplesMap = new Map(db.samples.map(s => [s.sample_id, s]));

// Generate HTML
let html = `<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>Buffer Database - NMR pH calibration</title>
  <link rel="stylesheet" href="./shared.css">
</head>
<body>
  <div class="page-container">
    <nav class="app-nav">
      <a href="../">Home</a>
      <a href="./">Docs</a>
      <a href="./buffers.html">Buffers</a>
      <a href="https://github.com/waudbylab/nmr-pH">GitHub</a>
    </nav>

    <header class="app-header">
      <h1>Buffer Database</h1>
      <p class="subtitle">Comprehensive NMR buffer chemical shift database</p>
      <span class="version">v${db.database_version}</span>
    </header>

    <section>
    <div class="meta">
      <strong>Database Version:</strong> ${db.database_version} ·
      <strong>Last Updated:</strong> ${db.last_updated} ·
      <strong>Buffers:</strong> ${db.buffers.length} ·
      <strong>Samples:</strong> ${db.samples.length}
    </div>
`;

// Add buffers as H2 sections
for (const buffer of db.buffers) {
  const sample = samplesMap.get(buffer.sample_id);
  const nuclei = Object.keys(buffer.chemical_shifts);

  html += `    <h2>${buffer.buffer_name}`;

  if (nuclei.length > 0) {
    html += `<span class="nuclei-badges">`;
    for (const nucleus of nuclei) {
      html += `<span class="nucleus-badge">${formatNucleus(nucleus)}</span>`;
    }
    html += `</span>`;
  }

  html += `</h2>

    <div class="buffer-meta">
      <strong>Buffer ID:</strong> <code>${buffer.buffer_id}</code><br>
      <strong>Family:</strong> ${buffer.buffer_family}<br>
      <strong>Solvent:</strong> ${sample ? formatSolvent(sample.solvent) : 'Unknown'}<br>
      <strong>Ionisation States:</strong> ${buffer.ionisation_states}<br>
      <strong>Sample:</strong> <code>${buffer.sample_id}</code> (T<sub>ref</sub> = ${sample?.reference_temperature_K || 298.15} K, I<sub>ref</sub> = ${sample?.reference_ionic_strength_M || 0} M)
    </div>

    <h3>Thermodynamic Parameters</h3>
    <div class="pka-section">
`;

  // pKa parameters - each on new line
  for (const pka of buffer.pKa_parameters) {
    html += `      <div class="pka-item">
        <strong>pKa${pka.pKa_index}:</strong> ${formatValue(pka.pKa)}
`;
    if (pka.dH_kJ_mol) {
      html += `        <span class="pka-line">ΔH = ${formatValue(pka.dH_kJ_mol, 1)} kJ/mol</span>
`;
    }
    if (pka.dCp_J_mol_K) {
      html += `        <span class="pka-line">ΔCp = ${formatValue(pka.dCp_J_mol_K, 0)} J/(mol·K)</span>
`;
    }
    if (pka.ion_size_angstrom) {
      html += `        <span class="pka-line">Ion size = ${pka.ion_size_angstrom} Å (for ionic strength corrections)</span>
`;
    }
    html += `        <span class="pka-line">Protonated charge: ${pka.protonated_charge > 0 ? '+' : ''}${pka.protonated_charge}</span>
      </div>
`;
  }

  html += `    </div>

    <h3>Chemical Shifts</h3>
`;

  // Chemical shifts table for each nucleus
  for (const [nucleus, resonances] of Object.entries(buffer.chemical_shifts)) {
    html += `    <p><strong>${formatNucleus(nucleus)}</strong></p>
    <div class="table-wrapper">
      <table>
        <thead>
          <tr>
            <th rowspan="2">Resonance</th>
            <th rowspan="2">Description</th>
`;

    // Add group headers for each ionisation state
    for (let i = 0; i < buffer.ionisation_states; i++) {
      html += `          <th colspan="3" class="group-header">State ${i}</th>
`;
    }

    html += `        </tr>
        <tr>
`;

    // Add sub-headers (δ, αT, αI) for each state
    for (let i = 0; i < buffer.ionisation_states; i++) {
      html += `          <th>δ (ppm)</th>
          <th>α<sub>T</sub> (ppm/K)</th>
          <th>α<sub>I</sub> (ppm/M)</th>
`;
    }

    html += `        </tr>
      </thead>
      <tbody>
`;

    for (const res of resonances) {
      html += `        <tr>
          <td><code>${res.resonance_id}</code></td>
          <td>${res.description || '-'}</td>
`;

      // Get shifts for each state with coefficients in separate columns
      for (let i = 0; i < buffer.ionisation_states; i++) {
        const stateShift = res.limiting_shifts.find(ls => ls.ionisation_state === i);
        if (stateShift) {
          // Chemical shift
          html += `          <td>${formatValue(stateShift.shift_ppm, 3)}</td>
`;
          // Temperature coefficient
          if (stateShift.temperature_coefficient_ppm_per_K) {
            html += `          <td class="coeff-cell">${formatValue(stateShift.temperature_coefficient_ppm_per_K, 4)}</td>
`;
          } else {
            html += `          <td class="empty-cell">-</td>
`;
          }
          // Ionic strength coefficient
          if (stateShift.ionic_strength_coefficient_ppm_per_M) {
            html += `          <td class="coeff-cell">${formatValue(stateShift.ionic_strength_coefficient_ppm_per_M, 2)}</td>
`;
          } else {
            html += `          <td class="empty-cell">-</td>
`;
          }
        } else {
          html += `          <td class="empty-cell">-</td>
          <td class="empty-cell">-</td>
          <td class="empty-cell">-</td>
`;
        }
      }

      html += `        </tr>
`;
    }

    html += `      </tbody>
      </table>
    </div>
`;
  }

  if (buffer.notes) {
    html += `    <p style="font-style: italic; color: var(--color-text-secondary);"><strong>Notes:</strong> ${buffer.notes}</p>
`;
  }
}

html += `    <div class="timestamp">
      Generated: ${new Date().toISOString().split('T')[0]}
    </div>
    </section>

    <footer class="app-footer">
      <p class="footer-affiliation">
        NMR pH calibration · <a href="https://waudbylab.org">Waudby Group</a> · UCL School of Pharmacy
      </p>
    </footer>
  </div>
</body>
</html>
`;

// Ensure docs directory exists
mkdirSync(join(rootDir, 'public/docs'), { recursive: true });

// Write output
const outputPath = join(rootDir, 'public/docs/buffers.html');
writeFileSync(outputPath, html);

console.log(`Generated: ${outputPath}`);
console.log(`  - ${db.samples.length} samples`);
console.log(`  - ${db.buffers.length} buffers`);
