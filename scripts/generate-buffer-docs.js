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
  <title>Buffers</title>
  <style>
    :root {
      --color-accent: #1a73e8;
      --color-text: #213547;
      --color-text-secondary: #666;
      --color-bg: #f5f5f5;
      --color-card: #fff;
      --color-border: #ddd;
      --color-code-bg: #f8f9fa;
    }

    * { box-sizing: border-box; }

    body {
      font-family: system-ui, -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
      line-height: 1.6;
      color: var(--color-text);
      background: var(--color-bg);
      margin: 0;
      padding: 1rem;
    }

    .container {
      max-width: 1100px;
      margin: 0 auto;
      background: var(--color-card);
      padding: 2rem;
      border-radius: 12px;
      box-shadow: 0 2px 8px rgba(0,0,0,0.1);
    }

    h1 {
      color: var(--color-accent);
      border-bottom: 2px solid var(--color-accent);
      padding-bottom: 0.5rem;
    }

    h2 {
      color: var(--color-text);
      margin-top: 2.5rem;
      padding-top: 1.5rem;
      border-top: 2px solid var(--color-border);
    }

    h3 {
      color: var(--color-text-secondary);
      margin-top: 1.5rem;
      font-size: 1.1rem;
    }

    a { color: var(--color-accent); text-decoration: none; }
    a:hover { text-decoration: underline; }

    code {
      background: var(--color-code-bg);
      padding: 0.2em 0.4em;
      border-radius: 4px;
      font-size: 0.9em;
      font-family: monospace;
    }

    table {
      width: 100%;
      border-collapse: collapse;
      margin: 1rem 0;
      font-size: 0.85rem;
      overflow-x: auto;
      display: block;
    }

    thead, tbody {
      display: table;
      width: 100%;
    }

    th, td {
      text-align: center;
      padding: 0.4rem 0.5rem;
      border: 1px solid var(--color-border);
      vertical-align: middle;
    }

    th {
      background: var(--color-code-bg);
      font-weight: 600;
    }

    th.group-header {
      background: #d0e8f2;
      font-weight: 700;
    }

    td:first-child, th:first-child {
      text-align: left;
    }

    td:nth-child(2), th:nth-child(2) {
      text-align: left;
    }

    tbody tr:nth-child(even) { background: #fafafa; }

    .coeff-cell {
      font-size: 0.8em;
      color: var(--color-text-secondary);
    }

    .empty-cell {
      color: var(--color-text-secondary);
    }

    .nav {
      margin-bottom: 2rem;
      padding: 1rem;
      background: var(--color-code-bg);
      border-radius: 8px;
    }

    .nav a {
      margin-right: 1.5rem;
      font-weight: 500;
    }

    .meta {
      background: var(--color-code-bg);
      padding: 1rem;
      border-radius: 8px;
      margin-bottom: 1.5rem;
      font-size: 0.9rem;
    }

    .meta strong { color: var(--color-text); }

    .nuclei-badges {
      display: inline-flex;
      gap: 0.4rem;
      margin-left: 0.75rem;
    }

    .nucleus-badge {
      background: var(--color-accent);
      color: white;
      padding: 0.2rem 0.4rem;
      border-radius: 3px;
      font-size: 0.75rem;
      font-weight: 500;
    }

    .buffer-meta {
      background: var(--color-code-bg);
      padding: 0.75rem 1rem;
      border-radius: 6px;
      margin: 1rem 0;
      font-size: 0.9rem;
    }

    .pka-section {
      background: #e8f5e9;
      padding: 1rem;
      border-radius: 6px;
      margin: 1rem 0;
    }

    .pka-item {
      margin-bottom: 0.75rem;
      line-height: 1.8;
    }

    .pka-item:last-child {
      margin-bottom: 0;
    }

    .pka-line {
      display: block;
      padding-left: 1rem;
      color: var(--color-text-secondary);
    }

    .timestamp {
      text-align: right;
      color: var(--color-text-secondary);
      font-size: 0.8rem;
      margin-top: 2rem;
    }

    sup { font-size: 0.7em; }
    sub { font-size: 0.7em; }

    .coeff {
      display: block;
      font-size: 0.85em;
      color: var(--color-text-secondary);
      margin-top: 0.15rem;
    }
  </style>
</head>
<body>
  <div class="container">
    <nav class="nav">
      <a href="../">Home</a>
      <a href="./">Docs</a>
      <a href="./buffers.html">Buffers</a>
      <a href="https://github.com/waudbygroup/nmr-pH">GitHub</a>
    </nav>

    <h1>Buffer Database</h1>

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

  for (const nucleus of nuclei) {
    html += `<span class="nucleus-badge">${formatNucleus(nucleus)}</span>`;
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

    <hr style="margin-top: 2rem;">
    <p style="color: var(--color-text-secondary); font-size: 0.875rem; text-align: center;">
      NMR pH calibration · <a href="https://waudbylab.org">Waudby Group</a> · UCL School of Pharmacy
    </p>
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
