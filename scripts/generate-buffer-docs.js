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
  <title>Buffer Database - NMR Buffer pH Estimation</title>
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
      margin-top: 2rem;
      padding-top: 1rem;
      border-top: 1px solid var(--color-border);
    }

    h3 { color: var(--color-accent); margin-top: 1.5rem; }
    h4 { color: var(--color-text-secondary); margin-top: 1rem; }

    a { color: var(--color-accent); }

    code {
      background: var(--color-code-bg);
      padding: 0.2em 0.4em;
      border-radius: 4px;
      font-size: 0.9em;
    }

    table {
      width: 100%;
      border-collapse: collapse;
      margin: 1rem 0;
      font-size: 0.9rem;
    }

    th, td {
      text-align: left;
      padding: 0.5rem 0.75rem;
      border: 1px solid var(--color-border);
    }

    th {
      background: var(--color-code-bg);
      font-weight: 600;
    }

    tr:nth-child(even) { background: #fafafa; }

    .nav {
      margin-bottom: 2rem;
      padding: 1rem;
      background: var(--color-code-bg);
      border-radius: 8px;
    }

    .nav a {
      margin-right: 1.5rem;
      text-decoration: none;
      font-weight: 500;
    }

    .nav a:hover { text-decoration: underline; }

    .meta {
      background: var(--color-code-bg);
      padding: 1rem;
      border-radius: 8px;
      margin-bottom: 1.5rem;
      font-size: 0.9rem;
    }

    .meta strong { color: var(--color-text); }

    .buffer-card {
      background: var(--color-card);
      border: 1px solid var(--color-border);
      border-radius: 8px;
      padding: 1.5rem;
      margin-bottom: 1.5rem;
    }

    .buffer-card h3 {
      margin-top: 0;
      border-bottom: 1px solid var(--color-border);
      padding-bottom: 0.5rem;
    }

    .nuclei-badges {
      display: flex;
      gap: 0.5rem;
      margin-bottom: 1rem;
    }

    .nucleus-badge {
      background: var(--color-accent);
      color: white;
      padding: 0.25rem 0.5rem;
      border-radius: 4px;
      font-size: 0.8rem;
      font-weight: 500;
    }

    .pka-info {
      background: #e8f5e9;
      padding: 0.75rem 1rem;
      border-radius: 6px;
      margin-bottom: 1rem;
    }

    .timestamp {
      text-align: right;
      color: var(--color-text-secondary);
      font-size: 0.8rem;
      margin-top: 2rem;
    }

    sup { font-size: 0.7em; }
  </style>
</head>
<body>
  <div class="container">
    <nav class="nav">
      <a href="../">App</a>
      <a href="./">Documentation</a>
      <a href="./buffers.html">Buffer Database</a>
      <a href="https://github.com/waudbygroup/nmr-pH">GitHub</a>
    </nav>

    <h1>Buffer Database</h1>

    <div class="meta">
      <strong>Database Version:</strong> ${db.database_version} ·
      <strong>Last Updated:</strong> ${db.last_updated} ·
      <strong>Buffers:</strong> ${db.buffers.length} ·
      <strong>Samples:</strong> ${db.samples.length}
    </div>

    <h2>Samples</h2>
    <table>
      <tr>
        <th>Sample ID</th>
        <th>Solvent</th>
        <th>Authors</th>
        <th>T<sub>ref</sub> (K)</th>
        <th>I<sub>ref</sub> (M)</th>
        <th>pH Range</th>
      </tr>
`;

// Add samples
for (const sample of db.samples) {
  const authors = sample.authors.map(a => a.name).join(', ');
  const phRange = sample.measurement_ranges?.pH
    ? `${sample.measurement_ranges.pH.min} - ${sample.measurement_ranges.pH.max}`
    : '-';

  html += `      <tr>
        <td><code>${sample.sample_id}</code></td>
        <td>${formatSolvent(sample.solvent)}</td>
        <td>${authors}</td>
        <td>${sample.reference_temperature_K}</td>
        <td>${sample.reference_ionic_strength_M}</td>
        <td>${phRange}</td>
      </tr>
`;
}

html += `    </table>

    <h2>Buffers</h2>
`;

// Add buffers
for (const buffer of db.buffers) {
  const sample = samplesMap.get(buffer.sample_id);
  const nuclei = Object.keys(buffer.chemical_shifts);

  html += `
    <div class="buffer-card">
      <h3>${buffer.buffer_name}</h3>
      <div class="nuclei-badges">
`;

  for (const nucleus of nuclei) {
    html += `        <span class="nucleus-badge">${formatNucleus(nucleus)}</span>\n`;
  }

  html += `      </div>

      <p>
        <strong>Buffer ID:</strong> <code>${buffer.buffer_id}</code><br>
        <strong>Family:</strong> ${buffer.buffer_family}<br>
        <strong>Solvent:</strong> ${sample ? formatSolvent(sample.solvent) : 'Unknown'}<br>
        <strong>Ionisation States:</strong> ${buffer.ionisation_states}
      </p>

      <div class="pka-info">
`;

  // pKa parameters
  for (const pka of buffer.pKa_parameters) {
    html += `        <strong>pKa${pka.pKa_index}:</strong> ${formatValue(pka.pKa)}`;
    if (pka.dH_kJ_mol) {
      html += ` · ΔH = ${formatValue(pka.dH_kJ_mol, 1)} kJ/mol`;
    }
    if (pka.dCp_J_mol_K) {
      html += ` · ΔCp = ${formatValue(pka.dCp_J_mol_K, 0)} J/(mol·K)`;
    }
    html += `<br>\n`;
  }

  html += `      </div>

      <h4>Chemical Shifts</h4>
`;

  // Chemical shifts table for each nucleus
  for (const [nucleus, resonances] of Object.entries(buffer.chemical_shifts)) {
    html += `
      <p><strong>${formatNucleus(nucleus)}</strong></p>
      <table>
        <tr>
          <th>Resonance</th>
          <th>Description</th>
`;

    // Add column headers for each ionisation state
    for (let i = 0; i < buffer.ionisation_states; i++) {
      html += `          <th>State ${i} (ppm)</th>\n`;
    }

    html += `        </tr>
`;

    for (const res of resonances) {
      html += `        <tr>
          <td><code>${res.resonance_id}</code></td>
          <td>${res.description || '-'}</td>
`;

      // Get shifts for each state
      for (let i = 0; i < buffer.ionisation_states; i++) {
        const stateShift = res.limiting_shifts.find(ls => ls.ionisation_state === i);
        if (stateShift) {
          html += `          <td>${formatValue(stateShift.shift_ppm, 3)}</td>\n`;
        } else {
          html += `          <td>-</td>\n`;
        }
      }

      html += `        </tr>
`;
    }

    html += `      </table>
`;
  }

  if (buffer.notes) {
    html += `      <p><em>${buffer.notes}</em></p>\n`;
  }

  html += `    </div>
`;
}

html += `
    <div class="timestamp">
      Generated: ${new Date().toISOString().split('T')[0]}
    </div>

    <hr style="margin-top: 2rem;">
    <p style="color: var(--color-text-secondary); font-size: 0.875rem; text-align: center;">
      NMR Buffer pH Estimation · <a href="https://waudbylab.org">Waudby Group</a> · UCL School of Pharmacy
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
