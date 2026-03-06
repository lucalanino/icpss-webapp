'use strict';

// ─── Static model parameters ─────────────────────────────────
// Feature order (28 features):
// [age, wbc, hb, plt, blasts,
//  sex1, sex0,
//  ASXL11, ASXL10, DNMT3A1, DNMT3A0, EZH21, EZH20,
//  RUNX11, RUNX10, SETBP11, SETBP10, STAG21, STAG20,
//  TET21, TET20, TP531, TP530, U2AF11, U2AF10,
//  karyo2, karyo1, karyo0]

const MEANS = [
  72.4402606806662, 17.024112961622,  11.1337436640116, 153.790731354091, 5.16364952932657,
  0.685734974656046, 0.314265025343954,
  0.370745836350471, 0.629254163649529, 0.0637219406227371, 0.936278059377263,
  0.0716871832005793, 0.928312816799421,
  0.166545981173063, 0.833454018826937, 0.0666183924692252, 0.933381607530775,
  0.0217233888486604, 0.97827661115134,
  0.598117306299783, 0.401882693700217, 0.0304127443881245, 0.969587255611875,
  0.054308472121651, 0.945691527878349,
  0.110065170166546, 0.0977552498189718, 0.792179580014482
];

const SDS = [
  9.65440207119458, 21.0172656255669, 2.31241255929585, 146.136278767461, 4.50193317625419,
  0.46439065498275, 0.46439065498275,
  0.483179484601256, 0.483179484601256, 0.244345427303944, 0.244345427303944,
  0.258062694414854, 0.258062694414854,
  0.372704981369746, 0.372704981369746, 0.249450276612498, 0.249450276612498,
  0.14583169357556, 0.14583169357556,
  0.490456091341119, 0.490456091341119, 0.171782354579012, 0.171782354579012,
  0.226707473827046, 0.226707473827046,
  0.313084345415057, 0.297091018691704, 0.405894556401112
];

const COEFS = [
  0.260086480013197,  0.0890091509135308, -0.291935519712112, -0.156677579825297, 0.170601826757536,
   0.0307305831941834, -0.0307305833167902,
   0.0259690176532307, -0.0259690176532307, 0.0530808223296164, -0.0530808222209811,
   0.0662332046869963, -0.0662332047176384,
   0.0719405588474033, -0.0719405588474033, 0.0330306864448873, -0.0330306864448873,
   0.0738044597259139, -0.0738044597778723,
  -0.0692759279246606,  0.0692759279246606, 0.09098535645263, -0.0909853563679199,
   0.0390218781545528, -0.0390218781721499,
   0.100817861018277, -0.00609948408468122, -0.0733007910926501
];

// Indices of demographics features: age=0, sex1=5, sex0=6
const DEMO_INDICES = [0, 5, 6];

// Gene pairs: [_1 index, _0 index] in the 28-feature vector
// prior p_mutated = MEANS[i1] for each gene
const GENE_PAIRS = [
  { name: 'ASXL1',  i1: 7,  i0: 8  },
  { name: 'DNMT3A', i1: 9,  i0: 10 },
  { name: 'EZH2',   i1: 11, i0: 12 },
  { name: 'RUNX1',  i1: 13, i0: 14 },
  { name: 'SETBP1', i1: 15, i0: 16 },
  { name: 'STAG2',  i1: 17, i0: 18 },
  { name: 'TET2',   i1: 19, i0: 20 },
  { name: 'TP53',   i1: 21, i0: 22 },
  { name: 'U2AF1',  i1: 23, i0: 24 },
];

const THRESHOLDS = [-0.5935335, 0.1467796, 0.5660161, 1.2584983];
const CLASSES = ['VL', 'L', 'I', 'H', 'VH'];

// Class color palette
const CLASS_COLORS = {
  VL: { bg: '#2166ac', text: '#ffffff' },
  L:  { bg: '#67a9cf', text: '#ffffff' },
  I:  { bg: '#f5a623', text: '#ffffff' },
  H:  { bg: '#d6604d', text: '#ffffff' },
  VH: { bg: '#b2182b', text: '#ffffff' }
};

// Full risk class names
const CLASS_NAMES = {
  VL: 'Very Low',
  L:  'Low',
  I:  'Intermediate',
  H:  'High',
  VH: 'Very High'
};

// Outcome data by risk class index (0=VL … 4=VH)
// Source: [Citation in press]
const OUTCOMES = [
  { surv5y: '71%', surv5y_ci: '65–77%',  aml3y: '6.5%', aml3y_ci: '4.3–9.3%' },
  { surv5y: '50%', surv5y_ci: '46–55%',  aml3y: '12%',  aml3y_ci: '9.9–14%'  },
  { surv5y: '29%', surv5y_ci: '23–36%',  aml3y: '15%',  aml3y_ci: '12–19%'   },
  { surv5y: '15%', surv5y_ci: '9.6–23%', aml3y: '17%',  aml3y_ci: '14–21%'   },
  { surv5y: '3.2%',surv5y_ci: '0.6–18%', aml3y: '27%',  aml3y_ci: '21–34%'   },
];

// ─── Encoding ─────────────────────────────────────────────────────────────────

/**
 * Build the 28-element encoded feature vector from raw patient inputs.
 * Missing values are represented as NaN.
 *
 * @param {Object} raw - Raw patient inputs:
 *   sex (1=Male, 0=Female), age, wbc, hb, plt, blasts,
 *   karyo (0/1/2), ASXL1, DNMT3A, EZH2, RUNX1, SETBP1, STAG2, TET2, TP53, U2AF1
 * @returns {number[]} 28-element array
 */
function encodePatient(raw) {
  const toNum = (v) => (v === null || v === undefined || v === '' ? NaN : Number(v));
  const toBin = (v, target) => {
    const n = toNum(v);
    if (isNaN(n)) return NaN;
    return n === target ? 1 : 0;
  };

  const sex = toNum(raw.sex);
  const age = toNum(raw.age);
  const wbc = toNum(raw.wbc);
  const hb  = toNum(raw.hb);
  const plt = toNum(raw.plt);
  const bla = toNum(raw.blasts);
  const kar = toNum(raw.karyo);

  return [
    age, wbc, hb, plt, bla,
    // sex: 1=Male, 0=Female → one-hot
    isNaN(sex) ? NaN : (sex === 1 ? 1 : 0),
    isNaN(sex) ? NaN : (sex === 0 ? 1 : 0),
    // Binary mutations: one-hot [mut1, mut0]
    toBin(raw.ASXL1,  1), toBin(raw.ASXL1,  0),
    toBin(raw.DNMT3A, 1), toBin(raw.DNMT3A, 0),
    toBin(raw.EZH2,   1), toBin(raw.EZH2,   0),
    toBin(raw.RUNX1,  1), toBin(raw.RUNX1,  0),
    toBin(raw.SETBP1, 1), toBin(raw.SETBP1, 0),
    toBin(raw.STAG2,  1), toBin(raw.STAG2,  0),
    toBin(raw.TET2,   1), toBin(raw.TET2,   0),
    toBin(raw.TP53,   1), toBin(raw.TP53,   0),
    toBin(raw.U2AF1,  1), toBin(raw.U2AF1,  0),
    // karyo: three-way one-hot [2, 1, 0]
    isNaN(kar) ? NaN : (kar === 2 ? 1 : 0),
    isNaN(kar) ? NaN : (kar === 1 ? 1 : 0),
    isNaN(kar) ? NaN : (kar === 0 ? 1 : 0),
  ];
}

// ─── Score computation ────────────────────────────────────────────────────────

/**
 * Compute the iCPSS score and risk class for an encoded patient vector.
 *
 * @param {number[]} encoded - 28-element array from encodePatient()
 * @param {boolean} includeDemographics - whether to include age/sex
 * @returns {{ score: number, riskClass: string }}
 */
function computeScore(encoded, includeDemographics) {
  // Standardize: z = (value - mean) / sd; NaN → 0
  let z = encoded.map((v, i) => {
    if (isNaN(v)) return 0;
    return (v - MEANS[i]) / SDS[i];
  });

  let coefs = COEFS;

  if (!includeDemographics) {
    // Remove age (0), sex1 (5), sex0 (6)
    z     = z.filter((_, i)     => !DEMO_INDICES.includes(i));
    coefs = coefs.filter((_, i) => !DEMO_INDICES.includes(i));
  }

  // Dot product
  const score = z.reduce((sum, zi, i) => sum + zi * coefs[i], 0);

  // Assign risk class
  const riskClass = assignClass(score);

  return { score, riskClass };
}

/**
 * Assign a risk class based on score and thresholds.
 * @param {number} score
 * @returns {string}
 */
function assignClass(score) {
  // CLASSES = ['VL','L','I','H','VH'], THRESHOLDS has 4 cutpoints
  // score < T[0] → VL, T[0]≤score<T[1] → L, ..., score≥T[3] → VH
  let idx = THRESHOLDS.length; // default to last class
  for (let i = THRESHOLDS.length - 1; i >= 0; i--) {
    if (score >= THRESHOLDS[i]) {
      idx = i + 1;
      break;
    }
    idx = i;
  }
  return CLASSES[idx];
}

// Karyo one-hot indices in the 28-feature vector: [karyo2, karyo1, karyo0]
const KARYO_INDICES = [25, 26, 27];

// ─── Confidence estimation ────────────────────────────────────────────────────

/**
 * Compute a probabilistic confidence score for cases where gene mutation
 * status or karyotype is unknown. Enumerates all combinations (2 states per
 * missing binary gene, 3 states for missing ternary karyotype), weights each
 * by training-set prior probabilities, and returns the probability that the
 * true risk class equals the computed class (shift = 0).
 *
 * @param {number[]} encoded - 28-element array from encodePatient()
 * @param {string}   baseClass - the class computed with NaN→0 imputation
 * @param {boolean}  includeDemographics
 * @returns {{ value: number, n: number }|null} null if all features known
 */
function computeConfidence(encoded, baseClass, includeDemographics) {
  const missingGenes = GENE_PAIRS.filter(g => isNaN(encoded[g.i1]));
  const karyoMissing = isNaN(encoded[KARYO_INDICES[0]]);

  const nMissing = missingGenes.length + (karyoMissing ? 1 : 0);
  if (nMissing === 0) return null;

  const N       = missingGenes.length;
  const baseIdx = CLASSES.indexOf(baseClass);

  // Pre-compute active feature indices (handle demographics exclusion once)
  const activeIdxs  = COEFS.map((_, i) => i).filter(i => includeDemographics || !DEMO_INDICES.includes(i));
  const activeCoefs = activeIdxs.map(i => COEFS[i]);

  // Baseline z-array: NaN → 0 (same imputation as computeScore)
  const baseZ = encoded.map((v, i) => isNaN(v) ? 0 : (v - MEANS[i]) / SDS[i]);

  // Karyo states to enumerate: 3 if missing, just [null] (no-op) if known
  // State 0 → karyo=0 (Low):  karyo_0=1, karyo_1=0, karyo_2=0
  // State 1 → karyo=1 (Int):  karyo_0=0, karyo_1=1, karyo_2=0
  // State 2 → karyo=2 (High): karyo_0=0, karyo_1=0, karyo_2=1
  const karyoStates = karyoMissing ? [0, 1, 2] : [null];

  let pShift = 0;

  for (const ks of karyoStates) {
    let karyoWeight = 1;
    // Pre-compute karyo z-values and weight for this state
    const karyoZ = [null, null, null];
    if (ks !== null) {
      for (let j = 0; j < 3; j++) {
        const idx      = KARYO_INDICES[j];
        const karyoVal = 2 - j;              // idx 25→val 2, idx 26→val 1, idx 27→val 0
        const rawVal   = (ks === karyoVal) ? 1 : 0;
        karyoZ[j]      = (rawVal - MEANS[idx]) / SDS[idx];
      }
      // Prior probability: MEANS[25]=P(karyo=2), MEANS[26]=P(karyo=1), MEANS[27]=P(karyo=0)
      karyoWeight = MEANS[KARYO_INDICES[2 - ks]];
    }

    for (let mask = 0; mask < (1 << N); mask++) {
      const z    = baseZ.slice();
      let weight = karyoWeight;

      // Apply karyo z-values if missing
      if (ks !== null) {
        for (let j = 0; j < 3; j++) z[KARYO_INDICES[j]] = karyoZ[j];
      }

      // Apply gene z-values
      for (let k = 0; k < N; k++) {
        const g     = missingGenes[k];
        const isMut = (mask >> k) & 1;
        const p     = MEANS[g.i1];

        if (isMut) {
          z[g.i1] = (1 - MEANS[g.i1]) / SDS[g.i1];
          z[g.i0] = (0 - MEANS[g.i0]) / SDS[g.i0];
          weight  *= p;
        } else {
          z[g.i1] = (0 - MEANS[g.i1]) / SDS[g.i1];
          z[g.i0] = (1 - MEANS[g.i0]) / SDS[g.i0];
          weight  *= (1 - p);
        }
      }

      const score  = activeIdxs.reduce((sum, i, k) => sum + z[i] * activeCoefs[k], 0);
      const clsIdx = CLASSES.indexOf(assignClass(score));

      if (Math.abs(clsIdx - baseIdx) >= 1) {
        pShift += weight;
      }
    }
  }

  return { value: 1 - pShift, n: nMissing };
}

// ─── Batch: xlsx parsing ──────────────────────────────────────────────────────

const BATCH_COLUMNS = ['id', 'sex', 'age', 'wbc', 'hb', 'plt',
  'blasts', 'karyo', 'ASXL1', 'DNMT3A', 'EZH2', 'RUNX1',
  'SETBP1', 'STAG2', 'TET2', 'TP53', 'U2AF1'];

/**
 * Parse an xlsx ArrayBuffer into an array of patient row objects.
 * @param {ArrayBuffer} buffer
 * @returns {{ patients: Object[], errors: string[] }}
 */
function parseXlsx(buffer) {
  const workbook = XLSX.read(buffer, { type: 'array' });
  const sheetName = workbook.SheetNames[0];
  const sheet = workbook.Sheets[sheetName];
  const rows = XLSX.utils.sheet_to_json(sheet, { defval: null, raw: false });

  if (rows.length === 0) {
    return { patients: [], errors: ['The uploaded file is empty or has no data rows.'] };
  }

  // Normalize column names (trim whitespace)
  const normalized = rows.map(row => {
    const out = {};
    for (const k of Object.keys(row)) {
      out[k.trim()] = row[k];
    }
    return out;
  });

  // Check required columns (id is optional)
  const firstRow = normalized[0];
  const missingCols = BATCH_COLUMNS.slice(1).filter(c => !(c in firstRow));
  const errors = [];
  if (missingCols.length > 0) {
    errors.push(`Missing columns: ${missingCols.join(', ')}`);
    return { patients: [], errors };
  }

  return { patients: normalized, errors: [] };
}

// Required fields that must be present for a row to be scored (no imputation)
const BATCH_REQUIRED = ['wbc', 'hb', 'plt', 'blasts'];

// Shared range rules — must match validateSingleForm in form helpers
const FIELD_RULES = [
  { key: 'wbc',    label: 'WBC',            min: 0.1, max: 500,  isInt: false },
  { key: 'hb',     label: 'Hemoglobin',     min: 4,   max: 20,   isInt: false },
  { key: 'plt',    label: 'Platelets',      min: 1,   max: 1500, isInt: true  },
  { key: 'blasts', label: 'BM Blast Count', min: 0,   max: 19,   isInt: true  },
];

/**
 * Validate a batch row against the same rules as single-patient mode.
 * Returns an array of issue strings (empty = valid).
 */
function validateBatchRow(row, includeDemographics) {
  const issues = [];

  for (const f of FIELD_RULES) {
    const v = row[f.key];
    if (v === null || v === undefined || String(v).trim() === '') continue; // missing handled separately
    const n = Number(v);
    if (isNaN(n)) { issues.push(`${f.label}: not a number`); continue; }
    if (n < f.min || n > f.max) { issues.push(`${f.label} out of range (${f.min}–${f.max})`); continue; }
    if (f.isInt && !Number.isInteger(n)) { issues.push(`${f.label}: must be a whole number`); }
  }

  // Karyo value must be 0, 1, or 2
  const karRaw = row.karyo;
  if (karRaw !== null && karRaw !== undefined && String(karRaw).trim() !== '') {
    const kar = Number(karRaw);
    if (isNaN(kar) || ![0, 1, 2].includes(kar)) {
      issues.push('Cytogenetic Risk must be 0 (Low), 1 (Intermediate), or 2 (High)');
    }
  }

  if (includeDemographics) {
    const ageRaw = row.age;
    if (ageRaw !== null && ageRaw !== undefined && String(ageRaw).trim() !== '') {
      const age = Number(ageRaw);
      if (isNaN(age) || age < 65 || age > 90) {
        issues.push('Age out of range for demographics adjustment (65–90)');
      }
    }
    const sexRaw = row.sex;
    if (sexRaw !== null && sexRaw !== undefined && String(sexRaw).trim() !== '') {
      const sex = Number(sexRaw);
      if (isNaN(sex) || ![0, 1].includes(sex)) {
        issues.push('Sex must be 0 (Female) or 1 (Male)');
      }
    }
  }

  return issues;
}

/**
 * Run batch computation on parsed patient rows.
 * Rows missing required fields or failing range validation are skipped.
 * @param {Object[]} patients
 * @param {boolean} includeDemographics
 * @returns {Object[]} result rows: { id, score, riskClass, skipped, skipReason }
 */
function computeBatch(patients, includeDemographics) {
  return patients.map((row, idx) => {
    const id = row.id !== null && row.id !== undefined && row.id !== ''
      ? row.id
      : `Patient_${idx + 1}`;

    const missing = BATCH_REQUIRED.filter(f => {
      const v = row[f];
      return v === null || v === undefined || String(v).trim() === '';
    });
    if (missing.length > 0) {
      return { id, score: null, riskClass: null, skipped: true,
               skipReason: `Missing required field(s): ${missing.join(', ')}` };
    }

    const issues = validateBatchRow(row, includeDemographics);
    if (issues.length > 0) {
      return { id, score: null, riskClass: null, skipped: true,
               skipReason: `Invalid value(s): ${issues.join('; ')}` };
    }

    const encoded    = encodePatient(row);
    const { score, riskClass } = computeScore(encoded, includeDemographics);
    const confidence = computeConfidence(encoded, riskClass, includeDemographics);
    return { id, score, riskClass, confidence, skipped: false, skipReason: '' };
  });
}

// ─── xlsx template generation ─────────────────────────────────────────────────

/**
 * Generate an xlsx template ArrayBuffer with the correct column headers.
 * @returns {ArrayBuffer}
 */
function generateTemplate() {
  const ws = XLSX.utils.aoa_to_sheet([BATCH_COLUMNS]);
  // Set column widths
  ws['!cols'] = BATCH_COLUMNS.map(c => ({ wch: Math.max(c.length + 2, 10) }));
  const wb = XLSX.utils.book_new();
  XLSX.utils.book_append_sheet(wb, ws, 'iCPSS_Template');
  return XLSX.write(wb, { type: 'array', bookType: 'xlsx' });
}

/**
 * Convert batch results to an xlsx ArrayBuffer for download.
 * @param {Object[]} results - array of { id, score, riskClass, skipped }
 * @returns {ArrayBuffer}
 */
function resultsToXlsx(results) {
  const header = ['ID', 'Class', 'Score', 'Confidence', 'Note'];
  const rows = results.map(r =>
    r.skipped
      ? [r.id, '', '', '', r.skipReason || 'Excluded']
      : [r.id, r.riskClass, r.score.toFixed(4),
         r.confidence !== null ? Math.round(r.confidence.value * 100) + '%' : '100%', '']
  );
  const ws = XLSX.utils.aoa_to_sheet([header, ...rows]);
  ws['!cols'] = [{ wch: 20 }, { wch: 8 }, { wch: 12 }, { wch: 14 }, { wch: 45 }];
  const wb = XLSX.utils.book_new();
  XLSX.utils.book_append_sheet(wb, ws, 'iCPSS_Results');
  return XLSX.write(wb, { type: 'array', bookType: 'xlsx' });
}

// ─── Helper: trigger file download ───────────────────────────────────────────

/**
 * Trigger a browser download of an ArrayBuffer as a file.
 * @param {ArrayBuffer} buffer
 * @param {string} filename
 * @param {string} mimeType
 */
function downloadBuffer(buffer, filename, mimeType) {
  const blob = new Blob([buffer], { type: mimeType });
  const url  = URL.createObjectURL(blob);
  const a    = document.createElement('a');
  a.href     = url;
  a.download = filename;
  document.body.appendChild(a);
  a.click();
  document.body.removeChild(a);
  URL.revokeObjectURL(url);
}

// ─── UI helpers ───────────────────────────────────────────────────────────────

/**
 * Render a result badge + score for single patient.
 * @param {{ score: number, riskClass: string }|null} result
 * @param {{ value: number, n: number }|null} confidence
 */
function renderSingleResult(result, confidence) {
  const container = document.getElementById('result-container');
  if (!result) {
    container.innerHTML = '';
    const sheet = document.getElementById('print-sheet');
    if (sheet) sheet.innerHTML = '';
    return;
  }
  const { score, riskClass } = result;
  const col = CLASS_COLORS[riskClass] || { bg: '#888', text: '#fff' };
  const fullName = CLASS_NAMES[riskClass] || riskClass;
  const classIdx = CLASSES.indexOf(riskClass);
  const outcome  = OUTCOMES[classIdx];

  // Build risk rail stops
  const railStops = CLASSES.map(c => {
    const isActive = c === riskClass;
    const stopCol  = CLASS_COLORS[c] || { bg: '#888' };
    const dotStyle = isActive ? ` style="background:${stopCol.bg}"` : '';
    return `<div class="rail-stop${isActive ? ' active' : ''}"><div class="rail-dot"${dotStyle}></div></div>`;
  }).join('');

  const railLabels = CLASSES.map(c => {
    const isActive = c === riskClass;
    const stopCol  = CLASS_COLORS[c] || { bg: '#888' };
    const lblStyle = isActive ? ` style="color:${stopCol.bg}"` : '';
    return `<span class="rail-lbl${isActive ? ' active' : ''}"${lblStyle}>${c}</span>`;
  }).join('');

  const outcomeBlock = outcome ? `
    <div class="outcome-stats">
      <div class="outcome-stat">
        <div class="outcome-value">${outcome.surv5y}</div>
        <div class="outcome-ci">(95% CI: ${outcome.surv5y_ci})</div>
        <div class="outcome-label">5-year overall survival</div>
      </div>
      <div class="outcome-stat">
        <div class="outcome-value">${outcome.aml3y}</div>
        <div class="outcome-ci">(95% CI: ${outcome.aml3y_ci})</div>
        <div class="outcome-label">3-year AML evolution</div>
      </div>
    </div>` : '';

  container.innerHTML = `
    <div class="result-section">
      <div class="result-class-name" style="color:${col.bg}">${fullName}</div>
      <div class="result-meta">Risk class: <strong>${riskClass}</strong> &ensp;·&ensp; Score: <strong>${score.toFixed(4)}</strong></div>
      ${confidence !== null ? `<div class="confidence-row${confidence.value < 0.80 ? ' confidence-unreliable' : ''}">Confidence: <strong>${Math.round(confidence.value * 100)}%</strong>${confidence.value < 0.80 ? ' <span class="confidence-warn">&#9888; Unreliable prediction</span>' : ''} <span class="confidence-note">(${confidence.n} unknown variable${confidence.n > 1 ? 's' : ''})</span></div>` : ''}
      <div class="risk-rail">
        <div class="risk-rail-track"></div>
        ${railStops}
      </div>
      <div class="rail-labels">${railLabels}</div>
      ${outcomeBlock}
      <div style="margin-top:18px;">
        <button type="button" class="btn-print" onclick="window.print()">Print / Save as PDF</button>
      </div>
    </div>`;
}

/**
 * Populate the hidden print sheet with full result + inputs for printing.
 */
function populatePrintSheet(result, raw, includeDemographics, confidence) {
  const sheet = document.getElementById('print-sheet');
  if (!sheet) return;

  const { score, riskClass } = result;
  const col      = CLASS_COLORS[riskClass] || { bg: '#888' };
  const fullName = CLASS_NAMES[riskClass] || riskClass;
  const classIdx = CLASSES.indexOf(riskClass);
  const outcome  = OUTCOMES[classIdx];

  const now     = new Date();
  const dateStr = now.toLocaleDateString('en-GB', { day: '2-digit', month: 'long', year: 'numeric' });
  const timeStr = now.toLocaleTimeString('en-GB', { hour: '2-digit', minute: '2-digit' });

  const railStops = CLASSES.map(c => {
    const isActive = c === riskClass;
    const sc = CLASS_COLORS[c] || { bg: '#888' };
    return `<div class="rail-stop${isActive ? ' active' : ''}"><div class="rail-dot"${isActive ? ` style="background:${sc.bg}"` : ''}></div></div>`;
  }).join('');

  const railLabels = CLASSES.map(c => {
    const isActive = c === riskClass;
    const sc = CLASS_COLORS[c] || { bg: '#888' };
    return `<span class="rail-lbl${isActive ? ' active' : ''}"${isActive ? ` style="color:${sc.bg}"` : ''}>${c}</span>`;
  }).join('');

  const outHtml = outcome ? `
    <div class="ps-outcomes">
      <div class="ps-stat">
        <div class="ps-stat-val">${outcome.surv5y}</div>
        <div class="ps-stat-ci">95% CI: ${outcome.surv5y_ci}</div>
        <div class="ps-stat-lbl">5-year overall survival</div>
      </div>
      <div class="ps-stat">
        <div class="ps-stat-val">${outcome.aml3y}</div>
        <div class="ps-stat-ci">95% CI: ${outcome.aml3y_ci}</div>
        <div class="ps-stat-lbl">3-year AML evolution</div>
      </div>
    </div>` : '';

  const karyoLabel = { '0': 'Low', '1': 'Intermediate', '2': 'High', '': 'Unknown', 'null': 'Unknown' };
  const mutGenes   = ['ASXL1','DNMT3A','EZH2','RUNX1','SETBP1','STAG2','TET2','TP53','U2AF1'];
  const mutPos     = mutGenes.filter(g => Number(raw[g]) === 1);
  const mutUnk     = mutGenes.filter(g => raw[g] === null || raw[g] === undefined);
  const demoStr    = includeDemographics
    ? `${raw.sex == 1 ? 'Male' : 'Female'}, age ${raw.age}`
    : 'Not included';

  sheet.innerHTML = `
    <div class="ps-watermark">For Investigational Use Only</div>
    <div class="ps-header">
      <div>
        <div class="ps-title"><strong>iCPSS</strong> Risk Score Calculator</div>
        <div class="ps-subtitle">Chronic Myelomonocytic Leukemia</div>
      </div>
      <div class="ps-timestamp">Computed: ${dateStr} &middot; ${timeStr}</div>
    </div>
    <hr class="ps-rule">
    <div class="ps-class-name" style="color:${col.bg}">${fullName}</div>
    <div class="ps-score-meta">Risk class: <strong>${riskClass}</strong> &ensp;&middot;&ensp; Score: <strong>${score.toFixed(4)}</strong>${confidence !== null ? ` &ensp;&middot;&ensp; Confidence: <strong>${Math.round(confidence.value * 100)}%</strong> (${confidence.n} unknown variable${confidence.n > 1 ? 's' : ''})` : ''}</div>
    <div class="risk-rail" style="margin:16px 0 4px"><div class="risk-rail-track"></div>${railStops}</div>
    <div class="rail-labels">${railLabels}</div>
    ${outHtml}
    <hr class="ps-rule">
    <div class="ps-section-label">Parameters Entered</div>
    <table class="ps-inputs-table">
      <tr>
        <td class="ps-il">WBC</td><td class="ps-iv">${raw.wbc} &times;10&#x2079;/L</td>
        <td class="ps-il">Hemoglobin</td><td class="ps-iv">${raw.hb} g/dL</td>
      </tr>
      <tr>
        <td class="ps-il">Platelets</td><td class="ps-iv">${raw.plt} &times;10&#x2079;/L</td>
        <td class="ps-il">BM Blast Count</td><td class="ps-iv">${raw.blasts}%</td>
      </tr>
      <tr>
        <td class="ps-il">Cytogenetic Risk</td><td class="ps-iv">${karyoLabel[String(raw.karyo)] || '&mdash;'}</td>
        <td class="ps-il">Demographics</td><td class="ps-iv">${demoStr}</td>
      </tr>
      <tr>
        <td class="ps-il">Mutations detected</td>
        <td class="ps-iv" colspan="3">${mutPos.length > 0 ? mutPos.join(', ') : 'None'}${mutUnk.length > 0 ? ` <span style="color:#8e939e;font-size:0.78rem">(not tested: ${mutUnk.join(', ')})</span>` : ''}</td>
      </tr>
    </table>
    <hr class="ps-rule">
    <div class="ps-ref">[Citation in press &mdash; reference to be added upon publication]</div>
    <div class="ps-disc"><strong>Investigational use only.</strong> This tool is intended for research purposes and to support &mdash; not replace &mdash; the clinical judgment of qualified healthcare professionals. It has not been reviewed or approved as a medical device by the FDA, EMA, or any other regulatory authority. Risk estimates are derived from a retrospective cohort and carry inherent uncertainty; they should be interpreted in the full clinical context of each individual patient. The authors assume no liability for clinical decisions made on the basis of this tool.</div>`;
}

/**
 * Render the batch results table.
 * @param {Object[]} results
 */
function renderBatchTable(results) {
  const container = document.getElementById('batch-results');
  if (!results || results.length === 0) {
    container.innerHTML = '<p class="empty-msg">No results yet.</p>';
    return;
  }
  const rows = results.map(r => {
    if (r.skipped) {
      return `<tr style="color:var(--text-mute,#8c7b6a)">
        <td>${escHtml(String(r.id))}</td>
        <td style="font-size:0.82rem">—</td>
        <td colspan="2" style="font-size:0.80rem;font-style:italic">${escHtml(r.skipReason || 'Excluded')}</td>
      </tr>`;
    }
    const col         = CLASS_COLORS[r.riskClass] || { bg: '#888', text: '#fff' };
    const confNull    = r.confidence === null;
    const confPct     = confNull ? 100 : Math.round(r.confidence.value * 100);
    const unreliable  = !confNull && r.confidence.value < 0.80;
    const confStr     = confPct + '%' + (unreliable ? ' &#9888;' : '');
    const confColor   = unreliable ? 'var(--warn,#c0392b)' : confNull ? 'var(--text-mute)' : 'var(--text-mid)';
    return `<tr${unreliable ? ' class="row-unreliable"' : ''}>
      <td>${escHtml(String(r.id))}</td>
      <td><span class="badge" style="background:${col.bg};color:${col.text}">${r.riskClass}</span></td>
      <td>${r.score.toFixed(4)}</td>
      <td style="color:${confColor}; font-variant-numeric:tabular-nums">${confStr}</td>
    </tr>`;
  }).join('');
  container.innerHTML = `
    <table class="results-table">
      <thead><tr><th>ID</th><th>Class</th><th>Score</th><th>Confidence</th></tr></thead>
      <tbody>${rows}</tbody>
    </table>`;
}

function escHtml(str) {
  return str.replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');
}

// ─── App initialization ───────────────────────────────────────────────────────

let batchResults = null;

document.addEventListener('DOMContentLoaded', () => {

  // ── Tab switching ─────────────────────────────────────────────────────────
  const tabs    = document.querySelectorAll('.tab-btn');
  const panels  = document.querySelectorAll('.tab-panel');

  tabs.forEach(btn => {
    btn.addEventListener('click', () => {
      tabs.forEach(t => t.classList.remove('active'));
      panels.forEach(p => p.classList.remove('active'));
      btn.classList.add('active');
      document.getElementById('panel-' + btn.dataset.tab).classList.add('active');
    });
  });

  // ── Single patient: compute on button click ───────────────────────────────
  const singleForm = document.getElementById('single-form');
  const singleBtn  = document.getElementById('btn-compute');
  const clearBtn   = document.getElementById('btn-clear');
  const demoToggle = document.getElementById('demo-toggle');

  singleBtn.addEventListener('click', () => {
    const raw = readSingleForm();
    const { errors, warnings } = validateSingleForm(raw, demoToggle.checked);
    if (errors.length > 0) {
      showErrors(errors, 'single-errors');
      clearWarnings('single-warnings');
      return;
    }
    clearErrors('single-errors');
    if (warnings.length > 0) {
      showWarnings(warnings, 'single-warnings');
    } else {
      clearWarnings('single-warnings');
    }
    const encoded    = encodePatient(raw);
    const result     = computeScore(encoded, demoToggle.checked);
    const confidence = computeConfidence(encoded, result.riskClass, demoToggle.checked);
    renderSingleResult(result, confidence);
    populatePrintSheet(result, raw, demoToggle.checked, confidence);
  });

  clearBtn.addEventListener('click', () => {
    singleForm.reset();
    // Sync demo section visibility after form reset (reset() doesn't fire 'change')
    document.getElementById('demo-toggle').dispatchEvent(new Event('change'));
    renderSingleResult(null);
    clearErrors();
    clearWarnings('single-warnings');
  });

  // ── Batch: template download ──────────────────────────────────────────────
  document.getElementById('btn-template').addEventListener('click', () => {
    const buf = generateTemplate();
    downloadBuffer(buf, 'icpss_template.xlsx',
      'application/vnd.openxmlformats-officedocument.spreadsheetml.sheet');
  });

  // ── Batch: file upload + auto-compute ────────────────────────────────────
  const batchFile = document.getElementById('batch-file');
  const batchDemo = document.getElementById('batch-demo-toggle');

  batchFile.addEventListener('change', () => {
    runBatch();
  });

  batchDemo.addEventListener('change', () => {
    // Re-compute if a file is already loaded
    if (batchFile.files.length > 0) {
      runBatch();
    }
  });

  function runBatch() {
    const file = batchFile.files[0];
    if (!file) return;
    clearErrors('batch-errors');
    clearWarnings('batch-warnings');
    document.getElementById('batch-results').innerHTML = '<p class="loading-msg">Computing…</p>';

    const reader = new FileReader();
    reader.onload = (e) => {
      const { patients, errors } = parseXlsx(e.target.result);
      if (errors.length > 0) {
        showErrors(errors, 'batch-errors');
        document.getElementById('batch-results').innerHTML = '';
        return;
      }
      batchResults = computeBatch(patients, batchDemo.checked);
      renderBatchTable(batchResults);
      document.getElementById('btn-download-results').disabled = false;
      const skipped = batchResults.filter(r => r.skipped).length;
      if (skipped > 0) {
        const s = skipped === 1;
        showWarnings([
          `${skipped} ${s ? 'row was' : 'rows were'} excluded from scoring due to missing or invalid data. ` +
          `${s ? 'It is' : 'They are'} flagged in the table and in the downloaded file.`
        ], 'batch-warnings');
      } else {
        clearWarnings('batch-warnings');
      }
    };
    reader.onerror = () => {
      showErrors(['Failed to read the file. Please try again.'], 'batch-errors');
      document.getElementById('batch-results').innerHTML = '';
    };
    reader.readAsArrayBuffer(file);
  }

  // ── Batch: download results ───────────────────────────────────────────────
  document.getElementById('btn-download-results').addEventListener('click', () => {
    if (!batchResults || batchResults.length === 0) return;
    const buf = resultsToXlsx(batchResults);
    const suffix = batchDemo.checked ? 'withDemo' : 'noDemo';
    downloadBuffer(buf, `icpss_results_${suffix}.xlsx`,
      'application/vnd.openxmlformats-officedocument.spreadsheetml.sheet');
  });
});

// ─── Form helpers ─────────────────────────────────────────────────────────────

function readSingleForm() {
  const g = (id) => document.getElementById(id);
  const sexEl   = document.querySelector('input[name="sex"]:checked');
  const karyoEl = document.querySelector('input[name="karyo"]:checked');

  // Read 3-state gene radio: Mut=1, WT=0, Unknown (value="")=null
  const geneVal = (name) => {
    const el = document.querySelector(`input[name="${name}"]:checked`);
    if (!el || el.value === '') return null;
    return Number(el.value);
  };

  return {
    sex:    sexEl   ? Number(sexEl.value) : null,
    age:    g('age').value,
    wbc:    g('wbc').value,
    hb:     g('hb').value,
    plt:    g('plt').value,
    blasts: g('blasts').value,
    karyo:  karyoEl ? karyoEl.value : null,
    ASXL1:  geneVal('ASXL1'),
    DNMT3A: geneVal('DNMT3A'),
    EZH2:   geneVal('EZH2'),
    RUNX1:  geneVal('RUNX1'),
    SETBP1: geneVal('SETBP1'),
    STAG2:  geneVal('STAG2'),
    TET2:   geneVal('TET2'),
    TP53:   geneVal('TP53'),
    U2AF1:  geneVal('U2AF1'),
  };
}

function validateSingleForm(raw, includeDemographics) {
  const errors = [];
  const warnings = [];

  // Required clinical fields
  const clinicalRequired = [
    { val: raw.wbc,   label: 'WBC'                     },
    { val: raw.hb,    label: 'Hemoglobin'              },
    { val: raw.plt,   label: 'Platelets'               },
    { val: raw.blasts, label: 'Bone Marrow Blast Count' },
  ];
  for (const f of clinicalRequired) {
    if (f.val === '' || f.val === null || f.val === undefined) {
      errors.push(`${f.label} is required.`);
    }
  }

  // Karyo is optional — when missing, confidence score reflects the uncertainty

  // Range + type validation — limits drawn from FIELD_RULES to stay in sync with batch mode
  const WARN_ABOVE = { wbc: 100, hb: 15 };
  const rangeFields = FIELD_RULES.map(f => ({
    val: raw[f.key], label: f.label, min: f.min, max: f.max, isInt: f.isInt,
    warnAbove: WARN_ABOVE[f.key] ?? null,
  }));
  for (const f of rangeFields) {
    if (f.val === '' || f.val === null || f.val === undefined) continue;
    const n = Number(f.val);
    if (isNaN(n)) {
      errors.push(`${f.label} must be a number.`);
      continue;
    }
    if (n < f.min || n > f.max) {
      errors.push(`${f.label} must be between ${f.min} and ${f.max}.`);
      continue;
    }
    if (f.isInt && !Number.isInteger(n)) {
      errors.push(`${f.label} must be a whole number.`);
      continue;
    }
    if (f.warnAbove !== null && n > f.warnAbove) {
      warnings.push(`${f.label} value (${n}) is unusually high — please verify.`);
    }
  }

  // Demographics validation (when toggle ON)
  if (includeDemographics) {
    if (raw.sex === null || raw.sex === undefined) {
      errors.push('Sex is required when demographics adjustment is enabled.');
    }
    if (raw.age === '' || raw.age === null || isNaN(Number(raw.age))) {
      errors.push('Age is required when demographics adjustment is enabled.');
    } else {
      const age = Number(raw.age);
      if (age < 65 || age > 90) {
        errors.push('Age must be between 65 and 90 for demographics-adjusted scoring.');
      }
    }
  }

  return { errors, warnings };
}

function showErrors(errs, containerId) {
  const id = containerId || 'single-errors';
  const el = document.getElementById(id);
  if (!el) return;
  el.innerHTML = errs.map(e => `<li>${escHtml(e)}</li>`).join('');
  el.style.display = 'block';
}

function clearErrors(containerId) {
  if (containerId) {
    const el = document.getElementById(containerId);
    if (el) { el.innerHTML = ''; el.style.display = 'none'; }
  } else {
    document.querySelectorAll('.error-list').forEach(el => {
      el.innerHTML = '';
      el.style.display = 'none';
    });
  }
}

function showWarnings(warns, containerId) {
  const el = document.getElementById(containerId);
  if (!el) return;
  el.innerHTML = warns.map(w => `<li>${escHtml(w)}</li>`).join('');
  el.style.display = 'block';
}

function clearWarnings(containerId) {
  const el = document.getElementById(containerId);
  if (el) { el.innerHTML = ''; el.style.display = 'none'; }
}
