'use strict';

// ─── Static model parameters (ported from R) ─────────────────────────────────
// Feature order (28 features):
// [age, wbc_cont, hb_cont, plt_cont, blast_cont,
//  sex1, sex0,
//  ASXL11, ASXL10, DNMT3A1, DNMT3A0, EZH21, EZH20,
//  RUNX11, RUNX10, SETBP11, SETBP10, STAG21, STAG20,
//  TET21, TET20, TP531, TP530, U2AF11, U2AF10,
//  cpss_karyo2, cpss_karyo1, cpss_karyo0]

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

const THRESHOLDS = [-0.5935335, 0.1467796, 0.5660161, 1.2584983];
const CLASSES = ['VL', 'L', 'I', 'H', 'VH'];

// Class color palette
const CLASS_COLORS = {
  VL: { bg: '#1a7f37', text: '#ffffff' },
  L:  { bg: '#2da44e', text: '#ffffff' },
  I:  { bg: '#f5a623', text: '#ffffff' },
  H:  { bg: '#cf222e', text: '#ffffff' },
  VH: { bg: '#82071e', text: '#ffffff' }
};

// ─── Encoding ─────────────────────────────────────────────────────────────────

/**
 * Build the 28-element encoded feature vector from raw patient inputs.
 * Missing values are represented as NaN.
 *
 * @param {Object} raw - Raw patient inputs:
 *   sex (1=Male, 0=Female), age, wbc_cont, hb_cont, plt_cont, blast_cont,
 *   cpss_karyo (0/1/2), ASXL1, DNMT3A, EZH2, RUNX1, SETBP1, STAG2, TET2, TP53, U2AF1
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
  const wbc = toNum(raw.wbc_cont);
  const hb  = toNum(raw.hb_cont);
  const plt = toNum(raw.plt_cont);
  const bla = toNum(raw.blast_cont);
  const kar = toNum(raw.cpss_karyo);

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
    // cpss_karyo: three-way one-hot [2, 1, 0]
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

// ─── Batch: xlsx parsing ──────────────────────────────────────────────────────

const BATCH_COLUMNS = ['id', 'sex', 'age', 'wbc_cont', 'hb_cont', 'plt_cont',
  'blast_cont', 'cpss_karyo', 'ASXL1', 'DNMT3A', 'EZH2', 'RUNX1',
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

/**
 * Run batch computation on parsed patient rows.
 * @param {Object[]} patients
 * @param {boolean} includeDemographics
 * @returns {Object[]} result rows: { id, score, riskClass }
 */
function computeBatch(patients, includeDemographics) {
  return patients.map((row, idx) => {
    const id = row.id !== null && row.id !== undefined && row.id !== ''
      ? row.id
      : `Patient_${idx + 1}`;
    const encoded = encodePatient(row);
    const { score, riskClass } = computeScore(encoded, includeDemographics);
    return { id, score, riskClass };
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
 * @param {Object[]} results - array of { id, score, riskClass }
 * @returns {ArrayBuffer}
 */
function resultsToXlsx(results) {
  const header = ['ID', 'Class', 'Score'];
  const rows = results.map(r => [r.id, r.riskClass, r.score.toFixed(4)]);
  const ws = XLSX.utils.aoa_to_sheet([header, ...rows]);
  ws['!cols'] = [{ wch: 20 }, { wch: 8 }, { wch: 12 }];
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
 */
function renderSingleResult(result) {
  const container = document.getElementById('result-container');
  if (!result) {
    container.innerHTML = '';
    return;
  }
  const { score, riskClass } = result;
  const col = CLASS_COLORS[riskClass] || { bg: '#888', text: '#fff' };
  container.innerHTML = `
    <div class="result-card">
      <div class="result-badge" style="background:${col.bg};color:${col.text}">
        ${riskClass}
      </div>
      <div class="result-score">Score: <strong>${score.toFixed(4)}</strong></div>
    </div>`;
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
    const col = CLASS_COLORS[r.riskClass] || { bg: '#888', text: '#fff' };
    return `<tr>
      <td>${escHtml(String(r.id))}</td>
      <td><span class="badge" style="background:${col.bg};color:${col.text}">${r.riskClass}</span></td>
      <td>${r.score.toFixed(4)}</td>
    </tr>`;
  }).join('');
  container.innerHTML = `
    <table class="results-table">
      <thead><tr><th>ID</th><th>Class</th><th>Score</th></tr></thead>
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
    clearErrors();
    const raw = readSingleForm();
    const errs = validateSingleForm(raw, demoToggle.checked);
    if (errs.length > 0) {
      showErrors(errs, 'single-errors');
      return;
    }
    const encoded = encodePatient(raw);
    const result  = computeScore(encoded, demoToggle.checked);
    renderSingleResult(result);
  });

  clearBtn.addEventListener('click', () => {
    singleForm.reset();
    renderSingleResult(null);
    clearErrors();
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
  const sexEl = document.querySelector('input[name="sex"]:checked');
  return {
    sex:        sexEl ? Number(sexEl.value) : null,
    age:        g('age').value,
    wbc_cont:   g('wbc').value,
    hb_cont:    g('hb').value,
    plt_cont:   g('plt').value,
    blast_cont: g('blasts').value,
    cpss_karyo: g('karyo').value !== '' ? g('karyo').value : null,
    ASXL1:      g('ASXL1').checked  ? 1 : 0,
    DNMT3A:     g('DNMT3A').checked ? 1 : 0,
    EZH2:       g('EZH2').checked   ? 1 : 0,
    RUNX1:      g('RUNX1').checked  ? 1 : 0,
    SETBP1:     g('SETBP1').checked ? 1 : 0,
    STAG2:      g('STAG2').checked  ? 1 : 0,
    TET2:       g('TET2').checked   ? 1 : 0,
    TP53:       g('TP53').checked   ? 1 : 0,
    U2AF1:      g('U2AF1').checked  ? 1 : 0,
  };
}

function validateSingleForm(raw, includeDemographics) {
  const errs = [];

  if (includeDemographics) {
    if (raw.sex === null || raw.sex === undefined) errs.push('Sex is required when demographics are included.');
    if (raw.age === '' || raw.age === null || isNaN(Number(raw.age))) {
      errs.push('Age is required and must be a number.');
    } else if (Number(raw.age) < 0 || Number(raw.age) > 120) {
      errs.push('Age must be between 0 and 120.');
    }
  }

  const numFields = [
    { val: raw.wbc_cont,   label: 'WBC',        min: 0,   max: 500  },
    { val: raw.hb_cont,    label: 'Hemoglobin',  min: 0,   max: 25   },
    { val: raw.plt_cont,   label: 'Platelets',   min: 0,   max: 3000 },
    { val: raw.blast_cont, label: 'BM Blasts',   min: 0,   max: 100  },
  ];
  for (const f of numFields) {
    if (f.val !== '' && f.val !== null && f.val !== undefined) {
      const n = Number(f.val);
      if (isNaN(n)) {
        errs.push(`${f.label} must be a number.`);
      } else if (n < f.min || n > f.max) {
        errs.push(`${f.label} must be between ${f.min} and ${f.max}.`);
      }
    }
  }

  if (raw.cpss_karyo === null || raw.cpss_karyo === '' || raw.cpss_karyo === undefined) {
    errs.push('CPSS Karyotype is required.');
  }

  return errs;
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
