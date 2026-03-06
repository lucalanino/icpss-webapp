// ── Extra UI wiring ─────────────────────────────────────────────────────────

const demoToggle  = document.getElementById('demo-toggle');
const demoSection = document.getElementById('demo-section');

function updateDemoSection() {
  demoSection.style.display = demoToggle.checked ? 'block' : 'none';
  if (!demoToggle.checked) {
    document.querySelectorAll('input[name="sex"]').forEach(r => r.checked = false);
    document.getElementById('age').value = '';
  }
}
demoToggle.addEventListener('change', updateDemoSection);
updateDemoSection();

// File name display
const batchFile     = document.getElementById('batch-file');
const batchFilename = document.getElementById('batch-filename');
batchFile.addEventListener('click', () => { batchFile.value = ''; });
batchFile.addEventListener('change', () => {
  batchFilename.textContent = batchFile.files.length > 0
    ? 'Selected: ' + batchFile.files[0].name : '';
});

// Drag-and-drop
const uploadArea = document.querySelector('.batch-upload-area');
uploadArea.addEventListener('dragover', e => {
  e.preventDefault();
  uploadArea.style.borderColor = 'var(--accent)';
});
uploadArea.addEventListener('dragleave', () => {
  uploadArea.style.borderColor = '';
});
uploadArea.addEventListener('drop', e => {
  e.preventDefault();
  uploadArea.style.borderColor = '';
  const files = e.dataTransfer.files;
  if (files.length > 0 && files[0].name.endsWith('.xlsx')) {
    const dt = new DataTransfer();
    dt.items.add(files[0]);
    batchFile.files = dt.files;
    batchFilename.textContent = 'Selected: ' + files[0].name;
    batchFile.dispatchEvent(new Event('change'));
  } else {
    alert('Please drop a .xlsx file.');
  }
});
