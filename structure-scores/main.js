// ---------- score formulas ---------------------------------------------
function tmScore(d, L) {
  const d0 = Math.max(1.24 * Math.cbrt(L - 15) - 1.8, 0.5);
  let s = 0;
  for (const di of d) s += 1 / (1 + (di / d0) ** 2);
  return s / L;
}
function pAtMost(d, L, thr) {
  let n = 0;
  for (const di of d) if (di <= thr) n++;
  return 100 * n / L;
}
function gdtTs(d, L) {
  return (pAtMost(d, L, 1) + pAtMost(d, L, 2)
        + pAtMost(d, L, 4) + pAtMost(d, L, 8)) / 4;
}
function gdtHa(d, L) {
  return (pAtMost(d, L, 0.5) + pAtMost(d, L, 1)
        + pAtMost(d, L, 2)   + pAtMost(d, L, 4)) / 4;
}
function rmsd(d, L) {
  let s = 0; for (const di of d) s += di * di;
  return Math.sqrt(s / L);
}

// ---------- state -------------------------------------------------------
let L = 12;                                          // protein length
let dist = Array.from({ length: L }, () => 0.8 + Math.random() * 0.6);

const MAX_DIST = 15;       // sliders + bar scale
const THRESHOLDS = [
  { v: 0.5, label: '0.5' }, { v: 1, label: '1' },
  { v: 2,   label: '2'   }, { v: 4, label: '4' },
  { v: 8,   label: '8 Å' },
];

// ---------- DOM refs ----------------------------------------------------
const $ = (id) => document.getElementById(id);
const Ls = $('L-slider'), Ln = $('L-num');
const sTm = $('s-tm'), sTs = $('s-gdtts'), sHa = $('s-gdtha'), sRm = $('s-rmsd');
const residuesEl = $('residues');
const formulaEl = $('formulas');

// ---------- bar colour by distance band -------------------------------
function bandColor(d) {
  if (d < 0.5)  return '#1eaa5c';
  if (d < 1)    return '#56c270';
  if (d < 2)    return '#c8d44e';
  if (d < 4)    return '#e6a34a';
  if (d < 8)    return '#d8625e';
  return '#5a4a52';
}

// ---------- build / rebuild residue rows -------------------------------
function buildRows() {
  residuesEl.innerHTML = '';
  for (let i = 0; i < L; i++) {
    const row = document.createElement('div');
    row.className = 'row';
    row.innerHTML = `
      <div class="name">res ${String(i + 1).padStart(2, '0')}</div>
      <input type="range" min="0" max="${MAX_DIST}" step="0.1" value="${dist[i]}">
      <div class="barwrap"><div class="bar"></div></div>
      <div class="dist"></div>
    `;
    // threshold lines
    const wrap = row.querySelector('.barwrap');
    THRESHOLDS.forEach((t, k) => {
      const line = document.createElement('div');
      line.className = 'threshold' + (i === 0 ? ' has-label' : '');
      line.style.left = (100 * t.v / MAX_DIST) + '%';
      if (i === 0) line.dataset.l = t.label;
      wrap.appendChild(line);
    });

    const slider = row.querySelector('input[type=range]');
    slider.addEventListener('input', (e) => {
      dist[i] = parseFloat(e.target.value);
      update();
    });
    residuesEl.appendChild(row);
  }
  update();
}

// ---------- redraw bars + score cards + formulas -----------------------
function update() {
  const rows = residuesEl.querySelectorAll('.row');
  rows.forEach((row, i) => {
    const d = dist[i];
    row.querySelector('.bar').style.width = (100 * d / MAX_DIST) + '%';
    row.querySelector('.bar').style.background = bandColor(d);
    row.querySelector('.dist').textContent = d.toFixed(2) + ' Å';
  });

  const tm = tmScore(dist, L);
  const ts = gdtTs(dist, L);
  const ha = gdtHa(dist, L);
  const rm = rmsd(dist, L);
  sTm.textContent = tm.toFixed(3);
  sTs.textContent = ts.toFixed(1);
  sHa.textContent = ha.toFixed(1);
  sRm.textContent = rm.toFixed(2) + ' Å';

  // Color the TM card by fold-similarity convention.
  sTm.style.color = tm >= 0.5 ? '#7eea9d' : tm >= 0.17 ? '#ffd479' : '#e08a82';

  // Detailed breakdown
  const d0 = Math.max(1.24 * Math.cbrt(L - 15) - 1.8, 0.5);
  const p05 = pAtMost(dist, L, 0.5);
  const p1  = pAtMost(dist, L, 1);
  const p2  = pAtMost(dist, L, 2);
  const p4  = pAtMost(dist, L, 4);
  const p8  = pAtMost(dist, L, 8);

  // per-residue TM contributions
  const tmTerms = dist.map(di => 1 / (1 + (di / d0) ** 2));
  const tmSum   = tmTerms.reduce((a, b) => a + b, 0);

  formulaEl.innerHTML =
`<span class="lbl">TM-score</span>
  L = ${L},  d₀ = max(1.24·∛(L-15) − 1.8, 0.5) = ${d0.toFixed(3)} Å
  Σᵢ 1 / (1 + (dᵢ/d₀)²) = ${tmSum.toFixed(4)}
  TM = ${tmSum.toFixed(4)} / ${L} = <b>${tm.toFixed(4)}</b>

<span class="lbl">GDT_TS</span>
  P₁ = ${p1.toFixed(1)}   P₂ = ${p2.toFixed(1)}   P₄ = ${p4.toFixed(1)}   P₈ = ${p8.toFixed(1)}
  GDT_TS = (${p1.toFixed(1)} + ${p2.toFixed(1)} + ${p4.toFixed(1)} + ${p8.toFixed(1)}) / 4 = <b>${ts.toFixed(2)}</b>

<span class="lbl">GDT_HA</span>
  P₀.₅ = ${p05.toFixed(1)}   P₁ = ${p1.toFixed(1)}   P₂ = ${p2.toFixed(1)}   P₄ = ${p4.toFixed(1)}
  GDT_HA = (${p05.toFixed(1)} + ${p1.toFixed(1)} + ${p2.toFixed(1)} + ${p4.toFixed(1)}) / 4 = <b>${ha.toFixed(2)}</b>

<span class="lbl">RMSD</span>
  Σᵢ dᵢ² = ${dist.reduce((a, b) => a + b * b, 0).toFixed(3)}
  RMSD = √(${dist.reduce((a, b) => a + b * b, 0).toFixed(3)} / ${L}) = <b>${rm.toFixed(3)} Å</b>`;
}

// ---------- L control ---------------------------------------------------
function setL(newL, regenDist = true) {
  L = Math.max(5, Math.min(30, newL | 0));
  Ls.value = L; Ln.value = L;
  if (regenDist) {
    if (dist.length < L) {
      while (dist.length < L) dist.push(1.0);
    } else {
      dist.length = L;
    }
  }
  buildRows();
}
Ls.addEventListener('input', (e) => setL(+e.target.value));
Ln.addEventListener('change', (e) => setL(+e.target.value));

// ---------- presets -----------------------------------------------------
const PRESETS = [
  {
    name: 'AlphaFold-like',
    desc: 'Most residues sub-Å, a few up to ~2 Å.',
    fn: (n) => Array.from({ length: n }, (_, i) =>
      0.3 + Math.random() * 0.7 + (i % 7 === 0 ? 1.5 * Math.random() : 0)),
  },
  {
    name: 'Mediocre',
    desc: 'Avg ~3 Å, mixed quality.',
    fn: (n) => Array.from({ length: n }, () => 1 + Math.random() * 5),
  },
  {
    name: 'Same fold, low accuracy',
    desc: 'Many 2–5 Å, a few outliers.',
    fn: (n) => Array.from({ length: n }, () => 1.5 + Math.random() * 4),
  },
  {
    name: 'Wrong fold',
    desc: 'Distances scattered, mostly > 8 Å.',
    fn: (n) => Array.from({ length: n }, () => 5 + Math.random() * 9),
  },
  {
    name: 'One bad outlier',
    desc: '11 perfect residues + 1 at 14 Å. RMSD blows up, TM/GDT stay good.',
    fn: (n) => {
      const a = Array.from({ length: n }, () => 0.3 + Math.random() * 0.4);
      a[a.length - 1] = 14;
      return a;
    },
  },
  {
    name: 'All identical',
    desc: 'Every residue at 1.5 Å — clean illustration.',
    fn: (n) => Array.from({ length: n }, () => 1.5),
  },
];

const presetsEl = $('presets');
PRESETS.forEach((p) => {
  const b = document.createElement('button');
  b.className = 'preset';
  b.textContent = p.name;
  b.title = p.desc;
  b.addEventListener('click', () => {
    dist = p.fn(L);
    buildRows();
  });
  presetsEl.appendChild(b);
});

// ---------- go ----------------------------------------------------------
buildRows();
