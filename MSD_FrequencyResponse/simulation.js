// Mass-Spring-Damper Frequency Response Simulator
// Ported from MSDForceExcitationSim.m (MATLAB App Designer)

'use strict';

// ── Constants ──────────────────────────────────────────────────────────────
const FPS = 60;
const DT  = 1 / FPS;
const COLORS = {
  r1: '#FFCB05', r1Dark: '#a07800',
  r2: '#00274C',
  spring: '#555555', damper: '#888888',
  arrow: '#c00000',
  text: '#1d1d1f', textDim: '#6e6e73',
  grid: '#e8e8ed', zeroline: '#c8c8cd',
  bg: '#ffffff'
};

// ── Canvas references ──────────────────────────────────────────────────────
const animCanvas     = document.getElementById('anim-canvas');
const posCanvas      = document.getElementById('pos-canvas');
const posforceCanvas = document.getElementById('posforce-canvas');
const magCanvas      = document.getElementById('mag-canvas');
const phaseCanvas    = document.getElementById('phase-canvas');
const animCtx     = animCanvas.getContext('2d');
const posCtx      = posCanvas.getContext('2d');
const posforceCtx = posforceCanvas.getContext('2d');
const magCtx      = magCanvas.getContext('2d');
const phaseCtx    = phaseCanvas.getContext('2d');

// ── Input helpers ──────────────────────────────────────────────────────────
const getNum = id => parseFloat(document.getElementById(id).value) || 0;

let waveform = 'sin'; // 'sin' (default) or 'cos'

function getParams() {
  return {
    m:    Math.min(100,   Math.max(0.01,  getNum('fr-m'))),
    k:    Math.min(50000, Math.max(0.1,   getNum('fr-k'))),
    c:    Math.min(10000, Math.max(0,     getNum('fr-c'))),
    F0:   Math.min(1000,  Math.max(0,     getNum('fr-F0'))),
    w:    Math.min(50,    Math.max(0,     getNum('fr-w'))),
    x0:   Math.min(10,    Math.max(-10,   getNum('fr-x0'))),
    v0:   Math.min(50,    Math.max(-50,   getNum('fr-v0'))),
    tend: Math.min(30,    Math.max(0.1,   getNum('fr-tend'))),
    waveform,
  };
}

// ── Physics ────────────────────────────────────────────────────────────────
function systemProps(p) {
  const wn   = Math.sqrt(p.k / p.m);
  const zeta = p.c / (2 * Math.sqrt(p.k * p.m));
  return { wn, zeta };
}

// Steady-state particular solution coefficients for F0*cos(w*t).
// x_p(t) = Xc*cos(w*t) + Xs*sin(w*t)
// At w=0 this collapses to Xc = F0/k (static deflection), giving a step response.
function particularCoeffs(p, wn, zeta) {
  const { F0, w, m } = p;
  if (F0 === 0) return { Xc: 0, Xs: 0 };
  const f0    = F0 / m;
  const delta = Math.pow(wn * wn - w * w, 2) + Math.pow(2 * zeta * w * wn, 2);
  if (delta < 1e-12) return { Xc: 0, Xs: 0 };
  const Xc = (wn * wn - w * w)    * f0 / delta;
  const Xs = (2 * zeta * w * wn)  * f0 / delta;
  return { Xc, Xs };
}

// Closed-form total response: x(t) = x_h(t) + x_p(t),
// where x_h satisfies x(0)=x0, x'(0)=v0 with the particular subtracted out.
//
// Particular coefficients (Xc, Xs) come from a cos-driven derivation. To
// support sin forcing without re-deriving, substitute cos→sin and sin→−cos
// in the particular solution. Equivalently:
//   cos: xp(t) = Xc cos + Xs sin,   xp(0) =  Xc,   xp'(0) =  Xs·w
//   sin: xp(t) = Xc sin − Xs cos,   xp(0) = −Xs,   xp'(0) =  Xc·w
function computeResponse(p) {
  const { wn, zeta } = systemProps(p);
  const { Xc, Xs }   = particularCoeffs(p, wn, zeta);
  const { x0, v0, tend, F0, w, waveform } = p;
  const isSin = waveform === 'sin';

  // Particular evaluated at t=0, used to pick homogeneous ICs
  const xp0 = isSin ? -Xs       :  Xc;
  const vp0 = isSin ?  Xc * w   :  Xs * w;
  const xh0 = x0 - xp0;
  const vh0 = v0 - vp0;

  const N = Math.ceil(tend * FPS) + 1;
  const t = new Float64Array(N);
  const x = new Float64Array(N);
  const F = new Float64Array(N);

  for (let i = 0; i < N; i++) {
    const ti = i * DT;
    t[i] = ti;

    // Homogeneous part (free-response solution with ICs xh0, vh0)
    let xh;
    if (zeta < 1 - 1e-6) {
      const wd = wn * Math.sqrt(1 - zeta * zeta);
      xh = Math.exp(-zeta * wn * ti) *
           (xh0 * Math.cos(wd * ti) + (zeta * wn * xh0 + vh0) / wd * Math.sin(wd * ti));
    } else if (zeta > 1 + 1e-6) {
      const sq = Math.sqrt(zeta * zeta - 1);
      const r1 = wn * (-zeta + sq);
      const r2 = wn * (-zeta - sq);
      const A1 = (vh0 + xh0 * wn * (zeta + sq)) / (2 * wn * sq);
      const A2 = (-vh0 - xh0 * wn * (zeta - sq)) / (2 * wn * sq);
      xh = A1 * Math.exp(r1 * ti) + A2 * Math.exp(r2 * ti);
    } else {
      xh = (xh0 + (vh0 + xh0 * wn) * ti) * Math.exp(-wn * ti);
    }

    const cw = Math.cos(w * ti);
    const sw = Math.sin(w * ti);
    const xp = isSin ? (Xc * sw - Xs * cw) : (Xc * cw + Xs * sw);

    x[i] = xh + xp;
    F[i] = isSin ? F0 * sw : F0 * cw;
  }

  // Operating point on FR curves
  const r = wn > 0 ? w / wn : 0;
  const M = (F0 > 0 && p.k > 0) ? Math.sqrt(Xc * Xc + Xs * Xs) * p.k / F0 : 0;
  const phi = (Math.abs(r) < 1e-9 && zeta === 0) ? 0
            : -Math.atan2(2 * zeta * r, 1 - r * r) * 180 / Math.PI;

  return { t, x, F, wn, zeta, Xc, Xs, r, M, phi };
}

// FR sweep curves: M(r), phi(r) for r in [0, 2.5]
function frSweep(zeta) {
  const N = 600;
  const rmax = 2.5;
  const rvect = new Float64Array(N);
  const M     = new Float64Array(N);
  const phi   = new Float64Array(N);
  for (let i = 0; i < N; i++) {
    const r = (i / (N - 1)) * rmax;
    rvect[i] = r;
    M[i] = 1 / Math.sqrt(Math.pow(1 - r * r, 2) + Math.pow(2 * zeta * r, 2));
    phi[i] = -Math.atan2(2 * zeta * r, 1 - r * r) * 180 / Math.PI;
  }
  return { rvect, M, phi };
}

// ── State ──────────────────────────────────────────────────────────────────
let resp     = null;
let sweep    = null;
let activeTab = 'pos';
let animRaf   = null;
let animStartTs = null;

// ── Derived label update ───────────────────────────────────────────────────
function fmt(n, digits) { return isFinite(n) ? n.toFixed(digits) : '—'; }

function updateDerivedLabels() {
  if (!resp) return;
  document.getElementById('fr-wn').textContent   = fmt(resp.wn, 3);
  document.getElementById('fr-zeta').textContent = fmt(resp.zeta, 4);
  document.getElementById('fr-M').textContent    = fmt(resp.M, 3);
  document.getElementById('fr-r').textContent    = fmt(resp.r, 3);
  document.getElementById('fr-phi').textContent  = fmt(resp.phi, 1);
  const type = resp.zeta < 1 - 1e-4 ? 'underdamped'
             : resp.zeta > 1 + 1e-4 ? 'overdamped'
             : 'critical';
  document.getElementById('fr-damp-type').textContent = type;
}

// ── Accessibility: hidden text descriptions ────────────────────────────────
function updateA11yDescriptions() {
  if (!resp) return;
  const p = getParams();
  const type = resp.zeta < 1 - 1e-4 ? 'underdamped'
             : resp.zeta > 1 + 1e-4 ? 'overdamped'
             : 'critically damped';

  document.getElementById('anim-description').textContent =
    `Mass-spring-damper with harmonic forcing. Natural frequency ${fmt(resp.wn, 3)} rad/s, ` +
    `damping ratio ${fmt(resp.zeta, 4)} (${type}). Applied force F(t) = ${fmt(p.F0, 2)} ${p.waveform}(${fmt(p.w, 2)}t) N. ` +
    `Initial position ${fmt(p.x0, 2)} m, initial velocity ${fmt(p.v0, 2)} m/s.`;

  const xMax = maxAbs(resp.x);
  document.getElementById('pos-description').textContent =
    `Position versus time. Time from 0 to ${fmt(p.tend, 1)} s. Peak displacement ${fmt(xMax, 3)} m.`;

  const fMax = maxAbs(resp.F);
  document.getElementById('posforce-description').textContent =
    `Position and applied force versus time. Position peak ${fmt(xMax, 3)} m. Force peak ${fmt(fMax, 2)} N.`;

  document.getElementById('mag-description').textContent =
    `Magnitude frequency response. Frequency ratio r = ${fmt(resp.r, 3)}. ` +
    `Dynamic magnification factor M = ${fmt(resp.M, 3)}.`;

  document.getElementById('phase-description').textContent =
    `Phase frequency response. Frequency ratio r = ${fmt(resp.r, 3)}. ` +
    `Phase shift ${fmt(resp.phi, 1)} degrees.`;
}

function maxAbs(arr) {
  let m = 0;
  for (let i = 0; i < arr.length; i++) {
    const a = Math.abs(arr[i]);
    if (a > m) m = a;
  }
  return m;
}

// ── CSV export ─────────────────────────────────────────────────────────────
function exportDataAsCSV() {
  if (!resp) return;
  const rows = ['time (s),position (m),force (N)'];
  for (let i = 0; i < resp.t.length; i++) {
    rows.push(`${resp.t[i].toFixed(6)},${resp.x[i].toFixed(6)},${resp.F[i].toFixed(6)}`);
  }
  const blob = new Blob([rows.join('\n')], { type: 'text/csv;charset=utf-8;' });
  const url  = URL.createObjectURL(blob);
  const a    = document.createElement('a');
  a.href     = url;
  a.download = 'msd-frequency-response-data.csv';
  a.style.visibility = 'hidden';
  document.body.appendChild(a);
  a.click();
  document.body.removeChild(a);
  URL.revokeObjectURL(url);
}

// ── Canvas resize ──────────────────────────────────────────────────────────
function fitCanvas(canvas) {
  const w = canvas.offsetWidth  * devicePixelRatio | 0;
  const h = canvas.offsetHeight * devicePixelRatio | 0;
  if (w === 0 || h === 0) return false;
  if (canvas.width !== w || canvas.height !== h) {
    canvas.width  = w;
    canvas.height = h;
    return true;
  }
  return false;
}

function resizeAll() {
  fitCanvas(animCanvas);
  // Only fit visible plot canvases (hidden panels have zero offset size)
  if (activeTab === 'pos')      fitCanvas(posCanvas);
  if (activeTab === 'posforce') fitCanvas(posforceCanvas);
  if (activeTab === 'mag')      fitCanvas(magCanvas);
  if (activeTab === 'phase')    fitCanvas(phaseCanvas);
}

// ── Full recompute + redraw ────────────────────────────────────────────────
function update() {
  if (animRaf) stopAnimation(true);
  const p = getParams();
  resp  = computeResponse(p);
  sweep = frSweep(resp.zeta);
  updateDerivedLabels();
  updateA11yDescriptions();
  resizeAll();
  drawAnimation(resp.x[0], 0);
  drawActivePlot();
}

function drawActivePlot(nowSec) {
  if (activeTab === 'pos')      drawPositionPlot(nowSec);
  if (activeTab === 'posforce') drawPositionForcePlot(nowSec);
  if (activeTab === 'mag')      drawMagPlot();
  if (activeTab === 'phase')    drawPhasePlot();
}

// ── Shared plot scaffolding helpers ────────────────────────────────────────
function clearCanvas(ctx, W, H) {
  ctx.clearRect(0, 0, W, H);
  ctx.fillStyle = COLORS.bg;
  ctx.fillRect(0, 0, W, H);
}

function drawGrid(ctx, pad, pw, ph, nx, ny) {
  ctx.strokeStyle = COLORS.grid;
  ctx.lineWidth   = 1;
  for (let i = 0; i <= nx; i++) {
    const cx = pad.l + i / nx * pw;
    ctx.beginPath(); ctx.moveTo(cx, pad.t); ctx.lineTo(cx, pad.t + ph); ctx.stroke();
  }
  for (let i = 0; i <= ny; i++) {
    const cy = pad.t + i / ny * ph;
    ctx.beginPath(); ctx.moveTo(pad.l, cy); ctx.lineTo(pad.l + pw, cy); ctx.stroke();
  }
}

function drawAxesBox(ctx, pad, pw, ph) {
  ctx.strokeStyle = '#c0c0c8';
  ctx.lineWidth   = 1;
  ctx.beginPath(); ctx.moveTo(pad.l, pad.t); ctx.lineTo(pad.l, pad.t + ph); ctx.stroke();
  ctx.beginPath(); ctx.moveTo(pad.l, pad.t + ph); ctx.lineTo(pad.l + pw, pad.t + ph); ctx.stroke();
}

// ── Position vs Time plot ──────────────────────────────────────────────────
function drawPositionPlot(nowSec) {
  if (!fitCanvas(posCanvas) && (posCanvas.width === 0)) return;
  const ctx = posCtx;
  const W   = posCanvas.width;
  const H   = posCanvas.height;
  const dpr = devicePixelRatio;
  clearCanvas(ctx, W, H);
  if (!resp) return;

  const pad = { t: 18 * dpr, r: 18 * dpr, b: 56 * dpr, l: 58 * dpr };
  const pw  = W - pad.l - pad.r;
  const ph  = H - pad.t - pad.b;

  let yMax = maxAbs(resp.x);
  yMax = yMax < 1e-9 ? 1 : yMax * 1.2;
  const tMax = resp.t[resp.t.length - 1];

  const toX = t  => pad.l + (t / tMax) * pw;
  const toY = xv => pad.t + (1 - (xv + yMax) / (2 * yMax)) * ph;

  drawGrid(ctx, pad, pw, ph, 5, 4);

  // Zero line
  ctx.strokeStyle = COLORS.zeroline;
  ctx.lineWidth   = 1.5;
  const cy0 = toY(0);
  ctx.beginPath(); ctx.moveTo(pad.l, cy0); ctx.lineTo(pad.l + pw, cy0); ctx.stroke();

  // Position curve
  ctx.strokeStyle = COLORS.r1Dark;
  ctx.lineWidth   = 1.6 * dpr;
  ctx.beginPath();
  for (let i = 0; i < resp.t.length; i++) {
    const cx = toX(resp.t[i]);
    const cy = toY(resp.x[i]);
    i === 0 ? ctx.moveTo(cx, cy) : ctx.lineTo(cx, cy);
  }
  ctx.stroke();

  // Now cursor
  if (nowSec !== undefined && nowSec >= 0) {
    const cx = toX(nowSec);
    ctx.strokeStyle = 'rgba(0,39,76,0.3)';
    ctx.lineWidth   = 1 * dpr;
    ctx.setLineDash([3 * dpr, 3 * dpr]);
    ctx.beginPath(); ctx.moveTo(cx, pad.t); ctx.lineTo(cx, pad.t + ph); ctx.stroke();
    ctx.setLineDash([]);
  }

  drawAxesBox(ctx, pad, pw, ph);

  // Tick labels
  ctx.fillStyle    = COLORS.textDim;
  ctx.font         = `${10.5 * dpr}px -apple-system, Helvetica, sans-serif`;
  ctx.textAlign    = 'center';
  ctx.textBaseline = 'top';
  for (let i = 0; i <= 5; i++) {
    const tv = (i / 5) * tMax;
    ctx.fillText(tv.toFixed(1), pad.l + i / 5 * pw, pad.t + ph + 5 * dpr);
  }
  ctx.textAlign    = 'right';
  ctx.textBaseline = 'middle';
  [-yMax, -yMax / 2, 0, yMax / 2, yMax].forEach(v => {
    ctx.fillText(v.toFixed(2), pad.l - 5 * dpr, toY(v));
  });

  // Axis labels
  ctx.fillStyle    = COLORS.textDim;
  ctx.font         = `${11 * dpr}px -apple-system, Helvetica, sans-serif`;
  ctx.textAlign    = 'center';
  ctx.textBaseline = 'top';
  ctx.fillText('Time (s)', pad.l + pw / 2, pad.t + ph + 22 * dpr);
  ctx.save();
  ctx.translate(14 * dpr, pad.t + ph / 2);
  ctx.rotate(-Math.PI / 2);
  ctx.textBaseline = 'top';
  ctx.fillText('Position (m)', 0, 0);
  ctx.restore();

  // Title-line: wn, zeta
  ctx.fillStyle = COLORS.text;
  ctx.font      = `${10.5 * dpr}px -apple-system, Helvetica, sans-serif`;
  ctx.textAlign = 'left';
  ctx.textBaseline = 'top';
  ctx.fillText(`ωₙ = ${resp.wn.toFixed(2)} rad/s    ζ = ${resp.zeta.toFixed(3)}`, pad.l + 4 * dpr, pad.t + 4 * dpr);
}

// ── Position & Force (dual y-axis) plot ────────────────────────────────────
function drawPositionForcePlot(nowSec) {
  if (!fitCanvas(posforceCanvas) && (posforceCanvas.width === 0)) return;
  const ctx = posforceCtx;
  const W   = posforceCanvas.width;
  const H   = posforceCanvas.height;
  const dpr = devicePixelRatio;
  clearCanvas(ctx, W, H);
  if (!resp) return;

  const pad = { t: 18 * dpr, r: 64 * dpr, b: 56 * dpr, l: 64 * dpr };
  const pw  = W - pad.l - pad.r;
  const ph  = H - pad.t - pad.b;

  let xMax = maxAbs(resp.x); xMax = xMax < 1e-9 ? 1 : xMax * 1.2;
  let fMax = maxAbs(resp.F); fMax = fMax < 1e-9 ? 1 : fMax * 1.2;
  const tMax = resp.t[resp.t.length - 1];

  const toX  = t  => pad.l + (t / tMax) * pw;
  const toYx = xv => pad.t + (1 - (xv + xMax) / (2 * xMax)) * ph;
  const toYf = fv => pad.t + (1 - (fv + fMax) / (2 * fMax)) * ph;

  drawGrid(ctx, pad, pw, ph, 5, 4);

  // Zero line
  ctx.strokeStyle = COLORS.zeroline;
  ctx.lineWidth   = 1.5;
  ctx.beginPath(); ctx.moveTo(pad.l, toYx(0)); ctx.lineTo(pad.l + pw, toYx(0)); ctx.stroke();

  // Force curve (right y-axis, blue)
  ctx.strokeStyle = COLORS.r2;
  ctx.lineWidth   = 1.4 * dpr;
  ctx.setLineDash([6 * dpr, 4 * dpr]);
  ctx.beginPath();
  for (let i = 0; i < resp.t.length; i++) {
    const cx = toX(resp.t[i]);
    const cy = toYf(resp.F[i]);
    i === 0 ? ctx.moveTo(cx, cy) : ctx.lineTo(cx, cy);
  }
  ctx.stroke();
  ctx.setLineDash([]);

  // Position curve (left y-axis, amber)
  ctx.strokeStyle = COLORS.r1Dark;
  ctx.lineWidth   = 1.8 * dpr;
  ctx.beginPath();
  for (let i = 0; i < resp.t.length; i++) {
    const cx = toX(resp.t[i]);
    const cy = toYx(resp.x[i]);
    i === 0 ? ctx.moveTo(cx, cy) : ctx.lineTo(cx, cy);
  }
  ctx.stroke();

  // Now cursor
  if (nowSec !== undefined && nowSec >= 0) {
    const cx = toX(nowSec);
    ctx.strokeStyle = 'rgba(0,39,76,0.3)';
    ctx.lineWidth   = 1 * dpr;
    ctx.setLineDash([3 * dpr, 3 * dpr]);
    ctx.beginPath(); ctx.moveTo(cx, pad.t); ctx.lineTo(cx, pad.t + ph); ctx.stroke();
    ctx.setLineDash([]);
  }

  drawAxesBox(ctx, pad, pw, ph);
  // Right axis line
  ctx.strokeStyle = '#c0c0c8';
  ctx.lineWidth   = 1;
  ctx.beginPath(); ctx.moveTo(pad.l + pw, pad.t); ctx.lineTo(pad.l + pw, pad.t + ph); ctx.stroke();

  // X tick labels
  ctx.fillStyle    = COLORS.textDim;
  ctx.font         = `${10.5 * dpr}px -apple-system, Helvetica, sans-serif`;
  ctx.textAlign    = 'center';
  ctx.textBaseline = 'top';
  for (let i = 0; i <= 5; i++) {
    const tv = (i / 5) * tMax;
    ctx.fillText(tv.toFixed(1), pad.l + i / 5 * pw, pad.t + ph + 5 * dpr);
  }
  // Left y ticks (position, amber)
  ctx.textBaseline = 'middle';
  ctx.textAlign    = 'right';
  ctx.fillStyle    = COLORS.r1Dark;
  [-xMax, -xMax / 2, 0, xMax / 2, xMax].forEach(v => {
    ctx.fillText(v.toFixed(2), pad.l - 5 * dpr, toYx(v));
  });
  // Right y ticks (force, blue)
  ctx.textAlign = 'left';
  ctx.fillStyle = COLORS.r2;
  [-fMax, -fMax / 2, 0, fMax / 2, fMax].forEach(v => {
    ctx.fillText(v.toFixed(1), pad.l + pw + 5 * dpr, toYf(v));
  });

  // Axis labels
  ctx.fillStyle    = COLORS.textDim;
  ctx.font         = `${11 * dpr}px -apple-system, Helvetica, sans-serif`;
  ctx.textAlign    = 'center';
  ctx.textBaseline = 'top';
  ctx.fillText('Time (s)', pad.l + pw / 2, pad.t + ph + 22 * dpr);

  // Left axis label (amber)
  ctx.save();
  ctx.translate(14 * dpr, pad.t + ph / 2);
  ctx.rotate(-Math.PI / 2);
  ctx.fillStyle = COLORS.r1Dark;
  ctx.textBaseline = 'top';
  ctx.fillText('Position (m)', 0, 0);
  ctx.restore();
  // Right axis label (blue)
  ctx.save();
  ctx.translate(W - 14 * dpr, pad.t + ph / 2);
  ctx.rotate(Math.PI / 2);
  ctx.fillStyle = COLORS.r2;
  ctx.textBaseline = 'top';
  ctx.fillText('Force (N)', 0, 0);
  ctx.restore();

  // Legend
  const legY = pad.t + 4 * dpr;
  ctx.font    = `${10 * dpr}px -apple-system, Helvetica, sans-serif`;
  ctx.textBaseline = 'middle';
  ctx.textAlign    = 'left';
  ctx.fillStyle    = COLORS.r1Dark;
  ctx.fillRect(pad.l + 4 * dpr, legY + 4 * dpr, 14 * dpr, 3 * dpr);
  ctx.fillStyle    = COLORS.text;
  ctx.fillText('Position', pad.l + 22 * dpr, legY + 5.5 * dpr);
  ctx.fillStyle    = COLORS.r2;
  ctx.fillRect(pad.l + 80 * dpr, legY + 4 * dpr, 14 * dpr, 3 * dpr);
  ctx.fillStyle    = COLORS.text;
  ctx.fillText('Force',    pad.l + 98 * dpr, legY + 5.5 * dpr);
}

// ── Magnitude FR plot ──────────────────────────────────────────────────────
function drawMagPlot() {
  if (!fitCanvas(magCanvas) && (magCanvas.width === 0)) return;
  const ctx = magCtx;
  const W   = magCanvas.width;
  const H   = magCanvas.height;
  const dpr = devicePixelRatio;
  clearCanvas(ctx, W, H);
  if (!resp || !sweep) return;

  const pad = { t: 18 * dpr, r: 18 * dpr, b: 56 * dpr, l: 58 * dpr };
  const pw  = W - pad.l - pad.r;
  const ph  = H - pad.t - pad.b;

  const rmax = 2.5;
  const Mmax = 10;

  const toX = r => pad.l + (r / rmax) * pw;
  const toY = m => pad.t + (1 - Math.min(m, Mmax) / Mmax) * ph;

  drawGrid(ctx, pad, pw, ph, 5, 5);

  // r=1 dashed reference line
  ctx.strokeStyle = '#888';
  ctx.lineWidth   = 1 * dpr;
  ctx.setLineDash([5 * dpr, 4 * dpr]);
  ctx.beginPath();
  ctx.moveTo(toX(1), pad.t); ctx.lineTo(toX(1), pad.t + ph);
  ctx.stroke();
  ctx.setLineDash([]);

  // Magnitude curve
  ctx.strokeStyle = COLORS.r1Dark;
  ctx.lineWidth   = 1.8 * dpr;
  ctx.beginPath();
  for (let i = 0; i < sweep.rvect.length; i++) {
    const cx = toX(sweep.rvect[i]);
    const cy = toY(sweep.M[i]);
    i === 0 ? ctx.moveTo(cx, cy) : ctx.lineTo(cx, cy);
  }
  ctx.stroke();

  // Operating-point marker (clamp to plot range)
  if (resp.r >= 0 && resp.r <= rmax) {
    const px = toX(resp.r);
    const py = toY(resp.M);
    ctx.fillStyle   = COLORS.r2;
    ctx.beginPath();
    ctx.arc(px, py, 6 * dpr, 0, 2 * Math.PI);
    ctx.fill();
    ctx.strokeStyle = '#fff';
    ctx.lineWidth   = 1.5 * dpr;
    ctx.stroke();
  }

  drawAxesBox(ctx, pad, pw, ph);

  // Tick labels
  ctx.fillStyle    = COLORS.textDim;
  ctx.font         = `${10.5 * dpr}px -apple-system, Helvetica, sans-serif`;
  ctx.textAlign    = 'center';
  ctx.textBaseline = 'top';
  [0, 0.5, 1, 1.5, 2, 2.5].forEach(r => {
    ctx.fillText(r.toFixed(1), toX(r), pad.t + ph + 5 * dpr);
  });
  ctx.textAlign    = 'right';
  ctx.textBaseline = 'middle';
  [0, 2, 4, 6, 8, 10].forEach(m => {
    ctx.fillText(m.toFixed(0), pad.l - 5 * dpr, toY(m));
  });

  // Axis labels
  ctx.fillStyle    = COLORS.textDim;
  ctx.font         = `${11 * dpr}px -apple-system, Helvetica, sans-serif`;
  ctx.textAlign    = 'center';
  ctx.textBaseline = 'top';
  ctx.fillText('Frequency Ratio  r = ω/ωₙ', pad.l + pw / 2, pad.t + ph + 22 * dpr);
  ctx.save();
  ctx.translate(14 * dpr, pad.t + ph / 2);
  ctx.rotate(-Math.PI / 2);
  ctx.textBaseline = 'top';
  ctx.fillText('Magnification Factor  M', 0, 0);
  ctx.restore();

  // Title-line readout
  ctx.fillStyle = COLORS.text;
  ctx.font      = `${10.5 * dpr}px -apple-system, Helvetica, sans-serif`;
  ctx.textAlign = 'left';
  ctx.textBaseline = 'top';
  ctx.fillText(`M = ${resp.M.toFixed(2)}    r = ${resp.r.toFixed(2)}`, pad.l + 4 * dpr, pad.t + 4 * dpr);
}

// ── Phase FR plot ──────────────────────────────────────────────────────────
function drawPhasePlot() {
  if (!fitCanvas(phaseCanvas) && (phaseCanvas.width === 0)) return;
  const ctx = phaseCtx;
  const W   = phaseCanvas.width;
  const H   = phaseCanvas.height;
  const dpr = devicePixelRatio;
  clearCanvas(ctx, W, H);
  if (!resp || !sweep) return;

  const pad = { t: 18 * dpr, r: 18 * dpr, b: 56 * dpr, l: 58 * dpr };
  const pw  = W - pad.l - pad.r;
  const ph  = H - pad.t - pad.b;

  const rmax = 2.5;
  // y range: -180 to 0 deg
  const toX = r => pad.l + (r / rmax) * pw;
  const toY = ph_deg => pad.t + (-ph_deg / 180) * ph; // -180→bottom, 0→top

  drawGrid(ctx, pad, pw, ph, 5, 4);

  // r=1 dashed
  ctx.strokeStyle = '#888';
  ctx.lineWidth   = 1 * dpr;
  ctx.setLineDash([5 * dpr, 4 * dpr]);
  ctx.beginPath();
  ctx.moveTo(toX(1), pad.t); ctx.lineTo(toX(1), pad.t + ph);
  ctx.stroke();
  ctx.setLineDash([]);

  // Phase curve
  ctx.strokeStyle = COLORS.r1Dark;
  ctx.lineWidth   = 1.8 * dpr;
  ctx.beginPath();
  for (let i = 0; i < sweep.rvect.length; i++) {
    const cx = toX(sweep.rvect[i]);
    const cy = toY(sweep.phi[i]);
    i === 0 ? ctx.moveTo(cx, cy) : ctx.lineTo(cx, cy);
  }
  ctx.stroke();

  // Operating-point marker
  if (resp.r >= 0 && resp.r <= rmax) {
    const px = toX(resp.r);
    const py = toY(Math.max(-180, Math.min(0, resp.phi)));
    ctx.fillStyle   = COLORS.r2;
    ctx.beginPath();
    ctx.arc(px, py, 6 * dpr, 0, 2 * Math.PI);
    ctx.fill();
    ctx.strokeStyle = '#fff';
    ctx.lineWidth   = 1.5 * dpr;
    ctx.stroke();
  }

  drawAxesBox(ctx, pad, pw, ph);

  // Tick labels
  ctx.fillStyle    = COLORS.textDim;
  ctx.font         = `${10.5 * dpr}px -apple-system, Helvetica, sans-serif`;
  ctx.textAlign    = 'center';
  ctx.textBaseline = 'top';
  [0, 0.5, 1, 1.5, 2, 2.5].forEach(r => {
    ctx.fillText(r.toFixed(1), toX(r), pad.t + ph + 5 * dpr);
  });
  ctx.textAlign    = 'right';
  ctx.textBaseline = 'middle';
  [0, -45, -90, -135, -180].forEach(p => {
    ctx.fillText(p.toFixed(0), pad.l - 5 * dpr, toY(p));
  });

  // Axis labels
  ctx.fillStyle    = COLORS.textDim;
  ctx.font         = `${11 * dpr}px -apple-system, Helvetica, sans-serif`;
  ctx.textAlign    = 'center';
  ctx.textBaseline = 'top';
  ctx.fillText('Frequency Ratio  r = ω/ωₙ', pad.l + pw / 2, pad.t + ph + 22 * dpr);
  ctx.save();
  ctx.translate(14 * dpr, pad.t + ph / 2);
  ctx.rotate(-Math.PI / 2);
  ctx.textBaseline = 'top';
  ctx.fillText('Phase Angle  φ (deg)', 0, 0);
  ctx.restore();

  // Title-line readout
  ctx.fillStyle = COLORS.text;
  ctx.font      = `${10.5 * dpr}px -apple-system, Helvetica, sans-serif`;
  ctx.textAlign = 'left';
  ctx.textBaseline = 'top';
  ctx.fillText(`φ = ${resp.phi.toFixed(1)}°    r = ${resp.r.toFixed(2)}`, pad.l + 4 * dpr, pad.t + 4 * dpr);
}

// ── Animation canvas ───────────────────────────────────────────────────────
function animLayout() {
  const W   = animCanvas.width;
  const H   = animCanvas.height;
  const dpr = devicePixelRatio;
  const pad = { t: 14 * dpr, r: 14 * dpr, b: 30 * dpr, l: 20 * dpr };
  const pw  = W - pad.l - pad.r;
  const ph  = H - pad.t - pad.b;

  const wallX  = pad.l + 0.05 * pw;
  // Equilibrium will float so the mass stays visible at maximum amplitude.
  // Compute extent of motion in physical units.
  const xMax = resp ? Math.max(maxAbs(resp.x), 1) : 1;
  // Reserve 1.6× xMax of canvas room to either side of equilibrium so the mass
  // never clips the wall or right edge.
  const physRange = Math.max(2 * xMax + 4, 6);
  const physW   = pw - 0.05 * pw;
  const scale   = physW / physRange;
  const eqX     = wallX + (xMax + 2) * scale; // some space between wall and equilibrium
  const massH   = Math.min(ph * 0.38, 44 * dpr);
  const massW   = massH;
  const midY    = pad.t + ph / 2;
  const massY   = midY - massH / 2;
  const wallTop = pad.t;
  const wallW   = 8 * dpr;
  const springY = midY - massH * 0.18;
  const damperY = midY + massH * 0.18;
  const connH   = massH * 0.12;

  return { W, H, dpr, pad, wallX, eqX, scale, massH, massW, midY, massY, wallTop, wallW, springY, damperY, connH };
}

function physToCanvasX(xPhys, L) {
  return L.eqX + xPhys * L.scale;
}

function springPath(x1, y, x2, nCoils, amp) {
  const pts = [];
  const leadFrac = 0.08;
  const len = x2 - x1;
  const leadLen = len * leadFrac;
  const zigStart = x1 + leadLen;
  const zigEnd   = x2 - leadLen;
  const zigLen   = zigEnd - zigStart;
  const nZig     = nCoils * 2;
  const spacing  = zigLen / nZig;

  pts.push([x1, y]);
  pts.push([zigStart, y]);
  for (let i = 0; i < nZig; i++) {
    const xi = zigStart + (i + 0.5) * spacing;
    const yi = y + (i % 2 === 0 ? amp : -amp);
    pts.push([xi, yi]);
  }
  pts.push([zigEnd, y]);
  pts.push([x2, y]);
  return pts;
}

function drawAnimation(xPhys, nowSec) {
  fitCanvas(animCanvas);
  const ctx = animCtx;
  const L   = animLayout();
  ctx.clearRect(0, 0, L.W, L.H);
  ctx.fillStyle = COLORS.bg;
  ctx.fillRect(0, 0, L.W, L.H);

  const massX  = physToCanvasX(xPhys, L);
  const massCX = massX + L.massW / 2;

  // Equilibrium dashed line — at mass center when at rest
  ctx.strokeStyle = '#bbb';
  ctx.lineWidth   = 1 * L.dpr;
  ctx.setLineDash([4 * L.dpr, 4 * L.dpr]);
  const eqLineX = L.eqX + L.massW / 2;
  ctx.beginPath();
  ctx.moveTo(eqLineX, L.pad.t);
  ctx.lineTo(eqLineX, L.H - L.pad.b);
  ctx.stroke();
  ctx.setLineDash([]);

  // Wall
  const wallBot = L.H - L.pad.b;
  ctx.fillStyle = '#e0e0e5';
  ctx.fillRect(L.wallX - L.wallW, L.wallTop, L.wallW, wallBot - L.wallTop);
  ctx.strokeStyle = '#b0b0b8';
  ctx.lineWidth   = 1.5 * L.dpr;
  ctx.beginPath();
  ctx.moveTo(L.wallX, L.wallTop); ctx.lineTo(L.wallX, wallBot);
  ctx.stroke();
  ctx.strokeStyle = '#c8c8d0';
  ctx.lineWidth   = 0.8 * L.dpr;
  const hatchStep = 8 * L.dpr;
  for (let y = L.wallTop; y < wallBot; y += hatchStep) {
    ctx.beginPath();
    ctx.moveTo(L.wallX - L.wallW, y);
    ctx.lineTo(L.wallX, y + hatchStep);
    ctx.stroke();
  }

  // Spring
  const springX1  = L.wallX;
  const springX2  = massX;
  const springLen = Math.max(springX2 - springX1, 4 * L.dpr);
  const amp       = Math.min(L.massH * 0.20, Math.max(5 * L.dpr, springLen * 0.08));
  const springPts = springPath(springX1, L.springY, springX2, 7, amp);

  ctx.strokeStyle = COLORS.spring;
  ctx.lineWidth   = 2 * L.dpr;
  ctx.beginPath();
  springPts.forEach(([sx, sy], i) => i === 0 ? ctx.moveTo(sx, sy) : ctx.lineTo(sx, sy));
  ctx.stroke();

  // Damper (piston)
  const dY      = L.damperY;
  const dX1     = L.wallX;
  const dX2     = massX;
  const dMid    = dX1 + (dX2 - dX1) * 0.6;
  const pistonH = L.connH;
  ctx.strokeStyle = COLORS.damper;
  ctx.lineWidth   = 2 * L.dpr;
  ctx.beginPath(); ctx.moveTo(dX1, dY); ctx.lineTo(dMid - pistonH * 1.2, dY); ctx.stroke();
  ctx.strokeRect(dMid - pistonH * 1.2, dY - pistonH, pistonH * 2.4, pistonH * 2);
  ctx.beginPath(); ctx.moveTo(dMid + pistonH * 1.2, dY); ctx.lineTo(dX2, dY); ctx.stroke();
  ctx.strokeStyle = '#b0b0b8';
  ctx.lineWidth   = 1 * L.dpr;
  ctx.beginPath(); ctx.moveTo(dMid - 1 * L.dpr, dY - pistonH); ctx.lineTo(dMid - 1 * L.dpr, dY + pistonH); ctx.stroke();

  // Mass block (maize)
  ctx.fillStyle   = 'rgba(255,203,5,0.22)';
  ctx.strokeStyle = COLORS.r1Dark;
  ctx.lineWidth   = 2 * L.dpr;
  ctx.fillRect(massX, L.massY, L.massW, L.massH);
  ctx.strokeRect(massX, L.massY, L.massW, L.massH);
  ctx.fillStyle    = COLORS.text;
  ctx.font         = `bold ${Math.max(10, 12 * L.dpr)}px -apple-system, Helvetica, sans-serif`;
  ctx.textAlign    = 'center';
  ctx.textBaseline = 'middle';
  ctx.fillText('m', massCX, L.midY);

  // Force arrow on the mass (scales with the active waveform at nowSec)
  const p = getParams();
  const wt = p.w * (nowSec || 0);
  const fNow = p.waveform === 'sin' ? Math.sin(wt) : Math.cos(wt);
  const arrowScale = (p.F0 > 0) ? fNow : 0;
  if (Math.abs(arrowScale) > 1e-3) {
    const arrowLenMax = L.massW * 0.9;
    const arrowLen    = arrowLenMax * arrowScale; // signed
    const arrowY      = L.midY;
    const startX      = massCX;
    const endX        = massCX + arrowLen;
    const headSize    = Math.min(Math.abs(arrowLen) * 0.4, 8 * L.dpr);
    const dir         = Math.sign(arrowLen);

    ctx.strokeStyle = COLORS.arrow;
    ctx.fillStyle   = COLORS.arrow;
    ctx.lineWidth   = 2.5 * L.dpr;
    // Shaft
    ctx.beginPath();
    ctx.moveTo(startX, arrowY);
    ctx.lineTo(endX - dir * headSize * 0.6, arrowY);
    ctx.stroke();
    // Head
    ctx.beginPath();
    ctx.moveTo(endX, arrowY);
    ctx.lineTo(endX - dir * headSize, arrowY - headSize * 0.6);
    ctx.lineTo(endX - dir * headSize, arrowY + headSize * 0.6);
    ctx.closePath();
    ctx.fill();
  }

  // Legend
  const legY  = L.H - L.pad.b / 2;
  const legFs = 9.5 * L.dpr;
  ctx.font         = `${legFs}px -apple-system, Helvetica, sans-serif`;
  ctx.textBaseline = 'middle';
  ctx.textAlign    = 'left';
  const legItems = [
    { color: COLORS.spring,  label: 'Spring' },
    { color: COLORS.damper,  label: 'Damper' },
    { color: COLORS.r1Dark,  label: 'Mass' },
    { color: COLORS.arrow,   label: 'Force' },
    { color: '#888888',      label: 'Equilibrium' },
  ];
  let lx = L.pad.l;
  legItems.forEach(item => {
    ctx.fillStyle = item.color;
    ctx.fillRect(lx, legY - 2 * L.dpr, 14 * L.dpr, 4 * L.dpr);
    ctx.fillStyle = COLORS.textDim;
    ctx.fillText(item.label, lx + 17 * L.dpr, legY);
    lx += ctx.measureText(item.label).width + 26 * L.dpr;
  });

  // x readout
  const xStr = `x = ${xPhys.toFixed(3)} m`;
  const fStr = `F = ${(p.F0 * fNow).toFixed(2)} N`;
  ctx.fillStyle    = COLORS.textDim;
  ctx.font         = `${9.5 * L.dpr}px -apple-system, Helvetica, sans-serif`;
  ctx.textAlign    = 'right';
  ctx.textBaseline = 'top';
  ctx.fillText(xStr, L.W - L.pad.r, L.pad.t);
  ctx.fillText(fStr, L.W - L.pad.r, L.pad.t + 14 * L.dpr);
}

// ── Animation playback ─────────────────────────────────────────────────────
const animBtn = document.getElementById('animate-btn');

function startAnimation() {
  if (animRaf) stopAnimation();
  if (!resp) return;
  animStartTs = null;
  animBtn.textContent = '⏹ Stop';
  animBtn.classList.add('running');
  animRaf = requestAnimationFrame(animFrame);
}

function stopAnimation(skipRedraw = false) {
  if (animRaf) { cancelAnimationFrame(animRaf); animRaf = null; }
  animBtn.textContent = '▶ Animate';
  animBtn.classList.remove('running');
  if (!skipRedraw && resp) {
    drawAnimation(resp.x[0], 0);
    drawActivePlot();
  }
}

function animFrame(ts) {
  if (!animStartTs) animStartTs = ts;
  const elap = (ts - animStartTs) / 1000;
  const tend = resp.t[resp.t.length - 1];
  if (elap >= tend) { stopAnimation(); return; }
  const i = Math.min(Math.round(elap * FPS), resp.x.length - 1);
  drawAnimation(resp.x[i], elap);
  if (activeTab === 'pos')      drawPositionPlot(elap);
  if (activeTab === 'posforce') drawPositionForcePlot(elap);
  // mag/phase plots are static during a run; no redraw needed
  animRaf = requestAnimationFrame(animFrame);
}

animBtn.addEventListener('click', () => {
  animRaf ? stopAnimation() : startAnimation();
});

document.getElementById('export-data-btn').addEventListener('click', exportDataAsCSV);

// ── Tab switching ──────────────────────────────────────────────────────────
const tabConfig = [
  { id: 'tab-pos',      panel: 'panel-pos',      key: 'pos' },
  { id: 'tab-posforce', panel: 'panel-posforce', key: 'posforce' },
  { id: 'tab-mag',      panel: 'panel-mag',      key: 'mag' },
  { id: 'tab-phase',    panel: 'panel-phase',    key: 'phase' },
];

function activateTab(key) {
  activeTab = key;
  tabConfig.forEach(tc => {
    const btn   = document.getElementById(tc.id);
    const panel = document.getElementById(tc.panel);
    const isActive = tc.key === key;
    btn.classList.toggle('active', isActive);
    btn.setAttribute('aria-selected', isActive ? 'true' : 'false');
    btn.setAttribute('tabindex',  isActive ? '0' : '-1');
    panel.classList.toggle('hidden', !isActive);
  });
  // Newly-revealed canvas needs to be sized and drawn
  resizeAll();
  drawActivePlot();
}

tabConfig.forEach((tc, idx) => {
  const btn = document.getElementById(tc.id);
  btn.addEventListener('click', () => activateTab(tc.key));
  btn.addEventListener('keydown', e => {
    if (e.key === 'ArrowLeft' || e.key === 'ArrowRight') {
      e.preventDefault();
      const dir = e.key === 'ArrowRight' ? 1 : -1;
      const nextIdx = (idx + dir + tabConfig.length) % tabConfig.length;
      const next = tabConfig[nextIdx];
      activateTab(next.key);
      document.getElementById(next.id).focus();
    } else if (e.key === 'Home') {
      e.preventDefault();
      activateTab(tabConfig[0].key);
      document.getElementById(tabConfig[0].id).focus();
    } else if (e.key === 'End') {
      e.preventDefault();
      activateTab(tabConfig[tabConfig.length - 1].key);
      document.getElementById(tabConfig[tabConfig.length - 1].id).focus();
    }
  });
});

// ── Input wiring ───────────────────────────────────────────────────────────
// ── Waveform toggle (sin / cos) ───────────────────────────────────────────
const waveformBtns = document.querySelectorAll('.seg-btn[data-waveform]');
waveformBtns.forEach(btn => {
  btn.addEventListener('click', () => {
    const next = btn.dataset.waveform;
    if (next === waveform) return;
    waveform = next;
    waveformBtns.forEach(b => {
      const isActive = b.dataset.waveform === waveform;
      b.classList.toggle('active', isActive);
      b.setAttribute('aria-checked', isActive ? 'true' : 'false');
    });
    update();
  });
});

function clampInputToBounds(el) {
  const v = parseFloat(el.value);
  if (!isFinite(v)) return;
  const min = parseFloat(el.min);
  const max = parseFloat(el.max);
  let clamped = v;
  if (isFinite(min) && clamped < min) clamped = min;
  if (isFinite(max) && clamped > max) clamped = max;
  if (clamped !== v) el.value = clamped;
}

document.querySelectorAll('input[type="number"]').forEach(el => {
  el.addEventListener('input', () => update());
  el.addEventListener('change', () => { clampInputToBounds(el); update(); });
});

// ── Resize observer ────────────────────────────────────────────────────────
const ro = new ResizeObserver(() => {
  resizeAll();
  drawAnimation(resp ? resp.x[0] : 0, 0);
  drawActivePlot();
});
[animCanvas, posCanvas, posforceCanvas, magCanvas, phaseCanvas].forEach(c => ro.observe(c));

// ── Init ───────────────────────────────────────────────────────────────────
update();
