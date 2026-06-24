// Mass-Spring-Damper Rotational Unbalance Simulator
// Adapted from MSD_FrequencyResponse. Physics derived in Lesson 20:
//   (m + m₀)·ẍ + c·ẋ + k·x = m₀·l·ω²·cos(ωt)
// In standard form:
//   ẍ + 2ζωₙẋ + ωₙ²x = E·ω²·cos(ωt)
//   ωₙ = √(k/(m+m₀))   ζ = c/(2(m+m₀)ωₙ)   E = m₀·l/(m+m₀)
//   M  = |X/E| = r²/√((1−r²)² + (2ζr)²)    φ = −atan2(2ζr, 1−r²)

'use strict';

// ── Constants ──────────────────────────────────────────────────────────────
const FPS = 60;
const DT  = 1 / FPS;
const COLORS = {
  r1: '#FFCB05', r1Dark: '#a07800',
  r2: '#00274C',
  spring: '#555555', damper: '#888888',
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

function getParams() {
  return {
    m:    Math.min(100,   Math.max(0.01,  getNum('ru-m'))),
    k:    Math.min(50000, Math.max(0.1,   getNum('ru-k'))),
    c:    Math.min(10000, Math.max(0,     getNum('ru-c'))),
    m0:   Math.min(10,    Math.max(0,     getNum('ru-m0'))),
    l:    Math.min(1,     Math.max(0,     getNum('ru-l'))),
    w:    Math.min(50,    Math.max(0,     getNum('ru-w'))),
    x0:   Math.min(10,    Math.max(-10,   getNum('ru-x0'))),
    v0:   Math.min(50,    Math.max(-50,   getNum('ru-v0'))),
    tend: Math.min(30,    Math.max(0.1,   getNum('ru-tend'))),
  };
}

// ── Physics ────────────────────────────────────────────────────────────────
function systemProps(p) {
  const Mtot = p.m + p.m0;
  const wn   = Math.sqrt(p.k / Mtot);
  const zeta = p.c / (2 * Mtot * wn);
  const E    = Mtot > 0 ? (p.m0 * p.l) / Mtot : 0;
  return { wn, zeta, E };
}

// Steady-state particular solution coefficients for E·ω²·cos(ωt) on the RHS
// of the standard-form EOM (mass already absorbed). xp(t) = Xc cos + Xs sin.
function particularCoeffs(p, wn, zeta, E) {
  const { w } = p;
  const A = E * w * w; // effective acceleration amplitude
  if (A === 0) return { Xc: 0, Xs: 0 };
  const delta = Math.pow(wn * wn - w * w, 2) + Math.pow(2 * zeta * w * wn, 2);
  if (delta < 1e-12) return { Xc: 0, Xs: 0 };
  const Xc = (wn * wn - w * w)   * A / delta;
  const Xs = (2 * zeta * w * wn) * A / delta;
  return { Xc, Xs };
}

// Closed-form total response: x(t) = x_h(t) + x_p(t).
//
// Arm angle starts vertical at t=0: θ(t) = π/2 + ω·t. The horizontal
// component of the centripetal force on the block is m₀·l·ω²·cos(θ),
// which becomes F(t) = −m₀·l·ω²·sin(ω·t). At t=0 there is no horizontal
// kick so the system avoids an artificial transient.
//
// Equivalently, the standard-form forcing is −E·ω²·sin(ωt). With cos-derived
// coefficients (Xc, Xs), the particular solution becomes
//   xp(t) = −Xc·sin(ωt) + Xs·cos(ωt),
//   xp(0) = Xs,   xp'(0) = −Xc·ω.
function computeResponse(p) {
  const { wn, zeta, E } = systemProps(p);
  const { Xc, Xs }      = particularCoeffs(p, wn, zeta, E);
  const { x0, v0, tend, m0, l, w } = p;

  const xh0 = x0 - Xs;
  const vh0 = v0 + Xc * w;

  const N = Math.ceil(tend * FPS) + 1;
  const t = new Float64Array(N);
  const x = new Float64Array(N);
  const F = new Float64Array(N);

  const Famp = m0 * l * w * w; // peak forcing amplitude (newtons)

  for (let i = 0; i < N; i++) {
    const ti = i * DT;
    t[i] = ti;

    // Homogeneous part
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
    x[i] = xh + (-Xc * sw + Xs * cw);
    F[i] = -Famp * sw;
  }

  // Operating point on FR curves
  const r = wn > 0 ? w / wn : 0;
  const denom = Math.sqrt(Math.pow(1 - r * r, 2) + Math.pow(2 * zeta * r, 2));
  const M = denom > 1e-12 ? (r * r) / denom : 0;
  const phi = (Math.abs(r) < 1e-9 && zeta === 0) ? 0
            : -Math.atan2(2 * zeta * r, 1 - r * r) * 180 / Math.PI;

  return { t, x, F, wn, zeta, E, Xc, Xs, r, M, phi, Famp };
}

// FR sweep curves: M(r) = r²/√(…), phi(r) for r in [0, 2.5]
function frSweep(zeta) {
  const N = 600;
  const rmax = 2.5;
  const rvect = new Float64Array(N);
  const M     = new Float64Array(N);
  const phi   = new Float64Array(N);
  for (let i = 0; i < N; i++) {
    const r = (i / (N - 1)) * rmax;
    rvect[i] = r;
    const denom = Math.sqrt(Math.pow(1 - r * r, 2) + Math.pow(2 * zeta * r, 2));
    M[i] = denom > 1e-12 ? (r * r) / denom : 0;
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
  document.getElementById('ru-wn').textContent   = fmt(resp.wn, 3);
  document.getElementById('ru-zeta').textContent = fmt(resp.zeta, 4);
  document.getElementById('ru-E').textContent    = fmt(resp.E, 4);
  document.getElementById('ru-M').textContent    = fmt(resp.M, 3);
  document.getElementById('ru-r').textContent    = fmt(resp.r, 3);
  document.getElementById('ru-phi').textContent  = fmt(resp.phi, 1);
  const type = resp.zeta < 1 - 1e-4 ? 'underdamped'
             : resp.zeta > 1 + 1e-4 ? 'overdamped'
             : 'critical';
  document.getElementById('ru-damp-type').textContent = type;
}

// ── Accessibility: hidden text descriptions ────────────────────────────────
function updateA11yDescriptions() {
  if (!resp) return;
  const p = getParams();
  const type = resp.zeta < 1 - 1e-4 ? 'underdamped'
             : resp.zeta > 1 + 1e-4 ? 'overdamped'
             : 'critically damped';

  document.getElementById('anim-description').textContent =
    `Mass-spring-damper with rotating unbalance. Natural frequency ${fmt(resp.wn, 3)} rad/s, ` +
    `damping ratio ${fmt(resp.zeta, 4)} (${type}). Unbalance m₀ = ${fmt(p.m0, 3)} kg at radius l = ${fmt(p.l, 3)} m, ` +
    `drive frequency ${fmt(p.w, 2)} rad/s, eccentricity E = ${fmt(resp.E, 4)} m. ` +
    `Peak forcing amplitude m₀·l·ω² = ${fmt(resp.Famp, 3)} N. ` +
    `Initial position ${fmt(p.x0, 2)} m, initial velocity ${fmt(p.v0, 2)} m/s.`;

  const xMax = maxAbs(resp.x);
  document.getElementById('pos-description').textContent =
    `Position versus time. Time from 0 to ${fmt(p.tend, 1)} s. Peak displacement ${fmt(xMax, 3)} m.`;

  const fMax = maxAbs(resp.F);
  document.getElementById('posforce-description').textContent =
    `Position and unbalance forcing versus time. Position peak ${fmt(xMax, 3)} m. Forcing peak ${fmt(fMax, 3)} N.`;

  document.getElementById('mag-description').textContent =
    `Displacement-to-eccentricity magnitude. Frequency ratio r = ${fmt(resp.r, 3)}. ` +
    `|X/E| = ${fmt(resp.M, 3)}.`;

  document.getElementById('phase-description').textContent =
    `Phase response. Frequency ratio r = ${fmt(resp.r, 3)}. ` +
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
  a.download = 'msd-rotational-unbalance-data.csv';
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

  ctx.strokeStyle = COLORS.zeroline;
  ctx.lineWidth   = 1.5;
  const cy0 = toY(0);
  ctx.beginPath(); ctx.moveTo(pad.l, cy0); ctx.lineTo(pad.l + pw, cy0); ctx.stroke();

  ctx.strokeStyle = COLORS.r1Dark;
  ctx.lineWidth   = 1.6 * dpr;
  ctx.beginPath();
  for (let i = 0; i < resp.t.length; i++) {
    const cx = toX(resp.t[i]);
    const cy = toY(resp.x[i]);
    i === 0 ? ctx.moveTo(cx, cy) : ctx.lineTo(cx, cy);
  }
  ctx.stroke();

  if (nowSec !== undefined && nowSec >= 0) {
    const cx = toX(nowSec);
    ctx.strokeStyle = 'rgba(0,39,76,0.3)';
    ctx.lineWidth   = 1 * dpr;
    ctx.setLineDash([3 * dpr, 3 * dpr]);
    ctx.beginPath(); ctx.moveTo(cx, pad.t); ctx.lineTo(cx, pad.t + ph); ctx.stroke();
    ctx.setLineDash([]);
  }

  drawAxesBox(ctx, pad, pw, ph);

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

  ctx.fillStyle = COLORS.text;
  ctx.font      = `${10.5 * dpr}px -apple-system, Helvetica, sans-serif`;
  ctx.textAlign = 'left';
  ctx.textBaseline = 'top';
  ctx.fillText(`ωₙ = ${resp.wn.toFixed(2)} rad/s    ζ = ${resp.zeta.toFixed(3)}`, pad.l + 4 * dpr, pad.t + 4 * dpr);
}

// ── Position & Forcing (dual y-axis) plot ──────────────────────────────────
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

  ctx.strokeStyle = COLORS.zeroline;
  ctx.lineWidth   = 1.5;
  ctx.beginPath(); ctx.moveTo(pad.l, toYx(0)); ctx.lineTo(pad.l + pw, toYx(0)); ctx.stroke();

  // Forcing curve (right y-axis, navy dashed)
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

  if (nowSec !== undefined && nowSec >= 0) {
    const cx = toX(nowSec);
    ctx.strokeStyle = 'rgba(0,39,76,0.3)';
    ctx.lineWidth   = 1 * dpr;
    ctx.setLineDash([3 * dpr, 3 * dpr]);
    ctx.beginPath(); ctx.moveTo(cx, pad.t); ctx.lineTo(cx, pad.t + ph); ctx.stroke();
    ctx.setLineDash([]);
  }

  drawAxesBox(ctx, pad, pw, ph);
  ctx.strokeStyle = '#c0c0c8';
  ctx.lineWidth   = 1;
  ctx.beginPath(); ctx.moveTo(pad.l + pw, pad.t); ctx.lineTo(pad.l + pw, pad.t + ph); ctx.stroke();

  ctx.fillStyle    = COLORS.textDim;
  ctx.font         = `${10.5 * dpr}px -apple-system, Helvetica, sans-serif`;
  ctx.textAlign    = 'center';
  ctx.textBaseline = 'top';
  for (let i = 0; i <= 5; i++) {
    const tv = (i / 5) * tMax;
    ctx.fillText(tv.toFixed(1), pad.l + i / 5 * pw, pad.t + ph + 5 * dpr);
  }
  ctx.textBaseline = 'middle';
  ctx.textAlign    = 'right';
  ctx.fillStyle    = COLORS.r1Dark;
  [-xMax, -xMax / 2, 0, xMax / 2, xMax].forEach(v => {
    ctx.fillText(v.toFixed(2), pad.l - 5 * dpr, toYx(v));
  });
  ctx.textAlign = 'left';
  ctx.fillStyle = COLORS.r2;
  [-fMax, -fMax / 2, 0, fMax / 2, fMax].forEach(v => {
    ctx.fillText(v.toFixed(2), pad.l + pw + 5 * dpr, toYf(v));
  });

  ctx.fillStyle    = COLORS.textDim;
  ctx.font         = `${11 * dpr}px -apple-system, Helvetica, sans-serif`;
  ctx.textAlign    = 'center';
  ctx.textBaseline = 'top';
  ctx.fillText('Time (s)', pad.l + pw / 2, pad.t + ph + 22 * dpr);

  ctx.save();
  ctx.translate(14 * dpr, pad.t + ph / 2);
  ctx.rotate(-Math.PI / 2);
  ctx.fillStyle = COLORS.r1Dark;
  ctx.textBaseline = 'top';
  ctx.fillText('Position (m)', 0, 0);
  ctx.restore();
  ctx.save();
  ctx.translate(W - 14 * dpr, pad.t + ph / 2);
  ctx.rotate(Math.PI / 2);
  ctx.fillStyle = COLORS.r2;
  ctx.textBaseline = 'top';
  ctx.fillText('Forcing (N)', 0, 0);
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
  ctx.fillText('Forcing', pad.l + 98 * dpr, legY + 5.5 * dpr);
}

// ── Magnitude FR plot (|X/E|) ──────────────────────────────────────────────
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
  const Mmax = 6;

  const toX = r => pad.l + (r / rmax) * pw;
  const toY = m => pad.t + (1 - Math.min(m, Mmax) / Mmax) * ph;

  drawGrid(ctx, pad, pw, ph, 5, 6);

  // r=1 dashed reference
  ctx.strokeStyle = '#888';
  ctx.lineWidth   = 1 * dpr;
  ctx.setLineDash([5 * dpr, 4 * dpr]);
  ctx.beginPath();
  ctx.moveTo(toX(1), pad.t); ctx.lineTo(toX(1), pad.t + ph);
  ctx.stroke();
  // M=1 dashed reference (the high-r asymptote)
  ctx.beginPath();
  ctx.moveTo(pad.l, toY(1)); ctx.lineTo(pad.l + pw, toY(1));
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

  ctx.fillStyle    = COLORS.textDim;
  ctx.font         = `${10.5 * dpr}px -apple-system, Helvetica, sans-serif`;
  ctx.textAlign    = 'center';
  ctx.textBaseline = 'top';
  [0, 0.5, 1, 1.5, 2, 2.5].forEach(r => {
    ctx.fillText(r.toFixed(1), toX(r), pad.t + ph + 5 * dpr);
  });
  ctx.textAlign    = 'right';
  ctx.textBaseline = 'middle';
  [0, 1, 2, 3, 4, 5, 6].forEach(m => {
    ctx.fillText(m.toFixed(0), pad.l - 5 * dpr, toY(m));
  });

  ctx.fillStyle    = COLORS.textDim;
  ctx.font         = `${11 * dpr}px -apple-system, Helvetica, sans-serif`;
  ctx.textAlign    = 'center';
  ctx.textBaseline = 'top';
  ctx.fillText('Frequency Ratio  r = ω/ωₙ', pad.l + pw / 2, pad.t + ph + 22 * dpr);
  ctx.save();
  ctx.translate(14 * dpr, pad.t + ph / 2);
  ctx.rotate(-Math.PI / 2);
  ctx.textBaseline = 'top';
  ctx.fillText('|X / E|', 0, 0);
  ctx.restore();

  ctx.fillStyle = COLORS.text;
  ctx.font      = `${10.5 * dpr}px -apple-system, Helvetica, sans-serif`;
  ctx.textAlign = 'left';
  ctx.textBaseline = 'top';
  ctx.fillText(`|X/E| = ${resp.M.toFixed(2)}    r = ${resp.r.toFixed(2)}`, pad.l + 4 * dpr, pad.t + 4 * dpr);
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
  const toX = r => pad.l + (r / rmax) * pw;
  const toY = ph_deg => pad.t + (-ph_deg / 180) * ph;

  drawGrid(ctx, pad, pw, ph, 5, 4);

  ctx.strokeStyle = '#888';
  ctx.lineWidth   = 1 * dpr;
  ctx.setLineDash([5 * dpr, 4 * dpr]);
  ctx.beginPath();
  ctx.moveTo(toX(1), pad.t); ctx.lineTo(toX(1), pad.t + ph);
  ctx.stroke();
  ctx.setLineDash([]);

  ctx.strokeStyle = COLORS.r1Dark;
  ctx.lineWidth   = 1.8 * dpr;
  ctx.beginPath();
  for (let i = 0; i < sweep.rvect.length; i++) {
    const cx = toX(sweep.rvect[i]);
    const cy = toY(sweep.phi[i]);
    i === 0 ? ctx.moveTo(cx, cy) : ctx.lineTo(cx, cy);
  }
  ctx.stroke();

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
  const xMax = resp ? Math.max(maxAbs(resp.x), 1) : 1;
  const physRange = Math.max(2 * xMax + 4, 6);
  const physW   = pw - 0.05 * pw;
  const scale   = physW / physRange;
  const eqX     = wallX + (xMax + 2) * scale;

  // Reserve top headroom for the rotating arm (radius l·scale + point-mass disk).
  // Arm tip swings up to l physical meters above the pivot, plus a small disk radius.
  const lPhys      = resp ? Math.max(getParams().l, 0) : 0;
  const armRadiusPx = lPhys * scale;
  const diskRadiusPx = 8 * dpr; // approximate; clamp during draw too
  const headroom    = armRadiusPx + diskRadiusPx + 4 * dpr;

  // Available block height after reserving headroom; clamp to a sensible band.
  const massHmax = ph * 0.42;
  const massH    = Math.min(massHmax, Math.max(28 * dpr, ph - headroom - 24 * dpr));
  const massW    = Math.min(massH, 50 * dpr);

  // Block sits with its top high enough to leave room for the headroom above.
  const blockTop = pad.t + headroom;
  const massY    = blockTop;
  const midY     = massY + massH / 2;
  const wallTop = pad.t;
  const wallW   = 8 * dpr;
  const springY = midY - massH * 0.18;
  const damperY = midY + massH * 0.18;
  const connH   = massH * 0.12;

  return { W, H, dpr, pad, wallX, eqX, scale, massH, massW, midY, massY,
           wallTop, wallW, springY, damperY, connH,
           armRadiusPx, headroom };
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

  // Equilibrium dashed line
  ctx.strokeStyle = '#bbb';
  ctx.lineWidth   = 1 * L.dpr;
  ctx.setLineDash([4 * L.dpr, 4 * L.dpr]);
  const eqLineX = L.eqX + L.massW / 2;
  ctx.beginPath();
  ctx.moveTo(eqLineX, L.pad.t);
  ctx.lineTo(eqLineX, L.H - L.pad.b);
  ctx.stroke();
  ctx.setLineDash([]);

  // Ground rails — top and bottom of the block, indicating that motion is
  // constrained to horizontal. Diagonal hatching shows the constrained side.
  const railGap   = 4 * L.dpr;
  const railThick = 3 * L.dpr;
  const blockTopY = L.massY;
  const blockBotY = L.massY + L.massH;
  const railLeft  = L.wallX;
  const railRight = L.W - L.pad.r;

  // Top rail (constraint surface above the block)
  const topRailY = blockTopY - railGap - railThick;
  ctx.fillStyle = '#e0e0e5';
  ctx.fillRect(railLeft, topRailY, railRight - railLeft, railThick);
  ctx.strokeStyle = '#b0b0b8';
  ctx.lineWidth   = 1.2 * L.dpr;
  ctx.beginPath();
  ctx.moveTo(railLeft,  topRailY + railThick);
  ctx.lineTo(railRight, topRailY + railThick);
  ctx.stroke();
  // Hatch above (away from the block)
  ctx.strokeStyle = '#c8c8d0';
  ctx.lineWidth   = 0.8 * L.dpr;
  const railHatchStep = 8 * L.dpr;
  const hatchLen      = 6 * L.dpr;
  for (let xx = railLeft; xx < railRight; xx += railHatchStep) {
    ctx.beginPath();
    ctx.moveTo(xx, topRailY);
    ctx.lineTo(xx - hatchLen, topRailY - hatchLen);
    ctx.stroke();
  }

  // Bottom rail (constraint surface below the block)
  const botRailY = blockBotY + railGap;
  ctx.fillStyle = '#e0e0e5';
  ctx.fillRect(railLeft, botRailY, railRight - railLeft, railThick);
  ctx.strokeStyle = '#b0b0b8';
  ctx.lineWidth   = 1.2 * L.dpr;
  ctx.beginPath();
  ctx.moveTo(railLeft,  botRailY);
  ctx.lineTo(railRight, botRailY);
  ctx.stroke();
  // Hatch below
  ctx.strokeStyle = '#c8c8d0';
  ctx.lineWidth   = 0.8 * L.dpr;
  for (let xx = railLeft; xx < railRight; xx += railHatchStep) {
    ctx.beginPath();
    ctx.moveTo(xx, botRailY + railThick);
    ctx.lineTo(xx - hatchLen, botRailY + railThick + hatchLen);
    ctx.stroke();
  }

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

  // Damper
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

  // Rotating-unbalance mechanism on top of the block.
  // Start vertical at t=0 so the initial horizontal force component is zero:
  // θ(t) = π/2 + ω·t (CCW). The arm tip is straight up at t=0.
  const p = getParams();
  const theta = Math.PI / 2 + p.w * (nowSec || 0);
  // Pivot anchored at top-center of the block (just inside the top edge).
  const pivotX = massCX;
  const pivotY = L.massY + 4 * L.dpr;
  const armLenPx = p.l * L.scale;
  // Tip in canvas pixels (canvas y grows downward, so subtract sin to put +y axis "up")
  const tipX = pivotX + armLenPx * Math.cos(theta);
  const tipY = pivotY - armLenPx * Math.sin(theta);

  // Arm
  if (armLenPx > 0.5) {
    ctx.strokeStyle = COLORS.r2;
    ctx.lineWidth   = 2 * L.dpr;
    ctx.beginPath();
    ctx.moveTo(pivotX, pivotY);
    ctx.lineTo(tipX, tipY);
    ctx.stroke();
  }

  // Pivot dot
  ctx.fillStyle = COLORS.r2;
  ctx.beginPath();
  ctx.arc(pivotX, pivotY, 3 * L.dpr, 0, 2 * Math.PI);
  ctx.fill();

  // Point mass m₀ — disk radius scales modestly with m₀ (clamped)
  const m0Radius = Math.max(4 * L.dpr, Math.min(12 * L.dpr, 4 * L.dpr + 6 * L.dpr * Math.sqrt(p.m0 / 0.5)));
  ctx.fillStyle   = 'rgba(0,39,76,0.22)';
  ctx.strokeStyle = COLORS.r2;
  ctx.lineWidth   = 1.6 * L.dpr;
  ctx.beginPath();
  ctx.arc(tipX, tipY, m0Radius, 0, 2 * Math.PI);
  ctx.fill();
  ctx.stroke();

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
    { color: COLORS.r2,      label: 'Unbalance' },
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

  // Corner readout
  const thetaDeg = ((theta * 180 / Math.PI) % 360 + 360) % 360;
  const xStr = `x = ${xPhys.toFixed(3)} m`;
  const tStr = `θ = ${thetaDeg.toFixed(0)}°`;
  ctx.fillStyle    = COLORS.textDim;
  ctx.font         = `${9.5 * L.dpr}px -apple-system, Helvetica, sans-serif`;
  ctx.textAlign    = 'right';
  ctx.textBaseline = 'top';
  ctx.fillText(xStr, L.W - L.pad.r, L.pad.t);
  ctx.fillText(tStr, L.W - L.pad.r, L.pad.t + 14 * L.dpr);
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
