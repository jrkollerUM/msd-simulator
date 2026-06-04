// Mass-Spring-Damper Simulator
// Ported from MSDSim.m (MATLAB App Designer)

'use strict';

// ── Constants ──────────────────────────────────────────────────────────────
const FPS = 60;
const DT  = 1 / FPS;
const COLORS = { r1: '#FFCB05', r2: '#00274C', spring: '#555555', damper: '#888888', wall: '#aaa', eq: '#bbb', text: '#1d1d1f', textDim: '#6e6e73', grid: '#e8e8ed', zeroline: '#c8c8cd', bg: '#ffffff' };

// ── Canvas references ──────────────────────────────────────────────────────
const animCanvas  = document.getElementById('anim-canvas');
const plotCanvas  = document.getElementById('plot-canvas');
const polesCanvas = document.getElementById('poles-canvas');
const animCtx     = animCanvas.getContext('2d');
const plotCtx     = plotCanvas.getContext('2d');
const polesCtx    = polesCanvas.getContext('2d');

// ── Input references ───────────────────────────────────────────────────────
const getNum = id => parseFloat(document.getElementById(id).value) || 0;
const setVal = (id, v) => { document.getElementById(id).value = v; };

function getParams(n) {
  const p = n === 1 ? 'r1' : 'r2';
  return {
    m:    Math.max(0.0001, getNum(`${p}-m`)),
    k:    Math.max(0.0001, getNum(`${p}-k`)),
    c:    Math.max(0,      getNum(`${p}-c`)),
    x0:   getNum(`${p}-x0`),
    v0:   getNum(`${p}-v0`),
    tend: Math.max(0.1,   getNum(`${p}-tend`)),
  };
}

// ── Physics ────────────────────────────────────────────────────────────────
function systemProps(p) {
  const wn   = Math.sqrt(p.k / p.m);
  const zeta = p.c / (2 * Math.sqrt(p.k * p.m));
  return { wn, zeta };
}

function computeResponse(p) {
  const { wn, zeta } = systemProps(p);
  const { x0, v0, tend } = p;
  const N = Math.ceil(tend * FPS) + 1;
  const t = new Float64Array(N);
  const x = new Float64Array(N);

  for (let i = 0; i < N; i++) {
    const ti = i * DT;
    t[i] = ti;
    if (zeta < 1 - 1e-6) {
      // Underdamped
      const wd = wn * Math.sqrt(1 - zeta * zeta);
      x[i] = Math.exp(-zeta * wn * ti) *
             (x0 * Math.cos(wd * ti) + (zeta * wn * x0 + v0) / wd * Math.sin(wd * ti));
    } else if (zeta > 1 + 1e-6) {
      // Overdamped
      const sq  = Math.sqrt(zeta * zeta - 1);
      const r1  = wn * (-zeta + sq);
      const r2  = wn * (-zeta - sq);
      const A1  = (v0 + x0 * wn * (zeta + sq)) / (2 * wn * sq);
      const A2  = (-v0 - x0 * wn * (zeta - sq)) / (2 * wn * sq);
      x[i] = A1 * Math.exp(r1 * ti) + A2 * Math.exp(r2 * ti);
    } else {
      // Critically damped
      x[i] = (x0 + (v0 + x0 * wn) * ti) * Math.exp(-wn * ti);
    }
  }
  return { t, x, wn, zeta };
}

function computeEnvelope(p, t) {
  const { wn, zeta } = systemProps(p);
  const { x0, v0 } = p;
  if (zeta >= 1 - 1e-6) return null;
  const wd = wn * Math.sqrt(1 - zeta * zeta);
  const C  = Math.sqrt(x0 * x0 + Math.pow((zeta * wn * x0 + v0) / wd, 2));
  const env = t.map(ti => C * Math.exp(-zeta * wn * ti));
  return env;
}

function computePoles(p) {
  const disc = p.c * p.c - 4 * p.m * p.k;
  if (disc >= 0) {
    return [
      { re: (-p.c + Math.sqrt(disc)) / (2 * p.m), im: 0 },
      { re: (-p.c - Math.sqrt(disc)) / (2 * p.m), im: 0 },
    ];
  }
  const re = -p.c / (2 * p.m);
  const im = Math.sqrt(-disc) / (2 * p.m);
  return [{ re, im }, { re, im: -im }];
}

// ── State ──────────────────────────────────────────────────────────────────
let resp1 = null, resp2 = null;
let r2Enabled = false;
let showEnv   = false;
let animWhich = 1;
let animRaf   = null;
let animStartTs = null;
let animData    = null;  // {t, x} for selected response

// ── Derived label update ───────────────────────────────────────────────────
function updateDerivedLabels() {
  function fmt(n, digits) { return isFinite(n) ? n.toFixed(digits) : '—'; }
  for (const n of [1, 2]) {
    if (n === 2 && !r2Enabled) continue;
    const p = getParams(n);
    const { wn, zeta } = systemProps(p);
    const pre = n === 1 ? 'r1' : 'r2';
    document.getElementById(`${pre}-wn`).textContent   = fmt(wn, 3);
    document.getElementById(`${pre}-zeta`).textContent = fmt(zeta, 4);
    const type = zeta < 1 - 1e-4 ? 'underdamped' : zeta > 1 + 1e-4 ? 'overdamped' : 'critical';
    document.getElementById(`${pre}-damp-type`).textContent = type;
  }
}

// ── Full recompute + redraw ────────────────────────────────────────────────
function update() {
  resp1 = computeResponse(getParams(1));
  resp2 = r2Enabled ? computeResponse(getParams(2)) : null;
  updateDerivedLabels();
  resizeAll();
  drawPlot();
  drawPoles();
  drawAnimation(resp1.x[0]);
}

// ── Canvas resize helpers ──────────────────────────────────────────────────
function fitCanvas(canvas) {
  const w = canvas.offsetWidth  * devicePixelRatio | 0;
  const h = canvas.offsetHeight * devicePixelRatio | 0;
  if (canvas.width !== w || canvas.height !== h) {
    canvas.width  = w;
    canvas.height = h;
    return true;
  }
  return false;
}

function resizeAll() {
  fitCanvas(animCanvas);
  fitCanvas(plotCanvas);
  fitCanvas(polesCanvas);
}

// ── Position-vs-Time plot ──────────────────────────────────────────────────
function drawPlot(nowSec) {
  const ctx = plotCtx;
  const W   = plotCanvas.width;
  const H   = plotCanvas.height;
  const dpr = devicePixelRatio;

  ctx.clearRect(0, 0, W, H);
  ctx.fillStyle = COLORS.bg;
  ctx.fillRect(0, 0, W, H);

  const pad = { t: 18 * dpr, r: 18 * dpr, b: 52 * dpr, l: 58 * dpr };
  const pw  = W - pad.l - pad.r;
  const ph  = H - pad.t - pad.b;

  if (!resp1) return;

  // Y-axis range
  const allX = [...resp1.x];
  if (resp2) allX.push(...resp2.x);
  let yMax = Math.max(...allX.map(Math.abs));
  yMax = yMax < 1e-9 ? 1 : yMax * 1.2;

  // X-axis range
  const tMax = Math.max(resp1.t[resp1.t.length - 1], resp2 ? resp2.t[resp2.t.length - 1] : 0);

  const toCanvasX = t  => pad.l + (t  / tMax)  * pw;
  const toCanvasY = xv => pad.t + (1 - (xv + yMax) / (2 * yMax)) * ph;

  // Grid
  ctx.strokeStyle = COLORS.grid;
  ctx.lineWidth   = 1;
  const nGridX = 5, nGridY = 4;
  for (let i = 0; i <= nGridX; i++) {
    const cx = pad.l + i / nGridX * pw;
    ctx.beginPath(); ctx.moveTo(cx, pad.t); ctx.lineTo(cx, pad.t + ph); ctx.stroke();
  }
  for (let i = 0; i <= nGridY; i++) {
    const cy = pad.t + i / nGridY * ph;
    ctx.beginPath(); ctx.moveTo(pad.l, cy); ctx.lineTo(pad.l + pw, cy); ctx.stroke();
  }

  // Zero line
  ctx.strokeStyle = COLORS.zeroline;
  ctx.lineWidth   = 1.5;
  const cy0 = toCanvasY(0);
  ctx.beginPath(); ctx.moveTo(pad.l, cy0); ctx.lineTo(pad.l + pw, cy0); ctx.stroke();

  // Draw response helper
  function drawCurve(resp, color, dash = []) {
    if (!resp) return;
    ctx.strokeStyle = color;
    ctx.lineWidth   = 1.5 * dpr;
    ctx.setLineDash(dash);
    ctx.beginPath();
    for (let i = 0; i < resp.t.length; i++) {
      const cx = toCanvasX(resp.t[i]);
      const cy = toCanvasY(resp.x[i]);
      i === 0 ? ctx.moveTo(cx, cy) : ctx.lineTo(cx, cy);
    }
    ctx.stroke();
    ctx.setLineDash([]);
  }

  // Envelopes
  if (showEnv) {
    const p1 = getParams(1);
    if (p1 && systemProps(p1).zeta < 1 - 1e-4) {
      const env = computeEnvelope(p1, [...resp1.t]);
      const envResp = { t: resp1.t, x: Float64Array.from(env) };
      const envRespNeg = { t: resp1.t, x: Float64Array.from(env.map(v => -v)) };
      drawCurve(envResp,    COLORS.r1, [4 * dpr, 4 * dpr]);
      drawCurve(envRespNeg, COLORS.r1, [4 * dpr, 4 * dpr]);
    }
    if (resp2) {
      const p2 = getParams(2);
      if (p2 && systemProps(p2).zeta < 1 - 1e-4) {
        const env = computeEnvelope(p2, [...resp2.t]);
        const envResp    = { t: resp2.t, x: Float64Array.from(env) };
        const envRespNeg = { t: resp2.t, x: Float64Array.from(env.map(v => -v)) };
        drawCurve(envResp,    COLORS.r2, [4 * dpr, 4 * dpr]);
        drawCurve(envRespNeg, COLORS.r2, [4 * dpr, 4 * dpr]);
      }
    }
  }

  drawCurve(resp1, COLORS.r1);
  if (resp2) drawCurve(resp2, COLORS.r2);

  // Now cursor
  if (nowSec !== undefined && nowSec >= 0) {
    const cx = toCanvasX(nowSec);
    ctx.strokeStyle = 'rgba(0,39,76,0.25)';
    ctx.lineWidth   = 1 * dpr;
    ctx.setLineDash([3 * dpr, 3 * dpr]);
    ctx.beginPath(); ctx.moveTo(cx, pad.t); ctx.lineTo(cx, pad.t + ph); ctx.stroke();
    ctx.setLineDash([]);
  }

  // Axes
  ctx.strokeStyle = '#c0c0c8';
  ctx.lineWidth   = 1;
  ctx.beginPath(); ctx.moveTo(pad.l, pad.t); ctx.lineTo(pad.l, pad.t + ph); ctx.stroke();
  ctx.beginPath(); ctx.moveTo(pad.l, pad.t + ph); ctx.lineTo(pad.l + pw, pad.t + ph); ctx.stroke();

  // Tick labels
  ctx.fillStyle  = COLORS.textDim;
  ctx.font       = `${10.5 * dpr}px -apple-system, Helvetica, sans-serif`;
  ctx.textAlign  = 'center';
  ctx.textBaseline = 'top';
  for (let i = 0; i <= nGridX; i++) {
    const t = i / nGridX * tMax;
    ctx.fillText(t.toFixed(1), pad.l + i / nGridX * pw, pad.t + ph + 5 * dpr);
  }
  ctx.textAlign    = 'right';
  ctx.textBaseline = 'middle';
  const yTicks = [-yMax, -yMax / 2, 0, yMax / 2, yMax];
  yTicks.forEach(v => {
    ctx.fillText(v.toFixed(2), pad.l - 5 * dpr, toCanvasY(v));
  });

  // Axis labels
  ctx.fillStyle    = COLORS.textDim;
  ctx.font         = `${11 * dpr}px -apple-system, Helvetica, sans-serif`;
  ctx.textAlign    = 'center';
  ctx.textBaseline = 'bottom';
  ctx.fillText('Time (s)', pad.l + pw / 2, H - 6 * dpr);
  ctx.save();
  ctx.translate(14 * dpr, pad.t + ph / 2);
  ctx.rotate(-Math.PI / 2);
  ctx.textBaseline = 'top';
  ctx.fillText('Position (m)', 0, 0);
  ctx.restore();

  // Legend
  const items = [{ label: `R1  ωₙ=${resp1.wn.toFixed(2)} ζ=${resp1.zeta.toFixed(3)}`, color: COLORS.r1 }];
  if (resp2) items.push({ label: `R2  ωₙ=${resp2.wn.toFixed(2)} ζ=${resp2.zeta.toFixed(3)}`, color: COLORS.r2 });
  let lx = pad.l + 8 * dpr;
  ctx.font = `${10 * dpr}px -apple-system, Helvetica, sans-serif`;
  items.forEach(item => {
    ctx.fillStyle = item.color;
    ctx.fillRect(lx, pad.t + 4 * dpr, 18 * dpr, 3 * dpr);
    ctx.fillStyle = COLORS.text;
    ctx.textAlign = 'left';
    ctx.textBaseline = 'middle';
    ctx.fillText(item.label, lx + 22 * dpr, pad.t + 5.5 * dpr);
    lx += ctx.measureText(item.label).width + 36 * dpr;
  });
}

// ── Pole-Zero Diagram ──────────────────────────────────────────────────────
function drawPoles() {
  const ctx  = polesCtx;
  const W    = polesCanvas.width;
  const H    = polesCanvas.height;
  const dpr  = devicePixelRatio;

  ctx.clearRect(0, 0, W, H);
  ctx.fillStyle = COLORS.bg;
  ctx.fillRect(0, 0, W, H);
  if (!resp1) return;

  const pad = { t: 18 * dpr, r: 18 * dpr, b: 52 * dpr, l: 58 * dpr };
  const pw  = W - pad.l - pad.r;
  const ph  = H - pad.t - pad.b;

  // Gather all poles for auto-scale
  const allPoles = computePoles(getParams(1));
  if (resp2) allPoles.push(...computePoles(getParams(2)));

  const maxRe = Math.max(...allPoles.map(p => Math.abs(p.re)));
  const maxIm = Math.max(...allPoles.map(p => Math.abs(p.im)));
  const rng   = Math.max(maxRe, maxIm, 1) * 1.4;

  const toX = re => pad.l + (re + rng) / (2 * rng) * pw;
  const toY = im => pad.t + (1 - (im + rng) / (2 * rng)) * ph;

  // Grid
  ctx.strokeStyle = COLORS.grid;
  ctx.lineWidth   = 1;
  for (let i = -1; i <= 1; i++) {
    const cx = toX(i * rng / 2);
    ctx.beginPath(); ctx.moveTo(cx, pad.t); ctx.lineTo(cx, pad.t + ph); ctx.stroke();
    const cy = toY(i * rng / 2);
    ctx.beginPath(); ctx.moveTo(pad.l, cy); ctx.lineTo(pad.l + pw, cy); ctx.stroke();
  }

  // Axes
  ctx.strokeStyle = '#c0c0c8';
  ctx.lineWidth   = 1;
  const cx0 = toX(0), cy0 = toY(0);
  ctx.beginPath(); ctx.moveTo(cx0, pad.t); ctx.lineTo(cx0, pad.t + ph); ctx.stroke();
  ctx.beginPath(); ctx.moveTo(pad.l, cy0); ctx.lineTo(pad.l + pw, cy0); ctx.stroke();

  // Stability region fill
  ctx.fillStyle = 'rgba(0,39,76,0.04)';
  ctx.fillRect(pad.l, pad.t, cx0 - pad.l, ph);

  // Draw poles
  function drawPolesFor(poles, color) {
    poles.forEach(p => {
      const px = toX(p.re), py = toY(p.im);
      const s  = 6 * dpr;
      ctx.strokeStyle = color;
      ctx.lineWidth   = 2 * dpr;
      ctx.beginPath();
      ctx.moveTo(px - s, py - s); ctx.lineTo(px + s, py + s);
      ctx.moveTo(px + s, py - s); ctx.lineTo(px - s, py + s);
      ctx.stroke();
    });
  }

  drawPolesFor(computePoles(getParams(1)), COLORS.r1);
  if (resp2) drawPolesFor(computePoles(getParams(2)), COLORS.r2);

  // Tick labels
  ctx.fillStyle    = COLORS.textDim;
  ctx.font         = `${10.5 * dpr}px -apple-system, Helvetica, sans-serif`;
  ctx.textAlign    = 'center';
  ctx.textBaseline = 'top';
  [-1, -0.5, 0, 0.5, 1].forEach(f => {
    const val = f * rng;
    const cx  = toX(val);
    ctx.fillText(val.toFixed(1), cx, pad.t + ph + 5 * dpr);
  });
  ctx.textAlign    = 'right';
  ctx.textBaseline = 'middle';
  [-1, -0.5, 0.5, 1].forEach(f => {
    const val = f * rng;
    const cy  = toY(val);
    ctx.fillText(val.toFixed(1), pad.l - 5 * dpr, cy);
  });

  // Axis labels
  ctx.fillStyle    = COLORS.textDim;
  ctx.font         = `${11 * dpr}px -apple-system, Helvetica, sans-serif`;
  ctx.textAlign    = 'center';
  ctx.textBaseline = 'bottom';
  ctx.fillText('σ  (real)', pad.l + pw / 2, H - 6 * dpr);
  ctx.save();
  ctx.translate(14 * dpr, pad.t + ph / 2);
  ctx.rotate(-Math.PI / 2);
  ctx.textBaseline = 'top';
  ctx.fillText('jω  (imag)', 0, 0);
  ctx.restore();
}

// ── Animation canvas ───────────────────────────────────────────────────────
// Physical world: x ranges from -10 to 10 m; wall at left
// Canvas maps physical coords to pixels

function animLayout() {
  const W   = animCanvas.width;
  const H   = animCanvas.height;
  const dpr = devicePixelRatio;
  const pad = { t: 12 * dpr, r: 12 * dpr, b: 28 * dpr, l: 20 * dpr };
  const pw  = W - pad.l - pad.r;
  const ph  = H - pad.t - pad.b;

  // Horizontal: wall at 5% from left; equilibrium at 40% from left; right edge at 95%
  const wallX   = pad.l + 0.05 * pw;
  const eqX     = pad.l + 0.40 * pw;    // equilibrium = x=0
  const physW   = (0.95 * pw) - (wallX - pad.l);  // usable pixels for motion
  const xRange  = 12;                    // physical range: -2 to +10 m visible
  const scale   = physW / xRange;       // px per meter

  const massH   = Math.min(ph * 0.38, 40 * dpr);
  const massW   = massH;
  const midY    = pad.t + ph / 2;
  const massY   = midY - massH / 2;     // top of mass rect
  const wallTop = pad.t;
  const wallBot = pad.t + ph;
  const wallW   = 8 * dpr;

  const springY  = midY - massH * 0.18;
  const damperY  = midY + massH * 0.18;
  const connH    = massH * 0.12;

  return { W, H, dpr, pad, wallX, eqX, scale, massH, massW, midY, massY, wallTop, wallBot, wallW, springY, damperY, connH };
}

function physToCanvasX(xPhys, layout) {
  return layout.eqX + xPhys * layout.scale;
}

// MATLAB-style spring: flat lead-in, fixed N zigzag coils, flat lead-out
// nCoils must be even; the zig/zag alternates ±amp from centerline y
function springPath(x1, y, x2, nCoils, amp) {
  const pts = [];
  const leadFrac = 0.08;  // fraction of length for flat lead at each end
  const len = x2 - x1;
  const leadLen = len * leadFrac;
  const zigStart = x1 + leadLen;
  const zigEnd   = x2 - leadLen;
  const zigLen   = zigEnd - zigStart;
  const nZig     = nCoils * 2;  // number of zig/zag segments

  pts.push([x1, y]);
  pts.push([zigStart, y]);
  for (let i = 0; i <= nZig; i++) {
    const xi = zigStart + (i / nZig) * zigLen;
    const yi = y + (i % 2 === 0 ? amp : -amp);
    pts.push([xi, yi]);
  }
  pts.push([zigEnd, y]);
  pts.push([x2, y]);
  return pts;
}

function drawAnimation(xPhys, nowSec) {
  displayX = xPhys;
  const ctx = animCtx;
  const L   = animLayout();
  ctx.clearRect(0, 0, L.W, L.H);
  ctx.fillStyle = COLORS.bg;
  ctx.fillRect(0, 0, L.W, L.H);

  const massX = physToCanvasX(xPhys, L);   // left edge of mass
  const massCX = massX + L.massW / 2;

  // Equilibrium dashed line
  ctx.strokeStyle = COLORS.eq;
  ctx.lineWidth   = 1 * L.dpr;
  ctx.setLineDash([4 * L.dpr, 4 * L.dpr]);
  ctx.beginPath();
  ctx.moveTo(L.eqX, L.pad.t);
  ctx.lineTo(L.eqX, L.pad.t + (L.H - L.pad.t - L.pad.b));
  ctx.stroke();
  ctx.setLineDash([]);

  // Wall — light gray fill + darker edge + diagonal hatching
  ctx.fillStyle = '#e0e0e5';
  ctx.fillRect(L.wallX - L.wallW, L.wallTop, L.wallW, L.H - L.wallTop - L.pad.b);
  ctx.strokeStyle = '#b0b0b8';
  ctx.lineWidth   = 1.5 * L.dpr;
  ctx.beginPath();
  ctx.moveTo(L.wallX, L.wallTop); ctx.lineTo(L.wallX, L.H - L.pad.b);
  ctx.stroke();
  // Hatching lines
  ctx.strokeStyle = '#c8c8d0';
  ctx.lineWidth   = 0.8 * L.dpr;
  const hatchStep = 8 * L.dpr;
  for (let y = L.wallTop; y < L.H - L.pad.b; y += hatchStep) {
    ctx.beginPath();
    ctx.moveTo(L.wallX - L.wallW, y);
    ctx.lineTo(L.wallX, y + hatchStep);
    ctx.stroke();
  }

  // Spring — fixed 7 coils, amplitude scales slightly with compression
  const springX1  = L.wallX;
  const springX2  = massX;
  const springLen = springX2 - springX1;
  const amp       = Math.min(L.massH * 0.20, Math.max(5 * L.dpr, springLen * 0.08));
  const springPts = springPath(springX1, L.springY, springX2, 7, amp);

  ctx.strokeStyle = COLORS.spring;
  ctx.lineWidth   = 2 * L.dpr;
  ctx.beginPath();
  springPts.forEach(([sx, sy], i) => i === 0 ? ctx.moveTo(sx, sy) : ctx.lineTo(sx, sy));
  ctx.stroke();

  // Damper (piston style)
  const dY       = L.damperY;
  const dX1      = L.wallX;
  const dX2      = massX;
  const dMid     = dX1 + (dX2 - dX1) * 0.6;
  const pistonH  = L.connH;
  ctx.strokeStyle = COLORS.damper;
  ctx.lineWidth   = 2 * L.dpr;
  // Rod from wall to piston end
  ctx.beginPath(); ctx.moveTo(dX1, dY); ctx.lineTo(dMid - pistonH * 1.2, dY); ctx.stroke();
  // Piston box
  ctx.strokeRect(dMid - pistonH * 1.2, dY - pistonH, pistonH * 2.4, pistonH * 2);
  // Rod from piston to mass
  ctx.beginPath(); ctx.moveTo(dMid + pistonH * 1.2, dY); ctx.lineTo(dX2, dY); ctx.stroke();
  // Cylinder line inside piston
  ctx.strokeStyle = '#b0b0b8';
  ctx.lineWidth   = 1 * L.dpr;
  ctx.beginPath(); ctx.moveTo(dMid - 1 * L.dpr, dY - pistonH); ctx.lineTo(dMid - 1 * L.dpr, dY + pistonH); ctx.stroke();

  // Mass block
  const activeColor = animWhich === 2 ? COLORS.r2 : COLORS.r1;
  const massFill = animWhich === 2 ? 'rgba(0,39,76,0.12)' : 'rgba(255,203,5,0.22)';
  ctx.fillStyle   = massFill;
  ctx.strokeStyle = activeColor;
  ctx.lineWidth   = 2 * L.dpr;
  ctx.fillRect(massX, L.massY, L.massW, L.massH);
  ctx.strokeRect(massX, L.massY, L.massW, L.massH);
  ctx.fillStyle    = COLORS.text;
  ctx.font         = `bold ${Math.max(10, 12 * L.dpr)}px -apple-system, Helvetica, sans-serif`;
  ctx.textAlign    = 'center';
  ctx.textBaseline = 'middle';
  ctx.fillText('m', massX + L.massW / 2, L.midY);

  // Legend bar
  const legY  = L.H - L.pad.b / 2;
  const legFs = 9.5 * L.dpr;
  ctx.font         = `${legFs}px -apple-system, Helvetica, sans-serif`;
  ctx.textBaseline = 'middle';
  ctx.textAlign    = 'left';
  const items = [
    { color: COLORS.spring, label: 'Spring' },
    { color: COLORS.damper, label: 'Damper' },
    { color: activeColor,   label: 'Mass' },
    { color: '#aaaaaa',     label: 'Equilibrium' },
  ];
  let lx = L.pad.l;
  items.forEach(item => {
    ctx.fillStyle = item.color;
    ctx.fillRect(lx, legY - 2 * L.dpr, 14 * L.dpr, 4 * L.dpr);
    ctx.fillStyle = COLORS.textDim;
    ctx.fillText(item.label, lx + 17 * L.dpr, legY);
    lx += ctx.measureText(item.label).width + 26 * L.dpr;
  });

  // x label
  const xStr = `x = ${xPhys.toFixed(3)} m`;
  ctx.fillStyle    = COLORS.textDim;
  ctx.font         = `${9.5 * L.dpr}px -apple-system, Helvetica, sans-serif`;
  ctx.textAlign    = 'right';
  ctx.textBaseline = 'top';
  ctx.fillText(xStr, L.W - L.pad.r, L.pad.t);
}

// ── Animation playback ─────────────────────────────────────────────────────
const animBtn = document.getElementById('animate-btn');

function startAnimation() {
  if (animRaf) stopAnimation();
  animData  = animWhich === 2 && resp2 ? resp2 : resp1;
  if (!animData) return;
  animStartTs = null;
  animBtn.textContent = '⏹ Stop';
  animBtn.classList.add('running');
  animRaf = requestAnimationFrame(animFrame);
}

function stopAnimation(skipRedraw = false) {
  if (animRaf) { cancelAnimationFrame(animRaf); animRaf = null; }
  animBtn.textContent = '▶ Animate';
  animBtn.classList.remove('running');
  if (!skipRedraw) {
    if (animData) drawAnimation(animData.x[0]);
    drawPlot();
  }
}

function animFrame(ts) {
  if (!animStartTs) animStartTs = ts;
  const elap = (ts - animStartTs) / 1000;
  const tend  = animData.t[animData.t.length - 1];
  if (elap >= tend) { stopAnimation(); return; }
  const i = Math.min(Math.round(elap * FPS), animData.x.length - 1);
  drawAnimation(animData.x[i]);
  drawPlot(elap);
  animRaf = requestAnimationFrame(animFrame);
}

animBtn.addEventListener('click', () => {
  animRaf ? stopAnimation() : startAnimation();
});

// ── Drag to set initial conditions ────────────────────────────────────────
let dragActive  = false;
let velBuf      = [];   // [{x, t}]
const VEL_BUF_N = 6;
let displayX    = 0;    // physical x currently rendered on the animation canvas

function canvasToPhysX(canvasX) {
  const L = animLayout();
  const dpr = devicePixelRatio;
  return (canvasX * dpr - L.eqX) / L.scale;
}

function massHitTest(canvasX, canvasY) {
  const L   = animLayout();
  const dpr = devicePixelRatio;
  const cx  = canvasX * dpr, cy = canvasY * dpr;
  const xPhys = dragActive ? currentDragX : displayX;
  const massLeft = physToCanvasX(xPhys, L);
  return cx >= massLeft && cx <= massLeft + L.massW &&
         cy >= L.massY   && cy <= L.massY + L.massH;
}

let currentDragX = 0;

animCanvas.addEventListener('mousedown', e => {
  const rect = animCanvas.getBoundingClientRect();
  const cx   = e.clientX - rect.left;
  if (animRaf) stopAnimation(true);
  dragActive = true;
  currentDragX = Math.max(-10, Math.min(10, canvasToPhysX(cx)));
  velBuf = [{ x: currentDragX, t: performance.now() }];
  drawAnimation(currentDragX);
  animCanvas.style.cursor = 'grabbing';
});

window.addEventListener('mousemove', e => {
  if (!dragActive) return;
  const rect = animCanvas.getBoundingClientRect();
  const cx   = e.clientX - rect.left;
  currentDragX = Math.max(-10, Math.min(10, canvasToPhysX(cx)));
  velBuf.push({ x: currentDragX, t: performance.now() });
  if (velBuf.length > VEL_BUF_N) velBuf.shift();
  drawAnimation(currentDragX);
});

window.addEventListener('mouseup', () => {
  if (!dragActive) return;
  dragActive = false;
  animCanvas.style.cursor = '';

  const x0 = currentDragX;
  let v0 = 0;
  if (velBuf.length >= 2) {
    const a = velBuf[0], b = velBuf[velBuf.length - 1];
    const dt = (b.t - a.t) / 1000;  // ms → s
    if (dt > 0.001) {
      const L   = animLayout();
      const dpr = devicePixelRatio;
      const dxPx = (b.x - a.x) * L.scale;
      v0 = Math.max(-50, Math.min(50, dxPx / dpr / dt));
      // convert px/s back to m/s
      v0 = (b.x - a.x) / dt;
      v0 = Math.max(-50, Math.min(50, v0));
    }
  }

  const prefix = animWhich === 2 ? 'r2' : 'r1';
  setVal(`${prefix}-x0`, x0.toFixed(3));
  setVal(`${prefix}-v0`, v0.toFixed(3));
  update();
  startAnimation();
});

// Touch support
animCanvas.addEventListener('touchstart', e => {
  e.preventDefault();
  const t = e.touches[0];
  const ev = { clientX: t.clientX, clientY: t.clientY };
  animCanvas.dispatchEvent(new MouseEvent('mousedown', ev));
}, { passive: false });

animCanvas.addEventListener('touchmove', e => {
  e.preventDefault();
  const t = e.touches[0];
  window.dispatchEvent(new MouseEvent('mousemove', { clientX: t.clientX, clientY: t.clientY }));
}, { passive: false });

animCanvas.addEventListener('touchend', e => {
  e.preventDefault();
  window.dispatchEvent(new MouseEvent('mouseup'));
}, { passive: false });

// ── Event wiring ──────────────────────────────────────────────────────────
const allInputs = document.querySelectorAll('input[type="number"], input[type="checkbox"], select');
allInputs.forEach(el => el.addEventListener('input', () => { if (!dragActive) update(); }));

document.getElementById('enable-r2').addEventListener('change', e => {
  r2Enabled = e.target.checked;
  const r2Panel = document.getElementById('r2-panel');
  const awRow   = document.getElementById('animate-which-row');
  r2Panel.querySelectorAll('input').forEach(i => { i.disabled = !r2Enabled; });
  awRow.style.display = r2Enabled ? 'flex' : 'none';
  if (!r2Enabled) {
    // Mirror r1 values into r2
    ['m', 'k', 'c', 'x0', 'v0', 'tend'].forEach(f => {
      setVal(`r2-${f}`, getNum(`r1-${f}`));
    });
    document.getElementById('r2-wn').textContent   = '—';
    document.getElementById('r2-zeta').textContent = '—';
    document.getElementById('r2-damp-type').textContent = '';
  }
  update();
});

document.getElementById('show-envelope').addEventListener('change', e => {
  showEnv = e.target.checked;
  update();
});

document.getElementById('animate-which').addEventListener('change', e => {
  animWhich = parseInt(e.target.value);
  if (animRaf) stopAnimation();
  update();
});

// ── Resize observer ────────────────────────────────────────────────────────
const ro = new ResizeObserver(() => {
  resizeAll();
  drawPlot(undefined);
  drawPoles();
  const xNow = animData ? animData.x[0] : (resp1 ? resp1.x[0] : 0);
  drawAnimation(xNow);
});
[animCanvas, plotCanvas, polesCanvas].forEach(c => ro.observe(c));

// ── Init ───────────────────────────────────────────────────────────────────
update();
