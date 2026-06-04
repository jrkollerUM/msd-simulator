# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Running the app

Open `index.html` directly in a browser — no build step, no server, no dependencies. Everything runs from `file://`.

## Architecture

This is a single-page static app with three files and no framework:

- **`index.html`** — layout only. Two-column structure: `#sidebar` (controls) + `#vis-panel` (three stacked `<canvas>` elements with IDs `anim-canvas`, `plot-canvas`, `poles-canvas`). All input IDs follow the pattern `r1-*` / `r2-*` (e.g. `r1-m`, `r2-x0`).
- **`style.css`** — dark theme using CSS custom properties in `:root`. All sizing uses `devicePixelRatio`-scaled values handled in JS, not CSS, so canvas CSS size ≠ canvas pixel size.
- **`simulation.js`** — everything else. No modules; runs in a single script scope.

### Data flow in `simulation.js`

All state is module-level: `resp1`, `resp2` (computed response objects `{t: Float64Array, x: Float64Array, wn, zeta}`), `r2Enabled`, `showEnv`, `animWhich`, `animRaf`, `animData`.

The central function is `update()`:
1. Calls `getParams(n)` to read inputs → `computeResponse(p)` (analytical solution) → stores in `resp1`/`resp2`
2. Calls `updateDerivedLabels()`, `resizeAll()`, `drawPlot()`, `drawPoles()`, `drawAnimation()`

Every `oninput` event on any control calls `update()`. There is no debouncing.

### Physics

Uses closed-form analytical solutions (no ODE solver). Three branches keyed on `zeta` with a ±1e-6 tolerance band around 1.0:
- `zeta < 1−1e-6` → underdamped (uses damped frequency `wd`)
- `zeta > 1+1e-6` → overdamped (two real exponentials)
- else → critically damped

Time vector is pre-sampled at 60 Hz (`FPS = 60`, `DT = 1/FPS`). `computeResponse` returns `Float64Array`s of length `ceil(tend * FPS) + 1`.

### Canvas rendering

All three canvases are redrawn from scratch on every `update()` call. Each draw function (`drawAnimation`, `drawPlot`, `drawPoles`) calls `animLayout()` / computes its own layout inline, accounting for `devicePixelRatio` on every measurement. No retained drawing objects.

`fitCanvas(canvas)` syncs `canvas.width/height` to the element's CSS size × DPR — called in `resizeAll()` before drawing and also in the `ResizeObserver` callback.

### Drag interaction

`mousedown` on `anim-canvas` starts a drag only if the click hits the mass bounding box (checked via `massHitTest`). Position samples are buffered in `velBuf` (last 6 `{x, t}` entries). On `mouseup`, v₀ is derived from the first/last buffer entries: `v0 = (b.x − a.x) / dt` in physical meters/second. The result is written back to the `r1-v0` or `r2-v0` input and `update()` is called. Touch events are forwarded as synthetic mouse events.

## Key constraints

- No build tooling, no npm, no modules (`import`/`export` not used).
- All canvas measurements must be multiplied by `devicePixelRatio` — raw pixel math without DPR will look blurry on HiDPI screens.
- Physical x range is ±10 m; velocity range is ±50 m/s. These limits are enforced on drag and in the input `min`/`max` attributes.
- The `eqX` coordinate in `animLayout()` is the canvas x-pixel for physical x=0 (equilibrium). `physToCanvasX(xPhys, L)` converts from meters to canvas pixels.
