# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Running the app

Open `index.html` directly in a browser — no build step, no server, no dependencies. Everything runs from `file://`.

## Architecture

This is a single-page static app with three files and no framework:

- **`index.html`** — layout only. Two-column structure: `#sidebar` (controls) + `#vis-panel` (three stacked `<canvas>` elements with IDs `anim-canvas`, `plot-canvas`, `poles-canvas`). All input IDs follow the pattern `r1-*` / `r2-*` (e.g. `r1-m`, `r2-x0`).
- **`style.css`** — minimalist white/light-gray theme using CSS custom properties in `:root` (UMich colors: R1=Maize `#FFCB05`, R2=Blue `#00274C`). All sizing uses `devicePixelRatio`-scaled values handled in JS, not CSS, so canvas CSS size ≠ canvas pixel size.
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

Physical parameter limits (enforced in both HTML and `getParams()`):
- **m**: 0.01–100 kg
- **k**: 0.1–50,000 N/m  
- **c**: 0–10,000 N·s/m
- **x₀**: ±10 m
- **v₀**: ±50 m/s
- **t_end**: 0.1–20 s

### Canvas rendering

All three canvases are redrawn from scratch on every `update()` call. Each draw function (`drawAnimation`, `drawPlot`, `drawPoles`) calls `animLayout()` / computes its own layout inline, accounting for `devicePixelRatio` on every measurement. No retained drawing objects.

`fitCanvas(canvas)` syncs `canvas.width/height` to the element's CSS size × DPR — called in `resizeAll()` before drawing and also in the `ResizeObserver` callback.

### Drag interaction

`mousedown` on `anim-canvas` stops any running animation, snaps the mass to the clicked x position, and initiates a drag. Position samples are buffered in `velBuf` (last 6 `{x, t}` entries). On `mouseup`, v₀ is derived from the first/last buffer entries: `v0 = (b.x − a.x) / dt` in physical meters/second. The result is written back to the `r1-v0` or `r2-v0` input, `update()` is called, and the animation auto-starts. `displayX` tracks the live animated position for smooth hit-testing mid-animation. Touch events are forwarded as synthetic mouse events.

## UI/UX behaviors

- **Parameter changes**: Any input change stops the animation and resets the mass to the new initial conditions.
- **Response 2 enable**: Checking "Enable Response 2" mirrors all R1 parameter values into R2 as a baseline.
- **Animation auto-start**: Releasing a drag automatically starts the animation from the new initial condition.
- **Live drag during animation**: Clicking anywhere on the animation canvas snaps the mass to that position and starts a drag, even if an animation is running.
- **Canvas labels**: "Characteristic Roots" (not "Pole-Zero Diagram"), x-axis "Real", y-axis "Imaginary".

## Spring drawing

The spring uses a fixed 7-coil MATLAB-style zigzag with flat leads at each end. Nodes are offset by half-spacing at entry/exit (`i + 0.5` instead of `i / nZig`) to eliminate abrupt vertical segments and create smooth angled transitions from the horizontal leads.

## Accessibility (WCAG 2.1 AA)

The app is being brought into WCAG 2.1 AA compliance. See `WCAG_2.1_AA_Compliance_Report.txt` for a full audit.

**Canvas accessibility**: Each of the three canvas elements (`anim-canvas`, `plot-canvas`, `poles-canvas`) has:
- `aria-label` describing the canvas purpose (e.g., "Animation of mass position over time")
- A hidden `<div role="region" aria-live="polite">` below the canvas with live-updated text descriptions of current state (damping type, key values, axis ranges)
- A "Download Data" button that exports current simulation results as CSV for screen reader / data inspection use

**Keyboard accessibility**:
- `anim-canvas` has `tabindex="0"` — reachable via Tab in natural order
- Arrow keys (↑/↓) adjust velocity by 0.5 m/s when canvas has focus; Shift+arrow = 0.1 m/s fine step
- Focus ring (3px `#00274C` outline) shown only during keyboard navigation via `:focus-visible`; mouse clicks do not trigger the outline (`:focus { outline: none }` suppresses it)
- Drag interaction fully preserved; `animCanvas.focus()` was intentionally removed from the mousedown handler to prevent focus ring appearing on drag

**Canvas legends**: Both the animation canvas component legend (Spring/Damper/Mass/Equilibrium) and the position plot response legend (R1/R2 ωₙ/ζ labels) are currently commented out in `simulation.js` due to insufficient color contrast between elements. Search for `// Legend` to locate them for future re-enabling with improved contrast.

**Related files**: `WCAG_2.1_AA_Compliance_Report.txt` documents issues and remediation priority.

## Key constraints

- No build tooling, no npm, no modules (`import`/`export` not used).
- All canvas measurements must be multiplied by `devicePixelRatio` — raw pixel math without DPR will look blurry on HiDPI screens.
- Physical x range is ±10 m; velocity range is ±50 m/s.
- The `eqX` coordinate in `animLayout()` is the canvas x-pixel for physical x=0 (equilibrium). `physToCanvasX(xPhys, L)` converts from meters to canvas pixels.
- Axis labels on time plot and characteristic roots plot are positioned 22 logical pixels below the plot baseline (after tick labels at +5*dpr); bottom padding is 72*dpr to ensure full visibility at all DPR values.
