# ME 360 Redesign — First-Pass Plan

**Status:** Draft v0.1 — first structural pass, not yet locked.
**Scope:** Full 16-week semester, both sections, flipped delivery.

---

## 1. Design Principles

1. **Pre-work carries derivations and motivation.** Each session has a paired 15–20 min video that handles (a) the *why* and (b) the key derivation steps. A worked example occasionally gets its own dedicated pre-work video, but is the exception — most worked examples happen in class. Students arrive primed on the concept, not expected to have already practiced it.
2. **Class time is guided problem solving in small groups.** Two-hour sessions structured as: 10 min recap/prework Q&A → 70 min group problem work (2–3 scaffolded problems) → 25 min instructor-led synthesis + extension → 15 min preview of next topic. Buffer for transitions baked in.
3. **No standalone "math review" block.** ODEs, Laplace, and linearization are taught *inline* the first time they are needed — students who need a refresher use the Canvas prereq module on their own schedule.
4. **Modeling depth ≥ controls depth.** Modeling unit is compressed but not gutted — mechanical, electrical, electromechanical remain core. Compression comes from removing the redundant front-loaded math review, not the modeling itself.
5. **Controls expands to ~⅓ of the course.** From 5 sessions to 9 sessions, covering classical design (root locus + light Bode) instead of high-level PID only.
6. **Portable to both sections.** All artifacts (videos, slides, worksheets, exams) live in a shared repository. Section-specific only in delivery style, not content.

---

## 2. Pre-Course & Day 1: Prereq Diagnostic

### Canvas Pre-Course Module (open ~1 week before class)
A self-paced review module with the prereq topics. Items are short (5–10 min each), modeled on the existing ME 240 prework format (video + PDF + optional MATLAB .mlx).

**Topics (all draw on ME 240 + Math 216):**
- **A1.** Solving 1st-order linear ODEs (separable, integrating factor)
- **A2.** Solving 2nd-order linear ODEs with constant coefficients (homogeneous + particular)
- **A3.** Newton's 2nd Law and free-body diagrams (multi-body, with springs/dampers)
- **A4.** Rigid body planar dynamics — Newton-Euler ΣF=mā_c, ΣM=Iα
- **A5.** Energy methods (work-energy theorem, conservative forces)
- **A6.** Single-DOF undamped/damped vibration — recognizing ω_n, ζ, time response form
- **A7.** Complex numbers, Euler's identity (sin/cos ↔ exponentials)
- **A8.** Vector calculus refresher (only if instructor judges weak — chain rule, derivatives in polar)

### Day 1 Diagnostic (first class session)
- 30-minute proctored diagnostic at the start of session 1, individual.
- Mixed format: ~12 short multiple-choice + ~3 short-answer (write the ODE / sketch the response).
- **Not graded for points.** Students get immediate Canvas feedback mapping each missed question to a specific review-module item.
- Aggregate results inform instructor on which inline reviews to emphasize.
- Rest of Day 1 session: course intro, motivating examples (cruise control, suspension, drone hover), and "what is system dynamics" framing.

---

## 3. Weekly Schedule (16 weeks, 32 sessions)

Format: **W#.S#** = Week, Session. Each session has a 15–20 min pre-work video (`PW`) and an in-class guided problem set (`IC`).

### Unit 1 — Foundations & Mechanical Modeling (Weeks 1–4, 8 sessions)

| Session | Topic | PW Video | In-Class Focus |
|---|---|---|---|
| W1.S1 | Course intro + prereq diagnostic | (none — Canvas pre-course module instead) | Diagnostic + motivating demos + course logistics |
| W1.S2 | What's a dynamic system? ODE → response intuition | "From physics to ODEs to response" | Classify 3 systems; identify inputs/outputs/states |
| W2.S1 | Linearization about an operating point | "Linearization in 15 min" | Linearize 3 nonlinear systems (pendulum, tank, motor) |
| W2.S2 | Translational mechanical systems | "Springs, dampers, masses — modeling rules" | Derive EOM for 2-mass spring-damper chain |
| W3.S1 | Rotational mechanical systems | "Rotational analogs + gears" | Derive EOM for gear train + flywheel |
| W3.S2 | Mixed translation/rotation, free-body diagrams | "Constraints and coupling" | Derive EOM for rack-and-pinion, pulley systems |
| W4.S1 | State-space representation | "State variables and why we use them" | Convert 3 EOMs to state-space form |
| W4.S2 | Modeling synthesis + computational solution (MATLAB) | "ode45 for state-space models" | Build & simulate a 2-DOF suspension model |

### Unit 2 — Laplace Domain & System Response (Weeks 5–7, 6 sessions)

| Session | Topic | PW Video | In-Class Focus |
|---|---|---|---|
| W5.S1 | Laplace transforms (just-in-time, mechanics only) | "Laplace as a tool, not a theory" | Apply Laplace to 3 ODEs; recognize standard pairs |
| W5.S2 | Transfer functions from EOMs | "EOMs → TFs in 4 steps" | Derive TFs for the systems built in Unit 1 |
| W6.S1 | **Exam 1** (Modeling + state-space + intro Laplace) | — | — |
| W6.S2 | 1st-order response (time constant, step response) | "τ and what it tells you" | Identify τ from data; predict response shape |
| W7.S1 | 2nd-order response — underdamped, ζ, ω_n, ω_d | "ζ, ω_n, and the response anatomy" | Match step responses to (ζ, ω_n); rise/settling time formulas |
| W7.S2 | 2nd-order response — over/critical damping; pole locations | "Reading the s-plane" | Map pole locations ↔ time response shapes |

### Unit 3 — Electrical, Electromechanical, & Frequency Response (Weeks 8–10, 6 sessions)

| Session | Topic | PW Video | In-Class Focus |
|---|---|---|---|
| W8.S1 | Electrical systems — RLC modeling | "Kirchhoff + impedance shortcuts" | Derive TF for 3 RLC circuits |
| W8.S2 | Electromechanical — DC motor | "DC motor: 2 equations, 1 system" | Build DC motor TF + state-space; identify parameters |
| W9.S1 | **(Spring break, no class)** | — | — |
| W9.S2 | **(Spring break, no class)** | — | — |
| W10.S1 | Frequency response concept — sinusoidal steady state | "Why systems care about frequency" | Compute |G(jω)| at specific frequencies; sketch by hand |
| W10.S2 | Bode plots — sketching from poles/zeros (light) | "Bode sketching rules in 15 min" | Sketch Bode for 4 TFs; read margins from given plots |

### Unit 4 — Control Systems (Weeks 11–15, 9 sessions)

| Session | Topic | PW Video | In-Class Focus |
|---|---|---|---|
| W11.S1 | Closed-loop systems & block diagram algebra | "Open vs. closed loop + block algebra" | Reduce 3 block diagrams; derive closed-loop TFs |
| W11.S2 | Performance specs: SS error, transient, sensitivity | "What do we actually want from a controller?" | Compute SS error for type 0/1/2 systems |
| W12.S1 | **Exam 2** (Laplace + response + freq response intro) | — | — |
| W12.S2 | Stability — Routh-Hurwitz, characteristic equation | "Stability without solving the polynomial" | Apply Routh to 3 systems; find stable gain ranges |
| W13.S1 | PID control — structure & effect of each term | "P, I, D — what each one fixes and breaks" | Tune P/PI/PID on simulated plant; observe trade-offs |
| W13.S2 | PID tuning methods — Z-N, manual, model-based | "Three ways to tune PID" | Tune a real plant model 3 ways; compare performance |
| W14.S1 | Root locus — sketching rules | "Root locus in 15 min" | Sketch RL for 4 plants; identify breakaway, asymptotes |
| W14.S2 | Root locus design — gain selection, lead/lag intro | "Designing with root locus" | Design controllers to meet ζ/ω_n specs using RL |
| W15.S1 | Stability margins from Bode (light) — GM, PM | "Gain & phase margin without the theory dump" | Read margins from Bode plots; relate to robustness |
| W15.S2 | Practical considerations — saturation, anti-windup, noise | "What breaks PID in the real world" | Diagnose & fix 3 broken PID implementations |

### Unit 5 — Synthesis & Final (Week 16)

| Session | Topic | PW Video | In-Class Focus |
|---|---|---|---|
| W16.S1 | Capstone case study (e.g., inverted pendulum or DC motor servo) | "End-to-end design walkthrough" | Group design project — full model → controller → sim |
| W16.S2 | Review for final | "What you should know — concept map" | Mixed Q&A + group practice problems |
| Finals week | **Final Exam** (cumulative) | — | — |

---

## 4. Comparison to Current Course

| Area | Current | Redesigned | Δ |
|---|---|---|---|
| Front-loaded math (ODE/Laplace/TF intro) | 5 lessons (L1–5) | 0 standalone + just-in-time inline | **−5** sessions reclaimed |
| Mechanical modeling | 4 lessons | 4 sessions | 0 |
| Electrical + Electromechanical | 3 lessons | 2 sessions | −1 |
| State-space | 2 lessons | 1 session | −1 |
| Time-domain response | 4 lessons | 3 sessions | −1 |
| Frequency response + Bode + Sys ID | 3 lessons | 2 sessions | −1 |
| **Controls** | **5 lessons** (PID at high level) | **9 sessions** (PID + Routh + RL + margins + practical) | **+4** |
| Capstone synthesis | 1 lesson | 1 session | 0 |
| Exams | 3 | 3 | 0 |

Net: enough room reclaimed from redundant math review + slight modeling compression to nearly double controls coverage, while introducing classical design tools (Routh, root locus, margins) that prepare students for ME 461 *or* leave non-ME-461 students with real control design ability.

---

## 5. Open Questions / Decisions Needed Before v1

1. **Sys ID / filtering content** — current course has it (L19–21). I dropped it from the schedule above to make room for controls. Acceptable, or should one session be added back?
2. **Capstone case study** — inverted pendulum (matches current course final lesson) or DC motor servo (more applied / industry-relevant)? Could rotate yearly.
3. **Diagnostic timing** — administer Day 1 during class as proposed, or push it to Canvas before Day 1? In-class costs 30 min but guarantees completion + lets students struggle freshly.
4. **Group composition** — fixed groups for the whole semester, rotating monthly, or per-session shuffled? Fixed groups build psych safety but risk dysfunction; rotating spreads exposure.
5. **Pre-work accountability** — auto-graded Canvas check after each video (low-stakes, ~5%/semester)? Without some accountability, flipped delivery often collapses.
6. **MATLAB vs. Python** — current course is MATLAB throughout. Keep, or open the door to Python? Strong recommendation: stay MATLAB for consistency with ME 240.
7. **Exam style change** — flipped courses sometimes shift to take-home or two-stage exams (individual then group). Worth considering or out of scope for v1?

---

## 6. Next Steps (Suggested)

1. **Get instructor feedback on this v0.1 structure** — particularly the controls expansion shape and the dropped Sys ID block.
2. **Draft the prereq diagnostic** (15 questions, with answer-to-review-module mapping).
3. **Storyboard one full week** as a worked example (videos + in-class problems + materials) to validate the time budget.
4. **Build the Canvas review module index** mapping each diagnostic question to a specific review item.
5. **Align with the other section's instructor** before locking v1.
