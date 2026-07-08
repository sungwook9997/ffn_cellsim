# FF overnight plan — drive the adherent cell through its MOTILITY DEVELOPMENT, and VISUALIZE the development (2026-07-09)

## ✅ OUTCOME (morning summary)

- **S1 (C1) DONE + committed `14e2cc5`** — leading-edge protrusion wired into the native `gpu_force_fn`; NATIVE
  (A5000, Nc=266k) reproduces the resting checkpoint (γ=0.151) → adheres flat → native disp≠0 (the disp=0 blocker
  is fixed).
- **C2 crawl mechanism VALIDATED + committed `6a2a54e`** — the cell CRAWLS: directed motility EMERGES from
  protrusion + clutch turnover (no explicit treadmill needed), **traction-driven** (clutches-OFF audit PASS,
  ratio 1e4–8e4×), **physiological ~45–60 nm/s** at coarse resolution (4.5 µm in 100 s).
- **Native crawl SPEED = OPEN ITEM (diagnosed)** — `FF_CRAWL_DIAGNOSIS_2026-07-09.md`. Root cause: the crawl drag
  is grid-dependent (Σγ∝Nc), so native v∝1/Nc; the physical fix Σγ=6πηR destabilizes the implicit solver → a
  focused solver-side numerics task (rigid-mode regularization / inertial term) for PI.
- **Development VISUALIZED** — `ff_development_storyboard.png` (resting→adhered+protruding[native]→crawling) +
  `ff_dev_crawl_demo_morph.html` (frame-animated crawl + COM trace, browser-verified) + dev-curves + native S1
  viewer. Native full-res HTMLs kept local (>100 MB); screenshots + curves committed.
- **Not reached (stretch, as planned)**: S3 (emergent polarization), S4 (FA↔collagen). Hard rules kept throughout
  (no magic numbers; F* 7 pN as-recorded; no gate loosened; physiological baseline).

---


**Starting point.** A3+A1 landed (`1106dac`): the validated resting checkpoint cell now adheres to the substrate
(`--from-resting` → γ=0.171 mN/m @ ΔP=40 Pa → STABLE ADHERED), and the substrate/FA-maturation/Piezo modules have
FF-native tests. The cell is ADHERENT but STATIC. This plan drives it through its **developmental program** —
adhere → polarize → protrude → translocate → remodel ECM — and makes each step **visible as an interactive record
of how the cell develops**. Direction (C1 → C2 → C3 → B) is the PI-confirmed `FF_EXPANSION_PLAN_2026-07-08` order;
this is EXECUTION, run autonomously overnight (no approval-waits; resolve-by-document, honor every hard rule).

## The developmental narrative (what the cell DOES, stage by stage)

| stage | the cell's development | mechanism landed | the "it develops" signal |
|---|---|---|---|
| **S0** ✅ | settles + grips the substrate (resting dome) | A3 `--from-resting` | STABLE ADHERED, 84/84 FA bound |
| **S1** (C1) | grows a LEADING EDGE — the front cortex bulges out | protrusion in the native `gpu_force_fn` | front cap displaces outward, cell polarizes in shape |
| **S2** (C2) | CRAWLS — net COM translocation across the dish | front-rear clutch treadmill (asymmetric traction) | COM drift > 0 (audit: clutches-OFF ≈ 0) |
| **S3** (C3) | polarity becomes EMERGENT, not imposed | `polarization_activegel` cap → `phat` | myosin rear-cap forms on its own; drift follows it |
| **S4** (B) | pulls on the COLLAGEN matrix, remodels it | FA↔Mikado fiber catch-slip clutch (live ECM DOF) | fibers recruit/align toward the cell (1/r, DCM⊗Kim) |

Each stage = **code (mechanism) → verify (analytic/audit gate) → VISUALIZE (interactive HTML of the developed cell)**.
The development is real physics emerging from the fine-grained machinery — never a scripted animation.

## HOW the development is visualized (the emphasized deliverable)

The guiding idea: **development = the FRAME axis.** The crawl driver already records `frames` (many timepoints) +
the `com` trajectory; `ff_crawl_viewer` already has a frame slider. So the cell's morphological evolution becomes an
interactive time-machine. Four visualization modes, all **browser-verified** (`browser_check.py`, real WebGL — never
bytes) and **full-res** (no filament/node downsampling, per PI 2026-07-07):

1. **Frame-animated developmental viewer** (the key new artifact) — one HTML per stage where the frame slider plays
   the cell DEVELOPING over real time: the front bulging (S1), the body translocating (S2), the ECM deforming (S4).
   Add a persistent **COM trajectory trace** (the path the cell walked) so the MOVEMENT is legible, not just shape.
   Cut-away default scene (cortex front hemisphere removed → nucleus/aster/substrate/FA/ECM inside visible).
2. **Field-overlay scenes** — each viewer carries σ_vm cortical stress + areal-strain scenes (turbo colorbar,
   `ff_virial_stress`), so the reader sees WHERE the developing cell is loaded (front-protrusion stress, rear-clutch
   traction footprint, ECM tension lines) as it evolves.
3. **Developmental storyboard** — one montage viewer (or figure) placing the key morphological states in sequence
   (t0 resting → adhered → protruding → crawling → remodeling) so the whole development reads at a glance.
4. **Developmental trajectory curves** — the observables vs time, reality-band-overlaid (per the viz-integrity rule):
   COM displacement, shape polarization (front-rear aspect / silhouette elongation), traction [nN], contact area
   (A/A0 top-down silhouette, NOT basal contact), bound-clutch fraction. These are the QUANTITATIVE development;
   HTML morphology stays primary (viz-interactive-html rule), the curves are the companion.

Naming: `ff_dev_S{n}_{stage}_morph.html` per stage + a capstone `ff_development_trajectory.html`. Committed as
persistent project state alongside the run npz (data npz is gitignored; the self-contained HTML carries the record).

## Execution tracks + constraints (honest about mac vs gbook)

- **Mac (CPU, explicit overdamped path) — fully autonomous, carries the VISUALIZATION.** The explicit crawl loop
  already applies protrusion (line 385), so the developmental morphology + frame-animation can be produced at
  DEV SCALE (a few hundred cortex filaments) tonight, entirely on the mac, every viewer browser-verified. This is
  where the "how it develops" visuals come from.
- **gbook (A5000, native large-dt cupy path) — launch + monitor, native motion.** C1's actual gap is that the
  NATIVE path (`gpu_force_fn`) drops protrusion → native disp=0. The fix is written + unit-checked on the mac, but
  NATIVE verification (COM drift at Nc≈266k) runs on the A5000. ssh has no python on PATH → full interpreter path;
  monitor by LOG FILE (not pgrep); code synced by rsync to `~/ff_scratch` (NOT the repo); repo-promote = PI. If
  gbook is unreachable overnight, the dev-scale mac track still delivers the full developmental story; native is a
  morning follow-up.
- **Autonomy rules (this run).** No approval-waits — resolve-by-document-and-continue. Honor every hard rule: no
  magic numbers (each new constant KB-anchored or surfaced), no gate-loosening, physiological baseline, A/A0 =
  top-down silhouette. **Seed-ensemble BEFORE any "spontaneous/emergent" claim** (per the piece-1/2 lesson — a
  single-seed crawl is not evidence). **Adversarially verify** each stage's motility claim (clutches-OFF audit +
  cross-seed sign stability) before calling it real. F* 7-vs-30 pN stays AS-RECORDED (no retune). Any absolute-force
  magnitude stays flagged (FA patch radius `a` still PI-gated).

## Order of work tonight

1. **S1 (C1)** — protrusion into `gpu_force_fn` (the native "one-line blocker"). Mac: verify explicit-path
   developmental viewer (front protrudes over frames). gbook: launch native crawl, monitor for COM drift.
2. **S2 (C2)** — front-rear clutch treadmill. Verify by the clutches-OFF audit + a small cross-seed ensemble
   (drift sign stable → real translocation, not a single-seed artifact). Frame-animated crawl viewer + COM trace.
3. **S3 (C3)** — wire the emergent polarization cap into `phat`. Viewer shows the myosin cap forming + drift following.
4. **S4 (B)** — FA↔collagen-I fiber clutch (live ECM DOF). Viewer shows the cell remodeling the Mikado matrix.
5. **Capstone** — assemble the developmental-trajectory viewer + the storyboard + the observable-vs-time curves.

Each stage commits on `dcm/main` with its browser-verified viewer; `make kb-check` before each commit; Dev Logs +
memory updated at the end. Stages are independent enough that if one stalls (e.g. C2 seed-instability, or a gbook
outage) the run continues to the next and the stall is documented for PI — never blocked.

## Definition of done (morning deliverable for PI)

- The cell is shown DEVELOPING: at minimum S1+S2 landed (protrusion in the native path + a crawling cell with net
  COM translocation verified by the OFF-audit + cross-seed), each with a browser-verified frame-animated viewer.
- A capstone developmental-trajectory viewer + the observable curves, so PI can SEE and READ how the cell developed.
- Honest ledger: what emerged vs what is still seed-unstable / gbook-pending / PI-gated, with no over-claim.
