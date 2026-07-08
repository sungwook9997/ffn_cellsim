# FF expansion plan — from the validated SUSPENDED cell to ADHERENT + MOTILE (2026-07-08)

**Starting point (the checkpoint):** `FF_CORTICAL_MECHANICS_STATE_2026-07-08.md` — a single MCF7 cell, full
compartment (native woven cortex + nucleus + membrane reservoir + pressure-borne turgor + biphasic drained cytoplasm
+ MT aster), validated on the CORTICAL-MECHANICS axis, but **SUSPENDED** (no substrate/adhesion) and STATIC. Viewer
record: `outputs/ff/figs/resting_full_compartment.html`.

This plan takes it to (A) an ADHERENT cell on a substrate — matching the PI's experiment (MCF7 on pV4D4 / collagen-I)
— and (C) a MOTILE cell that translocates. Each phase lists **code status** (from the 2026-07-08 machinery audit) and
**KB status** (from the 2026-07-08 TAG/DuckDB coverage audit — "is the KB information sufficient", PI's Q3).

## KB coverage summary (Q3 — is TAG info sufficient for the expansion?)

Direct DuckDB audit (162 KB-claims / 355 SourceEvidence). Per expansion topic:

| topic | KB-claims | SE | verdict |
|---|---|---|---|
| FA / integrin clutch (α5β1-FN, catch-slip, talin, vinculin, maturation) | 16 | 9 | **OK** (KB-2.1/2.2/2.4/2.5/2.6/2.7, High) |
| substrate stiffness (compliant, PAA, durotaxis, Bangasser-Odde) | 9 | 6 | **OK** (KB-1.5/1.21/1.14, High) |
| ECM / collagen fiber (Mikado) | 17 | 8 | **OK** (KB-1.x collagen mechanics) |
| traction force | 7 | 8 | **OK** (KB-2.12) |
| protrusion / lamellipodium / Arp2-3 / filopodium | 7 | 7 | **OK** (KB-3.6/3.7/3.8/3.12/3.21, H5-LP-1) |
| membrane / spreading (A/A0, reservoir) | 8 | 12 | **OK** |
| polarization / active gel | 2 | 3 | OK (thin but present) |
| **laminin (α6β1, PI Lam4)** | 1 | 0 | **THIN — but literature-limited, not a KB gap**: KB-PIV-5 documents that direct α6β1–laminin single-molecule kinetics are GENUINELY ABSENT in the literature → a grounded proxy is already registered. Nothing more to collect. |
| polymerization / barbed-end | 1 | 0 | OK (KB-3.6 ratchet; keyword under-counts) |
| **LINC / nucleus-cytoskeleton (nesprin/SUN)** | 0 | 1 | **GENUINE KB GAP** — needs new SE rows (Lombardi 2011, Kirby-Lammerding) before piece-5 |

**Q3 answer:** For the **adherent expansion the KB is SUFFICIENT** — FA/FN, substrate stiffness, ECM/collagen,
traction, protrusion, spreading are all richly covered (High confidence). Two items are NOT collectable/absent:
(a) direct α6β1-laminin kinetics (real literature gap; the PI's Lam4 path uses the registered proxy), and (b) **LINC
mechanocoupling — a genuine KB gap to fill** if the expansion includes nucleus positioning/mechanotransduction.

## Phase A — ADHERENT cell on a compliant substrate  [code: MOSTLY THERE · KB: OK]

Goal: the suspended cell settles on an elastic substrate, forms integrin catch-slip clutches, generates traction,
and spreads into a domed adherent cell — the PI's MCF7-on-substrate baseline.
- **Built (GPU-native, wired):** `fa_clutch_warp.py` (catch-slip clutch KMC, validated 3/3), `fa_maturation.py`
  (talin/vinculin/Piezo), `substrate.py` (compliant Winkler E=5 kPa, Bangasser-Odde series stiffness), spreading
  (`motility_warp` basal ratchet + area growth + volume closure). `ff_cell_on_substrate.py` already demos a settled,
  clutch-gripped, traction-generating cell.
- **Gaps to close:** (A1) **no FF-native tests** for `substrate.py` / `fa_maturation.py` / `piezo.py` — validation
  debt; add series-stiffness + rigid-limit + analytic-gate tests. (A2) **absolute magnitude KB-blocked**: the FA
  patch radius `a` (sets k_sub) is unregistered, and the α5β1-FN catch-slip peak is F*≈7 pN in code vs KB/Kong-2009
  30 pN — **surface to PI** before any traction-magnitude claim. (A3) wire the resting full-compartment builder to
  the adherent path (`--static --mature --spread`) so the CHECKPOINT cell (not a fresh one) adheres.
- **KB: sufficient** (FA/FN, substrate, traction, spreading all OK). Only the magnitude datum (`a`) + the 7-vs-30 pN
  reconciliation need a PI call.

## Phase B — Adhesion to COLLAGEN-I fibers (the PI's actual matrix)  [code: MISSING · KB: OK]

Goal: FA clutches bind MOVABLE collagen-I ECM fibers, not just a rigid plane.
- **Decisive gap:** the clutch anchor is a fixed point/plane; the Mikado collagen network (`ecm_mikado.py`) is a
  separate, unbound structure — the cell touches it only by repulsive excluded-volume. **No kernel forms a catch-slip
  clutch onto a live ECM fiber node.** So "adhere to collagen-I" is not possible today.
- **Step:** let `attach_fa_clutches` target Mikado fiber nodes and make the clutch spring anchor a live ECM DOF
  (Newton reaction onto the fiber), so traction remodels the ECM (ties into the validated DCM⊗ECM 1/r remodel work).
- **KB: OK** — collagen mechanics + integrin-FN + traction all covered; the DCM⊗Kim ECM remodel is already validated.

## Phase C — MOTILITY / translocation  [code: NOT SUFFICIENT · KB: OK except LINC]

Goal: the adherent cell crawls (net COM translocation). All 5 pieces exist as validated components, but the closed
loop is broken at integration.
- **C1 (the decisive one-line blocker):** `gpu_force_fn` — the ONLY native-scale, large-dt path — assembles
  `clutches` + `spread` but **never applies the leading-edge protrusion** (`ff_crawl_on_substrate.py:288-317` has no
  `protrude` branch). Add `leading_edge_push_kernel` + `protrusion_reaction_kernel` to it. This is the "real unlock".
- **C2:** front-rear **clutch treadmill** — polarity-dependent engage/rupture (front nascent, rear release) instead of
  uniform host-side rebind — so traction is asymmetric and converts protrusion into net drift.
- **C3:** wire **polarization (piece 2)** — feed `polarization_activegel`'s emergent cap into the polarity axis
  `phat` instead of hard-coding it (ideally an emergent 3-D myosin cap).
- **KB: OK** — traction, protrusion/lamellipodium, polarization all covered.

## Phase D — Nucleus mechanocoupling (LINC, piece 5)  [code: MISSING · KB: GAP]

Goal: the nucleus mechanically couples to the migration machinery (positioning, mechanotransduction, confined
migration). **No FF LINC module exists** (only archived HOOMD `linc.py`); the nucleus is a free bead-shell.
- **Step:** a Warp LINC coupling (nesprin/SUN springs cortex/MTOC ↔ nuclear surface). **KB gap first**: register the
  LINC SE rows (Lombardi 2011; Kirby-Lammerding) before building — this is the one genuine KB gap.

## Cross-cutting

- **Validation debt:** substrate / maturation / piezo have no FF tests. Add before any adherent result is called
  thesis-grade.
- **Magnitude grounding (PI calls):** FA patch radius `a`; α5β1-FN catch-slip F* 7 vs 30 pN.
- **Residual cortical-mechanics calibration** (from the state doc): over-stiff high-strain tail (6Z), drag ~2×
  factor, coherent↔linear baseline — refine but not blocking.

## Recommended order

1. **A3 + A1** — ✅ **DONE 2026-07-08** (PI-confirmed direction). `--from-resting` wires the validated resting
   checkpoint into the adherent driver (reproduces γ=0.171 mN/m @ ΔP=40 Pa → STABLE ADHERED); 25 FF-native tests
   added for substrate / FA-maturation / Piezo (all analytic gates PASS). F* 7-vs-30 pN: kept as-recorded per PI,
   flagged for later debugging. Adversarially reviewed (0 bugs). Details: `FF_ADHERENT_A3_A1_2026-07-08.md`.
   Still PI-gated before any absolute-traction claim: FA patch radius `a`.
   — wire the checkpoint cell to the adherent path + add the missing FF tests (fast, unblocks a
   defensible adherent baseline matching the PI's experiment).
2. **C1** — the one-line protrusion fix in `gpu_force_fn` (unblocks native-scale motion) + C2 clutch treadmill.
3. **B** — FA↔collagen-fiber binding (the PI's real matrix).
4. **D** — LINC (KB rows first, then code).
Magnitude PI-calls (A2) gate any absolute-force claim throughout.
