# Compartment Activation — Master Plan (2026-06-09)

The strategic roadmap from the **current state** (a validated, default-OFF
compartment platform) to a **working integrated physiological MCF7 cell** and
then a **multicell doublet**. This is the 전체 계획; the per-compartment 세부 계획
live in `COMPARTMENT_ACTIVATION_DETAILED_PLANS_2026-06-09.md`.

Governing contracts: `AGENTS.md`/`CLAUDE.md` hard rules, `H7_GATE_B_CONTRACT`,
`H7_SURFACE_MANIFOLD_EXPLICIT_CORTEX`, `COMPARTMENT_PLATFORM_PLAN_2026-06-08`,
`COMPARTMENT_PHYSICS_AUDIT_2026-06-09`.

---

## 0. Where we are (2026-06-09)

| Tier | Compartments | State |
|---|---|---|
| **Core** (always assembled) | cortex, crosslinkers, myosin | LIVE |
| **Baseline** (ON at setpoint) | cytoplasm (65.9 Pa·s), enclosed_volume (turgor 133 Pa), nucleus, membrane_surface | LIVE |
| **Optional, wired** (PI-gated activation) | fa, rigid_ligand_coating, substrate, lamellipodium, membrane_load, erm, turnover | LIVE |
| **Geometry** (force-free, not wired) | surface_manifold | EXPERIMENTAL |
| **Missing, authored** (default-OFF, NOT wired, raise-on-enable) | ventral_stress_fibers, linc, intermediate_filaments, microtubules, osmotic_regulation, membrane_reservoir, cadherin_junction, junctional_actin | EXPERIMENTAL/STUB |

The registry/recipe spine, 8 compartment modules, and a 2-pass adversarial physics
audit are landed and Lead-verified (commits `eee86a0`, `2928252`, `7990614`). No
compartment in the bottom three rows is enabled in production.

## 1. The two frontiers (do not conflate them)

The work splits into two orthogonal frontiers. The compartment platform is the
**breadth** frontier; it does NOT fix the **depth** frontier.

**A. DEPTH — cortical-γ closure (cortex-internal, pre-existing workstream).**
Gate-A established that resting active γ_soft ≈ 3.06e-3 mN/m sits ~114× below the
MCF7 datum (0.27 mN/m): the wall is **transmission/lever**, not generation
(myosin binds, loads, walks, steps — but contraction does not transmit into
spanning hoop tension on a non-condensing M-SHAKE backbone). **Gate-B** (the
compression-side buckling/condensation lever, band Hosseini [0.18, 0.40] mN/m) is
the open frontier — this is the in-flight WIP on this branch. **No compartment
below fixes this**; it is a property of the cortex backbone + integrator.

**B. BREADTH — the compartment platform (this plan).** Adds the rest of the
mechanobiology as new, separately-gated observables (traction, spread area,
junction tension, nuclear strain, volume regulation, blebbing). These are mostly
**orthogonal** to the cortical-γ number. The one coupling to watch: stress fibers
and junctional actin also generate/transmit tension, so their bond prefixes must
be denylisted from the cortical-γ estimator (no contamination).

→ The two frontiers proceed **in parallel**. Closing Gate-B is not a prerequisite
for activating, say, the FA-traction or cadherin-junction observables.

## 2. The universal activation procedure (every compartment)

A compartment graduates `EXPERIMENTAL → LIVE` through the same five-step gate.
This is the reusable contract; the detailed plans fill in the per-compartment
specifics.

1. **Parameter ratification (PI).** Every `None`/flagged constant gets a PI
   decision. The detailed plans propose literature-anchored candidates (with a
   SOLID / ORDER-ESTIMATE / NEEDS-MEASUREMENT confidence tag). PI ratifies, rejects,
   or commissions a measurement. No constant goes LIVE as an invented number.
2. **Wiring.** Add (a) a `configs/mcf7_baseline.yaml` `optional_subsystems` block,
   (b) the `resolve_*` call in `manifest.resolve_baseline`, (c) the
   `_extend_snapshot_*` / `attach_*` calls in `cell.py` (single-writer), (d) a
   `CellBuildOptions` flag, (e) flip the registry `CompartmentSpec.status` to LIVE
   and add `manifest_path`. The disabled path stays bit-identical.
3. **No-γ-contamination control (HARD).** Add the compartment's bond prefix to
   `cortex/cortical_tension.py` `ADHESION_BOND_TYPE_PREFIXES` (or a new non-cortical
   family), and add a test asserting the suspended cortical-γ is **unchanged**
   (<1%) ON vs OFF with myosin OFF + turgor OFF (the manifold's zero-tension-injection
   gate, generalized). Stress fibers additionally need the `sf_myosin_*` prefix split.
4. **Pairwise physiological gate.** Run the compartment ON in the FULL physiological
   baseline (cytoplasm + turgor + nucleus + membrane ON), measure its OBSERVABLE,
   compare to the literature BAND with mandatory CONTROLS (OFF/ON differencing,
   rigid/limit parity, the no-contamination check), and a written PASS / PARTIAL /
   REFUTE threshold — authored BEFORE the run, not loosened.
5. **GPU/CFL + recipe-level gate.** Honour the hot-path priority (P1 CPU debt noted +
   microbenched; native ForceCompute for the simple springs). Confirm the dt/CFL
   impact (critical for microtubules). Once the pairwise gate passes, the owning
   recipe (e.g. `adherent_passive`) gets a recipe-level observable gate.

## 3. Dependency DAG (what must precede what)

```
cortex/myosin/xlink (core) ─┬─ enclosed_volume ── osmotic_regulation
                            ├─ membrane_surface ─┬─ membrane_reservoir
                            │                    └─ erm (couples)
                            ├─ nucleus ── linc ──(couples to cortex / SF / MT)
                            ├─ fa ─┬─ rigid_ligand_coating (LIVE)
                            │      ├─ substrate (soft-gel variant)
                            │      └─ ventral_stress_fibers
                            ├─ lamellipodium ── membrane_load
                            ├─ intermediate_filaments   (no hard dep)
                            ├─ microtubules             (no hard dep; CFL risk)
                            └─ cadherin_junction ── junctional_actin   (2-cell)
```

## 4. Phased activation sequence

Ordered by (a) dependency, (b) parameter readiness, (c) cost/risk. Each phase is a
PI-gated checkpoint; do not batch-enable.

**Phase A — Adherent operating point (LIVE compartments).** No new physics; close
out what is already wired. (i) FA + rigid_ligand_coating → `adherent_passive`
recipe, with the equilibration prelude, gated on a **traction** observable. (ii)
lamellipodium + membrane_load → `adherent_active_spread`, gated on **spread
area / protrusion velocity**. (iii) erm, turnover as cortex-coupled modulators.
*Readiness: now (parameters already resolved); needs the traction gate authored.*

**Phase B — Ventral stress fibers.** First EXPERIMENTAL activation. Blocked on:
ratify `N_filaments` (→ μ_SF), split `sf_myosin_*` prefix, FA LIVE. Gate: single-SF
tension band (Kumar 2006 ~10-30 nN), basal-plane diagnostic, denylisted from γ.

**Phase C — Volume / membrane remodelling.** (i) osmotic_regulation — cheapest
(a scalar setpoint updater on the existing turgor force); blocked only on an MCF7
`Lp`. Gate: τ_RVD in seconds-minutes (Hoffmann 2009). (ii) membrane_reservoir —
blocked on `σ_crit_bleb` + `f_excess`. Gate: bleb nucleation threshold + tension
buffering.

**Phase D — Internal mechanics (confined_migration).** (i) LINC — blocked on
`k_linc` (nonlinear nesprin; propose effective stiffness or tabulated). (ii)
intermediate_filaments — blocked on the nonlinear strain-stiffening constitutive
law (tabulated potential). (iii) microtubules — blocked on `Y_stretch` + the
**dt/CFL impact** (stiff bending may force a global dt cut — quantify first).
Gates: nuclear strain under confinement; perinuclear cage stiffness; MT
compressive load-bearing.

**Phase E — Multicell.** (i) cadherin_junction — most runnable (catch-bond oracle
+ derived k_trans); needs the 2-cell build + partner search + the Iturri SE
registration. Gate: junction tension. (ii) junctional_actin — STUB; blocked on the
α-catenin catch set. Gate: junction-to-cortex force transmission.

**Phase F — Integration.** Promote validated subsets in `full_physiological`,
author recipe-level gates, and only then make integrated claims. **Never enable
full_physiological globally; no biological-closure claim until pairwise AND
recipe-level gates pass.**

## 5. Cross-cutting infrastructure (parallel to the phases)

- **Surface manifold (geometry-only).** Wire `surface_manifold` as broad-phase +
  local frames + the single soft normal-confinement (PI Magic-Number Block for
  `k_conf` ∈ [3e-6, 3e-5] N/m), with the mesh-resolution grid-invariance master
  gate. HARD: no manifold DOF in the γ budget. Needed by confinement (Phase D) and
  contact (Phase E).
- **GPU-main port.** P0 stack is GREEN. Clear the P1 CPU debt before any P1
  compartment runs at production scale: native ForceCompute for erm + substrate
  (simple springs), then the binder on-device nlist port. Track in
  `GPU_MAIN_PORT_*`.
- **Cortical-γ depth frontier (Gate-B).** Proceeds independently. The compartment
  platform does not wait on it.
- **KB / citation integrity.** Each ratified parameter → a `SourceEvidence` row +
  `refresh.sh` so the `verdict=OK` gate passes before the value is cited in a
  deliverable (cadherin Iturri + LINC f_rest re-anchor are already queued).

## 6. Definition of done (per compartment / overall)

- **Per compartment:** parameters ratified; wired LIVE with disabled-path
  bit-identity; no-γ-contamination test green; pairwise gate authored + PASS (or an
  honest PARTIAL/REFUTE recorded); GPU/CFL handled; recipe updated.
- **Overall (a working physiological cell):** `suspended_round` (with Gate-B
  resolved) + `adherent_active_spread` pass their recipe-level gates on the full
  baseline; the cell reproduces the MCF7 cortical tension, traction, and spread
  observables within band, from a physiological baseline, with every compartment at
  its setpoint and no γ contamination. Multicell adds the junction-tension gate.

---

*The per-compartment parameter candidates, gate designs, wiring steps, test plans,
and effort estimates are in `COMPARTMENT_ACTIVATION_DETAILED_PLANS_2026-06-09.md`.*
