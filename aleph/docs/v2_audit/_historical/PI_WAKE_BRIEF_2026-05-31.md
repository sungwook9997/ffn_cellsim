---
archived_on: 2026-07-28
superseded_by: aleph/docs/v2_audit/AC_EXECUTION_PLAN_2026-07-25.md
reason: >
  Written BEFORE the 2026-07-25 PI reframe, i.e. for a different objective — forward prediction
  and magnitude matching, rather than inferring per-cell-type parameters with the gate on
  mechanical connectedness. Archived, not deleted: its measurements and reasoning stand as a
  record of what was true then. Nothing in it may be quoted as current state; STATE.md is that.
  Selected mechanically: pre-reframe AND cited by no live file (code, STATE.md, CLAUDE.md,
  cell_engine/, gate_contracts/, tests, Makefile). Citations from run outputs and from other
  pre-reframe documents were not treated as protective.
---

# PI wake brief — 2026-05-31 autonomous session (12h /loop)

> Single boot-point. Branch `phase1/h3-cortex` (NOT pushed to ffn/foundation — your call).
> Full trail: `outputs/h3/production/AUTONOMOUS_LOG_2026-05-31.md`. Root-cause brief:
> `docs/v2_audit/KU35_FLOOR_ROOT_CAUSE_2026-05-31.md`. Notion progress page: 370120da…2732.

## Outcome in one line

**The months-long KU-3.5 cortical-tension floor is FULLY DIAGNOSED, and every blocker on the
path to measuring it is fixed + committed (all regression-green). The one remaining piece is a
core-physics myosin redesign → your sign-off (decision #1).**

## What landed (this session, all regression-green at each step)

| Commit | What |
|---|---|
| `75f464e` `c960406` `d33230b` `3d70e6b` | **G1 floor blocker (a):** FA contact-footprint seeding — the cortex shell was origin-centred, disjoint from the z=0 substrate, so ~0 clutch formed (9/52). Now seeds each FA under a south-cap cortex bead → clutch 9→**45** loaded. + per-bead fan-in/exclusion safety. Production-capable. |
| `f2f4e75` | **Pre-existing long-run crash FIXED.** HOOMD nlist has a hard 7-bonded-exclusions/particle cap; the uncapped xlink binder + the myosin cap being blind to other bonds let a cortex bead reach degree 8 → crash (the G1-OFF baseline crashed too — NOT a G1 bug). Both binders now share a per-bead degree budget. Verified: 2.5e6 steps crash-free. |
| `db1c756` | **EXTEND modules** (Track A `ecm/substrate.py`, H.8 `cell/membrane_surface.py`, H.9 `cell/nucleus.py`) scaffolded + adversarially verified, **default-off**. Nucleus net-force/COM-drift fixed (radius/direction decorrelation). |
| `206044b` | **H.8 membrane_surface integrated** into Cell.build (None-guard default-off) + the **KU-3.5 root-cause brief**. |

Also resolved: **G2** — myosin DOES step (`step_advances` 203→1396 over 2.5e6 steps); the
earlier `step_advances=0` was a too-short-sim artifact, not a defect.

## DECISION #1 (load-bearing) — myosin contractile mechanism redesign

**The KU-3.5 floor cause:** the myosin **binned-r0 ratchet is a LUMPED PROXY** (violating the
CLAUDE.md "no lumped mechanisms" hard rule). A "step" only relabels the head–actin bond to a
lower rest-length bin **on the same actin bead** — it does NOT transport actin material (the
AFINES `pos_a_end` grip-walk), so the spring relaxes out each tick, no strain accumulates →
γ_soft ~25,000× under band, r/r0=1.000 even after 1396 advances. (Force ceiling 0.33 pN <
F_stall 0.5 pN; bipolar sidedness absent. The method-of-planes measurement is faithful.)

**Proposed fix (NOT implemented — core physics, your call):** convert the proxy to a real motor —
advance the head's **attachment bead index** along its bound actin filament toward the minus end
(the fractional accumulator already exists), re-target the head–actin bond to the new downstream
bead while keeping the head in place → the bond is stretched → **sustained** pull; keep r0≈0
(idealized cross-bridge). + restore a force ceiling that reaches F_stall + enforce bipolar
sidedness. Detail + file:line in `KU35_FLOOR_ROOT_CAUSE_2026-05-31.md` §3. **Once this lands,
re-run the now-crash-free v4 production to confirm γ lifts into [0.35, 0.65].** (Note: a full v4
production run is ~40 h wall at the constrained dt — a dt co-tune / GPU is worth considering.)

**A concrete, implementable design brief is now ready: `KU35_GRIP_WALK_DESIGN_2026-05-31.md`** —
exact file:line change table (the single load-bearing edit is the kernel at `myosin.py:985-1022`; the
bond re-target primitive ALREADY exists via the per-tick `set_snapshot`), an **opt-in `stepping_mode`
(default = current `binned_r0`, for safe A/B)**, a 3-tier validation plan whose Tier-1 micro-diagnostic
proves sustained-tension-vs-relaxation in MINUTES (~1e5 steps, no 40 h run), and the PI design
decisions to ratify — esp. **(a) cortex filament polarity** (the cortex actin has NONE; minus-end must
be ASSIGNED — recommend label-only "minus = bead 0") and **(b) the bipolar sidedness binding rule**.
**Critical new finding in the brief:** at v0=0.2 µm/s + the FA-limited dt, walking ONE integer bead
takes ~3.7e6 steps (2.5 s sim) → the walk must use a CONTINUOUS sub-bead `pos_a_end` to make force
from step 1, AND **prior production runs (≤4e5 steps) were too short to show contraction even if the
mechanism had been right** — a run-length issue independent of the proxy. **Process: ratify the design
(esp. a+b), then it implements from the spec** (additive opt-in first, A/B, then make default).

## DECISION #2 — ratify or revert the autonomous PI-flags (all documented, all reversible)

- **Test-contract change:** `test_fa_physical_capture_radius_barely_clutches` → `_fully_clutches`
  (assertion inverted; it previously characterized the floor BUG, which G1 fixes). Revert = `git
  revert` the G1 commit chain. + 3 ligand-z=0 → contact-footprint assertion updates.
- **HOOMD library-capacity constants** (not physics; same class as the existing `MAX_HEADS_PER_BEAD`):
  `_MAX_CLUTCH_FANIN=1`, the south-cap depth, `_MAX_CORTEX_BEAD_DEGREE=6`. Tune if you prefer.
- **Doc-only:** `ecm/substrate.py effective_E_sub` diagnostic is N/m³ not Pa (relabel before any
  E_sub gate). **Nucleus** seeding leaves a ~4% residual net force (antipodal-pair seeding zeroes
  it exactly — flagged refinement).

## DECISION #3 — EXTEND scope (carried from the research handoff, your earlier boot prompt directed proceeding)

- substrate (Track A) + nucleus (H.9) Cell.build integration is the remaining default-off wiring
  (membrane done). Track C composite KU-3.5/KU-3.1 (γ_total = γ_membrane + γ_cortex) is a
  gate-contract change → you. In-vivo ECM presets gate on H.9 pore.
- Carried from `PI_DECISION_QUEUE_2026-05-30.md`: A2 (myosin stepping — now answered: artifact +
  the deeper proxy issue above), C1 integrin F*→30 pN, D1 ffn/foundation push (77+ commits), etc.

## What is NOT done (next session)
- The myosin redesign (decision #1) — gated on you.
- substrate + nucleus Cell.build integration (default-off plumbing).
- Figures for the findings (being generated this session into `outputs/h3/figs/`).
