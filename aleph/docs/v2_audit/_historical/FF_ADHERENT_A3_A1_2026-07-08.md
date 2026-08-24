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

# FF Phase-A wiring — A3 (checkpoint → adherent) + A1 (FF-native tests) (2026-07-08)

**Context.** The FF cortical-mechanics axis is validated on a SUSPENDED cell (`FF_CORTICAL_MECHANICS_STATE_2026-07-08.md`).
The expansion plan's recommended first step (`FF_EXPANSION_PLAN_2026-07-08.md`, PI-confirmed 2026-07-08) is **Phase A3 + A1**:
wire the validated resting checkpoint cell onto the adherent path, and close the FF-native validation debt for the
compliant-substrate / FA-maturation / Piezo modules. Both landed on `dcm/main`.

## A1 — FF-native validation tests (validation debt closed)

The compliant-substrate machinery (`ff/substrate.py`), the FA mechanosensor (`ff/fa_maturation.py`), and the Piezo
reporter (`ff/piezo.py`) were GPU-native + wired but had **zero FF tests** (the pre-existing `tests/test_substrate*.py`
are legacy HOOMD spheroid tests under `archive/hoomd_legacy`, not these Warp modules). Added three analytic-ground-truth
test files under `tests/ff/` (**25 tests, all PASS**; every expected value derived from first principles, not from the
module's own output):

- **`tests/ff/test_substrate.py`** (11) — flat-punch bridge `k_sub = 2·E·a/(1−ν²)` exact + `k_sub ∝ E`; the
  equilibrium-anchor kernel places the anchor at the clutch↔substrate series balance so the actin feels the
  **Bangasser-Odde series stiffness `k_series = k_int·k_sub/(k_int+k_sub)`** (derived independently, checked against the
  kernel); rigid limit `k_sub→∞ ⇒ anchor→rest` (fixed-pin recovered) + soft limit + unbound + Newton-pair reaction +
  CFL + input-validation guards.
- **`tests/ff/test_fa_maturation.py`** (8) — talin Bell parity `k_unfold(F)=k_u0·exp(F·dx/kBT)` (function + one-step Warp
  kernel); **Hill FA-growth fixed point at F=F_th** (the KB constraint `k_g0 = 2·k_d` puts the Hill half-max at
  threshold ⇒ `dA/dt=0` exactly at F_th; grows above, disassembles below) + kernel parity; vinculin steady state
  `N* = k_rec·p·N_max/(k_rec·p + k_diss)` + `k_int_eff = k_int·(1+α·N)`; clutch-load readout; force-gated disassembly.
- **`tests/ff/test_piezo.py`** (6) — the **UNIT-TRAP guard**: at the resting bilayer tension γ_mem≈10 pN/µm the KB-3.10
  logistic reads `P_open≈0.0347` (essentially CLOSED); the dropped-`×1000`-conversion bug (10 treated as pN/nm vs
  γ_half=5 pN/nm) would read `≈0.966` (wide open) — the failure the conversion prevents, checked against the ACTUAL
  `resolve_membrane()` resting tension; + logistic midpoint = 0.5; monotone/saturation; feedback gain = 0 (default-off).

⚠️ **DEBUG RECORD (PI 2026-07-08).** The FA maturation runs downstream of the `ff/fa_clutch` catch-slip bond, whose peak
is **F*≈7 pN** from the recorded Pereverzev params while KB/Kong-2009 report **~30 pN** — surfaced to PI, **kept as-recorded
(NOT retuned to a target)**. Recorded in `test_fa_maturation.py` docstring because the FA-growth threshold F_th=5 pN and the
clutch peak F* sit on the SAME force axis: when the absolute-traction magnitude is later reconciled (a PI call), both move
together. This is the debugging pointer PI asked to leave for the later magnitude pass.

## A3 — adherent run starts from the VALIDATED resting checkpoint (`--from-resting`)

The adherent driver `scripts/ff_crawl_on_substrate.py::build()` built a FRESH un-relaxed cortex (sparse myosin `//160`,
default membrane, no pre-relaxation) — an unphysical baseline that then adhered, violating the physiological-baseline HARD
rule. Added a **`--from-resting`** mode that reproduces the checkpoint recipe (`ff_resting_full_compartment.py`) and
**pre-relaxes to the resting turgor set-point before adhering**:

- dense myosin `n_myo = n_filaments//10`, membrane reservoir `f_excess=0.25`, nucleus `R_nuc=0.70R`, and a pre-relaxation
  to the interphase set-point (`pressure_setpoint=TURGOR_DP0=40 Pa`, biphasic drained solid `K_drained=300 Pa`) via the
  SAME `simulate_whole_cell_compression_on_device` the checkpoint uses.
- **Node-layout safe by construction:** both the checkpoint relaxation and `build()` produce the IDENTICAL
  `[cortex(Nc) ; MT_arms ; MTOC ; nucleus]` layout (both call `merge_aster_into_cortex` then `concat`), so the relaxed
  positions map back with no reconciliation; `merged.pos` is updated so run()'s V0 + face triangulation stay consistent.
- Substrate contact (`z_sub / basal cap / anchors / front`) is then derived from the RELAXED cortex geometry.
- **Backward-compatible:** `from_resting=False` keeps `n_myo_ratio=160`, `f_excess=0.0`, no relaxation, and derives contact
  from `pos_all[:Nc] == cx.net.pos` — identical to the original path.

**Verification (CPU dev-scale).** `--from-resting --static --mature --spread --fa-maturation --piezo`:
- pre-relaxed to the checkpoint state: **ΔP=40.3 Pa → γ=0.171 mN/m** (interphase FF_STAGE6V band 0.15≈0.17), V/V0=0.999,
  R_eq=8.48 µm — i.e. the adherent run starts FROM the validated resting physiological state, not a null baseline.
- adhered: **STABLE ADHERED**, V/V0=1.000, 84/84 FA bound, traction ~18 nN.
- Piezo live: membrane 10.0 pN/µm → P_open=0.0347 (rest≈closed) — matches the A1 gate.
- viewer `outputs/ff/figs/a3_adherent_from_resting_morph.html` **browser-verified** (WebGL OK, no JS errors): cortex woven
  actin + nucleus (3000 beads, 0.70R) + rigid substrate plane + 84/84 FA integrin clutches rendered full-res.

**Production path.** native = `--cortex-fil 38000 --from-resting --microtubules --implicit --device cuda:0` on the gbook
A5000 (pre-relax is `simulate_...`, the same native path the checkpoint uses). CPU is the small-N dev demo.

## Still PI-gated (before any ABSOLUTE-traction claim — not needed for A3/A1)

- **FA patch radius `a`** — sets the absolute `k_sub` magnitude (surfaced, KB-blocked, KB-PIV-5 datum family). Default =
  nascent-adhesion radius 50 nm; only the E-DEPENDENCE (biphasic/durotaxis shape) is gated, no magnitude claim.
- **α5β1-FN catch-slip peak F*** — code 7 pN (Pereverzev-recorded) vs KB/Kong-2009 30 pN. **PI 2026-07-08: keep 7 as-recorded,
  record it as an important later-debugging factor** (done: `test_fa_maturation.py` + `test_fa_clutch.py`).

## Note — pre-existing unrelated test failure (surfaced to PI, NOT mine)

`tests/ff/test_compartments_shared.py::test_resolve_nucleus_bridge_and_bands` FAILS on `dcm/main` independent of this work
(reproduces in isolation; this session touched only the 3 new test files + a 1-line import + the `build()`/`main()` wiring in
`ff_crawl_on_substrate.py`, none of which import `common/compartments.py`). Cause: the test asserts nucleus grid-invariance by
comparing an explicit `E_nuc=5000 Pa` cell against a DEFAULT-`E_nuc` cell, but the default `E_nuc` was changed to ~399 Pa in
the DCM fidelity-remediation session — so the two use different moduli and the intensive product differs. A stale DCM test
(test-contract), not an FF or A3/A1 defect. Flagged for the DCM/compartments owner.
