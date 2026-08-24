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

# Parameter-inconsistency reconciliation — NMII capture radius + L_p

**Date:** 2026-07-24 · **Branch:** `codex/ff-ac-codex` · **Env:** ffn_sim (Mac, CPU tests only) · **Scope:**
the two parameter INCONSISTENCIES flagged by `PARAM_PROVENANCE_AUDIT_2026-07-24.md` rows 12 (NMII capture) and
36 (L_p). Rigorous units-explicit reconciliation; a fix is applied only where the correct value is unambiguous
**and** the change keeps tests green.

---

## 1 — NMII head→actin capture radius

### Values found (file:line)

| Value | Location | Role | Provenance as written |
|---|---|---|---|
| **0.05 µm** | `ac/cell/assemble.py:155` `NMII_CAPTURE_UM` (was) | motor KMC `params.capture_radius` (`assemble.py:559`); default resting-seed bind reach (`assemble.py:793`); NG-1 native gate (`ac/motor/native_gates/ng1_two_filament_stall.py:183`) | "provisional — unused while unbound" |
| **0.210 µm** | `ff/hand_kmc.py:125` `NMIIA_MYOSIN.capture_radius_um` | FF-lane NMIIA hand attach gate | "head_actin_capture_perp 210 nm" |
| **0.6 µm** | ~20 seed-probe scripts (`ac_gate_a_*`, `ac_myosin_*_probe.py`, `ac_645_force_breakdown.py`, `ac_gate_a_straddle_closure.py`, …) via `resting_bound_myosin_capture_um=0.6` | forces resting-bound-myosin seeding in diagnostics | labeled `*_TEST` |

`ac/engine/cortex_motor_slice.py` (GATE-B) does **not** hardcode a capture: `CortexMotorParams.capture_radius`
is a `_require_param` field supplied by the caller. The "GATE-B used 0.6" is the caller (seed probe / slice
driver) passing 0.6, not a constant in the slice.

### Are they the same physical quantity?

**0.05 and 0.210 are the SAME quantity** — the perpendicular point-to-segment reach within which a myosin head
may bind an actin segment (gated in `ac/motor/segment_query.py`: "bind only if point-to-segment distance ≤
this"). Provenance of the correct value (`configs/phase1_h3.yaml:377–386`, KU-3.5 binding fix 2026-05-29,
`docs/KU35_myosin_binding_diagnosis.md`):

> heads sit ~head_rest_length = 200 nm laterally from the actin segment they bind. Eligibility = perpendicular
> distance to the segment LINE ≤ `head_actin_capture_perp` (the head's physical reach via its spring).
> `head_actin_capture_perp: 2.1e-7 m = head_rest_length + 10 nm slack`.
> **"Legacy 50 nm was for the wrong radial-offset placement that left heads ~824 nm from any actin."**

So **0.210 µm = r0_head (200 nm) + ~10 nm slack** is the sourced/derived-correct reach. `assemble.py`'s own
`NMII_HEAD_OFFSET_UM = r0_head = 0.200`, so 0.05 also **internally contradicted the same module's arm offset**.
The 0.05 is the explicitly-RETIRED legacy value, carried as a stale straggler on the `ac/cell` surface while
`ff/hand_kmc.py` already held the fixed 0.210.

**0.6 µm is the same quantity but a deliberate diagnostic loosening**, not a physiological reach. The seed
probes set it to force resting-bound-myosin heads to capture actin in the GATE-A/B convergence and force-budget
studies (all `*_TEST`-labeled). It is a binding-forcing knob, **left as-is** — it is not a physiological
constant and does not belong at 0.210.

### Verdict

**Single physiological value: 0.210 µm** (= r0_head + ~10 nm slack; KU-3.5 fix; matches `ff/hand_kmc.py`). The
0.05 was documented-legacy-wrong; the 0.6 is a distinct diagnostic override, not this quantity's setpoint.

### What I did — FIXED

- `ac/cell/assemble.py:155` `NMII_CAPTURE_UM` **0.05 → 0.210**, with a provenance comment cross-referencing
  `ff/hand_kmc.py` / `phase1_h3.yaml` / KU-3.5, noting it is NOT the 0.6 diagnostic reach.
- `ac/motor/native_gates/ng1_two_filament_stall.py:93` — updated the now-stale `< 2·NMII_CAPTURE_UM (0.05)`
  parenthetical to `≪ 2·NMII_CAPTURE_UM (0.210)`. Comment only; the `ACTIN_SPACING_UM = 0.004 ≪ 2·capture`
  invariant still holds (more so).

**Why this is safe, not a guess:**
- **NG-1 native gate (GPU, not runnable on Mac): binding is invariant.** `_build_two_filaments` places the two
  actin lines at `y = ±NMII_HEAD_OFFSET_UM = ±0.200`, i.e. exactly at the head arm reach, so a head's
  point-to-segment perpendicular distance is ≈ 0. Widening the capture 0.05 → 0.210 cannot un-bind or newly-bind
  anything in that 2-filament geometry (0 ≤ 0.05 ≤ 0.210), and the measured isometric-stall force comes from the
  already-bound heads — unchanged. This is deterministic from the geometry, not an empirical guess.
- **GATE-A resting baseline: unaffected** — heads are UNBOUND at t0; the constant does not enter.
- **GATE-B / seed probes: unaffected** — those pass their own `capture_radius` / `resting_bound_myosin_capture_um`
  (0.6), overriding the module default.
- The change only widens the default free-head bind reach used by the `myosin_activate.py` diagnostic KMC to its
  physiologically-correct value (more actin reachable = correct), which that tool reports honestly.

**CPU tests:** `pytest aleph/tests/ac/` → exit 0, **~986 passed / 84 skipped** (skips are CUDA-only kernel
gates), **0 failed / 0 error**. GREEN.

**Caveat (GPU re-validation):** NG-1 and the native GATE-B path run on the gbook A5000, not the dev Mac. The
NG-1 stall-force outcome is invariant by the geometry argument above, but a native NG-1 + GATE-B re-run should
still confirm on the next gbook session as routine hygiene (no result change expected).

---

## 2 — L_p (membrane hydraulic conductivity)

### Values found (file:line)

| Value | Location | Unit as written | Source |
|---|---|---|---|
| **1.6e-8** | `ac/cell/assemble.py:133` `L_P`; `ac/cell/ng4_payload.py:53` | µm/(s·Pa) | "gamma_floor code value; draft/MCF7" |
| **1.0e-12** | `ac/fluid/params_i0b1.yaml:141` `L_p.value` | m/(s·Pa) | Jung 2011 MCF7/AQP5; KB-DRAFT-3.B-26; Kedem-Katchalsky form |

Both feed the same law `s_water = L_p·(σ·ΔΠ_osm − ΔP)` — code `L_P` → `MembraneFluxBC` (`assemble.py:886`,
`ng4_payload.py:85`); the yaml row documents the same BC (`ac/fluid/boundary.py`).

### Units-explicit reconciliation

Convert the yaml value into the code's µm-based units:

```
1.0e-12 m/(s·Pa) × (1e6 µm / 1 m) = 1.0e-6 µm/(s·Pa)
```

Compare on a common basis:

```
yaml  = 1.0e-6  µm/(s·Pa)
code  = 1.6e-8  µm/(s·Pa)
ratio = 1.0e-6 / 1.6e-8 = 62.5×
```

### Verdict — GENUINELY OFF (~62.5×), not agree-after-units

The unit conversion is a clean ×1e6 and, **after** applying it, the two surfaces still differ by ~62.5× (the
code value is ~63× *smaller*). The audit's arithmetic (1e-12 m = 1e-6 µm ≠ 1.6e-8 µm) is confirmed. The yaml
itself already flags this: `reason_provisional: "…unit reconciliation open"`, `note: "Reconcile against code
value 1.6e-8 (gamma_floor)"`.

**Sourced/correct value = the yaml's 1.0e-12 m/(s·Pa) = 1.0e-6 µm/(s·Pa)** (Jung 2011 MCF7/AQP5,
Kedem-Katchalsky form, KB-DRAFT-3.B-26 — a real MCF7-specific anchor, and physically an ordinary cell L_p
magnitude). The code's 1.6e-8 µm/(s·Pa) (= 1.6e-14 m/(s·Pa)) is the older unsourced "gamma_floor" number and is
~63× too small.

### What I did — RECOMMEND (not fixed)

No code change. Reasons this is a recommendation, not an edit:
- It is a **62.5× change to a physics magnitude on the fluid runtime**, and the yaml row is explicitly
  **`owner: PI`, `provisional: true`, "unit reconciliation open"** — i.e. a PI-owned open reconciliation, which
  the gate discipline (no unilateral change of a PI-owned provisional constant) routes to PI, not to a silent
  in-place edit.
- **Not load-bearing at rest**, so leaving it does not corrupt GATE-A/B: at the resting baseline ΔΠ_osm = ΔP_hyd
  (zero net flux, `params_i0b1.yaml:174`), so `L_p` multiplies ~0. It bites only under perturbation (bleb /
  spreading / osmotic-response dynamics), which is post-GATE-A/B — so a fix carries dynamic-regime risk with no
  resting-baseline benefit.

**Recommended action for PI:** set the code surfaces to the sourced value —
`assemble.py:133` and `ng4_payload.py:53` `L_P = 1.6e-8 → 1.0e-6` µm/(s·Pa) — to match
`params_i0b1.yaml` (1.0e-12 m/(s·Pa), Jung 2011), and drop the yaml's `unit reconciliation open` flag. Best done
together with the membrane-flux / osmotic-response dynamic gate (where L_p first becomes load-bearing) so the
63× change is validated in the regime it affects, not in the rest state where it is inert.

---

## Two-line verdicts

1. **NMII capture radius — SINGLE physiological value (0.210 µm).** 0.05 and 0.210 are the same quantity;
   0.05 was documented-legacy-wrong (retired by KU-3.5), 0.210 = r0_head + 10 nm slack is sourced. **Fixed**
   `assemble.py` 0.05 → 0.210 (+ NG-1 comment); CPU tests green; NG-1/GATE-B binding invariant by geometry.
   The 0.6 in seed probes is a distinct `*_TEST` diagnostic override, left as-is.
2. **L_p — GENUINELY OFF (~62.5× after units).** 1e-12 m/(s·Pa) = 1e-6 µm/(s·Pa) ≠ code 1.6e-8 µm/(s·Pa).
   Sourced value is the yaml 1e-12 m/(s·Pa) (Jung 2011). **Recommend** aligning code to it (PI-owned, ×62.5,
   inert at rest) — not fixed here.
