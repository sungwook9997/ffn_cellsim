# H.4 — FA + motor-clutch in HOOMD

**Budget**: 3 weeks (Plan v2 said 2 wk; +1 wk for D5 Stam-Hocky minifilament)
**Branch**: `phase1/h4-fa-clutch` (cut from `v2/foundation`)
**Reference**: Plan v2 §3 Unit H.4, [PHASE_0_3_DECISIONS.md](../PHASE_0_3_DECISIONS.md) D2/D5/D6/D7, [AFINES_ALGORITHM_NOTES.md §3/§4/§6](../AFINES_ALGORITHM_NOTES.md)
**Owner**: Sub Session. Runs in parallel with Main’s H.2/H.3 once H.1 BAOAB freeze lands.
**Prereq**: H.1 BAOAB integrator (`ffn_sim/integrator/baoab.py`) merged to `v2/foundation`. Read-only for Sub.

## Goal

Build full-fidelity focal adhesion + motor-clutch dynamics: each FA is a
HOOMD bond group of dynamic integrin↔ligand bonds with D2 Pereverzev
catch-slip off-rate. Reproduce v1 Worker B's KU-2.x verdicts (per-clutch
force 5–20 pN, biphasic/saturating, F* lifetime peak, FA growth Hill) and
extend to the ECM substrate produced by H.1.

## Deliverables

| File | Content |
| --- | --- |
| `ffn_sim/bridge/fa.py` | FA particle group + integrin↔ligand bond setup. v2 `FocalAdhesion` is a HOOMD bond group + metadata dict (NOT the v1 archived `bridge/types.py` dataclass). |
| `ffn_sim/bridge/integrin_bonds.py` | D2 Pereverzev catch-slip dynamic bond updater for integrins. `k_off(F) = k_s·exp(F/F_s) + k_c·exp(−F/F_c)` per KU-2.5; KU-2.18 defaults. Python `hoomd.custom.Action`. |
| `ffn_sim/bridge/motor.py` | D5 Stam-Hocky myosin minifilament (shared module with `ffn_sim/cortex/myosin.py` — extract common implementation). |
| `ffn_sim/bridge/fa_growth.py` | FA growth emergent from clutch population dynamics. NO Hill wrapper — use direct count of bonded integrins as the FA size proxy. Hill function (KU-2.17) only invoked in `tests/validation/` for oracle comparison. |
| `ffn_sim/validation/pereverzev.py` | Pereverzev closed-form `k_off(F)` oracle for `tests/`. |
| `ffn_sim/configs/phase1_h4.yaml` | KU-2.x parameters: 20–50 nascent FAs per cell + 5–10 mature, N_total clutches per FA, Pereverzev/Bell-Evans rates, Hill stepping params. |
| `ffn_sim/tests/test_h4_topology.py` | Topology smoke. |
| `ffn_sim/tests/validation/test_ku24_motor_clutch.py` | KU-2.4 motor-clutch biphasic/saturating verdict — match v1 Worker B 0.34 % prominence saturating call within ± 0.1 %. |
| `ffn_sim/tests/validation/test_ku25_catch_peak.py` | KU-2.5 catch-bond lifetime peak at F* ≈ 30 pN. |
| `ffn_sim/tests/validation/test_ku217_fa_growth.py` | KU-2.17 FA growth threshold at `F_th^{per-FA} = 50 pN` with `N_engaged ≈ 10`. |
| `ffn_sim/outputs/h4/REPORT.md` | KU-2.x gate evidence, per-clutch force histogram, F-V curve, FA growth trace. |

## Implementation spec

### FA topology

- Per cell: 20–50 nascent FAs + 5–10 mature FAs.
- Each FA: cluster of integrin particles at a substrate position, type `integrin`, count `N_total = 50` per FA (KU-2.4).
- Substrate: provided by H.1 (`ffn_sim/ecm/`) once H.1 lands. For H.4 isolation testing, use a synthetic ligand particle plane at z = 0.
- Engaged clutches = active `md.bond.Harmonic` bonds between `integrin` and the nearest ligand/ECM particle.
- Bond params: `k_int^bare = 1e-3 N/m` (KU-2.7 baseline), `r0 = 0` (idealised attachment).
- Vinculin allostery (KU-2.7): `k_int^eff = k_int^bare · (1 + α·N_vin)` where N_vin is a per-FA dynamic count. Implement as a Python-side multiplier applied each step before HOOMD bond force evaluation.

### Integrin off-rate (D2 + KU-2.5)

- Pereverzev two-pathway: `k_off(F) = k_s·exp(F/F_s) + k_c·exp(−F/F_c)`.
- KU-2.18 defaults: `k_s = 0.5 s⁻¹, F_s = 30 pN, k_c = 0.4 s⁻¹, F_c = 7 pN`.
- `IntegrinBondUpdater(hoomd.custom.Action)` runs every N steps:
  - For each currently-engaged integrin bond, read `Force` from HOOMD bond data.
  - Compute Pereverzev `k_off(|F|)`.
  - Sample break probability `1 − exp(−k_off · Δt_batch)`.
  - On fire: remove bond, set `integrin.state = unbound`.
- Binding: for each unbound integrin within `R_FA = 1.5 μm` capture radius of a free ligand (Plan H.4), sample `k_on = 1 s⁻¹` per Δt_batch; on fire add bond.

### Motor stepping (D5 + D6)

- Shared `ffn_sim/bridge/motor.py` (also imported by `ffn_sim/cortex/myosin.py`) implements the Stam-Hocky minifilament.
- Per-head per-step: read force on the head's bond to actin, project onto filament direction, compute D6 Hill velocity, advance `pos_a_end`.
- See AFINES algorithm notes §4.3 for the kernel structure.

### Force-velocity (D6)

- Hill 1938, `v(F) = v0·(F_s − F)/(F_s + F/a)`, `a/F_s = 0.5` (Kovács 2003 NMII).
- NOT AFINES piecewise-linear stall.
- `v0 = 1 μm/s, F_s = 0.5 pN, a = 0.25 pN`.

### Talin / vinculin substeps (KU-2.6 / KU-2.7)

- Talin unfolding emerges from Bell-Evans on talin domain bond, not from a wrapper. Per-domain `k_unfold(F) = k_0·exp(F·x_β / kT)`.
- Vinculin recruitment ODE → per-step explicit Euler `dN_vin/dt = k_rec·(N_max − N_vin) − k_unbind·N_vin`. Apply `k_int^eff` allostery to bond stiffness each step.

### Excluded volume (D7)

- `md.pair.LJ` WCA repulsive-only between `integrin` × `ligand` and `integrin` × `actin_ecm` to prevent overlap.

## Validation acceptance

| Gate | Criterion | KU |
| --- | --- | --- |
| Per-clutch force | 5–20 pN at mature FA (v1 PASS) | KU-2.12 |
| Biphasic / saturating verdict | saturating, 0.34 % prominence ± 0.1 % vs v1 | KU-2.8 (v1: saturating) |
| Catch-bond lifetime peak | F* ≈ 30 pN ± 5 pN (Pereverzev maximum from closed-form) | KU-2.5 |
| FA growth Hill | `F_th^{per-FA} = 50 pN` triggers `N_engaged ≈ 10` | KU-2.17 |
| Pereverzev oracle | HOOMD-emergent k_off matches `ffn_sim/validation/pereverzev.py` closed-form within ± 5 % across F ∈ [0, 60] pN | KU-2.5 + D2 |
| Hill oracle | HOOMD-emergent v matches `v(F)` Hill form within ± 5 % | D6 |

## Open implementation questions

- Standalone vs ECM-coupled test: H.4 should ship a standalone version
  (synthetic substrate) by week 2, then integrate H.1's ECM by week 3.
- Multi-head minifilament cost: per minifilament 20 heads × dynamic
  binding to actin each step is expensive. Profile against the
  `XlinkUpdater` pattern; expect ~ms per minifilament per update call.
  With N_FA × N_motors small (∼50), tractable.
