# I2b excluded-volume — integration patch-notes (Session C, `ac/solid-ev`)

Per AC_PARALLEL_SESSIONS_2026-07-16.md §1.3: this track does **not** edit `ff/`. Every `ff/` change the
increment implies is written here as a documented patch-note + adapter; the **lead applies these in spine
order** (I2b lands after I2, before I3). Nothing below is applied in this worktree.

Deliverable of this track (all net-new under `ac/solid/`, CPU-green on the dev Mac):
- `wca_analytic.py` — pure-NumPy WCA oracle (ground truth).
- `steric_reference.py` — pure-NumPy brute-force + cell-list EV force (the algorithm oracle) + the CFL
  `k_EV` stiffness term.
- `steric_warp.py` — the Warp-CUDA `StericForce` (§1.4 primitive) + hash-grid WCA kernel (authored source;
  lead gates it natively).
- `params_i0b2b.yaml` — `sigma_EV` (GAP→PI) + `k_EV` (Magic-Number Block, my call).
- Tests `tests/ac/solid/{test_wca_oracle,test_steric_reference}.py` — 24 gates, all green.

---

## Patch-note 1 — SUBSUME the MT-tip↔cortex `soft_contact_kernel` (replace, not add)

**File / symbol:** `ff/network_warp.py :: soft_contact_kernel` (def L124) and its three launch sites
(L1232, L1324–1325, L1383), fed by `mtp_d` (built L942–946: each MT tip ↔ its nearest cortex node via
`cKDTree`, params `r_mt_c`, `k_mt_c`).

**Before:** the ONLY steric channel is a hand-wired *one pair per MT tip* — the outermost MT bead pushed
off its single nearest cortex node. Cortex↔cortex, SF↔SF, SF↔cortex, MT-shaft↔anything all interpenetrate
(ENGINE_ARCHITECTURE §3 gap; the mesh has no steric volume under compression).

**After:** `ac/solid/steric_warp.StericForce` computes WCA excluded volume over **all** filament nodes via
one shared `wp.HashGrid`. The MT-tip↔cortex case is then a **subset** of the general all-fiber EV (MT tip
and cortex node are different `fiber_id`s → they already repel through the general kernel).

**Double-count guard (REQUIRED):** the all-fiber EV must **SUBSUME** the old kernel, i.e. **delete** the
`soft_contact_kernel` launches (L1232/1324/1383) and the `mtp_d` build (L942–946) **in the same commit**
that wires in `StericForce` — do NOT run both, or the MT-tip↔cortex pair carries *two* repulsions. Adapter:
set `sigma_EV` so the general shell reproduces the intended `r_mt_c` contact range for the MT-tip↔cortex
pair (or, if MT and actin get different steric diameters per `params_i0b2b.yaml::sigma_EV.candidate_values`,
verify the tip↔cortex range is preserved). Regression: with EV parameters matched, the tip↔cortex force on
a confinement test must equal the old `soft_contact_kernel` result to tolerance before the old code is
removed (author this as a native parity check).

**Same-commit checklist:** remove `soft_contact_kernel` def + 3 launches + `mtp_d`; add `StericForce`
construction (fiber_id from the merged filament indexing) + one `accumulate` call per inner step.

---

## Patch-note 2 — add `k_EV` to the CFL `kmax`

**File / symbol:** `ff/network_warp.py` `kmax` assembly (L855–870; sibling build variants L1513, L1593,
L1706). Currently `dt_mu = 0.1 / kmax` with `kmax = max(bending κ/seg³, xl_k, turgor, nucleus, IF, membrane…)`.

**After:** fold in the steric stiffness — `kmax = max(kmax, k_EV_cfl)` where
`k_EV_cfl = ac.solid.steric_reference.cfl_stiffness_term(sigma_EV, epsilon, f_cap)`:
- without a cap this is the WCA well curvature `U''(r_min) = k_EV`;
- with the optional `force_cap_EV` it is the deeper-overlap stiffness `U''(r_eq(f_cap))` — the stiffness the
  integrator actually sees, never the divergent `r→0` core.

**Guard:** `k_EV` is the numerical repulsion scale from `params_i0b2b.yaml` (Magic-Number Block: derived from
the operating load × penetration tolerance, grid-invariant, NOT tuned to a crowding outcome). **CFL cost is
real** — with the physical `sigma=7 nm`, `k_EV` is large and throttles `dt`; this is a logged argument for
the `sigma_EV` discretization decision (PI), never a license to soften `k_EV` below the tolerance threshold
or to loosen a gate.

---

## Patch-note 3 — `StericForce.accumulate` in the inner force-assembly loop (§1.4)

**Contract:** AC_PARALLEL_SESSIONS §1.4 — each force primitive is `accumulate(state, out_force)` adding its
per-node contribution; the lead-owned integrator sums `MyosinForce`(I3) + `StericForce`(I2b) +
`PressureCoupling`(I1a) + the ported bending/WLC/xlink terms.

**Wire-in:** construct one `StericForce(fiber_id, sigma_EV, k_EV, device, active=active_mask, force_cap=…)`
at cell build (over the merged filament node indexing, with `active` = the dormant/live mask so parked nodes
are skipped), and call `steric.accumulate(state, f_d)` once per inner mechanical iteration, alongside the
other primitives, BEFORE the integrator advances. Own-row read-modify-write is race-free (each node's row is
written by exactly one thread) → composes with the other accumulators with no atomics.

**OFF == bit-identical:** if `StericForce` is simply not in the primitive list, forces are unchanged; and
when every pair is beyond `r_c` the kernel adds exactly zero (test_steric_reference::test_no_overlap_is_zero_force).
So EV-OFF and no-overlap are both bit-identical regressions — the required OFF gate.

---

## Patch-note 4 — reconcile with I7 LINC/plate soft-contact (single-channel invariant, hard-truth #6)

Hard-truth #6 lists **`LINC+soft_contact`** as one of the 6 double-count instances, guarded **at I7** ("land
the LINC wire-in + the soft-contact/plate term TOGETHER, per FSI_ALL_COMPARTMENT"). I2b's all-fiber EV must
be **the single steric channel**: when I7 adds its nucleus cap→LINC load path with a plate/soft-contact
term, that term must NOT re-implement filament↔filament (or filament↔nucleus-surface) steric that the I2b
EV already covers.

**Invariant to hold at I7 (flag for the lead):** exactly ONE steric channel acts on any given node pair.
Concretely — (a) filament↔filament steric = I2b `StericForce` only; (b) the I7 LINC tether is a *distinct*
force family (an attractive/elastic tether along the LINC bond), NOT a steric term, so LINC + I2b-EV on the
same node is two DIFFERENT families (not a double-count); (c) any nucleus-envelope↔filament *contact*
(soft-contact/plate) must be routed through the **same** EV mechanism (either extend `StericForce` to include
envelope surface nodes with a nucleus `fiber_id`, or explicitly hand that pair class to I7's term and exclude
it from I2b) — never both. Decision deferred to the I7 combined baseline (this note is the hand-off).

---

## Native-gate spec (lead runs on the gbook A5000, in spine order)

Authored CPU-green here; the lead runs these natively (the dev Mac cannot run Warp-CUDA):

1. **Warp ↔ brute-force parity (grid-invariance).** `StericForce._accumulate_pos` on a crowded native patch
   == `steric_reference.steric_force_bruteforce` to summation tolerance, for hash-grid cell sizes {r_c,
   1.5 r_c, 3 r_c}. Certifies the device HashGrid is a pure accelerator (Magic-Number-Block grid-invariance).
2. **OFF / zero-overlap == bit-identical.** EV not in the primitive list, and an all-separated config, both
   reproduce the pre-I2b forces bit-for-bit.
3. **Compressed-patch resists collapse (finite steric volume).** A locally compressed native filament patch
   develops a net outward (expansive) virial and does not collapse to zero volume — **PASSIVE cell** (no
   myosin), structural only.
4. **No fiber interpenetration on the native cell (crowding gate).** On the full-population native
   `--from-resting` PASSIVE cell (70,686 cortical F-actin + SF/MT/… at physiological density), a
   penetration-depth diagnostic (port the `dcm_neighbor_warp.penetration_depth_kernel` idea to filament
   nodes) reports `max_pen / mean_seg` below tolerance. **Field-render:** a full-res interactive 3-D
   cell-morphology **HTML** (peel/slab/cut) coloured by local steric penetration depth (per the §1.9 viz
   rule + PI 2026-07-02/07); browser-verify with `browser_check.py` — a WebGL grep is meaningless.
5. **CFL / dt ledger.** Report `k_EV`, its `kmax` contribution, the resulting `dt_mu`, and the wall-time
   delta from EV-ON. If `dt` collapses (physical `sigma=7 nm`), that is the `sigma_EV`-discretization
   FINDING for PI (coarse node-node vs segment-segment) — report, do not soften `k_EV`.
6. **k_EV magnitude close-out.** Measure the passive-cell operating per-node load `F_op`, set
   `k_EV = F_op / (tol·r_c)` with `tol ≤ 0.05`, confirm gates 3–4 hold, and REPORT the value + population/
   peak-GPU-byte ledger (never lower a biological density to fit memory or buy dt).

**Scope quarantine (verify finding, Session-C brief):** gates 3–4 pre-I3 certify a **PASSIVE** cell
(structural/regression). The load-bearing **no-interpenetration-under-contractile-load** verdict **DEFERS**
to a re-run after I3 / at I9 (myosin ON). Do not report the passive verdict as the load-bearing one.

---

## Open items surfaced to PI

- **`sigma_EV` discretization decision (I0-B2b GAP).** Physical filament diameter (actin 7 nm, MT 25 nm) is
  under-resolved by node-node EV at 0.5 µm node spacing AND drives a punishing CFL. Faithful options:
  (a) per-type physical `sigma` + **segment-segment (edge-edge) EV** (port `dcm_neighbor_warp`
  `edge_edge_contact_kernel`) — resolution-faithful; (b) node-node EV with coarse-grained `sigma ~ node
  spacing` — a declared numerical simplification. PI call; do not pick silently.
- **Crosslink-adjacent exclusion (I4 reconcile).** Same-filament pairs are excluded; crosslinked
  cross-filament nodes are not. If a crosslink rest length < `r_c`, EV + crosslink coexist on one node pair
  (two families, not a double-count) but produce a standing repulsion — reconcile at I4.
