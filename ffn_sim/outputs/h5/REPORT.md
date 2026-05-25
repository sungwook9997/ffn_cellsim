# H.5 — Lamellipodium 단계 1 (Bieling/Funk greenfield) closeout report

**Branch**: `phase1/h5-lamellipodium` (cut from `phase1/h3-cortex` @ `885db0c`)
**Status**: 🟨 단계 1 implementation freeze landed
**Brief**: [`ffn_sim/docs/briefs/H5_lamellipodium.md`](../../docs/briefs/H5_lamellipodium.md) — derived from PHASE_0_3_DECISIONS.md §D1 (brief was deferred per PHASE_0_CLOSEOUT.md row C until H.3 design choices landed).
**Spec ref**: PHASE_0_3_DECISIONS.md §D1 (greenfield Arp2/3, Bieling 2016 + Funk 2022 mechanisms).

## Scope (단계 1)

PI authorization "프롬프트 작성 후 H.5 그리고 마무리까지 진행바람". Single-iteration drop: brief + Sanity Gate + topology + three batched Updaters + tests + viz + KU-5.x skeletons + REPORT.

## Deliverables landed (단계 1)

| File | Status |
| --- | --- |
| `ffn_sim/docs/briefs/H5_lamellipodium.md` — derived from D1 spec | ✅ NEW |
| `ffn_sim/cell/lamellipodium.py` (~810 lines) — Sanity Gate §1–6 + ResolvedH5 + WaveMembranePin + 3 D1 Updaters + build_lamellipodium_simulation | ✅ NEW |
| `ffn_sim/cell/__init__.py` — exports | ✅ extended |
| `ffn_sim/configs/phase1_h5.yaml` — KU-5.x literal params + acceptance bands | ✅ NEW |
| `ffn_sim/tests/test_lamellipodium.py` (17 PASS / 1 SKIP) | ✅ NEW |
| `ffn_sim/tests/validation/test_ku5x_lamellipodium.py` — KU-5.1/5.2/5.3 production skeletons + oracle utilities (4 PASS / 3 SKIP opt-in `H5_KU5_PRODUCTION=1`) | ✅ NEW |
| `ffn_sim/scripts/h5_vis.py` — 5 figures (WAVE topology, Bell-Evans rate curves, Updater activity, Bieling oracle, pipeline summary) | ✅ NEW |
| `ffn_sim/outputs/h5/REPORT.md` + `figs/*.png` (this document + 5 figures) | ✅ |

## Architecture: cell/lamellipodium.py (D1 greenfield)

### Topology
- **WAVE membrane plane** at `y = +Y_max` (top of box, opposite cortex shell which sits at origin).
- `n_WAVE` WAVE particles (`wave_particle` type) at uniform random (x, z) on the plane.
- **Pinned to plane** via `WaveMembranePin(md.force.Custom)` — harmonic spring in y: `U = ½ k_wave_pin (y − Y_max)²`. CFL gate inherited from `cortex/erm.py` pattern: τ_pin = γ_b / k_wave_pin = 3.91 μs ≫ dt = 13 ns (CFL safe).
- **One mother actin seed per WAVE** at construction (D1 brief option (a) per H.5 brief open questions). Mother type `actin_lamel`; tangent = −ŷ (grows away from membrane into cytosol).

### D1 Bell-Evans rate Updaters (D2 batched, mirrors `crosslinkers.XlinkBondUpdater` pattern)

1. **`BarbedEndElongationUpdater`** — `k_elong(F) = k_elong⁰·exp(−F·δ_elong/kT)` with `k_elong⁰ = 11.6 /s`, `δ_elong = 2.7 nm` (Bieling 2016). On fire: appends actin bead at rest length ℓ₀ along the barbed-end tangent, demotes old tip to interior, promotes new bead to barbed-end status.

2. **`ArpBranchingUpdater`** — `k_b(F) = k_b⁰·(1 − 0.2·F/F_stall)` with `k_b⁰ = 0.037 /s` per WAVE (Bieling 2016). **Funk 2022 abortive emerges naturally**: when `F > F_stall/0.2` per WAVE, `(1 − 0.2·F/F_stall)` goes negative and is **clamped to zero** — no separate scalar override per D1 spec. On fire: inserts daughter actin at mother's barbed-end + ℓ₀ along a tangent rotated 72° from mother (random azimuthal in plane perpendicular to mother). Branch bond + branch angle (`md.angle.Harmonic` at t₀ = 72°, k = 100 pN·μm/rad² per D1 default).

3. **`CappingUpdater`** — `k_cap(F) = k_cap⁰·exp(−F·δ_cap·sin θ/kT)` with `k_cap⁰ = 3 /s`, `δ_cap = 0.3 nm` (Funk 2022). `sin θ` is the angle between barbed-end tangent and membrane normal (+ŷ). On fire: marks barbed end CAPPED (state mutation only; topology preserved).

### Force partition
Network load `F_total` set externally via `set_network_load(F)`. Per-barbed-end / per-WAVE allotment = `F_total / N_free_barbed_ends` (D1 spec).

### Snapshot-based dynamic topology growth
The two elongation/branching Updaters APPEND particles + bonds to the live HOOMD `Snapshot` each batch tick (mirrors `crosslinkers.XlinkBondUpdater._extend_snapshot` pattern). Helper `_extend_snapshot_with_new_actins` handles particle / bond / angle merging into a fresh Snapshot.

## Sanity Gate §1–6 verification

| § | Coverage |
| --- | --- |
| 1 Dimensional | 3 tests (batch_dt, Funk threshold, 72° branch angle) |
| 2 Boundary | 4 tests (n_WAVE=0, Y_max>L_box/2, negative k_elong, D2 CFL auto-shrink) |
| 3 Topology | 4 tests (counts, plane y-coord, mother ℓ₀ below WAVE, anchor pair tags) |
| 5 Sign/sense | 2 tests (mother tangent direction, 72° branch geometry) |
| 6 Measurement | 4 tests (sim construct, demo BAOAB no-NaN, elongation event count, Funk abortive clamp) |
| Production opt-in | 1 SKIP (`H5_PRODUCTION=1`) |

## Test results

| Suite | PASS / SKIP / FAIL |
| --- | --- |
| `test_lamellipodium.py` (단계 1 new) | 17 PASS / 1 SKIP |
| `test_ku5x_lamellipodium.py` (단계 1 new — oracle utils) | 4 PASS / 3 SKIP (production opt-in) |
| Main scope regression | **276 PASS / 16 SKIP / 0 FAIL** |

H.3 baseline (단계 8 closeout): 259 PASS / 15 SKIP. H.5 단계 1 net: **+17 PASS + +4 PASS + +1 + +3 SKIP = +21 PASS + +4 SKIP, 0 regressions on H.3 / H.2 / H.1 / BAOAB**.

## Figures

All written to `ffn_sim/outputs/h5/figs/`. Per CLAUDE.md visualization-integrity rules.

| Filename | Purpose |
| --- | --- |
| `fig_h5_wave_plane_topology.png` | 3D scatter of WAVE particles (red squares on top membrane plane) + mother actin seeds (green dots, ℓ₀ below each WAVE). Membrane plane reference surface overlay. |
| `fig_h5_bell_evans_rates.png` | k_elong(F), k_b(F), k_cap(F) analytical curves on log y-axis. Funk abortive threshold annotated. |
| `fig_h5_updater_activity.png` | Live BAOAB-time-series: free barbed-end count (top) + cumulative elongation/branching/capping events (bottom) over 6 × 200 BAOAB steps. |
| `fig_h5_force_velocity_oracle.png` | Bieling 2016 v(F) = v₀·exp(−F·δ/kT) closed-form with 1/e force and v₀/e annotations. KU-5.2 reference for ±30 % production-gate band. |
| `fig_h5_pipeline_summary.png` | Text panel: 단계 1 deliverables + demo cell diagnostics + open production gates. |

## Open / next iteration

ALL remaining items are **PI sign-off / dedicated multi-hour production session** paths (same class as H.3 단계 8 close):

1. **`Cell.build` integration** — extend `cell.py:Cell.build` with `with_lamellipodium` option + `p_lamellipodium` parameter; wires WaveMembranePin + three Updaters into the unified cortex+xlinks+myosin+lamellipodium sim. Deferred to a follow-up session for proper integration testing (multiple force computes interleaving).
2. **KU-5.1 dendritic density** — multi-hour steady-state BAOAB run.
3. **KU-5.2 Bieling force-velocity** — multi-hour parameter sweep vs imposed cortex-side load.
4. **KU-5.3 Funk abortive** — multi-hour branching-rate-vs-force sweep around the 500 Pa per-WAVE threshold.

H.5 STATUS: 🟨 implementation freeze (단계 1 complete). PI sign-off + dedicated production session needed for ✅ DONE ratification.

## Phase 1 Main scope final state

H.1 + H.2 + H.3 + H.5 all at 🟨 / ✅ implementation freeze. **276 PASS / 16 SKIP / 0 FAIL** Main scope test count. H.4 (Sub session) on `phase1/h4-fa-clutch` (week 1-2 committed 2026-05-25 by PI as `6eddb11` on `phase1/h3-cortex` — orthogonal to Main work).

Remaining for Phase 1 closeout:
- H.3 production gates (ERM CFL, L_p FULL, KU-3.x, 3-way 60s, variable-length L_p)
- H.5 production gates (KU-5.1/5.2/5.3)
- H.4 Sub session weeks 3-4 (Sub-owned)
- H.7 single-cell integration (depends on H.1+H.3+H.4+H.5 all DONE)
