# H.7 Gate-B result + session record (2026-06-08→09)

**Branch** `h7/full-cell-integration`. **Verdict: REFUTE** (buckling is NOT the cortical-
tension lever) — but the session's larger finding is a **mesh confound that reframes the
whole γ-floor saga**, plus the conclusion that **the active-myosin cortical tension is
still floored even on the corrected (connected) mesh**.

---

## 1. The headline result (connected-mesh Gate-B, dt 1×, no-nucleus)

Full ×40 native run, suspended/rounded MCF7, FA off, turgor on, s_grip developed to ~0.37,
nonconv=0 (numerically sound). Plateau structural cortical tension γ_struct [mN/m]:

| condition | γ_struct | s_grip | bond/ℓ₀ |
|---|---|---|---|
| rigid_myoON  | 0.0780 | 0.37 | 1.000 |
| rigid_myoOFF | 0.0816 | 0    | 1.000 |
| relaxed_myoON  | 0.00136 | 0.37 | 0.984 (condensed) |
| relaxed_myoOFF | 0.00401 | 0    | 0.993 |

**γ_active = γ_struct(myoON) − γ_struct(myoOFF):**
- rigid:   0.0780 − 0.0816 = **−0.0036 mN/m**
- relaxed: 0.00136 − 0.00401 = **−0.0027 mN/m**

Both ≈ 0 (slightly negative, within the ~0.02–0.03 γ_rigid noise). CONFIRM band was the
Hosseini MCF7 IQR [0.18, 0.40]; PARTIAL ≥0.03. → **REFUTE.**

### Interpretation (honest)
- **The ~0.08 mN/m γ_struct is ~100% PASSIVE** (turgor + structural Lagrange tension):
  myoOFF=0.082 ≈ myoON=0.078, so **myosin contributes ~0 to cortical tension** even with
  s_grip developed to 0.37. The mid-run "~0.1, near the band" signal was passive turgor,
  NOT active myosin tension (hollow).
- **The active-myosin γ-floor PERSISTS on the connected mesh.** γ_soft (myosin active MOP)
  stayed ~2e-4 mN/m throughout — **~1000× below the MCF7 datum (0.27)**, same as Gate-A.
  Connectivity rescued only the PASSIVE structural channel (the rigid-Lagrange γ now carries
  turgor because the network percolates), NOT the active transmission.
- **Buckling makes it worse:** relaxed mode releases the load-bearing compression constraint
  → the rigid-Lagrange channel collapses (γ_struct 0.08 → 0.001, ~20× drop). So buckling is
  definitively not the lever.
- **The real imbalance:** myosin active stress (~2e-4) is **~400× smaller than the passive
  turgor tension (~0.08)** in this model. Generation works (s_grip walks to 0.37); the
  local contraction does not convert into spanning cortical tension and is dwarfed by turgor.

**Net:** Gate-A's active-γ REFUTE HOLDS even on the connected mesh + with buckling permitted.
The lever is neither buckling nor (for the active channel) connectivity — it is the
**myosin-active-stress vs turgor force balance**, which is the next thing to investigate.

---

## 2. The mesh confound (the session's most consequential finding)

`connected_mesh` was a **default-False opt-in** added 2026-06-04/07 (the cortex rebuild's
percolated mesh) for backward-compat. The native-vs-cupy PARITY driver
(`h7_native_fullcell_go`) omitted it (harmless for a parity test); **Gate-A's driver
(`h7_gate_a_native`) copied that pattern → measured γ on the FRAGMENTED mesh** (0 inter-
filament bridge crosslinks, giant-component ~7%); the new Gate-B native driver inherited it
again. A disconnected cortex cannot transmit spanning tension, so **Gate-A's γ-floor REFUTE
was partly a fragmented-mesh artifact** — confirmed here: on the connected mesh the PASSIVE
structural tension rose ~30× (3e-3 → 0.08). (The ACTIVE floor, however, is real — see §1.)

Root cause: a silent default-False flag that is effectively REQUIRED for any transmission/γ
measurement, with no guard, propagated by copy-paste from a parity test into γ measurements.

### Guards added (commit 35e1144)
- `build_baseline_cell(constrained=True)` now (1) **requires `connected_mesh` explicitly**
  (no silent default) and (2) **refuses the fragmented mesh** unless `allow_fragmented_dev=True`.
- `cell.py` logs the seeded giant-component / z / homeless on every connected build + warns
  if giant < 0.5.
- Fixed callers: `test_cortical_tension`, `h7_gate_a_native` → connected_mesh=True (23 tests pass).
- ⚠️ **Open gap:** build CFL guards check nucleus/membrane/ERM/turgor but NOT the
  myosin/xlink springs — which is why dt 10× silently violated the myosin-backbone CFL.

---

## 3. What was built this session (infrastructure, validated)

- **Relaxed-constraint (unilateral) M-SHAKE** + **τ_bend-EMA Euler buckling gate** in the
  cupy integrator (`constrained_baoab.py`) — tension side inextensible, compression side
  releases above F_crit on the τ_bend-averaged sustained load. 12 sanity-gate tests.
  Commit 5567194.
- **Native CUDA port** of the relaxed mode (`shake_relaxed_kernel` + orchestrator + module.cc
  + wrapper) — bit-parity vs cupy (1.7e-13), rigid-limit superset exact. Commit 79c1861.
- **Native Gate-B driver** (`h7_gate_b_native.py`, checkpointed, --aggregate). Commit 8ce6f55.
- **Sanctioned no-nucleus build** (`allow_no_nucleus`) — nucleus proven decoupled from
  cortical γ (0% across the KU-3.B2.1 band; removal Δ within seed noise). Commit 2b9c440.

These are reusable, but the relaxed-mode machinery answered a question whose verdict is "no"
(buckling is not the lever) — the effort would have been avoided by auditing the mesh first.

---

## 4. The dt / performance investigation (resolved)

- The constrained dt (6.98e-7 s = `0.001·τ_bend`) is **NOT over-conservative** — the ratified
  reconciler (`dt_reconcile`, B1/KU-1.26, safety 0.1) shows the binding CFL is the **myosin
  minifilament backbone spring** (k_backbone=4.24e-5 N/m, τ=9.2µs → dt_max=9.2e-7 s). The
  actual dt is ~0.76× of that, essentially correct (1.3× headroom only).
- **s_grip development is ~step-bound** (dt 10× gave only ~1.5×/step, not 10×) AND dt 10×
  violates the myosin-backbone CFL → under-developed/inaccurate s_grip. So **dt does not
  help the myoON conditions**; the cost is intrinsic.
- **Why production is slow / scale reality:** the run cost = (myosin contraction physical
  time ~2.5 s, from v0=0.2µm/s NMIIA Kovács 2003 + Hill) ÷ (myosin-backbone-CFL dt ~7e-7 s)
  ≈ 3.5e6 steps/condition. This is physics, not inefficiency. **5000-cell × 80 h at the
  fine-grained dt is impossible on any hardware** (~10¹¹ steps/cell × 5000); that long-
  timescale multicellular regime is the center-based-CBM (Layer-2) domain, fed by single-
  cell parameters via the scale-bridge. The fine-grained tool is a µs–seconds, single-cell
  computational microscope.
- Mitigation for future production: develop the operating point ONCE + checkpoint, then run
  perturbations from it; a semi-implicit/multi-timescale treatment of the stiff myosin
  backbone spring is the legitimate ~10× lever (unlike the nucleus/dt detours).

---

## 5. Open items / next direction
- **PRIMARY:** investigate the **myosin-active-stress vs turgor balance** — why does myosin
  contribute ~0 to cortical tension (active stress ~400× < turgor)? Candidates: myosin
  force/density too low, turgor too high, or local-contraction→global-tension conversion
  failing in this architecture. This is the real γ-floor question, distinct from buckling/mesh.
- Gate-A should be formally re-run/superseded on the connected mesh (its fragmented result is
  confounded for the passive channel; the active REFUTE stands).
- Close the build-guard gap (myosin/xlink spring CFLs not checked at build).
- A strategy brainstorm (separate session) on using the fine-grained platform beyond
  "CBM parameter supplier" — prompt drafted.

## Artifacts
- `outputs/h7/production/h7_gateB_cm_SUMMARY.json` + `h7_gateB_cm_{rigid,relaxed}_myo{ON,OFF}.cond.json`
- diagnostics: `h7_nucleus_contribution_check.py`, `h7_no_nucleus_dt.py`, `h7_dt_push_test.py`
- native parity: `native/ffn_hoomd_plugin/test_constrained_relaxed_parity.py`
