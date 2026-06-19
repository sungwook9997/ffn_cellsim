# H.7 Native GPU optimization — sub-session summary (2026-06-07)

> ⚠️ **CORRECTION (2026-06-19, grounding-pass C1/C5).** Two numbers below drifted in prose and are corrected here. (1) **"44× integrator-only" is wrong** — the project's own same-day record (`H7_NATIVE_FULLCELL_GO_2026-06-07.md`) measured the same microbench at **39.47×** (native 5598 / cupy 142 steps/s); the only committed-JSON-verifiable figure is the **2.43× full-cell** (`h7_native_fullcell_go.json`). (2) The compartment-force **"parity 1e-15..1e-12" is unbacked** — the on-disk parity script compares native vs **CPU** (not cupy) at tolerance **rel 1e-8**, stdout-only, no committed artifact. Every speedup/parity number here is **GPU-only, gbook-measured, with NO committed build trace (.so/CMakeCache) and no CI** → unreproduced on dev. See the grounding-pass table + `project-rebuild-audit` memory.

Device-integrator sub-session, concurrent with the Lead's manifold/FA session. Goal:
make the constrained full-cell production fast enough for Gate-A/B (~2×10⁸ steps).
All deliverables are **additive + opt-in, GPU-only, CPU-fallback-safe**, in new files
or `native/` — no edits to the Lead's runtime/cell-physics files (zero collision).

## What landed (all validated + committed on `h7/full-cell-integration`)

| Lever | What | Result | Commits |
|---|---|---|---|
| **Integrator** | native C++/CUDA constrained L-M BAOAB (Fixman + M-SHAKE, 4 one-thread-per-chain kernels) | bit/tol parity; **~39.5× integrator-only** (⚠️ corrected from "44×" — see banner), **2.43× full-cell** (Lead-measured, JSON-verified), γ_soft+γ_rigid parity | cacb27b 0e2d490 24f12ad |
| **Compartment forces** | native `FFNRadialShellForce` ForceCompute (turgor/membrane/nucleus, sphere-equivalent radial laws, block-reduced) | parity 1e-15..1e-12; 3 forces **173 µs** (was 1860 cupy / 3700 cpu) | b659674 2a18203 (+ Lead wire a3a63ca) |
| **Gate-A verify** | everything-native on `build_baseline_cell` (real production path) | **6.0× vs cupy / ~4.6 d per 2e8**, nonconv=0, γ-parity (Lead's 5.86× agrees) | a8601a6 |
| **Binder enabler** | native `FFNAttachmentSpringForce` fixed-pool (head↔actin springs in device arrays; bind/unbind = `set_attachments` toggle, NO set_snapshot) | parity 3.4e-16 vs md.bond.Harmonic; toggle **17 µs** (vs **51.9 ms** set_snapshot); end-to-end **5.3×** vs set_snapshot binder | 3a26565 145ebb2 |

Refuted honestly: the **30× dt step-count lever** is stable but does NOT preserve γ
(rms T_bond 5.5× off) → not usable for the γ gates (14114c1). The **LJ** is the HOOMD
wall: Tree-nlist already optimal (Tree 1107 vs Cell 1340 µs), nucleus already excluded,
buffer/r_cut tuned — only a mesh-restricted excluded volume (research) remains (1cb0e5c).

## The everything-native step (verified, N=15,400, gate_a_profile b06ebd4)

```
1982 µs/step (505 steps/s, 6.0× vs cupy):
  binder updaters     1240 µs (62.6%)  global get/set_snapshot — set_snapshot = 51.9 ms/firing
  LJ excluded volume   366 µs (18.5%)  HOOMD wall (mesh = research)
  bond/angle/comp      236 µs (11.9%)  native
  integrator+triggers  139 µs ( 7.0%)  native
```

## The path forward (Lead's lane — the binder)

The **next ~2× lever is the binder** (62.6%): swap the myosin/xlink updaters from
`set_snapshot` bond mutation to the **`FFNAttachmentSpringForce` pool** — the enabler +
the pattern (`h7_binder_pool_demo`) are built + validated. Each binder firing becomes a
`set_attachments(head_tag, actin_tag, k, r0)` array toggle (17 µs) instead of a 51.9 ms
global rebuild; the force computes the springs every step (22 µs) on-device.

Projected (when the Lead wires it):
- **+ binder pool** (remove set_snapshot + the bond-list build) → **~8.4–16× / ~1.7–3.3 d**
  (the binding *decisions*, if vectorised like the demo's pool binder, are cheap).
- The binding-decision physics (Stam-Hocky myosin minifilament, catch-bond xlink) stays
  the Lead's lane; the spring force + the toggle are the provided generic enabler.

## Integration recipe (for the Lead)

```python
from ffn_hoomd_plugin import (
    NativeConstrainedBaoabUpdater,      # the integrator (swap onto a built constrained cell)
    NativeRadialShellForce,             # .from_{nucleus,membrane,turgor}(p, tag_range)  [Lead wired a3a63ca]
    NativeAttachmentSpringForce,        # pool_size=K; .set_attachments(head_tag, actin_tag, k, r0) per firing
)
```
All GPU-only with CPU fallback; harnesses: `scripts/h7_native_{gate_a_verify,gate_a_profile,
binder_pool_demo}.py`, `native/.../test_{shake,fixman,constrained,attachment_spring}_parity.py`.

## Status: my GPU-force/integrator lane is complete + verified. The remaining levers
are the binder wiring (Lead's myosin/xlink lane — enabler handed) and mesh-LJ (research).
