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

# Cortex percolation fix — plan (2026-07-23)

**Author:** analysis/diagnosis session (NOT the gbook production/main session — this is a handoff).
**Status of the code change:** `density_per_fil 1.0 → 20.0` is **edited in the working tree** (`ff/architecture_spec.py`
CORTEX spec) but **NOT committed** (unprompted-commit rule). Everything below C is for the main session to execute
on the A5000, single-commit-aligned.

---

## TL;DR

The resting-baseline convergence rollback and the "myosin heads float 0.13 µm off actin" finding are **both symptoms
of a structurally broken cortex**, not of myosin parameters or missing sourcing. The native dump shows the cortical
F-actin network is **fragmented into 11,430 pieces (giant component 79.6 %)** because crosslinks were placed at
**1 per filament** (`density_per_fil=1.0`), a ×40 coarse-graining relic that puts the network at **mean degree 2 —
below the ln(N)=11.2 percolation threshold**. The fix is a single sourced value: **`density_per_fil = 20`** (DERIVED
from the measured native mesh ξ≈30 nm), which raises mean degree to 40 → single spanning network (predicted giant
100 %). Cytosim is the construction blueprint (connectivity **emerges** from binding kinetics at physiological
crosslinker density, ~8–16 xl/filament) — but it only covers the cortex network; the membrane/turgor/fluid that make
the cell actually deform/bleb are separate components we already have and couple to via connectors.

---

## 1. Diagnosis — why the cortex is not trustworthy (native dump, `outputs/ac/cell_assembled/assembled_state.npz`)

Quantified structural audit of the committed native build (494,802 actin nodes / 70,686 fibers):

| # | Defect | Measured | Should be |
|---|---|---|---|
| 1 | **Zero-thickness shell** | every actin node at exactly R=7.400 µm, std=0 | ~0.2 µm thick (KB-3.1) |
| 2 | **Monodisperse length** | all fibers 3.0 µm / 7 nodes, identical | exponential (KB-3.18) |
| 3 | **Myosin heads splayed off actin** | head→actin-node median 0.133 µm; only 15.9 % bind within 0.05 µm capture | heads bound on the filament |
| 4 | **Crosslinks ultra-sparse** | 70,686 xl = **1/filament**, 24.7 % of nodes carry any xl | ~10–30/filament (physiological) |
| 5 | **KILLER — fragmented network** | **11,430 connected components**, giant 79.6 %, 9,478 fragments <10 nodes | 1 spanning component |

Defects 1, 2 are produced by `ff/weave.py` defaults (`overlap_free=False` → `cortex_thickness_um→0.0`;
`length_dist="mono"`). Defects 3, 4, 5 are the crosslink-density relic (below). The dump is **not stale** — it is
exactly what the current `weave()` path emits.

**These explain the resting blocker directly.** The 23d residual decomposition is per-fiber translation 35 % +
fiber internal 44 %: a fiber with 1 crosslink and unbound heads is a near-free floppy mode (near-zero operator
eigenvalue) → wrecks CG conditioning. "Convergence #2 conditioning stall" and "under-percolated network" are the
**same defect**; three preconditioner families were exhausted because you cannot precondition a structurally
pathological operator — you fix the structure.

## 2. Root cause — `density_per_fil=1.0` is a ×40 coarse-graining relic

- `ff/weave.py:154`: `n_xl = round(n_filaments * density_per_fil)`. Cortex default was `1.0` (`architecture_spec.py`
  CORTEX, note said "α-actinin **1:1**" — a convenience, not a sourced crosslinker density).
- Filament-connectivity graph mean degree = `2·n_xl/n_fil = 2.0`. Erdős–Rényi giant component at mean degree 2
  solves `S=1−e^(−2S)` → **S=0.7968**, matching the observed **79.6 %** exactly. The fragmentation is a
  **mathematical certainty of the density**, not a build bug.
- Single spanning component needs mean degree ≥ **ln(70,686)=11.2**; we were at 2.0 (1/5.6 of threshold).
- The repo's own docs already recorded that at the molecular reach ε=60 nm the network fragments (giant 1.5 %,
  z=0.20) and only the coarse `mesoscale_reach=√(A/n)` percolates (99 %) — so "spanning mesh 99 %" was **always a
  coarse-graining artifact** (`H7_SURFACE_MANIFOLD_EXPLICIT_CORTEX_2026-06-07.md`). At native it lands at 79.6 %.

**Feasibility — the network is NOT geometrically starved (native dump, cross-fiber near-pairs vs reach):**

| reach µm | cross-fiber pairs | /fiber | max mean-deg | giant%(ER) |
|---|---|---|---|---|
| 0.06 (molecular ε) | 2,011,475 | 28.5 | 54 | 100 |
| 0.10 (current reach) | 5,590,579 | 79.1 | 145 | 100 |
| — current build (cap 1.0) | 70,686 | 1.0 | 2.0 | 79.7 |

At the true molecular ε=60 nm there are 28.5 near-pairs/fiber available; `density_per_fil=1.0` **discards 96.5 % of
them**. Note: at native, `mesoscale_reach = √(4πR²/n) = 0.10 µm` ≈ filamin length scale (150 nm), so the reach is
already ~molecular at native and needs no change — the crutch was only pathological at ×40 (~1000 filaments → reach
0.84 µm).

## 3. Cytosim is the construction blueprint (and its scope limits)

Cytosim (Nédélec–Foethke 2007; our FF mechanics already port NF2007, same s·µm·pN units, and it is our bending-energy
parity oracle) builds a cortex by **self-organization**: filaments placed/grown, then **connectivity EMERGES from
crosslinker/motor binding kinetics** (Hand binds to the *closest point* of any filament within ε; Bell force-dependent
unbinding). Belmonte–Leptin–Nédélec 2017 used **12,000–64,000 connectors for 1,500–4,000 filaments = 8–16/filament**
(mean degree 16–32, well above threshold). Every one of our 5 defects is a symptom of having *placed a static network*
instead of letting it emerge.

**But Cytosim is a cytoskeleton engine, not a cell engine.** Its object model has **no membrane** (`Space` is a
confinement potential, usually rigid), no turgor/osmotic pressure, no poroelastic fluid/FSI, limited fiber-fiber
steric, no nucleus/LINC/IF/MT as first-class. The cell-shape-level behaviors (dent, crumple, **bleb**) are a coupled
**network ⊗ membrane ⊗ fluid ⊗ pressure** problem. We already have those components (Helfrich membrane, ERM clutch,
turgor, Biot poroelasticity, bleb machinery). So the architecture is: **build the cortex COMPONENT the Cytosim way,
keep coupling it to our membrane/fluid/nucleus components via `ac/engine/` connectors** — which is exactly what the
component/connector graph already targets. We cannot *run* Cytosim (CPU C++ vs our Warp-CUDA-only rule); we port the
paradigm, most of which already exists (`ac/weave/crosslink_kmc_warp.py`, Stam-Hocky heads w/ Hill+Bell, NF2007).

Emergent connectivity is also **forward-compatible with blebbing**: a placed static mesh cannot nucleate a new cortex
under a bleb, but a binding-kinetics cortex can reform — where Cytosim shines.

## 4. The ratified change (PI 2026-07-23)

`ff/architecture_spec.py` CORTEX: **`density_per_fil = 1.0 → 20.0`** (edited, uncommitted).

**Derivation (satisfies the Magic-Number Block — derived, grid-invariant, not tuned to a gate):**
`crosslinks/filament = L / δ_filamin = 3.0 µm / 150 nm = 20`, filamin being the main cortical crosslinker (KB-3.18)
at the measured native mesh ξ≈30 nm (Bovellan 2014, KB-3.1/3.18). Self-consistency: actin geometry (100 fil/µm²,
L=3 µm, h=200 nm) gives ρ=1500 µm⁻² → ξ=1/√ρ=26 nm ≈ measured 30 nm ✓. Three independent cross-checks converge:
ξ-derivation ~20; Cytosim/Belmonte 8–16; geometric ≤28.5 available at ε=60 nm.

**Verified prediction (non-CUDA, spec arithmetic):** n_xl 70,686 → **1,413,720**; mean degree 2.0 → **40**; ER giant
79.7 % → **100.00 %**. Geometrically suppliable (5.59M pairs available at reach 0.10 µm ≫ 1.41M needed).

**Honesty caveat:** 20 is the physiological value *regardless of whether it fixes convergence*. If D still fails,
that is a separate finding — do NOT un-source the density to chase the gate.

**Coupling caveat:** 20 is for the **L=3 µm monodisperse** build. When length → exponential/shorter (Tier 3-E), total
crosslinks are fixed by crosslinker concentration (not by how filaments are chopped), so `density_per_fil` must be
**re-derived** (shorter filaments ⇒ more filaments ⇒ fewer xl each at fixed mean degree).

## 5. Execution plan (main session, A5000, single-commit-aligned gbook)

**Gates:** density value PI-ratified ✓. gbook must be single-commit-aligned (Codex `SYNC NOW`) before the GO probe —
the handoff notes gbook is currently partial-synced (driver/inner_mechanics on the 30 pN side-branch → acceptance
unreliable).

- **C — rebuild + re-audit structure.** Commit the density edit. Rebuild native cortex (`dump_state.py`). Re-run the
  structural audit (§6): confirm **single spanning component** (mean degree ~40, giant ~100 %, fragments→0). No code
  beyond the density value; reach unchanged (already molecular at native).
- **D — resting convergence probe (THE decisive test).** `--from-resting --overlap-free-cortex
  --membrane-subdivisions 6` on the percolated cortex. Check whether the per-fiber-translation/floppy residual
  component (35 %+44 %) collapses and whether the 0.21 pN projected-force gate is reached (`inner_converged=True`,
  `outer_accepted=True`, no rollback). This validates or refutes "structure defect = convergence defect."
- **E — remaining structural fidelity (parallel, lower priority).** Make `overlap_free=True` the native default
  (restore 0.2 µm thickness; re-measure near-pairs on the overlap-free build — margin 54 vs threshold 11.2 is large
  so still percolates) + turn on the exponential length distribution (`length_dist="exponential"`, KB-3.18), then
  re-derive `density_per_fil` per the coupling caveat.
- **F — Cytosim-paradigm proper (after C/D validate).** Promote static force-free placement to emergent binding:
  wire `ac/weave/crosslink_kmc_warp.py` at physiological crosslinker density, and bind heads to the nearest actin
  point via kinetics (fixes defect #3, the 0.13 µm splay). This is a wiring job, not from-scratch.

## 6. Parity / test implications (flag before committing)

- **γ-floor bit-parity breaks.** The historical γ-floor cell used `density_per_fil=1.0`; density 20 changes every
  cortex build. `tests/ff/test_weave.py` asserts statistical parity with `gamma_floor.build_crosslinked_cortex` — it
  will need its expected n_xl updated, or the change gated behind a flag. To reproduce a legacy γ run, pin
  `density_per_fil=1.0` explicitly. Most other tests (`test_connected_mesh`, `test_cell_full`, …) use their own
  `dynamic_crosslinkers.n_xl` configs and are unaffected.
- **Cost.** n_xl ×20 (70,686 → 1.41M crosslink springs) raises build/step cost and memory. A5000 16 GB should absorb
  it (membrane subdiv-7 already hits 837k nodes / 752 MB), but log peak GPU bytes on the C rebuild.

## 7. Reproduce (structural audit — CPU, dev Mac)

```python
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
d = np.load("aleph/outputs/ac/cell_assembled/assembled_state.npz", allow_pickle=True)
pos = d['pos_post'].astype(float); n_actin = int(d['n_actin']); fo = d['fiber_offsets'].astype(int)
xl = d['xl'].astype(int)
a=[];b=[]
for i in range(fo.shape[0]-1):
    s,e=fo[i],fo[i+1]
    if e-s>=2: idx=np.arange(s,e-1); a.append(idx); b.append(idx+1)
sa=np.concatenate(a); sb=np.concatenate(b)
ea=np.concatenate([sa,xl[:,0]]); eb=np.concatenate([sb,xl[:,1]])
g=coo_matrix((np.ones(ea.shape[0]),(ea,eb)),shape=(n_actin,n_actin))
ncomp,lab=connected_components(g,directed=False); sizes=np.bincount(lab)
print("components", ncomp, "giant%", 100*sizes.max()/n_actin)   # baseline: 11430, 79.6%; target after C: ~1, ~100%
```
