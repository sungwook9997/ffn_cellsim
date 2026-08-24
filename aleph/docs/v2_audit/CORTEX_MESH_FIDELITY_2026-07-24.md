# Cortex mesh / crosslink-density fidelity — PI catch (2026-07-24)

**Status: HONEST GAP surfaced for PI decision. Not yet fixed.** The PI flagged that the cortex crosslink
density looks too low for a "complete" cortex. Verified — the concern is valid. This memo separates the two
coupled issues, sources what can be sourced, and gives a recommendation. Literature values below are from the
author's domain knowledge; **DOIs/exact numbers are marked VERIFY — they are NOT yet in the Contract-Graph**
(a `tag_query` for cortex mesh size / crosslink spacing / density-per-filament returned **0 rows** — a real KB
sourcing gap).

## 1. What the model actually builds (verified in code)

- Cortex: **70,686 filaments**, **length ≈ 3 µm**, **segment `seg_um = 0.5 µm`** ⇒ 7 nodes/filament ⇒
  **494,802 actin nodes** (`assemble.py::_cortex_region` → `CORTEX_REGION` → `architecture_spec.py::CORTEX`).
- Crosslinks: `CrosslinkerSpec(hand=ALPHA_ACTININ, density_per_fil=20.0)` ⇒ **n_xl = 70,686 × 20 = 1,413,720**
  (`aleph/laws/weave.py:154` `n_xl = round(n_filaments · density_per_fil)`).
- **Provenance of the `20`:** chosen for **network PERCOLATION**, not sourced biology — the code comment says
  verbatim *"density_per_fil=20 keeps the network single-spanning under the −44% near-pair loss, 2.1× headroom"*
  (`assemble.py:181`). It is a **numerical spanning threshold**, and per the physiological-baseline HARD rule it
  should be the sourced in-vivo value. **UNSOURCED.**

## 2. The two coupled issues (they are NOT the same knob)

**Issue A — crosslink DENSITY (the PI's wording).** 1,413,720 crosslinks over 70,686 × 3 µm of filament ≈
**6.6 crosslinks/µm ≈ 1 per 150 nm** of contour. That is in the ballpark of a physiological crosslink spacing
but on the *sparse* side, and it is a percolation number, not a sourced one.

**Issue B — mesh RESOLUTION (the deeper cause).** The **`seg_um = 0.5 µm` discretization caps the pore/mesh
size at ~0.5 µm** regardless of density: crosslinks can only land on nodes that are 500 nm apart, so the 20
crosslinks/filament cluster at 7 nodes (~3 links/node to other filaments) rather than tiling a fine mesh. So the
effective cortex mesh is **~500 nm, set by the discretization, not the density.**

## 3. Physiological target (SOURCED 2026-07-24 via web + review; primary DOIs to enter KB)

Confirmed against Chugh & Paluch 2018 "The actin cortex at a glance" (J Cell Sci 131:jcs186254, PMC6080608)
and the primaries it cites:

| Quantity | Physiological value | Source | In KB? |
|---|---|---|---|
| Cortical actin **mesh / pore size** | **~30–200 nm** (EM tomography 50–200 nm cell-type-dependent, Morone; ~30 nm control → ~100 nm on perturbation) | **Morone 2006 J Cell Biol**; Chugh & Paluch 2018 review | **NO (0 rows) — must add** |
| **Crosslink spacing** | **≈ the mesh, 30–200 nm** ("might correspond to the distance between crosslinkers", Chugh 2018) | Chugh & Paluch 2018 | **NO — must add** |
| Cortex **thickness** h_cortex | ~100–1000 nm; ~190–200 nm typical (model ~200 nm ✓) | Clark 2013 (KB-3.1/3.5); Chugh 2018 | partial |
| Cortex **filament length** | **SHORT: formin ~1 µm + Arp2/3 ~100 nm mix** | Bovellan 2014 Curr Biol; Fritzsche 2016; Chugh 2018 | **NO — must add** |
| Areal actin density | ~100 filaments/µm² (model baseline ✓) | model KB | yes |

**Two gaps, both confirmed:**
- **Mesh (Issue B):** model **~500 nm** (0.5 µm discretization) vs physiological **~50–100 nm** ⇒ **~5–10× coarser**.
- **Filament length:** model **3 µm** vs physiological **~1 µm (formin) / ~0.1 µm (Arp2/3)** ⇒ **~3–30× too long**, and
  the model has NO Arp2/3 short-filament population (nucleator="formin" only) — a missing architecture, not just a
  number. (Bovellan: cortex is ~2/3 formin, ~1/3 Arp2/3 by mass.)

## 4. What a physiological-mesh cortex would cost (derivation)

To resolve a ~100 nm mesh the segment length must drop from **0.5 µm → ~0.1 µm** (5×), i.e.
**~2.5 M actin nodes** (494,802 × 5) for the cortex alone, plus a crosslink density raised to keep ~100 nm
crosslink spacing on the finer contour (**density_per_fil ≈ 30–60**, from the sourced mesh once VERIFY-ed).
Raising `density_per_fil` alone (e.g. 20 → 40) on the current 0.5 µm mesh increases connectivity but does **NOT**
fix the mesh size — Issue B needs the finer discretization. Memory: ~2.5 M cortex nodes is well within an A5000
(the current 494 k + membrane subdiv-8 655 k already fit).

## 5. Honest caveat on the GATE-A / GATE-B milestones

GATE A (converged resting baseline) and GATE B (event-driven dynamic tension) were validated **on this
0.5 µm-mesh cortex**. Those results establish the **solver convergence and the dynamic MECHANISM** (tension
emerges from binding events) — they are correct at the coarse-grained resolution. They do **not** yet establish
**physiological mesh fidelity**; the γ magnitude and any mesh-sensitive quantity carry this coarse-graining as a
caveat (on top of the already-flagged PI-GAP rate constants). I under-flagged this when reporting the milestones
— the PI is right.

## 6. Recommendation (PI decision)

1. **Source the physiological cortex mesh + crosslinker density** (Morone/Bovellan/Chugh/Fritzsche) into the KB
   (it is a genuine 0-row gap), so `density_per_fil` and `seg_um` become sourced, not percolation/convenience.
2. **Refine the cortex discretization** `seg_um` 0.5 → ~0.1 µm (cortex ~2.5 M nodes) to reach the physiological
   mesh — a **significant** fidelity increment: re-run the GATE-A convergence and GATE-B dynamic slice at the
   finer mesh. This is the real "complete cortex."
3. Until then, treat the current cortex as an **explicitly coarse-grained** baseline (label it as such in the
   viz + reports), not a physiological-mesh cortex.

**This is a PI-scoped decision** (it re-opens the validated GATE-A/B at 5× the node count), so it is surfaced
rather than done unilaterally. Ready to execute step 1 (sourcing) + step 2 (a finer-mesh build + re-validation)
on approval.

## 7. Derived fine config + the configurable cortex geometry (2026-07-24)

The cortex geometry is now **configurable** (was hard-coded in `architecture_spec.CORTEX`). Three additive
`CellConfig` knobs in `aleph/components/incumbent/assemble.py` thread onto a COPY of `CORTEX_REGION.arch` via
`dataclasses.replace` (the module-level `CORTEX` constant is never mutated):

| knob | baseline (= today) | what it drives (downstream `ff.weave.weave`) |
|---|---|---|
| `cortex_seg_um` | 0.5 µm | mesh/pore resolution — `beads_per_filament = round(length/seg)+1` |
| `cortex_length_um` | 3.0 µm | representative filament contour length L |
| `cortex_density_per_fil` | 20.0 | crosslinks — `n_xl = round(n_fil·density)`; contour spacing = `length/density` |

**The defaults reproduce today's arch bit-for-bit** (verified: `_cortex_region(70686).arch` equals the old
replace-only-`n_filaments`+`R_um` result; test `test_default_cortex_region_is_bit_identical_to_legacy`).

### 7.1 Derivation (no magic numbers — from the sourced §3 mesh)

The single sourced quantity is the cortical **mesh ≈ crosslink spacing ≈ 50–100 nm** (Morone 2006; Bovellan
2014; Chugh & Paluch 2018, §3). Take the midpoint **mesh_target = 75 nm = 0.075 µm**. Two geometric identities
close the config for the R=7.5 µm, 70,686-filament, ~100/µm² cortex at the current representative L = 3.0 µm:

- **Mesh resolution = the discretization.** Crosslinks can only land on nodes `seg_um` apart, so the pore size
  is capped at `seg_um`. To resolve a 75 nm mesh ⇒ **`seg_um = mesh_target = 0.075 µm`.**
  Then `beads/fil = round(3.0/0.075)+1 = 41` ⇒ **nodes = 70,686 × 41 = 2,898,126.**
- **Crosslink contour spacing = mesh (Chugh 2018).** spacing `= length/density`, so
  `density = length / mesh_target = 3.0 / 0.075 =` **40** ⇒ `n_xl = round(70,686 × 40) =` **2,827,440.**

### 7.2 Recommended fine config + cost

```python
CellConfig(
    n_filaments=70686,
    cortex_seg_um=0.075,          # = 75 nm mesh_target (Morone/Bovellan/Chugh §3) — DERIVED, not tuned
    cortex_length_um=3.0,         # representative L kept (the short-filament fix is a separate architecture gap)
    cortex_density_per_fil=40.0,  # = length/mesh_target = 3.0/0.075 — holds ~75 nm crosslink spacing
)
```

| quantity | coarse baseline | **fine (0.075 µm)** | ratio |
|---|---|---|---|
| beads / filament | 7 | 41 | 5.86× |
| cortex actin nodes | 494,802 | **2,898,126** | 5.86× |
| crosslinks (`n_xl`) | 1,413,720 | **2,827,440** | 2.0× |
| mesh / pore size | ~500 nm | **~75 nm** | — |

**Rough peak bytes.** The mechanical actin state (pos vec3d 24 B + force 24 B + a few int32 topology arrays,
plus the crosslink arrays 2×int32+k+rest per xl) scales ~linearly with node + crosslink count, so the cortex
mechanical footprint grows ~**5.9×** on nodes / ~**2×** on crosslinks. For calibration the memo §4 anchor is
subdiv-6 full-native ≈ **592 MB** on the A5000 16 GB with the 494 k coarse cortex; the fine cortex adds ~2.4 M
nodes (~a few hundred MB for the actin arrays) — total **well within the A5000 16 GB** alongside membrane
subdiv-8 (655 k verts), nucleus, myosin heads, and the fluid grid. **A5000-feasible.**

### 7.3 Reaching ~50 nm (honest note)

The 50 nm lower bound needs **`seg_um = 0.05 µm` ⇒ beads/fil = 61 ⇒ 70,686 × 61 = 4,311,846 nodes (~8.7×)**,
with `density = 3.0/0.05 = 60 ⇒ n_xl = 4,241,160`. That is a **~9× actin-node** build; the mechanical state
alone is still ≲ 1–2 GB, but stacked with membrane subdiv-8 + nucleus + explicit myosin heads + the Biot field
grid it starts to **pressure the 16 GB** — it should be memory-accounted before committing, and if tight is the
case for a larger / multi-GPU device (per the hard rule: never lower biological density to fit memory). **The
practical fine config is `seg_um ≈ 0.075–0.1 µm` (75–100 nm mesh, comfortably A5000-feasible); 50 nm is a
stretch validation resolution, not the first fine build.**

### 7.4 Caveat carried forward

`cortex_length_um` stays 3.0 µm in the recommended config: the **filament-length gap (§3: physiological ~1 µm
formin + ~0.1 µm Arp2/3, and NO Arp2/3 short-filament population)** is a deeper *architecture* fix (a missing
nucleator population), not just this number. The fine-mesh config fixes **Issue B (resolution)**; the short-
filament architecture remains a separate PI item.

## 8. TAG/KB cross-check (2026-07-24) — important corrections

A TAG query against the Contract-Graph refines §3:

- **Native mesh is ~30 nm, not 50–100 nm.** KB-3.1 / KB-3.18 state verbatim: *"NATIVE mesh ~30 nm [Bovellan2014];
  the ~100 nm figure is the PERTURBED mDia1/Arp2-3-inhibited value, NOT native"* (30–50 nm native). So the model's
  ~500 nm mesh is **~10–16× coarser** than native — the gap is larger than the earlier ~5–10×. The physiological
  target segment length is **~0.03–0.05 µm**, not 0.075 µm.
- **The coarse-graining is DOCUMENTED, not hidden.** KB-3.18: *"PRODUCTION count N=70,686 cortical = 100/µm² ×
  4π(7.5µm)² (coarse REPRESENTATIVE filaments; magnitude order-correct vs native 30–50 nm mesh)"* + see
  `docs/v2_audit/CORTEX_COUNT_CROSSCHECK_2026-07-21.md`. So the 70,686 / 0.5 µm cortex is an **acknowledged
  representative coarse-graining** (magnitude order-correct), NOT a silent convenience bug like the turgor. The PI's
  ask is to push it toward native mesh fidelity — the right direction, from a documented baseline.
- **Cortex filament length is an open "PI A/B/C decision"** (KB-3.18), consistent with §2/§3: native cortical
  filaments are short (~hundreds nm) + high-turnover; the model's 3 µm is the "general linear-actin range", flagged
  as a PI choice, not a sourced cortex value.

**Cost of the NATIVE target** (30–50 nm ⇒ seg_um 0.03–0.05): ~7–17× the 494 k nodes (≈3.5–8.5 M cortex nodes),
which pressures the A5000 16 GB alongside membrane-subdiv8 + nucleus + heads + Biot grid — likely the
larger/multi-GPU case (never lower density per the HARD rule). So the practical path is a **resolution ladder**:
the derived seg_um=0.075 (2.9 M nodes, A5000-feasible) as the first fine build to test mesh-sensitivity, then the
30–50 nm native as the stretch/multi-GPU validation resolution. The turgor (§ audit) is separate and more urgent
because it's a *silent* HeLa proxy, not a documented coarse-graining.

## 9. NATIVE mesh-convergence result (2026-07-24) — the physiological mesh hits a solver+memory WALL

Ran the fine config natively (gbook A5000 16 GB) after gbook came back (root cause of the "down" was the Mac's
Tailscale being stopped, not gbook — fixed with `tailscale up --accept-routes`). The CLEAN turgor baseline
(no myosin, subdiv=8) at the FINE cortex (`seg_um=0.075`, `density=40`, 2,898,126 actin nodes / 3,569,158 total):

| mesh | actin nodes | pathA (matrix-free, weak) | pathB (explicit A_c, strong) |
|---|---|---|---|
| coarse ~500 nm | 494,802 | **converges 0.1398** ✓ | — |
| fine ~75 nm | 2,898,126 | **plateaus max\|PF\|=3.46** ✗ (it0 4.75→3.46, stuck t=0) | **OOM** ✗ (85 MB alloc fails; explicit A_c > 16 GB) |

**The physiological mesh is NOT a free config change — it is a real solver+memory wall:**
- Finer discretization worsens the conditioning (classic FEM: κ↑ as h↓), so the weak matrix-free solver **plateaus
  ~25× above the coarse** (3.46 vs 0.14) and the line search stalls (t=0).
- The strong solver (explicit crosslink-weighted BSR `A_c`) **does not fit 16 GB at 2.9 M nodes** (OOM).
- So the coarse mesh's clean convergence (0.14) is partly *because* it is a well-conditioned coarse-grained
  representation; refining toward native exposes the conditioning + memory cost.

**Implication (honest, answers the PI's mesh question):** reaching physiological mesh fidelity requires ONE of:
(a) a **memory-efficient strong solver** for the fine mesh (matrix-free geometric/algebraic multigrid, not an
explicit `A_c`) — a solver-development effort; (b) **more GPU memory** (multi-GPU / larger card — the HARD-rule
response to native memory pressure, never lower density); (c) a **resolution LADDER** — the finest mesh that both
fits AND converges on the A5000 (testing ~200 nm now). The coarse baseline's GATE-A/GATE-B results stand at their
resolution; the physiological-mesh target is gated on (a) or (b).

### §9 update — fine-mesh GATE-A convergence is IMPRACTICAL on the A5000 (all solvers)
Retested with the memory fix: `disable_fiber_block=True` (frees 2.85 GB) + `--membrane-subdiv 6`, fine 75 nm.
- pathA weak (OFF): plateaus max|PF|≈4.4 (even higher than the subdiv-8 3.46 — dropping the fiber block weakens
  the OFF preconditioner).
- pathB strong (explicit A_c): impractically SLOW — the 200 nm (1.13 M-node) ON run was stuck at it 0 for >24 min
  (>200 s/outer for the 800-iter block-Jacobi CG on the explicit A_c); killed.
- pathA + fq-coarse (matrix-free, fq-iters 200): also impractically slow — the fine (2.9 M) ON ran 26 min without
  completing outer iterations; killed.
**Verdict: the physiological ~75 nm mesh does not converge PRACTICALLY on one A5000 with the current solver stack
(weak=plateau, strong=OOM-or-too-slow). The resolution LADDER shows ~200 nm is the practical A5000 rung (weak
solver, max|PF|≈0.21, ~2 min). Reaching physiological mesh needs a genuinely fast+memory-efficient multigrid
(research) OR multi-GPU. The coarse/200 nm GATE-A/B results stand at their resolution.**

## 10. Fine GATE-B γ (dynamic) — mesh has a MODEST effect; coarse γ is roughly valid
Ran the DYNAMIC GATE-B (NMII binding events → emergent γ) at fine 75nm vs the coarse baseline (explicit inner
solve, fits ~2 GB — no OOM, unlike the static GATE-A). Steady state:

| | coarse ~500nm | fine ~75nm |
|---|---|---|
| bound% | ~85% | **98.2%** |
| γ_total | ~3.7 pN/µm | **4.41 pN/µm** (~19% higher) |
| max\|PF\| | ~0.5 | 4.19 (not equilibrium — same fine-mesh conditioning) |

**Answer to "does mesh change emergent tension": MODESTLY.** The finer mesh (75nm ≈ the ~50nm NMII capture reach)
puts an actin node within reach of nearly every straddle-placed head, so bound% rises 85→98% and γ rises
3.7→4.41 (~19%). So the emergent tension is set MOSTLY by the myosin mechanism (duty / f_stall / k_xb), with a
modest mesh sensitivity via the bound fraction — the **coarse GATE-B γ (3.7) is in the right ballpark, not a
gross coarse-mesh artifact.** Caveat: the fine GATE-B does not reach mechanical equilibrium (max|PF| 4.19 vs
coarse 0.5 — the same fine-mesh conditioning the static GATE-A exposed), so the fine γ 4.41 is a rough estimate
at a non-equilibrium state, not a converged value.
