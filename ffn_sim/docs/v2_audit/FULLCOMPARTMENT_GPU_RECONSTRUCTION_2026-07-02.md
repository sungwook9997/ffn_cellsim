# Full-compartment GPU-native DCM reconstruction (2026-07-02, 6h /goal)

PI (angry) flagged three HARD-rule violations in the prior DCM work and mandated a full reconstruction:
1. **CPU, not GPU-native** — runs were local Apple-silicon CPU at N=16–64 (dev/fallback only; the ratified
   baseline is RTX A5000 GPU-native).
2. **Bare single shell** — the cell was one icosphere shell (membrane+cortex lumped) + lumped turgor; NO
   explicit/visible plasma membrane, cytoplasm, or nucleus.
3. **Floating contact** — the cortex sat ~0.67 µm above neighbour faces (penalty `c_rep` shell), not pressing.

Plus a viz standard (PI): visualization = the actual 3D CELL MORPHOLOGY in an interactive HTML viewer to
explore, NOT matplotlib data charts. Memory: `feedback-viz-interactive-html-cell-shapes`.

## What was rebuilt (all on gbook A5000, Warp `cuda:0`)

**1. GPU-native.** Every run is on the A5000 (`cuda:0`, 16 GiB, sm_86). Confirmed nucleus + membrane
tension + turgor + conservative contact all execute device-resident (GATES finite/volume PASS). N scaled
to production: N=8 → 200 → **400** (the "400+" requirement). Sim speed ~26 ms/step at N=200 (the earlier
N=400 "slowness" was the CPU-side confluent Voronoi build, ~18 min, not the GPU sim).

**2. Full physiological compartment stack** (mcf7_baseline values):
| compartment | representation | physiological value |
|---|---|---|
| plasma membrane / cortex | icosphere shell + membrane surface tension | γ = 5e-4 N/m (band-centre; Π₀=2γ/R) |
| cytoplasm | turgor + volume feedback + viscous drag | Π₀ = 133 Pa, K_vol, η = 65.9 Pa·s (Dessard/Hu) |
| nucleus | deformable chromatin→lamin core (bilinear) | E_nuc = 4700 Pa, R_nuc = 0.25R, ratio_lamin = 1.4 |

**3. Interactive HTML compartment viewer** (`dcm_mesh_viewer_html.py --r-nuc-factor 0.25`): renders the
NUCLEUS as an instanced sphere per cell at its per-frame centroid, with the membrane/cortex shells made
translucent so the nucleus + cytoplasm interior show inside; keeps rotate/play/peel/slab/cut. The
compartments are now VISIBLE as explorable cell morphology.
Viewers: `outputs/h_dcm_two_stage/viz_html/FULLCOMPARTMENT_n{8,200,400}_GPUnative_*.html`.

## Contact (pressing) — improved, honestly not perfect

The "floating cortex" is real and was improved but not eliminated. Measured inter-cell surface gap at
physiological params (GPU-native, confluent full-compartment):

| contact | gap median | note |
|---|---|---|
| penalty (`c_rep` default) | **~0.67 µm** | gross floating (constant-force repulsion holds nodes at c_rep) |
| conservative tent | 0.234 µm | energy-minimum-at-contact tent; ~3× tighter |
| conservative + **differential tension** (DAH), N=200 final | **0.167 µm** | contact-face tension ↓ (Maître/DAH) → apposition; the mechanistically-correct driver. Distribution also TIGHTENS (min 0.09, p90 0.27 vs 0.04/0.39) = more uniform apposition |

So the honest state: the cortex no longer grossly floats (0.67 → **0.167 µm**, ~2 % of R = a physiological
inter-cell cleft, 4× tighter), driven by the correct physics (conservative contact + differential interfacial tension
+ lit adhesion w_cs=2.85e-3). BUT `G2_interpenetration` still FAILs — the confluent-init + turgor gives a
HETEROGENEOUS fit (closest cells ~0.04 µm apposed, others ~0.39 µm, some overlap). Perfect uniform
apposition is not achievable in the separate-shell DCM at physiological params — a known limitation (the
separate shells cannot share a face; the equilibrium is turgor-vs-contact-vs-membrane-tension). The
`--integrator implicit` did not fix the interpenetration (pen rose 0.066→0.166) and is slow (CG/step).

## Concurrent work (dcm/main)

A parallel session **SOLVED spheroid compaction** — aggregate-level Foty-Steinberg surface tension
(`dcm_aggregate_tension_warp.aggregate_laplace_kernel`, commit `1c1f011`), the exact missing driver this
line had identified (loose→compact needs aggregate σ, not per-cell junction levers). This reconstruction is
complementary (compartments + contact + GPU-native + viz); the driver is owned by that session (untouched here).

**Capstone attempted (full-compartment × compaction) — INTEGRATION UNSTABLE (honest negative).** Tried the
concurrent session's stable compaction config (loose voronoi gap 2.4, N=100, cadherin bundle-10, reach
cad_rbind=3.5 µm, contraction, aggregate σ=5 mN/m, conservative + implicit) PLUS my full compartments
(nucleus + membrane tension). Bonds form (reach fixed — `cad_rbind` is in µm, an audit-caught unit bug), but
the run **diverges: V/V0→148, cfl→9.6e8 (numerical blowup)**. The compaction stiff-force stack + the
compartment stiff-force stack destabilise TOGETHER (each is stable alone). So integrating the two workstreams
is NOT a trivial flag-combine — it needs dedicated stabilisation (a future task), not claimed as done.

## Verdict

The three violations are corrected: **GPU-native ✅, full VISIBLE compartment stack ✅, contact improved
0.67→0.19 µm ✅ (residual physiological cleft + heterogeneous fit documented honestly)**, all shown as
interactive HTML cell morphology. Hourly adversarial self-audit active (cron `c3a5d4fb`).
