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

So the honest state (**substantially reframed by audit#7** — see the correction note; my earlier "interpenetration
is an unfixed pathological limitation" was over-pessimistic and wrong on the mechanism). (1) **Gap: fixed.** The
cortex no longer grossly floats (0.67 → **0.156 µm**, ~2 % of R = a physiological inter-cell cleft), driven by the
correct physics (conservative contact + differential interfacial tension + lit adhesion w_cs=2.85e-3). (2) **The
`pen_frac`≈2.1–2.6 "interpenetration" is the CONFLUENT-INIT geometry, not a contact failure.** Facts that force
this reading:

- It is present **from t≈0** (my N=400 log: pen 1.997 at step 500, before any contact dynamics settle) — the
  confluent Voronoi cells *share interfaces* by construction, and `pen_frac = max_node_into_neighbour / mean_edge`
  (a worst-single-node metric calibrated for *separated* aggregating cells) reads a shared interface as "overlap."
- It is **robust to the contact method**: penalty/tent gives pen 2.1, and the concurrent session's `--ipc`
  log-barrier (CCD, penetration-free *by construction*) gives pen **2.6** on the same confluent full-compartment
  (`2df86bb`) — a *better* contact does not lower it, because the overlap is geometric (built into the init), not
  dynamic. `--ipc`'s CCD prevents *new* penetration; it cannot remove the confluent init's shared-interface seam.
- It is **robust to the integrator** (matched-N control): explicit pen 2.227 ≈ implicit 2.107 at N=400 (both peak
  3.24). (An earlier N=100 implicit 1.321 *looked* like a win but is a size effect — pen_frac rises with cell
  count / shared-interface area; the controlled comparison must be at fixed N, and there it's a wash.)
- The run is **healthy**: V/V0=1.000, porosity ~0.10 (space-filling), A/A0=0.999, faceted (asph 0.054). If cells
  massively over-occupied space you'd see porosity ≪0 / V/V0≫1; you don't. A few sharp-Voronoi-vertex nodes poke,
  the bulk is properly space-filling.

So the **G2 gate (pen<0.3) is an aggregation-regime gate** (built for separate cells that must not touch);
**confluent space-filling tissue inherently has pen>0.3**, and the concurrent session's validated answer
explicitly accepts pen~2.6 at V/V0=1.0 as *known healthy equilibrium overlap*. The pressing contact is therefore
**real** (gap fixed, cells apposed + space-filling), and the residual pen is the confluent seam, not a pathology.
*(What I got wrong earlier: I read the confluent-geometry pen as a contact-dynamics failure and claimed it needed
a to-be-built IPC. Corrected by audit#7 below.)*

## Concurrent work (dcm/main) — corrected chain (per `d13372d`, `774c77c`, `2df86bb`)

A parallel session probed the *dynamic* compaction path and its result went through two retractions I initially
mirrored — both now corrected:

- `1c1f011` "compaction SOLVED" (aggregate Foty-Steinberg σ, Rg −26%) → `88bd2c8` "premature, native N=400 blows
  up (V/V0=0.39, pen=123), needs IPC" → **`d13372d`: it was a 10⁶× UNIT BUG.** `dP_agg` was `2.0e6·σ/R` but DCM
  positions are in **metres**, so the Laplace pressure was 3.3e8 Pa, not ~330 Pa — a numerical CRUSH that faked
  both the "compaction" (cell-collapse/overlap, not densification) and the "pen=123." Fix: `2.0e6→2.0`.
- With the unit fixed, a **physical σ-sweep (1–20 mN/m, the full Foty-Steinberg range)** shows aggregate σ is
  indistinguishable from baseline — **NO dynamic compaction** (`774c77c`, honest negative), all clean (pen≈0.3 on
  a *loose* start with `--ipc`). Dynamic loose→compact is not force-achievable at physical magnitudes (junction
  levers and aggregate σ both too weak vs turgor-incompressible cells; amplifying beyond lit range = forbidden
  magic-number tuning).

**Two corrections to my own prior claims fall out of this:** (a) my "their pen=123 converges with my pen=2.1 = same
contact failure" was **wrong** — theirs was a unit-bug crush, mine is the confluent seam; different causes. (b) my
"needs a to-be-built true log-barrier IPC" was **wrong** — `--ipc` **already IS** a genuine Li-2020 log-barrier
(`nearest_face_ipc_kernel`, `dcm_contact_implicit_warp.py:129-151`); no new IPC is needed, and it gives the same
confluent pen (~2.6) because the overlap is geometric.

**Convergence (the real one):** the concurrent session's endorsed answer (`2df86bb`, "the validated answer") is
**confluent-init full-compartment N=400 (cortex+turgor+nucleus+membrane tension), `--ipc`, implicit → V/V0=1.000,
porosity ~0.10 compact, faceted, pen~2.6 healthy** — i.e. it *converged on this reconstruction's approach*. My
build (penalty contact, pen 2.1) and theirs (`--ipc`, pen 2.6) are equivalent; the confluent-init full-compartment
is the shared, validated spheroid. Dynamic self-assembly remains FF-mature territory.

**Capstone attempted (full-compartment × compaction) — INTEGRATION UNSTABLE (honest negative).** Tried the
concurrent session's stable compaction config (loose voronoi gap 2.4, N=100, cadherin bundle-10, reach
cad_rbind=3.5 µm, contraction, aggregate σ=5 mN/m, conservative + implicit) PLUS my full compartments
(nucleus + membrane tension). Bonds form (reach fixed — `cad_rbind` is in µm, an audit-caught unit bug), but
the run **diverges: V/V0→148, cfl→9.6e8 (numerical blowup)**. The compaction stiff-force stack + the
compartment stiff-force stack destabilise TOGETHER (each is stable alone). So integrating the two workstreams
is NOT a trivial flag-combine — it needs dedicated stabilisation (a future task), not claimed as done.

## Verdict

All three violations are corrected: **GPU-native ✅ (A5000, verified no CPU regression), full VISIBLE +
mechanically-active compartment stack ✅ (nucleus force-field + membrane/cortex shells + cytoplasm + turgor,
all force-kernels launch every step; Chrome-headless render-verified), contact ✅ — floating fixed
(0.67→0.156 µm) and pressing (cells apposed + space-filling, porosity ~0.10, V/V0=1.000)**. All shown as
interactive HTML cell morphology (render-proof + cross-section PNGs committed). Hourly adversarial self-audit
active (cron `c3a5d4fb`); it worked as intended — it caught and corrected *my own* over-claims in both
directions: audit#3 corrected "contact clean/faceted-good" (→ be precise about the pen metric), and **audit#7
corrected the opposite over-pessimism** — I had framed the `pen_frac`≈2.1 as an unfixed pathological
interpenetration needing a to-be-built IPC. It is neither: `--ipc` (a genuine Li-2020 log-barrier) already
exists and gives the *same* pen (~2.6), because the pen is the **confluent-init geometry** (shared Voronoi
interfaces, present from t≈0, robust to both contact method and integrator) at a healthy V/V0=1.0 — not a
contact failure. The G2 gate (pen<0.3) is an aggregation-regime gate; confluent space-filling tissue inherently
exceeds it. HONEST bottom line: the cells are faceted, compartmented (mechanically), GPU-native, and properly
space-filling; the residual pen is the confluent seam, and the concurrent session's independently-validated
answer (`2df86bb`, confluent full-compartment + `--ipc`) converged on exactly this.
