# Phase C — node-FACE contact 切り分け + audit cleanup (2026-06-23)

Branch `h7/compartment-platform`. Continues `PHASE_C_SPHEROID_AGGREGATION_RESEARCH_2026-06-22.md`
and the 9ec3ee1 "2-cell contact REVEALS cells overlap as rigid spheres" claim.

## ⭐ HEADLINE — the "rigid-sphere overlap, do NOT flatten" alarm was a MEASUREMENT/RENDER ARTIFACT

PI asked to first 切り分け (separate) the two failure modes for the open contact problem:
is the node-FACE contact (a) **failing to prevent cell-cell overlap**, or (b) **missing a
flattening mechanism** entirely? Built a proper 2-cell diagnostic (`scripts/twocell_overlap_diag.py`:
cross-section + oblateness + midplane-crossing + point-in-other-mesh) and a long-settle rep sweep
(`scripts/twocell_rep_sweep_long.py`). **Answer: NEITHER failure is real.**

**The 2-cell contact DOES flatten cells into a clean flat junction AND DOES prevent macroscopic
overlap.** The prior session's opposite conclusion came from (1) a **3D trisurf render** that draws
each cell as a smooth interpolated sphere — which *hides* the flat contact face — and (2) an
**under-settled run** (cells still barely-touching/round near the initial pack). The cross-section
view + the oblateness metric (both render-independent) disprove it. This is the **3rd
measurement/visualization artifact** in this lineage (cf. the retracted "spheroid SPREADS
A/A0→1.94" 7-cell-peeling artifact; the retracted "rigid marbles/interpenetration" flattening
misread). Pattern: trust cross-sections + scalar geometry metrics, never the 3D surface render or a
short run.

### Evidence — long-settle (10k) 2-cell rep sweep, adh=5e7, gap=2.05, γ=1e-4
| rep | NN/R | per-cell Psi | **oblate (axial/lateral)** | cross-midplane | deep/R | inside-other-mesh |
|---|---|---|---|---|---|---|
| 4e7 | 1.65 | 0.966 | **0.943** | 26/324 | 0.089 | 30 |
| 2e8 | 1.58 | 0.956 | **0.905** | 32/324 | 0.047 | 34 |
| 1e9 | 1.60 | 0.957 | **0.916** | 23/324 | 0.075 | 32 |
| 5e9 | 1.63 | 0.954 | **0.922** | 21/324 | 0.065 | 20 |

- **oblate < 1 at every rep** ⇒ cells go oblate toward the contact = real flattening (the flat
  junction is visible in `figs/twocell_flattening_VERDICT.png`; the smooth 3D render in 9ec3ee1's
  `twocell_contact_flattening.png` cannot show it).
- **crossing does NOT grow with softer rep** (4e7 → 26, 2e8 → 32) and **deep crossing is
  0.05–0.09 R ≈ 0.4–0.7 µm = ~20% of mean_edge** ⇒ this is **sub-edge interface interleaving** of two
  non-conforming flat meshes touching, NOT rigid-sphere interpenetration. `pen_frac` (node-into-
  nearest-face) reports ~0.05 here, which is actually *honest* for this regime — it only
  underestimates when a node tunnels DEEP (see latent risk below), which does not happen at 2-cell.
- per-cell Psi 0.95–0.97 is high because a 2-cell has only ONE contact face; Psi is N-dependent
  (an interior cell in the OLD n400 has ~12 faces → Psi 0.876). The right 2-cell metric is
  per-junction oblateness, which shows flattening.

### Latent risk to check at N≥400 (NOT a 2-cell problem)
`contact_grid_kernel`'s repulsion fires only for `sign<0 AND min_d < c_rep` (c_rep = 0.30·mean_edge
≈ 0.75 µm). A node pushed deeper than c_rep gets ZERO restoring force → free tunnelling; and the
grid query radius (c_adh+face_reach) caps detection range. At 2-cell the worst depth was 0.67 µm
(just under c_rep), so it never bites. Under N≥400 collective compression a node could exceed c_rep.
**Principled fix if it manifests** (legit lever #3 = fix a real bug, not outcome-tuning): excluded
volume on the INSIDE (sign<0) must be unconditional in depth — remove the `min_d<c_rep` cap on the
repel branch (keep it only as the OUTER soft-shell range on the sign≥0 adhesion side). Deferred until
the N=400 run shows whether deep_pen actually grows; do not pre-emptively change the contact law.

## Audit cleanup (PI 2026-06-23 confirmed issues)

### dead `dcm_coupling_host.py` — REMOVED
`DcmNodeCouplingHost` (node-NODE coupling) was an untracked (`??`) experiment, superseded by A1's
node-FACE re-enable (the `--coupling` flag uses `contact_grid_kernel`'s adh branch, a different
mechanism). Zero imports anywhere. Removed local + gbook.

### A1 `--coupling` double-counts adhesion + makes de-cohesion geometric — SURFACE TO PI
With `--cadherin --coupling` both on, there are **two** cell-cell adhesion channels:
(1) the node-FACE bilinear adhesion (`coh_adh=adh_strength`, geometric/distance-based, releases
when separation > c_adh) and (2) the Rakshit cadherin catch-bond ensemble (force-dependent). This
**double-counts** the adhesion energy AND the node-face channel's geometric release is a
**de-cohesion-by-geometry switch** — which violates the hard rule that de-cohesion must emerge from
catch-bond rupture (memory `feedback-junction-switch-fine-grained`).
**Proper unification (PI decision needed):** the node-FACE adhesion magnitude ω should be MODULATED
by the local cadherin catch-bond state (continuous node-face *geometry* for mesh-independence ×
catch-slip *kinetics* for the force-dependent strength/rupture), so there is ONE adhesion field whose
de-cohesion is emergent. This is a mechanism/contract change → not implemented unilaterally.
Interim: do not run `--coupling` together with `--cadherin` for production de-cohesion claims.

### E triplet-gate is TAUTOLOGICAL — SURFACE TO PI (gate change)
`measure_triplet_angle` returns `gaps.mean()` of the three junction opening angles, which **sum to
360° by construction**, so the mean is **identically 120° for ANY three cell centroids** — the gate
passes regardless of the actual configuration and measures nothing. A non-tautological test must use
either (a) `std(gaps)` ≈ 0 (symmetry of the centroid triangle) or, faithfully, (b) the actual
dihedral angle between the two cell-cell contact interfaces at the junction edge vs Young–Dupré
`cos(φ/2)=η/2`. Per the no-gate-loosening rule this is PI-authored; flagged, not silently rewritten.

### C5 ECM mechano-feedback — confirmed UNWIRED; recommend DEFER
`dcm_substrate_warp.py` has `chan_odde_traction_factor`, `e_sub_to_k_sub`,
`run_dcm_substrate_well_mechano_warp`, `dcm_substrate_well_strainstiff_kernel` (self-test PASS), but
the driver imports only the rigid-dish well/wetting kernels — the mechano-feedback is never launched.
Wiring it touches the substrate-WELL/wetting PROXY path that the proxy-free stack is deprecating, and
the backlog itself notes C5 "needs the deformable/3D substrate (C6 Mikado ECM Warp port) to be
meaningful." Recommend wiring C5 together with C6, not onto the rigid-dish proxy now. Surface to PI.

## Figures
- `figs/twocell_flattening_VERDICT.png` — 4-rep cross-sections, flat junction at every rep (the
  correction to 9ec3ee1).
- `figs/twocell_long_rep{4e+07,2e+08,1e+09,5e+09}.png` — per-rep cross-section + all-nodes diagnostic.
- `figs/n400_contact_valid.png` — (pending) N=400 top-down + equatorial slab vs OLD n400 baseline.

## NEXT
- N=400 aggregation settle (running, gbook A5000, ~50 min ETA) → confirm round all-touching contact
  spheroid at scale + whether deep_pen grows (latent-tunnel check). Compare to OLD n400
  (aniso 1.02, asph 0.013, contact 1.0, V/V0 1.12, Psi 0.876).
- PI decisions: A1 unification, E gate rewrite, C5 defer-vs-wire.
