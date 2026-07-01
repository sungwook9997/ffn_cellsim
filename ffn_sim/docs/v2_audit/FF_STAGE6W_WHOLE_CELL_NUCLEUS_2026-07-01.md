# FF Stage 6W — whole-cell AFM: shared compartments + the nucleus does not change the cortical-tension measurement

**Date:** 2026-07-01  **Engine:** FF (Warp, A5000)  **Branch:** dcm/main
**PI directives:** "nucleus mechanics나 다른 것이 구현된 상태로 돌리는 것 맞나?" + "진짜 전부 다 필요할텐데" +
"둘이 공유하면 안되나? dcm 것 그대로 써도 되면 그대로 갖고 오고, 정밀도에 차이가 있으면 새로이 구현" — build the
whole-cell compartments, reuse-first, share ONE implementation across the two Warp engines.
**Follows** FF_STAGE6V (the γ-floor resolved: cortical tension is osmotic-pressure-borne).

## One line

The whole cell (cortex + turgor + **stiff nucleus**) is now assembled from **shared engine-agnostic
compartments**, and a native virtual-AFM sweep shows the nucleus is **mechanically absent from the
cortical-tension measurement** (γ is nucleus-independent to 4 sig figs below the geometric contact strain),
engaging only under deep confinement — and even there it is negligible in the current over-stiff-cortex
regime. This closes the "are we running with nucleus mechanics?" question for the cortical-tension result:
cortex + turgor was the correct and sufficient model for it.

## Part 1 — shared compartments (`common/compartments.py`)

The full HOOMD compartment stack (`archive/hoomd_legacy/cell/`: nucleus, membrane_surface, cytoplasm, IF,
MT, SF, LINC) was **HOOMD-only** — not present in the going-forward Warp engines, so the FF AFM/γ work had
been running cortex + turgor only. The fix reuses the already-validated DCM kernel rather than rewriting:
`common/compartments.py` is a **verbatim, bit-identical port** of `dcm/radial_shell_warp` — one radial-shell
kernel with three laws, the engine-agnostic home both Warp engines import:

- **law 0 — nucleus** : bilinear about R0 (chromatin slope inside a lamin knee, +lamin slope past it).
- **law 1 — membrane**: S=4πR², γ=γ_mem+K_A·(S−A0)/A0, ΔP=2γ/R (inward Laplace).
- **law 2 — turgor**  : V=4/3πR³, ΔP=Π₀−K_vol·(V−V0)/V0 (outward osmotic).

Bit-parity guard (`tests/ff/test_compartments_shared.py`): `common.run_radial_shell_warp` == the DCM original
to `max|Δ| = 0`. `resolve_nucleus` bridges KU-3.B2 lit constants to FF units (1 Pa ≡ 1 pN/µm²): E_nuc=5 kPa
(band 1–10), ratio_lamin=3× (in-situ 1.4–5, NOT the 10× isolated value), knee_strain=0.10 →
k_chrom=4π·E·R/n_beads, k_lamin=(ratio−1)·k_chrom, d_knee=knee_strain·R, F_knee=k_chrom·d_knee (grid-invariant:
n_beads·k_chrom = 4π·E·R). **Cytoplasm (H.10) was already shared** — both engines run η=65.9 Pa·s (Dessard)
per-node drag. Follow-up (PI, /loop-idle): re-point `dcm/radial_shell_warp` at `common/` to drop the temporary
duplication.

## Part 2 — whole-cell AFM wiring (`simulate_whole_cell_compression_on_device`)

A NEW function (the validated `simulate_compressed_shell_on_device` cortex-only path is untouched → no
regression). It appends a Fibonacci-sphere nucleus bead cloud to the cortex node array and adds the shared
law-0 restoring force (`nucleus_shell_kernel`, verbatim law-0 math, accumulating). Two **physical**
cortex↔nucleus couplings:

1. **incompressible-cytoplasm displacement** — the nucleus occupies volume, so the compressible cytoplasm is
   V_cyto = V_hull − V_nuc; compressing the cell squeezes a smaller cytoplasm → ΔP rises faster (the nucleus
   is felt hydrostatically before contact). *In regulated-ΔP mode this coupling is off by construction* (ΔP is
   clamped), which is exactly the cortical-tension regime.
2. **direct compression** — once the plate half-gap drops below R_nuc (strain > 1−R_nuc/R_cell ≈ 0.667 for
   R_nuc/R_cell = 2.5/7.5 = 1/3), the plates press the nucleus beads directly and its stiff strain-stiffening
   law contributes.

Cortex kernels (bending/crosslink/myosin/turgor) act on cortex indices only (turgor launched `dim=Nc`); the
plate and integrator act on all nodes; centroid + convex hull use cortex nodes only.

## Part 3 — native sweep (N=70686, regulated ΔP=40 Pa interphase)

`gamma_apparent = ½·ΔP·R_eq` (liquid-drop), sweeping strain across the contact threshold, cortex-only vs
+nucleus:

| strain | F cortex-only [pN] | F +nucleus [pN] | γ cortex-only | γ +nucleus | contact |
|---|---|---|---|---|---|
| 0.050 | 2.07e6 | 2.03e6 | 0.1502 | 0.1502 | – |
| 0.150 | 1.92e7 | 1.88e7 | 0.1502 | 0.1502 | – |
| 0.300 | 7.72e7 | 7.56e7 | 0.1502 | 0.1502 | – |
| 0.450 | 1.74e8 | 1.70e8 | 0.1505 | 0.1505 | – |
| 0.600 | 3.09e8 | 3.02e8 | 0.1543 | 0.1543 | – |
| **0.667** | 3.81e8 | 3.73e8 | 0.1571 | 0.1572 | **contact** |
| 0.720 | 4.44e8 | 4.35e8 | 0.1607 | 0.1609 | contact |
| 0.780 | 5.21e8 | 5.10e8 | 0.1627 | 0.1629 | contact |

- **γ_apparent is nucleus-independent** (identical to 4 sig figs across arms) and **confinement-independent**
  (~0.150 mN/m over 0–60% strain; drifts to 0.163 at 78% as R_eq grows) — the liquid-drop signature, ≈ the
  Fischer-Friedrich interphase 0.17. *Caveat:* this γ is a **turgor-Laplace proxy** (½·ΔP·R_eq with regulated
  ΔP), NOT a force-fit of the indentation curve the way a real AFM extracts tension. The 0.150→0.163 drift is
  the single-radius proxy artifact (R_eq bulges ~8% at fixed ΔP), within Fischer-Friedrich's ±10% — i.e. the
  cortex is only *approximately* liquid-drop. The force-derived γ would disagree (see gap #1).
- `nucleus_contact` flips True at exactly strain **0.667** = 1−R_nuc/R_cell (Claim B).
- F_plate is a **constant 0.98×** with the nucleus at ALL strains (below AND above contact) — a flat ratio,
  the fingerprint of a setup offset, not a nucleus mechanical contribution (investigated next).

## Part 4 — clean isolation (the decisive test)

The flat 0.98 came from the nucleus changing V0 (=V_hull−V_nuc) → k_plate (∝1/V0, 3.8% higher) → dt_mu (via
the CFL) → slightly-less-converged at fixed 3000 steps. Re-run with **identical explicit k_plate + dt_mu
passed to both arms** (removing the coupling):

| strain | F cortex-only [pN] | F +nucleus [pN] | ΔF(nucleus) | contact |
|---|---|---|---|---|
| 0.300 | 7.720114e7 | 7.720114e7 | **−3.6e-7 pN (−0.000%)** | below |
| 0.667 | 3.811236e8 | 3.811236e8 | +2.1e-6 pN (+0.000%) | threshold |
| 0.780 | 5.209296e8 | 5.209770e8 | **+4.7e4 pN (+0.009%)** | past |

**Below contact the cortex sub-system is bit-identical with/without the nucleus** → the nucleus is genuinely
decoupled (Claim C confirmed: the 0.98 was purely numerical). **Past contact the nucleus adds a small
*positive* stiffening** (+0.009% at strain 0.78) — real (correct sign, physics restored) but **k_nuc-limited**
and swamped by the over-stiff cortex.

## Part 5 — interpretation (honest + regime-specific)

- **The cortical-tension measurement does not require the nucleus.** γ is the turgor Laplace tension (ΔP·R/2);
  the plate never reaches the nucleus at the strains where cortical tension is measured (<50%). So the
  cortex+turgor model used to resolve the γ-floor (6V) was correct and sufficient — the "are we running with
  nucleus mechanics?" concern does not change that result.
- **The nucleus IS wired and engages at deep strain** (>67%), as it must for whole-cell high-strain modulus
  and confined migration.
- **Why the nucleus is numerically negligible here is regime-specific.** Its shell stiffness (k_chrom+k_lamin
  ≈ 734 pN/µm) is ~10⁵× below the current cortex plate force (~10⁸ pN). That cortex force is itself ~10³× too
  high (FF_STAGE6V honest gap #1: the inextensible reshape makes the shell resist area change too strongly;
  the real cortex is soft/liquid-drop via turnover, nN-scale force). **HYPOTHESIS (untested — not run here):**
  in a realistic soft cortex (nN-scale force) the same ~734 pN/µm nucleus would be a comparable fraction and
  would plausibly produce the measurable deep-indentation stiffening AFM sees — directionally consistent with
  the known cell-mechanics picture (cortical tension at low deformation, nucleus-dominated response at high
  indentation; Caille 2002, Lammerding, Stephens 2017). This is a plausible extrapolation, NOT a result of this
  milestone — the soft-cortex counterfactual has not been run. What IS established: "nucleus negligible" is a
  statement about the current over-stiff-cortex regime, not a claim that the nucleus is mechanically
  unimportant.
- **Omitted couplings (honest gap).** The two cortex↔nucleus couplings modelled here (volume displacement +
  direct contact) are the *minimal physical set*. A real cell also couples the nucleus to cortex deformation via
  **LINC complexes / perinuclear actin / the microtubule cage**, which transmit sub-contact strain to the
  nucleus earlier than pure geometry predicts. We omit these; adding them (with the soft-cortex regime) is the
  path to a quantitative whole-cell deep-AFM response, and is PI-roadmap, not part of this milestone.

## Sanity gate

- **Measurement-protocol consistency:** γ extracted by the same liquid-drop route as the experiment, on the
  deformed (pressed) cell. ✓
- **No magic number:** nucleus constants are KU-3.B2 lit values via `resolve_nucleus` (band-guarded);
  pressure_setpoint = Fischer-Friedrich interphase ΔP; nothing tuned to a target. ✓
- **Grid-invariance:** nucleus bridge is intensive (n_beads·k_chrom = 4π·E·R); native N=70686. ✓
- **No regression:** new function; validated cortex path untouched; CPU tests + kb-check green. ✓
- **Isolation:** with fixed numerics the cortex is bit-identical below contact and the nucleus ΔF is +sign,
  k_nuc-scale, contact-gated — physically correct. ✓

## Adversarial verification

A 4-lens skeptic panel (numerics/convergence, code-correctness, physics/measurement-protocol,
literature/biology) each read the code + both datasets and tried to REFUTE claims A–E; a synthesizer
aggregated. **All five claims confirmed, none refutable** (conf 0.72–0.97). Numerics independently
reproduced the k_plate ratio to 6 sig figs (1.03846 = V0c/(V0c−V_nuc)) and the below-contact bit-identity at
machine epsilon (rel diff 4.6e-15). The panel's two carried caveats — γ is a turgor-Laplace *proxy* (not a
force-fit), and the soft-cortex corollary is an *untested extrapolation* — are folded into Parts 3 and 5 above.

## Figures
- `outputs/ff/figs/wholecell_afm_native.png` — Panel A: F_plate(strain) cortex-only vs +nucleus (overlap;
  clean-isolation ΔF annotated). Panel B: γ_apparent(strain) confinement- + nucleus-independent (interphase
  band overlaid). Panel C: xz cross-sections at strain 0.15/0.60/0.78 — cortex shell squashed between plates,
  nucleus untouched at 0.15/0.60, **compressed only at 0.78**.

## Files
- `common/compartments.py` — shared radial-shell kernel (bit-identical port) + `resolve_nucleus` (KU-3.B2).
- `ff/network_warp.py` — `simulate_whole_cell_compression_on_device` + `nucleus_shell_kernel` +
  `_seed_nucleus_cloud` (commit 8e34b39).
- `tests/ff/test_compartments_shared.py`, `tests/ff/test_compression.py` (whole-cell engagement + reduction).

Related: [[project-cell-mechanics-extend]], [[project-gamma-floor-likely-deficit]], FF_STAGE6V (pressure-borne
tension), FF_STAGE6U (virtual-AFM). Sources: Fischer-Friedrich 2014 (10.1038/srep06213); KU-3.B2 (E_nuc,
ratio_lamin); Dessard 2024 (η).
