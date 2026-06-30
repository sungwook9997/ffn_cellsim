# FF single-cell network vs Taeyoon Kim's BD actin-network model (Stage 6i)

**Date:** 2026-06-30 · **Engine:** `ffn_sim/ff/` · PI correction: the single-cell FF network must be
validated against a PEER fine-grained model (Taeyoon Kim), NOT fed into the multi-cell spheroid
A/A0 = a+b/R+c/R² law (that was a category error — A/A0 is the tissue-scale law; reverted).

## Peer oracle

Taeyoon Kim, *Simulation of Actin Cytoskeleton Structure and Rheology* (MIT MS thesis, 2007;
`references/181655768-MIT.pdf`): a 3D Brownian-dynamics cubic box of polymerizing actin filaments
cross-linked by ACPs, measuring network morphology (mesh size, connectivity, cross-linking angle) +
viscoelastic rheology. (Also Kim's Soft Matter 2021 `d0sm01911a.pdf` + Acta Biomater 2025
`1-s2.0-S1742706125004039` are in `references/`.)

FF is athermal mechanical-equilibrium (no Brownian), so the directly-comparable Kim observables are
STRUCTURAL (thermal G′/G″ from MSD are not FF-native; elastic modulus is a follow-up).

## Result — FF reproduces Kim's structural scaling

`ff/kim_network.py` builds a 3D cubic-box cross-linked actin network at concentration C_A (lit-anchored:
actin rise 2.7 nm/monomer, ρ_L = C_A·N_A·rise), ACP crosslinks at ratio R = C_ACP/C_A.

| observable | FF result | Kim / actin literature |
|---|---|---|
| mesh ξ(151 µM) | 64 nm (ρ_L^(−1/2)); measured 78 nm | ~50–100 nm cortex mesh (Kim; Schmidt 1989) |
| mesh-size scaling | **ξ ∝ C_A^(−0.48)** (measured, box network) | **ξ ∝ C_A^(−1/2)** (Schmidt/MacKintosh/Kim) |
| connectivity | **z = 2R** (exact, linear in ACP ratio) | connectivity rises with crosslinker ratio (Kim Fig 2-11) |

Figure: `outputs/ff/figs/kim_network_comparison.png`. Tests: `tests/ff/test_kim_network.py` (3).

The FF actin-network construction reproduces the semiflexible-network mesh-size scaling (ξ ∝ C_A^(−1/2))
and the linear connectivity-vs-crosslinker-ratio that Kim's BD model + the actin-network literature
establish — validating the FF single-cell network's STRUCTURAL physics against the peer model. No
tuning: C_A and R are controlled variables, the rise/ρ_L are lit-anchored.

## Mechanics — FF reproduces the elastic floppy→rigid transition (Stage 6i, `shear_modulus`)

Athermal simple-shear of the FF cross-linked box network (pin top/bottom at the affine-sheared
position, relax interior under actin-segment + crosslink springs + bending, read σ_xz/γ):

| connectivity z | G [pN/µm²] | regime |
|---|---|---|
| 0.2–0.6 | ~0 | FLOPPY (sub-threshold) |
| 1.0 | 0.18 | onset |
| 2.0 | 0.58 | rising |
| 3.0 | 4.2 | RIGID |
| 4.0 | 6.8 | rigid |

G≈0 below a connectivity threshold then rises steeply — the cross-linked-network elastic
**floppy→rigid transition** (Head, Levine & MacKintosh 2003 PRE 68; the regime structure Kim 2007
builds on). FF reproduces it. (`k_xl` sets the G magnitude scale; the TRANSITION/trend is the robust
Kim comparison — magnitude needs the matched ACP/actin stiffness.) Figure panel (C).

This validates the FF single-cell network's MECHANICS (not just structure) against the peer
cross-linked-network picture. Together: FF reproduces Kim/MacKintosh structure (ξ∝C_A^−1/2, z=2R) AND
mechanics (rigidity transition) — the correct single-cell validation.

## Mechanics MAGNITUDE — absolute G, lit-anchored stiffness (Stage 6c, task c — 2026-06-30)

The Stage-6i `shear_modulus` used a SINGLE knob `k_xl` for BOTH the actin backbone and the crosslink
junction (conflating them), so the absolute G was uncalibrated. Stage 6c separates them, each
SEPARATELY SOURCED (no tuning to outcome), and resolves the magnitude. Figure
`outputs/ff/figs/kim_shear_modulus.png`.

**Two stiffnesses, both lit-anchored:**
- **Actin backbone** — INEXTENSIBLE (FF-native NF2007 §5.3 reshape) by default, justified by the
  sourced axial modulus EA = 4.4e4 pN (Kojima, Ishijima & Yanagida 1994 PNAS 91:12962; ⇒ k_axial =
  EA/L_seg ≈ 1.8e5 pN/µm). A finite-EA spring mode is available for the sourced-stiffness regime.
- **Crosslink junction** `k_xl` — α-actinin **4.6e5 pN/µm** (455 pN/nm), filamin **8.2e5 pN/µm**
  (Ferrer 2008 PNAS 105:9221, AFM).

**Result (1 pN/µm² = 1 Pa exactly):**
1. **Crosslink-limited regime (k_xl ≪ EA/L_seg):** G is **LINEAR in k_xl** over 5 decades (G/k_xl
   constant) and matches the analytic affine form **G ≈ f_na·k_xl·ρ_L·ℓc** with a STABLE non-affine
   factor **f_na = 0.024** (< 1; the affine estimate is the upper bound, the network relaxes
   non-affinely — Head/Levine/MacKintosh). This is the robust, no-tuning analytic cross-check
   (analytic is primary per the oracle-is-crosscheck rule). Tests:
   `test_shear_modulus_crosslink_limited_linear`, `test_shear_modulus_rises_with_connectivity`.
2. **Sourced-stiffness regime (k_xl ~ EA/L_seg, finite EA):** at C_A=150–300 µM, z=3, the FF athermal
   **enthalpic** G is **~0.1–0.9 MPa** — cortex-stiff (the precise value is density- and
   non-affinity-dependent; reshape gives the rigid-actin upper bound, finite-EA ~½).
3. **The crosslinker-stiffness UNIT-SLIP (the headline):** the repo's `ALPHA_ACTININ.link_k =
   FILAMIN.link_k = 0.1 pN/µm` is ~**4.5e6× too soft** vs Ferrer (a 1000× pN/µm-vs-pN/nm slip
   COMPOUNDED with adopting the AFINES soft surrogate ~0.1 as physical). At the broken value G ≈
   0.19 Pa (≈10⁶× softer than the cortex); the sourced value lands in the cortex-stiff range. So the
   crosslinker stiffness spanned the *entire* physiological range — the elastic cortex was ~10⁶× too
   soft. (NOTE: this is the ELASTIC modulus axis; it does NOT change the active-γ floor, which is
   force-magnitude limited — [[project-gamma-floor-likely-deficit]].)

**Comparison branches.** FF athermal = the ENTHALPIC branch (in-vivo cortex ~0.1–1 kPa; Kim high-f
enthalpic 100–300 Pa). The dilute in-vitro Kim/Gardel low-f G' (0.07–1000 Pa) is the THERMAL/entropic
branch FF does not target. FF's sourced enthalpic G (~kPa–MPa, dense) is in the right enthalpic regime;
an exact matched absolute needs a finite-EA implicit solve at matched (C_A, R, prestress) — the follow-up.

**PI-gated (do NOT edit live constants unilaterally):** correcting `hand_kmc` `link_k` 0.1 → 4.6e5
(α-actinin) / 8.2e5 (filamin) changes a production runtime constant (feeds the γ-floor cortex + the
CFL) and the Ferrer stiffness datum is not yet in the KB → surface to PI + register SE rows first
(`references/SE_REGISTRATION_CANDIDATES_2026-06-30.md` §5). The `configs/phase1_h3.yaml` actin
intra-filament `k_xl=1e-7 N/m` is a separate clean 1000× slip → 1e-4 N/m, and its KU-3.19/Furuike
attribution is wrong (comment fix).

## Next (FF↔Kim, follow-up)

- Finite-EA IMPLICIT solver for the exact sourced-stiffness absolute G at matched (C_A, R, prestress)
  vs Kim's enthalpic high-f branch (the explicit solver converges at sourced stiffness only when
  k_xl ~ k_axial; the crosslink-limited absolute uses reshape).
- Cross-linking-angle distribution + bundle-vs-isotropic morphology (Kim's ACPc/ACPB) vs FF.
- (Thermal G′/G″ would need a thermostatted FF variant — out of the current athermal scope.)
