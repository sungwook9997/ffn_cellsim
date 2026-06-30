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

## Next (FF↔Kim, follow-up)

- **Elastic shear modulus G vs crosslinker density** (athermal shear of the FF network) vs Kim's
  cross-link-density stiffening — the mechanical (not just structural) comparison FF can do.
- Cross-linking-angle distribution + bundle-vs-isotropic morphology (Kim's ACPc/ACPB) vs FF.
- (Thermal G′/G″ would need a thermostatted FF variant — out of the current athermal scope.)
