"""Deformable-mesh nucleus (pillar P4) — the Active Cell foundation compartment upgrade (I2).

The DCM/FF nucleus is a radial-spring bead CLOUD (``nucleus_shell_kernel``): no lateral area
coupling, no bending, one deformation mode, no rupture. That is the "개살구" nucleus. This track
REBUILDS it as a genuinely deformable triangulated shell using the **membrane Helfrich machinery
already built** (``ff/membrane_surface.py``: dihedral bending + area tension), applied to the
nuclear envelope:

    triangulated lamina shell
      + Seung–Nelson dihedral bending  (κ̃ = 8π·κ_NE / Σ_ref, icosphere-calibrated — reused from membrane)
      + areal elasticity with the framework-#6 lamin split
            small strain → CHROMATIN (soft internal polymer net)
            baseline     → LAMIN-B    (uniform meshwork shell)
            large strain → LAMIN-A/C  (strain-stiffening; engages past the knee)
      + nucleoplasm VOLUME constraint (ν→½ incompressible;  E_vol = ½ K_vol (V−V0)²/V0)
      + nucleoplasm viscous drag
      + LINC tether stubs (nucleus↔cytoskeleton; nonlinear nesprin spring)
      + chromatin internal polymer net
      + EMERGENT envelope rupture past a sourced envelope-tension threshold (NOT tuned)

I2 (this module set) = the standalone deformable-mesh nucleus + its analytic gates. The I7 load path
(perinuclear actin cap → LINC → oblate flatten) is this track's NEXT increment (gated behind I3+I4);
only its standalone oracles (capstan ∮T·κ ds, oblate volume conservation, LINC force-extension) are
pre-authored here.

Analytic-first (§1.5): every gate is a pure-NumPy closed-form oracle + pytest FIRST (``*_analytic.py``,
ZERO Warp import — Warp-only-contract exempt), then the Warp kernel SOURCE (``envelope.py``). The dev
Mac authors kernel source + runs the host oracles only; the lead runs the native CUDA gates on gbook.

CRITICAL INTERFACE (§1.4): the fluid track OWNS ``ac/fluid/domain.py :: Domain.set_nucleus_boundary``.
This track ships ``DeformableNucleusMaskProvider`` (``mask_provider.py``) conforming to that signature —
the moving oblate mesh is I1a's inner relative-no-flux boundary. We code against the signature + a local
test-double; we do NOT edit ``ac/fluid/``. ff/ changes are staged in ``INTEGRATION.md``, never edited.

Sanity Gate (before first native execution, per the Sanity Gate Protocol):
  * dimensional analysis: κ̃ [pN·µm], areal moduli [pN/µm], K_vol [pN/µm²], K_linc [pN/µm];
  * boundary cases: flat hinge → zero bending force; zero areal strain → zero tension;
    sub-knee strain → chromatin tangent, super-knee → lamin-A/C tangent;
  * conservation: closed-mesh volume via the divergence-theorem signed-tet sum; Σf=0 per force;
  * numerical: FD-gradient is the sign arbiter of every force (force == −∂E/∂x to <1e-6·‖f‖);
    add κ̃, K_areal, K_vol, K_linc, K_chrom to kmax (CFL);
  * measurement consistency: sphere bending → 8πκ (Willmore); Young-Laplace ΔP = 2σ/R.

Parameter gate: ``params_i0b2.yaml`` (I0-B2). Magnitude gates (E_nuc-dependent tension, rupture strain,
oblate aspect) are INVALID until their blocking parameter closes; the STRUCTURAL gates (8πκ topology,
FD-gradient sign, volume conservation to machine precision, knee modulus-ratio crossover, Young-Laplace,
rupture-emerges-above-threshold) are magnitude-independent and pass now.

Runtime: NVIDIA Warp on CUDA GPU only (I0-A). No HOOMD, no CPU simulation path.
"""
