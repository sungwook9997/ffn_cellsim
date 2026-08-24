"""Fluid-first poroelastic substrate (pillar P1) — the Active Cell foundation layer.

The cell is a fluid-saturated porous continuum coupled to explicit filaments; the pore fluid is a
first-class always-solved DOF, not an optional channel. Physics level is inertialess Darcy/Biot
(intracellular Re ~ 1e-5..1e-10 — momentum inertia is negligible), NOT Navier-Stokes.

Increment I1a (this module set) = the CONSERVATIVE p/mass Biot substrate:

    S p_dot + alpha div(v_s) + div(q) = s_water          fluid-content balance   (S = 1/M)
    p = p_ext + p_bar + p_excess                          pressure split (p_bar = Pi_0 pre-stress)

solved on the LIVE membrane-minus-nucleus domain with conservative moving-boundary remap, an outer
membrane hydraulic-flux BC, and a live-nucleus relative no-flux BC. It runs under ONE outer physical
clock with a frozen-state inner mechanical solve (retiring the 0-D 6*pi*eta*R -> gamma_solid clock
lock in the same change). The mixed q-p / v_f velocity form is I1b; conservative solute transport is
I1c; neither lives here.

Reuse (NOT the held-face diffusion stencil): device arrays + Peskin spread/interp from
``ff/biot_fluid_warp.py``; live membrane geometry from ``ff/membrane_surface.py``.

Package map (I1a + I1b + I1c authored; Warp-CUDA source + host-numpy acceptance oracles):

    Warp-CUDA runtime SOURCE (gated on the A5000 by the lead; constructing requires CUDA):
      field_grid.py     conservative FV device arrays + a separate filament HashGrid
      biot_substrate.py I1a conservative p/mass update kernel + PressureCoupling.accumulate (§1.4)
      boundary.py       membrane hydraulic-flux BC + nucleus relative-no-flux BC
      domain.py         OWNS Domain + set_nucleus_boundary(provider) (§1.4); static-sphere placeholder
      scheduler.py      outer physical clock + inner mechanical solve (6piR->gamma_solid retirement)
      velocity.py       I1b Darcy discharge q + pore-fluid velocity v_f = v_s + q/phi
      transport.py      I1c conservative RAD monomer transport + MonomerField

    Host-numpy ACCEPTANCE ORACLES (dev-Mac green; the CPU verification of the identical stencil):
      consolidation_analytic / greens_analytic / manufactured   I1a closed forms (Terzaghi/Green/MMS)
      fv_reference.py   I1a discrete conservative-FV reference (constant-state/mass/content/SPD/convergence)
      ibm_reference.py  Peskin spread/interp adjoint (Interp = h^d Spread^T) + net-force projection
      darcy_analytic.py I1b Darcy kinematics (flux linearity, v_f-v_s=q/phi, radial conservation)
      transport_analytic.py / transport_reference.py  I1c FRAP + advection + machine-precision conservation

    Parameter ledgers: params_i0b1.yaml (I1a, closed) · params_i0b1b.yaml (phi GAP) · params_i0b1c.yaml
      (D_c draft, c0 GAP).   Hand-off: INTEGRATION.md (ff/ retirement patch-notes) · NATIVE_GATE_SPEC.md.

Parameter gate: ``params_i0b1.yaml`` (I0-B1, PI-ratified 2026-07-16). alpha=1.0 is ratified
(derived-by-limit); analytic gates still sweep alpha as an oracle. I1b/I1c native runs are blocked on the
phi / c0 GAPs (surfaced to PI); the analytic gates are parameter-agnostic and pass now.

Sanity Gate (before first native execution, per the Sanity Gate Protocol):
  * dimensional analysis of every term in the fluid-content balance;
  * boundary cases: impermeable limit conserves mass; constant-state preservation;
  * conservation: global fluid-content change == integrated membrane flux to machine precision;
  * numerical: CFL for the explicit consolidation step; pressure-work sign; net-force projection
    (closed-system COM drift ~ 0);
  * measurement consistency: Terzaghi/Green analytic oracle; manufactured moving-boundary solution.

Runtime: NVIDIA Warp on CUDA GPU only (I0-A). No HOOMD, no CPU simulation path.
"""
