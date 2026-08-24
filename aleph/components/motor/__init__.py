"""Head-resolved Stam-Hocky NMII minifilament (pillar P2) — the Active Cell active-force root.

The ONLY actomyosin contractile-motor primitive in the new engine: one head-resolved, multi-particle
bipolar minifilament from which cortical tension, stress fibers, and the perinuclear cap all reuse the
SAME motor (P2). There is NO lumped ``f_myo`` constant-force dipole and NO imposed ``N_side * F_head``
single-link substitute anywhere in the final engine — those are the RETIRED aggregate kernels
(``ff/network_warp.myosin_kernel`` + ``ff/myosin_linear.minifilament_kernel``), demoted to Warp
diagnostics/controls only (INTEGRATION.md).

Increment I3 (this module set) = the head-resolved minifilament:

    * explicit backbone beads (bipolar bare-zone + two anti-parallel head-bearing arms);
    * explicit heads on both sides, each a crossbridge spring (stiffness ``k_xb``) to its bound actin site;
    * per-head Hill force-velocity stepping   v(F) = v0 * (1 - F/F_s) / (1 + (F/F_s)/kappa);
    * per-head Bell slip detachment           p_off(f) = k_off0 * exp(|f| / f0);
    * all CUDA-resident, no per-step host state (I0-A GPU-only contract).

FV-form reconciliation (⚠ recorded, see params_i0b3.yaml::kappa_hill): I0-A (2026-07-16, HARD) mandates a
Hill force-velocity; a 2026-07-07 PI note ratified a LINEAR aggregate FV for non-muscle NMII. The general
Hill hyperbola above has the linear law as its exact ``kappa -> inf`` limit and Hill-1938 muscle at
``kappa = 0.25``; the curvature ``kappa`` is left an I0-B3 GAP for PI to select by isoform/assay. The Hill
MACHINERY is built (mandated); the curvature is NOT chosen here.

Honest scope (P2, do NOT overclaim): at the quasi-static inner equilibrium the inter-anchor sliding velocity
``v_slide -> 0`` so each engaged head sits at its stall force (a constant) — force-velocity bites only the
transient. The fine-grained content that SURVIVES at the settled state is (1) the DERIVED contractile
magnitude, (2) the Bell-EMERGENT engaged fraction (fewer heads bound under load, self-limiting), and (3)
topology reformation. The native gamma-floor magnitude is density-floored → report-not-tune (§6.2): NEVER
add heads/density/stiffness to close a floor.

Analytic modules (host NumPy, ZERO Warp — acceptance oracles, Warp-only-contract exempt by construction):
  * ``hill_fv_analytic``        — Hill-1938 force-velocity closed form + inverse + linear limit + work/power.
  * ``bell_kinetics_analytic``  — Bell slip + Pereverzev catch-slip off-rates, per-tick probabilities,
    steady-state engaged fraction (slip falls with load; NMII catch-slip rises to F* then falls, Kovacs 2007).
  * ``ensemble_stall_analytic`` — mean-field + stochastic (binomial) head population; ensemble stall EMERGES.
  * ``minifilament_topology``   — bipolar topology/count invariants + per-head Newton force-balance closure.
  * ``powerstroke_analytic``    — the power-stroke -> force coupling: the walked abscissa advances the actin
    attachment point (x_att = anchor + abscissa*walk_dir), so a stepping head makes a DIRECTED contractile
    force and force-velocity self-limits at F_stall. Mirrors the device crossbridge/load kernels bit-for-formula
    and is the LOCAL gate that catches a force<->power-stroke decoupling (which the un-launched kernels hide).

The Warp kernel SOURCE (device runtime, authored not executed on the dev Mac) lives in ``hand.py`` (the
per-head hand/KMC device API — I3 OWNS the §1.4 frozen interface) and ``minifilament_warp.py``
(``MyosinForce.accumulate`` + topology build). Downstream tracks (I4/I6/I8) CONSTRUCT hands via the
``hand.py`` API in their OWN presets.py and never edit the KMC core.

Sanity Gate (before first native execution, per the Sanity Gate Protocol):
  * dimensional analysis: v0 [um/s], F_s [pN], k_xb [pN/um], k_off0 [1/s], f0 [pN], kappa [-];
  * boundary cases: force-free at v=v0 (F=0); zero velocity at F=F_s (stall); OFF (no bound heads) => zero force;
  * conservation / sign: ATP-per-step work >= 0 below stall, efficiency eta = F*d/dG_ATP in [0,1) (2nd law);
  * numerical: per-head Newton closure residual -> 0; k_xb, k_on/k_off rates enter CFL (kmax);
  * sign-sense: minifilament PULLS its two actin anchors together (contractile, never extensile);
  * measurement consistency: single-head FV matches the Hill-1938 closed form; ensemble stall EMERGES from
    the Bell-governed bound population, it is NOT the imposed product N_side * F_head.

Runtime: NVIDIA Warp on CUDA GPU only (I0-A). No HOOMD (the archived Stam-Hocky myosin.py / motor.py are a
read-only geometry/kinetics port SPECIFICATION, never imported or executed), no CPU simulation path.
"""

from __future__ import annotations

from aleph.components.motor.bell_kinetics_analytic import (
    attach_probability,
    bell_off_rate,
    bell_f0_from_x_beta,
    catch_slip_off_rate,
    catch_slip_peak_force,
    detach_probability,
    engaged_fraction_catch_slip,
    engaged_fraction_steady,
)
from aleph.components.motor.ensemble_stall_analytic import (
    ensemble_stall_meanfield,
    ensemble_stall_stochastic,
)
from aleph.components.motor.hill_fv_analytic import (
    hill_force,
    hill_velocity,
    linear_velocity,
    power_per_head,
    step_work,
)
from aleph.components.motor.minifilament_topology import (
    MinifilamentTopology,
    head_newton_residual,
    working_stroke_strain,
)
from aleph.components.motor.powerstroke_analytic import (
    attachment_point,
    crossbridge_force,
    stall_abscissa,
    step_engaged_head,
)

__all__ = [
    "attach_probability",
    "bell_off_rate",
    "bell_f0_from_x_beta",
    "catch_slip_off_rate",
    "catch_slip_peak_force",
    "detach_probability",
    "engaged_fraction_catch_slip",
    "engaged_fraction_steady",
    "ensemble_stall_meanfield",
    "ensemble_stall_stochastic",
    "hill_force",
    "hill_velocity",
    "linear_velocity",
    "power_per_head",
    "step_work",
    "MinifilamentTopology",
    "head_newton_residual",
    "working_stroke_strain",
    "attachment_point",
    "crossbridge_force",
    "stall_abscissa",
    "step_engaged_head",
]
