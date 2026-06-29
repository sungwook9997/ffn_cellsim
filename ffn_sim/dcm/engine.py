"""Single import surface for the Phase-C Warp DCM engine.

Re-exports the validated entry points + per-force kernels so production code imports
ONE module. See ``ENGINE.md`` (this directory) for which entry point to use, the force
set, the hybrid/remesh architecture, and the parity contract + GPU verdict.

Entry points (each validated; pick per ENGINE.md):
  run_multicell  — canonical hybrid loop: N cells; per-cell turgor + cortex edges +
                   node-node cohesion + node-face contact (+ hash-grid, host remesh).
  run_hybrid     — single-cell hybrid (turgor + edges + remesh).
  run_crawl      — single cell + lamellipodial traction-tether (substrate crawl).
  run_diff       — end-to-end DIFFERENTIABLE loop (grads w.r.t. dP0, k_edge).

Per-force runners (parity-gated vs committed HOOMD refs; building blocks for custom loops):
  integrator : run_baoab_warp, run_shake_warp, run_fixman_warp
  cortex     : run_radial_shell_warp, run_harmonic_bond_warp, run_harmonic_angle_warp,
               run_wca_pair_warp
  DCM        : run_node_face_contact_warp, run_dcm_turgor_warp, run_dcm_cohesion_warp,
               run_dcm_substrate_well_warp, run_dcm_substrate_wetting_warp,
               run_lamellipodium_tether_warp

Importing this module pulls in ``warp`` (each submodule calls ``wp.init()``; idempotent).
"""

from __future__ import annotations

# canonical loops
from ffn_sim.dcm.dcm_warp_hybrid_multicell import run_multicell
from ffn_sim.dcm.dcm_warp_hybrid import run_hybrid
from ffn_sim.dcm.dcm_warp_crawl import run_crawl
from ffn_sim.dcm.dcm_warp_diff import run_diff

# integrator stack
from ffn_sim.dcm.baoab_warp import run_baoab_warp
from ffn_sim.dcm.shake_warp import run_shake_warp
from ffn_sim.dcm.fixman_warp import run_fixman_warp

# compartment / cortex forces
from ffn_sim.dcm.radial_shell_warp import run_radial_shell_warp
from ffn_sim.dcm.network_warp import (
    run_harmonic_bond_warp, run_harmonic_angle_warp, run_wca_pair_warp)

# DCM forces
from ffn_sim.dcm.dcm_contact_warp import run_node_face_contact_warp
from ffn_sim.dcm.dcm_turgor_warp import run_dcm_turgor_warp
from ffn_sim.dcm.dcm_cohesion_warp import run_dcm_cohesion_warp
from ffn_sim.dcm.dcm_substrate_warp import (
    run_dcm_substrate_well_warp, run_dcm_substrate_wetting_warp)
from ffn_sim.dcm.lamellipodium_warp import run_lamellipodium_tether_warp

#: the recommended top-level entry points, by use-case (see ENGINE.md)
ENTRY_POINTS = {
    "multicell": run_multicell,
    "single_cell": run_hybrid,
    "crawl": run_crawl,
    "differentiable": run_diff,
}

__all__ = [
    "run_multicell", "run_hybrid", "run_crawl", "run_diff", "ENTRY_POINTS",
    "run_baoab_warp", "run_shake_warp", "run_fixman_warp",
    "run_radial_shell_warp", "run_harmonic_bond_warp", "run_harmonic_angle_warp",
    "run_wca_pair_warp",
    "run_node_face_contact_warp", "run_dcm_turgor_warp", "run_dcm_cohesion_warp",
    "run_dcm_substrate_well_warp", "run_dcm_substrate_wetting_warp",
    "run_lamellipodium_tether_warp",
]
