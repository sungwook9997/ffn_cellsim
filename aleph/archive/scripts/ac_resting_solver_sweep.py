"""Resting-baseline solver sweep at the now-well-posed operating point (subdiv 6 + geometry fix + preload).

FIRE and rest-length preload both fail on the coupled stiff operator. Before building new multigrid, test
whether the EXISTING implicit/conditioned solvers (analytic_implicit with coarse modes, tournament, anderson,
contact_schwarz) can descend the ~2.22 pN cortex-reaction residual now that subdiv>=6 makes an equilibrium
exist and the geometry fix makes the ERM load path radial. Reports residual_candidate + convergence flags.
CUDA-gated -> gbook only.
"""
from __future__ import annotations

from aleph.components.incumbent.assemble import CellConfig
from aleph.components.incumbent.driver import run_from_resting


def run(solver: str, subdiv=6, n_filaments=8000, *, preload=True, cg_it=64, coarse_modes=0,
        coarse_it=8, n_inner=400):
    cfg = CellConfig(n_filaments=n_filaments, overlap_free_cortex=True,
                     erm_radial_pairing=True, membrane_subdivisions=subdiv)
    r = run_from_resting(
        cfg, n_inner=n_inner, inner_solver=solver, preload_erm_balance=preload,
        implicit_cg_max_iterations=cg_it, implicit_coarse_iterations=coarse_it,
        implicit_coarse_modes=coarse_modes,
    )
    print(f"  {solver:20s} cgit={cg_it:3d} coarse={coarse_modes} | "
          f"start={r.residual_start:8.4f} cand={r.residual_candidate:8.4f} end={r.residual_end:8.4f} "
          f"conv={r.inner_converged!s:5} accept={r.outer_accepted!s:5}")


if __name__ == "__main__":
    import sys
    nf = 70686 if "--native" in sys.argv else 8000
    print(f"=== solver sweep n={nf} subdiv=6 overlap_free+radial+ERM-preload ===")
    run("analytic_implicit", n_filaments=nf, coarse_modes=0)
    run("analytic_implicit", n_filaments=nf, coarse_modes=3)
    run("analytic_implicit", n_filaments=nf, coarse_modes=3, cg_it=128)
    run("anderson", n_filaments=nf)
    run("tournament", n_filaments=nf)
    run("contact_tournament", n_filaments=nf)
