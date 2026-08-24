"""Resting-baseline solver sweep INCLUDING the wired ERM-Schwarz candidate (2026-07-22d fork-2).

Runs run_from_resting at the well-posed operating point (subdiv 6 + geometry fix + ERM preload) with the
tournament solvers, comparing the ERM-Schwarz-augmented tournament against the Newton and WCA-contact
tournaments. Reports the driver's projected residual_start/candidate/end + convergence. CUDA-gated → gbook.
"""
from __future__ import annotations

from aleph.components.incumbent.assemble import CellConfig
from aleph.components.incumbent.driver import run_from_resting


def run(solver, nf, subdiv=6, n_inner=200, cg_it=96):
    cfg = CellConfig(n_filaments=nf, overlap_free_cortex=True, erm_radial_pairing=True,
                     membrane_subdivisions=subdiv)
    r = run_from_resting(cfg, n_inner=n_inner, inner_solver=solver, preload_erm_balance=True,
                         implicit_cg_max_iterations=cg_it, implicit_coarse_iterations=0,
                         implicit_coarse_modes=0)
    print(f"  {solver:20s} n_inner={n_inner} | start={r.residual_start:8.4f} "
          f"cand={r.residual_candidate:8.4f} end={r.residual_end:8.4f} "
          f"conv={r.inner_converged!s:5} accept={r.outer_accepted!s:5}", flush=True)


if __name__ == "__main__":
    import sys
    nf = 70686 if "--native" in sys.argv else 8000
    ni = 200
    for a in sys.argv:
        if a.startswith("--ninner="):
            ni = int(a.split("=")[1])
    print(f"=== ERM solver sweep n={nf} subdiv=6 overlap_free+radial+ERM-preload  n_inner={ni} ===", flush=True)
    run("erm_tournament", nf, n_inner=ni)
    run("erm_schwarz", nf, n_inner=ni)
    run("analytic_implicit", nf, n_inner=ni)
