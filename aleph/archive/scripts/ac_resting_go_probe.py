"""Resting-baseline GO probe — FIRE from-resting with the geometry fix + subdiv-6 membrane.

After (1) the inward-dispersion geometry fix (radial ERM tethers) and (2) recognizing the subdiv-3 membrane
puts 44 pN/node > f_rupt (11.4 pN) so no single-ERM equilibrium exists, this runs the actual GO gate at the
ratified dynamic membrane resolution (subdiv 6 -> 0.69 pN/node < f_rupt) where an equilibrium DOES exist.

Reports the StabilityReport GO fields: residual_start/candidate/end, inner_converged, outer_accepted.
CUDA-gated -> gbook only.
"""
from __future__ import annotations

import dataclasses

from aleph.components.incumbent.assemble import CellConfig
from aleph.components.incumbent.driver import run_from_resting


def go(n_filaments: int, subdiv: int, *, preload: bool, prestrain: float, n_inner: int = 200) -> None:
    cfg = CellConfig(
        n_filaments=n_filaments,
        overlap_free_cortex=True,
        erm_radial_pairing=True,
        membrane_subdivisions=subdiv,
    )
    rep = run_from_resting(
        cfg, n_inner=n_inner, inner_solver="fire",
        preload_erm_balance=preload, cortex_prestrain=prestrain,
    )
    d = dataclasses.asdict(rep) if dataclasses.is_dataclass(rep) else vars(rep)
    keys = ["inner_iters", "inner_attempts", "inner_converged", "inner_tolerance_um",
            "inner_max_displacement_um", "residual_start", "residual_candidate", "residual_end",
            "outer_accepted", "outer_rolled_back", "force_finite", "pos_finite"]
    print("=" * 96)
    print(f"n={n_filaments} subdiv={subdiv} preload={preload} prestrain={prestrain} n_inner={n_inner}")
    for k in keys:
        if k in d:
            print(f"    {k:28s}: {d[k]}")


if __name__ == "__main__":
    import sys
    if "--native" in sys.argv:
        go(70686, 6, preload=False, prestrain=0.0)
        go(70686, 6, preload=True, prestrain=0.0)
    else:
        go(8000, 3, preload=False, prestrain=0.0)   # control: no equilibrium (44pN/node > f_rupt)
        go(8000, 6, preload=False, prestrain=0.0)   # equilibrium exists; pure FIRE
        go(8000, 6, preload=True, prestrain=0.0)    # + ERM preload


def sweep_prestrain(n_filaments=8000, subdiv=6):
    """Coupled preload sweep: ERM preload + cortex pretension. Report residual_start only (n_inner=1)."""
    for ps in [0.0, 1e-5, 5e-5, 1e-4, 5e-4, 1e-3]:
        go(n_filaments, subdiv, preload=True, prestrain=ps, n_inner=1)
