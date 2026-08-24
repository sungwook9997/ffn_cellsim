"""GATE A blocker (B), consultation finding #2 (rotational invariance): does the directional crossbridge create a
net TORQUE, so no free-cell static equilibrium exists?

The directional crossbridge force ``f = k (d·w - r0) w`` uses a FROZEN world-space ``walk_dir`` w; such a noncentral
force is not automatically rotationally invariant. Each bound head applies +f at the head and -f at the actin
attachment, a couple ``(r_head - r_attach) × f`` that is nonzero whenever the head sits off the actin (defect#3). The
bipolar minifilament SHOULD cancel these couples (+ and - heads oppose), but that must be measured: if the summed
couple is nonzero, the free cell has a net torque, it would spin, and ``max|P F| → 0`` is unreachable regardless of
preconditioner. This reports, at the frozen seeded state, the net force ``ΣF`` and the net torque about the centre of
mass ``Σ (r_i - r_com) × F_i`` for the SEEDED cell vs a CLEAN (myosin-unbound) cell — isolating the myosin torque.
The scale reference is ``N · f_head · R`` (the torque if every head's couple added coherently). CUDA-gated → gbook.
"""
from __future__ import annotations

import sys

import numpy as np
import warp as wp

from aleph.components.incumbent.assemble import CellConfig, build_cell
from aleph.components.incumbent.driver import _accumulate_all


def net_force_torque(cell):
    _accumulate_all(cell, cell.pos_d, cell.f_d)
    wp.synchronize_device(cell.device)
    F = cell.f_d.numpy()
    r = cell.pos_d.numpy()
    com = r[:cell.n_actin].mean(axis=0)
    net_F = F.sum(axis=0)
    tau = np.cross(r - com, F).sum(axis=0)           # net torque about the actin COM
    return net_F, tau, np.linalg.norm(F, axis=1).max()


def build(native, fraction):
    n = 70686 if native else 1500
    seed = dict(resting_bound_myosin_fraction=fraction, resting_bound_myosin_force_pn=1.5,
                resting_bound_myosin_source="moment_TEST", resting_bound_myosin_capture_um=0.6) if fraction else {}
    return build_cell(CellConfig(n_filaments=n, overlap_free_cortex=True, erm_radial_pairing=True,
                                 membrane_subdivisions=6, **seed))


def main() -> None:
    native = "--native" in sys.argv
    R = 7.5
    print(f"=== GATE A total force/torque (native={native}) ===")
    for tag, frac in (("CLEAN (myosin unbound)", None), ("SEEDED (frac0.5 f1.5pN)", 0.5)):
        cell = build(native, frac)
        net_F, tau, fmax = net_force_torque(cell)
        # coherent-couple reference for the seeded case: n_bound heads × f_head × R
        nb = cell.ledger.get("resting_bound_myosin_n_bound", 0) if frac else 0
        ref = f"  [coherent-couple ref n·f·R ≈ {nb * 1.5 * R:.0f} pN·µm]" if nb else ""
        print(f"  {tag:<26} |ΣF| {np.linalg.norm(net_F):9.3f} pN   |Στ| {np.linalg.norm(tau):11.3f} pN·µm   "
              f"max|F| {fmax:7.3f}{ref}")
    print("Reading: if SEEDED |Στ| ≫ CLEAN and is an appreciable fraction of n·f·R, the directional crossbridge")
    print("leaves a NET TORQUE (frozen walk_dir not rotationally invariant) → no free-cell static equilibrium;")
    print("if SEEDED |Στ| ≈ CLEAN (both ~roundoff·N), the bipolar geometry cancels the couples and equilibrium exists.")


if __name__ == "__main__":
    main()
