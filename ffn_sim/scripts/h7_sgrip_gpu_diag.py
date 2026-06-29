"""H.7 — minimal repro: myosin grip_walk s_grip accumulates on CPU, NOT on GPU.

See ``docs/v2_audit/H7_GPU_MYOSIN_SGRIP_BUG_2026-06-07.md``. On a CPU device the
cortex-myosin grip-stretch ``_head_grip_s`` grows (~14 pm/tick); on a GPU device
it stays exactly 0.0 despite heads binding and a nonzero Hill ``v_step`` — which
makes the Gate-A experiment invalid on GPU. This is the device-contrast repro for
the fix (instrument ``cortex/myosin.py`` Step-3 walk write-back on GPU).

Usage:
    python -m ffn_sim.scripts.h7_sgrip_gpu_diag cpu
    python -m ffn_sim.scripts.h7_sgrip_gpu_diag gpu     # gbook (CUDA)
"""

from __future__ import annotations

import sys
from copy import deepcopy

import numpy as np

import hoomd

from ffn_sim.archive.hoomd_legacy.cell.manifest import build_baseline_cell, load_manifest


def main() -> int:
    dev_name = sys.argv[1] if len(sys.argv) > 1 else "cpu"
    dev = (hoomd.device.GPU(notice_level=0) if dev_name == "gpu"
           else hoomd.device.CPU(notice_level=0))
    m = deepcopy(load_manifest("mcf7_baseline.yaml"))
    m.setdefault("cortex_overrides", {}).setdefault("cortex", {})["n_filaments"] = 100
    cell = build_baseline_cell(
        manifest=m, device=dev, seed=1, constrained=True, equilibrate=True,
        equilibrate_steps=300, equilibrate_softstart_steps=80,
        allow_unpressurized_dev=True,
    )
    sim, ma = cell.simulation, cell.myosin_action
    print(f"device={dev_name} stepping={ma.stepping_mode} "
          f"v0={ma.p.v0_per_head:.2e} batch_dt={ma.p.batch_dt:.2e}", flush=True)
    for burst in range(4):
        sim.run(500)
        bnd = np.asarray(ma._head_bound_to_actin) >= 0
        sg = np.asarray(ma._head_grip_s)
        snap = sim.state.get_snapshot()
        bt = list(snap.bonds.types)
        tid = np.asarray(snap.bonds.typeid)
        aids = [i for i, n in enumerate(bt) if n.startswith("cortex_myosin_attach")]
        na = int(np.isin(tid, aids).sum()) if aids else 0
        print(f"  burst{burst} n_bound={int(bnd.sum())} attach_snap={na} "
              f"sgmax={sg.max():.3e} sgmean_bound="
              f"{(sg[bnd].mean() if bnd.any() else 0.0):.3e} "
              f"n_step_adv={ma._n_step_advances_total} n_bind={ma._n_bind_total}",
              flush=True)
    sgmax = float(np.asarray(ma._head_grip_s).max())
    verdict = ("s_grip ACCUMULATES (expected)" if sgmax > 0.0
               else "s_grip STUCK AT 0.0 (BUG — GPU grip_walk write-back)")
    print(f"\n[{dev_name}] final sgmax={sgmax:.3e} → {verdict}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
