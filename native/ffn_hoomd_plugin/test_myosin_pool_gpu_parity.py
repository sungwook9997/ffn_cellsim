"""GPU parity: native FFNAttachmentSpringForce (fed from MyosinAttachmentPool)
vs the EXACT production myosin attach-bond force. GPU-only; mirrors the CPU gate.

This is the GBOOK side of the Stage-2 binder gate. It cannot run on a CPU-only
box (``FFNAttachmentSpringForce`` requires a GPU). It reuses the SAME deterministic
binding-state setup as ``test_myosin_pool_cpu_parity.py`` and asserts that the
compiled CUDA pool force equals the production ``md.bond.Harmonic`` attach-bond
force to tight tolerance, AND equals the host ``reference_forces`` (kernel == CPU
mirror).

# GBOOK-VALIDATE: run on gbook (RTX A5000, HOOMD-blue 7.0.1) after building the
#   plugin. Exact commands:
#     cd ~/ffn_cellsim/native/ffn_hoomd_plugin
#     cmake -S . -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build -j
#     PYTHONPATH="$PWD/build:$PWD/python:$HOME/ffn_cellsim" \
#         python test_myosin_pool_gpu_parity.py
#   Expect: "myosin-pool-gpu-parity PASS" with rel < 1e-9 on both force checks.
"""

from __future__ import annotations

from copy import deepcopy
from pathlib import Path

import numpy as np
import yaml
import hoomd
import hoomd.md as md
import gsd.hoomd

from ffn_sim.cortex.myosin import (
    MyosinStepUpdater,
    cortex_myosin_attach_bin_names,
    cortex_myosin_attach_bin_rest_lengths,
    extend_state_with_cortex_myosin,
    generate_cortex_myosin_layout,
    resolve_cortex_myosin,
)

# On gbook the package __init__ imports the compiled _ffn_native fine.
from ffn_hoomd_plugin import MyosinAttachmentPool, NativeAttachmentSpringForce


_CFG = Path(__file__).resolve().parents[2] / "ffn_sim" / "configs" / "phase1_h3.yaml"


def _p_myo(n_motors, dt):
    cfg = deepcopy(yaml.safe_load(open(_CFG)))
    cfg["cortex"]["myosin"]["stepping_mode"] = "grip_walk"
    cfg["cortex"]["myosin"]["n_motors_per_cell"] = n_motors
    return resolve_cortex_myosin(cfg, dt=dt)


def _build_actin(n_fil, beads_per_fil, ell0, L):
    n_actin = n_fil * beads_per_fil
    pos = np.zeros((n_actin, 3))
    rng = np.random.default_rng(7)
    for f in range(n_fil):
        base = np.array([-0.5 * (beads_per_fil - 1) * ell0,
                         (f - 0.5 * (n_fil - 1)) * 0.6 * ell0,
                         (rng.random() - 0.5) * 0.2 * ell0])
        for j in range(beads_per_fil):
            pos[f * beads_per_fil + j] = base + j * ell0 * np.array([1.0, 0, 0])
    s = gsd.hoomd.Frame()
    s.particles.N = n_actin
    s.particles.types = ["actin_cortex"]
    s.particles.typeid = np.zeros(n_actin, np.uint32)
    s.particles.position = pos
    s.particles.mass = np.ones(n_actin)
    bonds = [(f * beads_per_fil + j, f * beads_per_fil + j + 1)
             for f in range(n_fil) for j in range(beads_per_fil - 1)]
    bg = np.asarray(bonds, np.uint32)
    s.bonds.N = int(bg.shape[0]); s.bonds.types = ["actin_backbone"]
    s.bonds.group = bg; s.bonds.typeid = np.zeros(bg.shape[0], np.uint32)
    s.configuration.box = [L, L, L, 0, 0, 0]
    return s, n_actin, pos


def _to_snap(frame):
    s = hoomd.Snapshot()
    s.configuration.box = list(frame.configuration.box)
    s.particles.N = int(frame.particles.N)
    s.particles.types = list(frame.particles.types)
    s.particles.typeid[:] = np.asarray(frame.particles.typeid)
    s.particles.position[:] = np.asarray(frame.particles.position)
    s.particles.mass[:] = np.asarray(frame.particles.mass)
    s.bonds.N = int(frame.bonds.N)
    s.bonds.types = list(frame.bonds.types)
    if int(frame.bonds.N) > 0:
        s.bonds.group[:] = np.asarray(frame.bonds.group)
        s.bonds.typeid[:] = np.asarray(frame.bonds.typeid)
    return s


def main() -> int:
    n_fil, beads_per_fil, L, ell0, dt, n_motors = 6, 5, 1.0e-5, 0.45e-6, 1.0e-9, 4
    p = _p_myo(n_motors, dt)
    p.head_actin_max_bind_dist = 1.5e-6
    p.head_actin_capture_perp = 1.5e-6

    base, n_actin, apos = _build_actin(n_fil, beads_per_fil, ell0, L)
    layout = generate_cortex_myosin_layout(
        p, R_cell=0.5 * L, motor_tag_start=n_actin, rng=np.random.default_rng(1))
    rng = np.random.default_rng(3)
    chosen = rng.choice(n_actin, size=n_motors, replace=False)
    for m in range(n_motors):
        shift = apos[chosen[m]] - layout.centers[m]
        layout.positions[m] += shift
        layout.centers[m] = apos[chosen[m]]
    frame = extend_state_with_cortex_myosin(base, layout, p)
    snap = _to_snap(frame)
    pos = np.asarray(snap.particles.position, float).copy()
    box_L = np.asarray(snap.configuration.box[:3], float)
    N = int(snap.particles.N)

    upd = MyosinStepUpdater(
        p_myo=p, layout=layout, kT=1.380649e-23 * 310.0, n_cortex_actin=n_actin,
        ell0_cortex=ell0, cortex_beads_per_filament=beads_per_fil)
    n_heads_total = 2 * p.n_heads_per_side * n_motors
    from scipy.spatial import cKDTree
    tree = cKDTree(pos[:n_actin])
    attach_bonds, attach_bins = [], []
    bin_width = p.head_actin_max_bind_dist / p.n_bins
    for h in range(n_heads_total):
        htag = upd._head_global_tag(h)
        d, j = tree.query(pos[htag])
        if d <= p.head_actin_max_bind_dist:
            upd._head_bound_to_actin[h] = int(j)
            fil, pj = upd._tag_to_fil_pos(int(j))
            upd._head_bound_filament[h] = fil
            upd._head_bound_bead_pos[h] = pj
            d_use = min(d, p.head_actin_max_bind_dist - 1e-12)
            attach_bonds.append((htag, int(j)))
            attach_bins.append(int(min(p.n_bins - 1, max(0, int(d_use / bin_width)))))
    attach_bonds = np.asarray(attach_bonds, np.int64)
    attach_bins = np.asarray(attach_bins, np.int64)
    n_engaged = int(attach_bonds.shape[0])
    print(f"[setup] engaged heads: {n_engaged}", flush=True)

    bin_names = cortex_myosin_attach_bin_names(p.n_bins)
    bin_r0 = cortex_myosin_attach_bin_rest_lengths(
        p.n_bins, p.head_actin_max_bind_dist, stepping_mode="grip_walk")
    r0_eps = float(bin_r0[0])

    # GROUND TRUTH: production md.bond.Harmonic over the attach bonds, on GPU.
    dev = hoomd.device.GPU(notice_level=0)
    gt = hoomd.Snapshot()
    gt.configuration.box = list(snap.configuration.box)
    gt.particles.N = N
    gt.particles.types = list(snap.particles.types)
    gt.particles.typeid[:] = np.asarray(snap.particles.typeid)
    gt.particles.position[:] = pos
    gt.particles.mass[:] = np.asarray(snap.particles.mass)
    gt.bonds.N = n_engaged
    gt.bonds.types = list(bin_names)
    gt.bonds.group[:] = attach_bonds.astype(np.uint32)
    gt.bonds.typeid[:] = attach_bins.astype(np.uint32)
    sg = hoomd.Simulation(device=dev, seed=1)
    sg.create_state_from_snapshot(gt)
    hb = md.bond.Harmonic()
    for nm, r0v in zip(bin_names, bin_r0):
        hb.params[nm] = dict(k=p.k_head_actin, r0=float(r0v))
    sg.operations.integrator = md.Integrator(dt=dt, forces=[hb], methods=[])
    sg.run(0)
    F_truth = np.asarray(hb.forces, float).copy()

    # NATIVE: pool -> FFNAttachmentSpringForce (no attach bonds in topology).
    sn = hoomd.Simulation(device=dev, seed=1)
    sn.create_state_from_snapshot(snap)  # snap has NO attach bonds at construction
    pool = MyosinAttachmentPool(pool_size=n_heads_total,
                                k_head_actin=p.k_head_actin, r0=r0_eps)
    pool.update_from_updater(upd)
    nat = NativeAttachmentSpringForce(pool_size=n_heads_total)
    sn.operations.integrator = md.Integrator(dt=dt, forces=[nat], methods=[])
    pool.push_to(nat)            # off-hot-path device-array write
    sn.run(0)
    F_nat = np.asarray(nat.forces, float).copy()

    # CPU reference too (kernel == host mirror).
    F_ref = pool.reference_forces(pos, box_L)

    sF = float(np.max(np.abs(F_truth))) + 1e-300
    d1 = float(np.max(np.abs(F_truth - F_nat))) / sF
    d2 = float(np.max(np.abs(F_ref - F_nat))) / sF
    ok = (d1 < 1e-9) and (d2 < 1e-9)
    print(f"[parity] native-vs-production rel={d1:.2e}  "
          f"native-vs-host-ref rel={d2:.2e}  {'PASS' if ok else 'FAIL'}", flush=True)
    print("myosin-pool-gpu-parity " + ("PASS" if ok else "FAIL"), flush=True)
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
