"""Parity + toggle + speed for FFNAttachmentSpringForce (fixed-pool binder enabler).

(1) Parity: native attachment springs vs an equivalent md.bond.Harmonic over the
    SAME K pairs (uniform k, r0) — force must match to ~eps.
(2) Toggle: set_attachments with k=0 on half the pool zeros those springs with NO
    set_snapshot / bond mutation; time the per-step force + the toggle.
(3) Context: the global set_snapshot it replaces was 51.9 ms/firing (gate_a_profile).

    cd ~/ffn_cellsim/native/ffn_hoomd_plugin
    PYTHONPATH="$PWD/build:$PWD/python:$HOME/ffn_cellsim" python test_attachment_spring.py
"""

from __future__ import annotations

import time

import numpy as np
import hoomd
import hoomd.md as md

from ffn_hoomd_plugin import NativeAttachmentSpringForce


def main() -> int:
    N, K = 6000, 3000
    L, k, r0 = 1.0e-5, 2.0e-3, 1.0e-7
    rng = np.random.default_rng(4)
    pos = (rng.random((N, 3)) - 0.5) * (0.5 * L)
    # K disjoint (head, actin) pairs from distinct particles
    perm = rng.permutation(N)[: 2 * K]
    head = perm[:K].astype(np.int32)
    actin = perm[K:].astype(np.int32)

    def _snap(with_bonds):
        s = hoomd.Snapshot()
        if s.communicator.rank == 0:
            s.configuration.box = [L, L, L, 0, 0, 0]
            s.particles.N = N
            s.particles.types = ["A"]
            s.particles.position[:] = pos
            s.particles.typeid[:] = np.zeros(N, np.int32)
            if with_bonds:
                s.bonds.N = K
                s.bonds.types = ["spring"]
                s.bonds.group[:] = np.stack([head, actin], axis=1)
                s.bonds.typeid[:] = np.zeros(K, np.int32)
        return s

    dev = hoomd.device.GPU(notice_level=0)

    # --- HOOMD harmonic-bond reference ---
    simb = hoomd.Simulation(device=dev, seed=1)
    simb.create_state_from_snapshot(_snap(True))
    hb = md.bond.Harmonic()
    hb.params["spring"] = dict(k=k, r0=r0)
    simb.operations.integrator = md.Integrator(dt=1e-9, forces=[hb], methods=[])
    simb.run(0)
    F_ref = np.asarray(hb.forces, np.float64).copy()

    # --- native attachment-spring ---
    simn = hoomd.Simulation(device=dev, seed=1)
    simn.create_state_from_snapshot(_snap(False))
    nat = NativeAttachmentSpringForce(pool_size=K)
    nat.set_attachments(head, actin, np.full(K, k), np.full(K, r0))
    simn.operations.integrator = md.Integrator(dt=1e-9, forces=[nat], methods=[])
    simn.run(0)
    F_nat = np.asarray(nat.forces, np.float64).copy()

    dF = float(np.max(np.abs(F_ref - F_nat)))
    sF = float(np.max(np.abs(F_ref))) + 1e-300
    okp = dF < 1e-9 * sF
    print(f"[parity] K={K} springs vs md.bond.Harmonic: max|nat-ref|={dF:.3e} rel={dF/sF:.2e}  "
          f"{'PASS' if okp else 'FAIL'}", flush=True)

    # --- toggle: k=0 on the first half (NO set_snapshot). Correctness checked
    #     element-wise vs a HOOMD bond reference built with ONLY the active half. ---
    on = slice(K // 2, K)
    simb2 = hoomd.Simulation(device=dev, seed=1)
    s2 = _snap(False)
    if s2.communicator.rank == 0:
        s2.bonds.N = K - K // 2
        s2.bonds.types = ["spring"]
        s2.bonds.group[:] = np.stack([head[on], actin[on]], axis=1)
        s2.bonds.typeid[:] = np.zeros(K - K // 2, np.int32)
    simb2.create_state_from_snapshot(s2)
    hb2 = md.bond.Harmonic(); hb2.params["spring"] = dict(k=k, r0=r0)
    simb2.operations.integrator = md.Integrator(dt=1e-9, forces=[hb2], methods=[])
    simb2.run(0)
    F_ref2 = np.asarray(hb2.forces, np.float64).copy()

    k_half = np.full(K, k); k_half[: K // 2] = 0.0
    t0 = time.perf_counter()
    nat.set_attachments(head, actin, k_half, np.full(K, r0))
    t_toggle = 1e6 * (time.perf_counter() - t0)
    simn.run(1)  # methods=[] -> positions fixed; forces recompute with the new pool
    F_half = np.asarray(nat.forces, np.float64).copy()
    dT = float(np.max(np.abs(F_ref2 - F_half)))
    okt = dT < 1e-9 * sF

    # --- per-step speed ---
    simn.run(40)
    t0 = time.perf_counter(); simn.run(600); t_step = 1e6 * (time.perf_counter() - t0) / 600
    print(f"[toggle] set_attachments (k=0 half) {t_toggle:.1f} us (NO set_snapshot; was 51,900 us)  "
          f"active-half vs HOOMD ref rel {dT/sF:.1e}  {'PASS' if okt else 'FAIL'}", flush=True)
    print(f"[speed]  per-step force {t_step:.1f} us/step (vs 51.9 ms/firing global set_snapshot)", flush=True)
    ok = okp and okt
    print("attachment-spring " + ("PASS" if ok else "FAIL"), flush=True)
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
