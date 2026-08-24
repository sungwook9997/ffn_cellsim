"""GATE A consultation gate: FD-verify the CORRECTED directional crossbridge tangent K = k w wᵀ.

The implicit operator's crossbridge tangent was fixed from the stale isotropic ``k I`` (central spring, ``d/|d|``
axis) to the exact rank-one ``K = k w wᵀ`` that matches the directional force ``f = k (d·w − r0) w``. The consultation
requires a post-change FD gate before trusting the solver. This isolates the crossbridge force
(``crossbridge_segment_kernel``) and its operator tangent (``add_segment_crossbridge_stiffness_kernel``) and checks,
at a frozen native seeded state, that the tangent equals the true force directional derivative:

    operator action  K·v  ==  −∂F/∂x · v  ==  −[F(x+εv) − F(x−εv)] / (2ε)

for probe fields v: random, rigid translation (⇒ 0), head perturbation ⊥ w (⇒ 0, the directional force ignores ⊥
motion), head perturbation ∥ w (⇒ nonzero, must match FD), plus a symmetry check vᵀKw = wᵀKv. CUDA-gated → gbook.
"""
from __future__ import annotations

import sys

import numpy as np
import warp as wp

from aleph.components.incumbent.assemble import CellConfig, build_cell
from aleph.components.incumbent.implicit_mechanics import add_segment_crossbridge_stiffness_kernel
from aleph.components.motor.segment_motor import crossbridge_segment_kernel

EPS = 1.0e-6


def main() -> None:
    native = "--native" in sys.argv
    n = 70686 if native else 1500
    cell = build_cell(CellConfig(
        n_filaments=n, overlap_free_cortex=True, erm_radial_pairing=True, membrane_subdivisions=6,
        resting_bound_myosin_fraction=0.5, resting_bound_myosin_force_pn=1.5,
        resting_bound_myosin_source="tangent_fd_TEST", resting_bound_myosin_capture_um=0.6,
    ))
    d = cell.device
    st = cell.myosin.segment_runtime.state
    hn = cell.myosin.head_node
    k_xb, r0_xb = cell.myosin.params.k_xb, cell.myosin.params.r0_xb
    nb = int(st["bound"].shape[0])
    pos0 = cell.pos_d.numpy()
    walk = st["walk_dir"].numpy()
    bound = st["bound"].numpy().astype(bool)
    head_idx = hn.numpy()

    def xb_force(pos_np):
        p = wp.array(np.ascontiguousarray(pos_np), dtype=wp.vec3d, device=d)
        f = wp.zeros(cell.n_total, dtype=wp.vec3d, device=d)
        wp.launch(crossbridge_segment_kernel, dim=nb,
                  inputs=[p, hn, st["bound"], st["seg_a"], st["seg_b"], st["bary_t"], st["abscissa"],
                          st["walk_dir"], k_xb, r0_xb, f], device=d)
        wp.synchronize_device(d)
        return f.numpy()

    def xb_tangent(v_np):
        v = wp.array(np.ascontiguousarray(v_np), dtype=wp.vec3d, device=d)
        out = wp.zeros(cell.n_total, dtype=wp.vec3d, device=d)
        wp.launch(add_segment_crossbridge_stiffness_kernel, dim=nb,
                  inputs=[cell.pos_d, v, hn, st["bound"], st["seg_a"], st["seg_b"], st["bary_t"], st["abscissa"],
                          st["walk_dir"], k_xb, r0_xb, out], device=d)
        wp.synchronize_device(d)
        return out.numpy()

    rng = np.random.default_rng(0)
    v_rand = rng.standard_normal((cell.n_total, 3))
    v_trans = np.tile([1.0, 0.7, -0.3], (cell.n_total, 1))
    v_perp = np.zeros((cell.n_total, 3)); v_par = np.zeros((cell.n_total, 3))
    hb = head_idx[bound]
    wb = walk[bound]
    perp = np.cross(wb, [0.0, 0.0, 1.0]); perp /= (np.linalg.norm(perp, axis=1, keepdims=True) + 1e-12)
    v_perp[hb] = perp
    v_par[hb] = wb

    print(f"=== crossbridge tangent FD gate (native={native}, n_bound={nb}) ===")
    print(f"{'probe field':<28}{'||K·v||':>12}{'||FD||':>12}{'rel err':>12}")
    for name, v in [("random (all nodes)", v_rand), ("rigid translation ⇒0", v_trans),
                    ("head ⊥ w ⇒0", v_perp), ("head ∥ w (must match)", v_par)]:
        Kv = xb_tangent(v)
        fd = (xb_force(pos0 + EPS * v) - xb_force(pos0 - EPS * v)) / (2.0 * EPS)   # ∂F/∂x · v
        target = -fd                                                              # operator action = −∂F/∂x·v
        rel = np.linalg.norm(Kv - target) / (np.linalg.norm(target) + 1e-12)
        print(f"{name:<28}{np.linalg.norm(Kv):>12.4f}{np.linalg.norm(fd):>12.4f}{rel:>12.3e}")

    # symmetry vᵀKw = wᵀKv
    v1, v2 = rng.standard_normal((cell.n_total, 3)), rng.standard_normal((cell.n_total, 3))
    s = abs(float((v1 * xb_tangent(v2)).sum()) - float((v2 * xb_tangent(v1)).sum()))
    print(f"symmetry |v1ᵀKv2 − v2ᵀKv1| = {s:.3e}  (PSD rank-one ⇒ ~roundoff)")
    print("GATE: rel err for random & ∥w near machine-eps, ⊥w & translation ⇒ ~0, symmetry ⇒ ~0  ⇒ tangent CORRECT")


if __name__ == "__main__":
    main()
