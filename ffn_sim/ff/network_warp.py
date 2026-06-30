"""GPU-native FF network force kernels (Warp) — the cortex/network forces, device-agnostic.

The whole point of Warp (vs Cytosim's C++ CPU core) is GPU-native + differentiable execution: the
same kernel runs on CPU (Mac dev) or CUDA (gbook A5000) by passing ``device``. The FF bending force
(``forces_warp``) and the implicit CG (``dcm.dcm_warp_implicit.device_cg``) are already Warp; this
module adds the remaining FF network forces — crosslinker/actin axial springs, myosin contractile
links, and the radial osmotic turgor — as Warp kernels (one thread per element, ``wp.atomic_add``,
float64), matching the DCM ``network_warp`` pattern. Together they let the FF cortex relaxation run
fully on-device (no numpy round-trip), so gbook's GPU is used via ``device="cuda"`` — NOT a separate
"GPU port".

Validation: each kernel is bit-parity-checked against the original numpy force in
``gamma_floor`` / ``network_contractility`` (`tests/ff/test_network_warp.py`).
"""

from __future__ import annotations

import numpy as np
import warp as wp

wp.init()


@wp.kernel
def link_spring_kernel(
    pos: wp.array(dtype=wp.vec3d),              # (N,) node positions
    links: wp.array(dtype=wp.int32, ndim=2),    # (L, 2) node-index pairs (i, j)
    k_arr: wp.array(dtype=wp.float64),          # (L,) stiffness [pN/µm]
    r0_arr: wp.array(dtype=wp.float64),         # (L,) rest length [µm]
    force: wp.array(dtype=wp.vec3d),            # (N,) out (atomic accumulate)
):
    """Hookean link force k·(L−r0)·û on i (toward j when stretched), −on j. Crosslinker + actin axial."""
    t = wp.tid()
    i = links[t, 0]
    j = links[t, 1]
    d = pos[j] - pos[i]
    L = wp.length(d)
    if L > wp.float64(1e-12):
        f = (k_arr[t] * (L - r0_arr[t]) / L) * d
        wp.atomic_add(force, i, f)
        wp.atomic_add(force, j, -f)


@wp.kernel
def myosin_kernel(
    pos: wp.array(dtype=wp.vec3d),
    links: wp.array(dtype=wp.int32, ndim=2),    # (M, 2) myosin link node pairs
    f_myo: wp.float64,                          # constant contractile force [pN]
    force: wp.array(dtype=wp.vec3d),
):
    """Myosin contractile force: constant magnitude f_myo pulling i,j together (toward each other)."""
    t = wp.tid()
    i = links[t, 0]
    j = links[t, 1]
    d = pos[j] - pos[i]
    L = wp.length(d)
    if L > wp.float64(1e-12):
        f = (f_myo / L) * d
        wp.atomic_add(force, i, f)              # on i toward j (shorten)
        wp.atomic_add(force, j, -f)


@wp.kernel
def turgor_kernel(
    pos: wp.array(dtype=wp.vec3d),
    centre: wp.vec3d,
    dP_area: wp.float64,                         # ΔP · (area/node) [pN]
    force: wp.array(dtype=wp.vec3d),
):
    """Radial outward osmotic turgor force ΔP·(area/node)·r̂ per node (Young-Laplace shell)."""
    i = wp.tid()
    d = pos[i] - centre
    r = wp.length(d)
    if r > wp.float64(1e-12):
        wp.atomic_add(force, i, (dP_area / r) * d)


@wp.kernel
def axpy_kernel(x: wp.array(dtype=wp.vec3d), step: wp.float64, F: wp.array(dtype=wp.vec3d)):
    """Explicit overdamped step x += step·F (step = dt/γ)."""
    i = wp.tid()
    x[i] = x[i] + step * F[i]


def link_spring_force_np(pos, links, k, r0, device="cpu"):
    """Convenience: Warp link-spring force as a numpy (N,3) array (for parity tests / numpy callers)."""
    N = pos.shape[0]
    pos_d = wp.array(np.ascontiguousarray(pos, dtype=np.float64), dtype=wp.vec3d, device=device)
    force_d = wp.zeros(N, dtype=wp.vec3d, device=device)
    if links.shape[0]:
        wp.launch(link_spring_kernel, dim=links.shape[0],
                  inputs=[pos_d, wp.array(np.ascontiguousarray(links, np.int32), dtype=wp.int32, device=device),
                          wp.array(np.ascontiguousarray(k, np.float64), dtype=wp.float64, device=device),
                          wp.array(np.ascontiguousarray(r0, np.float64), dtype=wp.float64, device=device),
                          force_d], device=device)
        wp.synchronize_device(device)
    return force_d.numpy().astype(np.float64)


def myosin_force_np(pos, links, f_myo, device="cpu"):
    """Convenience: Warp myosin contractile force as a numpy (N,3) array."""
    N = pos.shape[0]
    pos_d = wp.array(np.ascontiguousarray(pos, dtype=np.float64), dtype=wp.vec3d, device=device)
    force_d = wp.zeros(N, dtype=wp.vec3d, device=device)
    if links.shape[0]:
        wp.launch(myosin_kernel, dim=links.shape[0],
                  inputs=[pos_d, wp.array(np.ascontiguousarray(links, np.int32), dtype=wp.int32, device=device),
                          wp.float64(f_myo), force_d], device=device)
        wp.synchronize_device(device)
    return force_d.numpy().astype(np.float64)
