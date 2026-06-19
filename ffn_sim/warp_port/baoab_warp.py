"""Warp port of the overdamped Leimkuhler-Matthews BAOAB step (Phase B, B1).

Mirrors the FROZEN ``integrator.baoab.LeimkuhlerMatthewsBAOAB`` per-step update
(and its device sibling ``integrator.baoab_device.OverdampedBAOABDevice``) as a
single ``warp.kernel``, one thread per particle. This is the smallest piece of
the runtime and the first port gated by bit-parity against a committed HOOMD
reference fixture (see ``WARP_SPIKE_THEN_DECOHESION_PLAN_2026-06-19.md`` §B1).

Update rule (per particle, identical to the frozen Action)::

    r(t+Δt) = r(t) + (F/γ)·Δt + √(kT/(2 γ Δt))·(W_n + W_{n-1})·Δt

then wrap into the (possibly triclinic / Lees-Edwards-sheared) box, accumulate
the integer image-flag delta, and store ``W_n`` as next step's ``W_{n-1}``.

Parity contract (the anti-drift discipline)
-------------------------------------------
* **Float-op ordering matches the reference exactly.** The reference computes
  ``dr = F * (1/γ) * Δt + √(kT/(2γΔt)) * (W_n + W_{n-1}) * Δt``; this kernel
  reproduces that grouping term-by-term so the IEEE-754 result is bit-identical
  on the deterministic part.
* **The box wrap uses ``wp.rint`` (round-half-to-even)** — verified equal to
  numpy's ``np.round`` used by the frozen ``_wrap_into_box`` (``wp.round`` rounds
  half away from zero and would break bit-parity).
* **Noise is INJECTED, not generated here.** numpy's PCG64 ``standard_normal``
  is not reproducible inside a Warp kernel; for a clean parity test the SAME
  per-step noise the reference drew (``np.random.default_rng(seed)``) is passed
  in as a committed array. At kT=0 the noise term is multiplied by 0 so parity
  is bit-for-bit regardless. Generating the noise device-resident (splitmix64
  counter-RNG, as in ``native/.../baoab_kernel.cu``) is a separate follow-up;
  it does not change the integrator arithmetic this kernel verifies.

The kernel writes ``pos`` / ``image`` / ``prv`` in place across K launches with
no host transfer between steps, mirroring the GPU-resident intent of the port.
"""

from __future__ import annotations

import numpy as np

import warp as wp

wp.init()


@wp.kernel
def baoab_step_kernel(
    pos: wp.array(dtype=wp.vec3d),          # (N,) rw — particle positions
    force: wp.array(dtype=wp.vec3d),        # (N,) ro — net force F(r(t))
    image: wp.array(dtype=wp.vec3i),        # (N,) rw — periodic image flags
    gamma: wp.array(dtype=wp.float64),      # (N,) ro — per-particle Stokes drag
    bd_pref: wp.array(dtype=wp.float64),    # (N,) ro — √(kT/(2 γ Δt)) precomputed
    prv: wp.array(dtype=wp.vec3d),          # (N,) rw — W_{n-1}
    noise: wp.array(dtype=wp.vec3d, ndim=2),  # (K, N) ro — injected W_n per step
    k: wp.int32,                            # current step index into noise
    dt: wp.float64,
    Lx: wp.float64,
    Ly: wp.float64,
    Lz: wp.float64,
    xy: wp.float64,
    xz: wp.float64,
    yz: wp.float64,
):
    """One overdamped L-M BAOAB step for particle ``i`` (thread = particle).

    Reproduces ``LeimkuhlerMatthewsBAOAB.act`` arithmetic and ``_wrap_into_box``
    term-by-term so kT=0 is bit-for-bit and kT>0 differs only by float-op order.
    """
    i = wp.tid()

    W = noise[k, i]
    p = prv[i]
    f = force[i]

    inv_gamma = wp.float64(1.0) / gamma[i]
    pref = bd_pref[i]

    # dr = F·(1/γ)·Δt + √(kT/(2γΔt))·(W_n + W_{n-1})·Δt   (grouping matches ref)
    dr = f * inv_gamma * dt + pref * (W + p) * dt
    np_ = pos[i] + dr

    # --- minimum-image wrap into the upper-triangular box (frozen _wrap_into_box) ---
    rx = np_[0]
    ry = np_[1]
    rz = np_[2]

    fz = rz / Lz
    fy = (ry - yz * Lz * fz) / Ly
    fx = (rx - xy * Ly * fy - xz * Lz * fz) / Lx

    nx = wp.rint(fx)  # round-half-to-even == numpy np.round (verified)
    ny = wp.rint(fy)
    nz = wp.rint(fz)

    fx = fx - nx
    fy = fy - ny
    fz = fz - nz

    ox = Lx * fx + xy * Ly * fy + xz * Lz * fz
    oy = Ly * fy + yz * Lz * fz
    oz = Lz * fz

    pos[i] = wp.vec3d(ox, oy, oz)

    img = image[i]
    image[i] = wp.vec3i(
        img[0] + wp.int32(nx),
        img[1] + wp.int32(ny),
        img[2] + wp.int32(nz),
    )

    prv[i] = W


def run_baoab_warp(
    *,
    pos0: np.ndarray,        # (N, 3) float64 initial positions (tag-ordered)
    image0: np.ndarray,      # (N, 3) int32 initial image flags (tag-ordered)
    force: np.ndarray,       # (N, 3) float64 constant net force (tag-ordered)
    gamma: np.ndarray,       # (N,) float64 per-particle Stokes drag (tag-ordered)
    kT: float,
    dt: float,
    noise: np.ndarray,       # (K, N, 3) float64 injected W_n per step (tag-ordered)
    box: tuple[float, float, float, float, float, float],  # (Lx,Ly,Lz,xy,xz,yz)
    prv0: np.ndarray | None = None,  # (N, 3) float64 W_{-1}; default zeros (frozen init)
    device: str = "cpu",
) -> dict[str, np.ndarray]:
    """Run K overdamped L-M BAOAB steps in Warp with INJECTED noise.

    Returns a dict with the final ``pos`` (N,3 float64), ``image`` (N,3 int32),
    and ``prv`` (N,3 float64) — all tag-ordered (row i == tag i). The position
    / image / prv arrays live on the Warp device across all K launches with no
    host transfer between steps (the GPU-resident intent); only the final state
    is copied back to host here.

    The ``force`` is constant across steps — the parity fixture uses
    ``md.force.Constant`` (position-independent) so the integrator step can be
    verified in isolation, exactly as the reference computed it.
    """
    pos0 = np.ascontiguousarray(pos0, dtype=np.float64)
    image0 = np.ascontiguousarray(image0, dtype=np.int32)
    force = np.ascontiguousarray(force, dtype=np.float64)
    gamma = np.ascontiguousarray(gamma, dtype=np.float64)
    noise = np.ascontiguousarray(noise, dtype=np.float64)
    N = pos0.shape[0]
    K = noise.shape[0]
    if noise.shape[1] != N:
        raise ValueError(f"noise N mismatch: noise {noise.shape} vs pos {pos0.shape}")

    # bd_prefactor = √(kT/(2 γ Δt)) — precomputed exactly as the frozen Action
    # does at attach() time (numpy float64), then passed device-resident.
    bd_pref = np.sqrt(kT / (2.0 * gamma * dt))

    if prv0 is None:
        prv0 = np.zeros((N, 3), dtype=np.float64)
    prv0 = np.ascontiguousarray(prv0, dtype=np.float64)

    Lx, Ly, Lz, xy, xz, yz = box

    pos_d = wp.array(pos0, dtype=wp.vec3d, device=device)
    image_d = wp.array(image0, dtype=wp.vec3i, device=device)
    force_d = wp.array(force, dtype=wp.vec3d, device=device)
    gamma_d = wp.array(gamma, dtype=wp.float64, device=device)
    bd_pref_d = wp.array(bd_pref, dtype=wp.float64, device=device)
    prv_d = wp.array(prv0, dtype=wp.vec3d, device=device)
    noise_d = wp.array(noise, dtype=wp.vec3d, device=device)  # (K, N) vec3d

    for k in range(K):
        wp.launch(
            baoab_step_kernel,
            dim=N,
            inputs=[
                pos_d, force_d, image_d, gamma_d, bd_pref_d, prv_d, noise_d,
                wp.int32(k),
                wp.float64(dt),
                wp.float64(Lx), wp.float64(Ly), wp.float64(Lz),
                wp.float64(xy), wp.float64(xz), wp.float64(yz),
            ],
            device=device,
        )
    wp.synchronize_device(device)

    return {
        "pos": pos_d.numpy().astype(np.float64),
        "image": image_d.numpy().astype(np.int32),
        "prv": prv_d.numpy().astype(np.float64),
    }
