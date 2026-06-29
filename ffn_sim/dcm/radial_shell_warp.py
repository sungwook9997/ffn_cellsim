"""Warp port of the radial-shell compartment force (Phase B, B2).

Ports ``FFNRadialShellForce`` — the single shape shared by all three production
compartment ``md.force.Custom`` forces (``cell.nucleus.NucleusConfinement``,
``cell.membrane_surface.MembraneSurfaceTension``,
``cortex.enclosed_volume.EnclosedVolumePressure``) and the native CUDA
``radial_shell_force.cu``. The shape is:

  1. shell centroid  ``r_c = ⟨r_i⟩``           (mean over tagged shell beads)
  2. mean radius      ``R_mean = ⟨|r_i − r_c|⟩``
  3. per-bead radial force ``F_i = F_mag · n̂_i`` from one of three laws:

     law 0 — nucleus  : bilinear about R0 (pa=k_chrom, pb=k_lamin, pc=d_knee, pd=F_knee)
     law 1 — membrane : S=4πR², γ=pa+pb·(S-pc)/pc, ΔP=2γ/R, F=-(ΔP·S/n)·n̂  (pa=γ_mem,pb=K_A,pc=A0)
     law 2 — turgor   : V=4/3πR³, ΔP=pa-pb·(V-pc)/pc, F=+(ΔP·S/n)·n̂        (pa=ΔP0,pb=K_vol,pc=V0)

mirroring ``radial_shell_force.cu`` term-by-term (which itself was validated to
match the numpy ``set_forces`` originals to rel<1e-8). The per-bead force law is
the parallel, expensive part; the global reduction is a standard mean.

Parity contract (the anti-drift discipline)
-------------------------------------------
* ``reduce="host"`` (the GATED path): centroid + R_mean are computed in numpy,
  IDENTICAL to the reference's ``pos.mean(axis=0)`` / ``radii.mean()``, and only
  the per-bead force law runs in the Warp kernel. This isolates and verifies the
  force LAW — deterministic and reproducible. Achieves ~1e-13-or-tighter rel.
* ``reduce="warp"`` (the FULL-port diagnostic): centroid + R_mean are reduced in
  Warp via ``wp.atomic_add`` (mirroring the .cu's two reduction passes). Atomic
  order differs from numpy's pairwise sum, so this is graded only at the .cu's
  own reduction-order tolerance and reported as a diagnostic, not the committed
  bit-parity claim (mirrors "modulo atomic-reduction order" in the native test).
"""

from __future__ import annotations

import math

import numpy as np

import warp as wp

wp.init()

_PI = math.pi


@wp.kernel
def rsf_reduce_centroid(
    pos: wp.array(dtype=wp.vec3d),
    tag: wp.array(dtype=wp.uint32),
    acc: wp.array(dtype=wp.float64),  # [sum_x, sum_y, sum_z, count]
    t0: wp.uint32,
    t1: wp.uint32,
):
    i = wp.tid()
    t = tag[i]
    if t >= t0 and t < t1:
        p = pos[i]
        wp.atomic_add(acc, 0, p[0])
        wp.atomic_add(acc, 1, p[1])
        wp.atomic_add(acc, 2, p[2])
        wp.atomic_add(acc, 3, wp.float64(1.0))


@wp.kernel
def rsf_reduce_radius(
    pos: wp.array(dtype=wp.vec3d),
    tag: wp.array(dtype=wp.uint32),
    acc_r: wp.array(dtype=wp.float64),  # [sum_r]
    cx: wp.float64,
    cy: wp.float64,
    cz: wp.float64,
    t0: wp.uint32,
    t1: wp.uint32,
):
    i = wp.tid()
    t = tag[i]
    if t >= t0 and t < t1:
        p = pos[i]
        dx = p[0] - cx
        dy = p[1] - cy
        dz = p[2] - cz
        wp.atomic_add(acc_r, 0, wp.sqrt(dx * dx + dy * dy + dz * dz))


@wp.kernel
def rsf_apply(
    pos: wp.array(dtype=wp.vec3d),
    tag: wp.array(dtype=wp.uint32),
    force: wp.array(dtype=wp.vec3d),
    energy: wp.array(dtype=wp.float64),
    t0: wp.uint32,
    t1: wp.uint32,
    cx: wp.float64,
    cy: wp.float64,
    cz: wp.float64,
    R_mean: wp.float64,
    cnt: wp.float64,
    law: wp.int32,
    R0: wp.float64,
    pa: wp.float64,
    pb: wp.float64,
    pc: wp.float64,
    pd: wp.float64,
):
    """Per-bead radial-shell force (thread = bead). Mirrors set_forces / the .cu."""
    i = wp.tid()
    t = tag[i]
    if t < t0 or t >= t1:
        force[i] = wp.vec3d(wp.float64(0.0), wp.float64(0.0), wp.float64(0.0))
        energy[i] = wp.float64(0.0)
        return

    p = pos[i]
    dx = p[0] - cx
    dy = p[1] - cy
    dz = p[2] - cz
    r = wp.sqrt(dx * dx + dy * dy + dz * dz)
    rs = r
    if r <= wp.float64(0.0):
        rs = wp.float64(1.0)
    nx = dx / rs
    ny = dy / rs
    nz = dz / rs
    if r <= wp.float64(0.0):
        nx = wp.float64(0.0)
        ny = wp.float64(0.0)
        nz = wp.float64(0.0)

    PI = wp.float64(_PI)
    Fmag = wp.float64(0.0)
    U = wp.float64(0.0)

    if law == 0:
        # nucleus: bilinear about R0 (pa=k_chrom, pb=k_lamin, pc=d_knee, pd=F_knee)
        d = r - R0
        ad = wp.abs(d)
        sgn = wp.float64(0.0)
        if d > wp.float64(0.0):
            sgn = wp.float64(1.0)
        if d < wp.float64(0.0):
            sgn = wp.float64(-1.0)
        if ad <= pc:
            Fmag = -pa * d
            U = wp.float64(0.5) * pa * d * d
        else:
            e = ad - pc
            uknee = wp.float64(0.5) * pa * pc * pc
            Fmag = -sgn * (pd + (pa + pb) * e)
            U = uknee + pd * e + wp.float64(0.5) * (pa + pb) * e * e
    elif law == 1:
        # membrane: S=4πR² Laplace (pa=γ_mem, pb=K_A, pc=A0); inward
        S = wp.float64(4.0) * PI * R_mean * R_mean
        gamma_tot = pa + pb * (S - pc) / pc
        dP = wp.float64(2.0) * gamma_tot / R_mean
        A_i = S / cnt
        Fmag = -(dP * A_i)
        Uc = pa * (S - pc)
        if Uc < wp.float64(0.0):
            Uc = wp.float64(0.0)
        U = (Uc + wp.float64(0.5) * pb * (S - pc) * (S - pc) / pc) / cnt
    else:
        # turgor: V=4/3πR³ (pa=ΔP0, pb=K_vol, pc=V0); outward
        V = (wp.float64(4.0) / wp.float64(3.0)) * PI * R_mean * R_mean * R_mean
        S = wp.float64(4.0) * PI * R_mean * R_mean
        dP = pa - pb * (V - pc) / pc
        A_i = S / cnt
        Fmag = dP * A_i
        U = (wp.float64(0.5) * (pb / pc) * (V - pc) * (V - pc)) / cnt

    force[i] = wp.vec3d(Fmag * nx, Fmag * ny, Fmag * nz)
    energy[i] = U


def run_radial_shell_warp(
    *,
    pos: np.ndarray,        # (N, 3) float64
    tag: np.ndarray,        # (N,) uint32
    tag_range: tuple[int, int],
    law: int,
    R0: float,
    pa: float,
    pb: float,
    pc: float,
    pd: float,
    reduce: str = "host",
    device: str = "cpu",
) -> dict:
    """Run the radial-shell force in Warp. Returns force (N,3), energy (N,), and
    the reduction scalars used.

    ``reduce="host"`` computes centroid/R_mean in numpy (identical to the
    reference) — the gated force-law parity path. ``reduce="warp"`` reduces in
    Warp (atomic add) — the full-port diagnostic.
    """
    pos = np.ascontiguousarray(pos, dtype=np.float64)
    tag = np.ascontiguousarray(tag, dtype=np.uint32)
    N = pos.shape[0]
    t0, t1 = int(tag_range[0]), int(tag_range[1])

    pos_d = wp.array(pos, dtype=wp.vec3d, device=device)
    tag_d = wp.array(tag, dtype=wp.uint32, device=device)

    mask = (tag >= t0) & (tag < t1)
    cnt = float(mask.sum())

    if reduce == "host":
        sub = pos[mask]
        centroid = sub.mean(axis=0)
        radii = np.linalg.norm(sub - centroid, axis=1)
        R_mean = float(radii.mean())
        cx, cy, cz = (float(centroid[0]), float(centroid[1]), float(centroid[2]))
    elif reduce == "warp":
        acc = wp.zeros(4, dtype=wp.float64, device=device)
        wp.launch(rsf_reduce_centroid, dim=N,
                  inputs=[pos_d, tag_d, acc, wp.uint32(t0), wp.uint32(t1)],
                  device=device)
        a = acc.numpy()
        cnt_w = a[3] if a[3] > 0.0 else 1.0
        cx, cy, cz = a[0] / cnt_w, a[1] / cnt_w, a[2] / cnt_w
        acc_r = wp.zeros(1, dtype=wp.float64, device=device)
        wp.launch(rsf_reduce_radius, dim=N,
                  inputs=[pos_d, tag_d, acc_r,
                          wp.float64(cx), wp.float64(cy), wp.float64(cz),
                          wp.uint32(t0), wp.uint32(t1)],
                  device=device)
        R_mean = float(acc_r.numpy()[0] / cnt_w)
        cx, cy, cz = float(cx), float(cy), float(cz)
    else:
        raise ValueError(f"reduce must be 'host' or 'warp', got {reduce!r}")

    force_d = wp.zeros(N, dtype=wp.vec3d, device=device)
    energy_d = wp.zeros(N, dtype=wp.float64, device=device)
    wp.launch(
        rsf_apply, dim=N,
        inputs=[pos_d, tag_d, force_d, energy_d,
                wp.uint32(t0), wp.uint32(t1),
                wp.float64(cx), wp.float64(cy), wp.float64(cz),
                wp.float64(R_mean), wp.float64(cnt),
                wp.int32(law),
                wp.float64(R0), wp.float64(pa), wp.float64(pb),
                wp.float64(pc), wp.float64(pd)],
        device=device,
    )
    wp.synchronize_device(device)
    return {
        "force": force_d.numpy().astype(np.float64),
        "energy": energy_d.numpy().astype(np.float64),
        "centroid": np.array([cx, cy, cz], dtype=np.float64),
        "R_mean": R_mean,
        "cnt": cnt,
    }
