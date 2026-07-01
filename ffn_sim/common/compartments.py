"""Shared engine-agnostic cell-compartment forces (Warp) — used by BOTH ff/ and dcm/.

A "compartment" is a mechanical sub-element of the cell that accumulates a per-node radial force from a
shell centroid. The single shared shape (nucleus / membrane / mean-radius turgor) is the radial-shell
force, ported VERBATIM (bit-identical) from ``dcm/radial_shell_warp.py`` (itself validated to match the
HOOMD ``cell.nucleus.NucleusConfinement`` / ``cell.membrane_surface.MembraneSurfaceTension`` /
``cortex.enclosed_volume.EnclosedVolumePressure`` references to rel<1e-8, and the native CUDA
``radial_shell_force.cu``). Placing it in ``common/`` is the engine-agnostic home the two Warp engines
share (per PI 2026-07-01: implement compartments ONCE, both engines import).

  law 0 — nucleus  : bilinear about R0 (pa=k_chrom, pb=k_lamin, pc=d_knee, pd=F_knee)  [stiff strain-stiffening]
  law 1 — membrane : S=4πR², γ=pa+pb·(S-pc)/pc, ΔP=2γ/R, F=-(ΔP·S/n)·n̂  (pa=γ_mem,pb=K_A,pc=A0)  [inward]
  law 2 — turgor   : V=4/3πR³, ΔP=pa-pb·(V-pc)/pc, F=+(ΔP·S/n)·n̂        (pa=ΔP0,pb=K_vol,pc=V0)  [outward]

Units: the FF engine works in pN·µm·s and **1 Pa ≡ 1 pN/µm²** exactly, so the SI compartment formulae
(nucleus.py etc.) hold numerically with moduli in pN/µm² and lengths in µm.

⚠️ Follow-up (PI-coordinate, when the DCM /loop is idle): make ``dcm/radial_shell_warp.py`` a thin
re-export of this module to remove the temporary duplication (kept duplicated now to avoid touching the
live-/loop DCM files). The kernel here is a verbatim copy → bit-identical.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

import warp as wp

wp.init()

_PI = math.pi


@wp.kernel
def rsf_reduce_centroid(pos: wp.array(dtype=wp.vec3d), tag: wp.array(dtype=wp.uint32),
                        acc: wp.array(dtype=wp.float64), t0: wp.uint32, t1: wp.uint32):
    i = wp.tid()
    t = tag[i]
    if t >= t0 and t < t1:
        p = pos[i]
        wp.atomic_add(acc, 0, p[0]); wp.atomic_add(acc, 1, p[1])
        wp.atomic_add(acc, 2, p[2]); wp.atomic_add(acc, 3, wp.float64(1.0))


@wp.kernel
def rsf_reduce_radius(pos: wp.array(dtype=wp.vec3d), tag: wp.array(dtype=wp.uint32),
                      acc_r: wp.array(dtype=wp.float64), cx: wp.float64, cy: wp.float64, cz: wp.float64,
                      t0: wp.uint32, t1: wp.uint32):
    i = wp.tid()
    t = tag[i]
    if t >= t0 and t < t1:
        p = pos[i]
        dx = p[0] - cx; dy = p[1] - cy; dz = p[2] - cz
        wp.atomic_add(acc_r, 0, wp.sqrt(dx * dx + dy * dy + dz * dz))


@wp.kernel
def rsf_apply(pos: wp.array(dtype=wp.vec3d), tag: wp.array(dtype=wp.uint32),
              force: wp.array(dtype=wp.vec3d), energy: wp.array(dtype=wp.float64),
              t0: wp.uint32, t1: wp.uint32, cx: wp.float64, cy: wp.float64, cz: wp.float64,
              R_mean: wp.float64, cnt: wp.float64, law: wp.int32,
              R0: wp.float64, pa: wp.float64, pb: wp.float64, pc: wp.float64, pd: wp.float64):
    """Per-bead radial-shell force (thread = bead). Verbatim from dcm/radial_shell_warp.rsf_apply."""
    i = wp.tid()
    t = tag[i]
    if t < t0 or t >= t1:
        force[i] = wp.vec3d(wp.float64(0.0), wp.float64(0.0), wp.float64(0.0))
        energy[i] = wp.float64(0.0)
        return
    p = pos[i]
    dx = p[0] - cx; dy = p[1] - cy; dz = p[2] - cz
    r = wp.sqrt(dx * dx + dy * dy + dz * dz)
    rs = r
    if r <= wp.float64(0.0):
        rs = wp.float64(1.0)
    nx = dx / rs; ny = dy / rs; nz = dz / rs
    if r <= wp.float64(0.0):
        nx = wp.float64(0.0); ny = wp.float64(0.0); nz = wp.float64(0.0)
    PI = wp.float64(_PI)
    Fmag = wp.float64(0.0)
    U = wp.float64(0.0)
    if law == 0:
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
        V = (wp.float64(4.0) / wp.float64(3.0)) * PI * R_mean * R_mean * R_mean
        S = wp.float64(4.0) * PI * R_mean * R_mean
        dP = pa - pb * (V - pc) / pc
        A_i = S / cnt
        Fmag = dP * A_i
        U = (wp.float64(0.5) * (pb / pc) * (V - pc) * (V - pc)) / cnt
    force[i] = wp.vec3d(Fmag * nx, Fmag * ny, Fmag * nz)
    energy[i] = U


def run_radial_shell_warp(*, pos: np.ndarray, tag: np.ndarray, tag_range: tuple[int, int], law: int,
                          R0: float, pa: float, pb: float, pc: float, pd: float,
                          reduce: str = "host", device: str = "cpu") -> dict:
    """Radial-shell force in Warp (verbatim from dcm/radial_shell_warp.run_radial_shell_warp).

    ``reduce='host'`` = numpy centroid/R_mean (gated bit-parity path); ``reduce='warp'`` = Warp reduction.
    Returns force (N,3), energy (N,), centroid, R_mean, cnt.
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
        R_mean = float(np.linalg.norm(sub - centroid, axis=1).mean())
        cx, cy, cz = float(centroid[0]), float(centroid[1]), float(centroid[2])
    elif reduce == "warp":
        acc = wp.zeros(4, dtype=wp.float64, device=device)
        wp.launch(rsf_reduce_centroid, dim=N, inputs=[pos_d, tag_d, acc, wp.uint32(t0), wp.uint32(t1)], device=device)
        a = acc.numpy()
        cnt_w = a[3] if a[3] > 0.0 else 1.0
        cx, cy, cz = a[0] / cnt_w, a[1] / cnt_w, a[2] / cnt_w
        acc_r = wp.zeros(1, dtype=wp.float64, device=device)
        wp.launch(rsf_reduce_radius, dim=N, inputs=[pos_d, tag_d, acc_r, wp.float64(cx), wp.float64(cy),
                  wp.float64(cz), wp.uint32(t0), wp.uint32(t1)], device=device)
        R_mean = float(acc_r.numpy()[0] / cnt_w)
        cx, cy, cz = float(cx), float(cy), float(cz)
    else:
        raise ValueError(f"reduce must be 'host' or 'warp', got {reduce!r}")
    force_d = wp.zeros(N, dtype=wp.vec3d, device=device)
    energy_d = wp.zeros(N, dtype=wp.float64, device=device)
    wp.launch(rsf_apply, dim=N, inputs=[pos_d, tag_d, force_d, energy_d, wp.uint32(t0), wp.uint32(t1),
              wp.float64(cx), wp.float64(cy), wp.float64(cz), wp.float64(R_mean), wp.float64(cnt),
              wp.int32(law), wp.float64(R0), wp.float64(pa), wp.float64(pb), wp.float64(pc), wp.float64(pd)],
              device=device)
    wp.synchronize_device(device)
    return {"force": force_d.numpy().astype(np.float64), "energy": energy_d.numpy().astype(np.float64),
            "centroid": np.array([cx, cy, cz], dtype=np.float64), "R_mean": R_mean, "cnt": cnt}


# ---- lit-anchored compartment parameter resolvers (physiological-baseline; 1 Pa ≡ 1 pN/µm²) ----

@dataclass(frozen=True, slots=True)
class ResolvedNucleus:
    """Nucleus radial-shell (law 0) params in FF units (pN, µm). Ported from archive H.9 nucleus.py."""
    R_nuc_um: float
    n_beads: int
    k_chrom: float      # pN/µm  = 4π·E_nuc·R_nuc/n_beads
    k_lamin: float      # pN/µm  = (ratio_lamin−1)·k_chrom
    d_knee_um: float    # µm     = knee_strain·R_nuc
    F_knee_pN: float    # pN     = k_chrom·d_knee


def resolve_nucleus(*, R_nuc_um: float, n_beads: int, E_nuc_Pa: float = 5.0e3,
                    ratio_lamin: float = 3.0, knee_strain: float = 0.10) -> ResolvedNucleus:
    """Resolve the nucleus bilinear radial-shell stiffnesses (KU-3.B2; archive/cell/nucleus.py bridge).

    E_nuc = in-situ nuclear Young's modulus 5 kPa (band 1–10 kPa, KU-3.B2.1); ratio_lamin = nucleus:
    cytoplasm stiffness 3× (in-situ 1.4–5×, NOT the 10× isolated value); knee_strain = lamin engagement
    ~10% radial strain. Bridge (grid-invariant: n_beads·k_chrom = 4π·E_nuc·R_nuc is intensive):
        k_chrom = 4π·E_nuc·R_nuc/n_beads,  k_lamin = (ratio−1)·k_chrom,  d_knee = knee_strain·R_nuc,
        F_knee = k_chrom·d_knee.  (1 Pa ≡ 1 pN/µm² → E_nuc in pN/µm², lengths in µm → k in pN/µm.)
    """
    if not (1.0e3 <= E_nuc_Pa <= 1.0e4):
        raise ValueError(f"E_nuc {E_nuc_Pa} Pa outside KU-3.B2.1 in-situ band [1e3,1e4]")
    if not (1.4 <= ratio_lamin <= 5.0):
        raise ValueError(f"ratio_lamin {ratio_lamin} outside in-situ band [1.4,5.0] (10× is isolated)")
    if n_beads <= 0 or R_nuc_um <= 0.0:
        raise ValueError("R_nuc_um and n_beads must be positive")
    k_chrom = 4.0 * math.pi * E_nuc_Pa * R_nuc_um / n_beads
    k_lamin = (ratio_lamin - 1.0) * k_chrom
    d_knee = knee_strain * R_nuc_um
    return ResolvedNucleus(R_nuc_um=R_nuc_um, n_beads=n_beads, k_chrom=k_chrom, k_lamin=k_lamin,
                           d_knee_um=d_knee, F_knee_pN=k_chrom * d_knee)
