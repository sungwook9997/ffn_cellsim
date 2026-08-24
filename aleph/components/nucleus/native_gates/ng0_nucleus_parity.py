"""NG-0 (I2 nucleus) — Warp-CUDA <-> NumPy device-parity gate (run on the gbook A5000 by the lead).

Certifies that each device nuclear-envelope force kernel IS its CPU-verified per-node force formula. It
builds ONE small oblate icosphere mesh at a fixed subdivision, PERTURBS it with a fixed seed, and drives
BOTH the Warp mesh-force kernel and a pure-NumPy recompute of the SAME per-node force from the IDENTICAL
config, comparing per-node forces to round-off:

    NG-0.n1  helfrich_bending_kernel        vs  lamina_analytic.dihedral_bending_forces (the exact numpy
             (REUSED from ff/membrane_surface)   formula the kernel implements — kernel-on-CUDA vs numpy)
    NG-0.n2  lamina_areal_tension_kernel    vs  numpy per-face 3-regime area-gradient tension
    NG-0.n3  nucleoplasm_volume_force_kernel vs  p_vol * geometry.mesh_volume_gradient
    NG-0.n4  nucleoplasm_volume_reduce_kernel vs geometry.mesh_volume (divergence-theorem scalar)
    NG-0.n5  chromatin_wlc_kernel           vs  numpy Marko-Siggia WLC link tension
    NG-0.n6  nucleoplasm_viscosity_kernel   vs  numpy -gamma_nuc * v
    NG-0.n7  linc_tether_kernel             vs  numpy nonlinear tension-only tether

PASS iff the relative residual max|warp - ref| / max|ref| < 1e-10 (round-off). The device kernels use
``wp.atomic_add`` across faces/hinges/links, so the accumulation order differs from numpy's — for a handful
of f64 terms per node this reorders at the ~1e-15 level, far below 1e-10; a real formula bug would show as a
gross residual. STRUCTURAL / magnitude-independent: every modulus / kappa / k here is a test value (a
round-off parity verdict is param-agnostic; the physiological magnitudes are I0-B2 GAP-PI). Do NOT tune,
loosen the tolerance, or edit a kernel to force a pass — a genuine NG-0 failure is a real on-device bug.

I0-A: building the Warp arrays + launching the kernels REQUIRES a CUDA device, so this runner only executes
on the gbook A5000 (never the dev Mac; the Mac runs only the pure-NumPy oracles).

Run (from ~/ffn_ac_native on gbook):
    PYTHONPATH=. ~/miniconda3/envs/ffn_sim/bin/python aleph/components/nucleus/native_gates/ng0_nucleus_parity.py
"""

from __future__ import annotations

import sys

import numpy as np
import warp as wp

from aleph.components.nucleus.envelope import (
    chromatin_wlc_kernel,
    face_reference_areas,
    helfrich_bending_kernel,
    lamina_areal_tension_kernel,
    linc_tether_kernel,
    nucleoplasm_viscosity_kernel,
    nucleoplasm_volume_force_kernel,
    nucleoplasm_volume_reduce_kernel,
)
from aleph.components.nucleus.geometry import (
    build_hinges,
    build_oblate_mesh,
    mesh_volume,
    mesh_volume_gradient,
)
from aleph.components.nucleus.lamina_analytic import calibrate_kappa_tilde, dihedral_bending_forces

# ----------------------------------------------------------------------------------------------------
# Fixed test case (deterministic; STRUCTURAL — magnitudes irrelevant to a round-off parity verdict).
# ----------------------------------------------------------------------------------------------------
DEVICE: str | None = None  # resolve the current CUDA device; the hardware contract forbids fixed ordinals
SEED = 20260717
RTOL = 1.0e-10

R_EQ = 5.0              # equatorial radius [um]
ASPECT = 1.5            # oblate a/c
SUBDIV = 2              # 162 nodes / 320 faces / 480 hinges (small but generic)
PERTURB = 0.20         # random node displacement [um] -> generic non-flat / strained config

# test moduli / constants (round-off verdict is param-agnostic)
KAPPA_NE = 0.02         # envelope bending rigidity [pN.um]
K_SOFT = 10.0           # lamina sub-knee areal modulus [pN/um]
K_AC = 30.0             # lamina-A/C stiff areal modulus [pN/um]
KNEE = 0.02             # areal-strain knee
EPS_RUPT = 0.10         # areal-strain rupture threshold (sigma -> 0 above it)
P_VOL = 3.0             # nucleoplasm pressure [pN/um^2]
CONTOUR_L = 2.0         # chromatin WLC contour length [um]
PERSIST = 0.05          # chromatin persistence length [um]
KBT = 4.28e-3           # kBT at 310 K [pN.um]
GAMMA_NUC = 7.5         # nucleoplasm drag [pN.s/um]
K_LINC = 8.0            # LINC stiffness [pN/um]
STIFFENING = 2.0        # LINC nonlinear stiffening [1/um^2]


def _residual(warp_arr: np.ndarray, ref_arr: np.ndarray) -> tuple[float, float]:
    warp_arr = np.asarray(warp_arr, dtype=np.float64)
    ref_arr = np.asarray(ref_arr, dtype=np.float64)
    max_abs = float(np.max(np.abs(warp_arr - ref_arr)))
    scale = float(np.max(np.abs(ref_arr)))
    rel = max_abs / scale if scale > 0.0 else max_abs
    return max_abs, rel


def _report(name: str, max_abs: float, rel: float) -> bool:
    passed = rel < RTOL
    print(f"NG-0 {name}: max_abs={max_abs:.3e} rel={rel:.3e} -> {'PASS' if passed else 'FAIL'}")
    return passed


# ---- numpy mirrors of the per-node kernel formulas (SAME formula, host recompute) ------------------
def _lamina_sigma_np(eps: np.ndarray) -> np.ndarray:
    """Host mirror of envelope._lamina_sigma: soft below knee, stiff above, 0 past rupture."""
    sigma = np.where(eps <= KNEE, K_SOFT * eps, K_SOFT * KNEE + K_AC * (eps - KNEE))
    return np.where(eps > EPS_RUPT, 0.0, sigma)


def _lamina_force_np(pos: np.ndarray, faces: np.ndarray, a0: np.ndarray) -> np.ndarray:
    """Host mirror of lamina_areal_tension_kernel: per-face area-gradient tension f=-sigma*0.5*(nhat x e)."""
    f = np.zeros_like(pos)
    p0 = pos[faces[:, 0]]
    p1 = pos[faces[:, 1]]
    p2 = pos[faces[:, 2]]
    nrm = np.cross(p1 - p0, p2 - p0)
    a2 = np.linalg.norm(nrm, axis=1)
    area = 0.5 * a2
    eps = (area - a0) / a0
    sigma = _lamina_sigma_np(eps)
    nhat = nrm / a2[:, None]
    c = (-sigma * 0.5)[:, None]
    np.add.at(f, faces[:, 0], c * np.cross(nhat, p2 - p1))
    np.add.at(f, faces[:, 1], c * np.cross(nhat, p0 - p2))
    np.add.at(f, faces[:, 2], c * np.cross(nhat, p1 - p0))
    return f


def _chromatin_force_np(pos: np.ndarray, la: np.ndarray, lb: np.ndarray) -> np.ndarray:
    """Host mirror of chromatin_wlc_kernel (clamp r at 0.999999 exactly as the kernel)."""
    f = np.zeros_like(pos)
    d = pos[lb] - pos[la]
    x = np.linalg.norm(d, axis=1)
    r = np.minimum(x / CONTOUR_L, 0.999999)
    u = d / x[:, None]
    fmag = (KBT / PERSIST) * (r + 1.0 / (4.0 * (1.0 - r) ** 2) - 0.25)
    fv = fmag[:, None] * u
    np.add.at(f, la, fv)
    np.add.at(f, lb, -fv)
    return f


def _linc_force_np(npos: np.ndarray, apos: np.ndarray, n_idx: np.ndarray, a_idx: np.ndarray,
                   rest: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Host mirror of linc_tether_kernel: tension-only nonlinear tether f=k*delta*(1+stiff*delta^2)."""
    nf = np.zeros_like(npos)
    af = np.zeros_like(apos)
    d = apos[a_idx] - npos[n_idx]
    L = np.linalg.norm(d, axis=1)
    delta = L - rest
    active = delta > 0.0
    u = d / L[:, None]
    fmag = K_LINC * delta * (1.0 + STIFFENING * delta * delta)
    fv = np.where(active[:, None], fmag[:, None] * u, 0.0)
    np.add.at(nf, n_idx, fv)
    np.add.at(af, a_idx, -fv)
    return nf, af


def main() -> int:
    wp.init()
    dev = wp.get_device(DEVICE)
    if not dev.is_cuda:
        raise RuntimeError(f"NG-0 nucleus parity requires CUDA; CPU execution is forbidden (resolved {dev})")
    device = str(dev)

    # --- build the resting oblate mesh, then a fixed-seed perturbed config both paths share ---
    verts0, faces = build_oblate_mesh(R_EQ, ASPECT, SUBDIV)
    hinges = build_hinges(faces)                                   # (Nh,4) int64
    a0 = face_reference_areas(verts0, faces)                       # resting per-face areas
    Nv, Nf, Nh = verts0.shape[0], faces.shape[0], hinges.shape[0]
    kappa_tilde = calibrate_kappa_tilde(KAPPA_NE, verts0, hinges)

    rng = np.random.default_rng(SEED)
    pos_np = np.ascontiguousarray(verts0 + PERTURB * rng.standard_normal((Nv, 3)))
    vel_np = np.ascontiguousarray(rng.standard_normal((Nv, 3)) * 0.3)      # node velocities for viscosity
    apos_np = np.ascontiguousarray(pos_np + 0.6 * rng.standard_normal((Nv, 3)))  # LINC anchor cloud
    n_link = 200
    la = rng.integers(0, Nv, n_link).astype(np.int32)
    lb = rng.integers(0, Nv, n_link).astype(np.int32)
    keep = la != lb
    la, lb = la[keep], lb[keep]                                    # chromatin links (distinct endpoints)
    n_idx = rng.integers(0, Nv, 120).astype(np.int32)
    a_idx = rng.integers(0, Nv, 120).astype(np.int32)
    rest = (rng.random(120) * 1.2).astype(np.float64)              # mix of tension / slack (tension-only)

    print(f"NG-0 I2 nucleus parity gate | device={dev} | warp={wp.config.version} | "
          f"mesh: Nv={Nv} Nf={Nf} Nh={Nh} | chromatin_links={la.shape[0]} linc_links={n_idx.shape[0]} | "
          f"seed={SEED} rtol={RTOL:.0e}")
    print(f"    params(test): kappa_tilde={kappa_tilde:.4e} k_soft={K_SOFT} k_ac={K_AC} knee={KNEE} "
          f"eps_rupt={EPS_RUPT} p_vol={P_VOL} L={CONTOUR_L} lp={PERSIST} gamma={GAMMA_NUC} k_linc={K_LINC}")

    # --- device arrays (shared) ---
    pos = wp.array(pos_np, dtype=wp.vec3d, device=device)
    vel = wp.array(vel_np, dtype=wp.vec3d, device=device)
    apos = wp.array(apos_np, dtype=wp.vec3d, device=device)
    faces_wp = wp.array(np.ascontiguousarray(faces, dtype=np.int32), dtype=wp.int32, device=device)
    hinges_wp = wp.array(np.ascontiguousarray(hinges, dtype=np.int32), dtype=wp.int32, device=device)
    a0_wp = wp.array(np.ascontiguousarray(a0), dtype=wp.float64, device=device)
    ruptured_wp = wp.zeros(Nf, dtype=wp.int32, device=device)
    la_wp = wp.array(la, dtype=wp.int32, device=device)
    lb_wp = wp.array(lb, dtype=wp.int32, device=device)
    n_idx_wp = wp.array(n_idx, dtype=wp.int32, device=device)
    a_idx_wp = wp.array(a_idx, dtype=wp.int32, device=device)
    rest_wp = wp.array(rest, dtype=wp.float64, device=device)

    results: list[bool] = []

    # NG-0.n1 Helfrich bending (REUSED ff kernel) vs the exact numpy formula --------------------------
    force = wp.zeros(Nv, dtype=wp.vec3d, device=device)
    wp.launch(helfrich_bending_kernel, dim=Nh, inputs=[pos, hinges_wp, wp.float64(kappa_tilde)],
              outputs=[force], device=device)
    wp.synchronize()
    ref = dihedral_bending_forces(pos_np, hinges, kappa_tilde)
    results.append(_report("helfrich_bending_kernel", *_residual(force.numpy(), ref)))

    # NG-0.n2 lamina areal tension --------------------------------------------------------------------
    nforce = wp.zeros(Nv, dtype=wp.vec3d, device=device)
    wp.launch(lamina_areal_tension_kernel, dim=Nf,
              inputs=[pos, faces_wp, a0_wp, ruptured_wp, wp.float64(K_SOFT), wp.float64(K_AC),
                      wp.float64(KNEE), wp.float64(EPS_RUPT)], outputs=[nforce], device=device)
    wp.synchronize()
    ref = _lamina_force_np(pos_np, faces, a0)
    results.append(_report("lamina_areal_tension_kernel", *_residual(nforce.numpy(), ref)))

    # NG-0.n3 nucleoplasm volume force ----------------------------------------------------------------
    nforce = wp.zeros(Nv, dtype=wp.vec3d, device=device)
    wp.launch(nucleoplasm_volume_force_kernel, dim=Nf, inputs=[pos, faces_wp, wp.float64(P_VOL)],
              outputs=[nforce], device=device)
    wp.synchronize()
    ref = P_VOL * mesh_volume_gradient(pos_np, faces)
    results.append(_report("nucleoplasm_volume_force_kernel", *_residual(nforce.numpy(), ref)))

    # NG-0.n4 nucleoplasm volume reduce (scalar) ------------------------------------------------------
    vol_out = wp.zeros(1, dtype=wp.float64, device=device)
    wp.launch(nucleoplasm_volume_reduce_kernel, dim=Nf, inputs=[pos, faces_wp], outputs=[vol_out],
              device=device)
    wp.synchronize()
    ref_v = mesh_volume(pos_np, faces)
    results.append(_report("nucleoplasm_volume_reduce_kernel",
                           *_residual(np.array([vol_out.numpy()[0]]), np.array([ref_v]))))

    # NG-0.n5 chromatin WLC net -----------------------------------------------------------------------
    force = wp.zeros(Nv, dtype=wp.vec3d, device=device)
    wp.launch(chromatin_wlc_kernel, dim=la.shape[0],
              inputs=[pos, la_wp, lb_wp, wp.float64(CONTOUR_L), wp.float64(PERSIST), wp.float64(KBT)],
              outputs=[force], device=device)
    wp.synchronize()
    ref = _chromatin_force_np(pos_np, la, lb)
    results.append(_report("chromatin_wlc_kernel", *_residual(force.numpy(), ref)))

    # NG-0.n6 nucleoplasm viscosity -------------------------------------------------------------------
    nforce = wp.zeros(Nv, dtype=wp.vec3d, device=device)
    wp.launch(nucleoplasm_viscosity_kernel, dim=Nv, inputs=[vel, wp.float64(GAMMA_NUC)],
              outputs=[nforce], device=device)
    wp.synchronize()
    ref = -GAMMA_NUC * vel_np
    results.append(_report("nucleoplasm_viscosity_kernel", *_residual(nforce.numpy(), ref)))

    # NG-0.n7 LINC tether (two force outputs) ---------------------------------------------------------
    nforce = wp.zeros(Nv, dtype=wp.vec3d, device=device)
    aforce = wp.zeros(Nv, dtype=wp.vec3d, device=device)
    wp.launch(linc_tether_kernel, dim=n_idx.shape[0],
              inputs=[pos, apos, n_idx_wp, a_idx_wp, rest_wp, wp.float64(K_LINC), wp.float64(STIFFENING)],
              outputs=[nforce, aforce], device=device)
    wp.synchronize()
    ref_n, ref_a = _linc_force_np(pos_np, apos_np, n_idx, a_idx, rest)
    max_n, rel_n = _residual(nforce.numpy(), ref_n)
    max_a, rel_a = _residual(aforce.numpy(), ref_a)
    results.append(_report("linc_tether_kernel", max(max_n, max_a), max(rel_n, rel_a)))

    ok = all(results)
    print(f"NG-0 NUCLEUS OVERALL: {'PASS' if ok else 'FAIL'} "
          f"({sum(results)}/{len(results)} kernels round-off-faithful)")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
