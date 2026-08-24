"""NG-0 — Warp-CUDA <-> NumPy source-faithfulness parity gate (run on the A5000 by the lead).

For each of the three fluid-spine kernels, build ONE small fixed-seed case (16^3, random field, no-flux
box), advance ONE step on BOTH the Warp-CUDA device path and the pure-NumPy acceptance oracle from the
IDENTICAL initial state + identical parameters/masks/BCs, and report the absolute + relative parity
residual:

    NG-0.a  biot_pmass_update_kernel  vs  fv_reference.BiotFVReference.step
    NG-0.b  rad_transport_kernel       vs  transport_reference.RADTransportReference.step
    NG-0.c  darcy_discharge_kernel     vs  darcy_analytic.darcy_flux (applied to the kernel's own stencil grad)

PASS iff the relative residual max|warp - ref| / max|ref| < 1e-10 (round-off). This certifies the device
kernel IS the CPU-verified stencil. STRUCTURAL / parameter-agnostic: the I0-B1 closed params (c_v, mobility,
S) are used for realism, but any positive test value would give the same verdict (this is a round-off check,
not a magnitude verdict). phi / c_0 are test values (they do not enter the stepped stencils).

I0-A: constructing the FieldGrid / BiotSubstrate / DarcyVelocity / MonomerField REQUIRES a CUDA device, so
this runner only executes on the gbook A5000 (never the dev Mac). Do NOT tune, loosen the tolerance, or
"fix" a kernel to force a pass — a genuine NG-0 failure is a real on-device bug the CPU oracle could not
catch, and must be reported honestly.

Run (from ~/ffn_ac_native on gbook):
    PYTHONPATH=. ~/miniconda3/envs/ffn_sim/bin/python aleph/components/fluid/native_gates/ng0_parity.py
"""

from __future__ import annotations

import sys

import numpy as np
import warp as wp

from aleph.components.fluid import darcy_analytic
from aleph.components.fluid.biot_substrate import BiotSubstrate
from aleph.components.fluid.field_grid import FieldGrid
from aleph.components.fluid.fv_reference import BiotFVReference
from aleph.components.fluid.transport import MonomerField
from aleph.components.fluid.transport_reference import RADTransportReference
from aleph.components.fluid.velocity import DarcyVelocity

# ----------------------------------------------------------------------------------------------------
# Fixed test case (deterministic; STRUCTURAL — magnitudes irrelevant to a round-off parity verdict).
# ----------------------------------------------------------------------------------------------------
SHAPE = (16, 16, 16)
DX = 1.0            # um
DEVICE: str | None = None  # resolve the current CUDA device; the hardware contract forbids fixed ordinals
SEED = 20260716
RTOL = 1.0e-10      # round-off pass threshold (relative)

# I0-B1 closed params (used for realism; a round-off parity verdict is param-agnostic).
STORAGE_S = 1.0e-3          # 1/Pa   (S = 1/M)
C_V = 50.0                  # um^2/s (spec NG-4 anchor)
MOBILITY = C_V * STORAGE_S  # um^2/(Pa*s)   (c_v = mobility / S)
ALPHA = 1.0                 # Biot-Willis (I0-B1 ratified)
D_C = 4.0                   # um^2/s (I0-B1c draft; does not change a round-off verdict)
PHI = 0.7                   # test value; does not enter the stepped RAD stencil
RHO_F_B = (1.3, -2.1, 0.7)  # Pa/um; exercises the Darcy body-force term (both paths get it identically)


def _residual(warp_arr: np.ndarray, ref_arr: np.ndarray) -> tuple[float, float]:
    """Return (max_abs, rel) = (max|warp-ref|, max|warp-ref| / max|ref|) with a zero-ref guard."""
    warp_arr = np.asarray(warp_arr, dtype=np.float64)
    ref_arr = np.asarray(ref_arr, dtype=np.float64)
    max_abs = float(np.max(np.abs(warp_arr - ref_arr)))
    scale = float(np.max(np.abs(ref_arr)))
    rel = max_abs / scale if scale > 0.0 else max_abs
    return max_abs, rel


def _report(name: str, max_abs: float, rel: float) -> bool:
    passed = rel < RTOL
    verdict = "PASS" if passed else "FAIL"
    print(f"NG-0 {name}: max_abs={max_abs:.3e} rel={rel:.3e} -> {verdict}")
    return passed


# ----------------------------------------------------------------------------------------------------
# NG-0.a  Biot p/mass update  (biot_pmass_update_kernel vs BiotFVReference.step)
# ----------------------------------------------------------------------------------------------------
def gate_biot_pmass(rng: np.random.Generator) -> bool:
    """One explicit conservative p/mass step: device kernel vs the FV NumPy reference, same input."""
    p0 = rng.standard_normal(SHAPE).astype(np.float64)          # random p on a no-flux box
    s_water = rng.standard_normal(SHAPE).astype(np.float64) * 1e-2
    div_vs = rng.standard_normal(SHAPE).astype(np.float64) * 1e-3
    dt = 0.5 * (STORAGE_S * DX * DX / (2.0 * 3 * MOBILITY))     # 0.5 * CFL (identical for both paths)

    # --- Warp-CUDA path: mask is all-FLUID by default (a no-flux box) ---
    grid = FieldGrid(SHAPE, DX, device=DEVICE)
    with wp.ScopedDevice(grid.device):
        grid.p = wp.array(np.ascontiguousarray(p0), dtype=wp.float64)
        grid.s_water = wp.array(np.ascontiguousarray(s_water), dtype=wp.float64)
        grid.div_vs = wp.array(np.ascontiguousarray(div_vs), dtype=wp.float64)
    sub = BiotSubstrate(grid, mobility=MOBILITY, storage_S=STORAGE_S, alpha=ALPHA)
    sub.step(dt)                                                # ping-pong -> grid.p is the new field
    p_warp = grid.p.numpy()

    # --- NumPy reference: identical grid / p0 / params / all-fluid mask / no dirichlet ---
    ref = BiotFVReference(SHAPE, DX, MOBILITY, STORAGE_S)       # mask=None -> all FLUID
    p_ref = ref.step(p0, dt, solid_dilatation_rate=div_vs, s_water=s_water, alpha=ALPHA)

    max_abs, rel = _residual(p_warp, p_ref)
    return _report("biot_pmass_update", max_abs, rel)


# ----------------------------------------------------------------------------------------------------
# NG-0.b  RAD monomer transport  (rad_transport_kernel vs RADTransportReference.step)
# ----------------------------------------------------------------------------------------------------
def gate_rad_transport(rng: np.random.Generator) -> bool:
    """One explicit conservative RAD step: device kernel vs the RAD NumPy reference, same input."""
    u0 = (rng.random(SHAPE) + 0.5).astype(np.float64)                     # positive density u = phi c
    reaction = rng.standard_normal(SHAPE).astype(np.float64) * 1e-2
    v_f = (rng.standard_normal((*SHAPE, 3)) * 0.5).astype(np.float64)     # cell-centred pore velocity (um/s)
    v_max = float(np.max(np.abs(v_f)))
    dt = 0.5 / (2.0 * 3 * D_C / (DX * DX) + v_max / DX)                   # 0.5 * combined CFL (both paths)

    # --- Warp-CUDA path ---
    grid = FieldGrid(SHAPE, DX, device=DEVICE)                            # all-FLUID mask (no-flux box)
    field = MonomerField(grid, phi=PHI, d_c=D_C, c0=None)
    with wp.ScopedDevice(grid.device):
        field.u = wp.array(np.ascontiguousarray(u0), dtype=wp.float64)
        field.reaction = wp.array(np.ascontiguousarray(reaction), dtype=wp.float64)
        v_f_wp = wp.array(np.ascontiguousarray(v_f), dtype=wp.vec3d)
    field.step(dt, v_f_wp)                                                # ping-pong -> field.u is the new field
    u_warp = field.u.numpy()

    # --- NumPy reference: identical grid / u0 / v_f / reaction / all-fluid mask ---
    ref = RADTransportReference(SHAPE, DX, PHI, D_C)                      # mask=None -> all FLUID
    u_ref = ref.step(u0, dt, v_f=v_f, reaction=reaction)

    max_abs, rel = _residual(u_warp, u_ref)
    return _report("rad_transport", max_abs, rel)


# ----------------------------------------------------------------------------------------------------
# NG-0.c  Darcy discharge  (darcy_discharge_kernel vs darcy_analytic.darcy_flux)
# ----------------------------------------------------------------------------------------------------
def _kernel_stencil_grad(p: np.ndarray, dx: float) -> np.ndarray:
    """Replicate darcy_discharge_kernel's masked central-difference grad on an ALL-FLUID box (host).

    ``darcy_flux`` is a pure constitutive law q = -mobility (grad p - rho_f b); it does NOT compute a grad
    from a field. To compare the kernel (which bundles grad + constitutive) against it apples-to-apples we
    reproduce the kernel's EXACT stencil here and feed that grad through darcy_flux. On an all-fluid box:
      interior -> central (p[i+1]-p[i-1])/(2 dx);
      i=0 edge -> xm falls back to pc  =>  (p[1]-p[0])/(2 dx)   (the kernel's one-sided-against-centre form);
      i=nx-1   -> xp falls back to pc  =>  (p[nx-1]-p[nx-2])/(2 dx).
    """
    inv2dx = 1.0 / (2.0 * dx)
    grad = np.zeros((*p.shape, 3), dtype=np.float64)
    for axis in range(3):
        xm = p.copy()   # default neighbour = centre (kernel's no-neighbour fallback)
        xp = p.copy()
        lo = [slice(None)] * 3
        hi = [slice(None)] * 3
        lo[axis] = slice(1, None)     # cells with a low neighbour
        hi[axis] = slice(0, -1)       # cells with a high neighbour
        src_lo = [slice(None)] * 3
        src_hi = [slice(None)] * 3
        src_lo[axis] = slice(0, -1)   # the low neighbour p[i-1]
        src_hi[axis] = slice(1, None) # the high neighbour p[i+1]
        xm[tuple(lo)] = p[tuple(src_lo)]   # all-fluid: neighbour is always fluid
        xp[tuple(hi)] = p[tuple(src_hi)]
        grad[..., axis] = (xp - xm) * inv2dx
    return grad


def gate_darcy_discharge(rng: np.random.Generator) -> bool:
    """Cell-centred Darcy discharge: device kernel vs darcy_flux on the kernel's own stencil grad."""
    p0 = rng.standard_normal(SHAPE).astype(np.float64)

    # --- Warp-CUDA path ---
    grid = FieldGrid(SHAPE, DX, device=DEVICE)                # all-FLUID mask (no-flux box)
    with wp.ScopedDevice(grid.device):
        grid.p = wp.array(np.ascontiguousarray(p0), dtype=wp.float64)
    dv = DarcyVelocity(grid, mobility=MOBILITY, phi=PHI, rho_f_b=RHO_F_B)
    dv.discharge()
    q_warp = dv.q.numpy()                                     # (nx, ny, nz, 3)

    # --- NumPy reference: replicate the kernel's stencil grad, then apply the darcy_flux constitutive ---
    grad = _kernel_stencil_grad(p0, DX)
    q_ref = darcy_analytic.darcy_flux(grad, MOBILITY, np.asarray(RHO_F_B, dtype=np.float64))

    max_abs, rel = _residual(q_warp, q_ref)
    return _report("darcy_discharge", max_abs, rel)


def main() -> int:
    wp.init()
    dev = wp.get_device(DEVICE)
    print(f"NG-0 fluid-spine parity gate | device={dev} | warp={wp.config.version} | "
          f"case={SHAPE} dx={DX} seed={SEED} rtol={RTOL:.0e}")
    print(f"    params: mobility={MOBILITY:.4e} S={STORAGE_S:.4e} c_v={C_V} alpha={ALPHA} "
          f"D_c={D_C} phi={PHI} rho_f_b={RHO_F_B}")

    rng = np.random.default_rng(SEED)
    results = [
        gate_biot_pmass(rng),
        gate_rad_transport(rng),
        gate_darcy_discharge(rng),
    ]
    ok = all(results)
    print(f"NG-0 OVERALL: {'PASS' if ok else 'FAIL'} ({sum(results)}/{len(results)} kernels round-off-faithful)")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
