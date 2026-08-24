r"""Mask-aware affine pressure trace from the Eulerian field to a live immersed surface (Warp-CUDA).

A renormalized trilinear stencil preserves constants but not linear fields when OUTSIDE cells are removed near
the membrane.  Surface traction and Kedem-Katchalsky flux require the **interior trace** of pressure, so this
module fits the local affine polynomial

``p(x_support + dx*r) = c0 + c1*r_x + c2*r_y + c3*r_z``

to FLUID cell centres in the four-point Peskin support.  The support may be centred at the first interior
finite-volume location for conditioning while the fitted polynomial is evaluated at the actual membrane
quadrature point.  The weighted normal equations remain local to one CUDA thread.  An insufficient or
rank-deficient stencil returns ``valid=0``; callers latch that as a device failure rather than inventing a
fallback pressure.

Sanity Gate:
    * Reproduction: constant and affine fields are exact to float64 round-off whenever the masked design matrix
      has rank four; this is the native spatial-linear-trace gate.
    * Dimensions: coordinates are normalized by ``dx``, so only the fitted intercept carries pressure units.
    * Boundary: the compact Peskin support reads FLUID cells only; OUTSIDE/NUCLEUS values cannot leak in.
    * Numerical: the determinant test is scaled by the support weight and ``64*eps`` (at most 4^3 samples),
      a derived round-off bound rather than an empirical acceptance knob.
    * Residency: all samples, normal equations, inversion and validity state remain in the calling CUDA kernel.
"""

from __future__ import annotations

import warp as wp

__all__ = ["masked_affine_pressure_trace", "membrane_surface_pressure_trace", "peskin4"]

_FLUID = wp.constant(1)


@wp.func
def peskin4(r: wp.float64) -> wp.float64:
    """Four-point Peskin regularized delta weight on a dimensionless grid coordinate."""
    a = wp.abs(r)
    if a <= wp.float64(1.0):
        return (
            wp.float64(3.0) - wp.float64(2.0) * a
            + wp.sqrt(wp.float64(1.0) + wp.float64(4.0) * a - wp.float64(4.0) * a * a)
        ) / wp.float64(8.0)
    if a <= wp.float64(2.0):
        inner = wp.max(
            wp.float64(0.0),
            wp.float64(-7.0) + wp.float64(12.0) * a - wp.float64(4.0) * a * a,
        )
        return (
            wp.float64(5.0) - wp.float64(2.0) * a - wp.sqrt(inner)
        ) / wp.float64(8.0)
    return wp.float64(0.0)


@wp.func
def masked_affine_pressure_trace(
    p: wp.array3d(dtype=wp.float64),
    mask: wp.array3d(dtype=wp.int32),
    origin: wp.vec3d,
    dx: wp.float64,
    support_center: wp.vec3d,
    evaluation_point: wp.vec3d,
) -> wp.vec2d:
    """Return ``(p(evaluation_point), valid)`` from a local FLUID affine fit."""
    gx = (support_center[0] - origin[0]) / dx
    gy = (support_center[1] - origin[1]) / dx
    gz = (support_center[2] - origin[2]) / dx
    bi = int(wp.floor(gx))
    bj = int(wp.floor(gy))
    bk = int(wp.floor(gz))
    nx = mask.shape[0]
    ny = mask.shape[1]
    nz = mask.shape[2]
    normal = wp.mat44d()
    rhs = wp.vec4d(0.0, 0.0, 0.0, 0.0)
    weight_sum = wp.float64(0.0)

    for di in range(-1, 3):
        ii = bi + di
        if ii < 0 or ii >= nx:
            continue
        rx = (origin[0] + wp.float64(ii) * dx - support_center[0]) / dx
        wx = peskin4(rx)
        for dj in range(-1, 3):
            jj = bj + dj
            if jj < 0 or jj >= ny:
                continue
            ry = (origin[1] + wp.float64(jj) * dx - support_center[1]) / dx
            wy = peskin4(ry)
            for dk in range(-1, 3):
                kk = bk + dk
                if kk < 0 or kk >= nz or mask[ii, jj, kk] != _FLUID:
                    continue
                rz = (origin[2] + wp.float64(kk) * dx - support_center[2]) / dx
                wz = peskin4(rz)
                weight = wx * wy * wz
                basis = wp.vec4d(wp.float64(1.0), rx, ry, rz)
                normal += weight * wp.outer(basis, basis)
                rhs += weight * basis * p[ii, jj, kk]
                weight_sum += weight

    det = wp.determinant(normal)
    scale = wp.max(wp.float64(1.0), weight_sum)
    scale2 = scale * scale
    determinant_floor = wp.float64(64.0) * wp.float64(2.220446049250313e-16) * scale2 * scale2
    if wp.abs(det) <= determinant_floor:
        return wp.vec2d(0.0, 0.0)
    coefficients = wp.mul(wp.inverse(normal), rhs)
    eval_basis = wp.vec4d(
        wp.float64(1.0),
        (evaluation_point[0] - support_center[0]) / dx,
        (evaluation_point[1] - support_center[1]) / dx,
        (evaluation_point[2] - support_center[2]) / dx,
    )
    return wp.vec2d(wp.dot(coefficients, eval_basis), 1.0)


@wp.func
def membrane_surface_pressure_trace(
    p: wp.array3d(dtype=wp.float64),
    mask: wp.array3d(dtype=wp.int32),
    origin: wp.vec3d,
    dx: wp.float64,
    centroid: wp.vec3d,
    n_out: wp.vec3d,
) -> wp.vec2d:
    """Return the one-sided second-order surface pressure trace and validity flag."""
    support_half = centroid - wp.float64(0.5) * dx * n_out
    support_full = centroid - dx * n_out
    trace_half = masked_affine_pressure_trace(p, mask, origin, dx, support_half, support_half)
    trace_full = masked_affine_pressure_trace(p, mask, origin, dx, support_full, support_full)
    if trace_half[1] == wp.float64(0.0) or trace_full[1] == wp.float64(0.0):
        return wp.vec2d(0.0, 0.0)
    return wp.vec2d(wp.float64(2.0) * trace_half[0] - trace_full[0], 1.0)
