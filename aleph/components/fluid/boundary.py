"""Moving-domain boundary conditions (Warp-CUDA kernel SOURCE) — I1a.

Two boundaries move with the live cell:

* **Outer membrane — hydraulic-flux (Neumann) BC.** Water crosses the bilayer by the Kedem-Katchalsky
  law: the volume flux per unit membrane area is ``s = L_p (sigma dPi_osm - dP)``, ``dP = p_local - p_ext``
  (Starling/osmotic form; reflection coeff ``sigma`` defaults 1). Written into the dedicated
  ``s_membrane`` source channel by quadrature over the **live membrane triangles**.  Each triangle evaluates
  the affine interior pressure trace at its surface centroid and deposits ``s*A`` to a normalized FLUID
  stencil, divided by cell volume.  Therefore ``sum(s_membrane*dV) == sum_triangle(s*A)`` while using the true
  rotationally invariant triangle area, rather than the voxel/L1 area of FLUID-OUTSIDE Cartesian faces. This
  REPLACES the
  0-D lumped ``tau_osm`` whole-cell reservoir (``network_warp.py`` scalar turgor) with a spatial surface BC.

* **Nuclear envelope — relative no-flux BC.** No Darcy discharge crosses the nucleus surface RELATIVE TO
  THE SOLID (the envelope moves with the skeleton; pore fluid does not seep into the nucleoplasm on the
  poroelastic timescale). Enforced structurally by the mask: NUCLEUS faces carry zero flux in
  ``biot_pmass_update_kernel``. This class provides the semantics + a device diagnostic that certifies
  zero net flux across the nucleus surface (the native-gate readback).

Runtime CUDA-only (I0-A). Host acceptance: ``fv_reference`` drives the same membrane flux through
``darcy_divergence(membrane_flux=...)`` / ``boundary_flux_integral`` (the content==flux gate).
"""

from __future__ import annotations

import warp as wp  # noqa: E402  (Warp-CUDA runtime; HOOMD never imported)

from aleph.components.fluid.field_grid import FieldGrid
from aleph.components.fluid.surface_trace import membrane_surface_pressure_trace, peskin4

__all__ = [
    "live_triangle_membrane_flux_kernel",
    "membrane_flux_source_kernel",
    "MembraneFluxBC",
    "NucleusNoFluxBC",
]

_FLUID = wp.constant(1)
_OUTSIDE = wp.constant(0)
_NUCLEUS = wp.constant(2)


@wp.kernel
def live_triangle_membrane_flux_kernel(
    pos: wp.array(dtype=wp.vec3d),
    faces: wp.array(dtype=wp.int32, ndim=2),
    p: wp.array3d(dtype=wp.float64),
    mask: wp.array3d(dtype=wp.int32),
    origin: wp.vec3d,
    dx: wp.float64,
    L_p: wp.float64,
    d_pi_osm: wp.float64,
    p_ext: wp.float64,
    sigma_refl: wp.float64,
    unresolved_faces: wp.array(dtype=wp.int32),
    integrated_flux: wp.array(dtype=wp.float64),
    s_membrane: wp.array3d(dtype=wp.float64),
) -> None:
    """Deposit live-triangle Kedem-Katchalsky volume flux conservatively onto FLUID cells.

    For triangle ``t``, ``Q_t = L_p*(sigma*dPi-(p_inside-p_ext))*A_t`` [um^3/s].  A mask-aware affine trace
    obtains ``p_inside`` and normalized four-point Peskin weights distribute ``Q_t/dV`` to FLUID cells.
    Consequently the discrete source integral equals the true-area surface quadrature to reduction round-off.
    """
    t = wp.tid()
    i0 = faces[t, 0]
    i1 = faces[t, 1]
    i2 = faces[t, 2]
    x0 = pos[i0]
    x1 = pos[i1]
    x2 = pos[i2]
    area2_vec = wp.cross(x1 - x0, x2 - x0)
    area2 = wp.length(area2_vec)
    if area2 <= wp.float64(0.0):
        wp.atomic_add(unresolved_faces, 0, 1)
        return

    n_out = area2_vec / area2
    centroid = (x0 + x1 + x2) / wp.float64(3.0)
    xq = centroid - wp.float64(0.5) * dx * n_out
    gx = (xq[0] - origin[0]) / dx
    gy = (xq[1] - origin[1]) / dx
    gz = (xq[2] - origin[2]) / dx
    bi = int(wp.floor(gx))
    bj = int(wp.floor(gy))
    bk = int(wp.floor(gz))
    nx = mask.shape[0]
    ny = mask.shape[1]
    nz = mask.shape[2]

    trace = membrane_surface_pressure_trace(p, mask, origin, dx, centroid, n_out)
    if trace[1] == wp.float64(0.0):
        wp.atomic_add(unresolved_faces, 0, 1)
        return
    w_sum = wp.float64(0.0)
    for di in range(-1, 3):
        ii = bi + di
        if ii < 0 or ii >= nx:
            continue
        rx = (origin[0] + wp.float64(ii) * dx - xq[0]) / dx
        wx = peskin4(rx)
        for dj in range(-1, 3):
            jj = bj + dj
            if jj < 0 or jj >= ny:
                continue
            ry = (origin[1] + wp.float64(jj) * dx - xq[1]) / dx
            wy = peskin4(ry)
            for dk in range(-1, 3):
                kk = bk + dk
                if kk < 0 or kk >= nz or mask[ii, jj, kk] != _FLUID:
                    continue
                rz = (origin[2] + wp.float64(kk) * dx - xq[2]) / dx
                wz = peskin4(rz)
                w = wx * wy * wz
                w_sum += w

    if w_sum <= wp.float64(0.0):
        wp.atomic_add(unresolved_faces, 0, 1)
        return
    p_inside = trace[0]
    influx = L_p * (sigma_refl * d_pi_osm - (p_inside - p_ext))
    q_face = influx * wp.float64(0.5) * area2
    wp.atomic_add(integrated_flux, 0, q_face)
    inv_cell_volume = wp.float64(1.0) / (dx * dx * dx)
    for di in range(-1, 3):
        ii = bi + di
        if ii < 0 or ii >= nx:
            continue
        rx = (origin[0] + wp.float64(ii) * dx - xq[0]) / dx
        wx = peskin4(rx)
        for dj in range(-1, 3):
            jj = bj + dj
            if jj < 0 or jj >= ny:
                continue
            ry = (origin[1] + wp.float64(jj) * dx - xq[1]) / dx
            wy = peskin4(ry)
            for dk in range(-1, 3):
                kk = bk + dk
                if kk < 0 or kk >= nz or mask[ii, jj, kk] != _FLUID:
                    continue
                rz = (origin[2] + wp.float64(kk) * dx - xq[2]) / dx
                wz = peskin4(rz)
                w = wx * wy * wz
                wp.atomic_add(s_membrane, ii, jj, kk, q_face * w / w_sum * inv_cell_volume)


@wp.kernel
def membrane_flux_source_kernel(
    p: wp.array3d(dtype=wp.float64),
    mask: wp.array3d(dtype=wp.int32),
    L_p: wp.float64,
    d_pi_osm: wp.float64,
    p_ext: wp.float64,
    sigma_refl: wp.float64,
    dx: wp.float64,
    s_membrane: wp.array3d(dtype=wp.float64),
) -> None:
    """Write membrane hydraulic influx into its dedicated source channel.

    ``s_membrane[cell] += (#OUTSIDE faces) * L_p (sigma dPi_osm - (p - p_ext)) / dx`` (inward positive).
    The caller clears this channel before launch; the independent interior ``s_water`` is untouched.
    """
    i, j, k = wp.tid()
    nx = mask.shape[0]
    ny = mask.shape[1]
    nz = mask.shape[2]
    if mask[i, j, k] != _FLUID:
        return
    n_faces = 0
    if i > 0 and mask[i - 1, j, k] == _OUTSIDE:
        n_faces += 1
    if i < nx - 1 and mask[i + 1, j, k] == _OUTSIDE:
        n_faces += 1
    if j > 0 and mask[i, j - 1, k] == _OUTSIDE:
        n_faces += 1
    if j < ny - 1 and mask[i, j + 1, k] == _OUTSIDE:
        n_faces += 1
    if k > 0 and mask[i, j, k - 1] == _OUTSIDE:
        n_faces += 1
    if k < nz - 1 and mask[i, j, k + 1] == _OUTSIDE:
        n_faces += 1
    if n_faces == 0:
        return
    d_p = p[i, j, k] - p_ext
    influx = L_p * (sigma_refl * d_pi_osm - d_p)  # inward volume flux per unit membrane area
    s_membrane[i, j, k] += wp.float64(n_faces) * influx / dx


@wp.kernel
def _nucleus_suppressed_flux_kernel(
    p: wp.array3d(dtype=wp.float64),
    mask: wp.array3d(dtype=wp.int32),
    mobility: wp.float64,
    dx: wp.float64,
    out_abs_flux: wp.array(dtype=wp.float64),
) -> None:
    """Accumulate the WOULD-BE Darcy flux ``|mobility·(p_fluid − p_nucleus)/dx|`` across every FLUID↔NUCLEUS
    face — the discharge the relative-no-flux mask holds at zero.

    ⚠ Non-vacuity fix (I2 no-flux reconciliation). The previous probe computed ``mobility·(p−p)/dx`` — a
    self-minus-self difference that is IDENTICALLY 0 for ANY field, so it certified nothing (it could never
    detect a leak). This reads the REAL neighbour cell ``p[nucleus]`` (a DIFFERENT cell) so the value is a
    genuine measurement: 0 only when the pore-pressure field is uniform ACROSS the envelope (nothing to
    suppress), and > 0 under any transmembrane gradient — the magnitude the no-flux BC is actively holding at
    zero. Since ``biot_pmass_update_kernel`` skips non-FLUID neighbours, the APPLIED cross-nucleus flux stays 0
    (verified separately by the global content balance, NG-2); this probe proves that suppression is LOAD-BEARING
    (has teeth) rather than a uniform-field accident, so NG-5 is falsifiable.
    """
    i, j, k = wp.tid()
    nx = mask.shape[0]
    ny = mask.shape[1]
    nz = mask.shape[2]
    if mask[i, j, k] != _FLUID:
        return
    pc = p[i, j, k]
    flux = wp.float64(0.0)
    if i > 0 and mask[i - 1, j, k] == _NUCLEUS:
        flux += wp.abs(mobility * (pc - p[i - 1, j, k]) / dx)
    if i < nx - 1 and mask[i + 1, j, k] == _NUCLEUS:
        flux += wp.abs(mobility * (pc - p[i + 1, j, k]) / dx)
    if j > 0 and mask[i, j - 1, k] == _NUCLEUS:
        flux += wp.abs(mobility * (pc - p[i, j - 1, k]) / dx)
    if j < ny - 1 and mask[i, j + 1, k] == _NUCLEUS:
        flux += wp.abs(mobility * (pc - p[i, j + 1, k]) / dx)
    if k > 0 and mask[i, j, k - 1] == _NUCLEUS:
        flux += wp.abs(mobility * (pc - p[i, j, k - 1]) / dx)
    if k < nz - 1 and mask[i, j, k + 1] == _NUCLEUS:
        flux += wp.abs(mobility * (pc - p[i, j, k + 1]) / dx)
    wp.atomic_add(out_abs_flux, 0, flux)


class MembraneFluxBC:
    """Outer-membrane Kedem-Katchalsky water-flux BC folded into the field source.

    Args:
        grid: The device field grid.
        L_p: Membrane hydraulic conductivity [um/(s*Pa)] (I0-B1 draft, MCF7/AQP5).
        p_ext: External hydrostatic pressure [Pa] (0 baseline).
        sigma_refl: Reflection coefficient (1.0 = ideal semipermeable).
        pos: Authoritative global live-node positions.  When supplied with ``faces``, production uses true-area
            triangle quadrature.  Omitting both retains the voxel-face operator for isolated diagnostics only.
        faces: Global membrane triangle indices.
    """

    def __init__(
        self,
        grid: FieldGrid,
        *,
        L_p: float,
        p_ext: float = 0.0,
        sigma_refl: float = 1.0,
        pos: wp.array | None = None,
        faces: wp.array | None = None,
        n_faces: int = 0,
    ) -> None:
        self.grid = grid
        self.L_p = float(L_p)
        self.p_ext = float(p_ext)
        self.sigma_refl = float(sigma_refl)
        if (pos is None) != (faces is None):
            raise ValueError("pos and faces must be supplied together for live-triangle membrane flux")
        self.pos = pos
        self.faces = faces
        self.n_faces = int(n_faces)
        if self.faces is not None and self.n_faces <= 0:
            raise ValueError("n_faces must be positive for live-triangle membrane flux")
        self.discretization = "LIVE_TRIANGLE_CONSERVATIVE" if self.faces is not None else "VOXEL_DIAGNOSTIC_ONLY"
        with wp.ScopedDevice(grid.device):
            self.unresolved_faces_d = wp.zeros(1, dtype=wp.int32, device=grid.device)
            self.integrated_flux_d = wp.zeros(1, dtype=wp.float64, device=grid.device)

    def apply(self, d_pi_osm: float) -> None:
        """Rebuild the membrane-only source for the current osmotic difference [Pa]."""
        g = self.grid
        with wp.ScopedDevice(g.device):
            g.s_membrane.zero_()
            self.integrated_flux_d.zero_()
            if self.faces is not None and self.pos is not None:
                wp.launch(
                    live_triangle_membrane_flux_kernel,
                    dim=self.n_faces,
                    inputs=[
                        self.pos, self.faces, g.p, g.mask, wp.vec3d(*g.origin.tolist()), wp.float64(g.dx),
                        wp.float64(self.L_p), wp.float64(d_pi_osm), wp.float64(self.p_ext),
                        wp.float64(self.sigma_refl), self.unresolved_faces_d, self.integrated_flux_d,
                    ],
                    outputs=[g.s_membrane],
                    device=g.device,
                )
            else:
                wp.launch(
                    membrane_flux_source_kernel,
                    dim=g.shape,
                    inputs=[g.p, g.mask, wp.float64(self.L_p), wp.float64(d_pi_osm),
                            wp.float64(self.p_ext), wp.float64(self.sigma_refl), wp.float64(g.dx)],
                    outputs=[g.s_membrane],
                    device=g.device,
                )

    def reset_diagnostics(self) -> None:
        """Start a new explicit gate/run window; normal subcycles preserve earlier unresolved-face failures."""
        self.unresolved_faces_d.zero_()
        self.integrated_flux_d.zero_()


class NucleusNoFluxBC:
    """Relative no-flux at the (moving) nuclear envelope — enforced by the mask, audited by a real probe.

    The update kernel already skips NUCLEUS neighbours, so no Darcy discharge crosses the nucleus surface in
    the SOLID frame; envelope motion is carried by the conservative moving-domain remap (``domain.Domain``),
    not by a flux. The audit has two halves at the native gate (NG-5):

    * :meth:`suppressed_flux` — the WOULD-BE cross-envelope discharge the mask holds at zero (reads the real
      nucleus-neighbour pressure, so it is > 0 under any transmembrane gradient). This is the falsifiable,
      non-vacuous replacement for the old ``p−p`` probe that could only ever read 0.
    * the global content balance (NG-2) — total fluid content changes ONLY by the membrane flux term, with no
      nucleus contribution, so the APPLIED cross-envelope flux is genuinely 0 even when :meth:`suppressed_flux`
      is large (the mask is load-bearing, not a uniform-field accident).
    """

    def __init__(self, grid: FieldGrid, *, mobility: float) -> None:
        self.grid = grid
        self.mobility = float(mobility)

    def suppressed_flux(self) -> float:
        """Total would-be |Darcy discharge| across nucleus faces that the relative-no-flux mask holds at zero.

        0 iff the pore-pressure field is uniform across the envelope (nothing to suppress); > 0 under any
        transmembrane gradient (the BC is actively enforcing the no-flux invariant). Out-of-hot-loop host read.
        """
        g = self.grid
        with wp.ScopedDevice(g.device):
            probe = wp.zeros(1, dtype=wp.float64)
            wp.launch(_nucleus_suppressed_flux_kernel, dim=g.shape,
                      inputs=[g.p, g.mask, wp.float64(self.mobility), wp.float64(g.dx)],
                      outputs=[probe])
            return float(probe.numpy()[0])

    # Back-compat alias: the applied cross-envelope flux is structurally 0 (mask-skip). The meaningful,
    # non-vacuous audit is suppressed_flux()>0 (teeth) + the content balance (no leak); see the class docstring.
    def leaked_flux(self) -> float:
        """Deprecated name for :meth:`suppressed_flux` (the probe was reframed to be non-vacuous)."""
        return self.suppressed_flux()
