"""D0 — membrane shell + Biot-Darcy cytosol traction rig (Lane D, NON-AUTHORITATIVE DIAGNOSTIC).

Minimal scene::

    deformable spherical membrane (own node set, lipid-sheet mechanics)
            ↕  live pressure traction / water-flux boundary
    Biot-Darcy cytosol pressure field  (FieldGrid + BiotSubstrate)
            ↕
    localized interior pressure / volume perturbation (one-sided patch)

Read-only primitives reused
    * ``aleph.laws.membrane_surface``       — membrane mesh, area-tension, Helfrich bending, κ̃ calibration
    * ``aleph.components.incumbent.membrane_pressure`` — pressure → membrane-normal nodal traction (``MembranePressureTraction``)
    * ``aleph.components.fluid.field_grid``       — ``FieldGrid`` shared Eulerian state
    * ``aleph.components.fluid.biot_substrate``   — ``BiotSubstrate`` conservative p/mass solver
    * ``aleph.components.fluid.boundary``         — ``MembraneFluxBC`` Kedem-Katchalsky water flux
    * ``aleph.components.incumbent.fsi_coupling``      — ``SolidDilatationCoupling`` membrane-motion → fluid moving-boundary

Lane-local glue (this file only): an overdamped quasi-static membrane relaxation
kernel and a localized interior fluid-source perturbation kernel.

Pressure bookkeeping (kept as THREE distinct scalars — never summed into one):
    * ``dP_hyd``   resting hydrostatic difference (interior−exterior). 40 Pa is the
      **HeLa diagnostic hydrostatic proxy** (Fischer-Friedrich), *not* an MCF7
      production value. Seeds ``grid.p_bar`` and ``grid.p``.
    * ``dPi_osm``  resting osmotic difference that drives Kedem-Katchalsky water
      flux. At rest ``dPi_osm == dP_hyd`` ⇒ **zero net trans-membrane flux**.
    * ``p_excess`` the deviatoric Biot pore-pressure field ``grid.p − grid.p_bar``;
      starts at 0 and is what the perturbation excites.

Units: µm, pN, s (1 Pa == 1 pN/µm²). float64 fields, Warp-CUDA only.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np
import warp as wp

from aleph.laws.membrane_surface import (
    build_membrane_mesh,
    build_membrane_hinges,
    membrane_area_kernel,
    helfrich_bending_kernel,
    calibrate_kappa_tilde,
)
from aleph.components.incumbent.membrane_pressure import MembranePressureTraction
from aleph.components.fluid.field_grid import FieldGrid, FLUID, OUTSIDE
from aleph.components.fluid.biot_substrate import BiotSubstrate
from aleph.components.fluid.boundary import MembraneFluxBC
from aleph.components.incumbent.fsi_coupling import SolidDilatationCoupling

# --- Physiological / diagnostic parameters (engine µm-pN-s) --------------------------------------------
R_MEM_UM = 7.5                 # membrane radius [µm] (MCF7 geometry)
DP_HYD_REST_PA = 40.0          # resting hydrostatic difference [Pa] — HeLa DIAGNOSTIC PROXY, not MCF7 production
DPI_OSM_REST_PA = 40.0         # resting osmotic difference [Pa] — equals dP_hyd at rest ⇒ zero net KK flux
# Resting shell tension is DERIVED from Young-Laplace so the baseline is force-balanced (physiological
# operating point, PI 2026-06-04): 2·γ_eff/R = dP_hyd ⇒ γ_eff = dP_hyd·R/2. The bilayer alone (γ_mem≈10
# pN/µm, KB-3.B1.1) cannot balance a 40 Pa turgor — in the real cell the CORTEX carries it — so this
# diagnostic single-shell uses an *effective* pre-tension that lumps the cortical contribution. This is a
# DIAGNOSTIC lump (declared), not a production cortex; the perturbation is measured FROM this balanced rest.
KAPPA_M_PN_UM = 0.0828         # Helfrich bending κ_m [pN·µm] (20 k_BT; KB-3.B1.2)
BIOT_MOBILITY = 5.0e-3         # Darcy mobility k/µ [µm⁴/(pN·s)] (I0-B1)
BIOT_STORAGE_S = 1.0e-4        # storativity 1/M [µm²/pN] (I0-B1)
BIOT_ALPHA = 1.0               # Biot-Willis (I0-B1 ratified)
L_P = 1.6e-8                   # membrane hydraulic conductivity [µm/(s·Pa)] (engine value)
DX_UM = 0.5                    # grid spacing [µm]
GRID_PAD_UM = 3.0              # fluid grid padding beyond the membrane [µm]


# ============================ lane-local glue kernels ============================
@wp.kernel
def _overdamped_relax_kernel(
    pos: wp.array(dtype=wp.vec3d),
    force: wp.array(dtype=wp.vec3d),
    inv_gamma: wp.float64,
    dt: wp.float64,
    pos_new: wp.array(dtype=wp.vec3d),
    vel: wp.array(dtype=wp.vec3d),
):
    """Overdamped quasi-static membrane relaxation ``x_new = x + dt·f/γ``; ``v = (x_new−x)/dt``.

    ``inv_gamma`` (= 1/γ, a *numerical* relaxation mobility [µm/(pN·s)]) is chosen for stability, not a
    physical drag claim — this rig is a quasi-static diagnostic, so the trajectory is a relaxation path.
    """
    n = wp.tid()
    dx = force[n] * (inv_gamma * dt)
    pos_new[n] = pos[n] + dx
    vel[n] = dx / dt


@wp.kernel
def _sub_ref_kernel(force: wp.array(dtype=wp.vec3d), f_ref: wp.array(dtype=wp.vec3d)):
    """Subtract the resting reference force so the relaxed config is an exact zero-stress equilibrium."""
    n = wp.tid()
    force[n] = force[n] - f_ref[n]


@wp.kernel
def _sphere_source_kernel(
    origin: wp.vec3d,
    dx: wp.float64,
    mask: wp.array3d(dtype=wp.int32),
    centre: wp.vec3d,
    radius: wp.float64,
    rate: wp.float64,
    s_water: wp.array3d(dtype=wp.float64),
):
    """Localized interior fluid volumetric source [1/s] in a small sphere of FLUID cells (one-sided patch)."""
    i, j, k = wp.tid()
    if mask[i, j, k] != wp.int32(1):  # FLUID == 1
        return
    x = origin[0] + (wp.float64(i) + wp.float64(0.5)) * dx
    y = origin[1] + (wp.float64(j) + wp.float64(0.5)) * dx
    z = origin[2] + (wp.float64(k) + wp.float64(0.5)) * dx
    d = wp.vec3d(x - centre[0], y - centre[1], z - centre[2])
    if wp.length(d) <= radius:
        s_water[i, j, k] = rate


@wp.kernel
def _sphere_pbump_kernel(
    origin: wp.vec3d,
    dx: wp.float64,
    mask: wp.array3d(dtype=wp.int32),
    centre: wp.vec3d,
    radius: wp.float64,
    amplitude: wp.float64,
    p: wp.array3d(dtype=wp.float64),
):
    """Add a localized pressure excess [Pa] (smooth cosine bump) to FLUID cells in a one-sided patch."""
    i, j, k = wp.tid()
    if mask[i, j, k] != wp.int32(1):
        return
    x = origin[0] + (wp.float64(i) + wp.float64(0.5)) * dx
    y = origin[1] + (wp.float64(j) + wp.float64(0.5)) * dx
    z = origin[2] + (wp.float64(k) + wp.float64(0.5)) * dx
    d = wp.length(wp.vec3d(x - centre[0], y - centre[1], z - centre[2]))
    if d <= radius:
        w = wp.float64(0.5) * (wp.float64(1.0) + wp.cos(wp.float64(3.14159265358979) * d / radius))
        p[i, j, k] = p[i, j, k] + amplitude * w


@wp.kernel
def _sphere_phold_kernel(
    origin: wp.vec3d,
    dx: wp.float64,
    mask: wp.array3d(dtype=wp.int32),
    centre: wp.vec3d,
    radius: wp.float64,
    amplitude: wp.float64,
    p_bar: wp.array3d(dtype=wp.float64),
    p: wp.array3d(dtype=wp.float64),
):
    """CLAMP FLUID cells in a patch to ``p_bar + amplitude·bump`` [Pa] (held Dirichlet-like overpressure)."""
    i, j, k = wp.tid()
    if mask[i, j, k] != wp.int32(1):
        return
    x = origin[0] + (wp.float64(i) + wp.float64(0.5)) * dx
    y = origin[1] + (wp.float64(j) + wp.float64(0.5)) * dx
    z = origin[2] + (wp.float64(k) + wp.float64(0.5)) * dx
    d = wp.length(wp.vec3d(x - centre[0], y - centre[1], z - centre[2]))
    if d <= radius:
        w = wp.float64(0.5) * (wp.float64(1.0) + wp.cos(wp.float64(3.14159265358979) * d / radius))
        p[i, j, k] = p_bar[i, j, k] + amplitude * w


@wp.kernel
def _pv_work_kernel(
    p: wp.array3d(dtype=wp.float64),
    div_vs: wp.array3d(dtype=wp.float64),
    mask: wp.array3d(dtype=wp.int32),
    dV: wp.float64,
    out: wp.array(dtype=wp.float64),
):
    """Fluid boundary-work density integral ``∫_FLUID p·div(v_s) dV`` [pN·µm/s]."""
    i, j, k = wp.tid()
    if mask[i, j, k] != wp.int32(1):
        return
    wp.atomic_add(out, 0, p[i, j, k] * div_vs[i, j, k] * dV)


# ============================ D0 rig ============================
@dataclass
class MembraneCytosolRig:
    """Assembles a spherical membrane shell filled with a Biot-Darcy cytosol field.

    All heavy state is device-resident (Warp-CUDA). Construction requires a CUDA device.
    """

    device: str
    subdivisions: int = 3
    r_mem_um: float = R_MEM_UM
    # --- built state (filled in __post_init__) ---
    mesh: object = field(default=None, init=False)
    grid: FieldGrid = field(default=None, init=False)
    substrate: BiotSubstrate = field(default=None, init=False)
    traction: MembranePressureTraction = field(default=None, init=False)
    flux_bc: MembraneFluxBC = field(default=None, init=False)
    coupling: SolidDilatationCoupling = field(default=None, init=False)

    def __post_init__(self):
        dev = self.device
        # ---- membrane mesh (own node set) ----
        self.mesh = build_membrane_mesh(self.r_mem_um, subdivisions=self.subdivisions)
        self.verts0 = np.ascontiguousarray(self.mesh.verts, np.float64)
        self.faces_np = np.ascontiguousarray(self.mesh.faces, np.int32)
        self.hinges_np = np.ascontiguousarray(build_membrane_hinges(self.mesh.faces), np.int32)
        self.n_nodes = self.verts0.shape[0]
        self.n_faces = self.faces_np.shape[0]
        self.kappa_tilde = float(calibrate_kappa_tilde(KAPPA_M_PN_UM, self.verts0, self.hinges_np))
        # derived pre-tension balancing the resting turgor (Young-Laplace) → force-balanced baseline
        self.gamma_eff = 0.5 * DP_HYD_REST_PA * self.r_mem_um
        # node outward normals (resting sphere) + tributary areas (for node_volume weight)
        self._normals0 = self._vertex_normals(self.verts0, self.mesh.faces)
        self._node_area = self._vertex_areas(self.verts0, self.mesh.faces)

        # ---- device membrane arrays ----
        self.pos_d = wp.array(self.verts0, dtype=wp.vec3d, device=dev)
        self.pos_new_d = wp.zeros(self.n_nodes, dtype=wp.vec3d, device=dev)
        self.force_d = wp.zeros(self.n_nodes, dtype=wp.vec3d, device=dev)
        self.vel_d = wp.zeros(self.n_nodes, dtype=wp.vec3d, device=dev)
        self.faces_d = wp.array(self.faces_np, dtype=wp.int32, device=dev)
        self.hinges_d = wp.array(self.hinges_np, dtype=wp.int32, device=dev)
        # node_volume: boundary-layer voxel A_i·dx (representative fluid volume swept by node) [µm³]
        self.node_volume_d = wp.array(self._node_area * DX_UM, dtype=wp.float64, device=dev)
        self.active_d = wp.zeros(self.n_nodes, dtype=wp.int32, device=dev)  # all active (>=0)
        self.f_rest_d = wp.zeros(self.n_nodes, dtype=wp.vec3d, device=dev)  # resting zero-stress reference

        # ---- fluid grid ----
        half = self.r_mem_um + GRID_PAD_UM
        nx = int(round(2.0 * half / DX_UM))
        origin = (-half, -half, -half)
        self.grid = FieldGrid((nx, nx, nx), DX_UM, origin, device=dev)
        self._nx = nx
        self._origin = np.asarray(origin, np.float64)
        self._build_mask()
        self._seed_resting_pressure()

        self.substrate = BiotSubstrate(self.grid, mobility=BIOT_MOBILITY,
                                       storage_S=BIOT_STORAGE_S, alpha=BIOT_ALPHA)
        self.traction = MembranePressureTraction(grid=self.grid, faces_d=self.faces_d,
                                                 n_faces=self.n_faces, p_ext=0.0)
        self.flux_bc = MembraneFluxBC(self.grid, L_p=L_P, p_ext=0.0, sigma_refl=1.0,
                                      pos=self.pos_d, faces=self.faces_d, n_faces=self.n_faces)
        self.coupling = SolidDilatationCoupling(self.grid)
        self._pv_work_d = wp.zeros(1, dtype=wp.float64, device=dev)

    # ---------------- construction helpers (host) ----------------
    @staticmethod
    def _vertex_normals(verts, faces):
        nrm = np.zeros_like(verts)
        for tri in faces:
            a, b, c = verts[tri[0]], verts[tri[1]], verts[tri[2]]
            fn = np.cross(b - a, c - a)
            for vi in tri:
                nrm[vi] += fn
        # outward for a centred sphere: align with radial
        rad = verts - verts.mean(axis=0)
        flip = np.sum(nrm * rad, axis=1) < 0
        nrm[flip] *= -1.0
        norms = np.linalg.norm(nrm, axis=1, keepdims=True)
        norms[norms == 0] = 1.0
        return nrm / norms

    @staticmethod
    def _vertex_areas(verts, faces):
        area = np.zeros(verts.shape[0])
        for tri in faces:
            a, b, c = verts[tri[0]], verts[tri[1]], verts[tri[2]]
            tri_a = 0.5 * np.linalg.norm(np.cross(b - a, c - a))
            for vi in tri:
                area[vi] += tri_a / 3.0
        return np.ascontiguousarray(area, np.float64)

    def _build_mask(self):
        nx = self._nx
        idx = (np.arange(nx) + 0.5) * DX_UM + self._origin[0]
        X, Y, Z = np.meshgrid(idx, idx, idx, indexing="ij")
        centre = self.verts0.mean(axis=0)
        r = np.sqrt((X - centre[0]) ** 2 + (Y - centre[1]) ** 2 + (Z - centre[2]) ** 2)
        # FLUID out to a half-cell collar past the membrane so the one-sided interior pressure trace
        # (samples at centroid − n̂·dx/2 and −n̂·dx) is valid at EVERY membrane face → closed-shell
        # traction balances (G0.3). The collar cells hold pressure only; the physical fluid is inside.
        mask = np.where(r < (self.r_mem_um + 0.75 * DX_UM), FLUID, OUTSIDE).astype(np.int32)
        self.grid.set_mask(np.ascontiguousarray(mask))
        self._mask_np = mask
        self.n_fluid = int((mask == FLUID).sum())

    def _seed_resting_pressure(self):
        """Seed p_bar = p = dP_hyd (HeLa proxy) on FLUID cells; p_excess = 0."""
        p = np.where(self._mask_np == FLUID, DP_HYD_REST_PA, 0.0).astype(np.float64)
        self.grid.set_pressure(np.ascontiguousarray(p))
        with wp.ScopedDevice(self.device):
            self.grid.p_bar = wp.array(np.ascontiguousarray(p), dtype=wp.float64)

    # ---------------- membrane force assembly (device) ----------------
    def accumulate_membrane_force(self, pos_d=None):
        """Zero and rebuild the membrane nodal force: pressure traction + tension + bending."""
        pos = self.pos_d if pos_d is None else pos_d
        self.force_d.zero_()
        # pressure → outward normal traction (reads grid.p via interior surface trace)
        self.traction.accumulate(pos, self.force_d)
        # in-plane tension (area-elastic, taut sheet)
        wp.launch(membrane_area_kernel, dim=self.n_faces,
                  inputs=[pos, self.faces_d, wp.float64(self.gamma_eff), self.force_d],
                  device=self.device)
        # Helfrich bending
        wp.launch(helfrich_bending_kernel, dim=self.hinges_np.shape[0],
                  inputs=[pos, self.hinges_d, wp.float64(self.kappa_tilde), self.force_d],
                  device=self.device)
        # subtract the resting zero-stress reference (zero until relax_to_rest sets it)
        wp.launch(_sub_ref_kernel, dim=self.n_nodes,
                  inputs=[self.force_d, self.f_rest_d], device=self.device)

    def relax_to_rest(self, n_iter: int = 1200, step: float = 5.0e-4) -> float:
        """Settle the shell to its DISCRETE mechanical equilibrium at resting pressure (fluid frozen).

        The continuum sphere is the equilibrium, but the discrete area-tension and pressure-traction
        operators leave a small per-node residual at the 12 icosphere defect vertices. Relaxing removes
        it so ``disp = pos − verts0`` measures the perturbation response from a clean force-balanced rest.
        Returns the final max per-node residual force [pN]. Updates ``verts0``/normals/areas to the
        relaxed reference.
        """
        self._seed_resting_pressure()
        self.f_rest_d.zero_()  # relax against the raw (un-referenced) force field
        for _ in range(n_iter):
            self.accumulate_membrane_force()
            wp.launch(_overdamped_relax_kernel, dim=self.n_nodes,
                      inputs=[self.pos_d, self.force_d, wp.float64(step), wp.float64(1.0),
                              self.pos_new_d, self.vel_d], device=self.device)
            self.pos_d.assign(self.pos_new_d)
        self.verts0 = self.pos_d.numpy().copy()
        self._normals0 = self._vertex_normals(self.verts0, self.mesh.faces)
        self._node_area = self._vertex_areas(self.verts0, self.mesh.faces)
        self.node_volume_d = wp.array(self._node_area * DX_UM, dtype=wp.float64, device=self.device)
        # freeze the residual as the zero-stress reference: rest is now an exact equilibrium
        self.accumulate_membrane_force()             # f_rest still 0 → raw residual in force_d
        raw_resid = float(np.max(np.linalg.norm(self.force_d.numpy(), axis=1)))
        self.f_rest_d.assign(self.force_d)
        self.accumulate_membrane_force()             # now referenced → net ~0
        net_resid = float(np.max(np.linalg.norm(self.force_d.numpy(), axis=1)))
        self._rest_raw_residual = raw_resid
        return net_resid

    # ---------------- perturbation ----------------
    def set_patch_source(self, centre, radius, rate):
        """Localized interior volumetric fluid source [1/s] in a small FLUID sphere (one-sided patch)."""
        self.grid.s_water.zero_()
        wp.launch(_sphere_source_kernel, dim=(self._nx, self._nx, self._nx),
                  inputs=[wp.vec3d(*self._origin.tolist()), wp.float64(DX_UM), self.grid.mask,
                          wp.vec3d(*centre), wp.float64(radius), wp.float64(rate), self.grid.s_water],
                  device=self.device)

    def clear_source(self):
        self.grid.s_water.zero_()

    def add_pressure_bump(self, centre, radius, amplitude):
        """Add a one-shot localized pressure excess [Pa] (cosine bump) to the cytosol field (perturbation)."""
        wp.launch(_sphere_pbump_kernel, dim=(self._nx, self._nx, self._nx),
                  inputs=[wp.vec3d(*self._origin.tolist()), wp.float64(DX_UM), self.grid.mask,
                          wp.vec3d(*centre), wp.float64(radius), wp.float64(amplitude), self.grid.p],
                  device=self.device)

    def hold_pressure_bump(self, centre, radius, amplitude):
        """Clamp a patch to a held overpressure ``p_bar + amplitude·bump`` [Pa] (sustained perturbation)."""
        wp.launch(_sphere_phold_kernel, dim=(self._nx, self._nx, self._nx),
                  inputs=[wp.vec3d(*self._origin.tolist()), wp.float64(DX_UM), self.grid.mask,
                          wp.vec3d(*centre), wp.float64(radius), wp.float64(amplitude),
                          self.grid.p_bar, self.grid.p], device=self.device)

    # ---------------- outer physical step ----------------
    def outer_step(self, dt_phys: float, inv_gamma: float, dpi_osm: float = DPI_OSM_REST_PA,
                   move_membrane: bool = True, n_mech: int = 12,
                   hold_patch=None) -> dict:
        """One outer physical step.

        Sequence: (1) quasi-static membrane relaxation — ``n_mech`` stable overdamped sub-iterations at the
        current cytosol pressure (each within the explicit stiff-tension stability limit); (2) membrane
        velocity ``v_s = Δx/dt_phys`` → moving-boundary source ``div(v_s)``; (3) Biot sub-cycle advancing the
        pore pressure, optionally clamping a held patch overpressure (``hold_patch=(centre,radius,amp)``).
        Returns a small host diagnostic dict (post-loop readback only).
        """
        # Co-march the membrane INSIDE the Biot sub-cycle. The undrained cytosol (S=1e-4 ⇒ nearly
        # incompressible) turns any membrane motion into a large div(v_s)/S pressure response, so the
        # moving boundary must advance at the small fluid sub-step dt_sub to stay resolved and bounded.
        dt_cfl = self.substrate.cfl_dt(0.9)
        n_sub = max(1, int(math.ceil(dt_phys / dt_cfl)))
        dt_sub = dt_phys / n_sub
        w_membrane = 0.0
        for _ in range(n_sub):
            if hold_patch is not None:
                self.hold_pressure_bump(*hold_patch)
            self.grid.div_vs.zero_()
            if move_membrane:
                # one overdamped membrane increment at the current pressure, over dt_sub
                self.accumulate_membrane_force()
                wp.launch(_overdamped_relax_kernel, dim=self.n_nodes,
                          inputs=[self.pos_d, self.force_d, wp.float64(inv_gamma), wp.float64(dt_sub),
                                  self.pos_new_d, self.vel_d], device=self.device)
                # membrane motion → moving-boundary source div(v_s) BEFORE committing the move
                self.coupling.update(self.pos_d, self.vel_d, self.node_volume_d, self.active_d)
                # membrane nodal work rate Σ f_pressure·v_s accumulated over the sub-cycle
                f_p = wp.zeros(self.n_nodes, dtype=wp.vec3d, device=self.device)
                self.traction.accumulate(self.pos_d, f_p)
                w_membrane += float(np.sum(f_p.numpy() * self.vel_d.numpy())) * (dt_sub / dt_phys)
                self.pos_d.assign(self.pos_new_d)
            self.flux_bc.apply(dpi_osm)
            self.substrate.step(dt_sub)
        if hold_patch is not None:
            self.hold_pressure_bump(*hold_patch)
        # fluid boundary work ∫ p·div(v_s) dV over the last sub-step (representative rate)
        self._pv_work_d.zero_()
        wp.launch(_pv_work_kernel, dim=(self._nx, self._nx, self._nx),
                  inputs=[self.grid.p, self.grid.div_vs, self.grid.mask, wp.float64(DX_UM ** 3),
                          self._pv_work_d], device=self.device)
        w_fluid = float(self._pv_work_d.numpy()[0])
        return {
            "n_sub": n_sub,
            "dt_sub": dt_sub,
            "w_membrane_pn_um_s": w_membrane,
            "w_fluid_pn_um_s": w_fluid,
            "integrated_flux_um3_s": float(self.flux_bc.integrated_flux_d.numpy()[0]),
            "total_content": self.grid.total_content(BIOT_STORAGE_S),
        }

    # ---------------- readback helpers ----------------
    def membrane_state(self):
        pos = self.pos_d.numpy()
        disp = pos - self.verts0
        self.accumulate_membrane_force()
        # isolate the pressure traction alone for the vector field
        f_tot = wp.zeros(self.n_nodes, dtype=wp.vec3d, device=self.device)
        self.traction.accumulate(self.pos_d, f_tot)
        return {
            "pos": pos,
            "disp": disp,
            "traction": f_tot.numpy(),
            "normals0": self._normals0,
        }

    def pressure_midplane(self, axis: int = 2):
        p = self.grid.pressure_to_host()
        k = self._nx // 2
        if axis == 2:
            return p[:, :, k]
        if axis == 1:
            return p[:, k, :]
        return p[k, :, :]

    def patch_centre(self, frac: float = 0.78):
        """Sub-membrane +x patch centre — close enough to the +x membrane that a held overpressure
        directly loads that region (visible asymmetric traction), not just the deep interior."""
        c = self.verts0.mean(axis=0)
        return (c[0] + frac * self.r_mem_um, c[1], c[2])
