r"""Deformable-mesh NUCLEUS + MEMBRANE compartments for the assembled Active Cell (assembly milestone-2).

Two §1.4 force primitives the lead-owned integrator (``ac.cell.driver``) sums into the global node force,
parallel to ``MyosinForce`` / ``StericForce`` / ``PressureCoupling``:

  * :class:`NucleusCompartment` (P4) — a triangulated nuclear-envelope shell (the Helfrich machinery, NOT a
    radial bead ball). Force families, in the §3 order the nucleus INTEGRATION.md fixes: reused
    ``helfrich_bending_kernel`` (κ̃=8πκ_NE/Σ_ref) + per-face ``lamina_areal_tension_kernel`` (framework-#6
    lamin split) + device-resident nucleoplasm volume penalty (ν→½) + ``linc_tether_kernel``
    (nucleus↔cytoskeleton — the SINGLE load path). Device-predicated rupture commits once per OUTER step.
    Retires ``ff.network_warp.nucleus_shell_kernel`` by non-use (same commit that wires this in).

  * :class:`MembraneCompartment` (P4) — the independent plasma-membrane Helfrich sheet
    (``ff.membrane_surface``): reused ``helfrich_bending_kernel`` + in-plane ``membrane_area_kernel`` (γ_mem
    plateau; the K_A reservoir upturn stays default-OFF per master hard-truth #8) + commit-safe explicit ERM
    tethers (membrane↔cortex; the membrane RIDES on the cortex, force-free at rest). It is also the osmotic envelope
    for the fluid domain: a persistent CUDA ``wp.Mesh`` classifies the live membrane from build time onward;
    the analytic sphere provider is retained only for isolated diagnostics.

  * :class:`MeshNucleusMaskProvider` — the fluid-domain adapter (``NucleusMaskProvider.classify_nucleus``)
    that masks the LIVE deformable oblate nucleus mesh (about ``centre_nuc``) as the relative-no-flux inner
    inclusion (NOT a static sphere R_nuc). The production path delegates to a device-resident Warp mesh;
    the float64 host winding-number implementation remains an acceptance oracle only.

GAP magnitudes (I0-B2, PI-authored — do NOT choose, do NOT tune to a gate): lamin-A/C vs lamin-B split,
envelope rupture strain, nucleoplasm viscosity, k_linc ("8 pN is a TENSION not a stiffness"), MCF7 oblate
aspect, and the nuclear radius / N:C ratio. Every one is passed as a clearly-labeled provisional/TEST value.
At the RESTING baseline none enters the force balance: areal strain ε=0 ⇒ σ=0; V=V0 ⇒ p_vol=0; LINC/ERM rest
= the formation length ⇒ force-free — exactly the myosin-heads-unbound / σ_EV-inactive discipline. The
magnitude verdicts are HELD until I0-B2 closes; the STRUCTURAL gates (resting stability, no-flux, net force)
run now.

Runtime: Warp CUDA GPU only (I0-A). Built on the gbook A5000; the dev Mac import/syntax-checks only (the
FieldGrid / device arrays reject a non-CUDA device). Host geometry construction and LINC/ERM pairing are pure
NumPy; runtime winding classification is CUDA-resident, with the float64 host winding implementation used only
as an acceptance oracle.

Units: FF µm·pN·s. Areal moduli [pN/µm]; K_vol [pN/µm²]; k_linc/k_erm [pN/µm]; κ̃ [pN·µm].
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import numpy.typing as npt
import warp as wp

from aleph.components.incumbent.erm_tether import (
    ERMBellKinetics,
    erm_bell_force_kernel,
    erm_bell_kmc_kernel,
    erm_tether_force_kernel,
    erm_tether_rupture_kernel,
    increment_erm_epoch_if_accepted_kernel,
)
from aleph.components.incumbent.live_mesh_domain import LiveMeshNucleusMaskProvider
from aleph.components.incumbent.preload_contract import triangle_surface_area

# ── ac/nucleus track primitives — read-only ─────────────────────────────────────────────────────────
from aleph.components.nucleus.envelope import (
    NucleusMesh,
    build_nucleus,
    conditional_rupture_update_kernel,
    lamina_areal_tension_kernel,
    linc_tether_kernel,
    nucleoplasm_volume_reduce_kernel,
)
from aleph.components.nucleus.lamina_analytic import LaminaParams
from aleph.components.nucleus.mask_provider import solid_angle_winding_number

# ── reused ff/ Helfrich sheet (membrane) — read-only ────────────────────────────────────────────────
from aleph.laws.membrane_surface import (
    MembraneMesh,
    build_membrane_hinges,
    build_membrane_mesh,
    erm_rupture_force,
    helfrich_bending_kernel,
    membrane_area_kernel,
)
from aleph.laws.membrane_surface import (
    calibrate_kappa_tilde as membrane_calibrate_kappa_tilde,
)

__all__ = [
    "NUCLEUS_DEFAULTS",
    "MEMBRANE_DEFAULTS",
    "MeshNucleusMaskProvider",
    "NucleusCompartment",
    "MembraneCompartment",
    "build_nucleus_compartment",
    "build_membrane_compartment",
    "density_resolved_erm_pairs",
    "nucleoplasm_volume_force_arr_kernel",
]

# ─────────────────────────────────────────────────────────────────────────────────────────────────
# Provisional / TEST magnitudes — every one is an I0-B2 PI GAP (surfaced, never tuned to a gate).
# None enters the RESTING result (ε=0, V=V0, tethers force-free). See module docstring.
# ─────────────────────────────────────────────────────────────────────────────────────────────────
NUCLEUS_DEFAULTS = dict(
    r_eq_um=5.0,            # MCF7 nuclear equatorial radius [µm] — provisional geometry (N:C vol ~0.30); PI GAP
    aspect=1.0,            # resting SUSPENDED baseline = sphere (adherent oblate flatten is I7); aspect PI GAP
    subdivisions=3,        # icosphere level (642 nodes / 1280 faces) — continuum shell resolution
    # LaminaParams (framework-#6 lamin split) — ALL I0-B2 PI GAP, TEST values; ε=0 ⇒ σ=0 at rest:
    k_chrom=15.0,          # chromatin internal-net areal modulus [pN/µm] (soft) — TEST
    k_lamin_b=15.0,        # lamin-B baseline meshwork areal modulus [pN/µm] (soft) — TEST
    k_lamin_ac=150.0,      # lamin-A/C strain-stiffening areal modulus [pN/µm] (stiff, ratio 5×) — TEST
    knee_strain=0.10,      # lamin-A/C engagement strain (KB-3.B2.2 ~0.10) — TEST
    eps_rupture=0.50,      # EMERGENT envelope rupture strain — PI GAP, TEST (never tuned)
    kappa_ne=0.0856,       # envelope bending rigidity [pN·µm] (~20 kBT) — TEST
    k_vol=1.0e3,           # nucleoplasm bulk penalty [pN/µm²] (ν→½ incompressible) — TEST; p_vol=0 at V=V0
    k_linc=1.0e2,          # LINC tether stiffness [pN/µm] — PI GAP ("8 pN is a TENSION not a stiffness")
    linc_stiffening=0.0,   # nesprin nonlinear stiffening [1/µm²] — TEST (Hookean baseline)
    linc_reach_um=3.0,     # nucleus↔cortex LINC pairing reach [µm] (spans the cortex-nucleus gap at rest)
    linc_per_node=1,       # LINC tethers seeded per nucleus surface node (baseline coupling; cap→LINC is I7)
)

MEMBRANE_DEFAULTS = dict(
    subdivisions=3,        # icosphere level (642 nodes / 1280 faces) — matched to the cortex shell scale
    gamma_mem=10.0,        # in-plane bilayer tension γ_mem [pN/µm] (KB-3.B1.1) — LIT-anchored plateau
    kappa_m=0.0828,        # Helfrich bending κ_m [pN·µm] (20 kBT; KB-3.B1.2) — LIT-anchored
    k_erm=4.6e3,           # single-ezrin–F-actin linker stiffness [pN/µm] = 4.6 pN/nm (PI-ratified 2026-07-21).
                           # Provenance precision: Braunger 2014 JBC reports most-probable system stiffness
                           # 2.2 pN/nm; Braunger's event-statistics dissertation analysis derives the 4.6 pN/nm
                           # single-bond estimate with N≈1. KB-3.B1.6 registration is PENDING, not yet SoT.
                           # Replaces the retired provisional 1e2 (a bleb-PDE continuum value, ~46× too soft).
    erm_reach_um=1.0,      # membrane↔cortex ERM pairing reach [µm] (coincident shells at rest)
    erm_radial_pairing=False,  # pair each membrane node to its most-RADIAL cortex node (not Euclidean-nearest)
                           # so the tether transmits the radial turgor load between the shells. Default off
                           # (backward-compat); the resting physiological baseline turns it on.
    # The Notion Contract-Graph has no MCF7 ERM linker-density datum. ``None`` preserves the diagnostic
    # one-tether-per-membrane-node build, but that build is not a physiological preload. Production must
    # supply a sourced, PI-ratified density and source label.
    erm_density_per_um2=None,
    erm_density_source="",
    erm_bell_kinetics=None,  # no biological defaults: MCF7 k_on/k_off0/F0/capture are Contract-Graph GAPs
    with_area_tension=True,  # γ_mem plateau ON (physiological); K_A reservoir upturn stays OFF (hard-truth #8)
)


# ─────────────────────────────────────────────────────────────────────────────────────────────────
# Device-resident nucleoplasm volume pressure (keeps NG-6 zero-roundtrip: p_vol computed on device,
# never read to host inside the inner solve).
# ─────────────────────────────────────────────────────────────────────────────────────────────────
@wp.kernel
def _zero_scalar_kernel(a: wp.array(dtype=wp.float64)) -> None:
    a[wp.tid()] = wp.float64(0.0)


@wp.kernel
def _nucleoplasm_pressure_kernel(
    vol: wp.array(dtype=wp.float64), v0: wp.float64, k_vol: wp.float64,
    p_vol_out: wp.array(dtype=wp.float64),
) -> None:
    """``p_vol = −K_vol (V − V0)/V0`` [pN/µm²] on device (positive when V<V0 ⇒ inflate). One thread."""
    p_vol_out[0] = -k_vol * (vol[0] - v0) / v0


@wp.kernel
def nucleoplasm_volume_force_arr_kernel(
    npos: wp.array(dtype=wp.vec3d), faces: wp.array(dtype=wp.int32, ndim=2),
    p_vol_arr: wp.array(dtype=wp.float64), nforce: wp.array(dtype=wp.vec3d),
) -> None:
    """Incompressible nucleoplasm force with a DEVICE-array pressure (reads ``p_vol_arr[0]``, so the
    volume→pressure→force chain stays fully on-device — the NG-6 zero-roundtrip variant of
    ``envelope.nucleoplasm_volume_force_kernel``). Force = ``p_vol·∂V/∂x_i`` along the volume gradient."""
    t = wp.tid()
    pv = p_vol_arr[0]
    i0 = faces[t, 0]
    i1 = faces[t, 1]
    i2 = faces[t, 2]
    p0 = npos[i0]
    p1 = npos[i1]
    p2 = npos[i2]
    wp.atomic_add(nforce, i0, pv * wp.cross(p1, p2) / wp.float64(6.0))
    wp.atomic_add(nforce, i1, pv * wp.cross(p2, p0) / wp.float64(6.0))
    wp.atomic_add(nforce, i2, pv * wp.cross(p0, p1) / wp.float64(6.0))


# ─────────────────────────────────────────────────────────────────────────────────────────────────
# Fluid-domain adapter — the live deformable oblate mesh as the relative-no-flux inner boundary.
# ─────────────────────────────────────────────────────────────────────────────────────────────────
class MeshNucleusMaskProvider:
    """Adapt the nucleus mesh to ``domain.NucleusMaskProvider.classify_nucleus(grid, mask)``.

    The nucleus track ships ``mask_provider.DeformableNucleusMaskProvider`` (``inside_mask``/``boundary``/
    ``update`` on the ``NucleusBoundaryProvider`` protocol); the fluid ``Domain`` consumes the DIFFERENT
    ``classify_nucleus(grid, mask)`` device-side signature. This lead-owned glue bridges them WITHOUT editing
    either track file. ``_inside_grid`` evaluates the float64 host winding number for an independent oracle;
    production classification is a CUDA D2D refresh + device winding query after :meth:`bind_live`.
    """

    def __init__(self, verts: npt.NDArray[np.float64], faces: npt.NDArray[np.int64], *, aabb_margin_um: float = 1.0):
        self.faces = np.ascontiguousarray(faces, np.int64)
        self.aabb_margin = float(aabb_margin_um)
        self._live: LiveMeshNucleusMaskProvider | None = None
        self.update(verts)

    def update(self, verts: npt.NDArray[np.float64]) -> None:
        self.verts = np.ascontiguousarray(verts, np.float64)

    def bind_live(self, pos_d, node_off: int, n_verts: int) -> None:
        """Bind the authoritative GLOBAL CUDA position array for live, GPU-resident classification."""
        self._live = LiveMeshNucleusMaskProvider(
            pos_d, int(node_off), int(n_verts), self.faces, reference_verts=self.verts)

    def _inside_grid(self, grid) -> np.ndarray:
        """Host int32 (nx,ny,nz) occupancy: 1 where a cell centre is inside the current envelope, else 0.

        AABB-restricted: only cells within the nucleus bounding box (+ margin) are winding-number tested; the
        rest are trivially outside — keeps the O(N_cells·N_faces) winding number tractable out-of-hot-loop.
        """
        nx, ny, nz = grid.shape
        o = np.asarray(grid.origin, np.float64)
        dx = float(grid.dx)
        lo = self.verts.min(axis=0) - self.aabb_margin
        hi = self.verts.max(axis=0) + self.aabb_margin
        # index ranges of the AABB on the grid (clamped)
        ilo = np.clip(np.floor((lo - o) / dx).astype(int), 0, [nx, ny, nz])
        ihi = np.clip(np.ceil((hi - o) / dx).astype(int) + 1, 0, [nx, ny, nz])
        inside = np.zeros((nx, ny, nz), np.int32)
        if np.any(ihi <= ilo):
            return inside
        ii = np.arange(ilo[0], ihi[0])
        jj = np.arange(ilo[1], ihi[1])
        kk = np.arange(ilo[2], ihi[2])
        gi, gj, gk = np.meshgrid(ii, jj, kk, indexing="ij")
        pts = np.stack([o[0] + gi.ravel() * dx, o[1] + gj.ravel() * dx, o[2] + gk.ravel() * dx], axis=1)
        w = solid_angle_winding_number(pts, self.verts, self.faces)
        occ = (w >= 0.5).reshape(gi.shape)
        inside[ilo[0]:ihi[0], ilo[1]:ihi[1], ilo[2]:ihi[2]] = occ.astype(np.int32)
        return inside

    def classify_nucleus(self, grid, mask: wp.array) -> None:
        if self._live is None:
            raise RuntimeError("MeshNucleusMaskProvider must bind_live() before production classification")
        self._live.classify_nucleus(grid, mask)

    @property
    def query_failures_d(self) -> wp.array:
        """Device counter from the most recent production classification."""
        if self._live is None:
            raise RuntimeError("MeshNucleusMaskProvider has no live device binding")
        return self._live.query_failures_d

    @property
    def surface_ties_d(self) -> wp.array:
        """Device count of grid centres resolved by the documented closed-surface tie policy."""
        if self._live is None:
            raise RuntimeError("MeshNucleusMaskProvider has no live device binding")
        return self._live.surface_ties_d

    @property
    def query_failure_mask_d(self) -> wp.array:
        """Per-cell query failures from the most recent live classification (gate readback only)."""
        if self._live is None or self._live.query_failure_mask_d is None:
            raise RuntimeError("MeshNucleusMaskProvider has not classified a live device grid")
        return self._live.query_failure_mask_d

    @property
    def surface_tie_mask_d(self) -> wp.array:
        """Per-cell surface ties from the most recent live classification (gate readback only)."""
        if self._live is None or self._live.surface_tie_mask_d is None:
            raise RuntimeError("MeshNucleusMaskProvider has not classified a live device grid")
        return self._live.surface_tie_mask_d


# ─────────────────────────────────────────────────────────────────────────────────────────────────
# Nucleus compartment
# ─────────────────────────────────────────────────────────────────────────────────────────────────
@dataclass
class NucleusCompartment:
    """Deformable-mesh nucleus force primitive — ``accumulate(pos, f)`` (§1.4) + ``rupture_step(pos)``.

    Holds the GLOBAL-indexed device topology (envelope nodes live in ``pos[node_off : node_off+n_verts]``),
    the calibrated κ̃, per-face reference areas + rupture flags, the reference volume V0 + K_vol, and the LINC
    tether pairs (nucleus node ↔ cortex anchor). Built by :func:`build_nucleus_compartment`.
    """
    device: str
    node_off: int
    n_verts: int
    mesh: NucleusMesh
    # device arrays (GLOBAL indices)
    faces_d: wp.array          # (Nf,3) int32
    hinges_d: wp.array         # (Nh,4) int32
    a0_d: wp.array             # (Nf,) float64
    ruptured_d: wp.array       # (Nf,) int32
    linc_n_d: wp.array         # (Nl,) int32  nucleus node
    linc_a_d: wp.array         # (Nl,) int32  cortex anchor
    linc_rest_d: wp.array      # (Nl,) float64
    # device scalars for the on-device volume→pressure→force chain
    vol_d: wp.array            # (1,) float64
    p_vol_d: wp.array          # (1,) float64
    # magnitudes (provisional / GAP)
    kappa_tilde: float
    v0_ref: float
    k_vol: float
    k_soft: float
    k_ac: float
    knee: float
    eps_rupt: float
    k_linc: float
    linc_stiffening: float
    mean_edge_um: float
    n_faces: int = field(init=False)
    n_hinges: int = field(init=False)
    n_linc: int = field(init=False)

    def __post_init__(self):
        self.n_faces = int(self.faces_d.shape[0])
        self.n_hinges = int(self.hinges_d.shape[0])
        self.n_linc = int(self.linc_n_d.shape[0])

    def accumulate(self, pos: wp.array, f: wp.array) -> None:
        """Sum the nucleus force families into the GLOBAL node force ``f`` (§3 order; device-resident)."""
        d = self.device
        # 1. Helfrich bending (reused ff kernel; κ̃ = 8πκ_NE/Σ_ref)
        if self.n_hinges:
            wp.launch(helfrich_bending_kernel, dim=self.n_hinges,
                      inputs=[pos, self.hinges_d, wp.float64(self.kappa_tilde)], outputs=[f], device=d)
        # 2. per-face lamina areal tension (framework-#6 lamin split; ε=0 ⇒ σ=0 at rest)
        if self.n_faces:
            wp.launch(lamina_areal_tension_kernel, dim=self.n_faces,
                      inputs=[pos, self.faces_d, self.a0_d, self.ruptured_d,
                              wp.float64(self.k_soft), wp.float64(self.k_ac), wp.float64(self.knee),
                              wp.float64(self.eps_rupt)], outputs=[f], device=d)
            # 3. nucleoplasm incompressible volume penalty (device-resident p_vol; 0 at V=V0)
            wp.launch(_zero_scalar_kernel, dim=1, inputs=[self.vol_d], device=d)
            wp.launch(nucleoplasm_volume_reduce_kernel, dim=self.n_faces,
                      inputs=[pos, self.faces_d], outputs=[self.vol_d], device=d)
            wp.launch(_nucleoplasm_pressure_kernel, dim=1,
                      inputs=[self.vol_d, wp.float64(self.v0_ref), wp.float64(self.k_vol)],
                      outputs=[self.p_vol_d], device=d)
            wp.launch(nucleoplasm_volume_force_arr_kernel, dim=self.n_faces,
                      inputs=[pos, self.faces_d, self.p_vol_d], outputs=[f], device=d)
        # 4. LINC tether (the SINGLE nucleus↔cytoskeleton load path; force-free at rest). npos=apos=pos,
        #    nforce=aforce=f, indices GLOBAL ⇒ the Newton pair lands on the right global nodes.
        if self.n_linc:
            wp.launch(linc_tether_kernel, dim=self.n_linc,
                      inputs=[pos, pos, self.linc_n_d, self.linc_a_d, self.linc_rest_d,
                              wp.float64(self.k_linc), wp.float64(self.linc_stiffening)],
                      outputs=[f, f], device=d)

    def rupture_step(self, pos: wp.array, accepted: wp.array) -> None:
        """Commit per-face rupture once, predicated by the device outer-acceptance latch."""
        if self.n_faces:
            wp.launch(conditional_rupture_update_kernel, dim=self.n_faces,
                      inputs=[pos, self.faces_d, self.a0_d, wp.float64(self.eps_rupt), accepted,
                              self.ruptured_d],
                      device=self.device)


def mesh_mean_edge(verts: np.ndarray, faces: np.ndarray) -> float:
    """Mean triangle-edge length [µm] — the discretization length ℓ for the bending CFL term κ̃/ℓ³."""
    p0 = verts[faces[:, 0]]
    p1 = verts[faces[:, 1]]
    p2 = verts[faces[:, 2]]
    e = np.concatenate([np.linalg.norm(p1 - p0, axis=1), np.linalg.norm(p2 - p1, axis=1),
                        np.linalg.norm(p0 - p2, axis=1)])
    return float(e.mean()) if e.size else 1.0


def _pair_to_nearest(src_pts: np.ndarray, dst_pts: np.ndarray, reach: float) -> tuple[np.ndarray, np.ndarray]:
    """Nearest-``dst`` index for each ``src`` point within ``reach`` → (src_idx, dst_idx) local arrays."""
    if src_pts.shape[0] == 0 or dst_pts.shape[0] == 0:
        return np.zeros(0, np.int64), np.zeros(0, np.int64)
    try:
        from scipy.spatial import cKDTree
        d, j = cKDTree(dst_pts).query(src_pts, k=1)
    except Exception:  # noqa: BLE001 — numpy fallback (no scipy)
        dd = np.linalg.norm(dst_pts[None, :, :] - src_pts[:, None, :], axis=2)
        j = dd.argmin(axis=1)
        d = dd[np.arange(dd.shape[0]), j]
    keep = d <= reach
    si = np.nonzero(keep)[0].astype(np.int64)
    return si, np.asarray(j, np.int64)[keep]


def _pair_radial(
    src_pts: np.ndarray, dst_pts: np.ndarray, centre: np.ndarray, reach: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Radial pairing: each ``src`` (membrane) node → the ``dst`` (cortex) node whose OUTWARD RADIAL
    direction from ``centre`` best aligns with the src node's, so the tether opposes the radial pressure.

    The Euclidean-nearest cortex node is often laterally offset (the radial membrane↔cortex gap ~0.1 µm is
    comparable to the cortex node spacing), giving an off-radial tether that cannot transmit the radial
    turgor load between the two shells. Nearest-neighbour on the unit-sphere directions is monotone in the
    angular distance, so a KD-tree on the direction vectors yields the most-radial partner and scales to the
    full-native populations (no dense angle matrix). The ``reach`` guard uses the true Euclidean distance.
    """
    if src_pts.shape[0] == 0 or dst_pts.shape[0] == 0:
        return np.zeros(0, np.int64), np.zeros(0, np.int64)
    c = np.asarray(centre, dtype=np.float64)
    su = src_pts - c
    du = dst_pts - c
    su = su / np.linalg.norm(su, axis=1, keepdims=True)
    du = du / np.linalg.norm(du, axis=1, keepdims=True)
    try:
        from scipy.spatial import cKDTree
        _, j = cKDTree(du).query(su, k=1)
    except Exception:  # noqa: BLE001 — numpy fallback (no scipy)
        j = (su @ du.T).argmax(axis=1)
    j = np.asarray(j, np.int64)
    d = np.linalg.norm(src_pts - dst_pts[j], axis=1)
    keep = d <= reach
    si = np.nonzero(keep)[0].astype(np.int64)
    return si, j[keep]


def density_resolved_erm_pairs(
    membrane_verts: np.ndarray,
    cortex_pts: np.ndarray,
    *,
    surface_area_um2: float,
    density_per_um2: float,
    reach_um: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Create one explicit ERM state per density-derived linker.

    The cortex builder stores seven contiguous nodes per physiological F-actin filament. Sampling evenly in
    that deterministic node order selects the requested count without a tunable random seed; every selected
    cortex node is used at most once. Each linker terminates on the nearest membrane finite-element node, so
    multiple molecular states may load one continuum membrane degree of freedom while retaining independent
    cortex endpoints, rest lengths, and rupture states.

    This helper supplies mechanism, not a default density. The caller must provide a sourced, PI-ratified
    ``density_per_um2``; :data:`MEMBRANE_DEFAULTS` intentionally leaves it ``None``.
    """
    mem = np.asarray(membrane_verts, dtype=np.float64)
    cortex = np.asarray(cortex_pts, dtype=np.float64)
    if mem.ndim != 2 or mem.shape[1] != 3 or cortex.ndim != 2 or cortex.shape[1] != 3:
        raise ValueError("membrane_verts and cortex_pts must both be (N, 3)")
    if not np.isfinite(surface_area_um2) or surface_area_um2 <= 0.0:
        raise ValueError("surface_area_um2 must be finite and positive")
    if not np.isfinite(density_per_um2) or density_per_um2 < 0.0:
        raise ValueError("density_per_um2 must be finite and nonnegative")
    if not np.isfinite(reach_um) or reach_um <= 0.0:
        raise ValueError("reach_um must be finite and positive")

    n_target = int(np.floor(density_per_um2 * surface_area_um2 + 0.5))
    if n_target == 0:
        return np.zeros(0, np.int64), np.zeros(0, np.int64)
    if n_target > cortex.shape[0]:
        raise ValueError(
            f"density requests {n_target} explicit ERM linkers but only {cortex.shape[0]} cortex nodes exist; "
            "increase mechanistic cortex resolution rather than duplicating a cortex endpoint"
        )

    # Midpoint-stratified indices are unique whenever n_target <= n_cortex.
    cortex_idx = np.floor(
        (np.arange(n_target, dtype=np.float64) + 0.5) * cortex.shape[0] / n_target
    ).astype(np.int64)
    try:
        from scipy.spatial import cKDTree

        dist, mem_idx = cKDTree(mem).query(cortex[cortex_idx], k=1)
    except Exception:  # noqa: BLE001 - deterministic NumPy fallback for small test fixtures
        if mem.shape[0] * n_target > 5_000_000:
            raise RuntimeError("SciPy cKDTree is required for a native-density ERM build") from None
        delta = mem[None, :, :] - cortex[cortex_idx, None, :]
        dist2 = np.einsum("...i,...i->...", delta, delta)
        mem_idx = dist2.argmin(axis=1)
        dist = np.sqrt(dist2[np.arange(n_target), mem_idx])
    if np.any(np.asarray(dist) > reach_um):
        n_bad = int(np.count_nonzero(np.asarray(dist) > reach_um))
        raise ValueError(
            f"{n_bad} density-derived ERM linkers exceed the {reach_um:g} um physical reach; "
            "refine the membrane mesh instead of silently dropping explicit linkers"
        )
    return np.asarray(mem_idx, np.int64), cortex_idx


def build_nucleus_compartment(
    pos_actin: np.ndarray, node_off: int, device: str, *, params: dict | None = None, centre=(0.0, 0.0, 0.0),
) -> tuple[NucleusCompartment, np.ndarray, MeshNucleusMaskProvider]:
    """Build the resting deformable nucleus: mesh + force-family device arrays + LINC pairs + mask provider.

    Args:
        pos_actin: (N_actin,3) cortex node positions [µm] (the LINC anchor candidates).
        node_off: global index of the first nucleus envelope node (nucleus nodes are appended after it).
        device: CUDA device string.
        params: overrides for :data:`NUCLEUS_DEFAULTS`.
        centre: nucleus centre ``centre_nuc`` [µm].

    Returns:
        ``(compartment, verts, mask_provider)`` — ``verts`` (Nv,3) to append to the global node array.
    """
    p = {**NUCLEUS_DEFAULTS, **(params or {})}
    lam = LaminaParams(k_chrom=p["k_chrom"], k_lamin_b=p["k_lamin_b"], k_lamin_ac=p["k_lamin_ac"],
                       knee_strain=p["knee_strain"], eps_rupture=p["eps_rupture"], kappa_ne=p["kappa_ne"])
    mesh = build_nucleus(p["r_eq_um"], lam, aspect=p["aspect"], subdivisions=p["subdivisions"], centre=centre)
    n_verts = int(mesh.verts.shape[0])

    faces_g = (mesh.faces + node_off).astype(np.int32)
    hinges_g = (mesh.hinges + node_off).astype(np.int32)

    # LINC pairs: each nucleus surface node → nearest cortex actin node within reach (force-free rest length).
    # cap→LINC→nucleus is I7; here LINC is the baseline nucleus↔cortex coupling (the single load path, §3).
    nuc_local, cortex_idx = _pair_to_nearest(mesh.verts, pos_actin, p["linc_reach_um"])
    if nuc_local.size:
        linc_n = (nuc_local + node_off).astype(np.int32)
        linc_a = cortex_idx.astype(np.int32)
        linc_rest = np.linalg.norm(mesh.verts[nuc_local] - pos_actin[cortex_idx], axis=1).astype(np.float64)
    else:
        linc_n = np.zeros(0, np.int32)
        linc_a = np.zeros(0, np.int32)
        linc_rest = np.zeros(0, np.float64)

    with wp.ScopedDevice(device):
        comp = NucleusCompartment(
            device=device, node_off=node_off, n_verts=n_verts, mesh=mesh,
            faces_d=wp.array(faces_g, dtype=wp.int32),
            hinges_d=wp.array(hinges_g, dtype=wp.int32),
            a0_d=wp.array(np.ascontiguousarray(mesh.A0_face, np.float64), dtype=wp.float64),
            ruptured_d=wp.zeros(int(mesh.faces.shape[0]), dtype=wp.int32),
            linc_n_d=wp.array(linc_n, dtype=wp.int32),
            linc_a_d=wp.array(linc_a, dtype=wp.int32),
            linc_rest_d=wp.array(linc_rest, dtype=wp.float64),
            vol_d=wp.zeros(1, dtype=wp.float64),
            p_vol_d=wp.zeros(1, dtype=wp.float64),
            kappa_tilde=float(mesh.kappa_tilde), v0_ref=float(mesh.V0), k_vol=p["k_vol"],
            k_soft=lam.k_soft, k_ac=lam.k_lamin_ac, knee=lam.knee_strain, eps_rupt=lam.eps_rupture,
            k_linc=p["k_linc"], linc_stiffening=p["linc_stiffening"],
            mean_edge_um=mesh_mean_edge(mesh.verts, mesh.faces))
    provider = MeshNucleusMaskProvider(mesh.verts, mesh.faces)
    return comp, mesh.verts, provider


# ─────────────────────────────────────────────────────────────────────────────────────────────────
# Membrane compartment
# ─────────────────────────────────────────────────────────────────────────────────────────────────
@dataclass
class MembraneCompartment:
    """Plasma-membrane Helfrich sheet force primitive — ``accumulate(pos, f)`` (§1.4).

    Envelope nodes live in ``pos[node_off : node_off+n_verts]``. Force families: reused
    ``helfrich_bending_kernel`` + in-plane ``membrane_area_kernel`` (γ_mem plateau; K_A upturn OFF) +
    commit-safe ERM tethers (membrane↔cortex; force-free at rest — the membrane rides on the cortex).
    """
    device: str
    node_off: int
    n_verts: int
    mesh: MembraneMesh          # immutable host topology/reference geometry (not live per-step state)
    faces_d: wp.array          # (Nf,3) int32 GLOBAL
    hinges_d: wp.array         # (Nh,4) int32 GLOBAL
    erm_m_d: wp.array          # (Ne,) int32  membrane node (global)
    erm_c_d: wp.array          # (Ne,) int32  cortex node (global)
    erm_bound_d: wp.array      # (Ne,) int32
    erm_rest_d: wp.array       # (Ne,) float64
    kappa_tilde: float
    gamma_mem: float
    k_erm: float
    f_rupt: float
    erm_density_target_per_um2: float | None
    erm_density_source: str
    erm_pairing_mode: str
    erm_bell_kinetics: ERMBellKinetics | None
    erm_detach_events_d: wp.array  # (1,) int32 cumulative accepted Bell detachments
    erm_attach_events_d: wp.array  # (1,) int32 cumulative accepted Bell attachments
    erm_rng_epoch_d: wp.array  # (1,) int32; advances only on accepted physical transactions
    with_area_tension: bool
    A0: float
    mean_edge_um: float
    n_faces: int = field(init=False)
    n_hinges: int = field(init=False)
    n_erm: int = field(init=False)

    def __post_init__(self):
        self.n_faces = int(self.faces_d.shape[0])
        self.n_hinges = int(self.hinges_d.shape[0])
        self.n_erm = int(self.erm_m_d.shape[0])

    def accumulate(self, pos: wp.array, f: wp.array) -> None:
        d = self.device
        if self.n_hinges:
            wp.launch(helfrich_bending_kernel, dim=self.n_hinges,
                      inputs=[pos, self.hinges_d, wp.float64(self.kappa_tilde)], outputs=[f], device=d)
        if self.with_area_tension and self.n_faces:      # γ_mem plateau (K_A reservoir upturn OFF, hard-truth #8)
            wp.launch(membrane_area_kernel, dim=self.n_faces,
                      inputs=[pos, self.faces_d, wp.float64(self.gamma_mem)], outputs=[f], device=d)
        if self.n_erm:                                    # membrane↔cortex ERM (force-free at rest). npos=cpos=pos
            if self.erm_bell_kinetics is None:
                wp.launch(erm_tether_force_kernel, dim=self.n_erm,
                          inputs=[pos, self.erm_m_d, self.erm_c_d, self.erm_bound_d,
                                  wp.float64(self.k_erm), self.erm_rest_d, wp.float64(self.f_rupt)],
                          outputs=[f], device=d)
            else:
                wp.launch(erm_bell_force_kernel, dim=self.n_erm,
                          inputs=[pos, self.erm_m_d, self.erm_c_d, self.erm_bound_d,
                                  wp.float64(self.k_erm), self.erm_rest_d],
                          outputs=[f], device=d)

    def commit_kinetics(self, pos: wp.array, accepted: wp.array, tau: float, rng_seed: int) -> None:
        """Commit one diagnostic rupture or Bell on/off tick at an accepted outer physical-time boundary."""
        if self.n_erm:
            if self.erm_bell_kinetics is None:
                wp.launch(erm_tether_rupture_kernel, dim=self.n_erm,
                          inputs=[pos, self.erm_m_d, self.erm_c_d, self.erm_bound_d,
                                  wp.float64(self.k_erm), self.erm_rest_d, wp.float64(self.f_rupt), accepted],
                          device=self.device)
            else:
                kinetics = self.erm_bell_kinetics
                wp.launch(
                    erm_bell_kmc_kernel,
                    dim=self.n_erm,
                    inputs=[
                        pos, self.erm_m_d, self.erm_c_d, self.erm_bound_d, wp.float64(self.k_erm),
                        self.erm_rest_d, wp.float64(kinetics.k_on_s), wp.float64(kinetics.k_off0_s),
                        wp.float64(kinetics.bell_force_pn), wp.float64(kinetics.capture_radius_um),
                        wp.float64(tau), wp.int32(rng_seed & 0x7FFFFFFF), self.erm_rng_epoch_d, accepted,
                        self.erm_detach_events_d, self.erm_attach_events_d,
                    ],
                    device=self.device,
                )
                wp.launch(
                    increment_erm_epoch_if_accepted_kernel,
                    dim=1,
                    inputs=[accepted, self.erm_rng_epoch_d],
                    device=self.device,
                )

    def rupture_step(self, pos: wp.array, accepted: wp.array) -> None:
        """Compatibility entry point for the diagnostic hard-threshold gate only.

        Production Bell kinetics require a physical ``tau`` and must use :meth:`commit_kinetics`; silently
        treating a Bell population as instantaneous rupture would reintroduce the abstraction this split
        removes.
        """
        if self.erm_bell_kinetics is not None:
            raise RuntimeError("Bell-kinetic ERM must be committed with commit_kinetics(pos, accepted, tau, seed)")
        self.commit_kinetics(pos, accepted, 0.0, 0)


def build_membrane_compartment(
    pos_actin: np.ndarray, node_off: int, device: str, *, R_mem: float, params: dict | None = None,
    centre=(0.0, 0.0, 0.0),
) -> tuple[MembraneCompartment, np.ndarray]:
    """Build the resting plasma-membrane Helfrich sheet + ERM pairs to the cortex; return (compartment, verts)."""
    p = {**MEMBRANE_DEFAULTS, **(params or {})}
    mm = build_membrane_mesh(R_mem, subdivisions=p["subdivisions"], centre=centre)
    hinges = build_membrane_hinges(mm.faces)
    kt = membrane_calibrate_kappa_tilde(p["kappa_m"], mm.verts, hinges)
    f_rupt = erm_rupture_force(p["kappa_m"], p["gamma_mem"])
    n_verts = int(mm.verts.shape[0])

    faces_g = (mm.faces + node_off).astype(np.int32)
    hinges_g = (hinges + node_off).astype(np.int32)

    density = p["erm_density_per_um2"]
    if density is None:
        if p.get("erm_radial_pairing", False):
            mem_local, cortex_idx = _pair_radial(mm.verts, pos_actin, np.asarray(centre), p["erm_reach_um"])
            pairing_mode = "RADIAL_ONE_PER_MEMBRANE_NODE_DIAGNOSTIC"
        else:
            mem_local, cortex_idx = _pair_to_nearest(mm.verts, pos_actin, p["erm_reach_um"])
            pairing_mode = "UNSOURCED_ONE_PER_MEMBRANE_NODE_DIAGNOSTIC"
    else:
        mem_local, cortex_idx = density_resolved_erm_pairs(
            mm.verts,
            pos_actin,
            # Count against the represented mesh area rather than the analytic
            # sphere area.  This keeps the molecular areal density invariant
            # when the continuum membrane resolution changes.
            surface_area_um2=triangle_surface_area(mm.verts, mm.faces),
            density_per_um2=float(density),
            reach_um=float(p["erm_reach_um"]),
        )
        pairing_mode = "EXPLICIT_DENSITY_RESOLVED"
    if mem_local.size:
        erm_m = (mem_local + node_off).astype(np.int32)
        erm_c = cortex_idx.astype(np.int32)
        erm_rest = np.linalg.norm(mm.verts[mem_local] - pos_actin[cortex_idx], axis=1).astype(np.float64)
    else:
        erm_m = np.zeros(0, np.int32)
        erm_c = np.zeros(0, np.int32)
        erm_rest = np.zeros(0, np.float64)

    with wp.ScopedDevice(device):
        comp = MembraneCompartment(
            device=device, node_off=node_off, n_verts=n_verts, mesh=mm,
            faces_d=wp.array(faces_g, dtype=wp.int32),
            hinges_d=wp.array(hinges_g, dtype=wp.int32),
            erm_m_d=wp.array(erm_m, dtype=wp.int32),
            erm_c_d=wp.array(erm_c, dtype=wp.int32),
            erm_bound_d=wp.ones(erm_m.shape[0], dtype=wp.int32),
            erm_rest_d=wp.array(erm_rest, dtype=wp.float64),
            kappa_tilde=float(kt), gamma_mem=p["gamma_mem"], k_erm=p["k_erm"], f_rupt=float(f_rupt),
            erm_density_target_per_um2=None if density is None else float(density),
            erm_density_source=str(p["erm_density_source"]), erm_pairing_mode=pairing_mode,
            erm_bell_kinetics=p["erm_bell_kinetics"],
            erm_detach_events_d=wp.zeros(1, dtype=wp.int32),
            erm_attach_events_d=wp.zeros(1, dtype=wp.int32),
            erm_rng_epoch_d=wp.zeros(1, dtype=wp.int32),
            with_area_tension=bool(p["with_area_tension"]), A0=float(mm.A0),
            mean_edge_um=mesh_mean_edge(mm.verts, mm.faces))
    return comp, mm.verts
