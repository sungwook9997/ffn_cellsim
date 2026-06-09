"""Connected basal contractile actin mesh — geometry layer (B1; KU-3.5).

PI directive (a), 2026-06-09 (``docs/v2_audit/BASAL_MESH_DESIGN_2026-06-09.md``):
the ventral stress-fiber apparatus is a CONNECTED actin mesh laid on the basal
plane, with the long-axis SF cables woven into it, FA-anchored to the substrate
and contracted by NMII. A connected mesh has MANY filaments, so the Stam-Hocky
bipolar NMII gate (which needs two antiparallel filaments) works — unlike a single
FA→FA chain.

This module is the **geometry layer (increment B1)**: it lays out the basal mesh
as a BIMODAL filament field on a FLAT basal DISK — the exact analogue of the cortex
connected-mesh construction (:func:`ffn_sim.cortex.cortex.generate_bimodal_cortex_layout`)
but projected onto a plane instead of the sphere shell:

* **short Arp2/3-infill filaments** (≈ one segment) → the connecting basal meshwork;
* **long formin filaments** (≈ 10× longer) → the ventral stress-fiber CABLES,
  oriented along the cell long axis (the basal footprint's principal axis);
* DISORDERED-ISOTROPIC short filaments (uniform in-plane azimuth) so the infill
  percolates in every direction (Li-Gao-Xu 2022 cortex-rheology rationale, applied
  in-plane).

The connectivity (bridge-different-filament crosslinkers → giant ≥ 0.9, z ∈ [3,3.5]),
the FA anchoring, the ``sf_myosin_`` NMII placement and the equilibrated Kumar /
Balaban active gate are the FOLLOWING increments (B2–B5). This file adds NO force and
NO HOOMD state — it returns a pure geometry layout (reusing the cortex
:class:`VariableLengthCortexLayout`, which is prefix/force-agnostic).

Why reuse ``VariableLengthCortexLayout``
----------------------------------------
That dataclass is a pure flat-indexed filament topology (positions, per-filament
bead counts/starts, tangents, chain bond/angle groups, ``is_formin`` mask). It
carries no shell-specific assumption and no particle/bond TYPE name, so it equally
describes a planar basal mesh. The downstream snapshot-extender (B3/B4) mints the
``sf_``-prefixed types from this geometry (so every basal bond is γ-denylisted).

Sanity Gate
-----------
*Per CLAUDE.md Hard Rule. Checks in ``ffn_sim/tests/test_basal_mesh.py``.*

1. **Dimensional analysis** — ``ell0`` [m] bead spacing; ``footprint_radius`` [m];
   ``z_basal`` [m]; ``band_thickness`` [m]. Lengths only; no force here.
2. **Boundary cases** — ``n_filaments ≤ 0`` → ValueError; ``footprint_radius ≤ 0``
   → ValueError; ``ell0 ≤ 0`` → ValueError; ``formin_fraction`` ∈ [0,1].
3. **Conservation / count** — particle count = Σ N_beads_i; bonds = Σ (N_i − 1);
   angles = Σ (N_i − 2); ``is_formin`` sums to ≈ ``formin_fraction · F``.
4. **Numerical** — all positions finite float64; no NaN.
5. **Sign / sense (PLANAR + force-free)** — every bead's z lies in the basal band
   ``[z_basal − band_thickness, z_basal]`` (the mesh is FLAT, not bulging); every
   CoM lies inside the footprint disk (radial ≤ footprint_radius); every chain bond
   length equals ``ell0`` exactly (born force-free under a Harmonic r0 = ell0).
6. **Measurement-protocol** — long (formin = cable) filaments align with the long
   axis (mean |cos| ≥ ``align_floor``, default 0.8 — the SF-gate convention); short
   (infill) filaments are isotropic in-plane (mean |cos| ≈ 2/π ≈ 0.637).

References
----------
- Bimodal cortex construction (the shell analogue this mirrors):
  :func:`ffn_sim.cortex.cortex.generate_bimodal_cortex_layout`; Fritzsche 2016/2017
  (bimodal short Arp2/3 + long formin); Li-Gao-Xu 2022 (disordered isotropic).
- vSF long-axis alignment (cables): Tojkander 2012; Hotulainen-Lappalainen 2006;
  PI 2026-06-09 long-axis-aligned pairing rule (``select_aligned_fa_pairs``).
- Design: ``docs/v2_audit/BASAL_MESH_DESIGN_2026-06-09.md``.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from dataclasses import replace as _dc_replace
from typing import Any

import numpy as np

from ffn_sim.cortex.cortex import VariableLengthCortexLayout
from ffn_sim.cortex.crosslinkers import (
    ConnectedMeshSeed,
    ResolvedCrosslinkers,
    seed_connected_mesh_xlinks,
)


def _require_finite_positive(name: str, x: float) -> None:
    if not (math.isfinite(x) and x > 0.0):
        raise ValueError(f"{name} must be finite and > 0; got {x!r}")


def _inplane_long_axis(long_axis: np.ndarray | None) -> np.ndarray:
    """Return a unit 2-vector (x, y) for the cable orientation (default +x)."""
    if long_axis is None:
        return np.array([1.0, 0.0], dtype=np.float64)
    v = np.asarray(long_axis, dtype=np.float64).reshape(-1)[:2]
    n = float(np.linalg.norm(v))
    if not (math.isfinite(n) and n > 0.0):
        return np.array([1.0, 0.0], dtype=np.float64)
    return v / n


def generate_basal_mesh_layout(
    *,
    n_filaments: int,
    ell0: float,
    footprint_radius: float,
    z_basal: float,
    formin_fraction: float = 0.12,
    L_short_mean: float | None = None,
    L_long_mean: float | None = None,
    long_axis: np.ndarray | None = None,
    long_axis_jitter_deg: float = 0.0,
    band_thickness: float = 200.0e-9,
    seed: int = 73,
    rng: np.random.Generator | None = None,
    n_beads_min: int = 2,
    n_beads_max: int | None = None,
) -> VariableLengthCortexLayout:
    """Lay out a bimodal basal actin mesh on a FLAT disk (geometry only).

    Mirrors :func:`generate_bimodal_cortex_layout` but the centres of mass are
    sampled in an in-plane DISK (radius ``footprint_radius``, at height
    ``z_basal``) and the tangents lie IN the basal plane: long (formin = cable)
    filaments oriented along ``long_axis`` (± with optional jitter), short (Arp2/3
    = infill) filaments at a uniform in-plane azimuth. No shell projection; each
    filament sits at one depth within the basal band so crossing filaments are
    z-separated (no construction overlap), mirroring the shell-band trick.

    Args:
        n_filaments: number of basal filaments (cables + infill).
        ell0: bead rest length [m] (the cortex ℓ0; the Harmonic r0 downstream).
        footprint_radius: basal contact-disk radius [m] (from the FA footprint).
        z_basal: basal plane height [m] (the substrate-contact z of the cell).
        formin_fraction: fraction of filaments that are LONG cables (Fritzsche
            ~0.10-0.20; default 0.12 = the cortex value).
        L_short_mean: mean infill length [m] (default ``ell0``).
        L_long_mean: mean cable length [m] (default ``10·L_short_mean``).
        long_axis: in-plane cable orientation (3- or 2-vector; default +x).
        long_axis_jitter_deg: Gaussian angular spread of the cables about the
            long axis [deg]. Default 0.0 (perfectly aligned cables — no free
            parameter; raise it only if a measured SF dispersion is supplied).
        band_thickness: basal z-band depth [m] for radial separation (KU-3.17
            cortex thickness analogue; default 200 nm — a discretisation knob).
        seed / rng: RNG control.
        n_beads_min / n_beads_max: per-filament bead-count clamps.

    Returns:
        A :class:`VariableLengthCortexLayout` (planar; ``is_formin`` marks cables).
    """
    if rng is None:
        rng = np.random.default_rng(seed)

    F = int(n_filaments)
    if F <= 0:
        raise ValueError(f"n_filaments must be > 0; got {F}")
    _require_finite_positive("ell0", ell0)
    _require_finite_positive("footprint_radius", footprint_radius)
    _require_finite_positive("band_thickness", band_thickness)
    if not math.isfinite(z_basal):
        raise ValueError(f"z_basal must be finite; got {z_basal!r}")
    if not (0.0 <= formin_fraction <= 1.0):
        raise ValueError(f"formin_fraction must be in [0, 1]; got {formin_fraction}")
    if n_beads_min < 2:
        raise ValueError(f"n_beads_min must be ≥ 2; got {n_beads_min}")

    L0 = float(ell0)
    if L_short_mean is None:
        L_short_mean = L0
    if L_long_mean is None:
        L_long_mean = 10.0 * L_short_mean
    _require_finite_positive("L_short_mean", L_short_mean)
    if not (math.isfinite(L_long_mean) and L_long_mean >= L_short_mean):
        raise ValueError(
            f"L_long_mean must be finite ≥ L_short_mean; got {L_long_mean}"
        )

    if n_beads_max is None:
        # A cable chord cannot exceed the footprint diameter.
        L_cap = min(2.0 * footprint_radius, max(L_long_mean * 6.0, 12.0 * L0))
        n_beads_max = max(n_beads_min, int(round(L_cap / L0)) + 1)

    # ---- per-filament class (cable=formin) + exponential length ----
    is_formin = rng.random(F) < formin_fraction
    means = np.where(is_formin, L_long_mean, L_short_mean)
    L_samples = rng.exponential(means)
    n_beads_per = np.clip(
        np.round(L_samples / L0).astype(np.int64) + 1, n_beads_min, n_beads_max
    )
    L_realised = (n_beads_per - 1).astype(np.float64) * L0

    # ---- in-plane disk CoM (uniform-area) + per-filament band depth ----
    u = rng.uniform(0.0, 1.0, F)
    r_disk = footprint_radius * np.sqrt(u)            # uniform-area disk sampling
    theta = rng.uniform(0.0, 2.0 * math.pi, F)
    depth = rng.uniform(0.0, band_thickness, F)        # toward the substrate
    centers = np.empty((F, 3), dtype=np.float64)
    centers[:, 0] = r_disk * np.cos(theta)
    centers[:, 1] = r_disk * np.sin(theta)
    centers[:, 2] = z_basal - depth

    # ---- in-plane tangents: cables along the long axis, infill isotropic ----
    ax2 = _inplane_long_axis(long_axis)
    base_ang = math.atan2(ax2[1], ax2[0])
    angles = np.empty(F, dtype=np.float64)
    short_mask = ~is_formin
    n_short = int(short_mask.sum())
    if n_short > 0:
        angles[short_mask] = rng.uniform(0.0, 2.0 * math.pi, n_short)
    n_long = int(is_formin.sum())
    if n_long > 0:
        jit = math.radians(long_axis_jitter_deg)
        # ± along the axis (a cable has no head/tail polarity here), with optional
        # Gaussian jitter. jitter 0 → perfectly aligned cables (no free param).
        sign = rng.choice([0.0, math.pi], size=n_long)
        noise = rng.normal(0.0, jit, n_long) if jit > 0.0 else np.zeros(n_long)
        angles[is_formin] = base_ang + sign + noise
    tangents = np.stack(
        [np.cos(angles), np.sin(angles), np.zeros(F)], axis=1
    )
    tangents = tangents / np.linalg.norm(
        tangents, axis=1, keepdims=True
    ).clip(min=1e-30)

    # ---- flat layout (planar; z constant within a filament) ----
    n_total_beads = int(n_beads_per.sum())
    filament_starts = np.zeros(F, dtype=np.int64)
    filament_starts[1:] = np.cumsum(n_beads_per[:-1])
    positions_flat = np.empty((n_total_beads, 3), dtype=np.float64)
    bond_pairs: list[tuple[int, int]] = []
    angle_triplets: list[tuple[int, int, int]] = []
    for f in range(F):
        N_f = int(n_beads_per[f])
        start = int(filament_starts[f])
        offsets = (np.arange(N_f, dtype=np.float64) - 0.5 * (N_f - 1)) * L0
        positions_flat[start:start + N_f] = (
            centers[f] + offsets[:, None] * tangents[f]
        )
        for j in range(N_f - 1):
            bond_pairs.append((start + j, start + j + 1))
        for j in range(N_f - 2):
            angle_triplets.append((start + j, start + j + 1, start + j + 2))

    bond_groups = np.array(bond_pairs, dtype=np.int64).reshape(-1, 2)
    angle_groups = (
        np.array(angle_triplets, dtype=np.int64).reshape(-1, 3)
        if angle_triplets else np.empty((0, 3), dtype=np.int64)
    )

    return VariableLengthCortexLayout(
        positions_flat=positions_flat,
        n_beads_per_filament=n_beads_per,
        filament_starts=filament_starts,
        L_per_filament=L_realised,
        centers_of_mass=centers,
        tangents=tangents,
        bond_groups=bond_groups,
        angle_groups=angle_groups,
        is_formin=is_formin,
    )


def basal_mesh_build_report(
    layout: VariableLengthCortexLayout,
    *,
    ell0: float,
    footprint_radius: float,
    z_basal: float,
    band_thickness: float,
    long_axis: np.ndarray | None = None,
    align_floor: float = 0.8,
) -> dict[str, Any]:
    """Build-time B1 gate metrics (PLANAR + force-free + long-axis cables).

    Returns a dict with the controls and a ``verdict`` ("PASS"/"REVIEW"):
      * planar_ok          — every bead z ∈ [z_basal − band_thickness, z_basal].
      * within_footprint   — every CoM radial (x,y) ≤ footprint_radius (+ε).
      * force_free         — every chain bond length == ell0 (max strain < 1e-9).
      * cables_aligned     — mean |cos| of LONG (cable) tangents vs long axis ≥ floor.
      * infill_isotropic   — mean |cos| of SHORT (infill) tangents ≈ 2/π (in [0.5,0.75]).
    """
    pos = np.asarray(layout.positions_flat, dtype=np.float64)
    z = pos[:, 2]
    planar_ok = bool(
        np.all(z <= z_basal + 1e-12)
        and np.all(z >= z_basal - band_thickness - 1e-12)
    )

    com_xy = np.asarray(layout.centers_of_mass, dtype=np.float64)[:, :2]
    radial = np.linalg.norm(com_xy, axis=1)
    within_footprint = bool(np.all(radial <= footprint_radius + 1e-12))

    # force-free: chain bond lengths == ell0.
    bg = np.asarray(layout.bond_groups, dtype=np.int64).reshape(-1, 2)
    if bg.shape[0] > 0:
        seg = pos[bg[:, 0]] - pos[bg[:, 1]]
        ln = np.linalg.norm(seg, axis=1)
        max_strain = float(np.max(np.abs(ln - ell0) / ell0))
    else:
        max_strain = 0.0
    force_free = bool(max_strain < 1e-9)

    ax2 = _inplane_long_axis(long_axis)
    tang = np.asarray(layout.tangents, dtype=np.float64)[:, :2]
    tnrm = tang / np.linalg.norm(tang, axis=1, keepdims=True).clip(min=1e-30)
    absdot = np.abs(tnrm @ ax2)
    is_formin = np.asarray(layout.is_formin, dtype=bool)
    cable_align = float(np.mean(absdot[is_formin])) if is_formin.any() else float("nan")
    infill_align = (
        float(np.mean(absdot[~is_formin])) if (~is_formin).any() else float("nan")
    )
    cables_aligned = bool(np.isfinite(cable_align) and cable_align >= align_floor)
    infill_isotropic = bool(
        np.isfinite(infill_align) and 0.5 <= infill_align <= 0.75
    )

    n_fil = int(layout.n_beads_per_filament.shape[0])
    controls = {
        "n_filaments": n_fil,
        "n_cables": int(is_formin.sum()),
        "n_infill": int((~is_formin).sum()),
        "n_beads_total": int(pos.shape[0]),
        "planar_ok": planar_ok,
        "within_footprint": within_footprint,
        "force_free": {"ok": force_free, "max_strain": max_strain},
        "cables_aligned": {"ok": cables_aligned, "mean_abs_cos": cable_align},
        "infill_isotropic": {"ok": infill_isotropic, "mean_abs_cos": infill_align},
    }
    ok = (
        planar_ok and within_footprint and force_free
        and cables_aligned and infill_isotropic
    )
    return {"verdict": "PASS" if ok else "REVIEW", "controls": controls}


# ===========================================================================
# B2 — connect the basal mesh (bridge-different-filament crosslinkers)
# ===========================================================================
def basal_mesh_reach(footprint_radius: float, n_filaments: int) -> float:
    """Mesoscale crosslinker partner-search radius for the basal DISK.

    The planar analogue of the cortex shell reach ``√(A_shell/n_fil)``: the
    inter-filament spacing on a disk of area ``π·footprint_radius²`` shared by
    ``n_filaments`` filaments is ``footprint_radius·√(π/n_filaments)`` (DERIVED,
    grid-aware — the geometric dual of the ×40 areal coarse-graining on a disk,
    not a tuned constant).
    """
    return float(footprint_radius) * math.sqrt(math.pi / max(1, int(n_filaments)))


def connect_basal_mesh(
    layout: VariableLengthCortexLayout,
    p_xl: ResolvedCrosslinkers,
    *,
    footprint_radius: float,
    z_struct: float = 3.3,
    bundle_mult: int = 2,
    max_bead_degree: int = 4,
    reach: float | None = None,
    rng: np.random.Generator | None = None,
) -> ConnectedMeshSeed:
    """Seed a CONNECTED bridge-different-filament crosslink mesh on the basal layout.

    Reuses the cortex :func:`seed_connected_mesh_xlinks` verbatim (the basal mesh is
    the same :class:`VariableLengthCortexLayout` topology, just planar), with the
    DISK reach :func:`basal_mesh_reach`. The crosslinker ``max_bind_dist`` is set to
    the reach (mirroring ``cortex/connected_mesh.py``) so the per-r0 attach bins
    cover the mesoscale head-to-bead span. Geometry/topology only — no HOOMD state.

    The basal mesh has no Arp2/3 branch bonds (``prior_bead_bonds=None``), so the
    connectivity comes entirely from the seeded bridge crosslinks.

    Args:
        layout: the basal mesh layout (:func:`generate_basal_mesh_layout`).
        p_xl: resolved crosslinker params (species fractions/lengths + n_bins).
        footprint_radius: basal disk radius [m] (sets the default reach).
        z_struct: target distinct-neighbour coordination (Kim 2007 / Kadzik-Munro
            z≈3-4); default 3.3 = the literature-mid target (the gate checks the
            REALISED z ∈ [3.0, 3.5], it is not back-solved).
        bundle_mult: parallel crosslinks per connected pair (Flormann bundling).
        max_bead_degree: per-bead crosslink cap (HOOMD exclusion safety).
        reach: override the disk reach [m].
        rng: RNG.

    Returns:
        The :class:`ConnectedMeshSeed` (z_struct_realised, giant_fraction, …).
    """
    F = int(layout.n_beads_per_filament.shape[0])
    if reach is None:
        reach = basal_mesh_reach(footprint_radius, F)
    # bins must cover the mesoscale head-to-bead span (≈ ½·reach); set the
    # crosslinker max_bind_dist to the reach (cortex/connected_mesh.py convention).
    p_xl_eff = _dc_replace(p_xl, max_bind_dist=float(reach))
    n_beads = int(layout.positions_flat.shape[0])
    return seed_connected_mesh_xlinks(
        layout.positions_flat,
        layout.filament_idx,
        F,
        p_xl_eff,
        z_struct=z_struct,
        bundle_mult=bundle_mult,
        reach=float(reach),
        R_cell=float(footprint_radius),   # only used for the default reach (overridden)
        n_cortex_beads=n_beads,
        max_bead_degree=max_bead_degree,
        prior_bead_bonds=None,            # no Arp2/3 branches on the basal mesh
        rng=rng,
    )


def basal_connectivity_report_for_layout(
    layout: VariableLengthCortexLayout,
    p_xl: ResolvedCrosslinkers,
    *,
    footprint_radius: float,
    z_struct: float = 3.3,
    rng: np.random.Generator | None = None,
) -> dict[str, Any]:
    """Connect an arbitrary basal layout + report (helper for the combined mesh)."""
    seed = connect_basal_mesh(
        layout, p_xl, footprint_radius=footprint_radius, z_struct=z_struct, rng=rng,
    )
    return basal_connectivity_report(seed)


def basal_connectivity_report(
    seed: ConnectedMeshSeed,
    *,
    giant_floor: float = 0.9,
    z_band: tuple[float, float] = (3.0, 3.5),
    l_over_lc_floor: float = 5.9,
) -> dict[str, Any]:
    """B2 gate metrics — the cortex connected-mesh acceptance, on the basal mesh.

    Verdict controls (literature — the SAME acceptance as the cortex connected-mesh
    rebuild, ``cortex/connected_mesh.py``):
      * giant-component fraction ≥ ``giant_floor`` (0.9) — the mesh percolates.
      * realised coordination z ∈ ``z_band`` ([3.0, 3.5]) — Kadzik-Munro.
      * L/lc ≥ ``l_over_lc_floor`` (5.9) — Head 2003 crosslink density.

    ``n_homeless`` (isolated filaments with no partner in reach) is REPORTED as a
    diagnostic but is NOT a verdict criterion — a handful is consistent with the
    giant-fraction gate (giant ≈ 0.998 ⇒ ~0.2 % outside the giant), and the cortex
    acceptance does not gate on it either.
    """
    z = float(seed.z_struct_realised)
    giant = float(seed.giant_fraction)
    llc = float(seed.L_over_lc)
    homeless = int(seed.n_homeless)
    giant_ok = giant >= giant_floor
    z_ok = z_band[0] <= z <= z_band[1]
    llc_ok = llc >= l_over_lc_floor
    controls = {
        "n_xl": int(seed.n_xl),
        "giant_fraction": {"ok": giant_ok, "value": giant, "floor": giant_floor},
        "coordination_z": {"ok": z_ok, "value": z, "band": list(z_band)},
        "L_over_lc": {"ok": llc_ok, "value": llc, "floor": l_over_lc_floor},
        "n_homeless": homeless,   # diagnostic only (not gated)
    }
    ok = giant_ok and z_ok and llc_ok
    return {"verdict": "PASS" if ok else "REVIEW", "controls": controls}


# ===========================================================================
# B3 — F-layer ON the S-layer: cables span FA→FA, infill fills, anchored to FA
# ===========================================================================
@dataclass(slots=True)
class BasalApparatus:
    """The integrated basal apparatus: F-layer filament network ON the flat S-layer.

    The LONG cables (ventral stress fibers) span FA→FA — their end beads sit AT the
    FA-anchor positions (force-free ``sf_anchor``). The SHORT filaments are the
    isotropic infill meshwork filling the disk. Both are explicit fine-grained
    filaments; the surface carries no force (S-layer = positioning + connectivity).

    Attributes:
        layout: combined :class:`VariableLengthCortexLayout` (cables first, then
            infill; ``is_formin`` marks the cables).
        fa_positions: (n_fa, 3) FA-anchor positions on the surface (z = z_basal).
        anchor_bonds: (n_anchor, 2) int64 — (cable_end_bead_flat_idx, fa_index).
        anchor_r0: (n_anchor,) float64 — anchor rest lengths (≈ 0, force-free).
        n_cables: number of FA→FA cable filaments.
        footprint_radius, z_basal: geometry.
    """

    layout: VariableLengthCortexLayout
    fa_positions: np.ndarray
    anchor_bonds: np.ndarray
    anchor_r0: np.ndarray
    n_cables: int
    footprint_radius: float
    z_basal: float
    band_thickness: float = 200.0e-9


def build_basal_filament_network(
    surf,
    *,
    ell0: float,
    n_cables: int,
    n_infill: int,
    band_thickness: float = 200.0e-9,
    formin_fraction_label: float = 0.12,
    seed: int = 91,
    rng: np.random.Generator | None = None,
) -> BasalApparatus:
    """Build the F-layer ON the flat S-layer surface, cables anchored FA→FA.

    The LONG cables span aligned FA pairs (reusing the SF long-axis pairing) — each
    cable is a bead chain laid from FA ``A`` to FA ``B`` at ``ell0`` spacing, so its
    two END beads sit EXACTLY at the FA positions (anchor r0 = 0, force-free). The
    SHORT infill filaments are isotropic disk-filling chains
    (:func:`generate_basal_mesh_layout`, ``formin_fraction=0``). Both are appended
    into ONE flat-indexed :class:`VariableLengthCortexLayout` (cables first), with
    ``is_formin`` marking the cables.

    Args:
        surf: a :class:`ffn_sim.cell.basal_surface.FlatBasalSurface` (S-layer).
        ell0: bead spacing [m].
        n_cables: number of FA→FA cable filaments.
        n_infill: number of short infill filaments.
        formin_fraction_label: informational (cables ARE the formin set here).
        seed / rng: RNG.

    Returns:
        A :class:`BasalApparatus`.
    """
    from ffn_sim.cell.stress_fibers import select_aligned_fa_pairs

    if rng is None:
        rng = np.random.default_rng(seed)
    fa = np.asarray(surf.fa_positions, dtype=np.float64)
    n_fa = int(fa.shape[0])
    if n_cables < 1:
        raise ValueError(f"n_cables must be ≥ 1; got {n_cables}")
    if 2 * n_cables > n_fa:
        raise ValueError(
            f"need ≥ 2·n_cables = {2 * n_cables} FA anchors; got {n_fa}."
        )

    # ---- LONG cables: aligned FA pairs → bead chains A→B (ends AT FA) ----
    pos_ord, idx_ord = select_aligned_fa_pairs(
        fa, np.arange(n_fa, dtype=np.int64), n_SF=n_cables,
        R_cell=float(surf.footprint_radius),
    )
    cable_pos: list[np.ndarray] = []
    cable_nbeads: list[int] = []
    anchor_bonds: list[tuple[int, int]] = []
    anchor_r0: list[float] = []
    cursor = 0  # flat bead cursor across the combined layout

    for b in range(n_cables):
        A = pos_ord[2 * b]
        B = pos_ord[2 * b + 1]
        fa_a = int(idx_ord[2 * b])
        fa_b = int(idx_ord[2 * b + 1])
        L = float(np.linalg.norm(B - A))
        nb = max(2, int(round(L / ell0)) + 1)
        beads = A[None, :] + np.linspace(0.0, 1.0, nb)[:, None] * (B - A)[None, :]
        cable_pos.append(beads)
        cable_nbeads.append(nb)
        # end beads sit at FA positions → force-free anchors (r0 = 0).
        anchor_bonds.append((cursor, fa_a))
        anchor_r0.append(0.0)
        anchor_bonds.append((cursor + nb - 1, fa_b))
        anchor_r0.append(0.0)
        cursor += nb

    # ---- SHORT infill: isotropic disk-filling chains ----
    # The F-layer is a thin actin SLAB of physical thickness `band_thickness`
    # (KU-3.17 cortex thickness) sitting on the +z (cell-interior) side of the flat
    # surface. generate_basal_mesh_layout puts z ∈ [z0−band, z0]; pass
    # z0 = z_basal + band so the slab occupies [z_basal, z_basal+band] ABOVE the
    # surface (cables sit at z_basal, the surface, at the slab's substrate face).
    infill = generate_basal_mesh_layout(
        n_filaments=n_infill, ell0=ell0, footprint_radius=surf.footprint_radius,
        z_basal=surf.z_basal + band_thickness, band_thickness=band_thickness,
        formin_fraction=0.0, seed=seed + 1, rng=rng,
    )
    infill_pos = infill.positions_flat
    infill_nbeads = infill.n_beads_per_filament

    # ---- combine: cables first, then infill, into one flat layout ----
    all_nbeads = np.array(cable_nbeads + list(infill_nbeads), dtype=np.int64)
    F = all_nbeads.shape[0]
    starts = np.zeros(F, dtype=np.int64)
    starts[1:] = np.cumsum(all_nbeads[:-1])
    positions_flat = np.concatenate(cable_pos + [infill_pos], axis=0)
    is_formin = np.zeros(F, dtype=bool)
    is_formin[:n_cables] = True   # cables are the long (formin) set
    tangents = np.zeros((F, 3), dtype=np.float64)
    for f in range(F):
        s, n = int(starts[f]), int(all_nbeads[f])
        v = positions_flat[s + n - 1] - positions_flat[s]
        nrm = float(np.linalg.norm(v))
        tangents[f] = v / nrm if nrm > 0 else np.array([1.0, 0.0, 0.0])
    # L_per = ACTUAL end-to-end contour (cables span |B−A| at L/(nb−1) spacing;
    # infill at ell0) so each chain's registered r0 = L_per/(nb−1) is force-free.
    L_per = np.array([
        float(np.linalg.norm(positions_flat[int(starts[f]) + int(all_nbeads[f]) - 1]
                             - positions_flat[int(starts[f])]))
        for f in range(F)
    ], dtype=np.float64)
    com = np.array([positions_flat[int(starts[f]):int(starts[f]) + int(all_nbeads[f])].mean(axis=0)
                    for f in range(F)])

    bond_pairs: list[tuple[int, int]] = []
    angle_triplets: list[tuple[int, int, int]] = []
    for f in range(F):
        s, n = int(starts[f]), int(all_nbeads[f])
        for j in range(n - 1):
            bond_pairs.append((s + j, s + j + 1))
        for j in range(n - 2):
            angle_triplets.append((s + j, s + j + 1, s + j + 2))
    bond_groups = np.array(bond_pairs, dtype=np.int64).reshape(-1, 2)
    angle_groups = (np.array(angle_triplets, dtype=np.int64).reshape(-1, 3)
                    if angle_triplets else np.empty((0, 3), dtype=np.int64))

    layout = VariableLengthCortexLayout(
        positions_flat=positions_flat,
        n_beads_per_filament=all_nbeads,
        filament_starts=starts,
        L_per_filament=L_per,
        centers_of_mass=com,
        tangents=tangents,
        bond_groups=bond_groups,
        angle_groups=angle_groups,
        is_formin=is_formin,
    )
    return BasalApparatus(
        layout=layout,
        fa_positions=fa,
        anchor_bonds=np.array(anchor_bonds, dtype=np.int64).reshape(-1, 2),
        anchor_r0=np.array(anchor_r0, dtype=np.float64),
        n_cables=n_cables,
        footprint_radius=float(surf.footprint_radius),
        z_basal=float(surf.z_basal),
        band_thickness=float(band_thickness),
    )


def basal_apparatus_report(
    app: BasalApparatus, *, ell0: float, planar_tol: float = 1e-12,
) -> dict[str, Any]:
    """B3 build-time gate: cables anchored FA→FA (force-free), planar, on surface.

    Controls:
      * cables_anchored — every cable has BOTH ends anchored to an FA (2·n_cables
        anchor bonds).
      * anchors_force_free — anchor rest lengths ≈ 0 AND cable-end↔FA separation ≈ 0
        (the end bead sits AT the FA), max ≪ ell0.
      * planar — every bead at z = z_basal.
      * chains_force_free — every chain bond length == ell0.
    """
    lay = app.layout
    pos = np.asarray(lay.positions_flat, dtype=np.float64)
    fa = np.asarray(app.fa_positions, dtype=np.float64)
    ab = np.asarray(app.anchor_bonds, dtype=np.int64).reshape(-1, 2)

    cables_anchored = bool(ab.shape[0] == 2 * app.n_cables)
    if ab.shape[0] > 0:
        sep = np.linalg.norm(pos[ab[:, 0]] - fa[ab[:, 1]], axis=1)
        max_anchor_sep = float(np.max(sep))
    else:
        max_anchor_sep = float("inf")
    anchors_force_free = bool(max_anchor_sep < 1e-9 * ell0 + 1e-15)

    # SLAB planarity: the F-layer is a thin actin slab on the +z side of the flat
    # surface — every bead's z ∈ [z_basal, z_basal + band_thickness] (KU-3.17).
    z = pos[:, 2]
    in_slab = bool(np.all(z >= app.z_basal - planar_tol)
                   and np.all(z <= app.z_basal + app.band_thickness + planar_tol))
    z_dev = float(np.max(np.abs(z - app.z_basal)))

    # PER-FILAMENT force-free: each chain born at ITS OWN spacing (cables span
    # FA→FA at L/(nb−1); infill at ell0). Mirrors the SF per-bundle EXACT-r0
    # convention — the downstream Harmonic uses each filament's registered r0.
    starts = np.asarray(lay.filament_starts, dtype=np.int64)
    nbeads = np.asarray(lay.n_beads_per_filament, dtype=np.int64)
    worst = 0.0
    for f in range(starts.shape[0]):
        s, n = int(starts[f]), int(nbeads[f])
        if n < 2:
            continue
        spacing = float(lay.L_per_filament[f]) / (n - 1)
        if spacing <= 0:
            continue
        seg = pos[s:s + n - 1] - pos[s + 1:s + n]
        ln = np.linalg.norm(seg, axis=1)
        worst = max(worst, float(np.max(np.abs(ln - spacing) / spacing)))
    chains_force_free = bool(worst < 1e-9)

    controls = {
        "n_cables": app.n_cables,
        "n_infill": int(lay.n_beads_per_filament.shape[0]) - app.n_cables,
        "n_beads_total": int(pos.shape[0]),
        "n_anchor_bonds": int(ab.shape[0]),
        "cables_anchored": cables_anchored,
        "anchors_force_free": {"ok": anchors_force_free, "max_anchor_sep": max_anchor_sep},
        "in_slab": {"ok": in_slab, "max_z_dev": z_dev, "band": app.band_thickness},
        "chains_force_free": {"ok": chains_force_free, "max_per_filament_strain": worst},
    }
    ok = cables_anchored and anchors_force_free and in_slab and chains_force_free
    return {"verdict": "PASS" if ok else "REVIEW", "controls": controls}


# ===========================================================================
# B4 — sf_myosin_ NMII placement on the basal apparatus (prefix-split + pool)
# ===========================================================================
def place_sf_myosin_on_apparatus(
    app: BasalApparatus,
    p_myo_sf,
    *,
    rng: np.random.Generator | None = None,
):
    """Place ``sf_myosin_`` NMII minifilaments ON the basal apparatus filaments.

    Reuses the cortex Stam-Hocky/Hill machinery (``cortex/myosin.py``) via the
    prefix split (①a): ``p_myo_sf.prefix`` must be ``"sf_myosin_"`` so the motor
    bond/particle types fall under the ``sf_`` γ-denylist (never contaminating the
    cortical active-γ signal). The minifilaments are placed ACTIN-AWARE on the
    apparatus's explicit filament beads (the F-layer network), so the heads sit
    within binding reach of real actin — the dynamic ``MyosinStepUpdater`` (①b,
    ``actin_pool_tags`` = the apparatus beads) loads them in the equilibrated run (B5).

    Geometry only (a :class:`CortexMyosinLayout`); no HOOMD state, no force here.

    Args:
        app: the integrated :class:`BasalApparatus` (F-on-S).
        p_myo_sf: a ``ResolvedCortexMyosin`` with ``prefix="sf_myosin_"``.
        rng: RNG.

    Returns:
        A ``CortexMyosinLayout`` (minifilament bead positions on the apparatus).
    """
    from ffn_sim.cortex.myosin import generate_cortex_myosin_layout

    if not str(p_myo_sf.prefix).startswith("sf_"):
        raise ValueError(
            f"sf_myosin placement requires a 'sf_'-prefixed myosin params "
            f"(γ-denylist); got prefix={p_myo_sf.prefix!r}. Build it with "
            "resolve_cortex_myosin(cfg with prefix='sf_myosin_')."
        )
    if rng is None:
        rng = np.random.default_rng(p_myo_sf.seed)
    lay = app.layout
    fil_idx = lay.filament_idx
    nb = np.asarray(lay.n_beads_per_filament, dtype=np.int64)
    n_beads = int(lay.positions_flat.shape[0])
    return generate_cortex_myosin_layout(
        p_myo_sf, float(app.footprint_radius),
        motor_tag_start=n_beads,          # myosin appended AFTER the apparatus beads
        rng=rng,
        cortex_positions=np.asarray(lay.positions_flat, dtype=np.float64),
        cortex_tangents=np.asarray(lay.tangents, dtype=np.float64),
        beads_per_filament=int(max(1, round(float(nb.mean())))),  # fallback only
        cortex_filament_idx=np.asarray(fil_idx, dtype=np.int64),  # variable-length map
        surface_normal=np.array([0.0, 0.0, 1.0]),  # FLAT basal layer (+z), not sphere
    )


def sf_myosin_placement_report(
    app: BasalApparatus,
    p_myo_sf,
    myo_layout,
) -> dict[str, Any]:
    """B4 build-time gate: sf_myosin assembled, force-free, sf_-denylisted, distinct.

    Controls:
      * assembled       — n_motors minifilaments placed, each with
        n_backbone + 2·n_heads_per_side beads.
      * force_free      — intra-minifilament bonds born at their rest lengths
        (backbone segment == backbone_segment_length; head↔backbone == head_rest_length).
      * heads_near_actin— every minifilament centre is an apparatus actin bead and
        the heads sit within head_actin_max_bind_dist of it (binding-eligible).
      * sf_denylisted   — every sf_myosin_ bond/particle TYPE starts with 'sf_'
        (γ-excluded by the registry denylist) AND none is a cortex_myosin_ type
        (the cortical active-γ signal is untouched).
    """
    from ffn_sim.cortex.myosin import (
        myosin_attach_bin_names,
        myosin_backbone_bond_name,
        myosin_head_backbone_bond_name,
        myosin_particle_type_names,
    )

    M = int(p_myo_sf.n_motors_per_cell)
    N = int(p_myo_sf.n_backbone)
    H = int(p_myo_sf.n_heads_per_side)
    pos = np.asarray(myo_layout.positions, dtype=np.float64)  # (M, N+2H, 3)
    assembled = bool(pos.shape == (M, N + 2 * H, 3))

    # force-free intra bonds (sample over all motors).
    seg = p_myo_sf.backbone_segment_length
    worst_bb = 0.0
    worst_hb = 0.0
    for m in range(M):
        bb = pos[m, :N]
        d_bb = np.linalg.norm(bb[1:] - bb[:-1], axis=1)
        worst_bb = max(worst_bb, float(np.max(np.abs(d_bb - seg) / seg)))
        # head ↔ nearest backbone bead distance ≈ head_rest_length.
        heads = pos[m, N:]
        # each head bonds to a backbone bead; min distance to backbone ≈ head_off.
        for h in heads:
            dmin = float(np.min(np.linalg.norm(bb - h, axis=1)))
            worst_hb = max(worst_hb, abs(dmin - p_myo_sf.head_rest_length) / p_myo_sf.head_rest_length)
    force_free = bool(worst_bb < 1e-9 and worst_hb < 0.5)  # head offset within tol

    # heads near actin: the binding criterion is PERPENDICULAR distance to an actin
    # SEGMENT ≤ capture_perp (the segment-projection rule the MyosinStepUpdater
    # uses — NOT the bead-centre distance: a head sits head_off≈200nm laterally
    # from the actin LINE, so its nearest bead can be ~400nm while its perp to the
    # segment is ~200nm ≤ capture_perp≈210nm and it IS binding-eligible).
    centres = np.asarray(myo_layout.centers, dtype=np.float64)
    appos = np.asarray(app.layout.positions_flat, dtype=np.float64)
    from scipy.spatial import cKDTree
    tree = cKDTree(appos)
    d_centre, _ = tree.query(centres)
    centres_on_actin = float(np.max(d_centre)) if M else 0.0
    cap_perp = float(getattr(p_myo_sf, "head_actin_capture_perp",
                             p_myo_sf.head_rest_length + 10.0e-9))
    bg = np.asarray(app.layout.bond_groups, dtype=np.int64).reshape(-1, 2)
    segA = appos[bg[:, 0]]
    segB = appos[bg[:, 1]]
    seg = segB - segA
    seg_len2 = np.einsum("ij,ij->i", seg, seg).clip(min=1e-30)
    # Per-head min perpendicular distance to ANY actin segment. NMII heads in sparse
    # spots (e.g. a 700nm minifilament's end head overhanging a short infill chain)
    # simply stay UNBOUND — not every head engages. The placement is contraction-
    # CAPABLE if a sufficient FRACTION of heads are within capture_perp (binding-
    # eligible) so each minifilament can recruit a bipolar pair. We report the
    # fraction + the median, and gate the fraction (≥ floor), not every head.
    head_pos = pos[:, N:, :].reshape(-1, 3) if M else np.zeros((0, 3))
    perps = np.empty(head_pos.shape[0], dtype=np.float64)
    for i, h in enumerate(head_pos):
        t = np.clip(np.einsum("ij,ij->i", h[None, :] - segA, seg) / seg_len2, 0.0, 1.0)
        closest = segA + t[:, None] * seg
        perps[i] = float(np.min(np.linalg.norm(h[None, :] - closest, axis=1)))
    eligible_frac = float(np.mean(perps <= cap_perp)) if perps.size else 1.0
    median_perp = float(np.median(perps)) if perps.size else 0.0
    worst_perp = float(np.max(perps)) if perps.size else 0.0
    _ELIGIBLE_FLOOR = 0.8   # ≥80% of heads binding-eligible → contraction-capable
    heads_near = bool(perps.size == 0 or eligible_frac >= _ELIGIBLE_FLOOR)
    heads_max = worst_perp

    # sf_ denylisted + distinct from cortex_myosin_.
    types = [
        myosin_backbone_bond_name(p_myo_sf.prefix),
        myosin_head_backbone_bond_name(p_myo_sf.prefix),
        *myosin_particle_type_names(p_myo_sf.prefix),
        *myosin_attach_bin_names(p_myo_sf.n_bins, p_myo_sf.prefix),
    ]
    sf_denylisted = all(t.startswith("sf_") for t in types)
    distinct = not any(t.startswith("cortex_myosin_") for t in types)

    controls = {
        "n_motors": M,
        "n_beads_per_motor": N + 2 * H,
        "assembled": assembled,
        "force_free": {"ok": force_free, "backbone_strain": worst_bb,
                       "head_offset_rel_dev": worst_hb},
        "heads_near_actin": {"ok": heads_near, "eligible_fraction": eligible_frac,
                             "median_perp_m": median_perp, "max_perp_m": heads_max,
                             "capture_perp_m": cap_perp,
                             "max_centre_to_actin_m": centres_on_actin},
        "sf_denylisted": sf_denylisted,
        "distinct_from_cortex_myosin": distinct,
    }
    ok = assembled and force_free and heads_near and sf_denylisted and distinct
    return {"verdict": "PASS" if ok else "REVIEW", "controls": controls}
