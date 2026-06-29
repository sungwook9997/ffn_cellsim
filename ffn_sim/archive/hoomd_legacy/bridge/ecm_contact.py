"""Cell↔ECM CONTACT manifold — the geometric contact relation (GEOMETRY ONLY).

The (b) spatial-substrate layer beneath the FA-traction field
(``bridge/manifold_traction.py``). Where ``manifold_traction`` bins the EXPLICIT
integrin↔ligand bond *forces* onto cell-surface patches, this module establishes
the underlying *geometric contact relation* between the cell's basal surface
(the :class:`~ffn_sim.common.surface_manifold.SurfaceManifold` patches) and the
ECM (its substrate ligand sites): **which patches are in adhesive contact, the
cell-surface↔ECM correspondence on the shared manifold frame, and the contact
gap / area** — independent of whether any bond is currently engaged.

This is the shared-coordinate infrastructure the H.7 surface-manifold design
calls out (``H7_SURFACE_MANIFOLD_EXPLICIT_CORTEX_2026-06-07.md`` §"value is the
SPATIAL/contact infrastructure … FA/ECM contact patches, shared frames"): the
cell side (manifold patches + local frames) and the ECM side (ligand sites) are
co-registered in one coordinate system, so the traction field, the lamellipodium
region masks (``cortex/manifold_regions.py``), and any later ECM-deformation
readout all live on the same patches.

Division of labour (HARD — same contract as ``manifold_traction.py``)
--------------------------------------------------------------------
* The :class:`SurfaceManifold` supplies **GEOMETRY ONLY**: patches (triangles),
  centroids, outward normals ``n̂`` and tangent frames ``(e1, e2)``, areas, and
  the bead→patch map (:meth:`SurfaceManifold.nearest_patch`). No mechanics.
* The ECM side supplies only **POSITIONS** — the substrate ligand particle
  positions (type ``"ligand"`` / :data:`TYPE_LIGAND`), read from the runtime
  snapshot exactly as ``manifold_traction`` reads them. This module creates **no
  force, no bond, no spring, no γ** — it computes a geometric *relation*, not a
  dynamics. ``numpy`` + the manifold's geometry only (no HOOMD force/integrator
  construction).

What "contact" means here
-------------------------
A cell-surface patch is in **adhesive contact** with the ECM iff at least one
ligand sits within a physical contact reach ``contact_gap`` of that patch
(measured to the patch centroid — the broad-phase home-patch distance, the same
nearest-centroid metric ``nearest_patch`` uses). The contact footprint is then
the union of those patches; it is the geometric basal footprint the cell adheres
over, and it must CONTAIN the FA-engaged patches (a bond can only form where a
ligand is reachable). That containment is the cross-consistency gate
(:func:`assert_traction_within_contact`).

Each ligand is mapped to exactly one home patch (a partition of the ligand set
over patches), so the manifold is the shared index between the ECM ligand sites
and the cell-surface patches.

Sanity Gate (per CLAUDE.md Sanity Gate Protocol; checked in
``tests/test_ecm_contact.py``)
-----------------------------------------------------------------------------
* **Dimensional.** Positions, centroids, ``contact_gap`` are lengths [m]; the
  per-ligand contact distance [m]; patch areas [m²]; contact area [m²];
  ``n̂``/frames dimensionless units. No quantity enters a force budget.
* **Boundary.** Zero ligands → empty contact mask, ``contact_area = 0`` (no
  raise). A ligand exactly ``contact_gap`` from its home patch is in contact
  (``≤``). ``contact_gap`` must be finite > 0.
* **Conservation / partition.** Every ligand maps to exactly one home patch
  (Σ per-patch ligand counts == n_ligand); no ligand is created or dropped.
* **Sign / sense.** The per-patch contact normal is the manifold's OUTWARD ``n̂``;
  the signed normal offset ``(ligand − centroid)·n̂`` is reported (positive =
  ligand outward of the patch, toward the substrate, for a basal patch).
* **Numerical.** float64 throughout; positions read tag-ordered (a HOOMD
  snapshot is tag-ordered) exactly as ``manifold_traction``.
* **Measurement-protocol consistency.** Contact is measured with the SAME
  nearest-centroid metric the FA-traction binning uses, so the contact patch set
  is a superset of the FA basal patch set on the same manifold — the
  cross-consistency gate.
"""

from __future__ import annotations

from typing import Any

import numpy as np

from ffn_sim.common.surface_manifold import SurfaceManifold

__all__ = [
    "TYPE_LIGAND",
    "bin_ecm_contact",
    "measure_ecm_contact_manifold",
    "assert_traction_within_contact",
]

# Substrate ligand particle type — the SINGLE source of truth is
# ``ffn_sim.archive.hoomd_legacy.bridge.fa.TYPE_LIGAND`` (== ``cell.cell.FA_TYPE_SUBSTRATE_LIGAND``).
# Re-declared here (rather than importing the heavy fa/cell modules) so this pure
# geometry module stays import-light; the driver asserts it against fa.py.
TYPE_LIGAND: str = "ligand"


def bin_ecm_contact(
    manifold: SurfaceManifold,
    ligand_positions: np.ndarray,
    *,
    contact_gap: float,
) -> dict[str, Any]:
    """Bin ECM ligand sites onto cell-surface patches → the contact relation.

    Pure-numpy core (no HOOMD): given the manifold geometry and the ECM ligand
    positions, compute which patches are in adhesive contact, the ligand→patch
    partition, the per-patch contact gap, and the contact area.

    Args:
        manifold: A :class:`SurfaceManifold` fitted to the cell surface (GEOMETRY
            ONLY).
        ligand_positions: ``(n_lig, 3)`` substrate ligand positions [m].
        contact_gap: Physical contact reach [m]; a ligand within this distance of
            its home patch centroid is in contact with that patch. Caller-supplied
            (the physics reach — e.g. FA ``capture_radius`` + ``h_integrin`` plus
            the patch geometric slack); not a tuned constant.

    Returns:
        A dict with the per-patch + ECM-side contact relation::

            {
              # --- per-patch arrays (length n_tri) ---
              "patch_centroids": (n_tri, 3) [m],
              "patch_normals":   (n_tri, 3),       # outward n̂
              "patch_area_m2":   (n_tri,) [m²],
              "patch_n_ligands": (n_tri,) int,     # ligands homed on this patch
              "patch_min_gap_m": (n_tri,) [m],     # nearest ligand distance (inf if none)
              "patch_normal_offset_m": (n_tri,) [m],  # (nearest lig − centroid)·n̂ (signed)
              "contact_mask":    (n_tri,) bool,    # ≥1 ligand within contact_gap
              "contact_patch_ids": (n_contact,) int,
              # --- ECM side (ligand → patch partition) ---
              "ligand_home_patch": (n_lig,) int,   # home patch of each ligand
              "ligand_gap_m":      (n_lig,) [m],   # distance to home centroid
              "ligand_in_contact": (n_lig,) bool,  # gap ≤ contact_gap
              # --- totals ---
              "n_ligands": int,
              "n_ligands_in_contact": int,
              "n_contact_patches": int,
              "contact_area_m2": float,            # Σ area of contact patches
              "contact_area_fraction": float,      # / manifold total area
              "mean_contact_gap_m": float,         # mean over in-contact ligands
              "contact_gap_m": float,              # the reach used
              # --- partition gate ---
              "partition": {"passed": bool, "sum_patch_counts": int, "n_lig": int},
            }
    """
    if not (np.isfinite(contact_gap) and contact_gap > 0.0):
        raise ValueError(f"contact_gap must be finite > 0; got {contact_gap!r}")

    n_tri = manifold.n_tri
    centroids = np.asarray(manifold.tri_centroids, dtype=np.float64)
    normals = np.asarray(manifold.tri_normals, dtype=np.float64)
    areas = np.asarray(manifold.tri_areas, dtype=np.float64)
    total_area = float(np.sum(areas))

    lig = np.asarray(ligand_positions, dtype=np.float64)
    if lig.size == 0:
        lig = np.empty((0, 3), dtype=np.float64)
    if lig.ndim != 2 or lig.shape[1] != 3:
        raise ValueError(f"ligand_positions must be (n, 3); got {lig.shape}")
    n_lig = int(lig.shape[0])

    patch_n_ligands = np.zeros(n_tri, dtype=np.int64)
    patch_min_gap = np.full(n_tri, np.inf, dtype=np.float64)
    patch_normal_offset = np.zeros(n_tri, dtype=np.float64)
    contact_mask = np.zeros(n_tri, dtype=bool)

    if n_lig > 0:
        home = manifold.nearest_patch(lig)                       # (n_lig,)
        delta = lig - centroids[home]                            # (n_lig, 3)
        gap = np.linalg.norm(delta, axis=1)                      # (n_lig,)
        normal_off = np.einsum("ij,ij->i", delta, normals[home])  # signed
        in_contact = gap <= contact_gap

        np.add.at(patch_n_ligands, home, 1)
        # Per-patch nearest ligand gap + its signed normal offset.
        order = np.argsort(-gap)  # process far→near so the nearest overwrites
        patch_min_gap[home[order]] = gap[order]
        patch_normal_offset[home[order]] = normal_off[order]
        # Contact: a patch is in contact iff it has ≥1 IN-CONTACT ligand.
        contact_home = home[in_contact]
        contact_mask[contact_home] = True
    else:
        home = np.empty((0,), dtype=np.int64)
        gap = np.empty((0,), dtype=np.float64)
        in_contact = np.empty((0,), dtype=bool)

    contact_ids = np.flatnonzero(contact_mask)
    contact_area = float(np.sum(areas[contact_mask]))
    n_in_contact = int(np.count_nonzero(in_contact))
    mean_gap = float(np.mean(gap[in_contact])) if n_in_contact else 0.0

    sum_counts = int(np.sum(patch_n_ligands))
    partition_ok = bool(sum_counts == n_lig)

    return {
        "patch_centroids": centroids,
        "patch_normals": normals,
        "patch_area_m2": areas,
        "patch_n_ligands": patch_n_ligands,
        "patch_min_gap_m": patch_min_gap,
        "patch_normal_offset_m": patch_normal_offset,
        "contact_mask": contact_mask,
        "contact_patch_ids": contact_ids,
        "ligand_home_patch": home,
        "ligand_gap_m": gap,
        "ligand_in_contact": in_contact,
        "n_ligands": n_lig,
        "n_ligands_in_contact": n_in_contact,
        "n_contact_patches": int(contact_ids.size),
        "contact_area_m2": contact_area,
        "contact_area_fraction": (
            contact_area / total_area if total_area > 0.0 else 0.0
        ),
        "mean_contact_gap_m": mean_gap,
        "contact_gap_m": float(contact_gap),
        "partition": {
            "passed": partition_ok,
            "sum_patch_counts": sum_counts,
            "n_lig": n_lig,
        },
    }


def _read_ligand_positions(sim: Any, ligand_type: str) -> np.ndarray:
    """Read tag-ordered ligand positions from a HOOMD simulation.

    Mirrors ``manifold_traction._read_tag_ordered_positions_and_bonds`` (a HOOMD
    snapshot is tag-ordered; we re-index into tag order, then select the rows
    whose particle type is ``ligand_type``).

    Args:
        sim: A built HOOMD ``Simulation``.
        ligand_type: The substrate ligand particle-type name (``"ligand"``).

    Returns:
        ``(n_lig, 3)`` ligand positions [m] (empty if the type is absent or has
        no particles).
    """
    particle_types = list(sim.state.particle_types)
    try:
        lig_tid = particle_types.index(ligand_type)
    except ValueError:
        return np.empty((0, 3), dtype=np.float64)
    with sim.state.cpu_local_snapshot as s:
        pos = np.asarray(s.particles.position)
        tg = np.asarray(s.particles.tag)
        typeid = np.asarray(s.particles.typeid)
        inv = np.empty_like(tg)
        inv[tg] = np.arange(tg.size)
        pos_by_tag = pos[inv].copy()
        typeid_by_tag = typeid[inv].copy()
    return pos_by_tag[typeid_by_tag == lig_tid].copy()


def measure_ecm_contact_manifold(
    sim: Any,
    manifold: SurfaceManifold,
    *,
    contact_gap: float,
    ligand_type: str = TYPE_LIGAND,
) -> dict[str, Any]:
    """Measure the cell↔ECM contact relation from a built FA-adhered cell.

    Reads the substrate ligand positions from ``sim`` and bins them onto the
    manifold patches (:func:`bin_ecm_contact`). The manifold supplies geometry;
    the ECM supplies positions; no force is created.

    Args:
        sim: A built, FA-adhered HOOMD ``Simulation`` (e.g. ``cell.simulation``).
        manifold: A :class:`SurfaceManifold` fitted to the cell surface.
        contact_gap: Physical contact reach [m] (see :func:`bin_ecm_contact`).
        ligand_type: Substrate ligand particle-type name (``"ligand"``).

    Returns:
        The dict from :func:`bin_ecm_contact`.
    """
    lig = _read_ligand_positions(sim, ligand_type)
    return bin_ecm_contact(manifold, lig, contact_gap=contact_gap)


def assert_traction_within_contact(
    contact: dict[str, Any], traction: dict[str, Any]
) -> dict[str, Any]:
    """Cross-consistency gate: FA-engaged patches ⊆ geometric contact patches.

    An integrin↔ligand bond can only form where a ligand is reachable, so every
    patch carrying an engaged FA bond (``manifold_traction``'s
    ``basal_patch_mask``) MUST lie inside the geometric contact footprint
    (:func:`bin_ecm_contact`'s ``contact_mask``) measured on the SAME manifold.
    A violation means the two layers disagree (a bug or mismatched ``contact_gap``
    / manifold) — the caller should treat it as fatal.

    Args:
        contact: The dict from :func:`bin_ecm_contact` /
            :func:`measure_ecm_contact_manifold`.
        traction: The dict from
            :func:`ffn_sim.archive.hoomd_legacy.bridge.manifold_traction.measure_fa_traction_field`
            (must share the same manifold → same ``n_tri``).

    Returns:
        A gate dict ``{"passed": bool, "n_fa_patches": int,
        "n_fa_outside_contact": int, "fa_outside_ids": (k,) int}``.

    Raises:
        ValueError: if the two masks have different lengths (different manifolds).
    """
    contact_mask = np.asarray(contact["contact_mask"], dtype=bool)
    fa_mask = np.asarray(traction["basal_patch_mask"], dtype=bool)
    if contact_mask.shape != fa_mask.shape:
        raise ValueError(
            f"contact mask {contact_mask.shape} and FA mask {fa_mask.shape} have "
            "different lengths — they must be measured on the same manifold."
        )
    outside = fa_mask & ~contact_mask
    outside_ids = np.flatnonzero(outside)
    return {
        "passed": bool(outside_ids.size == 0),
        "n_fa_patches": int(np.count_nonzero(fa_mask)),
        "n_fa_outside_contact": int(outside_ids.size),
        "fa_outside_ids": outside_ids,
    }
