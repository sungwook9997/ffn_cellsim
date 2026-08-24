"""The nuclear envelope: the same closed-surface primitive as the plasma membrane, at its own radius.

WHY THIS IS A SEPARATE FILE AND NOT A SECOND ARGUMENT.  Geometrically a nuclear envelope is exactly
what :mod:`aleph.world.build.membrane` builds — a closed icosphere with faces, hinges and a discrete
enclosed volume — and this module adds no geometry at all.  What it adds is the two things that are
NOT shared: where the radius comes from, and the fact that the mesh spacing has no source.  Both are
provenance, both are population-specific, and folding them into a ``population="nucleus"`` argument
would put a sourced radius and an unsourced one behind the same parameter name.

THE RADIUS IS DERIVED FROM A MEASURED RATIO, NOT TYPED.  ``aleph.laws.cell_geometry`` carries the MCF7
nucleus:cell RADIUS ratio as a measured distribution — **0.68 ± 0.08**, from imaging flow cytometry at
n = 2164 and confirmed independently by photoacoustic sizing at n = 37 — so the nuclear radius is
``nc_ratio x R_cell`` and moves with the cell radius rather than beside it.  At the volume-anchored
MCF7 central radius that is 0.68 x 7.5 = **5.1 µm**.

⚠ **That conflicts with the incumbent, and the incumbent says so itself.**
``components/incumbent/compartments.py`` carries ``r_eq_um = 5.0`` labelled *"provisional geometry
(N:C vol ~0.30); PI GAP"*.  Two values for one quantity.  This module takes the ``laws/`` one because
it is sourced and because ``components/`` is frozen PORT SOURCE whose own comment declines to defend
its number — but the conflict is reported to the PI rather than closed here.  The 2% difference is
not the point; a sourced value and a self-declared gap sitting in one engine is.

⚠ **PI-GAP — THE ENVELOPE'S MESH SPACING HAS NO SOURCE IN THIS REPOSITORY.**  The membrane's
subdivision is derived from the 50-100 nm cortical / membrane-skeleton mesh (Morone 2006, Bovellan
2014, Chugh & Paluch 2018).  That band is a property of the actin-spectrin cortex, and the nuclear
envelope's structural mesh is the LAMIN meshwork, which is a different network with a different face
size.  No lamin meshwork spacing is registered anywhere in this repository.  Using the cortical band
for the envelope is therefore an **analogy, not a derivation**, and :data:`ENVELOPE_MESH_PI_GAP` is
the label that says so.  :func:`build_envelope` refuses to run without a ``mesh_provenance`` string, so
the label cannot be lost between the decision and the artifact — it is written into the run record by
:meth:`ClosedSurface.record`.  When a sourced lamin spacing arrives, the subdivision changes by itself.

WHAT IS DELIBERATELY ABSENT.  No lamina areal tension, no Helfrich bending, no nucleoplasm volume law,
no LINC tether, no chromatin.  PHASE 1 is geometry.  Note in particular that the envelope built here
is a SPHERE, not the oblate spheroid the incumbent's ``build_oblate_mesh`` produces: an aspect ratio is
a shape parameter with its own provenance, it is 1.0 in the incumbent's own default, and inventing one
here would be a physiological decision disguised as a builder default.

Sanity Gate (recorded per CLAUDE.md, before first execution):
  * dimensional — radii and mesh spacing [µm]; the N:C ratio is dimensionless.
  * boundary — the ratio must lie in ``(0, 1)``: a nucleus at or beyond the cell radius is not a
    nucleus, and the envelope must fit inside the membrane it will later couple to; a nuclear radius
    that does not is refused here rather than discovered as an inside-out cell.  Every refusal of
    :func:`~aleph.world.build.membrane.build_closed_surface` applies unchanged.
  * conservation/invariant — inherited from the shared builder: Euler at every level, every edge shared
    by exactly two faces, claims tiling the arena's live prefix.
  * CFL/precision — no integration; float64.
  * sign sense — inherited: outward winding, checked through a positive discrete enclosed volume.
  * measurement protocol — inherited; host readbacks are build-time only.

engine units: length µm, area µm², volume µm³.  Runtime: NVIDIA Warp on CUDA.  CPU-importable; the
builder requires a CUDA arena and raises otherwise.
"""

from __future__ import annotations

import math

from aleph.laws.cell_geometry import CellGeometryDist, resolve_cell_geometry
from aleph.world.arena import WorldArena
from aleph.world.build.membrane import ClosedSurface, build_closed_surface

__all__ = ["ENVELOPE_MESH_PI_GAP", "build_envelope", "nuclear_radius_um"]

#: The label a caller must carry when it reuses the cortical 50-100 nm mesh band for the envelope.
#: It is a sentence rather than a flag because it ends up verbatim in the run record, where a reader
#: who was not in this session has to be able to tell a derivation from an analogy.
ENVELOPE_MESH_PI_GAP = (
    "PI-GAP: no lamin meshwork spacing is registered in this repository. The value used is the "
    "cortical/membrane-skeleton mesh band (50-100 nm; Morone 2006, Bovellan 2014, Chugh & Paluch "
    "2018), which is a property of the actin-spectrin cortex and NOT of the nuclear lamina. This is "
    "an analogy, not a derivation, and the subdivision it produces inherits that status."
)


def nuclear_radius_um(geometry: CellGeometryDist | str = "mcf7") -> float:
    """Central nuclear radius [µm] = N:C radius ratio x cell radius, from sourced geometry.

    Reads ``aleph.laws.cell_geometry`` rather than carrying a number, so the nuclear radius moves with
    the cell radius: a swept ``R_cell`` cannot leave the nucleus behind at a value that was correct for
    a different cell.

    Args:
        geometry: a :class:`~aleph.laws.cell_geometry.CellGeometryDist` or a registry name/alias.

    Returns:
        The central nuclear radius [µm].

    Raises:
        ValueError: if the N:C ratio is not strictly inside ``(0, 1)`` — a nucleus at or past the cell
            radius is not a nucleus, and nothing downstream could place it inside the membrane.
    """
    dist = geometry if isinstance(geometry, CellGeometryDist) else resolve_cell_geometry(geometry)
    if not (math.isfinite(dist.nc_ratio) and 0.0 < dist.nc_ratio < 1.0):
        raise ValueError(
            f"{dist.name} carries an N:C radius ratio of {dist.nc_ratio}, which is not strictly inside "
            "(0, 1). A nucleus at or beyond the cell radius is not a nucleus; the envelope has to fit "
            "inside the membrane it will be coupled to."
        )
    return float(dist.R_nuc_um)


def build_envelope(
    arena: WorldArena,
    *,
    radius_um: float | None = None,
    mesh_um: float | None = None,
    mesh_provenance: str | None = None,
    centre_um: tuple[float, float, float] = (0.0, 0.0, 0.0),
    population: str = "nuclear_envelope",
) -> ClosedSurface:
    """The nuclear envelope, as a closed surface population.

    Args:
        arena: the world to claim from.
        radius_um: nuclear radius [µm].  **REQUIRED, no default.**  Get it from
            :func:`nuclear_radius_um` rather than typing one — that is what makes it derived.
        mesh_um: largest acceptable mean triangle edge [µm].  **REQUIRED, no default.**  The
            subdivision level follows from it.
        mesh_provenance: where ``mesh_um`` came from.  **REQUIRED, and unlike everywhere else in this
            package it may not be empty.**  The envelope's spacing is a live PI-GAP
            (:data:`ENVELOPE_MESH_PI_GAP`), and a gap that is not written into the artifact is a gap
            that gets quoted as a result later.  Pass :data:`ENVELOPE_MESH_PI_GAP` verbatim while the
            gap stands, and a citation once it does not.
        centre_um: nucleus centre [µm].  Concentric with the cell unless a caller says otherwise.
        population: the claim name.

    Returns:
        The built :class:`~aleph.world.build.membrane.ClosedSurface`.

    Raises:
        ValueError: on a missing radius, spacing or provenance, or any refusal from the shared builder.
    """
    if radius_um is None:
        raise ValueError(
            "radius_um has no default and must be declared. Derive it with nuclear_radius_um(), which "
            "reads the sourced N:C ratio, rather than typing a nuclear radius here."
        )
    if mesh_um is None:
        raise ValueError(
            "mesh_um has no default and must be declared. The subdivision level is DERIVED from it; a "
            "default spacing would be a typed subdivision level under a new name."
        )
    if not mesh_provenance or not mesh_provenance.strip():
        raise ValueError(
            "mesh_provenance must be a non-empty statement of where mesh_um came from. The nuclear "
            "envelope's mesh spacing is an OPEN PI-GAP: no lamin meshwork spacing is registered in "
            "this repository, so any value used is an analogy to the cortical mesh. Pass "
            "ENVELOPE_MESH_PI_GAP while that stands. A gap that is not written into the artifact is a "
            "gap that gets quoted as a result."
        )
    return build_closed_surface(
        arena, population, radius_um=radius_um, mesh_um=mesh_um,
        mesh_provenance=mesh_provenance, centre_um=centre_um,
    )


# ── self-check ──────────────────────────────────────────────────────────────────────────────────
def _demo() -> None:
    """Self-check: the derived radius, the level it implies, and the provenance that may not be empty."""
    from aleph.world.arena import Kind
    from aleph.world.build.membrane import icosphere_counts, mean_edge_um, subdivisions_for_mesh

    # The radius is the sourced ratio applied to the sourced cell radius, not a number typed here.
    dist = resolve_cell_geometry("mcf7")
    assert dist.R_cell_um == 7.5 and dist.nc_ratio == 0.68
    r_nuc = nuclear_radius_um()
    assert math.isclose(r_nuc, 0.68 * 7.5), "N:C x R_cell"
    assert math.isclose(r_nuc, 5.1)
    assert 0.0 < r_nuc < dist.R_cell_um, "the envelope has to fit inside the membrane"
    # It moves with the cell radius — the whole reason it is a ratio and not a length.
    assert nuclear_radius_um("ifc") > r_nuc, "a bigger cell carries a bigger nucleus"
    assert math.isclose(nuclear_radius_um("ifc"), 0.68 * 9.44)

    # The level follows from the spacing at THIS radius, and it is not the membrane's.
    level = subdivisions_for_mesh(r_nuc, 0.100)
    assert level == 6, "the nuclear radius reaches the band one level earlier than the cell radius"
    assert level != subdivisions_for_mesh(7.5, 0.100), "a derived level is radius-dependent"
    assert 0.050 <= mean_edge_um(r_nuc, level) <= 0.100, "inside the band it was derived from"
    assert icosphere_counts(level) == (40_962, 122_880, 81_920)

    # A ratio outside (0, 1) is not a nucleus.
    import dataclasses
    for bad in (0.0, 1.0, 1.2):
        try:
            nuclear_radius_um(dataclasses.replace(dist, nc_ratio=bad))
        except ValueError as exc:
            assert "not strictly inside" in str(exc)
        else:  # pragma: no cover
            raise AssertionError(f"an N:C ratio of {bad} must refuse")

    # Radius, spacing and provenance all refuse to default; provenance also refuses to be blank.
    arena = WorldArena(capacity={Kind.NODE: 10, Kind.FACE: 10, Kind.ANGLE4: 10})
    for kwargs, needle in (
        ({}, "radius_um has no default"),
        ({"radius_um": r_nuc}, "mesh_um has no default"),
        ({"radius_um": r_nuc, "mesh_um": 0.1}, "must be a non-empty statement"),
        ({"radius_um": r_nuc, "mesh_um": 0.1, "mesh_provenance": "   "}, "must be a non-empty statement"),
    ):
        try:
            build_envelope(arena, **kwargs)
        except ValueError as exc:
            assert needle in str(exc), f"{kwargs} raised the wrong refusal: {exc}"
        else:  # pragma: no cover
            raise AssertionError(f"build_envelope{kwargs} must refuse")

    # With a provenance it gets past the refusals and stops at the missing device, not before.
    try:
        build_envelope(arena, radius_um=r_nuc, mesh_um=0.1, mesh_provenance=ENVELOPE_MESH_PI_GAP)
    except RuntimeError as exc:
        assert "needs a CUDA arena" in str(exc)
    else:  # pragma: no cover
        raise AssertionError("a device=None arena must refuse to build")

    assert "PI-GAP" in ENVELOPE_MESH_PI_GAP and "lamin" in ENVELOPE_MESH_PI_GAP

    print(
        f"envelope self-check OK — N:C {dist.nc_ratio} x R_cell {dist.R_cell_um} um -> R_nuc "
        f"{r_nuc:.2f} um; mesh 0.100 um -> subdiv {level}: {icosphere_counts(level)[0]:,} verts / "
        f"{icosphere_counts(level)[2]:,} faces / {icosphere_counts(level)[1]:,} hinges "
        f"(mean edge {mean_edge_um(r_nuc, level) * 1e3:.1f} nm, PI-GAP)"
    )


if __name__ == "__main__":
    _demo()
