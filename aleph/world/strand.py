"""Building a strand into the arena: the first composition, and the first place "no constants" bites.

A STRAND is the first object in this engine — a node alone is inert, a segment is a constraint, an
angle is a force, and only their composition is a filament that can be named, counted and rendered.
This module builds one into a :class:`~aleph.world.arena.WorldArena` by claiming three contiguous
ranges and writing the topology between them::

    N nodes  +  (N-1) segments  +  (N-2) angle triples  +  polarity  +  a tip

WHAT IT DELIBERATELY DOES NOT DO.  It computes no force and assembles no operator.  A segment here
carries a rest length and a MATERIAL ARC COORDINATE and nothing else; whether the axial law is a spring
or a projected constraint is the solver's business, and this engine has already settled that one —
``cortex_state.CORTEX_CHANNELS_NOT_BOUND`` records axial inextensibility as "an NF2007 constraint
carried by the solver, not an accumulated force", with the note that binding a spring there would
double-count the backbone.  So the topology is built and the law is left where it lives.

THE ARC COORDINATE IS THE POINT.  Each segment records the rest-configuration arc interval it spans.
That is what a motor binds to — an address ``(segment, u)`` on a continuum — and it is why a node per
G-actin monomer is unnecessary: a binding site is a position along the material, not an allocation.
Discretising the filament at 30 nm and binding at 2.7 nm resolution are compatible statements, and the
arc coordinate is what makes them so.

WHY THERE ARE NO DEFAULTS HERE.  Under the PI's ruling of 2026-08-15, every value labelled
physiological is physiological FOR SOME PARTICULAR CELL and is not certain, so it is a declared axis
rather than a literal.  ``contour_um`` and ``seg_um`` therefore have no defaults and the build refuses
without them — the same discipline as ``MembraneAreaCard.capacity``, which has no default because
"supplying one would recreate ``RESERVOIR_STRAIN = 0.60`` … under a new name".  A default segment
length would recreate the ``seg_um = 0.5`` the provenance audit ranks #2 load-bearing with a source of
"none (0 KB rows)".

AND THE REALISED LENGTH IS REPORTED, NOT THE REQUESTED ONE.  A contour of 3.0 µm at a requested 0.05 µm
gives 60 segments exactly, but 3.0 at 0.055 gives 55 segments of 0.0545 µm.  The requested figure is
then not the one the physics sees.  Both are recorded and :attr:`Strand.seg_um_realised` is the one a
downstream law must read, because a rounding that changes a length by 1% is the kind of thing that is
invisible until it is quoted.

Sanity Gate (recorded per CLAUDE.md, before first execution):
  * dimensional — positions, contour and segment lengths [µm]; the arc coordinate is a length [µm]
    measured along the rest configuration; polarity is a sign and carries no unit.
  * boundary — a strand of fewer than 3 nodes has no angle triple and therefore no bending at all; it
    is rejected rather than built silently flexible, since a filament with no bending stiffness is a
    different physical object.  A non-finite or non-positive length is rejected.
  * conservation/invariant — the node, segment and angle ranges are claimed from one arena and the
    topology indexes only inside them, so a strand cannot address another population's nodes.  Asserted.
  * CFL/precision — no integration and no force; float64 throughout.  The realised segment length feeds
    the bending coefficient ``alpha = kappa / seg^3``, where a 1% length error is a 3% coefficient
    error, which is why the realised value is the reported one.
  * sign sense — ``polarity`` is +1 when the barbed end is the LAST node and -1 when it is the first;
    the tip node is derived from it rather than passed separately, so the two cannot disagree.
  * measurement protocol — host-side construction only; nothing is uploaded and no device is touched.

engine units: length µm.  Runtime: this module is pure host geometry and is CPU-importable.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import numpy.typing as npt

from aleph.world.arena import Claim, Kind, WorldArena

__all__ = ["Strand", "build_strand"]


@dataclass(frozen=True, slots=True)
class Strand:
    """One built filament: its claims, its topology, and the length the physics will actually see.

    Attributes:
        population: the population the ranges belong to.
        nodes / segments / angles: the three claims, contiguous within their kinds.
        position: ``(N, 3)`` rest positions [µm], in node order along the strand.
        seg_node: ``(N-1, 2)`` segment endpoints, as GLOBAL node indices into the arena.
        seg_rest_um: ``(N-1,)`` rest length of each segment [µm].
        seg_arc_um: ``(N-1, 2)`` the rest-configuration arc interval ``[s0, s1)`` each segment spans
            [µm] — the coordinate a motor binds at, and the reason monomer nodes are unnecessary.
        angle_idx: ``(N-2, 3)`` bending triples, as GLOBAL node indices.
        polarity: ``+1`` if the barbed end is the last node, ``-1`` if it is the first.
        tip_node: GLOBAL index of the barbed end, derived from ``polarity``.
        contour_um_requested / contour_um_realised: the asked-for and built contour lengths [µm].
        seg_um_requested / seg_um_realised: likewise.  **Downstream laws must read the realised one.**
    """

    population: str
    nodes: Claim
    segments: Claim
    angles: Claim
    position: npt.NDArray[np.float64]
    seg_node: npt.NDArray[np.int64]
    seg_rest_um: npt.NDArray[np.float64]
    seg_arc_um: npt.NDArray[np.float64]
    angle_idx: npt.NDArray[np.int64]
    polarity: int
    tip_node: int
    contour_um_requested: float
    contour_um_realised: float
    seg_um_requested: float
    seg_um_realised: float

    @property
    def n_nodes(self) -> int:
        """Number of nodes in this strand."""
        return self.nodes.count

    @property
    def seg_um_error(self) -> float:
        """Relative difference between the requested and realised segment length, dimensionless.

        Reported rather than checked: rounding a contour into a whole number of segments HAS to move
        one of the two lengths, and which one moved is a fact the artifact should carry rather than a
        failure to prevent.
        """
        return abs(self.seg_um_realised - self.seg_um_requested) / self.seg_um_requested


def _require_length(name: str, value: float | None) -> float:
    """Return a declared positive length, or refuse.

    ``None`` is refused rather than defaulted. A default here would be an unsourced constant governing
    the physics under a new name, which is exactly what the ratified ``capacity`` pattern exists to
    prevent.
    """
    if value is None:
        raise ValueError(
            f"{name} has no default and must be declared. Every value labelled physiological is "
            "physiological for some particular cell and is not certain, so it arrives as a swept axis "
            "with a scope — never as a literal chosen here."
        )
    v = float(value)
    if not np.isfinite(v) or v <= 0.0:
        raise ValueError(f"{name} must be finite and positive; got {value!r}")
    return v


def build_strand(
    arena: WorldArena,
    population: str,
    *,
    start: npt.ArrayLike,
    direction: npt.ArrayLike,
    contour_um: float | None = None,
    seg_um: float | None = None,
    polarity: int = +1,
) -> Strand:
    """Claim ranges for one strand and build its topology into ``arena``.

    Args:
        arena: the world to claim from.  Node, segment and angle ranges are taken from its tail.
        population: the name the ranges are attributed to, e.g. ``"cortex"``.
        start: ``(3,)`` position of the first node [µm].
        direction: ``(3,)`` direction the strand runs in; normalised here, so its magnitude is ignored.
        contour_um: rest contour length [µm].  **REQUIRED** — no default.
        seg_um: requested discretisation length [µm].  **REQUIRED** — no default.  The realised value
            is reported on the returned :class:`Strand` and is the one a downstream law must read.
        polarity: ``+1`` if the barbed end is the last node, ``-1`` if it is the first.

    Returns:
        The built :class:`Strand`.

    Raises:
        ValueError: on a missing or non-positive length, a degenerate direction, a polarity that is
            not ±1, or a discretisation giving fewer than 3 nodes (which would leave the filament with
            no bending triple at all, i.e. a different physical object).
    """
    contour = _require_length("contour_um", contour_um)
    seg = _require_length("seg_um", seg_um)
    if polarity not in (+1, -1):
        raise ValueError(f"polarity must be +1 or -1; got {polarity!r}")

    d = np.asarray(direction, np.float64).reshape(3)
    norm = float(np.linalg.norm(d))
    if not np.isfinite(norm) or norm <= 0.0:
        raise ValueError(f"direction is degenerate; got {direction!r}")
    d = d / norm
    origin = np.asarray(start, np.float64).reshape(3)

    n_seg = int(round(contour / seg))
    n_nodes = n_seg + 1
    if n_nodes < 3:
        raise ValueError(
            f"contour_um={contour} at seg_um={seg} gives {n_nodes} nodes; a strand needs at least 3 so "
            "it carries one bending triple. A filament with no bending stiffness is a different "
            "physical object, not a coarser one."
        )
    seg_realised = contour / n_seg

    nodes = arena.claim(population, Kind.NODE, n_nodes)
    segments = arena.claim(population, Kind.SEGMENT, n_seg)
    angles = arena.claim(population, Kind.ANGLE3, n_nodes - 2)

    arc = np.arange(n_nodes, dtype=np.float64) * seg_realised
    position = origin[None, :] + arc[:, None] * d[None, :]

    gid = nodes.lo + np.arange(n_nodes, dtype=np.int64)
    seg_node = np.stack([gid[:-1], gid[1:]], axis=1)
    seg_rest = np.full(n_seg, seg_realised, np.float64)
    seg_arc = np.stack([arc[:-1], arc[1:]], axis=1)
    angle_idx = np.stack([gid[:-2], gid[1:-1], gid[2:]], axis=1)

    tip_node = int(gid[-1]) if polarity > 0 else int(gid[0])

    strand = Strand(
        population=population, nodes=nodes, segments=segments, angles=angles,
        position=position, seg_node=seg_node, seg_rest_um=seg_rest, seg_arc_um=seg_arc,
        angle_idx=angle_idx, polarity=polarity, tip_node=tip_node,
        contour_um_requested=contour, contour_um_realised=float(arc[-1]),
        seg_um_requested=seg, seg_um_realised=seg_realised,
    )
    assert_inside_arena(strand)
    return strand


def assert_inside_arena(strand: Strand) -> None:
    """Assert every index the strand writes falls inside its own claimed node range.

    A strand that reaches outside its range would silently couple two populations through shared nodes
    — which is a weld, not a connection, and the one arrangement this engine's architecture forbids
    outright. Cheap, and it fails at build rather than at the first accepted step.

    Raises:
        AssertionError: if any segment or angle index leaves the node claim.
    """
    lo, hi = strand.nodes.lo, strand.nodes.hi
    for name, idx in (("segment", strand.seg_node), ("angle", strand.angle_idx)):
        if idx.size and (int(idx.min()) < lo or int(idx.max()) >= hi):
            raise AssertionError(
                f"{strand.population} {name} topology addresses [{int(idx.min())}, {int(idx.max())}] "
                f"outside its node claim [{lo}, {hi}). Sharing a node across populations is a weld, "
                "not a connection."
            )


def _demo() -> None:
    """Self-check: the composition rule, the refusals, and the realised-length report."""
    arena = WorldArena(capacity={Kind.NODE: 1_000, Kind.SEGMENT: 1_000, Kind.ANGLE3: 1_000})

    # A 3.0 µm filament at 0.05 µm: 60 segments exactly, so realised == requested.
    s = build_strand(arena, "cortex", start=(0.0, 0.0, 0.0), direction=(1.0, 0.0, 0.0),
                     contour_um=3.0, seg_um=0.05)
    assert s.n_nodes == 61 and s.segments.count == 60 and s.angles.count == 59, "N / N-1 / N-2"
    assert s.seg_um_error == 0.0
    assert np.isclose(s.contour_um_realised, 3.0)
    assert np.isclose(np.linalg.norm(s.position[-1] - s.position[0]), 3.0)
    assert s.tip_node == s.nodes.hi - 1, "polarity +1 puts the barbed end last"

    # 3.0 at 0.055 does NOT divide: 55 segments of 0.0545 µm. The realised value is what physics sees.
    r = build_strand(arena, "cortex", start=(0.0, 1.0, 0.0), direction=(1.0, 0.0, 0.0),
                     contour_um=3.0, seg_um=0.055)
    assert r.segments.count == 55
    assert r.seg_um_realised != r.seg_um_requested
    assert 0.0 < r.seg_um_error < 0.02, f"rounding moved the length by {r.seg_um_error:.4f}"
    assert np.allclose(r.seg_rest_um, r.seg_um_realised)

    # The arc coordinate tiles the contour with no gap — it is what a motor binds at.
    assert np.allclose(r.seg_arc_um[:-1, 1], r.seg_arc_um[1:, 0])
    assert np.isclose(r.seg_arc_um[0, 0], 0.0)
    assert np.isclose(r.seg_arc_um[-1, 1], r.contour_um_realised)

    # Claims are contiguous and attributed; the second strand starts where the first ended.
    assert r.nodes.lo == s.nodes.hi
    arena.assert_partitioned()
    assert arena.population_of(Kind.NODE, s.nodes.lo) == "cortex"

    # Both lengths REFUSE to default — the ratified pattern, applied at the first opportunity.
    for kwargs in ({"contour_um": 3.0}, {"seg_um": 0.05}, {}):
        try:
            build_strand(arena, "cortex", start=(0, 0, 0), direction=(1, 0, 0), **kwargs)
        except ValueError as exc:
            assert "no default" in str(exc)
        else:  # pragma: no cover
            raise AssertionError(f"build_strand{kwargs} must refuse")

    # Too coarse to carry a bending triple is a different object, not a coarser one.
    try:
        build_strand(arena, "cortex", start=(0, 0, 0), direction=(1, 0, 0),
                     contour_um=1.0, seg_um=1.0)
    except ValueError as exc:
        assert "at least 3" in str(exc)
    else:  # pragma: no cover
        raise AssertionError("a 2-node strand must refuse")

    print("strand self-check OK —",
          f"{s.n_nodes} nodes / {s.segments.count} segments / {s.angles.count} angles; "
          f"realised seg {r.seg_um_realised:.6f} µm vs requested {r.seg_um_requested}")


if __name__ == "__main__":
    _demo()
