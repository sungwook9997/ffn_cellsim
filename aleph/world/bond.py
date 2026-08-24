"""Bond families, and the count that refuses to be a discretisation artefact.

A BOND is the only primitive that can break. Segments, angles and faces are structure — they define what
an object IS — while a bond is a stateful, identified pairing with kinetics, and it is the only kind of
connection this architecture declares. Field interaction (steric, drag) needs no declaration because it
depends on position alone; a specific persistent pairing does, because "same world" cannot express
"still attached to THAT one".

THE DEFECT THIS MODULE IS SHAPED AROUND.  In the incumbent, the membrane–cortex ERM linkage is declared
as ONE connector between two components. Physically it is a POPULATION of N tethers whose count is a
physiological areal density — and because the declaration hid the population, nobody had to answer how
many. The answer turned out to be: ``erm_density_per_um2`` defaults to ``None``, the build falls back to
one tether per membrane mesh vertex, and the count is therefore set by the DISCRETISATION. At subdiv 6
that is 58/µm²; the "physiological 235/µm²" written in a comment is the subdiv-7 icosphere vertex
density, 231.8/µm², i.e. a mesh number wearing a physiological label.

So :class:`BondCount` is REQUIRED to construct a :class:`BondFamily`, and it refuses without a scope, a
provenance and a measured support. A family cannot exist without answering "how many, against what, for
which cell, on whose authority". That is the whole design; everything else here is arithmetic.

THE COMPONENT PAIR IS DERIVED, NOT DECLARED.  A bond stores two GLOBAL node indices and nothing else.
Which populations it joins is read off the arena's ID ranges at query time. Declaring the pair on the
edge is what allowed one declaration to stand in for a population, and it is also what forced four
identical NMII motor edges and an INTERNAL/INTER split that is bookkeeping rather than physics.

WHAT IS DELIBERATELY ABSENT.  No kinetic law, no attach/detach, no free list, no snapshot twins. Nothing
detaches yet, so a free list would be machinery for a case that does not exist; no transaction runs yet,
so snapshot twins would be state nothing restores. They arrive with the first family that has kinetics,
which is also when the right shape for them becomes measurable rather than guessable.

Sanity Gate (recorded per CLAUDE.md, before first execution):
  * dimensional — ``rest_um`` [µm]; ``stiffness_pn_per_um`` [pN/µm]; an areal density is [1/µm²],
    volumetric [1/µm³], per-filament dimensionless. The support carries the matching unit and the
    resolver multiplies them, so a mismatched basis is a unit error and is rejected by name.
  * boundary — a resolved count of zero is returned as zero and NOT rejected: a density low enough to
    give no bonds over the support is a physical answer, unlike a count that was never asked for. A
    count exceeding the arena's remaining BOND capacity raises with both numbers.
  * conservation/invariant — every bond index must lie inside the arena's live NODE prefix, so a bond
    cannot address an unclaimed node; asserted at build.
  * CFL/precision — no integration; float64. Note the recorded stiffness enters the CFL bound through
    the row sum of its endpoints, which is why it is stored per bond rather than per family.
  * sign sense — ``rest_um`` is a length and must be non-negative; a bond joining a node to itself is a
    self-loop with no direction and is rejected.
  * measurement protocol — host-side construction; no device is touched and nothing is uploaded.

engine units: length µm, stiffness pN/µm.  Runtime: pure host bookkeeping, CPU-importable.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import StrEnum
from typing import Literal

import numpy as np
import numpy.typing as npt

from aleph.units.provenance import PROVENANCE
from aleph.world.arena import Claim, Kind, WorldArena

__all__ = ["SourceClass", "Basis", "BondCount", "BondFamily", "build_bond_family"]

Basis = Literal["areal", "volumetric", "per_filament", "explicit"]

#: Unit the support must carry for each basis, so a mismatch is caught by name rather than by magnitude.
_SUPPORT_UNIT: dict[str, str] = {
    "areal": "um^2", "volumetric": "um^3", "per_filament": "filaments", "explicit": "none",
}


class SourceClass(StrEnum):
    """Where a value came from, extending the engine's LIVE vocabulary rather than re-spelling it.

    The first four are imported from :data:`aleph.units.provenance.PROVENANCE` and checked against
    it at import, for the reason ``driver.OMITTABLE_CHANNELS`` gives for sharing one vocabulary: a second
    spelling lets a typo pass silently.

    ``UNRATIFIED_PROXY`` is the fifth and it is the one the PI's 2026-08-15 ruling makes necessary. A
    value can be perfectly SOURCED and still be wrong for the run: Π₀ = 40 Pa is a real measurement — of
    HeLa. The sourcing is fine; the SCOPE is what was never recorded. So a value whose scope does not
    match the run's target cell is a proxy, not a source, however good the paper is.
    """

    SOURCED = "SOURCED"
    DERIVED = "DERIVED"
    CONVENIENCE = "CONVENIENCE"
    PI_GAP = "PI_GAP"
    UNRATIFIED_PROXY = "UNRATIFIED_PROXY"


_missing = set(PROVENANCE) - {c.value for c in SourceClass}
if _missing:  # pragma: no cover - a drift guard, not a branch
    raise ImportError(
        f"SourceClass has drifted from forces_manifest.PROVENANCE; missing {sorted(_missing)}. "
        "One vocabulary on purpose — a second spelling lets a typo pass silently."
    )


@dataclass(frozen=True, slots=True)
class BondCount:
    """How many bonds a family has CAPACITY for, against what, for which cell, on whose authority.

    Every field is required. There is no default anywhere in this class, and that is the point: the
    ERM count became a mesh artefact because "how many" had a fallback.

    ⚠ **CAPACITY, NEVER OCCUPANCY. The distinction was forced by the first kinetic family
    (2026-08-21) and it is not a nicety.** For the NMII head->actin crossbridge the two questions
    separate cleanly:

    * **capacity** — one site per head, ``2H`` per minifilament. Structural, static, answerable, and
      this is what a ``BondCount`` states.
    * **occupancy** — how many are attached right now. ``params_i0b3.yaml`` already forbids imposing
      it: the engaged fraction is EMERGENT from ``k_on``/``k_off`` and must not be set from a duty
      ratio.

    **An imposed occupancy is a latch wearing a measurement's label** — the ERM defect in its kinetic
    form. So a family that has kinetics reports its occupancy; it never accepts one.

    Attributes:
        basis: what the density is per — ``areal`` / ``volumetric`` / ``per_filament`` / ``explicit``.
        value: the density in that basis, or the count itself when the basis is ``explicit``.
        scope: the cell line, state, assay and temperature the value is true FOR. Required, because a
            sourced value out of scope is a proxy.
        source_class: see :class:`SourceClass`.
        provenance: the citation or the derivation. Required and non-empty.
    """

    basis: Basis
    value: float
    scope: str
    source_class: SourceClass
    provenance: str

    def __post_init__(self) -> None:
        if self.basis not in _SUPPORT_UNIT:
            raise ValueError(f"basis must be one of {sorted(_SUPPORT_UNIT)}; got {self.basis!r}")
        if not np.isfinite(self.value) or self.value < 0.0:
            raise ValueError(f"count value must be finite and nonnegative; got {self.value!r}")
        if not str(self.scope).strip():
            raise ValueError(
                "scope is required: a value labelled physiological is physiological FOR SOME PARTICULAR "
                "cell, and without saying which, a SOURCED value cannot be distinguished from a proxy."
            )
        if not str(self.provenance).strip():
            raise ValueError(
                "provenance is required and must be non-empty. A count with no stated origin is how a "
                "discretisation artefact acquires a physiological label."
            )

    @property
    def support_unit(self) -> str:
        """The unit the support must be given in for this basis."""
        return _SUPPORT_UNIT[self.basis]

    def resolve(self, support: float) -> int:
        """Resolve to a whole bond count against a MEASURED ``support``.

        Args:
            support: the measured area [µm²], volume [µm³] or filament count the density is taken
                against — measured, never analytic. The ERM builder that got this right counts against
                ``triangle_surface_area`` of the actual mesh rather than the analytic sphere, so the
                areal density stays invariant when the continuum resolution changes. Ignored, and must
                be 1, when the basis is ``explicit``.

        Returns:
            The resolved count, rounded half-up. Zero is a legitimate answer.

        Raises:
            ValueError: on a non-positive or non-finite support, or an ``explicit`` basis given a
                support other than 1.
        """
        if not np.isfinite(support) or support <= 0.0:
            raise ValueError(f"support must be finite and positive; got {support!r}")
        if self.basis == "explicit":
            if support != 1.0:
                raise ValueError("an explicit count takes no support; pass support=1.0")
            return int(round(self.value))
        return int(np.floor(self.value * float(support) + 0.5))


@dataclass(slots=True)
class BondFamily:
    """One kinetic family and its bond population, claimed out of an arena.

    Attributes:
        name: the family's name, e.g. ``"alpha_actinin"``.
        chemistry_card: the identity that actually distinguishes families. Two declared edges sharing a
            card are ONE family; the declared component pair never was the identity.
        count: the :class:`BondCount` the population was resolved from — kept so the artifact carries
            the provenance rather than only the number.
        support: the measured support the count was resolved against.
        bonds: the arena claim for this family's BOND range.
        node_i / node_j: ``(B,)`` GLOBAL node indices. The component pair is DERIVED from these.
        rest_um / stiffness_pn_per_um: ``(B,)`` per-bond rest length and stiffness.
    """

    name: str
    chemistry_card: str
    count: BondCount
    support: float
    bonds: Claim
    node_i: npt.NDArray[np.int64]
    node_j: npt.NDArray[np.int64]
    rest_um: npt.NDArray[np.float64]
    stiffness_pn_per_um: npt.NDArray[np.float64] = field(repr=False)
    #: ⚠ Whether ``node_j`` is TOPOLOGY or merely the AS-BUILT value of a partner that is STATE.
    #:
    #: Raised by the crossbridge family on 2026-08-21, and it is a real limit of this dataclass rather
    #: than a labelling preference. ``node_i``/``node_j`` are fixed at construction, which is correct
    #: for a crosslink or a tether — those join the same two nodes for as long as they exist. **A
    #: crossbridge does not.** The head detaches, the filament slides, and it re-attaches at a
    #: DIFFERENT site: its actin partner is state that the kinetics owns, not a build-time address.
    #:
    #: The port source has been bitten by exactly this and fixed it: ``components/motor/
    #: INTEGRATION.md:142`` records that the first I3 source *"pulled the head toward a FIXED anchor
    #: and never READ the abscissa"*, so the runtime realised a passive spring network plus KMC
    #: bookkeeping while the ACTIVE directed contraction — the whole point — was never generated, and
    #: the host oracles missed it because the kernels do not launch on the dev machine. Storing a
    #: static ``node_j`` here without saying so would re-introduce a defect this project has already
    #: paid for once.
    #:
    #: ``True`` therefore obliges ``live_partner_owner`` to name what holds the live value, so an
    #: artifact can never show ``node_j`` as if it were the partner in force at that step.
    partner_is_dynamic: bool = False
    #: Who owns the live partner when ``partner_is_dynamic``. Required in that case, empty otherwise.
    live_partner_owner: str = ""

    @property
    def n_bonds(self) -> int:
        """Number of bonds in this family."""
        return int(self.node_i.size)

    def component_pairs(self, arena: WorldArena) -> dict[tuple[str, str], int]:
        """Derive which population pairs this family actually joins, and how many bonds join each.

        This is the query that replaces a declared ``component_a``/``component_b``. A family may join
        more than one pair — that is a fact about the build, not a violation — and an INTERNAL bond
        simply reports the same population twice, with no separate declaration needed.
        """
        out: dict[tuple[str, str], int] = {}
        for i, j in zip(self.node_i.tolist(), self.node_j.tolist()):
            a = arena.population_of(Kind.NODE, i) or "<unclaimed>"
            b = arena.population_of(Kind.NODE, j) or "<unclaimed>"
            key = (a, b) if a <= b else (b, a)
            out[key] = out.get(key, 0) + 1
        return out

    def provenance_row(self) -> dict[str, object]:
        """The artifact row for this family — count, basis, support, scope and source, together.

        Together on purpose: a count separated from its scope is how "235/µm²" came to look like a
        measurement.
        """
        return {
            "family": self.name, "chemistry_card": self.chemistry_card,
            "n_bonds": self.n_bonds, "basis": self.count.basis, "density": self.count.value,
            "support": self.support, "support_unit": self.count.support_unit,
            "scope": self.count.scope, "source_class": str(self.count.source_class),
            "provenance": self.count.provenance,
            # CAPACITY, never occupancy — and the row says so rather than leaving n_bonds to be read
            # as "how many are attached", which is what a kinetic family's reader will assume.
            "n_bonds_is": "CAPACITY" if not self.partner_is_dynamic else
                          "CAPACITY; occupancy is emergent and lives with the kinetics",
            "partner_is_dynamic": self.partner_is_dynamic,
            "live_partner_owner": self.live_partner_owner or None,
            "node_j_is": ("the AS-BUILT partner only; the partner in force at any step is state, "
                          f"owned by {self.live_partner_owner}" if self.partner_is_dynamic
                          else "topology, fixed for the life of the bond"),
        }

    def __post_init__(self) -> None:
        """Refuse a dynamic-partner family that does not say where the live partner lives.

        Raises:
            ValueError: if ``partner_is_dynamic`` and ``live_partner_owner`` is empty, or if
                ``live_partner_owner`` is set on a static family.
        """
        if self.partner_is_dynamic and not str(self.live_partner_owner).strip():
            raise ValueError(
                f"family {self.name!r} declares a dynamic partner but names no owner for the live "
                "value. Without one, node_j reads as the partner in force at that step — which is "
                "how components/motor/INTEGRATION.md:142 records the first I3 source realising a "
                "passive spring network while the active contraction was never generated."
            )
        if not self.partner_is_dynamic and str(self.live_partner_owner).strip():
            raise ValueError(
                f"family {self.name!r} names a live-partner owner but does not declare the partner "
                "dynamic. One of the two is wrong and guessing which would settle it silently."
            )


def build_bond_family(
    arena: WorldArena,
    name: str,
    *,
    chemistry_card: str,
    count: BondCount,
    support: float,
    pairs: npt.ArrayLike,
    rest_um: npt.ArrayLike,
    stiffness_pn_per_um: float,
    partner_is_dynamic: bool = False,
    live_partner_owner: str = "",
) -> BondFamily:
    """Claim a BOND range and build a family into ``arena``.

    Args:
        arena: the world to claim from.
        name: the family name.
        chemistry_card: the card that identifies this family's chemistry.
        count: the resolved-from :class:`BondCount`. **Required** — a family cannot be built without
            answering how many, against what, for which cell.
        support: the MEASURED support to resolve ``count`` against.
        pairs: ``(P, 2)`` candidate node index pairs, GLOBAL. The first ``count.resolve(support)`` of
            them are taken, in the order given, so selection is the caller's deterministic choice and
            never a random draw here.
        rest_um: ``(P,)`` or scalar rest length per candidate pair.
        stiffness_pn_per_um: per-bond stiffness.

    Returns:
        The built :class:`BondFamily`.

    Raises:
        ValueError: if the resolved count exceeds the candidate pairs available, if a pair is a
            self-loop, or if a rest length or stiffness is negative or non-finite.
        AssertionError: if any bond addresses a node outside the arena's live NODE prefix.
    """
    resolved = count.resolve(support)
    p = np.asarray(pairs, np.int64).reshape(-1, 2)
    if resolved > p.shape[0]:
        raise ValueError(
            f"{name}: the count resolves to {resolved} bonds but only {p.shape[0]} candidate pairs "
            "exist. Supply more candidates rather than lowering the density — the density is a "
            "physiological axis and the candidate list is a build detail."
        )
    p = p[:resolved]

    r = np.asarray(rest_um, np.float64).reshape(-1)
    if r.size == 1:
        r = np.full(resolved, float(r[0]))
    else:
        r = r[:resolved]
    if r.size != resolved:
        raise ValueError(f"{name}: rest_um has {r.size} entries for {resolved} bonds")
    if not np.all(np.isfinite(r)) or np.any(r < 0.0):
        raise ValueError(f"{name}: rest_um must be finite and nonnegative")
    if not np.isfinite(stiffness_pn_per_um) or stiffness_pn_per_um < 0.0:
        raise ValueError(f"{name}: stiffness must be finite and nonnegative")
    if resolved and np.any(p[:, 0] == p[:, 1]):
        raise ValueError(f"{name}: a bond joining a node to itself has no direction")

    live = arena.n_live(Kind.NODE)
    if resolved and (int(p.min()) < 0 or int(p.max()) >= live):
        raise AssertionError(
            f"{name}: bond addresses node [{int(p.min())}, {int(p.max())}] outside the live NODE "
            f"prefix [0, {live}). A bond to an unclaimed node is a connection to nothing."
        )

    claim = arena.claim(name, Kind.BOND, resolved) if resolved else Claim(name, Kind.BOND, 0, 0)
    return BondFamily(
        name=name, chemistry_card=chemistry_card, count=count, support=float(support), bonds=claim,
        node_i=p[:, 0].copy(), node_j=p[:, 1].copy(), rest_um=r,
        stiffness_pn_per_um=np.full(resolved, float(stiffness_pn_per_um)),
        # ⚠ These were missing until 2026-08-21 and the first kinetic family had to reach them with
        # `dataclasses.replace`. That WORKS — replace re-runs `__post_init__`, so the pairing
        # validator still fires — but it means a builder could hand back a family whose dynamic
        # partner was never declared, and the validator would never see the omission because there
        # was nothing to omit. Session A ratcheted its own workaround rather than this gap; the gap
        # is closed here so the ratchet can come off.
        partner_is_dynamic=bool(partner_is_dynamic),
        live_partner_owner=str(live_partner_owner),
    )


def _demo() -> None:
    """Self-check: the count refuses to be a bare number, and the component pair is derived."""
    from aleph.world.strand import build_strand

    arena = WorldArena(capacity={Kind.NODE: 500, Kind.SEGMENT: 500, Kind.ANGLE3: 500, Kind.BOND: 500})
    a = build_strand(arena, "cortex", start=(0, 0, 0), direction=(1, 0, 0), contour_um=1.0, seg_um=0.1)
    b = build_strand(arena, "cortex", start=(0, 0.05, 0), direction=(1, 0, 0), contour_um=1.0, seg_um=0.1)

    # A count cannot exist without scope, provenance, or a legal basis.
    for kwargs, expect in (
        ({"scope": ""}, "scope is required"),
        ({"provenance": "  "}, "provenance is required"),
        ({"basis": "per_angstrom"}, "basis must be one of"),
        ({"value": -1.0}, "finite and nonnegative"),
    ):
        base = dict(basis="areal", value=10.0, scope="MCF7 interphase adherent 37C",
                    source_class=SourceClass.PI_GAP, provenance="no Contract-Graph datum")
        try:
            BondCount(**{**base, **kwargs})
        except ValueError as exc:
            assert expect in str(exc), (kwargs, exc)
        else:  # pragma: no cover
            raise AssertionError(f"BondCount({kwargs}) must refuse")

    # An explicit count takes no support; a density does, and it must be measured.
    per_fil = BondCount(basis="per_filament", value=20.0, scope="MCF7 interphase adherent 37C",
                        source_class=SourceClass.CONVENIENCE,
                        provenance="density_per_fil = 20 (spanning threshold, not sourced)")
    assert per_fil.resolve(2.0) == 40
    assert per_fil.support_unit == "filaments"
    try:
        per_fil.resolve(0.0)
    except ValueError as exc:
        assert "finite and positive" in str(exc)
    else:  # pragma: no cover
        raise AssertionError("a zero support must refuse")

    # A density low enough to give no bonds is an ANSWER, not an error.
    sparse = BondCount(basis="areal", value=1e-6, scope="MCF7", source_class=SourceClass.PI_GAP,
                       provenance="deliberately sparse, for the boundary case")
    assert sparse.resolve(1.0) == 0

    ai = a.nodes.lo + np.arange(a.n_nodes)
    bi = b.nodes.lo + np.arange(b.n_nodes)
    cand = np.stack([ai, bi], axis=1)

    fam = build_bond_family(
        arena, "alpha_actinin", chemistry_card="alpha_actinin_ferrer2008", count=per_fil, support=0.5,
        pairs=cand, rest_um=0.05, stiffness_pn_per_um=4.6e5)
    assert fam.n_bonds == 10, "20 per filament x 0.5 filaments of support"
    assert fam.bonds.count == fam.n_bonds

    # The pair is DERIVED. Both strands are cortex, so this reads as an internal bond with no
    # separate INTERNAL declaration anywhere.
    assert fam.component_pairs(arena) == {("cortex", "cortex"): 10}

    row = fam.provenance_row()
    assert row["scope"] and row["provenance"] and row["support_unit"] == "filaments"
    assert row["n_bonds"] == 10

    # A bond to a node nobody claimed is a connection to nothing.
    try:
        build_bond_family(arena, "ghost", chemistry_card="none", count=per_fil, support=0.05,
                          pairs=np.array([[0, arena.n_live(Kind.NODE) + 5]]), rest_um=0.05,
                          stiffness_pn_per_um=1.0)
    except AssertionError as exc:
        assert "connection to nothing" in str(exc)
    else:  # pragma: no cover
        raise AssertionError("an out-of-range bond must raise")

    print("bond self-check OK —", row)


def _demo_dynamic_partner() -> None:
    """The 2026-08-21 additions: capacity-vs-occupancy, and a partner that is state.

    Kept separate from :func:`_demo` so the original self-check stays readable as the thing it was
    written for, and so a failure here points at the kinetic-family additions specifically.
    """
    import numpy as _np

    from aleph.world.arena import Claim, Kind

    common = dict(
        name="nmii_head_actin", chemistry_card="nmii_head_actin_catchslip_kovacs2007",
        count=BondCount(basis="per_filament", value=2.0, scope="MCF7 cortical shell, resting",
                        source_class=SourceClass.PI_GAP, provenance="self-check"),
        support=1.0, bonds=Claim("nmii_head_actin", Kind.BOND, 0, 2),
        node_i=_np.array([0, 1], dtype=_np.int64), node_j=_np.array([10, 11], dtype=_np.int64),
        rest_um=_np.array([0.01, 0.01]), stiffness_pn_per_um=_np.array([1.0, 1.0]))

    # 1. A dynamic partner with no owner is refused. Without one, node_j reads as the partner in
    #    force at that step, which is the defect INTEGRATION.md:142 records being paid for once.
    try:
        BondFamily(**common, partner_is_dynamic=True)
    except ValueError as exc:
        assert "names no owner" in str(exc), exc
    else:
        raise AssertionError("a dynamic partner with no owner must be refused")

    # 2. ...and an owner without the flag is refused too. One of the two is wrong and guessing which
    #    would settle it silently — the failure mode this whole class exists to stop.
    try:
        BondFamily(**common, live_partner_owner="CrossbridgeState.bound")
    except ValueError as exc:
        assert "does not declare the partner dynamic" in str(exc), exc
    else:
        raise AssertionError("an owner without the flag must be refused")

    # 3. Declared properly, the artifact row cannot be misread.
    fam = BondFamily(**common, partner_is_dynamic=True,
                     live_partner_owner="aleph.world.families.nmii_cortex_crossbridge.CrossbridgeState.bound")
    row = fam.provenance_row()
    assert row["partner_is_dynamic"] is True
    assert "AS-BUILT" in row["node_j_is"] and "state" in row["node_j_is"]
    assert row["n_bonds_is"].startswith("CAPACITY")
    assert "occupancy is emergent" in row["n_bonds_is"]

    # 4. A static family says so, and says nothing about occupancy it does not have.
    stat = BondFamily(**common)
    assert stat.provenance_row()["node_j_is"] == "topology, fixed for the life of the bond"
    assert stat.provenance_row()["n_bonds_is"] == "CAPACITY"

    print("bond dynamic-partner self-check OK — capacity and occupancy cannot be confused in a row")


def _demo_builder_forwards_dynamic_partner() -> None:
    """`build_bond_family` carries the dynamic-partner declaration through, PI-adjacent 2026-08-21.

    Added because session A had to reach the fields with ``dataclasses.replace`` — which is correct
    and re-runs the validator, but leaves a builder able to hand back an UNDECLARED family, and the
    validator cannot object to an omission that was never expressible.
    """
    import numpy as _np

    from aleph.world.arena import Kind, WorldArena

    # BOND capacity for three families of two: the refusals still CLAIM before they validate,
    # so a tight capacity makes this self-check fail on the arena rather than on the point.
    arena = WorldArena(capacity={Kind.NODE: 8, Kind.BOND: 8}, device=None)
    arena.claim("p", Kind.NODE, 8)
    kw = dict(chemistry_card="c", support=1.0, pairs=_np.array([[0, 1], [2, 3]]),
              rest_um=_np.array([0.01, 0.01]), stiffness_pn_per_um=1.0,
              count=BondCount(basis="explicit", value=2.0, scope="self-check",
                              source_class=SourceClass.PI_GAP, provenance="self-check"))

    fam = build_bond_family(arena, "dyn", **kw, partner_is_dynamic=True,
                            live_partner_owner="X.bound")
    assert fam.partner_is_dynamic and fam.live_partner_owner == "X.bound"
    assert "AS-BUILT" in fam.provenance_row()["node_j_is"]

    # The pairing validator still fires through the builder, not only through replace().
    for bad, expect in (({"partner_is_dynamic": True}, "names no owner"),
                        ({"live_partner_owner": "X.bound"}, "does not declare the partner dynamic")):
        try:
            build_bond_family(arena, "bad", **kw, **bad)
        except ValueError as exc:
            assert expect in str(exc), (bad, exc)
        else:
            raise AssertionError(f"build_bond_family{bad} must be refused")

    print("bond builder self-check OK — the declaration survives the builder, and so does the pairing")


if __name__ == "__main__":
    # ⚠ At the END, and all three run. Until 2026-08-21 this block sat above the two demos below,
    # so `python -m aleph.world.bond` bound them but never called them — they existed only under
    # pytest. Found by the owner of `observe_gamma.py` after the same shape shipped there and broke
    # its documented invocation outright; an AST sweep of `world/` and `laws/` found this as the only
    # other instance. A module's own entry point must be able to reach its whole body, and running
    # only the first of three self-checks is the quiet version of not running them.
    _demo()
    _demo_dynamic_partner()
    _demo_builder_forwards_dynamic_partner()
