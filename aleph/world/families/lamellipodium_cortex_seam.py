r"""``lamellipodium_cortex_seam`` — the same family as ``sf_cortex_transient``, and still unanswerable.

VERDICT: **DUPLICATE**, and separately **BLOCKED**. Two findings, and neither replaces the other. This
module builds nothing: the seam is the *same physics* as ``sf_cortex_transient`` and therefore may not
become a second :class:`~aleph.world.bond.BondFamily`; and the family they share cannot be counted, so
the value it defers to is blocked as well. Folding a connector into another does not answer the
question — it moves it, and this module's job is to say where it moved to.

WHY IT IS A DUPLICATE, ARGUED FROM PHYSICS
------------------------------------------
1. **One force law.** Both edges resolve to ``engine.filament_crosslink.FilamentCrosslinkConnector`` —
   a state-gated Hookean central force between two barycentrically-interpolated segment points, in two
   components owning disjoint arrays. That module's own docstring: *"Five declared edges across three
   families are the SAME mechanics … The family label is biology; the force law is one."*
2. **One chemistry.** ``contracts.py:734`` (``sf_cortex_transient``) and ``contracts.py:798`` (this
   edge) declare the identical ``chemistry_card="transient_actin_crosslink"``, and a committed test
   asserts the two cards equal. ``bond.py`` makes the card the identity in as many words: *"Two
   declared edges sharing a card are ONE family; the declared component pair never was the identity."*
3. **One rate law.** The card's kinetics is α-actinin bind/unbind (Bell slip bond) over the SOURCED
   ``ALPHA_ACTININ`` record. Nothing in it is stress-fibre-specific: α-actinin crosslinks actin to
   actin, and a lamellipodial filament and a cortical filament are both actin.
4. **The only remaining difference is not stored.** The two edges differ solely in which population
   the actin-side endpoint lives in — and a bond holds two GLOBAL node indices and nothing else, with
   the pair derived from the arena's ID ranges at query time. A second family would therefore differ
   in a field that does not exist, which is the ``bond.py`` argument for why "membrane-side NMII" and
   "SF-side NMII" are not two families.

⚠ **THE FOLD IS NOT BY NAME, AND HERE IS THE CONTROL.** The sibling ``filopodium_cortex_root`` carries
the SAME ``ConnectorFamily.TRANSIENT_ACTIN`` and the same ``* -> cortex`` shape, and it is **NOT** a
duplicate: its card is ``formin_fascin_root_coupling``, a formin-nucleated bundle root, not an
α-actinin crosslink. Folding by family label would have swallowed it; folding by chemistry does not.
A rule that folds everything is not a discriminator, and this one refuses a case.

WHAT SURVIVES THE FOLD — the count, and it is the ERM class again
-----------------------------------------------------------------
The prompt for this lane asks whether this connector hides a population behind one declaration, as the
membrane–cortex ERM edge did. **It does, in the same class and a milder key.** The incumbent's count for
this chemistry is not a mesh vertex count; it is a PLACEMENT count. ``engine/sf_cortex_transient.py``
records it against itself, in capitals: *"THE COUNT IS PLACEMENT-DERIVED, NOT A SOURCED INVENTORY"* —
the joints are whatever geometry happened to fall inside the sourced α-actinin capture radius
ε = 0.06 µm, and at full native (job 90) the SF→cortex median separation was 2.59 µm, **43× that
radius**. The same module refuses a multiplicity above one because *"multiplicity > 1 needs a sourced
crosslinker-per-segment density; there is none"*, and ``alpha_actinin_kinetics.py`` closes with *"How
MANY crosslinks there should be is not here and is not sourced."*

So the ERM sentence transposes exactly: a mesh number wore a physiological label there; **a placement
number wears one here.** The capture radius is SOURCED and must never be widened to produce pairs —
widening a sourced constant until the geometry cooperates is tuning to an outcome — which leaves the
count with no source at all rather than a bad one.

For THIS edge the gap is compounded, because the lamellipodial side is itself a gap:
``world/build/lamellipodium.py`` calls its own count *"the worst of the five"* — 200 filaments is an
EXPLICIT number from ``laws/architecture_spec``, over its own 8 × 8 µm patch, and *"there is no
lamellipodial actin AREAL density in the KB at all."* PHASE 1 built exactly those 200 strands. A
per-filament density resolved against 200 filaments would multiply one unanswered question by another
and return a number that looks measured.

⚠ AND THE SEAM HAS NO GEOMETRY YET. ``build_lamellipodium`` takes its patch origin from the caller and
refuses to derive it, because where a leading edge sits depends on the adhesion state and *"deriving it
from a cell radius would be inventing a spread state."* Until a spread state is chosen, the
lamellipodium and the cortex stand in no determined spatial relation, so the number of pairs inside any
capture radius is a free consequence of an unmade decision — which is precisely the failure the
43× measurement above already demonstrated for the stress fibres.

⚠ THIS IS ALSO THE FIRST FAMILY WITH KINETICS — PI, NOT US
-----------------------------------------------------------
``transient_actin_crosslink`` binds and unbinds by construction; the incumbent refuses to bind the
mechanical edge without a rate law because *a transient crosslink that can bind and never unbind is a
permanent weld*, which the charter forbids by name. ``bond.py`` deliberately carries no kinetics, no
attach/detach, no free list and no snapshot twins, deferring them to *"the first family that has
kinetics"*. **That family is this one.** Its shape is a contract decision and is escalated, not built
here — recorded as ``kinetic-runtime-shape`` in :data:`BLOCKED_BY`.

WHAT WOULD CLOSE THIS
---------------------
Answers, in the order they unblock each other, all of them for the ONE folded family:
  1. the run's target cell line and state (nothing below is scoped without it);
  2. a transient-actin-crosslinker inventory for that cell — and the BASIS is part of the question: a
     seam is a boundary, so ``areal`` over a contact area and ``per_filament`` over the actin-side
     filaments are different physical claims, not two spellings;
  3. whether one density covers both the SF→cortex and lamellipodium→cortex candidate sets, since one
     :class:`~aleph.world.bond.BondCount` carries one value and one support — if the seam density
     differs from the SF density, the fold is right about the chemistry and wrong about the count, and
     that is a third finding rather than a reason to unfold;
  4. the spread state that puts the lamellipodial patch somewhere, so a support can be MEASURED;
  5. the kinetic-runtime shape, from PI.

Nothing here is a placeholder awaiting a default. Each is a question with no answer on disk today.

Sanity Gate (recorded per CLAUDE.md, before first execution)
  * dimensional — no quantity is formed. Were it buildable: ``rest_um`` [µm], ``stiffness_pn_per_um``
    [pN/µm], and a ``per_filament`` count dimensionless against a filament count. The two mechanical
    magnitudes ARE sourced (Ferrer 2008 PNAS AFM, PI-approved 2026-06-30); only the count is absent,
    which is why a partial provenance may not be reported as a provenance.
  * boundary — the refusal is unconditional and takes no arguments that could switch it off. It is not
    reachable-but-rare: every call raises, and the self-check asserts that.
  * conservation/invariant — nothing is claimed from the arena, so an arena passed in leaves with its
    live counts unchanged. The self-check asserts that too: a blocked family that consumed capacity
    would be a leak wearing a refusal.
  * CFL/precision — no integration, no array, no device. Host bookkeeping only.
  * sign sense — not applicable; no force law is written here, deliberately (see the kinetics note).
  * measurement protocol — no device is touched, nothing is uploaded, nothing is read back.

engine units: length µm, stiffness pN/µm.  Runtime: pure host bookkeeping, CPU-importable.
"""

from __future__ import annotations

from typing import NoReturn

from aleph.world.families import ConnectorGapError, FamilySpec

__all__ = ["CONNECTOR", "VALUE_CONNECTOR", "DUPLICATE_OF_ENTRY", "BLOCKED_BY", "SPEC",
           "build_lamellipodium_cortex_seam"]

#: The declared connector this module answers for, exactly as the device-run census spells it.
CONNECTOR: str = "lamellipodium_cortex_seam"

#: The connector this one is the same physics AS, and whose family carries both candidate sets.
VALUE_CONNECTOR: str = "sf_cortex_transient"

#: For ``families.DUPLICATE_OF``. Exposed rather than inserted: ``__init__.py`` is shared by fourteen
#: concurrent sessions and file separation is the only protection they have, so Lead merges this.
DUPLICATE_OF_ENTRY: tuple[str, str] = (CONNECTOR, VALUE_CONNECTOR)

#: The unanswered questions, by name, for Lead's decision queue. Ordered as they unblock each other.
BLOCKED_BY: tuple[str, ...] = (
    "target-cell-scope",
    "transient-actin-crosslinker-inventory",
    "count-basis-areal-vs-per-filament",
    "one-density-for-both-candidate-sets",
    "lamellipodial-spread-state-for-a-measurable-support",
    "kinetic-runtime-shape",
)

SPEC = FamilySpec(
    connector=CONNECTOR,
    populations=("lamellipodium", "cortex"),
    basis="per_filament",
    blocked_by=BLOCKED_BY,
    source=(
        "DUPLICATE of sf_cortex_transient by chemistry card 'transient_actin_crosslink' "
        "(contracts.py:734 vs :798), one force law (filament_crosslink), one rate law "
        "(alpha_actinin_kinetics). Mechanics SOURCED — ALPHA_ACTININ.link_k 4.6e5 pN/um and "
        "capture_radius_um 0.06 (Ferrer 2008 PNAS AFM, PI-approved 2026-06-30). COUNT: none. The "
        "incumbent's is placement-derived inside the capture radius and says so; the lamellipodial "
        "side is a PI-GAP of its own (no areal actin density in the KB). Basis 'per_filament' is the "
        "one the incumbent's refusal names, NOT a ratified choice — see count-basis-areal-vs-per-filament."
    ),
)


def build_lamellipodium_cortex_seam(*_args: object, **_kwargs: object) -> NoReturn:
    """Refuse to build, twice over: this is a duplicate family AND an uncounted one.

    Takes and ignores any arguments on purpose — there is no argument that could make it succeed, and a
    signature that looked configurable would invite someone to pass the missing density in.

    Raises:
        ConnectorGapError: always. The refusal is this module's product.
    """
    raise ConnectorGapError(
        CONNECTOR,
        list(BLOCKED_BY),
        note=(
            f"and before any of that: this edge is the SAME PHYSICS as {VALUE_CONNECTOR} (identical "
            "chemistry card, force law and rate law; the population pair is derived from arena ID "
            "ranges and is not stored), so it must NOT become a second BondFamily. The count question "
            "moves to the folded family and is not answered by the fold. See this module's docstring."
        ),
    )


def _demo() -> None:
    """Self-check: the refusal fires, the fold is recorded, and nothing is consumed to produce it."""
    assert not SPEC.buildable, "a spec with unanswered questions may never report buildable"
    assert SPEC.blocked_by == BLOCKED_BY
    assert SPEC.connector == CONNECTOR

    # The fold points somewhere else, and never at itself: a self-referential duplicate is a family
    # that quietly builds while claiming to defer.
    assert DUPLICATE_OF_ENTRY[0] == CONNECTOR
    assert DUPLICATE_OF_ENTRY[1] != CONNECTOR, "a connector cannot be a duplicate of itself"

    # Every call refuses, with or without arguments, and names every open question.
    for args, kwargs in (((), {}), ((object(),), {"density_per_filament": 20.0, "support": 200.0})):
        try:
            build_lamellipodium_cortex_seam(*args, **kwargs)
        except ConnectorGapError as exc:
            assert exc.connector == CONNECTOR
            assert set(exc.missing) == set(BLOCKED_BY), exc.missing
            assert VALUE_CONNECTOR in str(exc), "the refusal must name what it defers to"
        else:  # pragma: no cover - the refusal is unconditional
            raise AssertionError(f"{CONNECTOR} must refuse; it built with {kwargs!r}")

    # A blocked family that consumed arena capacity would be a leak wearing a refusal. Only run this
    # where an arena can be made without a device; the refusal above is the part that must always hold.
    try:
        from aleph.world.arena import Kind, WorldArena
    except Exception:  # pragma: no cover - import environment, not a branch of this module
        pass
    else:
        arena = WorldArena(capacity={Kind.NODE: 8, Kind.BOND: 8})
        before = {k: arena.n_live(k) for k in (Kind.NODE, Kind.BOND)}
        try:
            build_lamellipodium_cortex_seam(arena)
        except ConnectorGapError:
            pass
        assert {k: arena.n_live(k) for k in (Kind.NODE, Kind.BOND)} == before, "refusal must not claim"

    print(f"{CONNECTOR} self-check OK — DUPLICATE_OF {VALUE_CONNECTOR}; "
          f"{len(BLOCKED_BY)} unanswered: {', '.join(BLOCKED_BY)}")


if __name__ == "__main__":
    _demo()
