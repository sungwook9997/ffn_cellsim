r"""``if_sf_plectin`` — the plectin IF↔actin crosslink. **BLOCKED: it does not build.**

The declaration is ``engine/contracts.py:764``::

    ConnectorContract("if_sf_plectin", ConnectorFamily.PLECTIN,
                      "intermediate_filament", "sf_arc", kinetics=True, commit_on_accept=True,
                      chemistry_card="plectin_if_actin")

and the committed device-run census records it ``NOT_BOUND`` / ``n_calls: 0`` / *"no runtime object in
this composed world"*.  So this module is a first implementation, and its product is a REFUSAL: five
questions have no answer, and :class:`~aleph.world.bond.BondCount` exists precisely so that inventing
one is not available.

IT IS NOT A DUPLICATE, AND THE ARGUMENT IS FROM PHYSICS.
    ``engine/filament_crosslink.py`` establishes that FIVE declared edges — ``sf_cortex_transient``,
    ``lamellipodium_cortex_seam``, ``filopodium_cortex_root``, ``if_sf_plectin``,
    ``mt_sf_spectraplakin`` — are the SAME MECHANICS: one state-gated Hookean central force between two
    interpolated segment points in two disjoint arrays.  **That is not a reason to fold them here.**
    ``bond.py`` puts family identity on the ``chemistry_card``, not on the force law and not on the
    declared component pair, and it is right to: the five edges declare FOUR cards, because chemistry
    sets *when a bond forms and breaks*, which is the state the family owns.  ``plectin_if_actin`` is
    the only edge carrying its card.  Folding by shared force law would be the same mistake as folding
    by shared name, made one level deeper — it would put plectin's on/off kinetics and Arp2/3-era
    transient-actin kinetics into one population that cannot then be given two rate laws.

    THE NEAR MISS, recorded because it is the trap this package warns about.  ``if_nucleus_linc``
    carries the card ``nesprin3_plectin_if_linc`` — the word *plectin* is in it, and it joins the SAME
    ``intermediate_filament`` population.  It is a different family: nesprin-3 at the nuclear envelope
    is a LINC joint with its own (also unsourced) rate law and its own partner, and the ``PLECTIN``
    vs ``LINC`` families are declared apart.  Judged by physics; the shared substring is noise.

    A THIRD plectin population EXISTS and is not declared as a connector at all.
    ``laws/intermediate_filaments.py:143`` builds tangential IF↔IF plectin crosslinks inside the cage.
    Under ``bond.py``'s derived-pair rule those could share this family's arena range and simply report
    ``("intermediate_filament", "intermediate_filament")`` from
    :meth:`~aleph.world.bond.BondFamily.component_pairs` — but plectin's IF-IF and IF-actin binding
    are different domains and would be a different card, so this is flagged for the PI rather than
    assumed either way.  It is not in ``DUPLICATE_OF``: that map is for DECLARED connectors.

⚠ **THE ERM DEFECT IS PRESENT HERE TWICE, NOT ONCE.**
    ``bond.py`` names the class of error: *a mesh number wearing a physiological label*.  For ERM the
    count fell back to one tether per membrane vertex.  For plectin BOTH the count and the pairing rule
    are discretisation numbers today:

    * **count** — ``laws/intermediate_filaments.py:143`` is *"one per bead to the nearest bead on a
      DIFFERENT spoke within reach"*.  One per BEAD.  The bead count is ``span / l_seg``, so the
      crosslink count is a function of ``seg_um`` and of nothing biological.  :func:`_demo` runs this
      rule at two resolutions and shows the count change by the resolution ratio, for the same cell.
    * **reach** — ``KB-DRAFT-3.B-10`` records ``crosslink_reach = ~l_seg``, and the law uses
      ``1.5 * l_seg_um``.  The reach is the SEGMENT LENGTH, so the candidate-pair set is set by the
      discretisation as well.  A count resolved against a candidate list that is itself a mesh artefact
      is not saved by having a density.

WHAT THE KNOWLEDGE BASE ACTUALLY HOLDS (queried, not recalled — ``outputs/tag_kb/kb.duckdb``).
    One row mentions this chemistry: ``KB-DRAFT-3.B-10`` *"Cytolinker crosslinks with correct identity
    (plectin IF-IF, IF-actin, IF-MT)"*, ``status: draft``, alias ``AUDIT-2026-07-15 gap:``, assumption
    *"DRAFT (NOT-IN-TAG). PROVISIONAL id — PI assigns final KB number at ingest"*.  Its whole
    ``value_range_si`` is ``ratio_xl (crosslink/backbone) = 0.1 (placeholder — PI-anchor needed)`` and
    ``crosslink_reach = ~l_seg``.  **A stiffness ratio and a reach — no count, and no density.**  The
    same row's own subtopic states the IF compartment has *"ZERO materialized KnowledgeClaim,
    SourceEvidence, ModelContract or ValidationGate rows"*.  ``PARAM-k_xl`` is ``default: TBD``,
    ``status: draft``.  Its two citations (Svitkina/Verkhovsky/Borisy 1996 JCB 135(4):991; Wiche 2015
    Curr Opin Cell Biol 32:21) are a fibroblast-lamella EM study and a review — neither is a count, and
    neither is in the target cell or state, so even a number read off them would be
    :attr:`~aleph.world.bond.SourceClass.UNRATIFIED_PROXY`, not ``SOURCED``.

⚠ **THE KINETICS QUESTION IS THE ONE TO ESCALATE, AND IT BLOCKS EVEN A COUNTED FAMILY.**
    The contract declares ``kinetics=True, commit_on_accept=True``.  ``bond.py`` has NO kinetics on
    purpose and defers them to *"the first family that has kinetics"*.  ``filament_crosslink.py`` is
    explicit about what a plectin edge without a rate law is: *"a transient crosslink with no rate law
    is not a transient crosslink; it is a permanent weld, which is the one thing the charter says a
    connector may never be"* — and its runtime RAISES on ``propose_events`` for exactly this edge
    rather than proposing nothing (``test_filament_crosslink.py:126``).  So a mechanics-only
    ``BondFamily`` here would be a weld even if every count question were answered.  Per this
    package's rule, that is carried to the PI and not built.

POPULATION NAMING, recorded rather than reconciled.  The contract's endpoint B is ``sf_arc``; the arena
population built by ``world/build/stress_fiber.py`` and claimed in PHASE 1 is ``stress_fiber`` (1,000
strands / 201,000 nodes).  There is exactly one SF population in the arena, so the endpoint is not
ambiguous, and :data:`SPEC` records the arena's name as ``FamilySpec`` requires.  The rename is not a
blocker; it is noted so nobody re-derives it.

Sanity Gate (recorded per CLAUDE.md, before first execution).  Nothing in this module executes physics —
the only reachable path is the refusal — so the gate records what a build WOULD have to pass, and the
one thing that does run:
  * dimensional — a plectin count is dimensionless; its density is [1/µm³] on a volumetric basis or
    dimensionless per filament.  ``rest_um`` [µm], ``stiffness_pn_per_um`` [pN/µm].  The only magnitude
    the KB offers, ``ratio_xl = 0.1``, is DIMENSIONLESS and multiplies ``k_bb`` — it is not a count and
    cannot be converted into one.
  * boundary — a resolved count of zero would be a legitimate answer (``bond.py``); a count of zero
    reached because nobody asked is what this module refuses instead.
  * conservation/invariant — the force law is central and equal-and-opposite by construction in the
    runtime that will carry it; nothing is asserted here because nothing is built here.
  * CFL/precision — no integration; float64.  ``_demo``'s bead arithmetic is integer.
  * sign sense — a plectin crosslink is a two-sided Hookean tie, unlike a steric or ratchet contact:
    it may pull AND push, so no unilateral gate belongs on it.  Recorded now so the eventual build does
    not acquire one by copying a contact edge.
  * measurement protocol — host-side; no device is touched, nothing is uploaded, nothing is read back.

engine units: length µm, stiffness pN/µm.  Runtime: pure host bookkeeping, CPU-importable.
"""

from __future__ import annotations

from aleph.world.arena import WorldArena
from aleph.world.bond import BondFamily
from aleph.world.families import ConnectorGapError, FamilySpec

__all__ = ["CHEMISTRY_CARD", "SPEC", "build_if_sf_plectin"]

#: The identity that distinguishes this family — not the force law, which four other edges share, and
#: not the declared component pair, which ``bond.py`` derives from arena ID ranges at query time.
CHEMISTRY_CARD = "plectin_if_actin"

SPEC = FamilySpec(
    connector="if_sf_plectin",
    populations=("intermediate_filament", "stress_fiber"),
    # The basis is itself one of the unanswered questions; ``per_filament`` is recorded so the census
    # has an entry, and the question is named in ``blocked_by`` so it is not mistaken for a ruling.
    basis="per_filament",
    blocked_by=(
        "count: how many plectin IF-actin crosslinks a cell has. No KnowledgeClaim, Parameter or "
        "ModelContract row in kb.duckdb carries one. KB-DRAFT-3.B-10 is a draft AUDIT gap row "
        "(NOT-IN-TAG, provisional id) whose only value is ratio_xl=0.1, a stiffness ratio flagged "
        "'placeholder - PI-anchor needed', and PARAM-k_xl is default TBD.",
        "basis: whether the count is per IF filament, per um^3 of IF-SF overlap, or per um of IF "
        "contour. bond.py has no linear basis, and no source picks one, so this is a modelling "
        "decision before it is a lookup.",
        "reach: the IF-SF separation within which a pair is crosslinkable, i.e. what the candidate "
        "pair list is. KB-DRAFT-3.B-10 gives '~l_seg' and laws/intermediate_filaments.py:143 uses "
        "1.5*l_seg_um - the SEGMENT LENGTH. The candidate set is therefore set by seg_um, so a "
        "density resolved against it inherits the discretisation it was meant to be invariant to.",
        "scope: the two citations (Svitkina 1996, fibroblast lamella EM; Wiche 2015, review) are "
        "neither the target cell nor its state, so any number taken from them is UNRATIFIED_PROXY "
        "and not SOURCED - the distinction the PI's 2026-08-15 ruling exists for.",
        "kinetics: the contract declares kinetics=True/commit_on_accept=True. bond.py has no kinetics "
        "by design and defers them to 'the first family that has kinetics'. A mechanics-only plectin "
        "family is a PERMANENT WELD, which the charter forbids a connector to be, so this is a PI "
        "escalation and blocks the build even if every count question above were answered.",
    ),
)


def build_if_sf_plectin(arena: WorldArena, **_: object) -> BondFamily:
    """Refuse to build ``if_sf_plectin``, naming what is unanswered.

    The signature is the one the eventual builder will have, so that a caller wiring PHASE 3 gets the
    refusal at the call site rather than at an import.

    Args:
        arena: the world that would be claimed from. Unused: nothing is claimed.
        **_: swallowed. A ``count`` supplied here would not unblock the family — four other questions
            remain, and the fifth makes a counted-but-kinetics-less family a weld.

    Raises:
        ConnectorGapError: always, with :attr:`SPEC`'s five questions.
    """
    raise ConnectorGapError(
        SPEC.connector,
        list(SPEC.blocked_by),
        note=(
            "Not a duplicate: chemistry_card 'plectin_if_actin' is unique among the five edges that "
            "share the crosslink force law. The kinetics question is the one to rule on first - it "
            "blocks the build independently of the count."
        ),
    )


def _bead_rule_count(span_um: float, seg_um: float, n_spokes: int) -> int:
    """Reproduce the incumbent's crosslink rule, so the artefact can be MEASURED rather than asserted.

    ``laws/intermediate_filaments.py:143`` lays one crosslink per bead to the nearest bead on a
    different spoke within ``1.5 * l_seg``. Beads per spoke is ``span / seg_um``, so this returns
    ``n_spokes * n_beads`` — a count with no biological input at all. It exists to be shown wrong.
    """
    n_beads = int(span_um / seg_um) + 1
    return n_spokes * n_beads


def _demo() -> None:
    """Self-check: the module refuses, and the rule it refuses to adopt is resolution-dependent."""
    # 1. The spec is a census entry that says BLOCKED, and every question is named.
    assert not SPEC.buildable
    assert len(SPEC.blocked_by) == 5, SPEC.blocked_by
    assert SPEC.populations == ("intermediate_filament", "stress_fiber")
    assert all(q.split(":")[0] in
               {"count", "basis", "reach", "scope", "kinetics"} for q in SPEC.blocked_by)

    # 2. The refusal is typed, so a partially-built family cannot be mistaken for a built one, and a
    #    caller catching NotImplementedError still gets the questions.
    try:
        build_if_sf_plectin(None, count="a number someone typed")  # type: ignore[arg-type]
    except ConnectorGapError as exc:
        assert exc.connector == "if_sf_plectin"
        assert len(exc.missing) == 5
        assert "kinetics" in str(exc) and "PERMANENT WELD" in str(exc)
        assert isinstance(exc, NotImplementedError)
    else:  # pragma: no cover
        raise AssertionError("if_sf_plectin must refuse to build")

    # 3. THE FINDING. The same cell, the same biology, two mesh resolutions — and the incumbent's rule
    #    returns counts in the ratio of the resolutions. A number that halves when seg_um doubles is a
    #    property of the mesh, not of plectin. This is the ERM defect, in this connector, executable.
    span, n_spokes = 5.0, 60
    coarse = _bead_rule_count(span, 0.50, n_spokes)
    fine = _bead_rule_count(span, 0.25, n_spokes)
    assert coarse == 660 and fine == 1260, (coarse, fine)
    assert fine > 1.9 * coarse - n_spokes, "halving seg_um must roughly double the count"
    # ...and neither number is the answer to "how many plectin crosslinks does this cell have".

    print(f"if_sf_plectin: BLOCKED on {len(SPEC.blocked_by)} — "
          f"{', '.join(q.split(':')[0] for q in SPEC.blocked_by)}. "
          f"bead-rule count is {coarse} at seg=0.50 um and {fine} at seg=0.25 um, same cell.")


if __name__ == "__main__":
    _demo()
