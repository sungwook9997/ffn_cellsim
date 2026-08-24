r"""``sf_cortex_transient`` — the transient actin crosslink between a stress fibre and the cortex.

**This module builds nothing, and that refusal is its product.** It is NOT a duplicate: under the
criterion converged on 2026-08-20 (:data:`aleph.world.families.DUPLICATE_CRITERION` — a family's
identity is its CHEMISTRY CARD) this edge is the CANONICAL member of the card
``transient_actin_crosslink``, and ``lamellipodium_cortex_seam`` folds INTO it.  So this is where that
card's count question lives, and the count cannot be answered.

────────────────────────────────────────────────────────────────────────────────────────────────
1. IDENTITY — read off the contract table, not asserted.

``contracts.py`` is the declaration of record and it is read, never re-spelled:

  * ``:730-735``  ``sf_cortex_transient``       TRANSIENT_ACTIN  sf_arc↔cortex        ``transient_actin_crosslink``
  * ``:797-803``  ``lamellipodium_cortex_seam`` TRANSIENT_ACTIN  lamellipodium↔cortex ``transient_actin_crosslink``
  * ``:843-849``  ``filopodium_cortex_root``    TRANSIENT_ACTIN  filopodium↔cortex    ``formin_fascin_root_coupling``
  * ``:905-909``  ``dorsal_arc_crosslink``      TRANSIENT_ACTIN  sf_arc↔sf_arc        **chemistry_card is None**

**One family, two declared pairs.**  D13 (``lamellipodium_cortex_seam``) reached this from its side and
this module confirms it from the contract table: same card, so ONE family whose
``BondFamily.component_pairs`` reports ``{("cortex","stress_fiber"): n, ("cortex","lamellipodium"): m}``
at query time.  A bond stores two global node indices; the pair was never the identity.

**The control that makes the fold non-vacuous** is ``filopodium_cortex_root``: same
``ConnectorFamily.TRANSIENT_ACTIN``, same runtime class, same cortex endpoint — and a DIFFERENT card.
It does not fold.  Neither do ``if_sf_plectin`` (``plectin_if_actin``) or ``mt_sf_spectraplakin``
(``spectraplakin_mt_actin``), which share only ``FilamentCrosslinkConnector``.  Folding on a shared
runtime class would be folding by implementation — one step worse than folding by name.

⚠ **``dorsal_arc_crosslink`` HAS NO DECLARED CARD, AND IT IS ALREADY ``DEVICE_RUN``.**  On physics it is
a candidate member of this one: same ``ConnectorFamily``, actin on both ends, and its stiffness is the
same symbol this edge's builder spends — ``sf_mechanics.py:100`` sets
``ALPHA_ACTININ_K_PN_PER_UM = float(ALPHA_ACTININ.link_k)`` and both take their stiffness from it
(``sf_mechanics.py:255``, ``sf_cortex_transient.py:267``).  But its ``chemistry_card`` is ``None``, and
under the criterion an undeclared identity is not a match — it is an unanswered question.  **Asserting a
card for another session's connector is exactly the invention this package forbids**, so it is filed as
a PI item (§3b) rather than folded.  It matters for the count: if that card is this one, its population
is part of this family and the number below has to cover it.

⚠ **AND THE SF-SIDE PARTNER POPULATION MAY NOT EXIST.**  The contract's endpoint is ``sf_arc`` —
transverse and dorsal ARCS — while PHASE 1 built ``stress_fiber``, the VENTRAL FA-to-FA bundles of
``world/build/stress_fiber.py``.  Those are different structures, not a rename, and D8 and D9 hit the
same class from their own side within the hour (``test_families_census.py``, xfail-strict).  ``SPEC``
below names ``stress_fiber`` because it is the only SF-like population the arena has, **and that choice
is recorded as a blocker rather than made silently**: if the arcs are a population still to be built,
this edge's count is being asked against the wrong one, and if ``sf_arc`` is stale then the ventral
bundles inherit an endpoint role written for arcs.  Which it is, is the PI's call.

────────────────────────────────────────────────────────────────────────────────────────────────
2. THE COUNT — checked against the ERM pattern, and it is the same class of defect.

The ERM case was *"a mesh number wearing a physiological label"*: the count fell out of the icosphere
subdivision.  This one is not a mesh number.  It is a **placement number wearing a chemistry label**,
and the disguise is better, because the constant it wears is genuinely SOURCED.

  * ``pair_sf_to_cortex_segments`` takes every SF segment with a cortex segment inside 0.06 µm, one
    each.  The RADIUS is sourced (Ferrer 2008's ε = 60 nm).  The COUNT is not: it is however many pairs
    this particular placement of two independent builders happens to produce.  The engine module says so
    itself, in capitals, and records the native measurement — SF→cortex nearest distance **median
    2.59 µm, 43× the capture radius** (job 90).  A support that is ~0 because of where filaments were
    put is not a physiological zero, and the two are indistinguishable in the artifact.
  * The candidate pool is set by construction too: ``sf_population.py:555`` designates the entire
    ``interior`` of every fibre as ``cortex_sites``.
  * The card's other declared member answers no better, and neither does the INTERNAL edge:
    ``dorsal_arc_crosslink``'s joints are built at ``sf_population.py:573-581`` as one per dorsal free
    end paired to the nearest arc apex, **with no radius test at all** — its count is
    ``len(dorsal_free_nodes)``, a topology number.  It reached ``DEVICE_RUN`` without being asked.

**The basis is DERIVED and is not among the gaps.**  A transient actin crosslink is a two-headed
actin-actin linker, so its population scales with apposed filament pairs — ``per_filament``, resolved
against a MEASURED count of apposed SF filaments.  Not areal (no membrane is involved) and not
volumetric (what is counted is the BOUND linker, not the cytoplasmic pool).  What is absent is the
VALUE, its SCOPE and its AUTHORITY.  ``PARAM_PROVENANCE_AUDIT_2026-07-24.md`` rows 24-25 source
``link_k`` and ``k_off0`` and are silent on multiplicity; ``pair_sf_to_cortex_segments`` refuses
multiplicity > 1 with the words *"there is none"*.  Searched and absent, not unsearched.

The support is blocked one level down as well: ``world/build/stress_fiber.py`` records that **how many
stress fibres a cell has is absent from the KB entirely** (*"the spec describes ONE fibre"*), so the
apposed-filament support cannot be measured against a population that is itself a PI-GAP.

────────────────────────────────────────────────────────────────────────────────────────────────
3. WHAT GOES TO THE PI, AND WHAT THIS MODULE MAY NOT DO.

  a. **This is a candidate for the first family with kinetics.**  ``bond.py`` defers attach/detach, free
     lists and snapshot twins to *"the first family that has kinetics, which is also when the right
     shape for them becomes measurable rather than guessable"*.  The rates ARE sourced (α-actinin k_on
     10 /s, k_off0 0.066 /s, Bell f0; Ferrer 2008), so this family qualifies on the datum.  **It is not
     built here.**  Choosing the kinetic state shape is a PI decision, and inventing one inside a
     connector module bakes it in before it is measured.
  b. **``dorsal_arc_crosslink``: an undeclared card AND a contract/runtime contradiction.**  Its
     contract declares ``kinetics=True`` (``contracts.py:906``, field 5 of ``ConnectorContract``) while
     its runtime ``propose_events`` is a documented no-op — ``sf_mechanics.py:511``, *"declares no
     kinetics, so it advances no irreversible state"*.  A transient crosslink that binds and never
     unbinds is a permanent weld, which ``CLAUDE.md`` forbids by name and which
     ``sf_cortex_transient.py`` itself names when it refuses to bind without a rate law.  One of the two
     is wrong and this session may not pick which.
  c. **What crosslinker does ``transient_actin_crosslink`` denote?**  The card names none, while the
     builder spends α-actinin's sourced constants.  Folding ``lamellipodium_cortex_seam`` in therefore
     hands the lamellipodial rear seam α-actinin's stiffness and rates — a physiological claim nobody
     has made.  The fold is right on the criterion; the card's CONTENT is a PI question.
  d. **The multiplicity itself** — value, scope and authority, per §2.

⚠ The capture radius is never widened to produce more pairs.  Widening a sourced constant until the
geometry cooperates is tuning to an outcome, which ``CLAUDE.md`` forbids; if the count is too low that
is a PLACEMENT question for whoever owns the SF inventory.

⚠ This module writes no entry into :data:`aleph.world.families.DUPLICATE_OF`.  That dict lives in
``families/__init__.py``, which fourteen concurrent sessions share and this session does not own; the
fold is exposed as :data:`FOLDS_IN` for Lead to merge.

Sanity Gate (recorded per CLAUDE.md, before first execution):
  * dimensional — the count basis is ``per_filament`` [dimensionless] against a support in filaments;
    the constants this module CITES are pN/µm (stiffness) and µm (capture radius), and none of them is
    re-declared here as a runtime value — citing is not re-spelling.
  * boundary — the only path through this module raises. There is no zero-count branch to define,
    because a zero resolved from an absent density is indistinguishable from a physiological zero, and
    ``bond.py`` returns the latter as a legitimate answer.
  * conservation/invariant — nothing is claimed from the arena, so no BOND capacity is consumed and
    ``assert_partitioned`` is unaffected by importing or calling this module. Asserted in ``_demo``.
  * CFL/precision — no integration, no device, no float64 array.
  * sign sense — not applicable; no length, stiffness or force is produced.
  * measurement protocol — host-side declaration only. Nothing is uploaded, nothing is read back, and
    the module must never import :mod:`aleph.engine`, which is frozen PORT SOURCE
    (``test_layer_directions``). Every engine fact above is a CITATION with a file and a line.

engine units: length µm, stiffness pN/µm.  Runtime: pure host declaration, CPU-importable.
"""

from __future__ import annotations

from typing import NoReturn

from aleph.world.arena import WorldArena
from aleph.world.families import ConnectorGapError, FamilySpec

__all__ = [
    "CONNECTOR", "CHEMISTRY_CARD", "FOLDS_IN", "CARD_UNDECLARED", "SPEC",
    "build_sf_cortex_transient_family",
]

#: The declared connector name, exactly as the committed device-run census spells it.
CONNECTOR: str = "sf_cortex_transient"

#: The identity that distinguishes families, read from ``engine/contracts.py:734`` — not chosen here.
CHEMISTRY_CARD: str = "transient_actin_crosslink"

#: Connectors that carry this card and therefore fold INTO this family rather than becoming their own.
#: Confirmed against the contract table (§1); claimed independently by D13 from its own side. Exposed
#: for Lead to merge into ``families.DUPLICATE_OF``, which this session does not own.
FOLDS_IN: tuple[str, ...] = ("lamellipodium_cortex_seam",)

#: Connectors whose ``chemistry_card`` is ``None`` and which are physics candidates for this card. NOT
#: folded — an undeclared identity is an open question, and asserting a card for someone else's
#: connector is the invention this package exists to prevent. See §1 and §3b.
CARD_UNDECLARED: tuple[str, ...] = ("dorsal_arc_crosslink",)

SPEC = FamilySpec(
    connector=CONNECTOR,
    populations=("cortex", "stress_fiber"),
    basis="per_filament",
    blocked_by=(
        "count-value: no transient-actin-crosslink multiplicity exists for ANY population pair. "
        "PARAM_PROVENANCE_AUDIT rows 24-25 source link_k and k_off0 and are silent on how many; "
        "pair_sf_to_cortex_segments refuses multiplicity > 1 with 'there is none'. The engine's "
        "operative count is placement-derived - whatever falls inside the sourced 0.06 um radius - and "
        "the native measurement is median 2.59 um = 43x that radius, so it is ~0 by placement, not by "
        "physiology, and the artifact cannot tell the two apart.",

        "count-scope: no cell line, state, assay or temperature can be attached to a count that does "
        "not exist. Ferrer 2008 is a single-molecule AFM measurement of the bond and carries no "
        "per-cell inventory; the run target is MCF7.",

        "count-authority: no PI ruling. STATE.md (e) carries six open items and a transient actin "
        "crosslink density is not among them. This module is the surface that puts it there.",

        "support: the per_filament support cannot be MEASURED, because the SF population it would be "
        "measured against is itself a PI-GAP - world/build/stress_fiber.py records that how many "
        "stress fibres a cell has is absent from the KB entirely ('the spec describes ONE fibre').",

        "kinetics: this connector declares kinetics=True and bond.py deliberately has no kinetic state, "
        "deferring it to 'the first family that has kinetics'. The rates are sourced, so this family "
        "qualifies on the datum - which makes choosing the state shape a PI decision, not a connector "
        "module's.",

        "partner-population: the contract's SF endpoint is `sf_arc` (transverse/dorsal ARCS) while "
        "PHASE 1 built `stress_fiber` (ventral, FA-to-FA) - different structures, not a rename. SPEC "
        "names `stress_fiber` because it is the only SF-like population the arena has, and that is "
        "recorded here rather than chosen silently: either the arcs are an unbuilt population and this "
        "count is asked against the wrong one, or the contract endpoint is stale. PI call; D8 and D9 "
        "reached it from their own side (test_families_census.py, xfail-strict).",

        "identity-open: two questions about what this card DENOTES. (i) dorsal_arc_crosslink is "
        "ConnectorFamily.TRANSIENT_ACTIN with chemistry_card=None (contracts.py:905-909) and is already "
        "DEVICE_RUN; if that card is this one, its population belongs to this family and the count must "
        "cover it. Its contract also declares kinetics=True while sf_mechanics.py:511 no-ops "
        "propose_events - a permanent weld. (ii) the card names no crosslinker while the builder spends "
        "alpha-actinin's constants, so folding lamellipodium_cortex_seam in hands the lamellipodial "
        "seam alpha-actinin's stiffness and rates, which nobody has claimed physiologically.",
    ),
    source=(
        "SOURCED and sufficient for the FORCE: alpha-actinin link_k 4.6e5 pN/um and capture radius "
        "0.06 um (Ferrer 2008 PNAS AFM, PI-approved 2026-06-30); rates k_on 10/s, k_off0 0.066/s, Bell "
        "f0. SOURCED and ABSENT for the COUNT: nothing, in any of the four bases. That distinction is "
        "the finding - a family whose stiffness is impeccable and whose population is a placement "
        "artifact is the ERM defect with a better disguise."
    ),
)


def build_sf_cortex_transient_family(
    arena: WorldArena, *, support_filaments: float | None = None
) -> NoReturn:
    """Refuse to build, naming every unanswered question.

    The signature is the one this family WOULD take — an arena to claim a BOND range out of, and the
    MEASURED count of apposed SF filaments to resolve a ``per_filament`` density against. Both
    arguments are accepted and neither is used, so a caller who already has an arena and a support
    still cannot get a family: what is missing is the density, not the plumbing.

    Args:
        arena: the world the BOND range would be claimed from. Untouched — no capacity is consumed.
        support_filaments: the measured apposed-SF-filament count the density would resolve against.
            Cannot be supplied honestly today; see ``SPEC.blocked_by[3]``.

    Raises:
        ConnectorGapError: always, with seven unanswered questions.
    """
    raise ConnectorGapError(
        CONNECTOR,
        list(SPEC.blocked_by),
        note=(
            f"card {CHEMISTRY_CARD!r} — canonical member; {', '.join(FOLDS_IN)} folds in, so this "
            f"refusal blocks that edge too. Card undeclared and unresolved for: "
            f"{', '.join(CARD_UNDECLARED)}. Do not widen the 0.06 um capture radius to raise the count."
        ),
    )


def _demo() -> None:
    """Self-check: the module refuses, the refusal carries every question, and it claims no capacity."""
    assert SPEC.connector == CONNECTOR
    assert not SPEC.buildable, "a spec with no sourced count must not report itself buildable"
    assert len(SPEC.blocked_by) == 7, SPEC.blocked_by
    assert SPEC.basis == "per_filament", "a two-headed actin-actin linker scales with apposed pairs"
    # Arena population names. NOT a rename of the contract's `sf_arc`: see blocked_by 'partner-
    # population'. This assertion pins the side this module took so a PI ruling breaks it loudly.
    assert SPEC.populations == ("cortex", "stress_fiber"), SPEC.populations

    # A fold points at something OTHER than this connector, and the two lists must not overlap:
    # a card cannot be both matched and undeclared for the same edge.
    assert CONNECTOR not in FOLDS_IN and CONNECTOR not in CARD_UNDECLARED
    assert not set(FOLDS_IN) & set(CARD_UNDECLARED)

    # The builder refuses whatever it is handed, including a real arena and a plausible support.
    from aleph.world.arena import Kind

    arena = WorldArena(capacity={Kind.NODE: 8, Kind.BOND: 8})
    for kwargs in ({}, {"support_filaments": 1.0}, {"support_filaments": 0.0}):
        try:
            build_sf_cortex_transient_family(arena, **kwargs)
        except ConnectorGapError as exc:
            assert exc.connector == CONNECTOR
            assert len(exc.missing) == 7
            assert "count-value" in str(exc) and "partner-population" in str(exc)
            assert CHEMISTRY_CARD in str(exc) and FOLDS_IN[0] in str(exc)
        else:  # pragma: no cover
            raise AssertionError(f"build must refuse; got a family for {kwargs}")
    assert arena.n_live(Kind.BOND) == 0, "a refusal must not consume arena capacity"

    # Every connector this module names must be one the committed census actually spells — the cheapest
    # guard against a fold or a gap that points at a connector nobody declared.
    import json
    import pathlib

    census = pathlib.Path("aleph/outputs/ac/connector_devicerun/native_record.json")
    if census.exists():
        names = {c["connector"] for c in json.loads(census.read_text())["connectors"]}
        for n in (CONNECTOR, *FOLDS_IN, *CARD_UNDECLARED):
            assert n in names, f"{n} is not a declared connector"

    print(f"{CONNECTOR} self-check OK — {SPEC!r}, card={CHEMISTRY_CARD!r}, folds_in={FOLDS_IN}")


if __name__ == "__main__":
    _demo()
