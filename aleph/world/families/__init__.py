"""Bond families, one module per declared connector.

WHY ONE FILE EACH.  Twenty-eight of the thirty cell-internal connectors have no runtime object at all
— ``NOT_BOUND`` in the committed device-run census, which is weaker than a stub — so this is mostly
first implementation rather than porting, and it is being done by parallel sessions. The ownership lane
cannot separate them: a lane is per-lane, not per session instance, so twenty-eight sessions inside
``world`` would have no protection from each other at all. **File separation is the protection.** One
session owns one module here and reads everything else.

WHAT A MODULE IN HERE OWES.  A :class:`~aleph.world.bond.BondFamily` cannot be constructed without a
:class:`~aleph.world.bond.BondCount`, and that refuses without a scope, a provenance and a measured
support — *how many, against what, for which cell, on whose authority*. A module that cannot answer
those does not get to guess: it raises :class:`ConnectorGapError` naming what is missing, and that
refusal is the module's product rather than its failure. Lead collects them into the decision queue.

THE DEFECT THIS SHAPE EXISTS FOR.  The incumbent declared the membrane-cortex ERM linkage as ONE
connector between two components. Physically it is a POPULATION of tethers whose count is an areal
density, and because the declaration hid the population **nobody had to answer how many**. The answer
turned out to be the icosphere vertex count — 642 at subdivision 3, which is what every resting run of
2026-08-20 actually measured. ``bond.py`` names that class: *"a mesh number wearing a physiological
label"*.

⚠ AND SOME OF THE THIRTY ARE THE SAME PHYSICS. Four NMII motor edges differ only in which population
they land on; six ``*_cytosol_transfer`` edges likewise. A module whose connector is a duplicate says so
and does NOT build a second family — the component pair is derived from the arena's ID ranges at query
time, which is the whole reason ``bond.py`` stores two global node indices and nothing else. Judge by
physics, never by name.
"""

from __future__ import annotations

__all__ = ["ConnectorGapError", "FamilySpec", "DUPLICATE_OF", "DUPLICATE_CRITERION"]


class ConnectorGapError(NotImplementedError):
    """A family cannot be built because an answer is missing, and inventing one is forbidden.

    Raised with the questions that are unanswered, not with a placeholder value. The charter's rule is
    that a physiological value which is unknown is surfaced to the PI; this is that surface, made into a
    type so a partially-built family cannot be mistaken for a built one.

    Args:
        connector: the declared connector name this module owns.
        missing: what could not be answered, one string per question.
        note: where the gap is recorded, or what would close it.
    """

    def __init__(self, connector: str, missing: list[str], note: str = "") -> None:
        self.connector = connector
        self.missing = list(missing)
        self.note = note
        detail = "; ".join(self.missing)
        super().__init__(
            f"{connector}: cannot build a BondFamily — {len(self.missing)} unanswered: {detail}"
            + (f". {note}" if note else "")
        )


#: Connectors that are the SAME PHYSICS as another and must not build a second family. A module whose
#: name appears as a key does not build; it records the duplication and defers to the value. Empty at
#: creation on purpose — every entry must be argued from physics in the owning module's docstring, and
#: an entry added without that argument is exactly the name-based folding this package forbids.
DUPLICATE_OF: dict[str, str] = {
    # Merged by Lead from each module's own DUPLICATE_OF_ENTRY. A module does not write this dict
    # itself — fourteen sessions share the file — so it exposes its verdict and Lead merges it here.
    # The argument for each lives in the owning module's docstring, which is where it must stay.
    "lamellipodium_cortex_seam": "sf_cortex_transient",       # card transient_actin_crosslink (D13)
    "filopodium_membrane_tip": "lamellipodium_membrane_contact",  # card actin_membrane_brownian_ratchet_contact (D12/D14)
    "sf_cortex_transient": "dorsal_arc_crosslink",            # (D10)
}

#: THE CRITERION, converged 2026-08-20 by three sessions independently and recorded here so the
#: remaining eleven do not re-derive it.
#:
#:     A family's identity is its CHEMISTRY CARD. A shared runtime class is not a duplicate, and a
#:     shared force law is not a duplicate.
#:
#: Each of the three arrived at it from a different direction, which is why it is written down:
#:
#: * `if_nucleus_linc` (D4) — the three LINC edges share ONE runtime class and one joint SoA, which
#:   looks like a duplicate. But `_EDGE_CHEMISTRY` gives them three different cards
#:   (`nesprin_actin_linc` / `microtubule_motor_linc` / `nesprin3_plectin_if_linc`). Verdict: the
#:   shared runtime is a shared CONTAINER, not shared physics. Folding them would hide three
#:   chemistries under one name — the mirror of the ERM defect, where one name hid a population.
#: * `lamellipodium_cortex_seam` (D13) — genuinely IS a duplicate of `sf_cortex_transient`: both
#:   declare `chemistry_card="transient_actin_crosslink"` and `filament_crosslink.py` says five edges
#:   are "the SAME mechanics". Its control is what makes the call trustworthy: the sibling
#:   `filopodium_cortex_root` sits in the same `ConnectorFamily.TRANSIENT_ACTIN` but carries a
#:   different card, `formin_fascin_root_coupling` — so folding by card is not folding by name.
#: * `if_sf_plectin` (D9) — `plectin_if_actin` is unique among the five edges; only the force law is
#:   shared. Not a duplicate.
#:
#: ⚠ AND FOLDING DOES NOT DISSOLVE THE COUNT. D13's second finding: `sf_cortex_transient.py` records
#: of itself that "THE COUNT IS PLACEMENT-DERIVED, NOT A SOURCED INVENTORY" — the number is however
#: many partners happened to fall inside a capture radius. That is the ERM defect in another dress: a
#: PLACEMENT number wearing a physiological label instead of a mesh number wearing one. A duplicate
#: verdict answers "is this one family or two"; it never answers "how many".
DUPLICATE_CRITERION = "chemistry card; a shared runtime class or force law is not a duplicate"


class FamilySpec:
    """What one connector module declares about itself, whether or not it can build.

    Every module in this package exposes ``SPEC`` of this type, so Lead can enumerate the whole census —
    built and blocked alike — without importing anything that raises.

    Args:
        connector: the declared connector name, exactly as it appears in the device-run census.
        populations: the two population names the bond joins, as the arena knows them. Recorded for the
            census only; the runtime pair is DERIVED from ID ranges at query time and is never stored on
            the bond.
        basis: what the count is measured against — ``"areal"``, ``"volumetric"``, ``"per_filament"``
            or ``"explicit"``.
        blocked_by: the unanswered questions. Empty means the module can build.
        source: provenance for the count, or the reason there is none.
    """

    __slots__ = ("connector", "populations", "basis", "blocked_by", "source")

    def __init__(self, connector: str, populations: tuple[str, str], basis: str,
                 blocked_by: tuple[str, ...] = (), source: str = "") -> None:
        if not connector:
            raise ValueError("a family spec needs the connector name it answers for")
        if len(populations) != 2:
            raise ValueError(f"{connector}: a bond joins exactly two populations, got {populations!r}")
        if basis not in ("areal", "volumetric", "per_filament", "explicit"):
            raise ValueError(f"{connector}: basis {basis!r} is not a recognised count basis")
        if not blocked_by and not source.strip():
            raise ValueError(
                f"{connector}: a family that claims to be buildable must name its provenance — an "
                "unblocked spec with no source is the defect this package exists to prevent"
            )
        self.connector = connector
        self.populations = populations
        self.basis = basis
        self.blocked_by = tuple(blocked_by)
        self.source = source

    @property
    def buildable(self) -> bool:
        """True when nothing is unanswered. Says nothing about whether the answer is RIGHT."""
        return not self.blocked_by

    def __repr__(self) -> str:
        state = "buildable" if self.buildable else f"blocked({len(self.blocked_by)})"
        return f"FamilySpec({self.connector!r}, {self.populations!r}, {self.basis!r}, {state})"
