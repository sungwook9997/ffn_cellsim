r"""``mt_sf_spectraplakin`` — the MT↔actin cytolinker family, which REFUSES to be built.

WHAT THIS CONNECTOR IS.  Spectraplakins (MACF1/ACF7) are ~500 kDa cytolinkers with an N-terminal
actin-binding CH pair and a C-terminal GAS2 microtubule-binding domain, so one molecule is a mechanical
bond between a microtubule lattice site and an actin material point.  ``engine/contracts.py:770``
declares it ``ConnectorFamily.SPECTRAPLAKIN``, ``microtubule -> sf_arc``, ``kinetics=True``,
``commit_on_accept=True``, chemistry card ``spectraplakin_mt_actin``.

WHAT ITS RUNTIME STATE IS.  ``connector_devicerun/native_record.json`` (build ``73a93883``, FULL NATIVE,
RTX 4090-1) records it ``NOT_BOUND`` — *"no runtime object in this composed world"*, ``n_calls: 0``,
``n_launch_kernels: 0``.  ``sf_mechanics.py:121`` calls it ``SEAMED``.  So this is first
implementation, not a port, and the artifact it produces is a REFUSAL.

---

IT IS NOT A DUPLICATE OF ``if_sf_plectin`` (D9), AND THE ARGUMENT IS NOT THE NAME.

``engine/filament_crosslink.py`` folds five declared edges — including this one and ``if_sf_plectin`` —
onto ONE force law, and it is right to: *"The family label is biology; the force law is one."*  A
shared force law is a shared FORCE LAW.  It is not a shared family, and the distinction is the whole of
``bond.py``'s design:

  * ``bond.py`` states the identity outright — ``chemistry_card`` is *"the identity that actually
    distinguishes families.  Two declared edges sharing a card are ONE family; the declared component
    pair never was the identity."*  The cards differ: ``spectraplakin_mt_actin`` vs ``plectin_if_actin``,
    and ``filament_crosslink.py`` says so itself — *"What differs across the five is CHEMISTRY, not
    mechanics."*
  * The duplicate case this package was written for is the four NMII motor edges: **one molecule**, one
    chemistry, four declarations differing only in which population they land on — which the arena
    derives from ID ranges anyway.  Here the MOLECULE differs.  Spectraplakin's partner is a
    microtubule via GAS2; plectin's is an intermediate filament.  Different partner, different reach,
    different on/off rates, and a count that is per-MT rather than per-IF.
  * Folding them would be the ERM defect in mirror image.  There, one declaration hid a POPULATION and
    nobody had to answer how many.  Here, one family would hide two CHEMISTRIES and nobody would have to
    answer which count was which — the same failure, one level up.

So ``DUPLICATE_OF`` gets no entry from this module, and that too is a judgement rather than an omission.
The same argument decides ``if_sf_plectin`` (D9) the same way; it is written here once.

⚠ **A separate naming discrepancy, recorded rather than resolved.**  The contract's endpoint is
``sf_arc`` — transverse/dorsal arcs.  What PHASE 1 actually stood up is population ``stress_fiber``:
``build/stress_fiber.py`` builds VENTRAL fibres, FA-to-FA on the basal plane.  Arcs and ventral fibres
are different structures with different MT apposition, so the population this family would bond to is
itself unsettled.  It is folded into the support question below rather than answered here.

---

THE FOUR QUESTIONS, AND WHY THREE OF THEM CHAIN.

``BondCount`` asks *how many · against what · for which cell · on whose authority*.  Answering them in
that order is impossible here, because the first three are one chain:

  1. **reach** (``rest_um``) — the crosslink's zero-force separation.  ABSENT.  The nearest KB datum is
     the sibling cytolinker claim ``KB-DRAFT-3.B-10``, which offers ``crosslink_reach = ~l_seg um`` and
     labels it *"placeholder — PI-anchor needed"*.  A reach set to the segment length is a
     DISCRETISATION, which is precisely the class of error ``bond.py`` exists to refuse: *"a mesh number
     wearing a physiological label"*.
  2. **support** — what the density is measured against.  It is measurable IN PRINCIPLE and only once
     (1) is answered: both populations stand in the arena today, so a neighbour query at the reach
     yields the MT-to-SF apposition (length, or candidate-pair count) directly from the built geometry.
     Without a reach there is no query, so the support cannot be measured — and per ``bond.py`` it must
     be *measured, never analytic*.
  3. **count / basis** — ABSENT, and the KB says so in the claim that covers exactly this connector.
     ``KB-DRAFT-3-11`` *"MT↔actin / MT↔cortex distributed crosslinks (MAPs, +TIPs, spectraplakins)"*
     carries ``MT-actin crosslink stiffness = NOT-IN-KB pN/µm`` and no count in any basis.  Whether the
     right basis is per-microtubule, per µm of apposition, or volumetric is not a modelling preference:
     it is what the missing measurement would have been reported in.
  4. **authority** — the claim is ``status: draft``, ``assumptions: AUDIT-2026-07-15 DRAFT (NOT-IN-TAG).
     PROVISIONAL id — PI assigns final KB number at ingest``, and its ACF7 source (Kodama et al. 2003
     Cell 115:343) is annotated **ABSENT** in the claim's own citation string.  ``source_audit`` holds no
     row for it, so under CLAUDE.md — *"before citing a KB source in a deliverable, confirm its
     ``source_audit.verdict`` is OK"* — it may not be cited as a source at all.  And Kodama is mouse
     ACF7; the target scope is MCF7 epithelial, so even a recovered value would enter as
     ``UNRATIFIED_PROXY``, not ``SOURCED``.

**One PI ruling — the reach — converts (2) from unanswerable into a measurement.**  That is why the
questions are listed in this order and not by severity.

⚠ **AND THIS FAMILY IS KINETIC, WHICH IS A SECOND ESCALATION.**  The contract declares
``kinetics=True``/``commit_on_accept=True``, and ``filament_crosslink.py`` is explicit about what a
kinetic edge without a rate law becomes: *"A transient crosslink with no rate law is not a transient
crosslink; it is a permanent weld, which is the one thing the charter says a connector may never be."*
``bond.py`` deliberately holds no kinetics — *"They arrive with the first family that has kinetics"* —
and no spectraplakin on/off-rate card exists.  So even a fully answered ``BondCount`` would not license
building this family here.  That is a PI decision and this module does not take it.

Sanity Gate (recorded per CLAUDE.md, before first execution):
  * dimensional — nothing is computed.  Had it built: ``rest_um`` [µm], ``stiffness_pn_per_um`` [pN/µm],
    and the count in whichever basis the missing measurement was reported in.  The basis recorded on
    :data:`SPEC` is the CANDIDATE, not a ruling — ``FamilySpec`` requires the field so the census can be
    enumerated without importing anything that raises, and ``blocked_by`` names it as unanswered.
  * boundary — there is no boundary case, because there is no admissible input: the builder raises on
    every call, including one that supplies a count, since a supplied count would be an invented one.
  * conservation/invariant — none reached.  The force law it would bind (equal-and-opposite central,
    ``filament_crosslink.py``'s Sanity Gate) is the incumbent's and is not re-derived here.
  * CFL/precision — no integration, no device, float64 throughout the (absent) arithmetic.
  * sign sense — not reached.
  * measurement protocol — pure host bookkeeping; no device is touched, nothing is uploaded, and no
    number is produced.  Nothing here is a CPU result because nothing here is a result.

engine units: length µm, stiffness pN/µm.  Runtime: pure host bookkeeping, CPU-importable.
"""

from __future__ import annotations

from typing import Never

from aleph.world.families import ConnectorGapError, FamilySpec

__all__ = ["CONNECTOR", "SPEC", "build_mt_sf_spectraplakin"]

CONNECTOR = "mt_sf_spectraplakin"

#: The unanswered questions, in the order the chain resolves: a reach makes the support measurable, and
#: only a measured support can carry a density. Names, not values — a value here would be the invention
#: this package forbids.
BLOCKED_BY: tuple[str, ...] = (
    "rest_um: the spectraplakin crosslink reach is absent; the only KB figure is the sibling "
    "cytolinker claim KB-DRAFT-3.B-10's 'crosslink_reach = ~l_seg um', self-labelled 'placeholder — "
    "PI-anchor needed', and a reach equal to the segment length is a discretisation, not a length",
    "support: the MT-to-SF apposition is measurable from the standing PHASE 1 arena by a neighbour "
    "query — but only at a reach, so it is blocked ON rest_um rather than absent from the world",
    "count + basis: no spectraplakin count exists in any basis. KB-DRAFT-3-11 covers exactly this "
    "connector and carries no count; whether it is per-microtubule, per um of apposition or volumetric "
    "is what the missing measurement would have been reported in",
    "stiffness_pn_per_um: KB-DRAFT-3-11 records 'MT-actin crosslink stiffness = NOT-IN-KB pN/um'",
    "scope + authority: KB-DRAFT-3-11 is status=draft, NOT-IN-TAG, PROVISIONAL id, and annotates its "
    "own ACF7 source (Kodama 2003 Cell 115:343) as ABSENT with no source_audit row. Kodama is mouse "
    "ACF7 against a target scope of MCF7 epithelial, so a recovered value enters as UNRATIFIED_PROXY",
    "kinetics: the contract declares kinetics=True/commit_on_accept=True and no spectraplakin on/off "
    "rate card exists. bond.py holds no kinetics by design; building one here is a PI decision",
    "endpoint population: the contract's endpoint is 'sf_arc' (transverse/dorsal arcs) but PHASE 1 "
    "stood up 'stress_fiber' (ventral, FA-to-FA, basal plane). Different structures, different MT "
    "apposition — which population this family bonds to is itself unsettled",
)

#: What this module declares about itself, buildable or not. ``populations`` uses the arena's names, per
#: :class:`~aleph.world.families.FamilySpec`; the contract's ``sf_arc`` is the discrepancy above.
SPEC = FamilySpec(
    connector=CONNECTOR,
    populations=("microtubule", "stress_fiber"),
    basis="per_filament",  # CANDIDATE — per microtubule. Named in BLOCKED_BY as unanswered.
    blocked_by=BLOCKED_BY,
    source="",  # there is none; that is the finding
)


def build_mt_sf_spectraplakin(*_args: object, **_kwargs: object) -> Never:
    """Refuse to build the family, naming what is unanswered.

    The signature takes and ignores everything on purpose: there is no argument list that would make
    this buildable, because supplying a count IS the invention. It raises before touching an arena.

    Raises:
        ConnectorGapError: always, carrying :data:`BLOCKED_BY`.
    """
    raise ConnectorGapError(
        CONNECTOR,
        list(BLOCKED_BY),
        note="Not a duplicate of if_sf_plectin — same force law, different chemistry card; see the "
             "module docstring. Answering rest_um first makes the support measurable from the "
             "standing PHASE 1 arena, which is the cheapest ruling that moves this.",
    )


def _demo() -> None:
    """Self-check: the module refuses, the refusal names its questions, and it claims no duplicate."""
    assert not SPEC.buildable, "a spec with unanswered questions must not report buildable"
    assert SPEC.connector == CONNECTOR
    assert SPEC.populations == ("microtubule", "stress_fiber")
    assert len(SPEC.blocked_by) == len(BLOCKED_BY) == 7

    # Every entry is a NAMED question, not a value: the name comes first, before a colon.
    for q in SPEC.blocked_by:
        assert ":" in q and q.split(":", 1)[0].strip(), q

    # The build refuses whatever it is handed — including a caller who brought a count.
    for args, kwargs in ((), {}), ((object(),), {"count": 1000, "support": 12.5}):
        try:
            build_mt_sf_spectraplakin(*args, **kwargs)
        except ConnectorGapError as exc:
            assert exc.connector == CONNECTOR
            assert len(exc.missing) == len(BLOCKED_BY)
            assert "rest_um" in str(exc) and "if_sf_plectin" in str(exc)
        else:  # pragma: no cover
            raise AssertionError("this family must not build")

    # The duplicate judgement is NEGATIVE, and a negative judgement leaves the table alone.
    from aleph.world.families import DUPLICATE_OF

    assert CONNECTOR not in DUPLICATE_OF, "this connector is not a duplicate; see the docstring"

    print(f"{CONNECTOR} self-check OK — blocked on {len(SPEC.blocked_by)}, no family built")


if __name__ == "__main__":
    _demo()
