"""``lamellipodium_membrane_contact`` — a CONSTRAINT for certain, and a bond only under one reading.

⚠ **RULED 2026-08-21, PI decision 1**: this edge is a non-penetration CONSTRAINT and ``world/`` gets a
destination for it — :mod:`aleph.world.contact`. This module declares :data:`CONTACT` on the shared
:data:`~aleph.world.contact.GAP_LAW_POINT_SURFACE`.

⚠⚠ **AND HERE THE RULING AND THE FROZEN SOURCE DO NOT LINE UP — SURFACED, NOT PAPERED OVER.** The
decision names this edge as one of *"the three cardless CONTACT edges"*. **It is not cardless.**
``contracts.py:791`` gives it ``chemistry_card="actin_membrane_brownian_ratchet_contact"`` and
``kinetics=True``; ``connector_joints.py:7-13`` splits the five CONTACT edges 3/2 on exactly that, and
the third CARDLESS edge is ``membrane_ecm_contact``, not this one.

**The ruling nevertheless holds for this edge, and for a better reason than a correction.** All five
CONTACT edges compute the SAME unilateral gap — the source says so in its own words — and under
:data:`aleph.world.contact.CONTACT_IDENTITY` that shared Phi is what makes them one constraint family.
A growing barbed end may push a membrane and may **never tow it inward**, under EVERY reading of the
ratchet. So the constraint half of this edge is unconditional and lands here today.

**What the card buys is a SECOND object, not a different one.** ``families.DUPLICATE_CRITERION`` says
a shared force law is not a duplicate, and that is precisely this case: the constraint is shared, the
chemistry is not. If §1 below resolves to the *tethered* ratchet, this edge is a constraint **plus** a
bond population; if it resolves to the *elastic* ratchet, it is a constraint and the card is vestigial.
Either way the constraint is already true, so :data:`CONTACT` is declared and the bond refusal below
stands verbatim, unweakened by the ruling. **One edge, two shapes** — which is why the decision's
three-way grouping and the source's 3/2 split can both be right about it.

The bond half of this module builds nothing. Under :mod:`aleph.world.families` that is a product, not a failure: the
four questions :class:`~aleph.world.bond.BondCount` asks — *how many, against what, for which cell, on
whose authority* — have five unanswered pieces between them, and every candidate answer available
today is either a geometry choice nobody has made, a density of the wrong THING, a density in a basis
that does not exist, or a keratocyte measurement wearing this cell's name.

WHAT THE FROZEN PORT SOURCE DECLARES (read, never imported — ``test_layer_directions`` forbids the
import). ``engine/contracts.py:791`` declares ``ConnectorFamily.CONTACT``, ``lamellipodium`` <->
``membrane``, ``kinetics=True``, ``commit_on_accept=True``, endpoint roles *"growing barbed-end
material point"* and *"membrane contact quadrature"*, chemistry card
``actin_membrane_brownian_ratchet_contact``. ``engine/connector_joints.py`` resolves it to
``JointKind.BROWNIAN_RATCHET_CONTACT`` and gives it a unilateral gap law that may push and may never
pull. ``connector_devicerun/native_record.json`` records the composed native run's verdict:
``NOT_BOUND``, ``n_calls: 0``, *"no runtime object in this composed world"*. So this is first
implementation, not a port.

**1. THE PRIOR QUESTION — IS THIS A BOND?**  ``bond.py`` draws the line: *"Field interaction (steric,
drag) needs no declaration because it depends on position alone; a specific persistent pairing does."*
Three of the five declared CONTACT edges carry ``kinetics=False`` and are therefore field interaction
by that criterion — they need no family. This edge and ``filopodium_membrane_tip`` carry
``kinetics=True``, and the kinetics is the ONLY thing that could make either a bond. Which it is
depends on a model choice nobody has recorded:

* the **elastic** Brownian ratchet — a filament tip bends and pushes when a thermal fluctuation opens a
  gap, and is never attached. Position alone. **Not a bond**, and this whole module is then a
  mis-declaration to be closed rather than a family to be built.
* the **tethered** ratchet — a subpopulation is bound to the membrane through attachment proteins while
  the rest push. **That subpopulation is a bond population**, and it is the only reading under which
  ``BondCount`` even applies.

The paper that separates them is *Mogilner & Oster 2003, "Force Generation by Actin Polymerization II:
The Elastic Ratchet and Tethered Filaments"*. It is registered in this repository's KB —
``source_audit.citation_key = MogilnerOster2003_BJ``, ``verdict = OK`` — and it has **zero ingested
chunks**. Registered, audited, and unreadable: the same defect class
``PER_CELL_STRUCTURE_COUNTS_2026-08-20.md`` §4 records for ``KB-3.8``'s ``r_filo``, where a
``verified/High`` row could not be checked against its own citations. **The corpus cannot decide
whether this connector is a bond**, so this module does not decide it either.

**2. HOW MANY — the one verified density is not this population's count, for three separate reasons.**
``KB-3.7`` (``verified``, ``confidence: High``) reads *"fil density ~100 per µm leading edge"*.

* *It counts the wrong objects.* 100/µm is a **filament** density; a bond here is a **contacting
  barbed end**. Mueller et al. 2017 — the source under that row — separates the two explicitly, reports
  filament density and barbed-end density as distinct quantities, states that a barbed end *"could
  either represent an actively growing or a capped filament"*, and confines elongation to a strip of
  width ``w`` behind the edge. The bond population is the growing, membrane-adjacent subset of the
  barbed ends, and **no fraction for it is stated anywhere in the corpus.**
* *It is in a basis that does not exist.* 100/µm is a LINEAR density, per µm of leading edge.
  ``BondCount.basis`` offers ``areal`` / ``volumetric`` / ``per_filament`` / ``explicit`` and nothing
  linear. Converting it to ``areal`` needs the ~200 nm sheet thickness as a second datum AND a MEASURED
  leading-edge membrane strip to resolve against — and the arena's membrane is an icosphere with no
  leading edge on it, so the only support in reach would once again be a mesh feature. That is
  precisely the ERM shape this package exists to refuse, and dressing it in an areal basis would hide
  it better rather than fix it.
* *Its scale factor is unratified.* The count follows from a leading-edge width, which is **PI decision
  1**, still open: 5–20 µm gives 500–2,000 filaments against ``KB-3.7``, while the standing build
  carries an explicit 200 — i.e. **2.5–10x under**, which is a fact about a placeholder and not a
  measurement.

**3. FOR WHICH CELL — the sourcing is fine and the SCOPE is missing.**  ``KB-3.7`` carries no cell line
in any field. Its five citations all audit ``OK`` (``Krause2014_NRMCB``, ``Mueller2017_Cell``,
``Kage2017_NatCommun``, ``Bieling2016_Cell``, ``MalyBorisy2001_PNAS``), and the density traces to
Mueller 2017's electron tomography of **fish keratocyte** lamellipodia. Under the PI's 2026-08-15
ruling that is :attr:`~aleph.world.bond.SourceClass.UNRATIFIED_PROXY`, not ``SOURCED`` — the Π₀ = 40 Pa
case verbatim, where a real measurement of the wrong cell passed as a physiological value because
nobody recorded the scope.

**4. DUPLICATE — this edge and ``filopodium_membrane_tip`` are ONE family, argued from physics.**
``bond.py``: *"chemistry_card: the identity that actually distinguishes families. Two declared edges
sharing a card are ONE family; the declared component pair never was the identity."* Both edges declare
``actin_membrane_brownian_ratchet_contact``; both resolve to ``JointKind.BROWNIAN_RATCHET_CONTACT``;
the port source's own test asserts that those two and no others are the ratchet kind, read from
``kinetics`` rather than from the card. The physics is one statement — *a growing actin barbed end
pushes on a membrane it cannot pull* — and it does not consult whether the filament arrived in a
branched sheet or a bundle. What genuinely differs (branched vs bundled architecture, the curvature and
tension the membrane opposes with) belongs to the two POPULATIONS and to the membrane's own law, not to
the pairing. And since a bond stores two global node indices with the population derived from arena ID
ranges at query time, one family covers both without a second declaration.

That is :data:`~aleph.world.families.DUPLICATE_CRITERION` applied verbatim — *a family's identity is its
chemistry card; a shared runtime class or force law is not a duplicate* — and here the card is shared,
not merely the law. ⚠ **A duplicate verdict answers "one family or two"; it never answers "how many",**
so §2 survives the fold intact.

**Which of the two is canonical is Lead's call and Lead has made it (2026-08-20): this name.** The
argument is the card's own mechanism — ``KB-3.6`` is registered as *"Brownian ratchet polymerization
velocity (**lamellipodium** protrusion)"*, so the chemistry is named for this population and the
filopodial tip is the same ratchet pushing the same membrane, not a second chemistry. The direction
matters operationally rather than physically: :data:`aleph.world.families.DUPLICATE_OF` maps *pusher ->
canonical*, so exactly ONE side of a pair may write an entry — if both did, the pair would cycle and
nothing would be built. **This module is the canonical side and therefore writes no entry**; D12
(``filopodium_membrane_tip``) writes it. :data:`DUPLICATE_OF_PROPOSED` here is the census record of that
ruling, not a second entry.

⚠ And being canonical is not being buildable. It settles who is responsible for answering §1-§3; it
does not answer them, and this module still refuses.

**5. WHAT IS DELIBERATELY NOT BUILT — the kinetic law.**  ``bond.py`` defers kinetics to *"the first
family that has kinetics"*, and this edge would be it: ``kinetics=True`` and ``commit_on_accept=True``,
with ``KB-3.6`` (``verified``) already holding the rate law it would want — ``v(F) = v0(1 - F/F_stall)^w``,
``w ~ 2``, ``F_stall`` per filament ~2–5 pN, ``δ = 2.7 nm``. Writing it here would put attach/detach,
a free list and snapshot twins into the tree on the strength of a family that cannot state its own
population. **Escalated to the PI, not built.**

Sanity Gate (recorded per CLAUDE.md, before first execution — for a module that refuses, the gate is
about the refusal being total and checkable):
  * dimensional — the gate's own finding. ``KB-3.7``'s density is [1/µm]; ``_SUPPORT_UNIT`` admits
    [µm²], [µm³], [filaments] and none. A mismatched basis is a unit error by name in ``BondCount``,
    and this one is caught here instead, before a number reaches it.
  * boundary — there is no count to bound. A resolved zero would be a physical answer; a count that was
    never asked for is not, and that is the state this connector is in.
  * conservation/invariant — nothing is claimed from the arena, so no BOND range is consumed and
    ``assert_partitioned`` is unaffected by importing or calling this module. Asserted in ``_demo``.
  * CFL/precision — no integration, no stiffness recorded, nothing enters a CFL bound.
  * sign sense — recorded for whoever unblocks this: the law is UNILATERAL. It may push and may never
    pull, gated on the gap, structurally and not by a stiffness value. A contact that pulls is an
    adhesion nobody declared.
  * measurement protocol — host-side, no device touched, nothing uploaded, nothing built. The module
    is importable on a CUDA-less machine, which is what lets Lead enumerate :data:`SPEC` for the whole
    census without importing anything that raises.

engine units: length µm, areal density 1/µm². Runtime: pure host bookkeeping, CPU-importable.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, NoReturn

from aleph.world.contact import GAP_LAW_POINT_SURFACE, ContactScope
from aleph.world.families import ConnectorGapError, FamilySpec

if TYPE_CHECKING:  # the census must not drag the arena in just to read a spec
    from aleph.world.arena import WorldArena

__all__ = ["CONNECTOR", "SPEC", "CONTACT", "DUPLICATE_OF_PROPOSED",
           "build_lamellipodium_membrane_contact"]

CONNECTOR = "lamellipodium_membrane_contact"

#: The pair this connector is the CANONICAL half of, argued from physics in §4 and ruled by Lead
#: 2026-08-20. Read *pusher -> canonical*, the same direction as
#: :data:`aleph.world.families.DUPLICATE_OF` — but the ENTRY there is D12's to write, not this
#: module's: only one side of a pair may push, or the two would cycle and neither would build. This is
#: the census record of the ruling, and a module may not edit the shared package file in any case.
DUPLICATE_OF_PROPOSED: dict[str, str] = {"filopodium_membrane_tip": CONNECTOR}

#: The unanswered questions, one per entry, in the order they must be answered: whether there is a
#: family at all comes before how big it is.
BLOCKED_BY: tuple[str, ...] = (
    "is-this-a-bond: elastic ratchet (position-only, NOT a bond) vs tethered ratchet (a persistent "
    "pairing, a bond). MogilnerOster2003_BJ is registered with source_audit verdict OK and has ZERO "
    "ingested chunks, so the corpus cannot decide it. Answer this one first — the other four are "
    "moot if the answer is 'elastic'",
    "contacting-barbed-end-fraction: KB-3.7's 100/um counts FILAMENTS; a bond is a growing barbed end "
    "in contact. Mueller2017 separates filament density from barbed-end density and states no "
    "fraction, so the bond population is not the filament population",
    "count-basis: KB-3.7's density is LINEAR (1/um of leading edge) and BondCount has no linear "
    "basis. An areal conversion needs the ~200 nm sheet thickness AND a MEASURED leading-edge "
    "membrane strip, which an icosphere membrane does not have — the ERM defect in a new dress",
    "leading-edge-width-um: PI decision 1, open. 5-20 um gives 500-2,000 against KB-3.7; the standing "
    "build carries an explicit 200, i.e. 2.5-10x under",
    "scope: KB-3.7 names no cell line and its density is keratocyte (Mueller2017, all five citations "
    "audit OK). Under PI 2026-08-15 that is UNRATIFIED_PROXY for the target cell, not SOURCED",
)

SPEC = FamilySpec(
    connector=CONNECTOR,
    populations=("lamellipodium", "membrane"),
    # Declared for the census as the basis a density would have to reach before it could be sourced
    # (`build/lamellipodium.py` says the same of the population's own count). That it does not fit
    # KB-3.7's linear density is itself one of the blockers above, not a detail this hides.
    basis="areal",
    blocked_by=BLOCKED_BY,
)

#: **The constraint half, and it is NOT contingent on any of** :data:`BLOCKED_BY`. A barbed end may
#: push the membrane and may never tow it inward, under both readings of §1, so this scope is true
#: today whether or not a bond population is ever declared beside it.
#:
#: ⚠ Its blocker is a DIFFERENT question from the five above and shares none of them: not *how many*
#: (a constraint has no population — the active set is an output of the broad phase), not *what
#: fraction of barbed ends* and not *what leading-edge width*. It is the actin/membrane contact
#: distance, and the ruling did not supply one either.
CONTACT = ContactScope(
    connector=CONNECTOR,
    point_population="lamellipodium",  # growing barbed ends
    surface_population="membrane",     # owns FACEs (build/membrane.py)
    blocked_by=(
        "contact distance d0 [um] for an actin barbed end against the membrane, WITH its scope. "
        "ABSENT: connector_joints.py defines no magnitude by design. Note it is NOT the ratchet's "
        "delta = 2.7 nm (KB-3.6) — that is a monomer half-length, the step a polymerisation event "
        "advances by, and reading it as the separation at which a contact force vanishes would be the "
        "unit-shaped mistake this package exists to catch.",
    ),
)

#: Recorded for whoever unblocks this, so the next session does not re-derive it: the kinetic law this
#: family would need already exists as KB-3.6 (verified) — v(F) = v0 (1 - F/F_stall)^w, w ~ 2,
#: F_stall ~2-5 pN per filament, delta = 2.7 nm. It is NOT implemented here; see §5 of the docstring.
KINETIC_LAW_IF_UNBLOCKED = "KB-3.6 (verified) — PI escalation, not an implementation licence"


def build_lamellipodium_membrane_contact(arena: WorldArena) -> NoReturn:
    """Refuse to build, naming what is unanswered.

    The signature is the one a builder will have when the gaps close, so unblocking this module is an
    edit to its body rather than to its callers. Nothing is claimed from ``arena``; it is accepted and
    left untouched.

    Args:
        arena: the world a family would be claimed out of. Unused, deliberately.

    Raises:
        ConnectorGapError: always, with the five unanswered questions.
    """
    del arena  # nothing is claimed: a refusal must not consume a BOND range
    raise ConnectorGapError(
        CONNECTOR,
        list(BLOCKED_BY),
        note=(
            "The CONSTRAINT half of this edge is already declared as CONTACT (PI decision 1, "
            "2026-08-21) and does not wait on any of these: a barbed end may push the membrane and "
            "may never tow it, under both readings of the ratchet. What is refused here is the BOND "
            "half, which exists only under the TETHERED reading. "
            "Do not invent a count. Also: this edge and filopodium_membrane_tip share the chemistry "
            "card actin_membrane_brownian_ratchet_contact and are ONE family (see DUPLICATE_OF_PROPOSED) "
            "— building both would be the second family this package forbids. The kinetic law is "
            "escalated to the PI, not built: bond.py defers kinetics to the first family that has it, "
            "and this would be that family."
        ),
    )


def _demo() -> None:
    """Self-check: the refusal is total, it names its questions, and it consumes nothing."""
    from aleph.world.arena import Kind, WorldArena

    assert not SPEC.buildable, "a spec with unanswered questions may never read as buildable"
    assert len(SPEC.blocked_by) == 5
    assert SPEC.populations == ("lamellipodium", "membrane")

    # The bond/field question comes first: an 'elastic' answer closes the connector instead of sizing it.
    assert SPEC.blocked_by[0].startswith("is-this-a-bond"), "the prior question must lead the queue"

    arena = WorldArena(capacity={Kind.NODE: 8, Kind.SEGMENT: 8, Kind.ANGLE3: 8, Kind.BOND: 8})
    before = arena.n_live(Kind.BOND)
    try:
        build_lamellipodium_membrane_contact(arena)
    except ConnectorGapError as exc:
        assert exc.connector == CONNECTOR
        assert len(exc.missing) == 5
        assert "leading-edge-width-um" in "; ".join(exc.missing)
        assert "filopodium_membrane_tip" in exc.note, "the duplication must travel with the refusal"
    else:  # pragma: no cover
        raise AssertionError("this connector cannot be built and must say so by type")

    # A refusal that consumed a claim would leave a family half-built and countable.
    assert arena.n_live(Kind.BOND) == before == 0
    arena.assert_partitioned()

    # The fold is one-directional and does not name itself, or the census would chase its own tail.
    assert DUPLICATE_OF_PROPOSED == {"filopodium_membrane_tip": CONNECTOR}
    assert CONNECTOR not in DUPLICATE_OF_PROPOSED

    # THE CONSTRAINT HALF is declared and does NOT inherit the bond's blockers. If it ever does, the
    # two shapes have been collapsed into one and the ruling has been undone by an edit.
    from aleph.world.contact import ContactFamily

    # ROLES, not order — and this scope is why roles were needed: it was declared point-first
    # while the other two were declared surface-first, and nothing could notice.
    assert CONTACT.point_population == "lamellipodium"
    assert CONTACT.surface_population == "membrane" and CONTACT.d0 is None
    assert not set(CONTACT.blocked_by) & set(BLOCKED_BY), (
        "the constraint is true under both readings of the ratchet; it may not wait on the bond's "
        "questions, and the only thing it waits on is a contact distance"
    )
    assert "2.7 nm" in CONTACT.blocked_by[0], "the delta-is-not-d0 trap must travel with the blocker"
    fam = ContactFamily(gap_law=GAP_LAW_POINT_SURFACE, scopes=(CONTACT,))
    assert fam.is_complementarity and not fam.enforceable

    print(f"{CONNECTOR}: constraint DECLARED on {GAP_LAW_POINT_SURFACE} (blocked on d0 only); "
          f"bond REFUSED, blocked on {len(SPEC.blocked_by)}; canonical half of the pair "
          "with filopodium_membrane_tip, which writes the DUPLICATE_OF entry")


if __name__ == "__main__":
    _demo()
