"""``membrane_cortex_contact`` — the compressive membrane/cortex channel. A CONSTRAINT, and NOT D1.

⚠ **RULED 2026-08-21, PI decision 1** (``docs/v2_audit/PI_DECISIONS_2026-08-21.md``): this edge is a
non-penetration CONSTRAINT, not a bond. Everything below was written before that ruling and is kept
because the ruling agreed with it — the argument is what produced the decision, so deleting it would
delete the evidence and leave only the verdict.

**What changed: the destination.** The refusal below closed with *"``aleph/world/`` has no home for a
non-penetration field"*, and that is now false: :mod:`aleph.world.contact` is the home. So this module
no longer refuses into a void. It declares :data:`CONTACT` — a :class:`~aleph.world.contact.ContactScope`
on the shared :data:`~aleph.world.contact.GAP_LAW_POINT_SURFACE` — and **the three count questions
below are not deferred, they are DISSOLVED**: a constraint has no population, so *how many* and
*against what* were questions about the wrong object. See :data:`CLOSED_BY_RULING`.

**What did not change: the two magnitudes.** ``d0`` is still registered nowhere, and the ERM
ownership-invariant conflict (⚠ below) is still a decision nobody has taken. Those survive verbatim.

TWO VERDICTS, BOTH REFUSALS OF A BOND. This module builds no :class:`~aleph.world.bond.BondFamily`. It
refuses twice, for two different reasons, and the second reason is the one that matters.

VERDICT 1 — it is NOT a duplicate of ``membrane_erm_cortex`` (D1). Under the criterion Lead recorded in
``families/__init__.py`` — *a family's identity is its CHEMISTRY CARD* — this is decided before the
mechanics are even read: D1 declares ``chemistry_card="ezrin_membrane_f_actin"`` and this edge declares
none at all (``contracts.py:702`` vs ``:981``, read off the built architecture, not off the source).
Two edges cannot be one family when one of them has no family identity.

The mechanics say the same thing independently, which is why both are recorded. The two edges join the
same population pair and are still different laws, because they occupy DISJOINT SIGN DOMAINS of the
same gap:

  * ERM is a molecular tether. ``erm_cortex_connector.py:33-34`` records in its own Sanity Gate that
    "a compressed / broken / resting tether is force-free (a molecular linker is a tether, not a
    strut)". It carries gap > 0 and nothing else.
  * CONTACT is a unilateral non-penetration law. ``connector_joints.py:175-181`` is force-free unless
    ``gap < 0``. It carries gap < 0 and nothing else.

Folding them into one family would produce a bilateral spring over the whole gap axis, and
``connector_joints.py:27-30`` states why that is not a parameter choice but a different object: "a
contact that pulls is an adhesion nobody declared. The sign is therefore structural, gated on the gap,
and cannot be turned off by a stiffness value." ``contracts.py:975-980`` reaches the same conclusion
from the load-path side and is the ratified one — PI decision D3 option C, 2026-07-28 — making this
edge the turgor transmission path precisely BECAUSE the ERM tether structurally cannot supply it. So
``DUPLICATE_OF`` gets no entry from here, and there is no second family either.

⚠ FLAG FOR LEAD — the two port sources contradict each other, and this module cannot settle it.
``erm_cortex_connector.py:32`` asserts as an ownership invariant that "the tether is the ONLY
membrane<->cortex mechanical coupling", which the 2026-07-28 ratification of this edge falsifies. The
ERM line is the older of the two and is presumed stale, but it is an ASSERTED INVARIANT in a frozen
port source, so which one the arena inherits is a decision, not a reading.

VERDICT 2 — it is not a BOND AT ALL, so ``BondCount``'s four questions have no answers to look for.
The missing chemistry card is not an omission to fill in; it is the finding. ``BondFamily`` takes
``chemistry_card`` as a required field precisely because that is what a family IS, and ``bond.py``
draws the line the same way: a bond is "a stateful, identified pairing with kinetics", while "field
interaction (steric, drag) needs no declaration because it depends on position alone".

That the cardless contacts are exactly the position-only ones is MEASURED, and the measurement is
SCOPED. Within ``ConnectorFamily.CONTACT`` — the five edges ``connector_joints.py:7-13`` lists — the
three properties coincide exactly, card iff kinetics iff not ``remap_on_accept``, splitting 2 / 3::

    lamellipodium_membrane_contact  card=actin_membrane_brownian_ratchet_...   kinetics=True   remap=False
    filopodium_membrane_tip         card=actin_membrane_brownian_ratchet_...   kinetics=True   remap=False
    membrane_cortex_contact         card=None                                  kinetics=False  remap=True
    membrane_ecm_contact            card=None                                  kinetics=False  remap=True
    nucleus_cortex_contact          card=None                                  kinetics=False  remap=True

(``membrane_erm_cortex``, for the D1 comparison, sits outside this family with
``card=ezrin_membrane_f_actin``, ``kinetics=True``, ``remap=False``.)

⚠ The coincidence is NOT architecture-wide and this module does not claim it is. 15 of the 38 declared
edges carry no card; 11 of them are ``immersed_transfer`` / ``fluid_boundary`` / ``environment_boundary``
edges, a different kind of field coupling and not this module's question. What IS exact architecture-wide
is the CONJUNCTION: cardless **and** ``remap_on_accept=True`` picks out those three contact edges and
nothing else in the cell.

So ``ConnectorFamily.CONTACT`` is not the boundary — ``connector_joints.py:19-22`` already says the
family label is not the physics boundary, and this is which side of it each edge lands on. The two
kinetic contacts are Brownian-ratchet chemistry and are bonds. These three have no chemistry, their
pair set is re-discovered by a broad phase on every accepted step, and nothing about a pair survives
the step that found it. There is no identity to keep, so there is nothing for a BOND to be.

⚠ UNRELATED FINDING, surfaced by taking that table and recorded rather than dropped:
``dorsal_arc_crosslink`` declares ``kinetics=True`` with NO chemistry card — the only edge in the
architecture that does. Under Lead's criterion a kinetic bond without a card has no family identity, so
whoever owns that edge has a question this module cannot answer for them.

⚠ FOR LEAD — the verdict is a CLASS verdict, and two of the three are other sessions' files.
``nucleus_cortex_contact`` is D6 and ``membrane_ecm_contact`` is one of the two unassigned FA-side
edges. If PI rules this edge out of the family census, the same ruling disposes of both, and it should
be taken once rather than three times.

THE ERM DEFECT WEARING A SECOND COSTUME — this is why the refusal is the product. The incumbent runtime
does not honour ``remap_on_accept``: ``connector_joints.py:47-54`` records, as a KNOWN GAP, that it
"holds a FIXED-CAPACITY pair list resolved at build". A frozen pair list has a count, and the endpoints
say what that count would be — ``contracts.py:982`` gives the membrane side as "live membrane contact
quadrature" and ``cortex_population.py:364-366`` gives the cortex side as material NODES. One pair per
membrane quadrature point IS the icosphere vertex count: 642 at subdiv 3, the same number every resting
run of 2026-08-20 measured. Building a family here would hand that number a ``BondCount``, a scope and a
provenance line, and it would then be a mesh number wearing a physiological label — the exact defect
``bond.py`` is shaped around, re-entered through a different connector. So the count is not merely
unknown here; asking for one is the wrong question, and answering it is how the defect returns.

WHAT WOULD CLOSE THIS. Not a density. Either (a) PI rules this edge is a field interaction and it
leaves the family census entirely, becoming a broad-phase contact the arena declares nothing for — the
reading this module argues for — or (b) PI rules it is a bond population, in which case the count needs
a physiological basis that is NOT the membrane discretisation, plus the two magnitudes nothing in the
port source defines (below). ``basis`` is recorded as ``explicit`` on the SPEC because that is the only
one of the four whose support unit is ``none``; it is a placeholder for "no support exists", not a
claim that an explicit count is known.

**PI took (a).** The SPEC stays blocked, permanently and by ruling rather than by an open question, and
``basis`` stays ``explicit`` as the same placeholder — the census enumerates this edge, and enumerating
it must not require inventing a basis for a count that does not exist.

Sanity Gate (recorded per CLAUDE.md, before first execution):
  * dimensional — nothing is computed. The two magnitudes this edge WOULD need are a contact distance
    [µm] and a stiffness [pN/µm] (``connector_joints.py:56-58``, where ``rest_um`` on a contact joint
    is the contact distance, not a rest length). Neither exists: ``connector_joints.py:43-46`` defines
    no magnitude by design, and a repo-wide search finds a ``contact_distance_um`` only on the AFM
    indenter, which is a different contact.
  * boundary — there is no count to be zero. A resolved zero would be an answer about a density; this
    module never reaches a density.
  * conservation/invariant — no arena state is touched, no claim is made, no capacity is consumed.
    :func:`build` raises before it reads the arena, so a failed import or a caught exception cannot
    leave a half-claimed BOND range behind.
  * CFL/precision — no integration, no device, no float.
  * sign sense — the sign IS the finding: gap < 0 here, gap > 0 for D1, and that disjointness is what
    refutes the duplication.
  * measurement protocol — host only; nothing is uploaded and no device is initialised. The card /
    kinetics / remap table above was read off ``reference_cell_architecture()`` as built rather than
    grepped out of the source, so a default that never reaches the object cannot be mistaken for a
    declaration. That read happened in the investigating session; ``world/`` may not import
    ``engine/``, so this module records the result and the demo cannot re-derive it.

engine units: length µm, stiffness pN/µm.  Runtime: pure host bookkeeping, CPU-importable.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, NoReturn

from aleph.world.contact import GAP_LAW_POINT_SURFACE, ContactScope
from aleph.world.families import ConnectorGapError, FamilySpec

if TYPE_CHECKING:  # pragma: no cover - annotation only; world/ never imports engine/
    from aleph.world.arena import WorldArena

__all__ = ["CONNECTOR", "SPEC", "CONTACT", "BLOCKED_BY", "CLOSED_BY_RULING", "build"]

CONNECTOR = "membrane_cortex_contact"

#: Questions that were open before PI decision 1 and are now CLOSED — recorded rather than deleted,
#: because each was closed by being about the wrong object and that is worth more than its answer.
#: Anyone re-opening a count question for this edge should read these first.
CLOSED_BY_RULING: tuple[str, ...] = (
    "IS THIS A BOND? — RULED: no. A non-penetration constraint (PI decision 1, 2026-08-21). Every "
    "line of evidence this module recorded held: no chemistry card, remap_on_accept=True, a pair set "
    "re-discovered per accepted step, and a unilateral law bond.py's bilateral primitive cannot state.",
    "how many — DISSOLVED, not answered. A constraint has no population: its active set is an OUTPUT "
    "of the broad phase at a configuration and changes every step. The icosphere-vertex count this "
    "module warned about is therefore not a count that is missing, it is a count that must never be "
    "asked for, and contact.ContactScope has nowhere to put one.",
    "against what — DISSOLVED with it. A support is what a DENSITY is resolved against; there is no "
    "density.",
)

#: What survives the ruling. Two entries, and neither is a count.
#:
#: ⚠ The first is why the SPEC stays blocked forever rather than becoming buildable: this module will
#: never build a :class:`~aleph.world.bond.BondFamily`, and the census must read that as "no family
#: here" and not as "a family nobody has got round to".
BLOCKED_BY: tuple[str, ...] = (
    "NOT A BOND, BY RULING — PI decision 1, 2026-08-21. This module builds no BondFamily and never "
    "will; its destination is aleph.world.contact and it declares CONTACT in this module. No "
    "BondCount question is well posed for it and none is queued. See CLOSED_BY_RULING.",
    "contact distance d0 [um] — the physiological membrane/cortex gap at which the force vanishes, "
    "WITH its scope. ABSENT and unchanged by the ruling: connector_joints.py defines no magnitude by "
    "design, the only contact_distance_um in the repo belongs to the AFM indenter, and the 0.030 um "
    "in the port source's tests is a fixture. This is the one physical magnitude a contact has.",
    "WHICH SOURCE STANDS — erm_cortex_connector.py:32 asserts the ERM tether is the ONLY "
    "membrane<->cortex mechanical coupling, which the 2026-07-28 ratification of this edge falsifies. "
    "An asserted invariant in a frozen port source is a decision to retire, not a line to ignore. "
    "Untouched by decision 1, which ruled on the SHAPE of this edge and not on whether it exists.",
)

#: Declared for the BOND census. Blocked by ruling, and the blocking is the point.
SPEC = FamilySpec(
    connector=CONNECTOR,
    populations=("membrane", "cortex"),
    basis="explicit",
    blocked_by=BLOCKED_BY,
)

#: **The destination.** One scope on the shared gap law — the membrane/cortex surface pair, tested by
#: a broad phase, carrying no pair list and no count. ``d0`` is ``None`` because it is unsourced, and
#: ``None`` is the census entry rather than a placeholder value.
#:
#: The FORMULATION is not declared here and must not be: whether the family is enforced as
#: ``0 <= Phi ⊥ gamma >= 0`` or as a regularised penalty is a property of the whole family and of the
#: solve, not of one surface pair. A scope that chose one would be one connector deciding for four.
#: ⚠ If the regularised form is taken, this edge's stiffness is NOT a free knob: it sets the turgor
#: transmission this edge was ratified to carry (PI D3 option C, 2026-07-28), so the approximation
#: error lands directly on a load path, which is where it must be reported.
#:
#: ⚠ **AND THE CORPUS HOLDS THE WRONG NUMBER IN ABUNDANCE.** A KB sweep on 2026-08-21 (12,810
#: chunks) found NO membrane/cortex separation and no ``parameter`` or ``knowledge_claim`` row for
#: one. What it did find, in five papers, is the cortex **THICKNESS**: 100–200 nm
#: (``Chugh2017_NatCellBiol``, ``Fritzsche2016_SciAdv``, ``Fritzsche2017_NatCommun``,
#: ``Lchtefeld2024_NatMethods``, ``mobility-of-molecular-motors``). That is a different quantity —
#: how thick the layer is, not how far it sits from the bilayer — and it is dimensionally correct,
#: abundantly sourced and plausible-looking, which is precisely what makes it the value most likely
#: to be promoted by someone who needs a number to make a run start. It is named in the blocker so
#: the trap travels with the gap rather than being rediscovered.
CONTACT = ContactScope(
    connector=CONNECTOR,
    point_population="cortex",      # filament NODEs
    surface_population="membrane",  # owns FACEs (build/membrane.py)
    blocked_by=BLOCKED_BY[1:] + (
        "⚠ NOT the cortex THICKNESS. 100-200 nm is well sourced in this corpus (Chugh2017, "
        "Fritzsche2016/2017, Lchtefeld2024) and is how thick the layer is, not how far it sits from "
        "the bilayer. A KB sweep 2026-08-21 over 12,810 chunks found no membrane/cortex separation "
        "and no registered parameter for one; the thickness is what is there instead.",
    ),
)


def build(arena: WorldArena) -> NoReturn:
    """Refuse to build a BOND, and name the destination that now exists.

    Args:
        arena: the world a family would have been claimed from. Never read — the refusal precedes it,
            so no capacity is consumed and no partial claim can be left behind. A constraint would not
            claim from it either: it is an address resolved from live geometry, not an allocation.

    Raises:
        ConnectorGapError: always. See the module docstring for both verdicts.
    """
    raise ConnectorGapError(
        CONNECTOR,
        list(BLOCKED_BY),
        note=(
            "RULED a constraint, not a bond (PI decision 1, 2026-08-21) — the destination is "
            "aleph.world.contact.ContactFamily and this module's scope is CONTACT. Still NOT a "
            "duplicate of membrane_erm_cortex: ERM carries gap > 0 (a tether, force-free in "
            "compression) and this edge carries gap < 0 (unilateral non-penetration) — disjoint sign "
            "domains, so folding them would make a bilateral spring, which connector_joints.py:27-30 "
            "calls an adhesion nobody declared."
        ),
    )


def _demo() -> None:
    """Self-check: the bond refusal survives the ruling, the destination exists, no count appears."""
    from aleph.world.families import DUPLICATE_OF

    assert not SPEC.buildable, "this module must not claim to be buildable"
    assert SPEC.connector == CONNECTOR
    assert SPEC.populations == ("membrane", "cortex")

    # The duplication was REFUTED, not deferred: no entry may appear for this connector.
    assert CONNECTOR not in DUPLICATE_OF, (
        "D2 is not a duplicate of D1 — disjoint gap sign domains; see the module docstring"
    )

    # build() refuses before touching the arena, so None is a legal argument and proves it.
    try:
        build(None)  # type: ignore[arg-type]
    except ConnectorGapError as exc:
        assert exc.connector == CONNECTOR
        assert "NOT A BOND, BY RULING" in exc.missing[0]
        assert "aleph.world.contact" in exc.note and "gap > 0" in exc.note
    else:  # pragma: no cover
        raise AssertionError("membrane_cortex_contact must refuse to build a bond")

    # THE DESTINATION. Blocked on d0 and the ERM conflict — and on nothing that is a count.
    assert CONTACT.connector == CONNECTOR
    # ROLES, not order: the membrane owns FACEs and is the surface; the cortex supplies nodes.
    assert CONTACT.surface_population == "membrane" and CONTACT.point_population == "cortex"
    assert not CONTACT.enforceable and CONTACT.d0 is None
    joined = " ".join(CONTACT.blocked_by)
    assert "contact distance d0" in joined
    assert "how many" not in joined and "against what" not in joined, (
        "a count question has re-entered the live queue; it was DISSOLVED, not deferred — "
        f"see CLOSED_BY_RULING ({len(CLOSED_BY_RULING)} entries)"
    )

    # The scope must be enforceable under the SAME gap law the other contact edges use, or the
    # identity criterion has been broken by declaring a private law for one surface pair.
    from aleph.world.contact import ContactFamily

    fam = ContactFamily(gap_law=GAP_LAW_POINT_SURFACE, scopes=(CONTACT,))
    assert fam.is_complementarity and not fam.enforceable

    print(f"{CONNECTOR} self-check OK — bond: refused by ruling ({len(BLOCKED_BY)} entries, "
          f"{len(CLOSED_BY_RULING)} closed); constraint: scope declared on {GAP_LAW_POINT_SURFACE}")


if __name__ == "__main__":
    _demo()
