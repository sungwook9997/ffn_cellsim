"""``nucleus_cortex_contact`` — STERIC, and now it has somewhere to be.

**The question this module was opened to answer:** is the nucleus–cortex interaction a steric field or
a bond? **Answer: steric** — and ⚠ **RATIFIED 2026-08-21 by PI decision 1**
(``docs/v2_audit/PI_DECISIONS_2026-08-21.md``), which rules this edge a non-penetration CONSTRAINT and
closes D6. Everything below predates the ruling and is kept because it produced it.

**What changed: the destination.** The refusal closed with *"``aleph/world/`` has no home for a
non-penetration field … the destination for this connector is a design decision in ``laws/``"*. It is
now :mod:`aleph.world.contact`, and this module declares :data:`CONTACT` on the shared
:data:`~aleph.world.contact.GAP_LAW_POINT_SURFACE`. The two count questions are **DISSOLVED rather
than answered** — see :data:`CLOSED_BY_RULING` — and with them the mesh-artefact trap the ⚠ below
measures, because there is no longer any object a per-vertex count could be attached to. The check
that measures it STAYS in :func:`_demo`, as a ratchet: if a vertex count is ever pinned and called a
density, it still fails.

**What did not change:** ``d0``. No perinuclear contact gap is registered for any cell line, the
0.030 µm in the port source's tests is a fixture, and the ruling did not supply one.

WHAT THE PORT SOURCE ACTUALLY DECLARES (read-only, never imported from here).
``engine/contracts.py:960`` declares ``nucleus_cortex_contact`` as ``ConnectorFamily.CONTACT``,
``nucleus``↔``cortex``, ``kinetics=False``, ``commit_on_accept=False``, with ``generation_required``,
``remap_on_accept`` and ``blocks_sleep_refine`` all True. ``engine/load_path.py:122-124`` classes it
``STERIC_CONTACT`` — *"unilateral steric non-penetration; force-free at and beyond the physiological
gap"* — alongside ``membrane_cortex_contact`` and ``membrane_ecm_contact``, and
``tests/ac/engine/test_connector_joints.py:227`` pins exactly those three as the steric set. The PI note
that introduced the edge (``contracts.py:950-953``) states the whole class is *"non-kinetic, all
force-free at the physiological gap, so none introduces a sourced magnitude beyond a numerical exclusion
stiffness derived from the existing WCA / derived-mobility pattern."*

WHY THAT SETTLES IT, IN THIS ARCHITECTURE'S OWN WORDS. ``bond.py`` opens by drawing exactly this line:
*"Field interaction (steric, drag) needs no declaration because it depends on position alone; a specific
persistent pairing does, because 'same world' cannot express 'still attached to THAT one'."* A contact
pair here has no identity to persist — ``generation_required`` + ``remap_on_accept`` mean the pair list
is regenerated from proximity and re-addressed on every accepted step. ``arena.py:76-78`` reached the
same conclusion from the allocation side when it closed the ``Kind`` set: *"``SITE`` is an address
resolved from live geometry, not an allocation."* A steric contact is a SITE. A bond is an allocation.
Nothing in this cell can be both.

⚠ IS THE COUNT A MESH NUMBER WEARING A PHYSIOLOGICAL LABEL? **Yes — and worse than the ERM case that
this package exists because of.** The only way to give this connector a bond population is one bond per
nuclear-envelope contact quadrature point, i.e. per envelope mesh vertex, which is precisely how the ERM
count became the icosphere vertex count. ``_demo`` computes it live: at the two envelope spacings this
repository can defend today the vertex count moves 40,962 → 163,842 and the implied areal density moves
125/µm² → 501/µm², a factor of four from a spacing that is **itself an open PI-GAP** — ``build/envelope.py``
records ``ENVELOPE_MESH_PI_GAP``, no lamin meshwork spacing is registered anywhere here, and the value in
use is an analogy to the *cortical* mesh. ERM at least borrowed a sourced band. This would borrow an
admitted gap and then report the result as a count.

⚠ DUPLICATE JUDGEMENT: **NOT ASSESSABLE, and the reason is the determination again.** The package's
ratified :data:`~aleph.world.families.DUPLICATE_CRITERION` is *"chemistry card; a shared runtime class or
force law is not a duplicate"*. Applied here it returns nothing to compare: all three steric edges carry
``chemistry_card=None`` in ``reference_cell_architecture()``. They do share a runtime class
(``ContractJointConnector``) and a force law (one unilateral kernel), and under the criterion neither of
those is a duplicate — the D4 LINC case is exactly that trap, three chemistries under one container.

But the ``None`` is not an oversight to be filled in. **A chemistry card names a chemistry, and a steric
exclusion has no molecule** — nothing binds, nothing unbinds, and the magnitude is a numerical exclusion
stiffness by the PI's own note. So the criterion for being a family cannot even be evaluated on this
connector, which is another independent line of evidence for STERIC: a bond family is identified by its
chemistry, and this edge has no chemistry to be identified by.

**No ``DUPLICATE_OF`` entry is written**, and not because the question is close. The key must point at a
family that exists; under the determination none of the three is one, and pointing at nothing would be
the name-based folding this package forbids, run in reverse. If the PI overrules the determination, the
first thing that ruling must supply is a chemistry card for a contact — at which point the duplicate
question becomes askable for the first time.

**✅ AND THE QUESTION IS NOW ASKABLE, WITHOUT A CARD.** The PI did not overrule; the PI ruled the other
way, and :data:`aleph.world.contact.CONTACT_IDENTITY` is the criterion that replaces the card for
things that have none:

    *A constraint's identity is its GAP FUNCTION — which primitive kinds are measured and by what
    metric. The surface pair is SCOPE, not identity; the force law is downstream and decides nothing.*

Applied here it returns a verdict where the bond criterion returned nothing to compare. All five
declared CONTACT edges compute the SAME Phi — ``connector_joints.py:7-13`` calls them *"the same
unilateral gap law"* in the source's own words — so they are **ONE constraint family with several
scopes**, not several families. This edge contributes the ``nuclear_envelope``/``cortex`` scope.

⚠ That is the same merge ``bond.py`` performs with a shared chemistry card, and it rejects the same
thing: the declared component pair, which is bookkeeping for a bond and scope for a constraint, and
identity for neither. It does NOT contradict ``families.DUPLICATE_CRITERION``'s *"a shared force law is
not a duplicate"* — that clause is about a shared IMPLEMENTATION downstream of a chemistry. Here there
is no chemistry, and Phi is not the force law: it is the geometry the force law is downstream OF.

⚠ **And folding does not dissolve d0.** ``DUPLICATE_CRITERION``'s own closing warning transfers
exactly: a duplicate verdict answers *"is this one family or two"* and never *"how many"* — here, never
*"how far apart"*. One family, four scopes, and the perinuclear gap is still unregistered.

WHERE THIS LEAVES A REAL GAP FOR LEAD — ✅ **CLOSED 2026-08-21.** ``aleph/world/`` had **no home for a
non-penetration field**: ``bond.py`` declares pairings, ``laws_bind.py`` binds per-primitive kernels,
and neither is a neighbour query; ``build/cortex.py:47-49`` notes the incumbent's WCA overlap
relaxation was deliberately not ported. :mod:`aleph.world.contact` is now that home. ⚠ It is a TYPE and
a CONTRACT, not a neighbour query: the broad phase itself is still unwritten, and this module's scope
declares where one would look, not how it looks.

Sanity Gate (recorded per CLAUDE.md, before first execution):
  * dimensional — nothing is computed. Had a family been built: an areal contact density is [1/µm²]
    against a support in [µm²], the contact gap ``d0`` is [µm] and the exclusion stiffness [pN/µm]. The
    gap and the stiffness are BOTH unsourced (below), so the check is vacuous and is recorded as vacuous
    rather than as passed.
  * boundary — the module's only executable path is a refusal, so its boundary case is that the refusal
    is total: :func:`build` raises for every input because it takes none.
  * conservation/invariant — none held here. The invariant that matters is negative and is checked by
    ``_demo``: no ``BondFamily`` is constructed and no arena claim is taken, so a blocked connector
    cannot consume BOND capacity that a real family will need.
  * CFL/precision — no integration, no device, float64 where arithmetic occurs (the mesh check only).
  * sign sense — the physics being declined is unilateral: overlap repels, separation is force-free, and
    a contact that pulls is an undeclared adhesion. That asymmetry is the second reason this is not a
    bond, which is bilateral by construction in ``bond.py`` (one ``rest_um``, one stiffness, no side).
  * measurement protocol — host-side only; nothing is uploaded and no device is touched.

Device-run census (``outputs/ac/connector_devicerun/native_record.json``): ``verdict: NOT_BOUND``,
``n_calls: 0``, ``force_channel: null`` — *"no runtime object in this composed world"*. Consistent with
the above: the edge has never run, and no measurement is being retracted by declining to build it.

engine units: length µm, stiffness pN/µm.  Runtime: pure host bookkeeping, CPU-importable.
"""

from __future__ import annotations

import math

from aleph.world.bond import BondFamily
from aleph.world.contact import GAP_LAW_POINT_SURFACE, ContactScope
from aleph.world.families import DUPLICATE_CRITERION, ConnectorGapError, FamilySpec

__all__ = ["CONNECTOR", "SPEC", "CONTACT", "BLOCKED_BY", "CLOSED_BY_RULING", "build"]

CONNECTOR = "nucleus_cortex_contact"

#: Closed by PI decision 1, 2026-08-21. Kept rather than deleted: two of the three were closed by
#: being about the wrong object, which is a stronger result than an answer would have been.
CLOSED_BY_RULING: tuple[str, ...] = (
    "PRIMARY — is a steric contact a bond at all? RULED: no, a non-penetration constraint (PI "
    "decision 1, closing D6). Every line of evidence this module recorded was upheld: "
    "load_path.py:122-124 STERIC_CONTACT; bond.py's exclusion of position-only field interaction; "
    "arena.py:76-78, a SITE is 'an address resolved from live geometry, not an allocation'. The "
    "destination is aleph.world.contact, declared as CONTACT below.",
    "how many — DISSOLVED. The pair list is regenerated from proximity every step "
    "(generation_required=True, remap_on_accept=True), so no persistent population exists for a count "
    "to be OF. That was already the finding; the ruling makes it structural — ContactScope has no "
    "field a count fits in.",
    "against what — DISSOLVED with it. The nucleus-cortex contact area is a STATE variable that "
    "changes as the cell deforms, and a support is what a DENSITY resolves against. There is none.",
)

#: What survives. ``d0`` is the only physical magnitude a contact has, and it is absent.
BLOCKED_BY: tuple[str, ...] = (
    "NOT A BOND, BY RULING — PI decision 1, 2026-08-21, closing D6. This module builds no BondFamily "
    "and never will; it declares CONTACT instead. See CLOSED_BY_RULING for what the ruling dissolved.",
    "which cell — no perinuclear contact gap d0 is registered for MCF7 or any other line. "
    "connector_joints.py:193 takes contact_distance_um as a caller argument and no caller in the "
    "repository supplies one; the 0.030 um that appears in tests is a fixture. Unchanged by the "
    "ruling, and now the ONLY physical magnitude this edge needs.",
    "on whose authority — contracts.py:950-953 states the magnitude is 'a numerical exclusion "
    "stiffness derived from the existing WCA / derived-mobility pattern'. That description is honest "
    "about a NUMERICAL parameter and dishonest about a physiological one, and contact.py keeps them "
    "apart: the complementarity form has no stiffness at all (gamma is a multiplier), and a penalty "
    "stiffness is declared as the regularisation it is. Which form the family takes is not this "
    "connector's call — it is one property of the whole family and of the solve.",
)

#: ``basis`` is required by :class:`FamilySpec` and no basis is right, because the count is not the
#: quantity in question. ``areal`` is recorded as the only one that could apply IF the PRIMARY question
#: is decided against this module — an areal density over the envelope surface — so the census entry is
#: enumerable. It asserts no value and no measurement.
SPEC = FamilySpec(
    connector=CONNECTOR,
    populations=("nuclear_envelope", "cortex"),
    basis="areal",
    blocked_by=BLOCKED_BY,
    source="",
)

#: **The destination.** The ``nuclear_envelope``/``cortex`` scope of the one gap law, carrying no pair
#: list and no count — which is exactly what makes the mesh-artefact trap unreachable rather than
#: merely unvisited. ``nuclear_envelope`` is the ratified spelling (PI decision 8, same sitting).
#:
#: ⚠ **THE ONE PAPER IN CORPUS THAT MODELS THIS CONTACT HAS NO d0 AT ALL, AND THAT IS THE FINDING.**
#: A KB sweep on 2026-08-21 (12,810 chunks) found exactly one: ``Fang2016_PhysRevE`` p4 gives a
#: cortex/envelope repulsion ``U = k(π/2·δ − tan(π/2·δ))`` for ``−1 ≤ δ < 0``, with ``k = 5e-14 J``
#: and
#:
#:     ``δij = (|r^C_i − r^N_j| − |r^C_0,i − r^N_0,j|) / |r^C_0,i − r^N_0,j|``
#:
#: — a *changed distance ratio* against the INITIAL configuration. So its gap is measured from
#: whatever the build happened to produce, which means the mesh sets the reference and the physics
#: inherits the discretisation: the ERM defect, arrived at from the law rather than from the count.
#: Its coefficient is a numerical potential constant, not a separation. ⚠ And
#: ``source_audit.verdict = CHECK`` (SE261), so under CLAUDE.md it may not be cited in a deliverable
#: as it stands. Recorded as the shape to AVOID, not as a source.
CONTACT = ContactScope(
    connector=CONNECTOR,
    point_population="cortex",                # filament NODEs
    surface_population="nuclear_envelope",    # NODE + FACE + ANGLE4 (build/envelope.py)
    blocked_by=BLOCKED_BY[1:] + (
        "⚠ NOT a gap measured against the initial configuration. Fang2016_PhysRevE p4 — the only "
        "cortex/envelope repulsion in this corpus — defines its gap as a RATIO of change against "
        "|r_0|, so the build sets the reference and the mesh sets the physics. Its 5e-14 J is a "
        "potential coefficient, not a separation, and source_audit says CHECK (SE261). A d0 must be "
        "an absolute separation with a cell scope, or it is the ERM defect reached via the law.",
    ),
)


def build() -> BondFamily:
    """Refuse to build a BOND, and name the destination that now exists. Never returns.

    Takes no arena and no parameters on purpose: a signature shaped for a build that the determination
    says should not happen is speculative scaffolding, and it would let a caller believe the only thing
    missing is an argument. That holds a fortiori after the ruling.

    Raises:
        ConnectorGapError: always.
    """
    raise ConnectorGapError(
        CONNECTOR,
        list(BLOCKED_BY),
        note=(
            "RULED a constraint, not a bond (PI decision 1, 2026-08-21, closing D6) — the destination "
            "is aleph.world.contact.ContactFamily and this module's scope is CONTACT. One contact per "
            "envelope mesh vertex would have reproduced the ERM defect on a mesh spacing that is "
            "itself an open PI-GAP (build/envelope.py: ENVELOPE_MESH_PI_GAP); under the ruling there "
            "is no count for it to become."
        ),
    )


# ── self-check ──────────────────────────────────────────────────────────────────────────────────
def _demo() -> None:
    """Self-check: the refusal is total, nothing is claimed, and the mesh-artefact check has teeth."""
    from aleph.world.build.membrane import icosphere_counts, subdivisions_for_mesh
    from aleph.world.contact import CONTACT_IDENTITY, ContactFamily

    # 1. The refusal names every open question and does not quietly become buildable.
    assert not SPEC.buildable
    # The duplicate judgement above was made under the ratified criterion; if that criterion is
    # re-ratified on something other than the chemistry card, this module's judgement must be re-made.
    assert "chemistry card" in DUPLICATE_CRITERION, DUPLICATE_CRITERION
    try:
        build()
    except ConnectorGapError as exc:
        assert exc.connector == CONNECTOR
        assert "NOT A BOND, BY RULING" in exc.missing[0]
        assert "aleph.world.contact" in exc.note and "ENVELOPE_MESH_PI_GAP" in exc.note
    else:  # pragma: no cover - build() has no return path
        raise AssertionError("a steric contact must not build a BondFamily")

    # 1b. THE DESTINATION, and that no count question came back with it.
    # ROLES, not order: build/envelope.py claims NODE + FACE + ANGLE4, so the envelope is the
    # surface. ⚠ The LAMINA could not be: build/lamina.py claims Kind.NODE and nothing else.
    assert CONTACT.surface_population == "nuclear_envelope"
    assert CONTACT.point_population == "cortex" and CONTACT.d0 is None
    joined = " ".join(CONTACT.blocked_by)
    assert "how many" not in joined and "against what" not in joined, (
        f"a count question re-entered the live queue; it was DISSOLVED — see CLOSED_BY_RULING "
        f"({len(CLOSED_BY_RULING)} entries)"
    )
    # The bond criterion could not be evaluated here; the constraint criterion can, and the merge it
    # licenses is what makes this ONE scope of a shared law rather than a private law for one pair.
    assert "gap function" in CONTACT_IDENTITY, CONTACT_IDENTITY
    fam = ContactFamily(gap_law=GAP_LAW_POINT_SURFACE, scopes=(CONTACT,))
    assert fam.is_complementarity and not fam.enforceable

    # 2. THE MESH-ARTEFACT CHECK. One bond per envelope contact quadrature point is one bond per mesh
    #    vertex, and the vertex count is set by a spacing this repository records as an open PI-GAP.
    #    If the density were physiological it would be invariant under that spacing. It is not: it
    #    quadruples. This assertion fails the moment someone pins a vertex count and calls it a count.
    r_nuc = 0.68 * 7.5  # N:C ratio x MCF7 central radius, per build/envelope.py — quoted, not decided.
    area = 4.0 * math.pi * r_nuc * r_nuc  # analytic, and labelled so: bond.py requires a MEASURED
    #                                       support, which is a further reason no count is available.
    counts = [icosphere_counts(subdivisions_for_mesh(r_nuc, m))[0] for m in (0.1, 0.05)]
    assert counts == [40962, 163842], counts
    coarse, fine = (v / area for v in counts)
    assert fine / coarse > 3.9, (
        "the per-vertex 'density' is supposed to be a discretisation artefact; if this ratio has "
        "become ~1 the count has been pinned and the ERM defect has been reproduced"
    )

    print(
        f"{CONNECTOR}: STERIC — a CONSTRAINT, not a family. blocked_by={len(SPEC.blocked_by)}; "
        f"per-vertex 'density' {coarse:.1f} -> {fine:.1f} /um^2 across two admissible envelope spacings"
    )


if __name__ == "__main__":
    _demo()
