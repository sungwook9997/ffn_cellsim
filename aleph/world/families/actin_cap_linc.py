r"""``actin_cap_linc`` — the perinuclear actin cap on the nuclear envelope, and why it does not build.

WHAT THIS CONNECTOR IS.  The perinuclear actin cap is a set of contractile actomyosin bundles that arch
OVER the apical face of the nucleus and terminate on the nuclear envelope, where each end is held by a
LINC complex: a KASH-domain nesprin in the outer nuclear membrane, bridged across the perinuclear space
to a SUN trimer in the inner membrane, which in turn grips the lamin meshwork.  Mechanically it is the
one path by which an actomyosin bundle can put a direct, non-steric load on the nucleus.  It is
therefore a BOND in this architecture's sense — a stateful, identified pairing that can break — and not
a field interaction.

**THIS MODULE BUILDS NOTHING, AND THAT REFUSAL IS ITS PRODUCT.**  :data:`SPEC` carries six unanswered
questions and :func:`build_actin_cap_linc` raises
:class:`~aleph.world.families.ConnectorGapError` naming them.  What is here besides the refusal is the
one decision that is genuinely this module's and is not blocked: how a candidate pair is chosen, and
from which side the count is allowed to come.

----------------------------------------------------------------------------------------------------
THE ERM TEST, APPLIED.  Is this connector a mesh number wearing a physiological label?
----------------------------------------------------------------------------------------------------

``bond.py`` records the incumbent ERM defect: one declared connector standing for a population, so
nobody had to answer how many, and the answer turned out to be the icosphere vertex count.  Run the
same test here.  In the frozen port source, ``engine/sf_population.py`` places the cap as
``for k in range(n_cap)`` with **``n_cap: int = 2``** — a signature default — and gives each cap bundle
``linc_local = [outer_a, outer_b]``, both outer sarcomere ends.  So the incumbent's answer to *how many
actin-cap LINC bonds does a cell have* is **four**, and four is ``2 x 2``: a typed-in fibre count times
"both ends of a chain".  It is not the ERM defect exactly — the count is not set by the ENVELOPE mesh,
so it does not move when the nuclear subdivision changes — but it is the same class of number.  It is
a BUILDER PARAMETER wearing a physiological label, and it may not be carried forward.

⚠ Note which way round the pairing runs, because it is the entire defence.  Candidates here are
generated **from the filament side**: one candidate per cap-bundle endpoint, each paired to the
envelope vertex nearest it (:func:`linc_candidate_pairs`).  Generating them from the envelope side —
"a LINC at every nuclear vertex" — is precisely how ERM acquired a mesh-set population, and it would
be indistinguishable from a physiological answer at any single subdivision.

----------------------------------------------------------------------------------------------------
DUPLICATE JUDGEMENT: **NOT** a duplicate of ``mt_nucleus_linc`` or ``if_nucleus_linc``
----------------------------------------------------------------------------------------------------

The three edges are declared as one ``ConnectorFamily.LINC`` in the port source and the round-2 brief
calls them "LINC 3형제", so the duplication question is live and must be answered from physics rather
than from the shared name.  It is answered NO, on three grounds:

1. **The load-bearing element is not shared.**  A ``BondFamily`` bond is a spring between two global
   nodes with a rest length and a stiffness.  The SUN trimer and the perinuclear bridge ARE common to
   all three, but they are the short, stiff end of the path; the compliance and the rest length are set
   by the cytoplasmic KASH protein, and those are three different molecules whose contour lengths
   differ by more than an order of magnitude.  Nesprin-2 giant is a spectrin-repeat rod hundreds of
   nanometres long that unfolds under load; nesprin-3 is short and reaches vimentin only through
   plectin, which inserts a second compliance in series.  Two springs whose rest lengths differ by 10x
   are not one family with a shared card.
2. **One of the three is not a passive spring at all.**  The MT edge is a motor attachment — the
   nuclear KASH protein recruits kinesin/dynein, and the force it delivers depends on motor state, not
   on extension.  Folding it in would put a motor and a tether under one rest length.
3. **The package criterion decides it.**  ``families.DUPLICATE_CRITERION`` — *"chemistry card; a
   shared runtime class or force law is not a duplicate"* — and ``bond.py``'s *"two declared edges
   sharing a card are ONE family"*.  The three do not share a card even in the port source, which
   spells them ``nesprin_actin_linc`` / ``microtubule_motor_linc`` / ``nesprin3_plectin_if_linc``.
   ⚠ The three edges DO share one runtime class and one joint SoA in the port source, which is what
   makes them look folded; that is a shared CONTAINER.  Session D4 (``if_nucleus_linc``) reached the
   same verdict from that direction, independently, and it is recorded in ``DUPLICATE_CRITERION``.

The genuine NMII precedent points the other way and is worth stating so the asymmetry is not mistaken
for inconsistency: "membrane-side NMII" and "SF-side NMII" ARE one family, because the crossbridge
chemistry is identical and only the population it lands on differs — and the population is derived from
arena ID ranges, never stored.  Here the chemistry itself differs.  ``DUPLICATE_OF`` therefore takes no
entry from this module — and note that a NO here buys nothing towards the count, since
``DUPLICATE_CRITERION``'s own closing warning is that a duplicate verdict never answers *how many*.

⚠ **BUT A REAL COUPLING BETWEEN THE THREE IS UNDECLARED, AND IT IS NOT A DUPLICATION.**  If a SUN
trimer can be occupied by one KASH protein at a time, the three families COMPETE for one finite
population of SUN sites on the inner membrane.  Three independent ``BondFamily`` populations, each
sizing itself against the same envelope, would over-subscribe that resource silently — structurally the
same defect as PI decision card **B1** (one NMII head bound by two motor connectors double-counts its
force).  Nothing in the arena expresses a shared site budget today.  This is raised, not resolved:
question 6 of :data:`BLOCKED_BY`.

----------------------------------------------------------------------------------------------------
WHY THE COUNT CANNOT BE ANSWERED HERE
----------------------------------------------------------------------------------------------------

``BondCount`` demands *how many, against what, for which cell, on whose authority*.  Taken in order:

* **How many** — the knowledge base has the row and the row has no number.  ``KB-DRAFT-3.B-04``
  ("LINC complex — nucleus<->cytoskeleton force transfer") records the LINC surface density in its own
  words as ``condition/lamin-A-dependent per um^2``, at ``status: draft``.  A density that is declared
  condition-dependent and left unvalued is not a value with wide error bars; it is an axis nobody has
  declared.
* **Against what** — areal (per µm² of envelope) and per-filament (per cap-bundle end) are both
  defensible, and they are not interconvertible without the cap-BUNDLE count, which is itself absent:
  ``world/build/stress_fiber.py`` states that *"how many stress fibres a cell has is absent entirely"*.
* **Which cell** — every LINC and actin-cap datum reachable from here is fibroblast or epithelial.  The
  target is MCF7, and the actin cap is reported disrupted or absent in carcinoma lines with altered
  lamin A/C.  The physiologically correct count for MCF7 may be **zero**, and zero is a decision, not a
  default — ``BondCount.resolve`` returns a zero count without complaint precisely so that answer can
  be given deliberately.
* **On whose authority** — nothing is ratified.  ``forces_manifest.yaml`` marks
  ``nucleus_linc_stiffness`` ``pi_undeclared`` with the reason *"the sourced 8 pN LINC datum is a
  TENSION, not a stiffness (Dejardin)"*, and PI decision card **C1** (2026-08-11) records that no
  validated LINC on/off-rate card exists and that binding must be surfaced to the PI first.

⚠ **AND THE FILAMENT-SIDE POPULATION IS NOT IN THE ARENA.**  The round-2 brief opens this connector on
the ground that both populations exist.  They do not.  ``world/build/stress_fiber.py`` builds VENTRAL
bundles only — FA-to-FA, in the basal plane, with an axis refused if it is parallel to the substrate
normal — and a perinuclear cap arches over the nucleus.  There is no cap endpoint in the arena to bond
to, so this connector is blocked ahead of its count as well as by it.  Reported to Lead rather than
worked around: :func:`linc_candidate_pairs` therefore takes endpoints as an argument and never assumes
a population.

⚠ **NO KINETIC LAW IS WRITTEN HERE, AND THIS FAMILY IS A CANDIDATE FOR THE FIRST THAT HAS ONE.**
``bond.py`` defers attach/detach, free lists and snapshot twins to *"the first family that has
kinetics"*.  The port source gives this edge ``{bind, unbind}`` event channels and a source-gated
kinetics object that refuses to propose without a rate card.  Per the charter that is a PI question, so
it is escalated and not built.

Sanity Gate (recorded per CLAUDE.md, before first execution):
  * dimensional — ``rest_um`` [µm]; ``stiffness_pn_per_um`` [pN/µm]; ``support_um2`` [µm²] and the
    count's ``areal`` basis carry the matching unit, which ``BondCount`` checks by name.  Endpoint
    positions are [µm] and must be in the same frame as the envelope's, since the pairing is metric.
  * boundary — a resolved count of zero is a legitimate answer and is returned, not raised (it is the
    live possibility for MCF7).  A count exceeding the candidate endpoints raises, with both numbers,
    rather than silently shortening the population.  An envelope with no vertices is refused.
  * conservation/invariant — every bond addresses a node inside the arena's live NODE prefix
    (asserted by :func:`~aleph.world.bond.build_bond_family`), and every candidate pair joins one
    filament-side node to one envelope vertex, so no pair is a self-loop and none is envelope-internal.
    The candidate count is the number of ENDPOINTS supplied, never the number of envelope vertices —
    asserted in :func:`_demo`, because that is the ERM failure made numerical.
  * CFL/precision — no integration; float64.  ``rest_um`` MUST be a molecular contour/resting length,
    never the built separation ``|x_i - x_j|``: setting the rest length to the distance the builder
    happened to produce turns a hundreds-of-nanometre extensible tether into a zero-force weld, and it
    is a weld that no gate downstream can see.  This module therefore refuses to default it.
  * sign sense — a tether pulls the nucleus toward the fibre and the fibre toward the nucleus with
    equal and opposite force; the bond stores no direction, and the pair order (filament, envelope) is
    bookkeeping only.  ``rest_um`` is a length and must be non-negative.
  * measurement protocol — host-side construction; no device is touched and nothing is uploaded.  The
    support must be the MEASURED envelope area ``Surface.area0_total_um2``, summed over the mesh's own
    triangles, never ``4 pi R^2`` — the discrete surface differs from the analytic one and substituting
    the analytic value is what lets an areal density drift with resolution.

engine units: length µm, stiffness pN/µm, area µm².  Runtime: pure host bookkeeping, CPU-importable.
"""

from __future__ import annotations

import numpy as np
import numpy.typing as npt

from aleph.world.bond import BondCount, BondFamily, SourceClass, build_bond_family
from aleph.world.families import ConnectorGapError, FamilySpec
from aleph.world.surface import Surface

__all__ = [
    "CONNECTOR", "CHEMISTRY_CARD", "BLOCKED_BY", "SPEC",
    "linc_candidate_pairs", "build_actin_cap_linc",
]

#: The declared connector name, exactly as the device-run census spells it.
CONNECTOR = "actin_cap_linc"

#: The chemistry that identifies this family. Distinct from the MT and IF LINC cards on purpose — see
#: the duplicate judgement in the module docstring; the shared SUN bridge is not the compliant element.
CHEMISTRY_CARD = "nesprin2g_sun_actin_cap_linc"

#: The questions that must be answered before a family exists. Named, so Lead can aggregate them into
#: the decision queue without reading prose, and so a partially-answered set is visibly partial.
BLOCKED_BY: tuple[str, ...] = (
    "count.value — HOW MANY: no value exists. KB-DRAFT-3.B-04 records the LINC surface density in its "
    "own words as 'condition/lamin-A-dependent per um^2', status draft. The port source's answer is "
    "n_cap=2 bundles x 2 ends = 4, which is a builder signature default, not a measurement.",

    "count.basis — AGAINST WHAT: areal (per um^2 of nuclear envelope) or per_filament (per cap-bundle "
    "end) is undeclared. They are not interconvertible, because the cap-BUNDLE count is itself absent "
    "(world/build/stress_fiber.py: 'how many stress fibres a cell has is absent entirely').",

    "count.scope — WHICH CELL: every reachable LINC / actin-cap datum is fibroblast or epithelial; the "
    "target is MCF7 and the actin cap is reported disrupted in carcinoma lines with altered lamin A/C. "
    "Zero is a live and legitimate answer, and it must be given deliberately rather than defaulted to.",

    "count.provenance / rest_um / stiffness_pn_per_um — WHOSE AUTHORITY: forces_manifest marks "
    "nucleus_linc_stiffness pi_undeclared ('the sourced 8 pN LINC datum is a TENSION, not a stiffness "
    "(Dejardin)'), and PI card C1 records that no validated LINC on/off-rate card exists. The nesprin-2 "
    "giant resting length has no registered value either, and it may NOT be taken from the built "
    "separation — that substitution turns the tether into a weld.",

    "population — THE FILAMENT SIDE DOES NOT EXIST: world/build/stress_fiber.py builds VENTRAL bundles "
    "only (basal plane, axis refused parallel to the substrate normal). No perinuclear cap population "
    "is in the arena, so there is no endpoint to bond to. Blocked ahead of the count, not only by it.",

    "SUN-site exclusivity — UNDECLARED COUPLING: if one SUN trimer holds one KASH protein at a time, "
    "actin_cap_linc / mt_nucleus_linc / if_nucleus_linc compete for one finite site population on the "
    "inner nuclear membrane. Three families sizing independently over-subscribe it silently — the same "
    "defect as PI card B1 for NMII heads. Nothing in the arena expresses a shared site budget.",
)

#: This module's census entry. ``buildable`` is False and stays False until :data:`BLOCKED_BY` empties.
#: The basis is recorded as ``areal`` because that is what KB-DRAFT-3.B-04's own wording implies; it is
#: nonetheless listed above as unanswered, since the per-filament reading is equally live and the row
#: that would settle it carries no value.
SPEC = FamilySpec(
    connector=CONNECTOR,
    populations=("stress_fiber", "nuclear_envelope"),
    basis="areal",
    blocked_by=BLOCKED_BY,
    source="",
)


def linc_candidate_pairs(
    endpoint_nodes: npt.ArrayLike,
    endpoint_pos: npt.ArrayLike,
    envelope: Surface,
) -> npt.NDArray[np.int64]:
    """Candidate ``(K, 2)`` global node pairs, one per FILAMENT-side endpoint, nearest envelope vertex.

    The direction of generation is the point of this function.  One candidate per cap endpoint means
    the candidate count is a property of the ACTIN population; generating from the envelope side — a
    LINC per nuclear vertex — would make the count a property of the nuclear subdivision, which is the
    ERM defect, and it would look identical to a physiological answer at any one resolution.

    Selection is deterministic and total: the pairs come back in endpoint order, so
    :func:`~aleph.world.bond.build_bond_family` taking a prefix of them is the caller's ordering choice
    and never a draw.  Several endpoints may land on one envelope vertex; that is a fact about the
    geometry, and it is not deduplicated here because a vertex is a mesh sample, not a LINC socket.

    Args:
        endpoint_nodes: ``(K,)`` GLOBAL node indices of the filament-side endpoints.
        endpoint_pos: ``(K, 3)`` positions of those endpoints [µm], in the envelope's frame.
        envelope: the built nuclear envelope surface.

    Returns:
        ``(K, 2)`` GLOBAL index pairs, column 0 the filament node, column 1 the envelope vertex.

    Raises:
        ValueError: if the two endpoint arrays disagree in length, or the envelope has no vertices.
    """
    nodes = np.asarray(endpoint_nodes, np.int64).reshape(-1)
    pos = np.asarray(endpoint_pos, np.float64).reshape(-1, 3)
    if nodes.size != pos.shape[0]:
        raise ValueError(f"{nodes.size} endpoint nodes against {pos.shape[0]} positions")
    if envelope.n_vertices == 0:
        raise ValueError("the envelope has no vertices; there is nothing to anchor to")
    if nodes.size == 0:
        return np.zeros((0, 2), np.int64)

    d2 = ((pos[:, None, :] - envelope.position[None, :, :]) ** 2).sum(axis=2)
    nearest = envelope.nodes.lo + np.argmin(d2, axis=1)
    return np.stack([nodes, nearest.astype(np.int64)], axis=1)


def build_actin_cap_linc(
    arena: object,
    *,
    envelope: Surface,
    endpoint_nodes: npt.ArrayLike,
    endpoint_pos: npt.ArrayLike,
    count: BondCount | None = None,
    support_um2: float | None = None,
    rest_um: float | None = None,
    stiffness_pn_per_um: float | None = None,
) -> BondFamily:
    """Build the actin-cap LINC family — **or refuse, which is what it does today**.

    Every physiological argument is keyword-only and defaults to ``None``, so the refusal is the
    default path and supplying a value is a deliberate act by a caller who has the PI's answer.

    Args:
        arena: the world to claim the BOND range from.
        envelope: the built nuclear envelope; its measured ``area0_total_um2`` is the only admissible
            support for an ``areal`` count.
        endpoint_nodes: ``(K,)`` GLOBAL filament-side endpoint node indices.
        endpoint_pos: ``(K, 3)`` positions of those endpoints [µm].
        count: the PI-ratified :class:`~aleph.world.bond.BondCount`.  **Required to build.**
        support_um2: the MEASURED envelope area [µm²] the density resolves against.  Pass
            ``envelope.area0_total_um2``, never ``4 pi R^2``.
        rest_um: the tether's resting length [µm].  **Required**, and it is a molecular length — never
            the built separation, which would weld the nucleus to the fibre at zero force.
        stiffness_pn_per_um: the tether stiffness [pN/µm].  **Required**; the 8 pN LINC datum is a
            tension and may not be substituted for it.

    Returns:
        The built :class:`~aleph.world.bond.BondFamily`.

    Raises:
        ConnectorGapError: whenever any required physiological answer is absent — the present state.
    """
    given = {
        "count": count, "support_um2": support_um2,
        "rest_um": rest_um, "stiffness_pn_per_um": stiffness_pn_per_um,
    }
    missing = [name for name, value in given.items() if value is None]
    if missing:
        raise ConnectorGapError(
            CONNECTOR,
            [f"{name} was not supplied" for name in missing] + list(BLOCKED_BY),
            note=("These are PI questions, not build details. Supplying a value here would repeat the "
                  "ERM defect with a different label — see this module's docstring."),
        )

    pairs = linc_candidate_pairs(endpoint_nodes, endpoint_pos, envelope)
    return build_bond_family(
        arena, CONNECTOR, chemistry_card=CHEMISTRY_CARD, count=count, support=float(support_um2),
        pairs=pairs, rest_um=float(rest_um), stiffness_pn_per_um=float(stiffness_pn_per_um),
    )


def _demo() -> None:
    """Self-check: the refusal holds, the count comes from the actin side, and the pair is derived."""
    from aleph.world.arena import Kind, WorldArena
    from aleph.world.families import DUPLICATE_OF
    from aleph.world.surface import build_surface, icosphere

    # The census entry says blocked, and says so by NAME rather than by a flag.
    assert not SPEC.buildable
    assert len(SPEC.blocked_by) == 6
    assert any(q.startswith("count.value") for q in SPEC.blocked_by)
    assert any("DOES NOT EXIST" in q for q in SPEC.blocked_by)

    # The LINC trio are NOT folded: the compliant element differs, so the cards differ.
    assert CONNECTOR not in DUPLICATE_OF, "the duplicate judgement is NO — see the module docstring"

    arena = WorldArena(capacity={Kind.NODE: 4_000, Kind.FACE: 4_000, Kind.ANGLE4: 8_000, Kind.BOND: 500})
    v, f = icosphere(2)
    env = build_surface(arena, "nuclear_envelope", vertices=v, faces=f, radius_um=5.1)

    # Four cap endpoints, arching over the envelope. Positions are a self-check fixture, not a build.
    ends = arena.claim("cap_stub", Kind.NODE, 4)
    ep_nodes = ends.lo + np.arange(4)
    ep_pos = np.array([[-5.4, 0.0, 1.0], [5.4, 0.0, 1.0], [0.0, -5.4, 1.0], [0.0, 5.4, 1.0]])

    # THE ERM TEST, MADE NUMERICAL: the candidate count follows the ACTIN side, not the mesh.
    pairs = linc_candidate_pairs(ep_nodes, ep_pos, env)
    assert pairs.shape == (4, 2), "one candidate per endpoint"
    assert env.n_vertices == 162 and pairs.shape[0] != env.n_vertices
    finer = build_surface(arena, "envelope_fine", vertices=icosphere(3)[0], faces=icosphere(3)[1],
                          radius_um=5.1)
    assert linc_candidate_pairs(ep_nodes, ep_pos, finer).shape[0] == 4, \
        "refining the envelope must not change how many LINC bonds a cell has"

    # Nearest is nearest, checked against brute force, and every pair joins the two sides.
    for k in range(4):
        d = np.linalg.norm(env.position - ep_pos[k], axis=1)
        assert pairs[k, 1] == env.nodes.lo + int(np.argmin(d))
    assert np.all(pairs[:, 0] != pairs[:, 1])
    assert np.all((pairs[:, 1] >= env.nodes.lo) & (pairs[:, 1] < env.nodes.lo + env.n_vertices))

    # A mismatched candidate request is refused by name rather than broadcast into a wrong answer.
    try:
        linc_candidate_pairs(ep_nodes[:2], ep_pos, env)
    except ValueError as exc:
        assert "endpoint nodes against" in str(exc), exc
    else:  # pragma: no cover
        raise AssertionError("a mismatched candidate request must refuse")

    # No endpoints is zero candidates, not an error: a cell with no cap is a physical configuration.
    assert linc_candidate_pairs([], np.zeros((0, 3)), env).shape == (0, 2)

    # THE PRODUCT: the builder refuses, and the refusal carries the questions.
    try:
        build_actin_cap_linc(arena, envelope=env, endpoint_nodes=ep_nodes, endpoint_pos=ep_pos)
    except ConnectorGapError as exc:
        assert exc.connector == CONNECTOR
        assert len(exc.missing) == 4 + len(BLOCKED_BY)
        assert "count was not supplied" in exc.missing
        assert "rest_um was not supplied" in exc.missing
    else:  # pragma: no cover
        raise AssertionError("actin_cap_linc has no count and must refuse to build")

    # Partial answers do not unblock it either.
    try:
        build_actin_cap_linc(arena, envelope=env, endpoint_nodes=ep_nodes, endpoint_pos=ep_pos,
                             rest_um=0.8, stiffness_pn_per_um=100.0)
    except ConnectorGapError as exc:
        assert "count was not supplied" in exc.missing and "rest_um was not supplied" not in exc.missing
    else:  # pragma: no cover
        raise AssertionError("a partly-answered family must still refuse")

    # And once a count EXISTS, the family builds and its component pair is DERIVED, never declared.
    # The values below are a self-check fixture and are labelled as one; they are not an answer.
    fixture = BondCount(
        basis="areal", value=4.0 / env.area0_total_um2,
        scope="SELF-CHECK FIXTURE - not a cell, not a condition, not a temperature",
        source_class=SourceClass.PI_GAP,
        provenance="fixture proving the builder works once the PI answers; NOT a physiological density")
    fam = build_actin_cap_linc(
        arena, envelope=env, endpoint_nodes=ep_nodes, endpoint_pos=ep_pos, count=fixture,
        support_um2=env.area0_total_um2, rest_um=0.8, stiffness_pn_per_um=100.0)
    assert fam.n_bonds == 4 and fam.chemistry_card == CHEMISTRY_CARD
    assert fam.component_pairs(arena) == {("cap_stub", "nuclear_envelope"): 4}, \
        "the population pair is read off arena ID ranges, not stored on the bond"
    assert fam.provenance_row()["source_class"] == "PI_GAP"

    # A density low enough to give no bonds is an ANSWER — the live possibility for MCF7.
    zero = BondCount(basis="areal", value=1e-12, scope="SELF-CHECK FIXTURE",
                     source_class=fixture.source_class, provenance="the zero-count boundary case")
    assert build_actin_cap_linc(arena, envelope=env, endpoint_nodes=ep_nodes, endpoint_pos=ep_pos,
                                count=zero, support_um2=env.area0_total_um2, rest_um=0.8,
                                stiffness_pn_per_um=100.0).n_bonds == 0

    print(f"actin_cap_linc self-check OK — BLOCKED on {len(BLOCKED_BY)} PI questions, 0 bonds built")


if __name__ == "__main__":
    _demo()
