r"""``if_nucleus_linc`` — the keratin cage's anchor on the nuclear envelope, and the count nobody has.

WHAT THIS CONNECTOR IS.  A LINC complex spanning the nuclear envelope: SUN1/2 in the inner membrane,
bound to the lamin meshwork on the nucleoplasmic side; nesprin-3 in the outer membrane, whose KASH
domain crosses the perinuclear space and whose cytoplasmic end binds **plectin**, which in turn binds
the intermediate filament.  Mechanically it is a passive elastic tether joining one IF spoke's inner
bead to one nuclear-envelope vertex.  It is the nuclear end of the load path
:mod:`aleph.world.build.intermediate_filament` exists to carry — that module builds the spokes and
stops, saying in its own docstring that "those nearest-neighbour searches are the LINC and cortex-anchor
BONDS — PHASE 3".  This is PHASE 3 for the nuclear end.

────────────────────────────────────────────────────────────────────────────────────────────────────
NOT A DUPLICATE, AND HERE IS THE ARGUMENT.

``actin_cap_linc``, ``mt_nucleus_linc`` and ``if_nucleus_linc`` share ONE runtime class
(``engine.linc_connector.LincRigConnectorAdapter``), one fixed-capacity joint SoA, and one module whose
first line is *"Graph-owned LINC joint population shared by actin-cap, MT, and IF edges"*.  Every
surface signal says duplicate.  They are not, for two reasons:

  1. **Three different chemistry cards.**  ``linc_connector._EDGE_CHEMISTRY`` is
     ``nesprin_actin_linc`` / ``microtubule_motor_linc`` / ``nesprin3_plectin_if_linc`` — three
     distinct molecular complexes, not one complex landing on three populations.  ``bond.py`` makes the
     chemistry card the family identity, and ``families.DUPLICATE_CRITERION`` records the rule three
     sessions converged on independently: *a shared runtime class is not a duplicate*.  The shared SoA
     is a shared CONTAINER.
  2. **One of the three is not even passive.**  ``mt_nucleus_linc`` couples through nesprin-4 / KASH5
     to kinesin-1 and dynein: it is a MOTOR attachment with a duty cycle, a stall force and a
     direction.  ``if_nucleus_linc`` is a passive spring through a cytolinker.  A family that generates
     force and a family that stores it cannot share a kinetic law.

Contrast the case the package warns about.  "Membrane-side NMII" and "SF-side NMII" ARE one family:
same head, same crossbridge cycle, same card — only the partner population differs, and the partner is
DERIVED from arena ID ranges rather than declared.  Here the molecule itself differs.  Folding these
three would hide three chemistries under one name, which is the ERM defect run backwards: there, one
name hid a population; here, one name would hide three chemistries.

Verdict: **NOT in** :data:`~aleph.world.families.DUPLICATE_OF`.  This module builds its own family.

────────────────────────────────────────────────────────────────────────────────────────────────────
⚠ THE COUNT IS THE ERM DEFECT AGAIN, ONE LAYER DEEPER — AND THIS TIME THE BOND IS NOT EVEN THE RIGHT
BOND.

``components/incumbent/compartments.py:116`` reads::

    linc_per_node=1,   # LINC tethers seeded per nucleus surface node (baseline coupling; ...)

so the LINC population size IS the nuclear mesh vertex count.  At the subdivision every resting run of
2026-08-20 used that is **642** — the same icosphere number ``bond.py`` records for ERM, from the same
mesh, wearing the same physiological label.  Nobody had to answer how many, for the same reason: one
declared connector stood in for a population.

Two further facts make this worse than ERM rather than equal to it:

  * **The incumbent's runtime LINC is not this connector.**  ``compartments.py:500`` pairs
    ``mesh.verts`` against ``pos_actin`` — nucleus to CORTEX ACTIN, inside a ``linc_reach_um = 3.0``
    search.  There is no intermediate filament on either end.  So ``if_nucleus_linc`` has no runtime
    object (``NOT_BOUND`` in ``connector_devicerun/native_record.json``) **and** no ancestor whose
    count could be inherited.  This is first implementation, not a port.
  * **The engine layer never answers it either.**  ``linc_connector._pack_linc_population`` sets
    ``capacity = len(specs)`` — the count is whatever the caller handed in.  The refusal has to happen
    here or it happens nowhere.

────────────────────────────────────────────────────────────────────────────────────────────────────
WHAT IS BLOCKED, AND WHAT IS NOT.

BLOCKED — see :data:`SPEC` for the exact wording Lead collects:

  * **the areal density** [1/µm² of nuclear envelope].  Neither the repository nor the Contract-Graph
    carries one: ``knowledge_claim`` and ``parameter`` both return **0 rows** for
    ``linc|nesprin|plectin|sun1|sun2|kash|lamin``, and of the **42** BM25 chunks that mention
    nesprin/plectin/LINC, **none** carries a per-area or per-nucleus number (queried 2026-08-20
    against ``outputs/tag_kb/kb.duckdb``; TAG-empty is not proof of absence, so this is reported as
    "not in this KB" rather than "does not exist").
  * **the nesprin-3 share.**  Even given a total LINC density, what separates this family's count from
    its two siblings' is the isoform split — nesprin-1/2 (actin) vs nesprin-3 (IF) vs nesprin-4/KASH5
    (MT).  Absent.  Without it, one density cannot be divided three ways except by invention.
  * **the stiffness.**  ``k_linc = 1.0e2`` pN/µm is labelled ``PI GAP`` in the incumbent and its own
    comment warns *"8 pN is a TENSION not a stiffness"* — i.e. the number behind it is not the quantity
    the slot wants.  No default here.
  * **the kinetic card.**  ``SourceGatedLincKinetics.propose_candidate`` refuses with *"no validated
    on/off-rate card exists"*; PI decision card C1 (2026-08-11) lists all three LINC edges under it.
  * **the scope, on top of all three.**  This cell is MCF7 epithelial, whose IF cage is keratin
    K8/K18; the nesprin-3/plectin literature is overwhelmingly vimentin in mesenchymal cells.  A number
    borrowed across that boundary is :attr:`~aleph.world.bond.SourceClass.UNRATIFIED_PROXY`, not
    ``SOURCED``, which is exactly the distinction that class was added for.  Whoever answers the
    density must also say which of the two it is.
  * **the per-filament escape is closed too.**  "One LINC per spoke" would make the count independent
    of the envelope mesh, but it makes it the CAGE count instead, and ``n_spokes`` is itself a
    provisional PI GAP (``build.intermediate_filament``, incumbent ``n_fil = 60``).  Trading a mesh
    number for a gap number is not an answer.

NOT BLOCKED, and implemented here:

  * **the pairing.**  Each spoke's inner bead is placed at ``R_nuc + gap`` by construction, on a
    Fibonacci direction, so its nearest envelope vertex is unambiguous and needs **no capture radius**.
    The incumbent's ``linc_reach_um = 3.0`` exists only because its undirected nucleus↔cortex search
    had to guess; that guess is not inherited.
  * **the rest length.**  Taken as the built separation, so the family is unstrained at ``t0``.  That
    is not a convenience: it is the ``r0_bind`` attach-unstrained rule already ratified in
    ``STATE.md`` (b) (``bc3fef26``), and it is what the incumbent does (``linc_rest = norm(...)``).
    ⚠ Whether a physiological LINC carries a resting PRESTRAIN is unknown and unsourced; recording
    zero is the ratified default, not a measurement.

WHAT IS DELIBERATELY ABSENT.  **No kinetic law.**  ``bond.py`` defers attach/detach, free lists and
snapshot twins to *"the first family that has kinetics"*, and this connector would be it — the engine
declares ``LINC_EVENT_CHANNELS = {bind, unbind}`` and ``commit_on_accept``.  Per the standing
instruction that is a PI escalation, not a thing to build: the shape of a free list is measurable once
one family detaches and guessable before.  This module builds a static population and says so.

Sanity Gate (recorded per CLAUDE.md, before first execution):
  * dimensional — positions and ``rest_um`` [µm]; ``stiffness_pn_per_um`` [pN/µm]; the density is
    [1/µm²] against a nuclear-envelope area support [µm²], which is what
    :attr:`~aleph.world.bond.BondCount.support_unit` checks by name.
  * boundary — zero spokes or zero envelope vertices refuses rather than returning an empty family: a
    cage with no nuclear anchor is recorded as "present" either way, which is the fact an empty claim
    destroys.  A resolved count of zero is a legitimate answer and is passed through to ``bond.py``.
    A count exceeding the candidates raises ``bond.py``'s message, which says to supply more candidates
    rather than lower the density.
  * conservation/invariant — **the candidate count is ``n_spokes x ranks`` and does not depend on the
    envelope vertex count.**  That is the anti-ERM invariant and :func:`_demo` asserts it directly by
    refining the envelope one level.  Every bond index lies in the arena's live NODE prefix, asserted
    by ``build_bond_family``.
  * CFL/precision — no integration; float64.  The stiffness is stored per bond because it enters the
    CFL bound through the row sum at its endpoints.
  * sign sense — the spoke's INNER bead is the nuclear end (``build.intermediate_filament``: "bead 0
    is the inner (nuclear) end"), so the caller supplies inner ends and a bond to an outer end would be
    the cage anchored to the cortex, which is a different connector.  ``rest_um >= 0`` by construction.
  * measurement protocol — host-side construction, pure bookkeeping, no device touched and nothing
    uploaded; positions arrive as host arrays read back BETWEEN accepted steps, never inside one.

engine units: length µm, stiffness pN/µm.  Runtime: pure host bookkeeping, CPU-importable.
"""

from __future__ import annotations

import numpy as np
import numpy.typing as npt

from aleph.world.bond import BondCount, BondFamily, build_bond_family
from aleph.world.families import ConnectorGapError, FamilySpec

__all__ = ["CONNECTOR", "CHEMISTRY_CARD", "SPEC", "candidate_pairs", "build_if_nucleus_linc"]

CONNECTOR = "if_nucleus_linc"

#: The card that IS this family's identity, taken verbatim from ``engine.linc_connector._EDGE_CHEMISTRY``
#: rather than re-spelled, for the reason ``SourceClass`` gives for sharing one vocabulary.
CHEMISTRY_CARD = "nesprin3_plectin_if_linc"

SPEC = FamilySpec(
    connector=CONNECTOR,
    populations=("intermediate_filament", "nucleus"),
    basis="areal",
    blocked_by=(
        "linc_areal_density_per_um2 — how many LINC complexes per µm² of nuclear envelope? "
        "0 rows in kb.duckdb knowledge_claim and parameter for linc|nesprin|plectin|sun|kash|lamin; "
        "0 of 42 matching BM25 chunks carry a per-area or per-nucleus number. The incumbent's answer "
        "is linc_per_node=1 x nucleus mesh vertices = 642 at subdiv 3, i.e. a mesh number.",
        "nesprin3_share_of_linc — what fraction of that density is nesprin-3/plectin (IF) rather than "
        "nesprin-1/2 (actin) or nesprin-4/KASH5 (MT)? Without it one density cannot be split across "
        "the three LINC families except by invention.",
        "scope_keratin_vs_vimentin — is a nesprin-3/plectin number measured on mesenchymal vimentin "
        "SOURCED for MCF7 keratin K8/K18, or UNRATIFIED_PROXY? Whoever supplies the density must say.",
        "k_if_nucleus_linc_pn_per_um — the incumbent's k_linc = 1.0e2 pN/µm is labelled PI GAP and its "
        "own comment says the source figure, 8 pN, is a TENSION not a stiffness.",
        "if_nucleus_linc_on_off_rates — PI card C1 2026-08-11: SourceGatedLincKinetics refuses with "
        "'no validated on/off-rate card exists'. ⚠ ESCALATION, not a datum: this would be the FIRST "
        "family with kinetics, which bond.py defers by design. PI decides the shape, this module does "
        "not build it.",
    ),
)


def candidate_pairs(
    *,
    spoke_inner_ids: npt.ArrayLike,
    spoke_inner_xyz: npt.ArrayLike,
    envelope_ids: npt.ArrayLike,
    envelope_xyz: npt.ArrayLike,
    ranks: int = 1,
) -> npt.NDArray[np.int64]:
    """Candidate ``(IF bead, envelope vertex)`` pairs, ordered so selection is never a random draw.

    For each spoke inner bead, the ``ranks`` nearest envelope vertices, emitted **rank-major**: every
    spoke's nearest vertex first, in spoke order, then every spoke's second-nearest, and so on.  A
    resolved count of ``n_spokes`` or fewer therefore takes at most one bond per spoke, and a larger
    count spreads over the cage before it doubles up on any one spoke.

    **This ordering is the anti-ERM property.**  The number of candidates is ``n_spokes x ranks`` and
    is independent of how finely the envelope is meshed, so refining the nuclear surface changes where
    the bonds land and never how many exist.

    Args:
        spoke_inner_ids: ``(S,)`` GLOBAL node index of each spoke's inner (nuclear) bead.
        spoke_inner_xyz: ``(S, 3)`` their positions [µm], read back between accepted steps.
        envelope_ids: ``(V,)`` GLOBAL node indices of the nuclear envelope vertices.
        envelope_xyz: ``(V, 3)`` their positions [µm].
        ranks: how many nearest envelope vertices to offer per spoke.  **1 offers exactly one candidate
            per spoke**; raise it when a density resolves to more bonds than there are spokes, which is
            what ``build_bond_family`` tells you to do rather than lowering the density.

    Returns:
        ``(S * ranks, 2)`` int64 pairs, ``[:, 0]`` the IF bead and ``[:, 1]`` the envelope vertex.

    Raises:
        ValueError: on an empty spoke or envelope list, mismatched id/position lengths, a non-positive
            ``ranks``, or ``ranks`` exceeding the number of envelope vertices.
    """
    s_id = np.asarray(spoke_inner_ids, np.int64).reshape(-1)
    v_id = np.asarray(envelope_ids, np.int64).reshape(-1)
    s_xyz = np.asarray(spoke_inner_xyz, np.float64).reshape(-1, 3)
    v_xyz = np.asarray(envelope_xyz, np.float64).reshape(-1, 3)

    if s_id.size == 0 or v_id.size == 0:
        raise ValueError(
            f"{CONNECTOR}: {s_id.size} spoke inner ends and {v_id.size} envelope vertices. A cage with "
            "no nuclear anchor is not a weak load path, it is no load path — and the composed cell "
            "would record the nucleus as coupled either way."
        )
    if s_id.size != s_xyz.shape[0] or v_id.size != v_xyz.shape[0]:
        raise ValueError(
            f"{CONNECTOR}: {s_id.size} spoke ids against {s_xyz.shape[0]} positions, and {v_id.size} "
            f"envelope ids against {v_xyz.shape[0]}. An id and its position must be the same index."
        )
    if ranks < 1:
        raise ValueError(f"{CONNECTOR}: ranks must be at least 1; got {ranks}")
    if ranks > v_id.size:
        raise ValueError(
            f"{CONNECTOR}: ranks={ranks} exceeds the {v_id.size} envelope vertices available"
        )

    # ponytail: brute-force (S x V) distance scan. S is the cage count (incumbent 60) and V the
    # envelope vertex count (642 at subdiv 3), so this is ~4e4 rows. Swap in a KD-tree if either grows
    # by two orders of magnitude; nothing about the ordering above changes if it does.
    d2 = ((s_xyz[:, None, :] - v_xyz[None, :, :]) ** 2).sum(axis=2)
    order = np.argsort(d2, axis=1, kind="stable")[:, :ranks]  # stable: ties break by envelope index
    pairs = np.empty((s_id.size * ranks, 2), np.int64)
    for r in range(ranks):
        pairs[r * s_id.size:(r + 1) * s_id.size, 0] = s_id
        pairs[r * s_id.size:(r + 1) * s_id.size, 1] = v_id[order[:, r]]
    return pairs


def build_if_nucleus_linc(
    arena,
    *,
    count: BondCount | None = None,
    support: float | None = None,
    stiffness_pn_per_um: float | None = None,
    spoke_inner_ids: npt.ArrayLike,
    spoke_inner_xyz: npt.ArrayLike,
    envelope_ids: npt.ArrayLike,
    envelope_xyz: npt.ArrayLike,
    ranks: int = 1,
    name: str = CONNECTOR,
) -> BondFamily:
    """Build the nesprin-3/plectin LINC family, or refuse by naming what is unanswered.

    The geometry is implemented and the numbers are not, which is the honest shape: this refuses today
    and builds the moment the PI supplies a density and a stiffness, with no further code.

    Args:
        arena: the world to claim the BOND range from.
        count: how many LINC bonds, against what, for which cell, on whose authority.  **No default** —
            see :data:`SPEC` for why there cannot be one.
        support: the MEASURED nuclear-envelope area [µm²] the density resolves against — the area of
            the mesh that exists, ``ClosedSurface.area0_total_um2``, never the analytic sphere.  That is
            what keeps the count invariant when the envelope resolution changes.  Pass ``1.0`` only for
            an ``explicit`` count.
        stiffness_pn_per_um: per-bond stiffness [pN/µm].  **No default** — see :data:`SPEC`.
        spoke_inner_ids / spoke_inner_xyz: the IF cage's inner (nuclear) beads and their positions.
        envelope_ids / envelope_xyz: the nuclear envelope vertices and their positions.
        ranks: passed to :func:`candidate_pairs`.
        name: the family name claimed in the arena.

    Returns:
        The built :class:`~aleph.world.bond.BondFamily`, unstrained at ``t0``.

    Raises:
        ConnectorGapError: if ``count``, ``support`` or ``stiffness_pn_per_um`` is absent.  **This is
            the expected outcome as of 2026-08-20** and is the module's product, not its failure.
        ValueError: from :func:`candidate_pairs` or ``build_bond_family``.
    """
    missing = []
    if count is None:
        missing.append(SPEC.blocked_by[0])
        missing.append(SPEC.blocked_by[1])
        missing.append(SPEC.blocked_by[2])
    if support is None and count is not None:
        missing.append(
            "support — the MEASURED nuclear-envelope area [µm²] the density resolves against. It is "
            "not optional and it is not the analytic sphere: measuring against the mesh that exists is "
            "what stops the count from moving with the resolution."
        )
    if stiffness_pn_per_um is None:
        missing.append(SPEC.blocked_by[3])
    if missing:
        raise ConnectorGapError(
            CONNECTOR, missing,
            note="No kinetic law is built either — see SPEC.blocked_by[-1]: this would be the first "
                 "family with kinetics and bond.py defers that to a PI decision. The pairing, the rest "
                 "length and the arena claim are implemented and need nothing further.",
        )

    pairs = candidate_pairs(
        spoke_inner_ids=spoke_inner_ids, spoke_inner_xyz=spoke_inner_xyz,
        envelope_ids=envelope_ids, envelope_xyz=envelope_xyz, ranks=ranks)

    # Attach-unstrained: the rest length IS the built separation (STATE.md (b), `r0_bind` bc3fef26).
    s_xyz = np.asarray(spoke_inner_xyz, np.float64).reshape(-1, 3)
    v_xyz = np.asarray(envelope_xyz, np.float64).reshape(-1, 3)
    s_pos = {int(i): k for k, i in enumerate(np.asarray(spoke_inner_ids, np.int64).reshape(-1))}
    v_pos = {int(i): k for k, i in enumerate(np.asarray(envelope_ids, np.int64).reshape(-1))}
    rest = np.linalg.norm(
        np.stack([s_xyz[s_pos[int(a)]] - v_xyz[v_pos[int(b)]] for a, b in pairs]), axis=1)

    return build_bond_family(
        arena, name, chemistry_card=CHEMISTRY_CARD, count=count, support=float(support),
        pairs=pairs, rest_um=rest, stiffness_pn_per_um=float(stiffness_pn_per_um))


def _demo() -> None:
    """Self-check: the refusal is the product, and the count does not move with the envelope mesh."""
    from aleph.world.arena import Kind, WorldArena
    from aleph.world.surface import build_surface, icosphere

    assert not SPEC.buildable and len(SPEC.blocked_by) == 5, SPEC
    assert CONNECTOR not in __import__(
        "aleph.world.families", fromlist=["DUPLICATE_OF"]).DUPLICATE_OF, \
        "the three LINC edges carry three chemistry cards; see this module's docstring"

    def _world(subdiv: int, n_spokes: int = 6):
        """One arena holding a nuclear envelope and a radial IF cage, host-side, no device."""
        verts, faces = icosphere(subdiv)
        arena = WorldArena(capacity={Kind.NODE: 20_000, Kind.FACE: 40_000, Kind.ANGLE4: 60_000, Kind.BOND: 5_000})
        env = build_surface(arena, "nucleus", vertices=verts, faces=faces, radius_um=5.0)
        # The cage's inner beads: Fibonacci directions at R_nuc + gap, the placement
        # `build.intermediate_filament` makes. Claimed as their own population so the derived pair reads.
        k = np.arange(n_spokes) + 0.5
        cosp = 1.0 - 2.0 * k / n_spokes
        sinp = np.sqrt(np.maximum(0.0, 1.0 - cosp**2))
        az = 10.166407384630519 * k
        inner = 5.15 * np.stack([sinp * np.cos(az), sinp * np.sin(az), cosp], axis=1)
        claim = arena.claim("intermediate_filament", Kind.NODE, n_spokes)
        ids = claim.lo + np.arange(n_spokes, dtype=np.int64)
        return arena, ids, inner, env

    # subdiv 3 = 642 envelope vertices, the resolution every resting run of 2026-08-20 used.
    arena, s_ids, s_xyz, env = _world(subdiv=3)
    kw = dict(spoke_inner_ids=s_ids, spoke_inner_xyz=s_xyz,
              envelope_ids=env.nodes.lo + np.arange(env.n_vertices), envelope_xyz=env.position)

    # 1. With no density and no stiffness it refuses, and the refusal NAMES the questions.
    try:
        build_if_nucleus_linc(arena, **kw)
    except ConnectorGapError as exc:
        assert exc.connector == CONNECTOR
        assert len(exc.missing) == 4, exc.missing          # 3 count questions + the stiffness
        assert "linc_areal_density_per_um2" in str(exc) and "TENSION not a stiffness" in str(exc)
        assert "first family with kinetics" in str(exc)
    else:  # pragma: no cover
        raise AssertionError("a family with no density and no stiffness must refuse")

    # 2. THE ANTI-ERM INVARIANT: refining the envelope changes where bonds land, never how many exist.
    coarse = candidate_pairs(**kw)
    arena2, s2, x2, env2 = _world(subdiv=4)
    assert env2.n_vertices > env.n_vertices, (env.n_vertices, env2.n_vertices)
    fine = candidate_pairs(spoke_inner_ids=s2, spoke_inner_xyz=x2,
                           envelope_ids=env2.nodes.lo + np.arange(env2.n_vertices),
                           envelope_xyz=env2.position)
    assert coarse.shape == fine.shape == (6, 2), (coarse.shape, fine.shape)
    assert not np.array_equal(coarse[:, 1], fine[:, 1] - (env2.nodes.lo - env.nodes.lo)), \
        "the mesh refined, so at least one anchor vertex must have moved"

    # 3. Each spoke gets its OWN nearest vertex, and it is on the spoke's own side of the nucleus.
    assert len(set(coarse[:, 0].tolist())) == 6
    for a, b in coarse:
        u = s_xyz[int(a) - int(s_ids[0])]
        v = env.position[int(b) - env.nodes.lo]
        assert float(u @ v) > 0.0, "the nearest envelope vertex must not be on the far side"

    # 4. Given a PI answer it builds, the pair is DERIVED, and the family is unstrained at t0.
    pi_count = BondCount(
        basis="explicit", value=6.0,
        scope="_demo only — NOT a physiological count; see SPEC.blocked_by",
        source_class="PI_GAP", provenance="a stand-in so the build path is exercised; the real count "
                                          "is the first entry of SPEC.blocked_by")
    fam = build_if_nucleus_linc(arena, count=pi_count, support=1.0, stiffness_pn_per_um=1.0e2, **kw)
    assert fam.n_bonds == 6 and fam.chemistry_card == CHEMISTRY_CARD
    assert fam.component_pairs(arena) == {("intermediate_filament", "nucleus"): 6}
    # Unstrained at t0 means the rest length IS the built separation, to the digit — not a nominal
    # gap. It exceeds gap_um = 0.15 because the nearest MESH vertex is not on the spoke's radial line,
    # which is a fact about the envelope's resolution and is why rest is measured, never assumed.
    for k, (a, b) in enumerate(zip(fam.node_i.tolist(), fam.node_j.tolist())):
        u = s_xyz[a - int(s_ids[0])]
        v = env.position[b - env.nodes.lo]
        assert abs(float(np.linalg.norm(u - v)) - float(fam.rest_um[k])) < 1e-12
    assert float(fam.rest_um.min()) >= 0.15 - 1e-9, fam.rest_um
    assert fam.provenance_row()["scope"] and fam.provenance_row()["support_unit"] == "none"

    # 5. A count larger than the candidates says to supply more candidates, not to lower the density.
    over = BondCount(basis="explicit", value=12.0, scope="_demo", source_class="PI_GAP",
                     provenance="deliberately over the 6 candidates rank=1 offers")
    try:
        build_if_nucleus_linc(arena, count=over, support=1.0, stiffness_pn_per_um=1.0e2, **kw)
    except ValueError as exc:
        assert "Supply more candidates" in str(exc), exc
    else:  # pragma: no cover
        raise AssertionError("an over-subscribed count must refuse")
    wide = candidate_pairs(**{**kw, "ranks": 2})
    assert wide.shape == (12, 2) and np.array_equal(wide[:6], coarse), "rank-major, nearest first"

    # 6. An empty cage is refused rather than recorded as a present-but-empty coupling.
    try:
        candidate_pairs(**{**kw, "spoke_inner_ids": np.zeros(0, np.int64),
                           "spoke_inner_xyz": np.zeros((0, 3))})
    except ValueError as exc:
        assert "no load path" in str(exc)
    else:  # pragma: no cover
        raise AssertionError("an empty cage must refuse")

    print(f"if_nucleus_linc self-check OK — BLOCKED on {len(SPEC.blocked_by)}: "
          + "; ".join(q.split(" — ")[0] for q in SPEC.blocked_by))


if __name__ == "__main__":
    _demo()
