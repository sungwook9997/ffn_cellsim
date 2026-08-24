r"""``filopodium_cortex_root`` — the bundle root, and the number nobody has ever written down.

WHAT THE CONNECTOR IS.  A filopodium is a parallel formin bundle standing off the cell surface
(:mod:`aleph.world.build.filopodium`); bead 0 of every filament in that bundle is its ROOT, sitting at
``root_R_um`` on the polar cap.  This edge is what holds that root end in the cortical actin network.
The incumbent declares it ``TRANSIENT_ACTIN``, ``filopodium -> cortex``, chemistry card
``formin_fascin_root_coupling``, endpoint roles *"filopodium bundle root"* / *"cortical actin material
point"* (``engine/contracts.py:843``), and its cortex endpoint addresses the NODE domain, not the
segment domain (``engine/cortex_population.py:352``) — a root anchors to a material point, it does not
walk a polar coordinate the way a motor head does.

**It has no runtime object.**  It is one of the 28 ``NOT_BOUND`` entries in
``outputs/ac/connector_devicerun/native_record.json``, so nothing here is a port.

WHY IT IS NOT A DUPLICATE, ARGUED FROM PHYSICS.  Three declared edges bind the SAME runtime,
``FilamentCrosslinkConnector`` — ``sf_cortex_transient``, ``lamellipodium_cortex_seam`` and this one —
because a state-gated Hookean central force between two interpolated points is one force law and
writing it three times would give the adjoint three places to drift.  **Sharing a force law is not
being the same family.**  The first two both declare card ``transient_actin_crosslink``: same
crosslinker, same rate law, two declarations standing for one population that the arena addresses by ID
range.  *Those two are duplicates of each other.*  This edge declares a different card, and the card is
what decides when a bond EXISTS — the property a :class:`~aleph.world.bond.BondFamily` is the carrier
for.  A bundle root held by fascin at the base of a formin bundle has a different hand, a different
on-rate and a different failure mode from a transient α-actinin-class crosslink between two networks
that happen to touch.  Folding them would put two rate laws behind one, which is the mirror image of
the mistake ``DUPLICATE_OF`` exists to prevent, so this module builds its own family and adds no entry.

That verdict was reached here independently and then met
:data:`~aleph.world.families.DUPLICATE_CRITERION`, which three other sessions converged on the same
day — *a family's identity is its chemistry card; a shared runtime class or force law is not a
duplicate*.  D13 uses THIS edge as the control that makes its own fold trustworthy: the two edges it
folds sit in the same ``ConnectorFamily.TRANSIENT_ACTIN`` as this one, so folding by card is
demonstrably not folding by name.

⚠ THE COUNT.  :class:`~aleph.world.bond.BondCount` asks four questions.  Two are answered here and two
are not, and the two that are not are not answerable from anything in this repository.

* *Against what* — **answered: ``per_filament``, support = the MEASURED ``n_strands`` of the built
  filopodium population.**  Every filament in a bundle has exactly one root end, so the root population
  is per-filament by construction, and the support is counted off what actually stands rather than off
  an analytic bundle.  This is the axis the ERM defect turned on and it is worth stating why it holds
  here: ``n_strands = n_fingers x filaments_per_bundle`` and **neither factor is a discretisation**.
  Refining ``seg_um`` doubles ``nodes_per_strand`` and leaves this count untouched — asserted in
  :func:`_demo`, because that assertion is the whole difference between a density and a mesh artefact.
* *How many* — **UNANSWERED.**  How many cortical attachments ONE root end makes has no value anywhere.
  ``architecture_spec.FILOPODIUM`` gives filaments per bundle (20, band 10–30) and the intra-bundle
  crosslink density (3.0/filament, with FILAMIN standing in because *"fascin kinetics are PI-gated"*)
  and says nothing whatever about the root.  Writing 1 would be a build convenience wearing a
  physiological label: the two root architectures in the literature — convergent elongation, where the
  root is embedded over a length in the network and crosslinked along it, versus formin tip-nucleation,
  where a single pointed end sits in the cortex — give different numbers AND different embedded
  lengths, and which one this cell uses is not recorded either.
* *For which cell* — **UNANSWERED.**  No line, state, assay or temperature.
* *On whose authority* — **UNANSWERED.**  ``formin_fascin_root_coupling`` appears exactly three times in
  the repository, all three of them declarations (``contracts.py``, ``protrusion.py``, and a test
  asserting the string), and **carries no parameter anywhere**: no rest length, no stiffness, no rate.
  The incumbent never had to answer this because its builder takes the joint specs from its caller
  (``filament_crosslink.build_filament_crosslink_connector``), so the count lived outside the runtime.
  There is no ``ValidationGate`` or ``ModelContract`` row, and those are PI-authored only.

So :func:`build_filopodium_cortex_root` raises :class:`~aleph.world.families.ConnectorGapError` naming
those questions, and :data:`SPEC` carries them in ``blocked_by`` for Lead's decision queue.  It builds
the moment the PI supplies the three missing answers, and not before.  The pairing arithmetic is here
and is exact, so what is blocked is the physiology and only the physiology.

⚠ ONE DEPENDENCY THAT IS NOT MINE.  The number of filopodia per cell is *absent from the contract
graph* — ``build/filopodium.py`` says so — so the support this family resolves against does not exist
yet either.  That gap belongs to the filopodium population, not to this edge, and is named here only so
the queue does not double-count it.

WHAT IS DELIBERATELY ABSENT.  No kinetic law.  ``bond.py`` defers attach/detach to *"the first family
that has kinetics"*, and this edge declares ``kinetic=True`` in the incumbent — so if it turns out to
be that first family, that is a PI item and not something to build here.

Sanity Gate (recorded per CLAUDE.md, before first execution):
  * dimensional — ``rest_um`` [µm], ``stiffness_pn_per_um`` [pN/µm]; the basis is ``per_filament`` so
    the support is a dimensionless filament count and :meth:`BondCount.resolve` multiplies a pure
    number by a pure number.  ``BondCount`` rejects a support carrying any other basis by name.
  * boundary — a count resolving to zero is a legal answer and returns an empty family (``bond.py``);
    a build with no answers refuses with :class:`ConnectorGapError` rather than a placeholder; more
    bonds than candidate pairs refuses in ``build_bond_family`` with both numbers.
  * conservation/invariant — every root index is derived from the population's own NODE claim, so it
    lies inside it by construction; asserted against the claim bounds before any bond is made, and
    ``build_bond_family`` re-checks against the arena's live NODE prefix.  Partner-major pair ordering
    keeps the attachments BALANCED across roots: with ``build_bond_family`` taking the first ``N``
    candidates in order, root-major ordering would give every bond to the first few roots.
  * CFL/precision — no integration; float64.  Stiffness is stored per bond because it enters the CFL
    bound through the row sum at its endpoints.
  * sign sense — bead 0 is the root and the last bead is the tip (``build/filopodium.py`` sets
    ``polarity=+1`` for exactly this reason); rooting the tip instead would invert what the bundle
    pushes against.  ``rest_um`` is a length and is non-negative.
  * measurement protocol — host-side index arithmetic only.  No device is touched, nothing is read
    back, and no position is inspected: which cortex node each root pairs with is the caller's
    deterministic choice, per ``build_bond_family``'s contract on ``pairs``.

engine units: length µm, stiffness pN/µm.  Runtime: pure host bookkeeping, CPU-importable.
"""

from __future__ import annotations

import numpy as np
import numpy.typing as npt

from aleph.world.arena import WorldArena
from aleph.world.bond import BondCount, BondFamily, build_bond_family
from aleph.world.families import ConnectorGapError, FamilySpec

__all__ = ["SPEC", "CONNECTOR", "CHEMISTRY_CARD", "root_nodes", "candidate_pairs",
           "build_filopodium_cortex_root"]

CONNECTOR = "filopodium_cortex_root"

#: The card the incumbent declares, kept as the family identity rather than the component pair. Two
#: edges sharing a card are one family; this one is shared with nothing.
CHEMISTRY_CARD = "formin_fascin_root_coupling"

#: Bead 0 of every bundle filament. Not a convention: ``build/filopodium.py`` lays the bundle outward
#: from the surface with uniform polarity, so bead 0 IS the root and the last bead is the barbed tip.
ROOT_BEAD = 0

SPEC = FamilySpec(
    connector=CONNECTOR,
    populations=("filopodium", "cortex"),
    basis="per_filament",
    blocked_by=(
        "how many cortical attachments ONE root end makes — no value exists; architecture_spec."
        "FILOPODIUM gives filaments per bundle and intra-bundle crosslinks and nothing for the root, "
        "and 1 would be a build convenience, not a measurement",
        "which root architecture this cell has — convergent elongation (root embedded over a length, "
        "crosslinked along it) vs formin tip-nucleation (one pointed-end anchor) give different "
        "numbers and different embedded lengths; unrecorded either way",
        "scope: no cell line, state, assay or temperature is attached to the root coupling anywhere",
        "the rest length and stiffness of the coupling — `formin_fascin_root_coupling` appears 3x in "
        "the repository, all declarations, and carries no parameter of any kind",
        "authority: no ValidationGate/ModelContract row, and fascin kinetics are already PI-gated in "
        "architecture_spec.FILOPODIUM ('fascin = PI-gated', SE section 5)",
    ),
    source="none found. Searched: laws/architecture_spec.FILOPODIUM, engine/{contracts,protrusion,"
           "filament_crosslink,cortex_population}.py, outputs/ac/topology/graph.json and the "
           "connector_devicerun census. The card is a name with no numbers behind it. INHERITED and "
           "not counted here: filopodia per cell is absent from the contract graph (build/filopodium.py"
           "), so this family's support does not exist yet either.",
)


def root_nodes(inventory: dict[str, object]) -> npt.NDArray[np.int64]:
    """The GLOBAL node index of every bundle filament's root bead.

    Args:
        inventory: the census returned by :func:`aleph.world.build.filopodium.build_filopodia`. A
            planned-but-unbuilt inventory has no ``claims`` and is refused — a root index into a
            population that was never claimed addresses somebody else's nodes.

    Returns:
        ``(n_strands,)`` global node indices, in strand order.

    Raises:
        KeyError: if the inventory was never built.
        AssertionError: if the derived indices leave the population's own NODE claim.
    """
    lo, hi = inventory["claims"]["node"]  # type: ignore[index]
    n_strands = int(inventory["n_strands"])  # type: ignore[arg-type]
    n_per = int(inventory["nodes_per_strand"])  # type: ignore[arg-type]
    roots = int(lo) + np.arange(n_strands, dtype=np.int64) * n_per + ROOT_BEAD
    assert n_strands * n_per == int(hi) - int(lo), (
        f"{CONNECTOR}: the inventory says {n_strands} strands x {n_per} beads but its NODE claim spans "
        f"{int(hi) - int(lo)}. The root index is arithmetic on that claim, so a mismatch means the "
        "roots would land in another population."
    )
    return roots


def candidate_pairs(
    roots: npt.ArrayLike, cortex_partners: npt.ArrayLike
) -> npt.NDArray[np.int64]:
    """Interleave roots with their candidate cortex partners, PARTNER-major.

    Args:
        roots: ``(R,)`` global root node indices.
        cortex_partners: ``(R, m)`` global cortex node indices — the ``m`` candidate material points
            offered for each root, in the caller's own preference order. Choosing them is the caller's
            deterministic job (``build_bond_family``'s contract on ``pairs``); this module owns the
            index arithmetic and touches no position.

    Returns:
        ``(R*m, 2)`` pairs ordered ``[all roots' 1st partner, all roots' 2nd, ...]``. That order is the
        point: ``build_bond_family`` takes the first ``N`` candidates, so root-major ordering would give
        every bond to the first roots and leave the rest of the bundle unattached.

    Raises:
        ValueError: if the partner table's first axis does not match ``roots``.
    """
    r = np.asarray(roots, np.int64).reshape(-1)
    p = np.asarray(cortex_partners, np.int64)
    if p.ndim == 1:
        p = p.reshape(-1, 1)
    if p.shape[0] != r.size:
        raise ValueError(
            f"{CONNECTOR}: {r.size} roots but a partner table of shape {p.shape}; one row per root"
        )
    return np.stack([np.tile(r, p.shape[1]), p.T.reshape(-1)], axis=1)


def build_filopodium_cortex_root(
    arena: WorldArena,
    inventory: dict[str, object],
    *,
    cortex_partners: npt.ArrayLike,
    count: BondCount | None = None,
    rest_um: float | None = None,
    stiffness_pn_per_um: float | None = None,
) -> BondFamily:
    """Build the root family, or refuse by name.

    Args:
        arena: the world holding both populations.
        inventory: the built filopodium census (see :func:`root_nodes`).
        cortex_partners: ``(R, m)`` candidate cortex node indices per root; see
            :func:`candidate_pairs`.
        count: how many attachments per root end, ``per_filament``, with its scope and authority.
            **No value exists** — see the module docstring.
        rest_um: the coupling's rest length [µm]. **No value exists.**
        stiffness_pn_per_um: the coupling's stiffness [pN/µm]. **No value exists.**

    Returns:
        The built :class:`~aleph.world.bond.BondFamily`.

    Raises:
        ConnectorGapError: whenever any of the three physiological answers is missing. This is the
            module's product today, not its failure.
        ValueError: if ``count`` is given on a basis other than ``per_filament`` — the support this
            family resolves against is a filament count, and a density in any other basis multiplied
            by it is a unit error.
    """
    missing = [q for q, given in (
        ("how many attachments per root end (`count`)", count),
        ("the coupling's rest length (`rest_um`)", rest_um),
        ("the coupling's stiffness (`stiffness_pn_per_um`)", stiffness_pn_per_um),
    ) if given is None]
    if missing:
        raise ConnectorGapError(
            CONNECTOR, missing + list(SPEC.blocked_by),
            "Not buildable from anything in this repository; PI-authored answers only. The pairing "
            "arithmetic is implemented and exact — only the physiology is missing.",
        )
    assert count is not None and rest_um is not None and stiffness_pn_per_um is not None
    if count.basis != "per_filament":
        raise ValueError(
            f"{CONNECTOR}: the support is a filament count, so the basis must be per_filament; got "
            f"{count.basis!r}. A count on another basis multiplied by a filament count is a unit error."
        )

    roots = root_nodes(inventory)
    return build_bond_family(
        arena, CONNECTOR, chemistry_card=CHEMISTRY_CARD, count=count, support=float(roots.size),
        pairs=candidate_pairs(roots, cortex_partners), rest_um=rest_um,
        stiffness_pn_per_um=stiffness_pn_per_um,
    )


def _demo() -> None:
    """Self-check: the refusal is the product, the count is seg-invariant, and the pairing is balanced."""
    from aleph.world.arena import Kind
    from aleph.world.build.filopodium import plan_filopodia

    per_cell = BondCount(
        basis="explicit", value=4.0, scope="demo only — NOT a cell", source_class="PI_GAP",
        provenance="filopodia per cell is absent from the contract graph; 4 here is a fixture size")
    per_bundle = BondCount(
        basis="explicit", value=5.0, scope="demo only — NOT a cell", source_class="PI_GAP",
        provenance="architecture_spec gives 20 in a band of 10-30; 5 here is a fixture size")
    geom = {"contour_um": 0.4, "spacing_um": 0.008, "cap_deg": 60.0, "root_R_um": 7.5}

    # THE anti-ERM assertion. Refining the discretisation must not move the population this family's
    # count is taken against — if it did, the count would be a mesh number wearing a physiological
    # label, which is the defect `bond.py` is shaped around.
    coarse = plan_filopodia(count=per_cell, bundle_count=per_bundle, seg_um=0.1, **geom)
    fine = plan_filopodia(count=per_cell, bundle_count=per_bundle, seg_um=0.05, **geom)
    assert coarse["n_strands"] == fine["n_strands"] == 20
    assert fine["nodes_per_strand"] == 2 * coarse["nodes_per_strand"] - 1, "the mesh DID refine"

    # Stand the two populations in an arena. Only the ID ranges matter here; no geometry is read.
    arena = WorldArena(capacity={Kind.NODE: 400, Kind.BOND: 400})
    n_per = int(coarse["nodes_per_strand"])
    n_strands = int(coarse["n_strands"])
    filo = arena.claim("filopodium", Kind.NODE, n_strands * n_per)
    cortex = arena.claim("cortex", Kind.NODE, 60)
    built = dict(coarse) | {"claims": {"node": (filo.lo, filo.hi)}}

    roots = root_nodes(built)
    assert roots.size == n_strands
    assert roots[0] == filo.lo and int(roots.max()) + n_per - 1 == filo.hi - 1
    assert np.all(np.diff(roots) == n_per), "one root per strand, at bead 0"

    # Three cortex partners per root, and the pair order is partner-major so a truncated take stays
    # balanced across roots rather than piling onto the first ones.
    partners = cortex.lo + (np.arange(n_strands * 3, dtype=np.int64) % 60).reshape(n_strands, 3)
    pairs = candidate_pairs(roots, partners)
    assert pairs.shape == (n_strands * 3, 2)
    assert np.array_equal(pairs[:n_strands, 0], roots), "the first R pairs are one per root"
    try:
        candidate_pairs(roots, partners[:-1])
    except ValueError as exc:
        assert "one row per root" in str(exc)
    else:  # pragma: no cover
        raise AssertionError("a partner table that does not match the roots must refuse")

    # The refusal, which is what this module actually delivers today.
    for kwargs in ({}, {"rest_um": 0.03}, {"rest_um": 0.03, "stiffness_pn_per_um": 100.0}):
        try:
            build_filopodium_cortex_root(arena, built, cortex_partners=partners, **kwargs)
        except ConnectorGapError as exc:
            assert exc.connector == CONNECTOR
            assert len(exc.missing) == 3 - len(kwargs) + len(SPEC.blocked_by)
        else:  # pragma: no cover
            raise AssertionError("a family with no physiological answer must refuse to build")
    assert arena.n_live(Kind.BOND) == 0, "a refused family claims nothing"
    assert not SPEC.buildable and len(SPEC.blocked_by) == 5

    # And it DOES build the moment the answers arrive — with a count that is nobody's ruling, so this
    # stays a fixture and may not be quoted as an inventory.
    hypothetical = BondCount(
        basis="per_filament", value=2.0, scope="DEMO FIXTURE — no cell, no authority",
        source_class="PI_GAP", provenance="not a value: exercises the path the PI's answer will take")
    fam = build_filopodium_cortex_root(
        arena, built, cortex_partners=partners, count=hypothetical, rest_um=0.03,
        stiffness_pn_per_um=100.0)
    assert fam.n_bonds == 2 * n_strands
    assert np.array_equal(np.sort(np.unique(fam.node_i, return_counts=True)[1]),
                          np.full(n_strands, 2)), "every root got exactly 2, none got 4"
    # The component pair is DERIVED — this module never declared filopodium->cortex to the bond.
    assert fam.component_pairs(arena) == {("cortex", "filopodium"): 2 * n_strands}
    assert fam.provenance_row()["support_unit"] == "filaments"

    try:
        build_filopodium_cortex_root(
            arena, built, cortex_partners=partners, rest_um=0.03, stiffness_pn_per_um=100.0,
            count=BondCount(basis="areal", value=2.0, scope="d", source_class="PI_GAP", provenance="d"))
    except ValueError as exc:
        assert "must be per_filament" in str(exc)
    else:  # pragma: no cover
        raise AssertionError("an areal density against a filament support must refuse")

    print(f"{CONNECTOR} self-check OK — BLOCKED on {len(SPEC.blocked_by)} questions; pairing exact for "
          f"{n_strands} roots, seg-invariant")


if __name__ == "__main__":
    _demo()
