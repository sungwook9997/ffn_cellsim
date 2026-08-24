r"""``mt_cortex_capture`` — cortical dynein capturing a microtubule. **It refuses, and that is the point.**

WHAT THE CONNECTOR IS.  ``engine/contracts.py:751`` declares it ``ConnectorFamily.MOTOR``, microtubule
<-> cortex, bidirectional with adjoint transfer, ``kinetics=True``, ``chemistry_card=
"microtubule_cortical_dynein"``; endpoint role A is *"dynamic microtubule plus-end or lattice motor
site"* and role B *"cortical dynein capture site"*.  ``engine/load_path.py:128`` gives it its own
``JointKind.CORTICAL_DYNEIN_CAPTURE = 7``: *"a bound cortical-dynein crossbridge on a microtubule
lattice site"*.  In the committed device-run census
(``outputs/ac/connector_devicerun/native_record.json``) it is one of the 36 ``unbound_connectors`` — no
runtime object at all.  So this module is a FIRST IMPLEMENTATION, not a port, and it is written against
a declaration rather than against running code.

NOT A DUPLICATE, BY THE RATIFIED CRITERION.  ``families/__init__.DUPLICATE_CRITERION`` says a family's
identity is its CHEMISTRY CARD, and that a shared runtime class or a shared force law is not a
duplicate.  ``microtubule_cortical_dynein`` occurs exactly ONCE across the 38 declared connectors — the
adjacent MT edges carry ``microtubule_motor_linc`` (D3) and ``spectraplakin_mt_actin`` (D8), and the
four NMII motor edges all carry ``nmii_head_actin_hill_bell``.  The temptation here is the force law:
``connector_joints.py`` records that *"motor, actin_anchor and fa_clutch are three families and ONE
mechanics"* — a state-gated bilateral Hookean central force that
``load_path._segment_joint_force_kernel`` already is.  Under the criterion that shared law is a shared
CONTAINER, and folding on it would hide a processive minus-end-directed motor inside a passive
crosslinker.  Two further contract differences confirm the split rather than merely permitting it: the
three LINC edges declare ``generation_required`` / ``remap_on_accept`` / ``blocks_sleep_refine`` and
capture declares none of them, and capture is the only MT edge with its own ``JointKind``.
**Verdict: its own family. Nothing is added to ``DUPLICATE_OF``.**

⚠ WHY IT CANNOT BE BUILT TODAY — the four questions, each searched before it was declared unanswerable.

*How many.*  Nothing anywhere states a number of cortical dynein capture sites, or of bound MT-cortex
crossbridges, per cell.  The Contract-Graph ``parameter`` table returns **0 rows** for
``dynein``/``capture``/``mt_cortex``/``microtubule``.  ``paper_chunks`` (12,810 chunks) has 75 chunks
mentioning dynein and **3** matching *"cortical dynein"*, all three from one axonal MT-bundle
simulation paper about kinesin-bundled extensile bundles in a neurite — a different compartment, a
different cell, and no areal density in any of them.  ``PI_GAP_EVIDENCE_CARDS_2026-07-25.md`` has no
card for this connector at all: its E-series cards (E3/E4/E5) are ERM.  The one adjacent literature
lead in ``references/clipped/`` — cortical dynein tracking during asymmetric division — is a *C.
elegans* one-cell embryo at mitosis, which is a scope refusal (``SourceClass.UNRATIFIED_PROXY``), not a
source.

*Against what* — and this one is a finding rather than a gap.  A "cortical dynein site" reads as an
AREAL density on the cortex, and **the cortex population has no measurable area**: ``world/`` builds it
as 4,197,654 NODE + SEGMENT + ANGLE3 and claims no FACE (``outputs/ac/world_phase1/phase1_native.json``).
The faces belong to ``membrane``.  So an areal basis would have to resolve its support against a THIRD
population's mesh — a membrane triangle area standing in for the support of a bond that does not touch
the membrane, which is the ERM defect's shape exactly.  Of the four bases, only ``per_filament`` has a
support this arena can measure today (the microtubule STRAND claim), and ``SPEC.basis`` records that as
a PROPOSAL, not an answer.

*Which cell.*  Unanswerable while the support is: the microtubule count itself is a PI-GAP.
``PER_CELL_STRUCTURE_COUNTS_2026-08-20.md`` §3 records *"no per-cell count anywhere in the corpus"* and
proposes 250-600/epithelial cell as a declared TEST AXIS, with PubMed evidence that an adenocarcinoma
line differs in nucleating capacity.  ``world/build/microtubule.py`` ships 500 strands from a
``BondCount`` whose ``source_class`` is ``PI_GAP``.  A per-filament capture count resolved against 500
would inherit an unratified support and read as a measurement.

*On whose authority.*  Nobody's.  There is no ratified source, no KB row and no PI card.

⚠ A SECOND ANSWER IS MISSING THAT IS NOT A NUMBER: **which MT sites are candidates at all.**  The
contract says *"plus-end OR lattice motor site"*.  Plus-end-only makes the population at most one per
microtubule and the natural basis ``per_filament``; a lattice reading makes every one of the 73,500 arm
beads a candidate and the basis a per-length density.  These are different models, not different
parameters, so ``sites`` is a required argument here with no default and the module pre-decides nothing.

⚠ NO KINETIC LAW IS BUILT HERE, DELIBERATELY.  ``bond.py`` defers attach/detach, the free list and the
snapshot twins to *"the first family that has kinetics"*.  This connector declares ``kinetics=True`` and
``commit_on_accept=True``, and ``connector_joints.py`` records that a kinetic connector proposing
nothing *"is indistinguishable from a permanent weld, which is the one thing the charter says a
connector may never be"* — so a silent no-op would be worse than a refusal.  Whether ``mt_cortex_capture``
is the family that opens that machinery is a Lead/PI call across the fourteen concurrent sessions, not
this module's; it is raised, not built.  It is worth flagging that a processive motor is the WRONG first
customer for it: dynein walks, so its bound state carries a velocity and a stall force as well as an
on/off rate, and the passive slip-bond shape that a crosslinker family would motivate does not cover it.

⚠ AND THE CAPTURE RADIUS NEVER SETS THE COUNT.  ``__init__.DUPLICATE_CRITERION`` records D13's finding
that ``sf_cortex_transient``'s count is *"however many partners happened to fall inside a capture
radius"* — a placement number wearing a physiological label.  Here the radius only FILTERS candidates;
``BondCount.resolve`` decides how many of them become bonds, and if the filter yields fewer candidates
than the count resolves to, ``build_bond_family`` raises rather than quietly returning the placement
number.  The radius is itself unsourced (the ERM analogue is PI-GAP card E4, *"candidate lead: NOT
FOUND"*), so it too is a required argument.

Sanity Gate (recorded per CLAUDE.md, before first execution):
  * dimensional — ``capture_radius_um`` and ``rest_um`` [µm], ``stiffness_pn_per_um`` [pN/µm], positions
    [µm]; the count is [1/filament] under the proposed basis and dimensionless once resolved.  Nothing
    here carries a time or a rate, because nothing here is kinetic.
  * boundary — zero candidates inside the radius is returned as an empty pair array and NOT an error:
    a cell whose microtubules do not reach the cortex is a physical configuration.  A count that then
    resolves above zero is refused by ``build_bond_family``, which is the correct place for it.
  * conservation/invariant — every returned pair joins one MT node to one cortex node, asserted against
    the arena's claims; each cortical site is used at most once, because one dynein anchor holding two
    microtubules is a pairing nobody declared and a BOND is *"a stateful, identified pairing"*.  The
    MTOC is excluded from the candidate set: it belongs to no arm and is not a lattice site.
  * CFL/precision — no integration; float64 throughout.  The stiffness is stored per bond, so it enters
    the CFL bound through its endpoints' row sums exactly as ``bond.py`` requires.
  * sign sense — ``rest_um`` is a length and non-negative; the aster radiates outward so the plus-end is
    the LAST node of each strand (``world/build/microtubule.py`` sign-sense clause), which is the
    convention ``microtubule_sites`` reads.
  * measurement protocol — pure host bookkeeping and geometry.  Positions are passed IN as host arrays
    rather than read from the arena, so the readback is the caller's explicit between-steps transfer and
    this module cannot become a hidden device round-trip.  A device-resident broad phase is PHASE 4.

engine units: length µm, stiffness pN/µm.  Runtime: pure host bookkeeping, CPU-importable.
"""

from __future__ import annotations

from typing import Literal

import numpy as np
import numpy.typing as npt

from aleph.world.arena import Claim, Kind, WorldArena
from aleph.world.bond import BondCount, BondFamily, SourceClass, build_bond_family
from aleph.world.families import ConnectorGapError, FamilySpec

__all__ = [
    "CONNECTOR", "CHEMISTRY_CARD", "SPEC", "BLOCKED_BY",
    "build", "microtubule_sites", "pair_to_nearest_cortex", "build_capture_family",
]

CONNECTOR = "mt_cortex_capture"
CHEMISTRY_CARD = "microtubule_cortical_dynein"
MT_POPULATION = "microtubule"
CORTEX_POPULATION = "cortex"

#: Which microtubule nodes may be captured. No default: the contract says "plus-end OR lattice motor
#: site" and those are two different models.
CandidateSites = Literal["plus_end", "lattice"]

#: The unanswered questions, in ``BondCount``'s own order, plus the three this connector adds. Lead
#: collects these into the decision queue; the module docstring carries what was searched for each.
BLOCKED_BY: tuple[str, ...] = (
    "how many — no count of cortical dynein capture sites or bound MT-cortex crossbridges exists: "
    "Contract-Graph `parameter` returns 0 rows for dynein/capture/mt_cortex/microtubule, and of "
    "12,810 paper_chunks the 3 'cortical dynein' hits are one axonal MT-bundle paper with no density",
    "against what — the proposed `per_filament` basis is the only one this arena can measure, because "
    "the cortex population claims no FACE and so has no area; an areal cortical density would have to "
    "resolve against the MEMBRANE's triangle area, which is the ERM defect's shape",
    "which cell — the per-filament support is itself a PI-GAP: PER_CELL_STRUCTURE_COUNTS_2026-08-20 §3 "
    "finds no per-cell MT count in the corpus and proposes 250-600/epithelial cell as a TEST AXIS, "
    "which world/build/microtubule.py ships as 500 under source_class=PI_GAP",
    "whose authority — no ratified source, no KB row, and no PI_GAP evidence card for this connector",
    "which MT sites are candidates — the contract says 'plus-end OR lattice motor site'; plus-end-only "
    "and lattice are different models with different bases, not two settings of one parameter",
    "capture radius — unsourced; the ERM analogue is card E4, 'candidate lead: NOT FOUND'",
    "kinetics — this edge declares kinetics=True and bond.py defers attach/detach to 'the first family "
    "that has kinetics'. A processive motor is the wrong first customer: dynein's bound state carries "
    "a walk velocity and a stall force, not only an on/off rate. PI call, not this module's",
)

SPEC = FamilySpec(
    CONNECTOR,
    (MT_POPULATION, CORTEX_POPULATION),
    "per_filament",
    blocked_by=BLOCKED_BY,
    source="",
)

_NOTE = (
    "SPEC.basis='per_filament' is this module's PROPOSAL, not an answer — it is recorded because "
    "FamilySpec requires one of four bases and 'unknown' is not among them. See the module docstring "
    "for what was searched. Nothing is added to DUPLICATE_OF: chemistry card "
    f"{CHEMISTRY_CARD!r} occurs exactly once across the 38 declared connectors."
)


def build(arena: WorldArena) -> BondFamily:
    """The census entry point. It raises today, and the refusal is this module's product.

    Args:
        arena: the world the family would be claimed out of. Unused — the refusal is prior to any
            geometry, and taking the argument keeps the signature the day an answer arrives.

    Raises:
        ConnectorGapError: always, naming every unanswered question in :data:`BLOCKED_BY`.
    """
    raise ConnectorGapError(CONNECTOR, list(BLOCKED_BY), note=_NOTE)


def _sole_claim(arena: WorldArena, population: str, kind: Kind) -> Claim:
    """The one claim a population holds of ``kind``, or a refusal naming how many it actually holds."""
    claims = arena.claims(population, kind)
    if len(claims) != 1:
        raise ValueError(
            f"{CONNECTOR}: expected exactly one {kind} claim for population {population!r}, found "
            f"{len(claims)}. This module reads the arithmetic topology world/build/microtubule.py "
            "builds (one contiguous range per population); a split population needs explicit arrays."
        )
    return claims[0]


def microtubule_sites(
    arena: WorldArena, sites: CandidateSites, *, population: str = MT_POPULATION
) -> npt.NDArray[np.int64]:
    """GLOBAL node indices of the microtubule sites a cortical dynein may capture.

    The topology is arithmetic, not an array: ``world/build/microtubule.py`` lays strand ``k`` bead
    ``i`` at ``lo + k*n_per + i`` and writes the MTOC last, so ``n_per`` is derived from the NODE and
    STRAND claims rather than passed in — a mismatch means the population was not built by that module
    and is refused rather than indexed on a guess.

    Args:
        arena: the world holding the population.
        sites: ``"plus_end"`` for the tip of each arm — the last bead, since the aster radiates outward
            — or ``"lattice"`` for every arm bead. The MTOC is in neither: it belongs to no arm.
        population: the microtubule population name.

    Returns:
        ``(S,)`` global node indices, ascending.

    Raises:
        ValueError: on an unrecognised ``sites``, or when the node and strand claims do not satisfy
            ``n_nodes == n_strand * n_per + 1``.
    """
    if sites not in ("plus_end", "lattice"):
        raise ValueError(
            f"{CONNECTOR}: sites must be 'plus_end' or 'lattice'; got {sites!r}. There is no default — "
            "the contract's 'plus-end OR lattice motor site' is two models, and picking one silently "
            "is the modelling decision this package exists to refuse."
        )
    nodes = _sole_claim(arena, population, Kind.NODE)
    strands = _sole_claim(arena, population, Kind.STRAND)
    n_strand = strands.count
    n_per, rem = divmod(nodes.count - 1, n_strand)
    if rem or n_per < 1:
        raise ValueError(
            f"{CONNECTOR}: {population} has {nodes.count} nodes over {n_strand} strands, which is not "
            "n_strand * nodes_per_strand + 1 (the arms plus one MTOC). This module cannot infer the "
            "topology of a population it did not recognise."
        )
    k = np.arange(n_strand, dtype=np.int64)
    if sites == "plus_end":
        return nodes.lo + k * n_per + (n_per - 1)
    return nodes.lo + np.arange(n_strand * n_per, dtype=np.int64)


def pair_to_nearest_cortex(
    mt_index: npt.ArrayLike,
    mt_position_um: npt.ArrayLike,
    cortex_claim: Claim,
    cortex_position_um: npt.ArrayLike,
    *,
    capture_radius_um: float,
) -> npt.NDArray[np.int64]:
    """Candidate pairs: each MT site to its nearest cortex node within ``capture_radius_um``.

    **This is a FILTER, never a count.** It says which pairings are geometrically possible; how many of
    them become bonds is :meth:`~aleph.world.bond.BondCount.resolve`'s answer, taken from the head of
    this list in the ascending MT order below. A count resolving above the number of candidates is a
    refusal in ``build_bond_family``, not a silent truncation to the placement number.

    One cortical site is used at most once: where two microtubules are nearest the same cortex node,
    the lower MT index keeps it and the other is dropped. One dynein anchor holding two microtubules is
    a pairing nobody declared.

    Args:
        mt_index: ``(S,)`` global node indices from :func:`microtubule_sites`, ascending.
        mt_position_um: ``(S, 3)`` host positions [µm] for those indices, in the same order.
        cortex_claim: the cortex NODE claim, so the returned indices are GLOBAL.
        cortex_position_um: ``(C, 3)`` host positions [µm] over ``cortex_claim``'s range, in range order
            — i.e. exactly what ``arena.download_nodes(cortex_claim, "position")`` returns.
        capture_radius_um: the reach beyond which a plus-end is not captured. Required and unsourced;
            see :data:`BLOCKED_BY`.

    Returns:
        ``(P, 2)`` global ``[mt_node, cortex_node]`` pairs, deterministic and ascending in the MT index.
        Empty when nothing is in reach — a physical answer, not an error.

    Raises:
        ValueError: on a non-finite or non-positive radius, a shape mismatch, or a cortex position
            block whose length is not ``cortex_claim.count``.
    """
    from scipy.spatial import cKDTree

    if not np.isfinite(capture_radius_um) or capture_radius_um <= 0.0:
        raise ValueError(f"{CONNECTOR}: capture_radius_um must be finite and positive")
    mt = np.asarray(mt_index, np.int64).reshape(-1)
    mt_xyz = np.asarray(mt_position_um, np.float64).reshape(-1, 3)
    cx_xyz = np.asarray(cortex_position_um, np.float64).reshape(-1, 3)
    if mt_xyz.shape[0] != mt.size:
        raise ValueError(f"{CONNECTOR}: {mt.size} MT sites but {mt_xyz.shape[0]} positions")
    if cx_xyz.shape[0] != cortex_claim.count:
        raise ValueError(
            f"{CONNECTOR}: {cx_xyz.shape[0]} cortex positions for a claim of {cortex_claim.count}. "
            "Pass the whole claimed range, in range order."
        )
    if mt.size == 0:
        return np.empty((0, 2), np.int64)

    dist, local = cKDTree(cx_xyz).query(mt_xyz, distance_upper_bound=float(capture_radius_um))
    in_reach = np.isfinite(dist)
    mt, local = mt[in_reach], local[in_reach].astype(np.int64)
    # One cortical site, one microtubule: np.unique keeps the FIRST occurrence, and the MT order is
    # ascending, so the tie goes to the lower MT index and the result stays deterministic.
    _, first = np.unique(local, return_index=True)
    keep = np.sort(first)
    return np.stack([mt[keep], cortex_claim.lo + local[keep]], axis=1)


def build_capture_family(
    arena: WorldArena,
    *,
    count: BondCount,
    support: float,
    pairs: npt.ArrayLike,
    rest_um: npt.ArrayLike,
    stiffness_pn_per_um: float,
) -> BondFamily:
    """Build the capture family once the four questions have been ANSWERED. Nothing here defaults.

    Kept separate from :func:`build` on purpose: this is the path the day a PI ruling lands, and it is
    exercised by the module self-check with a count that says ``CONVENIENCE`` in its own provenance, so
    the geometry is tested without a physiological claim being made anywhere.

    Args:
        arena: the world to claim the BOND range from.
        count: the ratified :class:`~aleph.world.bond.BondCount`.
        support: the MEASURED support to resolve it against — under the proposed ``per_filament`` basis,
            the microtubule STRAND claim count.
        pairs: ``(P, 2)`` candidates from :func:`pair_to_nearest_cortex`.
        rest_um: per-pair or scalar rest length [µm] — the bound crossbridge's unstrained extension.
        stiffness_pn_per_um: the crossbridge stiffness [pN/µm].

    Returns:
        The built :class:`~aleph.world.bond.BondFamily`, named for the connector and carrying
        :data:`CHEMISTRY_CARD`.

    Raises:
        ValueError: from ``build_bond_family``, or if the built bonds do not in fact join the
            microtubule and cortex populations — a pair list assembled from the wrong ranges would
            otherwise build a perfectly valid family of the wrong physics.
    """
    family = build_bond_family(
        arena, CONNECTOR, chemistry_card=CHEMISTRY_CARD, count=count, support=support,
        pairs=pairs, rest_um=rest_um, stiffness_pn_per_um=stiffness_pn_per_um,
    )
    joined = set(family.component_pairs(arena))
    expected = {tuple(sorted((MT_POPULATION, CORTEX_POPULATION)))}
    if family.n_bonds and joined != expected:
        raise ValueError(
            f"{CONNECTOR}: the built bonds join {sorted(joined)}, not {sorted(expected)}. The component "
            "pair is DERIVED from the arena's ID ranges, so this is the pair list being wrong, not the "
            "derivation."
        )
    return family


def _demo() -> None:
    """Self-check: the module refuses by default, and the geometry it would use is correct.

    No physiological number is asserted anywhere below. The count carries ``CONVENIENCE`` and says so in
    its own provenance, which is what lets the pairing be tested without a claim being made.
    """
    # 1. The product: a refusal that names every unanswered question.
    arena = WorldArena(capacity={Kind.NODE: 400, Kind.STRAND: 40, Kind.BOND: 100})
    try:
        build(arena)
    except ConnectorGapError as exc:
        assert exc.connector == CONNECTOR
        assert len(exc.missing) == len(BLOCKED_BY) == 7
        assert "how many" in str(exc) and "kinetics" in str(exc)
    else:  # pragma: no cover
        raise AssertionError("mt_cortex_capture must refuse: no count exists")
    assert not SPEC.buildable and SPEC.populations == (MT_POPULATION, CORTEX_POPULATION)

    # 2. A cortex shell on the unit sphere, and an aster of 4 arms x 5 beads + 1 MTOC, laid the way
    #    world/build/microtubule.py lays it: bead i of arm k at lo + k*n_per + i, tip last.
    n_cortex, n_strand, n_per = 200, 4, 5
    cortex = arena.claim(CORTEX_POPULATION, Kind.NODE, n_cortex)
    mt_nodes = arena.claim(MT_POPULATION, Kind.NODE, n_strand * n_per + 1)
    arena.claim(MT_POPULATION, Kind.STRAND, n_strand)

    golden = np.pi * (3.0 - np.sqrt(5.0))
    t = np.arange(n_cortex, dtype=np.float64)
    z = 1.0 - 2.0 * (t + 0.5) / n_cortex
    r = np.sqrt(np.maximum(0.0, 1.0 - z * z))
    cortex_xyz = np.stack([r * np.cos(golden * t), r * np.sin(golden * t), z], axis=1)

    arms = np.array([[1.0, 0, 0], [-1.0, 0, 0], [0, 1.0, 0], [0, -1.0, 0]], np.float64)
    # Arm 3 is SHORT: its tip stops at 0.5, well inside the shell, so it must not be captured.
    reach = np.array([1.0, 1.0, 1.0, 0.5])
    mt_xyz = np.concatenate(
        [np.outer(np.arange(1, n_per + 1) / n_per * reach[k], arms[k]) for k in range(n_strand)]
        + [np.zeros((1, 3))]  # the MTOC
    )

    # 3. Plus-ends are the arm TIPS and the MTOC is in neither candidate set.
    tips = microtubule_sites(arena, "plus_end")
    lattice = microtubule_sites(arena, "lattice")
    assert tips.tolist() == [mt_nodes.lo + k * n_per + n_per - 1 for k in range(n_strand)]
    assert lattice.size == n_strand * n_per
    mtoc = mt_nodes.hi - 1
    assert mtoc not in set(tips.tolist()) | set(lattice.tolist()), "the MTOC is not a lattice site"
    for bad in ("tip", "", None):
        try:
            microtubule_sites(arena, bad)  # type: ignore[arg-type]
        except ValueError as exc:
            assert "no default" in str(exc)
        else:  # pragma: no cover
            raise AssertionError(f"sites={bad!r} must refuse")

    # 4. The radius FILTERS. Three long arms reach the shell; the short one does not.
    local = tips - mt_nodes.lo
    cand = pair_to_nearest_cortex(tips, mt_xyz[local], cortex, cortex_xyz, capture_radius_um=0.2)
    assert cand.shape == (3, 2), f"three arms in reach, got {cand.shape[0]}"
    assert np.all(cand[:, 0] == tips[:3]) and np.all(np.diff(cand[:, 0]) > 0), "ascending, deterministic"
    assert np.all((cand[:, 1] >= cortex.lo) & (cand[:, 1] < cortex.hi))
    # Nothing in reach is an ANSWER, not an error.
    assert pair_to_nearest_cortex(
        tips, mt_xyz[local], cortex, cortex_xyz, capture_radius_um=1e-6).shape == (0, 2)
    # One cortical site is used at most once: two arms aimed at the same point give ONE pair.
    twin = np.array([cortex_xyz[0] * 0.99, cortex_xyz[0] * 0.99], np.float64)
    assert pair_to_nearest_cortex(
        tips[:2], twin, cortex, cortex_xyz, capture_radius_um=0.2).shape == (1, 2)

    # 5. The count governs, the placement does not. Two per filament over 1.5 filaments of support is
    #    3 bonds, and the candidate list happens to hold exactly 3 — so a count of 4 must REFUSE
    #    rather than fall back to the 3 the geometry offered.
    per_fil = BondCount(
        basis="per_filament", value=2.0, scope="NOT PHYSIOLOGICAL — module self-check only",
        source_class=SourceClass.CONVENIENCE,
        provenance="a number chosen to exercise the resolver; no cell has been claimed to have it")
    fam = build_capture_family(
        arena, count=per_fil, support=1.5, pairs=cand, rest_um=0.02, stiffness_pn_per_um=50.0)
    assert fam.n_bonds == 3 and fam.chemistry_card == CHEMISTRY_CARD
    assert fam.component_pairs(arena) == {(CORTEX_POPULATION, MT_POPULATION): 3}
    row = fam.provenance_row()
    assert row["support_unit"] == "filaments" and row["source_class"] == "CONVENIENCE"

    greedy = BondCount(basis="explicit", value=4.0, scope="NOT PHYSIOLOGICAL — module self-check only",
                       source_class=SourceClass.CONVENIENCE,
                       provenance="one more than the geometry offers, to prove the refusal")
    try:
        build_capture_family(arena, count=greedy, support=1.0, pairs=cand, rest_um=0.02,
                             stiffness_pn_per_um=50.0)
    except ValueError as exc:
        assert "only 3 candidate pairs" in str(exc), exc
    else:  # pragma: no cover
        raise AssertionError("a count above the candidates must refuse, never truncate to placement")

    print(f"mt_cortex_capture self-check OK — REFUSES with {len(BLOCKED_BY)} unanswered; {row}")


if __name__ == "__main__":
    _demo()
