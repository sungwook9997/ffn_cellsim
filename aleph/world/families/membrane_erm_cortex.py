"""``membrane_erm_cortex`` — the ERM tether population, on a density axis the PI called TEMPORARY.

This is the connector :mod:`aleph.world.bond` was written ABOUT. Its docstring names this edge as the
defect the whole primitive is shaped around: the incumbent declares ONE connector between two
components where the physics is a POPULATION of ezrin/radixin/moesin tethers whose count is an areal
density, and because the declaration hid the population **nobody had to answer how many**. So this
module is the test of whether that shape actually refuses, and it did: it refused for four weeks, and
PI decision 3 of 2026-08-21 answered it with an axis that says out loud that it expires.

WHAT IS ANSWERED, AND IT IS MOST OF IT.  Three of :class:`~aleph.world.bond.BondCount`'s four questions
have answers, and the mechanism is not in doubt:

* *against what* — the MEASURED mesh area, ``ClosedSurface.area0_total_um2``, never the analytic sphere.
  ``bond.py`` singles this out as the one thing the incumbent's ERM builder got right, and it is why an
  areal density here stays invariant when the continuum membrane resolution changes.
* *which cell* — MCF7 interphase adherent 37 °C, the target of every axis in ``world/build``.
* *whose authority, for the FORCE* — ``k_erm`` = 4.6e3 pN/µm is `KB-3.B1.6`, PI-ratified 2026-07-21
  (Braunger 2014 JBC 289:9833, ``source_audit`` verdict OK, checked in ``kb.duckdb`` 2026-08-20;
  colloidal-probe AFM series-spring deconvolution of the single ezrin–F-actin bond). The per-molecule
  unbinding force ≈ 50 pN is activation-INDEPENDENT, and PIP2 sets the bond NUMBER rather than the
  per-bond force — which is exactly why the number cannot be inferred from the force.
* *how many* — **a DECLARED, TEMPORARY axis as of PI decision 3**, and nothing stronger; see below.
  Nothing in the Contract-Graph answers it, and no mesh may be allowed to.

**WHAT CHANGED - PI decision 3, 2026-08-21** (`docs/v2_audit/PI_DECISIONS_2026-08-21.md` §3): the
density becomes a **DECLARED AXIS - a band and a test point - and the PI's own word for it is
TEMPORARY.**  *"This axis expires when a sourced areal density arrives, and any result standing on it
carries that."*  So the word is stamped into :data:`ERM_DENSITY_AXIS` and travels in every
``provenance_row`` this module produces; a run that quotes a number resting on it is quoting an axis
with an expiry date, and the artifact says so without anyone having to remember.

⚠ **WHAT THE AXIS IS NOT - 231.8/µm² IS NOT RATIFIED, AND IS REFUSED BY NAME.**  That number is
``subdiv 7 = 163,842`` vertices over the measured shell of a 7.5 µm cell, ``4*pi*(7.5)² = 706.9 µm²``,
and it reproduces to the fourth digit: **163842 / 706.858 = 231.789**.  It is a mesh artefact wearing a
physiological label - the ORIGINAL ERM defect and the reason the whole class has a name - and the
incumbent's comment calling it "physiological 235/µm²" is that artefact being read back as a datum.
⚠ **And at the other end, the three L0 resting runs of 2026-08-20 ran ``subdiv 3`` = 642 vertices =
0.908/µm²**, two and a half orders below even the diagnostic-only zebrafish proxy, and below the
capacity minimum this module's band starts at.  One-per-vertex is not a density with an unfortunate
value; it is not a density, at any subdivision.  :func:`refuse_vertex_density` refuses the whole class
STRUCTURALLY - it compares the RESOLVED count against the membrane's own vertex count, so it catches
subdiv 3, 6, 7 and every level nobody has run yet, rather than blacklisting three magnitudes.

**WHERE THE BAND'S TWO ENDS COME FROM, AND WHAT THEY DO AND DO NOT SUPPORT.**  Decision 6 of the same
sitting refused a band whose endpoints could not carry one, so these are stated the same way:

* **Lower - the capacity condition, and it is genuinely a lower bound**, which is the one thing a
  band's lower end is entitled to be.  ``VG-ERM-capacity`` records a back-computed minimum of
  **3.2649/µm²**, and that value was computed with ``erm_rupture_force = 11.4347 pN`` - which the SAME
  card corrects as the CONTINUUM tube-extraction scale (`KB-3.B1.4`), **not** the per-ERM force.  With
  the sourced single-bond figure (~50 pN, `KB-3.B1.6`, Braunger 2014) the same condition gives
  **0.7467/µm²**, and that is the lower end.  ⚠ It assumes a bound fraction of 1; the bound fraction is
  one of the two items still blocked below, and it can only RAISE this end.
* **Upper - 600/µm², zebrafish embryonic ectoderm**, the only DIRECT ERM areal density anywhere in the
  record.  Other-species embryo, which the PI ruled diagnostic-only, so it stands here as a band end
  and never as the value.
* **Test point - the band's geometric centre, 21.166/µm².**  A declared convention over a 2.9-decade
  band, not a measurement, and saying so is the point: it is WHY the axis is temporary.

⚠ **AND THE AXIS IS A CAPACITY AXIS, NOT AN ENGAGED ONE.**  ``bond.py`` separated the two on
2026-08-21 (`f99400a6`): a :class:`~aleph.world.bond.BondCount` says how many tethers the structure
HAS, and how many are attached at a given step is occupancy, which is emergent from ``k_on``/``k_off``
and may never be imposed — ``params_i0b3.yaml`` already forbade imposing an engaged fraction from a
duty ratio.  ERM tethers attach and shed; ρ_ERM is the population, the **bound/active fraction** is the
occupancy, and they are two of the bundle's three items precisely because they are different
questions.  ``provenance_row`` carries ``n_bonds_is`` so a reader cannot take one for the other, and
``partner_is_dynamic`` stays **False** here for the same reason the family is blocked: nothing detaches
yet.  ⚠ It has to flip to ``True`` in the same change that lands ``erm_bell_kinetics`` — once the KMC
sheds a tether, ``node_j`` is state and the family owes a ``live_partner_owner``.

**WHAT IS STILL BLOCKED, AND WHY THE MODULE STILL REFUSES BY DEFAULT.**  ``VG-ERM-capacity`` approves
the bundle as ONE unit - *"areal density + bound/active fraction + Bell {k_on, k_off=1.3/s, F0, capture
radius}"*.  Decision 3 answers the density and **only** the density.  The bound fraction and the Bell
kinetics are untouched, so ``count`` still has no default and a caller that passes nothing still gets
:class:`~aleph.world.families.ConnectorGapError` naming the two that remain.  What decision 3 makes
possible is a caller that passes :func:`erm_density_axis` **explicitly**: a static tether population on
a declared temporary axis, with no kinetics - which is what PHASE 3 needs and all it may have.

NOT A DUPLICATE, AND THE PI ALREADY SAID SO IN WRITING.  ``membrane_cortex_contact`` joins the SAME
population pair, which under this package's :data:`~aleph.world.families.DUPLICATE_CRITERION` proves
nothing either way — identity is the CHEMISTRY CARD. Here the two do not merely carry different cards:
the contact edge has no chemistry to card. ``engine/contracts.py:981`` declares it
``ConnectorFamily.CONTACT`` with ``kinetics=False``, and the contract's own note, recording PI decision
D3 option C of 2026-07-28, settles it outright: *"The two edges are complementary, not redundant: ERM
is a tension tether that is force-free in compression, so this CONTACT edge is the only one that can
carry an outward pressure."* It is the turgor transmission path — a position-only steric interaction,
which by ``bond.py``'s own rule needs no bond declaration at all, since "field interaction (steric,
drag) needs no declaration because it depends on position alone". This connector is its opposite on
every axis: unilaterally TENSILE, a persistent IDENTIFIED pairing, stateful, and kinetic. Not entered
in ``DUPLICATE_OF``; the shared population pair is not evidence that it should be.

⚠ THIS IS THE FIRST FAMILY WITH KINETICS, AND IT IS NOT MINE TO WRITE.  ``bond.py`` defers attach /
detach / free lists / snapshot twins to "the first family that has kinetics". ERM is that family: the
incumbent already drives a Bell-SLIP KMC over it (``erm_bell_kmc_pair_kernel``), where the off-rate
rises with tensile load and the shedding IS the mechanistic bleb onset. **This module writes no kinetic
law.** It is raised to the PI, together with the count it is bundled with, rather than built.

WHAT WOULD CLOSE THIS, in the PI's own priority order: (1) MCF7 quantitative proteomics × cortical /
membrane localization fraction; (2) another human epithelial DIRECT areal density; (3) other-species
embryo, diagnostic-only. Registering any of them is a ``ModelContract``/``ValidationGate`` change and is
PI-authored, never auto-created.

CONVENTIONS DECLARED, NOT DERIVED — recorded here so they are visible rather than discovered later:

* ``rest_um`` is the AS-BUILT pair separation, so the family is force-free in the built configuration.
  That is the incumbent's convention and it is defensible, but it is a choice with a consequence: the
  resting state carries ZERO ERM prestress by construction. It is also not the molecule — real ezrin is
  ~20 nm, while the arena's built radial gap is ``h_cortex/2`` = **0.100 µm**, five times that, because
  the cortex shell is placed half a cortical thickness inside the membrane. A rest length that is a
  geometry artefact is the same class of defect as a count that is a mesh artefact; it is smaller only
  because ``rest_um`` is per-bond state rather than the population size.
* Cortex endpoints are midpoint-stratified in arena node order (the incumbent's rule, which needs no
  seed and reuses no cortex node). That is a uniform sample over the SHELL only because the cortex
  builder draws filament frames isotropically — an inherited property, stated so it can be rechecked
  if that builder changes.

Sanity Gate (recorded per CLAUDE.md, before first execution):
  * dimensional — density [1/µm²] × support [µm²] → count [1]; ``rest_um`` [µm]; ``k_erm`` [pN/µm].
    The basis is ``areal`` and :class:`~aleph.world.bond.BondCount` rejects a support carrying any
    other unit by name.
  * boundary — a resolved count of zero is a legal answer and is returned as an empty family; a count
    exceeding the available cortex nodes raises with both numbers rather than reusing an endpoint, and
    a pair separated by more than the reach raises rather than being silently dropped.
  * conservation/invariant — every endpoint is a GLOBAL arena node id inside the live NODE prefix
    (asserted by ``build_bond_family``); each cortex node is used at most once, so no cortex degree of
    freedom carries two molecular states; the force this family will later feed is the two-array
    ``-f``/``+f`` pair, which is Newton's 3rd law across never-merged owner arrays.
  * CFL/precision — no integration here; float64 throughout. ``k_erm`` = 4.6e3 pN/µm enters the CFL
    bound through the row sum at its endpoints, which is why ``bond.py`` stores stiffness per bond.
  * sign sense — the tether is unilateral: tensile only, force-free when compressed or at rest.
    ``rest_um`` is a length and is non-negative by construction (a norm).
  * measurement protocol — host-side construction from downloaded positions; no device is touched and
    nothing is uploaded. Positions are passed in rather than read from the arena, so this module is
    CPU-importable and testable without a card.

engine units: length µm, stiffness pN/µm.  Runtime: pure host bookkeeping, CPU-importable.
"""

from __future__ import annotations

import math

import numpy as np
import numpy.typing as npt

from aleph.world.bond import BondCount, BondFamily, SourceClass, build_bond_family
from aleph.world.families import ConnectorGapError, FamilySpec

__all__ = [
    "CONNECTOR", "CHEMISTRY_CARD", "K_ERM_PN_PER_UM", "SPEC",
    "ERM_DENSITY_BAND_PER_UM2", "ERM_DENSITY_TEST_POINT_PER_UM2", "ERM_DENSITY_AXIS",
    "ERM_AXIS_EXPIRES_WHEN", "erm_density_axis", "refuse_vertex_density",
    "erm_candidate_pairs", "build_membrane_erm_cortex",
]

CONNECTOR = "membrane_erm_cortex"

#: The chemistry, which is what identifies a family — never the declared component pair.
CHEMISTRY_CARD = "ezrin_braunger2014"

#: Single ezrin–F-actin linker stiffness [pN/µm]. KB-3.B1.6, PI-ratified 2026-07-21; Braunger 2014 JBC
#: 289(14):9833 (source_audit OK). Retires the provisional 100 pN/µm, a bleb-PDE continuum value ~46x
#: too soft. ⚠ ``laws/network_warp.py:817`` still carries ``k_erm_pN_um`` defaulting to 50.0, labelled
#: there as a CFL convenience and NOT single-molecule ezrin — that is a different path, and the
#: disagreement is reported here rather than silently reconciled.
K_ERM_PN_PER_UM = 4.6e3

#: The three names the PI bundled as ONE contract, kept whole because they were approved whole.
#: ⚠ ``rho_ERM_per_um2`` is no longer in :data:`SPEC.blocked_by` — PI decision 3 answered it with the
#: temporary axis below — but it stays listed HERE, because the bundle is the contract unit and the
#: axis is not a source. What came off the blocked list is the blank, not the requirement.
_BUNDLE: tuple[str, ...] = (
    "rho_ERM_per_um2",
    "erm_bound_fraction",
    "erm_bell_kinetics",
)

#: What is still unanswered after decision 3, and therefore what a caller passing no count is told.
_STILL_BLOCKED: tuple[str, ...] = ("erm_bound_fraction", "erm_bell_kinetics")

#: The registered ``VG-ERM-capacity`` back-computed minimum [1/µm²] — necessary-capacity, a LOWER
#: BOUND, and explicitly disqualified from being adopted as the value. Recorded at
#: ``tests/ac/cell/test_preload_contract.py:23`` and grid-invariant by that file's own third test.
VG_ERM_CAPACITY_MIN_RECORDED_PER_UM2 = 3.2649139201533086

#: The force that minimum was computed with [pN] — and ⚠ it is the CONTINUUM tube/tether-extraction
#: scale ``f_t = 2*pi*sqrt(2*kappa_m*(T_m + gamma_MCA))`` (`KB-3.B1.4`), **not** a single-molecule
#: force. AC_DECISION_CARDS 2026-07-22 P2 Card 2 §1 corrects exactly this mislabelling.
F_CONTINUUM_TETHER_PN = 11.43470677829702

#: The single ezrin–F-actin unbinding force [pN] — `KB-3.B1.6`, Braunger 2014, activation-INDEPENDENT.
#: This is the force a per-tether capacity condition is entitled to use.
F_SINGLE_ERM_UNBINDING_PN = 50.0

#: Zebrafish embryonic ectoderm areal density [1/µm²] — the ONLY direct ERM areal density in the
#: record, other-species embryo, ruled diagnostic-only. Stands as a band END, never as the value.
ZEBRAFISH_ERM_PER_UM2 = 600.0

#: The declared axis [1/µm²]. The lower end is the capacity condition re-evaluated with the SOURCED
#: single-bond force (the recorded minimum scales as 1/f), at bound fraction 1 — the bound fraction is
#: still blocked and can only raise it. The upper end is the zebrafish proxy.
ERM_DENSITY_BAND_PER_UM2: tuple[float, float] = (
    VG_ERM_CAPACITY_MIN_RECORDED_PER_UM2 * F_CONTINUUM_TETHER_PN / F_SINGLE_ERM_UNBINDING_PN,
    ZEBRAFISH_ERM_PER_UM2,
)

#: The test point [1/µm²]: the band's GEOMETRIC centre, which is the scale-free choice over a band
#: spanning 2.9 decades. A declared convention, not a measurement — and that is why the axis expires.
ERM_DENSITY_TEST_POINT_PER_UM2 = math.sqrt(
    ERM_DENSITY_BAND_PER_UM2[0] * ERM_DENSITY_BAND_PER_UM2[1])

#: The PI's own word, in the artifact rather than in a session note.
ERM_AXIS_EXPIRES_WHEN = (
    "TEMPORARY (PI decision 3, 2026-08-21, the PI's own word): this axis expires the moment a sourced "
    "areal density arrives, and any result standing on it carries that. Closing order, PI's priority: "
    "(1) MCF7 quantitative proteomics x cortical/membrane localization fraction, (2) another human "
    "epithelial DIRECT areal density, (3) other-species embryo, diagnostic-only. Registering one is a "
    "ModelContract/ValidationGate change and is PI-authored, never auto-created."
)

#: The whole axis as one row, for a driver that records what it ran on.
ERM_DENSITY_AXIS: dict[str, object] = {
    "axis": "rho_ERM_per_um2",
    "status": "TEMPORARY",
    "band_per_um2": ERM_DENSITY_BAND_PER_UM2,
    "test_point_per_um2": ERM_DENSITY_TEST_POINT_PER_UM2,
    "band_low_basis": (
        f"VG-ERM-capacity recorded minimum {VG_ERM_CAPACITY_MIN_RECORDED_PER_UM2:.4f}/um^2 was computed "
        f"with f={F_CONTINUUM_TETHER_PN:.4f} pN, the CONTINUUM tube-extraction scale (KB-3.B1.4) which "
        f"P2 Card 2 section 1 corrects as NOT the per-ERM force. rho_min scales as 1/f, so at the "
        f"sourced single-bond {F_SINGLE_ERM_UNBINDING_PN:.0f} pN (KB-3.B1.6, Braunger 2014) the same "
        f"necessary condition gives {ERM_DENSITY_BAND_PER_UM2[0]:.4f}/um^2, at bound fraction 1. The "
        "bound fraction is still blocked and can only RAISE this end. NECESSARY-only: a lower bound is "
        "what a band's lower end is entitled to be, and nothing more."
    ),
    "band_high_basis": (
        f"zebrafish embryonic ectoderm {ZEBRAFISH_ERM_PER_UM2:.0f}/um^2 — the only DIRECT ERM areal "
        "density in the record; other-species embryo, PI-ruled diagnostic-only, so it is a band END."
    ),
    "test_point_basis": (
        "the band's GEOMETRIC centre over 2.9 decades. A declared convention, NOT a measurement, and "
        "not the value for any cell — which is what makes the axis temporary rather than provisional."
    ),
    "not_ratified": (
        "231.789/um^2 (= subdiv 7, 163,842 vertices / 706.858 um^2) is a MESH ARTEFACT and is NOT this "
        "axis, whatever the incumbent's 'physiological 235/um^2' comment says; 0.908/um^2 (= subdiv 3, "
        "642 vertices) is the same artefact at the resolution the 2026-08-20 L0 runs actually ran, and "
        "it sits BELOW this band. One-per-vertex is not a density at any subdivision."
    ),
    "count_is": (
        "CAPACITY, never occupancy (bond.py, f99400a6, 2026-08-21). rho_ERM is how many tethers the "
        "membrane-cortex interface HAS; how many are ATTACHED at a step is emergent from the Bell "
        "kinetics and may not be imposed - params_i0b3.yaml forbids imposing an engaged fraction from "
        "a duty ratio. The bound/active fraction is a SEPARATE bundle item and is still blocked."
    ),
    "expires_when": ERM_AXIS_EXPIRES_WHEN,
}

_GAP_NOTE = (
    "VG-ERM-capacity (PI-authored, AC_DECISION_CARDS 2026-07-22 P2 Card 2) approves the bundle as ONE "
    "unit: areal density + bound/active fraction + Bell {k_on, k_off=1.3/s, F0, capture radius}. PI "
    "decision 3 of 2026-08-21 answers the DENSITY only, and as a TEMPORARY declared axis - band "
    f"[{ERM_DENSITY_BAND_PER_UM2[0]:.4f}, {ERM_DENSITY_BAND_PER_UM2[1]:.0f}]/um^2, test point "
    f"{ERM_DENSITY_TEST_POINT_PER_UM2:.3f}/um^2 - not as a source. The bound fraction and the Bell "
    "kinetics are untouched, so this family still refuses a caller who passes no count. Pass "
    "erm_density_axis() EXPLICITLY for a static, kinetics-free tether population on the temporary "
    "axis. k_erm=4.6e3 pN/um and k_off0=1.3/s ARE sourced. ERM is also the FIRST family with kinetics "
    "and bond.py defers that shape on purpose."
)

SPEC = FamilySpec(
    connector=CONNECTOR,
    populations=("cortex", "membrane"),
    basis="areal",
    blocked_by=_STILL_BLOCKED,
    source=(
        "FORCE sourced, COUNT on a TEMPORARY declared axis (PI decision 3, 2026-08-21). k_erm 4.6e3 "
        "pN/um = KB-3.B1.6 PI-ratified 2026-07-21 (Braunger 2014 JBC, source_audit OK). rho_ERM band "
        f"[{ERM_DENSITY_BAND_PER_UM2[0]:.4f}, {ERM_DENSITY_BAND_PER_UM2[1]:.0f}]/um^2 with test point "
        f"{ERM_DENSITY_TEST_POINT_PER_UM2:.3f}/um^2, EXPIRING on a sourced areal density. NOT the "
        "icosphere vertex density (231.789/um^2 at subdiv 7, 0.908/um^2 at the subdiv 3 the 2026-08-20 "
        "L0 runs used) - that class is refused structurally by refuse_vertex_density()."
    ),
)


def erm_density_axis(value_per_um2: float | None = None) -> BondCount:
    """The declared TEMPORARY density axis as a :class:`BondCount`, at the test point or a band value.

    Args:
        value_per_um2: where on the axis to sit [1/µm²]. ``None`` takes
            :data:`ERM_DENSITY_TEST_POINT_PER_UM2`. Any value outside
            :data:`ERM_DENSITY_BAND_PER_UM2` is refused with both ends named — a swept axis whose
            sweep can leave its own band is not an axis.

    Returns:
        The count, carrying ``UNRATIFIED_PROXY`` and the word TEMPORARY in its provenance, so
        ``BondFamily.provenance_row`` puts both in the artifact without a driver remembering to.

    Raises:
        ValueError: on a non-finite value or one outside the declared band.
    """
    lo, hi = ERM_DENSITY_BAND_PER_UM2
    v = float(ERM_DENSITY_TEST_POINT_PER_UM2 if value_per_um2 is None else value_per_um2)
    if not math.isfinite(v) or not (lo <= v <= hi):
        raise ValueError(
            f"rho_ERM={value_per_um2!r} is outside the declared band [{lo:.4f}, {hi:.0f}]/um^2. The "
            "band's ends are a necessary-capacity lower bound and a diagnostic-only zebrafish proxy; "
            "a value outside them is not a point on this axis, and the two numbers that most often "
            "arrive here - 231.789 (subdiv 7) and 0.908 (subdiv 3) - are mesh artefacts rather than "
            "densities. Widening the band is a PI call, not a call-site one."
        )
    return BondCount(
        basis="areal", value=v,
        scope="MCF7 interphase adherent 37C — TEMPORARY declared axis, NOT a measurement of this cell",
        source_class=SourceClass.UNRATIFIED_PROXY,
        provenance=(
            f"TEMPORARY declared axis (PI decision 3, 2026-08-21): band "
            f"[{lo:.4f}, {hi:.0f}]/um^2, test point {ERM_DENSITY_TEST_POINT_PER_UM2:.3f}/um^2, this "
            f"run at {v:.3f}/um^2. {ERM_AXIS_EXPIRES_WHEN}"
        ),
    )


def refuse_vertex_density(resolved_count: int, n_membrane_vertices: int) -> None:
    """Refuse a count that is the membrane's VERTEX COUNT — the ERM defect, at any subdivision.

    The original defect was not the magnitude 231.8; it was ``one tether per membrane mesh vertex``,
    which produces a different wrong number at every resolution and looks physiological at each. So
    this compares the resolved count against the mesh rather than blacklisting magnitudes: it fires at
    subdiv 3, 6, 7 and at levels nobody has run yet.

    ⚠ It also fires on a legitimate density that happens to resolve to exactly the vertex count. That
    is deliberate and it is cheap — the axis moves by a hair and the coincidence is gone — while the
    defect it catches has cost this engine three run sets.

    Args:
        resolved_count: how many tethers the count resolved to.
        n_membrane_vertices: vertices in the membrane mesh the density was resolved against.

    Raises:
        ValueError: when the two are equal.
    """
    if resolved_count and resolved_count == int(n_membrane_vertices):
        raise ValueError(
            f"{CONNECTOR}: the density resolves to exactly {resolved_count} tethers against a mesh of "
            f"{n_membrane_vertices} vertices, i.e. ONE TETHER PER VERTEX. That is not a density with an "
            "unfortunate value — it is the discretisation, and it is the defect bond.py was written "
            "about: it yields 231.789/um^2 at subdiv 7 and 0.908/um^2 at subdiv 3, each of which has "
            "been read as physiological in this repository. Change the density, not the mesh."
        )


def erm_candidate_pairs(
    n_pairs: int,
    *,
    membrane_pos_um: npt.ArrayLike,
    membrane_lo: int,
    cortex_pos_um: npt.ArrayLike,
    cortex_lo: int,
    reach_um: float,
) -> tuple[npt.NDArray[np.int64], npt.NDArray[np.float64]]:
    """Pick ``n_pairs`` cortex endpoints and their membrane partners, as GLOBAL arena node ids.

    The cortex endpoints are midpoint-stratified in arena node order — ``floor((i + 0.5) * N / n)`` —
    which is unique for every ``n <= N``, needs no seed, and uses no cortex node twice. Each endpoint
    then terminates on the NEAREST membrane vertex, so several molecular states may load one continuum
    membrane degree of freedom while keeping independent cortex endpoints and rest lengths.

    Args:
        n_pairs: how many pairs to build. Comes from a resolved :class:`BondCount`, never typed.
        membrane_pos_um: ``(M, 3)`` host membrane vertex positions [µm], in arena order.
        membrane_lo: the membrane NODE claim's low index, added to make the ids global.
        cortex_pos_um: ``(N, 3)`` host cortex node positions [µm], in arena order.
        cortex_lo: the cortex NODE claim's low index.
        reach_um: the physical pairing reach [µm]. A pair longer than this is a build error.

    Returns:
        ``(pairs, rest_um)`` — ``(n_pairs, 2)`` global ids ordered ``(membrane, cortex)``, and the
        as-built separations [µm].

    Raises:
        ValueError: on a malformed position array, a non-positive reach, a count exceeding the
            available cortex nodes, or any pair exceeding ``reach_um``.
    """
    mem = np.asarray(membrane_pos_um, np.float64)
    cor = np.asarray(cortex_pos_um, np.float64)
    if mem.ndim != 2 or mem.shape[1] != 3 or cor.ndim != 2 or cor.shape[1] != 3:
        raise ValueError("membrane_pos_um and cortex_pos_um must both be (N, 3)")
    if not np.isfinite(reach_um) or reach_um <= 0.0:
        raise ValueError(f"reach_um must be finite and positive; got {reach_um!r}")
    n = int(n_pairs)
    if n < 0:
        raise ValueError(f"n_pairs must be nonnegative; got {n_pairs!r}")
    if n == 0:
        return np.zeros((0, 2), np.int64), np.zeros(0, np.float64)
    if n > cor.shape[0]:
        raise ValueError(
            f"{CONNECTOR}: the density asks for {n} tethers but the cortex offers {cor.shape[0]} "
            "nodes. Raise the mechanistic cortex resolution rather than reusing a cortex endpoint — a "
            "shared endpoint is two molecular states pretending to be one."
        )

    from scipy.spatial import cKDTree

    ci = np.floor((np.arange(n, dtype=np.float64) + 0.5) * cor.shape[0] / n).astype(np.int64)
    dist, mi = cKDTree(mem).query(cor[ci], k=1)
    dist = np.asarray(dist, np.float64)
    if np.any(dist > reach_um):
        n_bad = int(np.count_nonzero(dist > reach_um))
        raise ValueError(
            f"{CONNECTOR}: {n_bad} of {n} tethers exceed the {reach_um:g} um reach (max "
            f"{dist.max():.4g} um). Refine the membrane mesh rather than dropping explicit linkers — a "
            "dropped linker lowers the density silently, which is the defect this module exists for."
        )
    pairs = np.stack([np.asarray(mi, np.int64) + int(membrane_lo), ci + int(cortex_lo)], axis=1)
    return pairs, dist


def build_membrane_erm_cortex(
    arena,
    *,
    membrane_pos_um: npt.ArrayLike,
    membrane_lo: int,
    membrane_area_um2: float,
    cortex_pos_um: npt.ArrayLike,
    cortex_lo: int,
    reach_um: float,
    count: BondCount | None = None,
) -> BondFamily:
    """Build the ERM tether family, or refuse by name.

    Args:
        arena: the :class:`~aleph.world.arena.WorldArena` to claim the BOND range from.
        membrane_pos_um: ``(M, 3)`` host membrane vertex positions [µm].
        membrane_lo: the membrane NODE claim's low index.
        membrane_area_um2: the MEASURED total mesh area [µm²] — ``ClosedSurface.area0_total_um2``,
            never ``4 pi R^2``. Passing the analytic sphere would make the density resolution-dependent
            again, in the opposite direction.
        cortex_pos_um: ``(N, 3)`` host cortex node positions [µm].
        cortex_lo: the cortex NODE claim's low index.
        reach_um: the physical pairing reach [µm].
        count: the areal :class:`BondCount`. ``None`` raises :class:`ConnectorGapError` naming what is
            still unanswered rather than falling back to anything — there is no fallback in this module
            to fall back TO. Pass :func:`erm_density_axis` EXPLICITLY to build on the temporary axis;
            an axis a caller has to name is not a default, which is the whole difference.

    Returns:
        The built :class:`BondFamily`.

    Raises:
        ConnectorGapError: when no count is given — the bundle's other two items are still blocked.
        ValueError: if ``count`` is not an areal basis, if the count resolves to one tether per
            membrane vertex, or on any geometric build error.
    """
    if count is None:
        raise ConnectorGapError(CONNECTOR, list(_STILL_BLOCKED), _GAP_NOTE)
    if count.basis != "areal":
        raise ValueError(
            f"{CONNECTOR}: ERM is an areal population — PIP2 sets the number of attachments per unit "
            f"membrane area (KB-3.B1.6). Got basis {count.basis!r}."
        )
    resolved = count.resolve(float(membrane_area_um2))
    refuse_vertex_density(resolved, np.asarray(membrane_pos_um, np.float64).reshape(-1, 3).shape[0])
    pairs, rest = erm_candidate_pairs(
        resolved, membrane_pos_um=membrane_pos_um, membrane_lo=membrane_lo,
        cortex_pos_um=cortex_pos_um, cortex_lo=cortex_lo, reach_um=reach_um)
    return build_bond_family(
        arena, CONNECTOR, chemistry_card=CHEMISTRY_CARD, count=count,
        support=float(membrane_area_um2), pairs=pairs, rest_um=rest,
        stiffness_pn_per_um=K_ERM_PN_PER_UM,
    )


def _demo() -> None:
    """Self-check: the axis is declared and temporary, the mesh artefact is refused, and the two
    remaining bundle items still block a caller who brings nothing."""
    from aleph.world.arena import Kind, WorldArena
    from aleph.world.strand import build_strand

    # The spec is honest about what decision 3 did and did NOT close: the density came off the blocked
    # list, the bound fraction and the Bell kinetics did not.
    assert not SPEC.buildable and SPEC.blocked_by == _STILL_BLOCKED
    assert "rho_ERM_per_um2" not in SPEC.blocked_by and "rho_ERM_per_um2" in _BUNDLE
    assert SPEC.populations == ("cortex", "membrane") and SPEC.basis == "areal"

    # THE AXIS. Band ends derived from their own cited inputs, test point at the geometric centre, and
    # the PI's word in the row rather than in a commit message.
    lo, hi = ERM_DENSITY_BAND_PER_UM2
    assert abs(lo - 0.7466666666666667) < 1e-12 and hi == 600.0
    assert abs(ERM_DENSITY_TEST_POINT_PER_UM2 - (lo * hi) ** 0.5) < 1e-12
    assert "TEMPORARY" in ERM_AXIS_EXPIRES_WHEN and "expires" in ERM_AXIS_EXPIRES_WHEN
    assert ERM_DENSITY_AXIS["status"] == "TEMPORARY"

    # ⚠ AND THE TWO NUMBERS THAT MAY NOT BE RATIFIED. 231.789 is subdiv 7 over the measured shell and
    # 0.908 is subdiv 3; the first is inside the band by arithmetic and is refused on IDENTITY, the
    # second is below it. Reproduced here so the claim is checked rather than asserted in prose.
    area_7_5 = 4.0 * np.pi * 7.5 ** 2
    assert abs(163_842 / area_7_5 - 231.789) < 1e-3, "the subdiv-7 vertex density, to the digit"
    assert abs(642 / area_7_5 - 0.908) < 1e-3, "and the subdiv 3 the 2026-08-20 L0 runs ran"
    # ⚠ AND THE BAND DOES NOT EXCLUDE EITHER OF THEM. 231.789 lands inside it, and 0.908 lands inside
    # it too — just above the force-corrected capacity floor and 3.6x under the recorded gate minimum.
    # THAT is why the guard is structural rather than a range check: a band is a claim about
    # MAGNITUDE, and the ERM defect is a claim about IDENTITY. A range check would have admitted both.
    assert lo < 163_842 / area_7_5 < hi
    assert lo < 642 / area_7_5 < VG_ERM_CAPACITY_MIN_RECORDED_PER_UM2

    arena = WorldArena(capacity={Kind.NODE: 4_000, Kind.SEGMENT: 4_000,
                                 Kind.ANGLE3: 4_000, Kind.BOND: 4_000})
    # Two concentric shells 0.1 um apart — the arena's own membrane/cortex offset, h_cortex/2.
    rng = np.random.default_rng(0)
    u = rng.normal(size=(600, 3))
    u /= np.linalg.norm(u, axis=1, keepdims=True)
    mem_pos, cor_pos = 7.5 * u, 7.4 * u[::-1]

    mem = build_strand(arena, "membrane", start=(0, 0, 0), direction=(1, 0, 0),
                       contour_um=5.99, seg_um=0.01)
    cor = build_strand(arena, "cortex", start=(0, 1, 0), direction=(1, 0, 0),
                       contour_um=5.99, seg_um=0.01)
    assert mem.n_nodes == cor.n_nodes == 600, (mem.n_nodes, cor.n_nodes)
    area = 4.0 * np.pi * 7.5 ** 2
    fixture = dict(membrane_pos_um=mem_pos, membrane_lo=mem.nodes.lo, membrane_area_um2=area,
                   cortex_pos_um=cor_pos, cortex_lo=cor.nodes.lo, reach_um=0.5)

    # A caller who brings nothing is still refused, and the refusal now names the TWO that remain
    # rather than three — the queue reads what is actually open.
    try:
        build_membrane_erm_cortex(arena, **fixture)
    except ConnectorGapError as exc:
        assert exc.connector == CONNECTOR
        assert tuple(exc.missing) == _STILL_BLOCKED, exc.missing
        assert "TEMPORARY declared axis" in str(exc) and "erm_density_axis()" in str(exc)
    else:  # pragma: no cover
        raise AssertionError("the bundle's other two items still block a defaulted build")
    assert arena.n_live(Kind.BOND) == 0, "a refused family must claim nothing"

    # THE AXIS, PASSED EXPLICITLY — which is what decision 3 makes possible. The fixture's 600-vertex
    # shell cannot hold the test point over a real MCF7 area, so the axis is exercised at its LOW end
    # here and the band arithmetic is checked separately above.
    axis = erm_density_axis(lo)
    fam = build_membrane_erm_cortex(arena, count=axis, **fixture)
    assert fam.n_bonds == int(np.floor(lo * area + 0.5)) == 528

    # THE WORD TRAVELS. A driver does not have to remember to record it: it is in the row.
    row = fam.provenance_row()
    assert "TEMPORARY" in row["provenance"] and "expires" in row["provenance"]
    assert row["source_class"] == "UNRATIFIED_PROXY" and row["support_unit"] == "um^2"
    # CAPACITY, not occupancy — and the partner is static only because nothing detaches YET.
    assert str(row["n_bonds_is"]).startswith("CAPACITY")
    assert fam.partner_is_dynamic is False and row["live_partner_owner"] is None
    assert "CAPACITY, never occupancy" in str(ERM_DENSITY_AXIS["count_is"])
    assert row["support"] == area and row["density"] == lo

    # The mechanism behind the axis is unchanged and still correct.
    assert fam.component_pairs(arena) == {("cortex", "membrane"): 528}
    assert np.unique(fam.node_j).size == fam.n_bonds
    assert np.allclose(fam.rest_um, 0.1, atol=5e-3), (fam.rest_um.min(), fam.rest_um.max())
    assert np.all(fam.stiffness_pn_per_um == K_ERM_PN_PER_UM)

    # ⚠ THE MESH ARTEFACT, REFUSED BY IDENTITY RATHER THAN BY MAGNITUDE. One tether per vertex at the
    # fixture's own resolution is 600/um^2-equivalent here and 231.789 at subdiv 7 — the guard fires on
    # both because it compares the COUNT to the MESH, which is what the defect actually is.
    one_per_vertex = erm_density_axis(mem_pos.shape[0] / area)
    try:
        build_membrane_erm_cortex(arena, count=one_per_vertex, **fixture)
    except ValueError as exc:
        assert "ONE TETHER PER VERTEX" in str(exc) and "231.789" in str(exc)
    else:  # pragma: no cover
        raise AssertionError("a one-per-vertex density must be refused at any subdivision")

    # And the axis refuses to leave its own band — in both directions, and by name.
    for outside in (10.0 * hi, 0.5 * lo, 0.0):
        try:
            erm_density_axis(outside)
        except ValueError as exc:
            assert "outside the declared band" in str(exc) and "231.789" in str(exc)
        else:  # pragma: no cover
            raise AssertionError(f"{outside} is outside the band and must refuse")

    # A density low enough to give no tethers is an ANSWER; a reach that cannot be met is not.
    empty = BondCount(basis="areal", value=1e-9, scope="boundary case",
                      source_class=SourceClass.PI_GAP, provenance="deliberately sparse")
    assert build_membrane_erm_cortex(arena, count=empty, **fixture).n_bonds == 0
    try:
        erm_candidate_pairs(10, membrane_pos_um=mem_pos, membrane_lo=mem.nodes.lo,
                            cortex_pos_um=cor_pos, cortex_lo=cor.nodes.lo, reach_um=0.001)
    except ValueError as exc:
        assert "exceed the 0.001 um reach" in str(exc)
    else:  # pragma: no cover
        raise AssertionError("an unreachable pairing must raise rather than drop linkers")

    # More tethers than cortex nodes is a resolution failure, never a shared endpoint.
    try:
        erm_candidate_pairs(cor_pos.shape[0] + 1, membrane_pos_um=mem_pos, membrane_lo=mem.nodes.lo,
                            cortex_pos_um=cor_pos, cortex_lo=cor.nodes.lo, reach_um=0.5)
    except ValueError as exc:
        assert "reusing a cortex endpoint" in str(exc)
    else:  # pragma: no cover
        raise AssertionError("an over-subscribed cortex must raise")

    # A volumetric count is a category error for a surface population.
    try:
        build_membrane_erm_cortex(
            arena, count=BondCount(basis="volumetric", value=1.0, scope="x",
                                   source_class=SourceClass.PI_GAP,
                                   provenance="wrong basis on purpose"), **fixture)
    except ValueError as exc:
        assert "areal population" in str(exc)
    else:  # pragma: no cover
        raise AssertionError("a volumetric ERM density must be refused by name")

    print(f"{CONNECTOR} self-check OK — rho_ERM is a TEMPORARY declared axis, band "
          f"[{lo:.4f}, {hi:.0f}]/um^2, test point {ERM_DENSITY_TEST_POINT_PER_UM2:.3f}/um^2; "
          f"231.789 (subdiv 7) and 0.908 (subdiv 3) both REFUSED; still BLOCKED on "
          f"{len(_STILL_BLOCKED)}: {', '.join(_STILL_BLOCKED)}")


if __name__ == "__main__":
    _demo()
