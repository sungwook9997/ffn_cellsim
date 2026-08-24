"""The first kinetic bond family, and the citations it stands on, checked rather than only written.

Each module's own ``_demo()`` is the primary check (the ``world`` convention) and
``test_families_census.py`` already runs every family's. This file adds the two kinds of thing a
``_demo()`` cannot do:

1. **It checks the claims this family makes ABOUT OTHER FILES.** ``nmii_cortex_crossbridge`` asserts
   in its docstring that the declared chemistry card names two retired laws, that ``hand_kmc``'s NMIIA
   preset is the retired Bell one, that ``components/motor/hand.py`` keeps a pure-Bell detach kernel
   for other consumers, and that today's step path cannot accept anything. Every one of those is a
   citation, and a citation nobody checks is a comment. These tests fail when the cited file changes —
   which is the point: the failure message says to update the docstring, not to restore the defect.
2. **It puts the Warp codegen type-check in the collector.** The dev machine has no card, so the only
   thing standing between this module and Lead's GPU host is that the kernel compiles. That belongs in
   pytest rather than in a shell command somebody has to remember.

⚠ WHAT IS NOT HERE.  No number this family produces is checked against a physical band, because it
produces none. Existence, sign, refusal and the commit gate are the whole content — see the module
docstring's "no magnitude is claimed anywhere".
"""

from __future__ import annotations

import ast
import pathlib

import numpy as np
import pytest

from aleph.laws import crossbridge_kmc as law
from aleph.world.arena import Kind, WorldArena
from aleph.world.bond import BondCount, SourceClass, build_bond_family
from aleph.world.build.cortex import StrandPopulation
from aleph.world.families import ConnectorGapError
from aleph.world.families import nmii_cortex_crossbridge as xb

_PACKAGE = pathlib.Path(__file__).resolve().parents[2]


# ── the two self-checks, collected ──────────────────────────────────────────────────────────────
def test_law_demo_passes() -> None:
    """``laws/crossbridge_kmc``'s own self-check is the primary check for the rate law."""
    law._demo()


def test_family_demo_passes() -> None:
    """``families/nmii_cortex_crossbridge``'s own self-check is the primary check for the family."""
    xb._demo()


def test_the_kernel_type_checks_for_cuda() -> None:
    """The kernel emits CUDA C++ without a device. This is the whole dev-machine handoff gate.

    CPU execution is forbidden and there is no card here, so "it compiles" is the strongest statement
    this machine may make about a kernel — and it is a real one: codegen type-checks every expression.
    It produces no number, so it is not a CPU result under any reading.
    """
    assert law._codegen_check() > 0


# ── the physics claim: a SIGN, not a magnitude ──────────────────────────────────────────────────
def test_the_off_rate_falls_with_load_then_rises() -> None:
    """THE claim of the rate law: catch below ``F*``, slip above it — the opposite of Bell at low load.

    Checked on a PROVISIONAL parameter set, which is legitimate precisely because the property is a
    sign: it holds for every admissible set, so no magnitude is being asserted by exhibiting one.
    """
    cs = _provisional()
    f_star = cs.peak_force_pn
    f = np.linspace(0.0, 4.0 * f_star, 2001)
    p = np.asarray(cs.off_rate(f), np.float64)
    assert np.all(np.diff(p[f < 0.98 * f_star]) < 0.0), "catch branch must SLOW detachment under load"
    assert np.all(np.diff(p[f > 1.02 * f_star]) > 0.0), "slip branch must ACCELERATE it past F*"


def test_a_slip_dominated_set_is_refused_at_construction() -> None:
    """A parameter set with no interior ``F*`` is a Bell bond wearing a catch bond's name."""
    with pytest.raises(ValueError, match="SLIP-DOMINATED"):
        law.CrossbridgeKinetics(k_catch0=0.01, x_catch_um=1.0e-3, k_slip0=1.0, x_slip_um=0.6e-3,
                                k_on=50.0, capture_radius_um=0.210, scope="x", provenance="x")


@pytest.mark.parametrize("field", sorted(law.NMII_CROSSBRIDGE_GAPS))
def test_every_kinetic_magnitude_refuses_to_default(field: str) -> None:
    """Six PI-GAPs, six refusals, each naming what is unresolved about that one."""
    with pytest.raises(ValueError, match="PI-GAP"):
        law.CrossbridgeKinetics(**{**_provisional_kwargs(), field: None})


# ── the citations this family stands on ─────────────────────────────────────────────────────────
def test_the_declared_chemistry_card_still_names_two_retired_laws() -> None:
    """``engine/contracts.py`` gives all four NMII motor edges ``nmii_head_actin_hill_bell``.

    Hill was superseded 2026-07-07 (linear FV, PI-ratified) and re-scoped as muscle-only by PI
    decision 7 of 2026-08-21; Bell was superseded for this head by the 2026-07-23 catch-slip
    correction. If this test fails the disagreement has CLOSED — update the module docstring and
    :data:`DECLARED_CHEMISTRY_CARD`, do not restore the old card.
    """
    contracts = (_PACKAGE / "engine" / "contracts.py").read_text()
    assert xb.DECLARED_CHEMISTRY_CARD in contracts, (
        "the declared card is no longer in engine/contracts.py — the disagreement this family reports "
        "has changed shape; re-read it before quoting the docstring.")
    assert contracts.count(f'chemistry_card="{xb.DECLARED_CHEMISTRY_CARD}"') == 4, (
        "four NMII motor edges carried this card; the count changed. The duplicate verdict in "
        "DUPLICATE_OF_ENTRY names three siblings and is derived from that count.")


def test_the_named_nmii_preset_in_laws_is_still_the_retired_bell_law() -> None:
    """``hand_kmc.NMIIA_MYOSIN`` is ``catch_slip=False`` — the law the 2026-07-23 fix retired.

    Reported by this family rather than edited: it is another lane's file with other consumers. When
    it is corrected this test fails, and that failure is the signal to drop the warning from the
    module docstring.
    """
    from aleph.laws.hand_kmc import NMIIA_MYOSIN

    assert NMIIA_MYOSIN.catch_slip is False, (
        "hand_kmc's NMIIA preset became catch-slip. Good — remove the warning from "
        "laws/crossbridge_kmc.py's docstring, which reports it as still-Bell.")


def test_the_pure_bell_detach_kernel_is_retained_for_its_other_consumers() -> None:
    """``components/motor/hand.py`` keeps BOTH kernels, and this bond's law is the catch-slip one."""
    hand = (_PACKAGE / "components" / "motor" / "hand.py").read_text()
    assert "def step_detach_kernel(" in hand, "the pure Bell kernel is what ERM/alpha-actinin use"
    assert "def step_detach_catch_slip_kernel(" in hand, "the NMII head's kernel"
    assert "Kovacs" in hand and "2026-07-23" in hand, "the fidelity correction this family follows"


def test_the_step_path_cannot_accept_anything_today() -> None:
    """The census claims ``commit_reachable_today: False``. This is that claim, checked.

    ``world/step.py`` defaults to :class:`UndefinedAcceptance`, whose verdict is
    ``ACCEPTANCE_UNDEFINED`` and whose ``is_accepted`` is False by design (PI decision 18). So nothing
    this family proposes can be applied on the current path — and when a real predicate lands, this
    test fails and the census row has to be re-derived rather than left stale.
    """
    from aleph.world.step import ACCEPTANCE_UNDEFINED, StepVerdict, UndefinedAcceptance

    verdict = StepVerdict(
        step_index=0, verdict=UndefinedAcceptance().decide(None), predicate_name="undefined",
        predicate_tests="nothing", dt_phys_s=0.0, residual_pn=None, tol_pn=None, n_events=0,
        wall_s=0.0)
    assert verdict.verdict == ACCEPTANCE_UNDEFINED
    assert verdict.is_accepted is False, "an unjudged step must never read as accepted"


# ── the refusals that are this module's product ─────────────────────────────────────────────────
def test_a_uniform_polarity_cortex_refuses_by_type_not_by_attribute_error() -> None:
    """A pending build configuration must be distinguishable from a broken module, by type."""
    arena, inventory, mf_pos, cor_pos, cortex_of = _fixture()
    with pytest.raises(ConnectorGapError) as excinfo:
        xb.plan_crossbridge_stations(
            nmii_inventory=inventory, nmii_pos_um=mf_pos, cortex=cortex_of(None),
            cortex_pos_um=cor_pos, reach_um=0.3)
    err = excinfo.value
    assert not isinstance(err, AttributeError)
    assert err.missing == ["cortex.polarity_per_strand"]
    assert "PI decision 4" in str(err)
    assert arena.n_live(Kind.BOND) == 0, "a refused family claims nothing"


def test_both_head_rows_on_one_filament_is_refused() -> None:
    """A bipolar minifilament straddling ONE filament pulls it against itself."""
    _, _, _, _, cortex_of = _fixture()
    cortex = cortex_of(np.array([+1, -1] * 3, np.int8))
    with pytest.raises(ValueError, match="against itself"):
        xb.assert_bipolar_station(cortex, 2, 2)


def test_two_same_polarity_filaments_are_refused() -> None:
    """Two filaments running the same way TRANSLATE the pair instead of contracting it."""
    _, _, _, _, cortex_of = _fixture()
    cortex = cortex_of(np.array([+1, -1] * 3, np.int8))
    with pytest.raises(ValueError, match="TRANSLATES the pair"):
        xb.assert_bipolar_station(cortex, 0, 2)
    xb.assert_bipolar_station(cortex, 0, 1)  # the only legal arrangement


def test_the_head_count_axis_has_a_test_point_and_refuses_a_band() -> None:
    """PI decision 6 ratified the axis and the test point 30. The band was NOT ratified."""
    assert xb.N_HEADS_PER_SIDE_TEST_POINT == 30
    with pytest.raises(TypeError, match="NOT ratified"):
        xb.crossbridge_count([10, 30], scope="x")
    count = xb.crossbridge_count(30, scope="MCF7 interphase adherent 37C")
    assert "NO BAND" in count.provenance
    assert count.resolve(442.0) == 26_520


def test_capacity_is_the_count_and_occupancy_never_is() -> None:
    """The kinetic split: capacity is structural and declared; occupancy is emergent and measured."""
    arena, inventory, mf_pos, cor_pos, cortex_of = _fixture()
    cortex = cortex_of(np.array([+1, -1] * 3, np.int8))
    h = inventory["n_heads_per_side"]
    family, state, census = xb.build_nmii_cortex_crossbridge(
        arena, nmii_inventory=inventory, nmii_pos_um=mf_pos, cortex=cortex, cortex_pos_um=cor_pos,
        reach_um=0.3, count=xb.crossbridge_count(h, scope="fixture", provenance="fixture"),
        k_xb_pn_per_um=1.0, r0_xb_um=0.0, kinetics=_provisional())

    assert family.n_bonds == inventory["n_minifilaments"] * 2 * h, "one site per head, exactly"
    assert np.unique(family.node_i).size == family.n_bonds, "a head carries ONE molecular state"
    assert state.bound_fraction == 0.0, "a bound t0 would be an imposed duty ratio"
    assert census["commit_reachable_today"] is False
    assert census["partner_is_dynamic"] is True, "declared through bond.py's own field, f99400a6"
    assert census["live_partner_owner"].endswith("CrossbridgeState.bound")
    assert family.partner_is_dynamic and family.live_partner_owner
    assert family.provenance_row()["n_bonds_is"].startswith("CAPACITY")

    # A count that is not the head count is a disagreement, not a rounding.
    with pytest.raises(ValueError, match="One site per head"):
        xb.build_nmii_cortex_crossbridge(
            arena, nmii_inventory=inventory, nmii_pos_um=mf_pos, cortex=cortex,
            cortex_pos_um=cor_pos, reach_um=0.3,
            count=xb.crossbridge_count(h + 1, scope="x", provenance="x"),
            k_xb_pn_per_um=1.0, r0_xb_um=0.0, kinetics=_provisional())

    # An areal density is a category error for a per-minifilament population.
    with pytest.raises(ValueError, match="PER MINIFILAMENT"):
        xb.build_nmii_cortex_crossbridge(
            arena, nmii_inventory=inventory, nmii_pos_um=mf_pos, cortex=cortex,
            cortex_pos_um=cor_pos, reach_um=0.3, k_xb_pn_per_um=1.0, r0_xb_um=0.0,
            kinetics=_provisional(),
            count=BondCount(basis="areal", value=1.0, scope="x", source_class=SourceClass.PI_GAP,
                            provenance="wrong basis on purpose"))


def test_a_kinetic_transition_commits_only_on_an_accepted_step() -> None:
    """The transaction, made structural: a refused commit changes nothing and says why."""
    state = xb.CrossbridgeState.all_free(8)
    proposal = np.ones(8, np.int8)
    with pytest.raises(RuntimeError, match="not accepted"):
        state.commit(proposal, accepted=False)
    assert state.bound_fraction == 0.0 and state.n_commits == 0, "a refusal must change nothing"
    assert state.commit(proposal, accepted=True) == 8
    assert state.bound_fraction == 1.0 and state.n_commits == 1


# ── fixtures ────────────────────────────────────────────────────────────────────────────────────
def _provisional_kwargs() -> dict:
    """The PROVISIONAL catch-slip set, labelled as such wherever it appears. NOT a run configuration."""
    return dict(
        k_catch0=0.35, x_catch_um=1.0e-3, k_slip0=0.35, x_slip_um=0.6e-3, k_on=50.0,
        capture_radius_um=0.210, scope="PROVISIONAL — test fixture, NOT a run configuration",
        provenance="ac_gate_b_cortex_motor_native.CATCH_SLIP_PROVISIONAL; every value an I0-B3 PI-GAP")


def _provisional() -> law.CrossbridgeKinetics:
    return law.CrossbridgeKinetics(**_provisional_kwargs())


def _fixture():
    """Six cortex filaments on two parallel lines and two minifilaments straddling them.

    Host geometry only — no device is touched, which is what makes this collectable on the dev
    machine at all. The cortex population is hand-constructed rather than built, because
    ``build_strand_population`` is a CUDA kernel and CPU execution is forbidden.
    """
    n_str, n_per_str, h, n_bb, n_mf = 6, 5, 3, 4, 2
    arena = WorldArena(capacity={Kind.NODE: 400, Kind.SEGMENT: 400, Kind.ANGLE3: 400, Kind.BOND: 400})
    cnodes = arena.claim("cortex", Kind.NODE, n_str * n_per_str)
    cor_pos = np.zeros((n_str * n_per_str, 3))
    for s in range(n_str):
        sl = slice(s * n_per_str, (s + 1) * n_per_str)
        cor_pos[sl, 0] = np.arange(n_per_str) * 0.1
        cor_pos[sl, 1] = 0.0 if s % 2 == 0 else 0.4
        cor_pos[sl, 2] = (s // 2) * 1.0

    n_per = n_bb + 2 * h
    nmii_nodes = arena.claim("nmii", Kind.NODE, n_mf * n_per)
    mf_pos = np.zeros((n_mf * n_per, 3))
    for m in range(n_mf):
        b = m * n_per
        mf_pos[b:b + n_per, 0] = 0.2
        mf_pos[b:b + n_per, 2] = m * 1.0
        mf_pos[b:b + n_bb, 1] = 0.2
        mf_pos[b + n_bb:b + n_bb + h, 1] = 0.05
        mf_pos[b + n_bb + h:b + n_per, 1] = 0.35

    inventory = {"n_minifilaments": n_mf, "n_bb": n_bb, "n_heads_per_side": h,
                 "nodes_per_minifilament": n_per,
                 "claims": {"node": (nmii_nodes.lo, nmii_nodes.hi)}}

    def cortex_of(polarity_per_strand):
        return StrandPopulation(
            population="cortex", nodes=cnodes, segments=cnodes, angles=cnodes,
            n_strands=n_str, nodes_per_strand=n_per_str,
            seg_node=None, seg_rest_um=None, seg_arc_um=None, angle_idx=None,
            radius_um=7.4, thickness_um=0.2, centre_um=(0.0, 0.0, 0.0),
            areal_density_um2=100.0, density_provenance="test fixture",
            contour_um_requested=0.4, contour_um_realised=0.4, seg_um_requested=0.1,
            seg_um_realised=0.1, seed=0, polarity=+1, topology_bytes=0,
            polarity_per_strand=polarity_per_strand,
            polarity_basis="test fixture; alternating by construction")

    return arena, inventory, mf_pos, cor_pos, cortex_of


# ── the dynamic-partner declaration, and the merge argument, re-derived ─────────────────────────
def test_the_pairing_validator_fires_on_the_path_this_family_actually_takes() -> None:
    """A half-declared dynamic partner is refused THROUGH THE BUILDER, not only on the dataclass.

    ⚠ Kept after its sibling ratchet was deleted, and deliberately: the family's route to these two
    fields changed once already (``dataclasses.replace`` until ``f99400a6`` forwarded them through
    ``build_bond_family``), and a validator is only worth having on the path in use. This test
    follows the path rather than the dataclass, so the next time the route moves it moves with it.

    The failure it prevents is not hypothetical — ``bond.py`` cites it: a family that declares a
    dynamic partner and names no owner lets ``node_j`` read as the partner in force at that step,
    which is how ``components/motor/INTEGRATION.md:142`` records the first I3 source realising a
    passive spring network while the active contraction was never generated.
    """
    arena, inventory, mf_pos, cor_pos, cortex_of = _fixture()
    common = dict(chemistry_card=xb.CHEMISTRY_CARD, support=1.0, pairs=np.array([[0, 1]]),
                  rest_um=0.0, stiffness_pn_per_um=1.0,
                  count=BondCount(basis="explicit", value=1.0, scope="fixture",
                                  source_class=SourceClass.PI_GAP, provenance="fixture"))

    with pytest.raises(ValueError, match="names no owner"):
        build_bond_family(arena, "half_declared_a", partner_is_dynamic=True,
                          live_partner_owner="", **common)
    with pytest.raises(ValueError, match="does not declare the partner dynamic"):
        build_bond_family(arena, "half_declared_b", partner_is_dynamic=False,
                          live_partner_owner="somewhere", **common)

    # And the family itself goes through that same builder, fully declared.
    family, _, census = xb.build_nmii_cortex_crossbridge(
        arena, nmii_inventory=inventory, nmii_pos_um=mf_pos,
        cortex=cortex_of(np.array([+1, -1] * 3, np.int8)), cortex_pos_um=cor_pos, reach_um=0.3,
        count=xb.crossbridge_count(inventory["n_heads_per_side"], scope="fixture",
                                   provenance="fixture"),
        k_xb_pn_per_um=1.0, r0_xb_um=0.0, kinetics=_provisional())
    assert family.partner_is_dynamic and family.live_partner_owner.endswith(
        "nmii_cortex_crossbridge.CrossbridgeState.bound")
    assert census["live_partner_owner"] == family.live_partner_owner


def _nmii_motor_contracts() -> list[ast.Call]:
    """Every ``ConnectorContract`` call in ``engine/contracts.py`` carrying this family's card.

    Parsed, not grepped, and NOT truncated — the merge argument is a claim about ALL of them, and
    Lead's own warning of 2026-08-21 is that a ``grep | head`` reported as exhaustive produced two
    defects in one hour.
    """
    tree = ast.parse((_PACKAGE / "engine" / "contracts.py").read_text())
    out = []
    for node in ast.walk(tree):
        if not (isinstance(node, ast.Call) and getattr(node.func, "id", "") == "ConnectorContract"):
            continue
        card = next((k.value.value for k in node.keywords
                     if k.arg == "chemistry_card" and isinstance(k.value, ast.Constant)), None)
        if card == xb.DECLARED_CHEMISTRY_CARD:
            out.append(node)
    return out


def test_the_four_nmii_motor_edges_differ_in_exactly_one_declared_field() -> None:
    """The merge argument, re-derived from the source: same everything, different target population.

    This is what makes :data:`DUPLICATE_OF_ENTRY` a physics verdict rather than a name verdict. If it
    fails, the four have diverged and the merge has to be re-argued before Lead applies it.
    """
    calls = _nmii_motor_contracts()
    assert len(calls) == 4, f"expected four NMII motor edges, parsed {len(calls)}"

    names, targets, tails = [], [], set()
    for call in calls:
        args = [a.value for a in call.args if isinstance(a, ast.Constant)]
        # positional: name, <ConnectorFamily.MOTOR>, source, target, kinetics, commit_on_accept
        family_attr = ast.unparse(call.args[1])
        names.append(args[0])
        targets.append(args[2])
        kw = {k.arg: ast.unparse(k.value) for k in call.keywords}
        tails.add((family_attr, args[1], args[3], args[4],
                   kw.get("endpoint_role_a"), kw.get("endpoint_role_b"),
                   kw.get("generation_required"), kw.get("remap_on_accept"),
                   kw.get("blocks_sleep_refine")))

    assert len(tails) == 1, f"the four edges are NOT identical apart from the target: {tails}"
    only = next(iter(tails))
    assert only[0] == "ConnectorFamily.MOTOR" and only[1] == "nmii"
    assert only[2] is True and only[3] is True, "all four are kinetic and commit on accept"
    assert len(set(targets)) == 4, f"four distinct target populations, got {targets}"
    assert set(names) == set(xb.DUPLICATE_OF_ENTRY) | {xb.CONNECTOR}, (
        "DUPLICATE_OF_ENTRY must name exactly the three siblings of this canonical side")


def test_the_control_a_motor_edge_that_is_not_this_card() -> None:
    """Folding by CARD is not folding by family — the control that makes the verdict trustworthy.

    ``mt_cortex_capture`` is ``ConnectorFamily.MOTOR`` too and carries a different card, so the merge
    separates a myosin crossbridge from a dynein capture site rather than absorbing it. Same shape as
    D13's ``filopodium_cortex_root`` control.
    """
    src = (_PACKAGE / "engine" / "contracts.py").read_text()
    assert src.count("ConnectorFamily.MOTOR") == 5, "five MOTOR edges; the count is the argument"
    assert len(_nmii_motor_contracts()) == 4, "four of them are this card"
    assert 'chemistry_card="microtubule_cortical_dynein"' in src, (
        "the fifth MOTOR edge's distinct card is the control; if it changed, re-argue the merge")


def test_the_contract_already_declared_the_partner_live() -> None:
    """Independent corroboration of the ``bond.py`` finding, from before this family existed.

    All four declare ``endpoint_role_b="live polar actin material coordinate"`` — *live*, and a
    *coordinate* rather than a node. The dynamic partner was in the contract all along and invisible
    in ``BondFamily`` until ``f99400a6``.
    """
    for call in _nmii_motor_contracts():
        role_b = next(k.value.value for k in call.keywords if k.arg == "endpoint_role_b")
        assert "live" in role_b and "coordinate" in role_b, role_b


# ── the pre-registration: written before the values, and checkable without them ─────────────────
def test_there_is_no_partial_constant_set_that_runs() -> None:
    """The correction this write-up produced: not five, not six — eight, or the geometry census.

    Lane A told Lead the runnable minimum was five (``k_xb`` + the four catch-slip). Tension needs
    ``r0_xb`` too, which is six, and ``CrossbridgeKinetics`` refuses without ``k_on`` and
    ``capture_radius_um``, which is eight. This test is the arithmetic, so the claim in the PI queue
    cannot drift back to a smaller number.
    """
    six = dict(k_catch0=0.35, x_catch_um=1.0e-3, k_slip0=0.35, x_slip_um=0.6e-3,
               scope="x", provenance="x")
    for absent in ("k_on", "capture_radius_um"):
        supplied = {**six, "k_on": 50.0, "capture_radius_um": 0.210, absent: None}
        with pytest.raises(ValueError, match="PI-GAP"):
            law.CrossbridgeKinetics(**supplied)

    assert len(xb.BLOCKED_BY) == 8
    assert set(law.NMII_CROSSBRIDGE_GAPS) | {"k_xb_pn_per_um", "r0_xb_um"} == set(xb.BLOCKED_BY)


def test_the_detach_only_tier_is_refused_by_the_all_free_initial_state() -> None:
    """Why there is no six-constant tier, checked rather than only argued.

    A detach-only measurement needs an initial bound set. ``CrossbridgeState`` is born all-free
    because a bound ``t0`` is an imposed duty ratio, and attachment is the only non-imposed source of
    one — so skipping ``k_on`` does not buy a cheaper experiment, it buys one whose initial condition
    is what ``params_i0b3.yaml`` forbids imposing.
    """
    state = xb.CrossbridgeState.all_free(16)
    assert state.bound_fraction == 0.0, "the only legal initial state"
    assert "tier_0_station_census" in xb.MEASUREMENT_TIERS
    assert "tier_1_attach_detach_cycle" in xb.MEASUREMENT_TIERS
    assert len(xb.MEASUREMENT_TIERS) == 2, "the tiers are zero and eight; there is nothing between"


def test_the_tiers_name_their_constants_and_tier_zero_names_none() -> None:
    """A tier that could not be enumerated is not a queue entry, only a paragraph."""
    assert xb.MEASUREMENT_TIERS["tier_0_station_census"]["requires"] == ()
    assert xb.MEASUREMENT_TIERS["tier_1_attach_detach_cycle"]["requires"] == xb.BLOCKED_BY
    for tier in xb.MEASUREMENT_TIERS.values():
        assert str(tier["buys"]).strip() and str(tier["blocked_on"]).strip()


def test_tier_zero_is_reachable_with_no_constant_at_all() -> None:
    """The claim that the station census needs none of the eight, exercised.

    ``plan_crossbridge_stations`` takes positions, a polarity-carrying population and a declared
    reach — and no stiffness, no rate, no force. If this ever starts needing one, the tier-0 entry in
    the PI queue is wrong and this test is what says so.
    """
    from inspect import signature

    _, inventory, mf_pos, cor_pos, cortex_of = _fixture()
    for fn in (xb.station_census, xb.plan_crossbridge_stations):
        params = set(signature(fn).parameters)
        assert params.isdisjoint(set(xb.BLOCKED_BY)), (
            f"{fn.__name__} acquired a PI-GAP argument: {params & set(xb.BLOCKED_BY)}")

    plan = xb.plan_crossbridge_stations(
        nmii_inventory=inventory, nmii_pos_um=mf_pos,
        cortex=cortex_of(np.array([+1, -1] * 3, np.int8)), cortex_pos_um=cor_pos, reach_um=0.3)
    assert plan["n_stations"] == inventory["n_minifilaments"]


def test_what_no_tier_buys_is_recorded_beside_what_they_do() -> None:
    """The half that decides what may be quoted must travel with the half that sells the ruling."""
    joined = " ".join(xb.NOT_BOUGHT_BY_ANY_TIER)
    assert "gamma" in joined and "does not STEP" in joined, "no tension: this family never steps"
    assert "PI decision 18" in joined, "no magnitude while acceptance is undefined"
    assert "86x" in joined, "the tolerance is loosened by the events under test"
    assert "TEMPORARY" in joined, "a declared axis is not a sourced number"


# ── TIER 0: the census reports where the builder refuses ────────────────────────────────────────
def _census_fixture():
    """Three minifilaments, one per outcome: stationed, no cortex in reach, no anti-parallel partner.

    Built so the two FAILURE reasons are distinguishable, because they have different owners — a
    placement problem belongs to ``build/nmii.py`` and a polarity problem to ``build/cortex.py``.
    """
    n_per_str, h, n_bb, n_mf = 5, 2, 4, 3
    n_per = n_bb + 2 * h
    arena = WorldArena(capacity={Kind.NODE: 400, Kind.SEGMENT: 400, Kind.ANGLE3: 400, Kind.BOND: 400})

    # Four filaments: a +/- pair at z=0 (a legal station), and a +/+ pair at z=10 (no partner).
    pol = np.array([+1, -1, +1, +1], np.int8)
    zs = [0.0, 0.0, 10.0, 10.0]
    cnodes = arena.claim("cortex", Kind.NODE, 4 * n_per_str)
    cor = np.zeros((4 * n_per_str, 3))
    for s_i in range(4):
        sl = slice(s_i * n_per_str, (s_i + 1) * n_per_str)
        cor[sl, 0] = np.arange(n_per_str) * 0.1
        cor[sl, 1] = 0.0 if s_i % 2 == 0 else 0.4
        cor[sl, 2] = zs[s_i]

    nmii_nodes = arena.claim("nmii", Kind.NODE, n_mf * n_per)
    mf = np.zeros((n_mf * n_per, 3))
    for m, z in enumerate((0.0, 500.0, 10.0)):   # stationed / far away / same-polarity neighbourhood
        b = m * n_per
        mf[b:b + n_per, 0] = 0.2
        mf[b:b + n_per, 2] = z
        mf[b:b + n_bb, 1] = 0.2
        mf[b + n_bb:b + n_bb + h, 1] = 0.05
        mf[b + n_bb + h:b + n_per, 1] = 0.35

    inventory = {"n_minifilaments": n_mf, "n_bb": n_bb, "n_heads_per_side": h,
                 "nodes_per_minifilament": n_per,
                 "claims": {"node": (nmii_nodes.lo, nmii_nodes.hi)}}
    cortex = StrandPopulation(
        population="cortex", nodes=cnodes, segments=cnodes, angles=cnodes,
        n_strands=4, nodes_per_strand=n_per_str,
        seg_node=None, seg_rest_um=None, seg_arc_um=None, angle_idx=None,
        radius_um=7.4, thickness_um=0.2, centre_um=(0.0, 0.0, 0.0),
        areal_density_um2=100.0, density_provenance="test fixture",
        contour_um_requested=0.4, contour_um_realised=0.4, seg_um_requested=0.1,
        seg_um_realised=0.1, seed=0, polarity=+1, topology_bytes=0,
        polarity_per_strand=pol, polarity_basis="test fixture")
    return arena, inventory, mf, cor, cortex


def test_the_census_reports_where_the_builder_refuses() -> None:
    """Same arrangement, opposite obligations: the builder must raise, the census must count."""
    _, inventory, mf, cor, cortex = _census_fixture()
    kw = dict(nmii_inventory=inventory, nmii_pos_um=mf, cortex=cortex, cortex_pos_um=cor,
              reach_um=0.3)

    census = xb.station_census(**kw)
    assert census["n_minifilaments"] == 3
    assert census["n_stationed"] == 1 and census["n_unstationed"] == 2

    with pytest.raises(ValueError, match="NO legal bipolar station"):
        xb.plan_crossbridge_stations(**kw)


def test_the_two_failure_reasons_are_counted_separately_because_they_have_different_owners() -> None:
    """A placement problem is build/nmii.py's; a polarity problem is build/cortex.py's.

    Summing them into one "unstationed" number would send the fix to the wrong builder — which is
    why the census returns both and never a total.
    """
    _, inventory, mf, cor, cortex = _census_fixture()
    census = xb.station_census(nmii_inventory=inventory, nmii_pos_um=mf, cortex=cortex,
                               cortex_pos_um=cor, reach_um=0.3)
    assert census["n_no_cortex_in_reach"] == 1, "the minifilament placed 500 um away"
    assert census["n_no_antiparallel_partner"] == 1, "the one among two +1 filaments"
    assert census["n_no_cortex_in_reach"] + census["n_no_antiparallel_partner"] == \
        census["n_unstationed"]


def test_an_unstationed_minifilament_is_minus_one_not_filament_zero() -> None:
    """A missing station must not be readable as a station on filament 0."""
    _, inventory, mf, cor, cortex = _census_fixture()
    census = xb.station_census(nmii_inventory=inventory, nmii_pos_um=mf, cortex=cortex,
                               cortex_pos_um=cor, reach_um=0.3)
    for m in census["unstationed"]:
        assert tuple(census["station"][m]) == (-1, -1)


def test_the_count_is_monotone_in_reach_which_is_why_one_reach_is_not_a_result() -> None:
    """Reach is a declared axis. A single point hides it, and picking the flattering point is banned.

    Monotonicity is what makes that concrete: the census can be driven to any count between its
    floor and its ceiling by the choice of one argument, so the CURVE is the result.
    """
    _, inventory, mf, cor, cortex = _census_fixture()
    counts = [xb.station_census(nmii_inventory=inventory, nmii_pos_um=mf, cortex=cortex,
                                cortex_pos_um=cor, reach_um=r)["n_stationed"]
              for r in (0.05, 0.3, 5.0, 1000.0)]
    assert counts == sorted(counts), f"stationed count must not fall as reach grows: {counts}"
    assert counts[0] < counts[-1], "the axis must actually move the answer, or it is not an axis"


def test_the_census_claims_no_magnitude_and_uses_no_gap_constant() -> None:
    """The two fields a reader checks before quoting anything out of this record."""
    _, inventory, mf, cor, cortex = _census_fixture()
    census = xb.station_census(nmii_inventory=inventory, nmii_pos_um=mf, cortex=cortex,
                               cortex_pos_um=cor, reach_um=0.3)
    assert census["no_pi_gap_constant_used"] is True
    assert census["magnitude_claimed"] is None
    assert census["cortex_polarity"]["n_plus"] + census["cortex_polarity"]["n_minus"] == 4


def test_the_census_reports_how_far_the_actin_is_without_being_asked() -> None:
    """One query answers what the first sweep needed nine reaches to locate.

    ``n_no_cortex_in_reach`` says the actin is further away than THIS reach. It does not say how far,
    so the 2026-08-21 tier-0 run had to recover the distance from where the count moved along a
    ladder. The distance is independent of ``reach_um`` and the instrument already holds it.
    """
    _, inventory, mf, cor, cortex = _census_fixture()
    near = [xb.station_census(nmii_inventory=inventory, nmii_pos_um=mf, cortex=cortex,
                              cortex_pos_um=cor, reach_um=r)["nearest_cortex_node_um"]
            for r in (0.05, 1000.0)]
    assert near[0] == near[1], "the distance to the actin cannot depend on the reach we asked about"
    assert near[0]["min"] <= near[0]["median"] <= near[0]["max"]
    # The fixture puts one minifilament 500 um away, so the spread must survive into the record
    # rather than being collapsed to a mean that hides it.
    assert near[0]["max"] > 100.0 and near[0]["min"] < 1.0


def test_a_population_that_can_never_attach_is_refused_not_built() -> None:
    """The 2026-08-21 TIER 0 defect, made unrepresentable.

    NMII shell [6.850, 7.050] and cortex shell [7.300, 7.500] are disjoint by 0.250 um against a
    capture-radius proxy of 0.210 um. The build succeeded and reported ``n_outside: 0``, and this
    family would have produced 26,520 rows that read as a motor. It now refuses, using the caller's
    OWN declared capture radius — no threshold is invented here — and testing the same condition the
    kernel's attach branch tests.
    """
    arena, inventory, mf, cor, cortex = _census_fixture()
    # Keep only the minifilament that has a legal station, then push the cortex radially away by
    # more than the declared capture radius — the shape of the real defect, at fixture scale.
    inv = {**inventory, "n_minifilaments": 1,
           "claims": {"node": (inventory["claims"]["node"][0],
                               inventory["claims"]["node"][0] + inventory["nodes_per_minifilament"])}}
    mf1 = mf[:inventory["nodes_per_minifilament"]]
    kw = dict(nmii_inventory=inv, nmii_pos_um=mf1, cortex=cortex, cortex_pos_um=cor,
              count=xb.crossbridge_count(inv["n_heads_per_side"], scope="fixture",
                                         provenance="fixture"),
              k_xb_pn_per_um=1.0, r0_xb_um=0.0)

    # In reach, inside capture: builds, and says so.
    near = xb.build_nmii_cortex_crossbridge(arena, reach_um=0.3, kinetics=_provisional(), **kw)
    assert near[2]["n_beyond_capture_at_t0"] == 0

    # Same geometry, a capture radius smaller than every built separation: refused.
    tiny = law.CrossbridgeKinetics(**{**_provisional_kwargs(), "capture_radius_um": 1e-4})
    with pytest.raises(ValueError, match="beyond the declared capture radius"):
        xb.build_nmii_cortex_crossbridge(arena, reach_um=0.3, kinetics=tiny, **kw)


def test_the_refusal_sends_the_reader_to_placement_not_to_rate_constants() -> None:
    """The message has to name the right owner, or it produces a search in the wrong module."""
    arena, inventory, mf, cor, cortex = _census_fixture()
    inv = {**inventory, "n_minifilaments": 1,
           "claims": {"node": (inventory["claims"]["node"][0],
                               inventory["claims"]["node"][0] + inventory["nodes_per_minifilament"])}}
    tiny = law.CrossbridgeKinetics(**{**_provisional_kwargs(), "capture_radius_um": 1e-4})
    with pytest.raises(ValueError) as excinfo:
        xb.build_nmii_cortex_crossbridge(
            arena, nmii_inventory=inv, nmii_pos_um=mf[:inventory["nodes_per_minifilament"]],
            cortex=cortex, cortex_pos_um=cor, reach_um=0.3, kinetics=tiny,
            count=xb.crossbridge_count(inv["n_heads_per_side"], scope="f", provenance="f"),
            k_xb_pn_per_um=1.0, r0_xb_um=0.0)
    msg = str(excinfo.value)
    assert "PLACEMENT result, not a kinetics one" in msg
    assert "inert for EVERY value of the eight PI-GAPs" in msg
    assert "station_census" in msg, "point at the tool that answers it without building anything"
