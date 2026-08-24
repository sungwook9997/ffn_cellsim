"""Host-side contracts of the native-spectrum driver, exercised without a device.

The driver's physics runs only on CUDA, but three of its decisions are pure host logic and each one
has already been a defect class in this repo:

* a family is "covered" because its kernel LAUNCHED — the WCA steric kernel launches over every node
  and contributes exactly nothing when no pair is inside its cutoff, which is the designed state of an
  ``overlap_free_cortex`` build;
* an attribute toggled for an ablation is not restored, so every later probe silently measures a
  different operator;
* a stability verdict quoted without the step it belongs to.

These are testable here, and they are what would make a green native run wrong rather than absent.
"""

from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pytest

from aleph.scripts import ac_observe_native_spectrum as driver


def _fake_cell(**overrides: object) -> SimpleNamespace:
    """Return a cell-shaped stand-in carrying only what the census and verdict read."""
    myosin = SimpleNamespace(
        backbone_bonds=SimpleNamespace(shape=(12,)),
        head_bonds=SimpleNamespace(shape=(20,)),
        backbone_angles=SimpleNamespace(shape=(8,)),
        head_arm_angles=SimpleNamespace(shape=(20,)),
        k_theta_backbone=1.0,
        k_theta_arm=0.0,                     # gated OFF: must be counted as zero, not as 20
        state={"bound": SimpleNamespace(shape=(5,))},
    )
    base = {
        "n_total": 100, "n_actin": 80, "n_fibers": 10, "n_tri": 30, "n_xl": 40,
        "steric": SimpleNamespace(n=100),
        "myosin": myosin,
        # `n_verts` is what a resized compartment shows up as; the real MembraneCompartment and
        # NucleusCompartment both carry it, and the census now records it so a run at a reduced
        # subdivision cannot read as a full-population one.
        "nucleus": SimpleNamespace(n_linc=7, n_verts=42),
        "membrane": SimpleNamespace(n_erm=9, n_faces=60, n_verts=32),
        "kmax": 2.0e6, "mechanical_kmax": 1.5e6, "dt_mu": 5.0e-8,
        "device": "device-under-test",
    }
    base.update(overrides)
    return SimpleNamespace(**base)


def _fake_solver(n: int = 100, *, mem_edges: int | None = 90) -> SimpleNamespace:
    edges = None if mem_edges is None else SimpleNamespace(shape=(mem_edges,))
    return SimpleNamespace(n=n, mem_edges_d=edges)


def test_the_launch_census_counts_a_gated_off_family_as_zero() -> None:
    """``k_theta_arm = 0`` means the kernel is not launched; counting its array size would over-state."""
    census = driver._stiffness_launch_census(_fake_cell(), _fake_solver())
    assert census["myosin_backbone_angle"] == 8
    assert census["myosin_arm_angle"] == 0
    assert census["membrane_area_edges"] == 90
    assert census["steric_wca"] == 100

    absent = driver._stiffness_launch_census(
        _fake_cell(steric=None, myosin=None, nucleus=None, membrane=None),
        _fake_solver(mem_edges=None))
    assert absent["steric_wca"] == 0 and absent["myosin_crossbridge"] == 0
    assert absent["membrane_erm"] == 0 and absent["membrane_area_edges"] == 0


def test_every_declared_family_names_an_attribute_that_actually_gates_a_launch() -> None:
    """The ablation is only exact if each toggle is the attribute the operator reads at launch time."""
    cell, solver = _fake_cell(), _fake_solver()
    for name, owner, attribute, _inert in driver._FAMILY_TOGGLES:
        target = cell if owner == "cell" else solver
        assert hasattr(target, attribute), f"{name} toggles {owner}.{attribute}, which does not exist"


def test_disabling_a_family_restores_the_attribute_even_when_the_probe_raises() -> None:
    """An un-restored toggle would make every LATER family's difference wrong, silently."""
    cell, solver = _fake_cell(), _fake_solver()
    original = cell.steric
    with driver._family_disabled(cell, solver, "cell", "steric", None):
        assert cell.steric is None
    assert cell.steric is original

    with pytest.raises(RuntimeError):
        with driver._family_disabled(cell, solver, "solver", "mem_edges_d", None):
            assert solver.mem_edges_d is None
            raise RuntimeError("probe blew up mid-ablation")
    assert solver.mem_edges_d is not None


def test_the_stability_verdict_is_the_product_the_integrator_is_actually_bounded_by() -> None:
    """``pos += dt_mu·P F`` is forward Euler; it is stable exactly while ``dt·λ_max < 2``."""
    cell = _fake_cell(dt_mu=5.0e-8, kmax=2.0e6)
    verdict = driver._explicit_step_verdict(cell, 1.0e7)

    assert verdict["stability_product_base"] == pytest.approx(0.5)
    assert verdict["base_step_stable"] is True
    # FIRE multiplies the step by up to 10; the product is reported, and deliberately NOT adjudicated
    # against this limit, because FIRE is damped MD and its bound scales differently
    assert verdict["stability_product_at_fire_max"] == pytest.approx(5.0)
    assert "NOT ADJUDICATED" in verdict["fire_max_verdict"]
    assert "fire_max_step_stable" not in verdict
    assert verdict["lambda_max_over_kmax"] == pytest.approx(5.0)
    assert verdict["stable_explicit_mobility_step_um_per_pN"] == pytest.approx(2.0e-7)


def test_a_soft_operator_leaves_headroom_and_says_how_much() -> None:
    cell = _fake_cell(dt_mu=1.0e-8, kmax=1.0e7)
    verdict = driver._explicit_step_verdict(cell, 1.0e7)
    assert verdict["stability_product_base"] == pytest.approx(0.1)
    assert verdict["headroom_factor_at_base"] == pytest.approx(20.0)
    assert verdict["stability_product_at_fire_max"] == pytest.approx(1.0)


def test_the_population_census_reports_dof_not_just_nodes() -> None:
    """The 2026-07-28 retraction happened because no record said how big the system was."""
    census = driver._population_census(_fake_cell(), _fake_solver())
    assert census["n_nodes_total"] == 100
    assert census["n_dof"] == 300
    assert census["n_erm_tethers"] == 9
    assert census["n_membrane_faces"] == 60
    # The per-compartment SIZES, without which a subdivision knob is invisible in the record: on
    # 2026-07-29 two runs 1.73x apart in node count both stamped the same "fraction of native".
    assert census["membrane_vertices"] == 32
    assert census["nucleus_vertices"] == 42
    assert set(census["stiffness_launch_census"]) >= {
        "bending", "crosslink", "steric_wca", "membrane_erm", "nucleus_linc"}


# ── The verdict: what it takes to be allowed to believe a matrix-free spectrum ────────────────────────
def _certification(
    *, conservative: bool = True, family_conservative: bool = True,
    unexercised_family_asymmetric: bool = False,
) -> dict[str, object]:
    """Return a certification record shaped like the driver's, with the failure modes selectable."""
    families = {
        "crosslink": {
            "exercised": True,
            "symmetry": {"conservative": family_conservative},
        },
        "steric_wca": {
            # a family that contributed nothing must not be able to decide the gate either way
            "exercised": False,
            "symmetry": None if not unexercised_family_asymmetric else {"conservative": False},
        },
    }
    return {
        "assembled_symmetry": {
            "conservative": conservative,
            "max_abs_asymmetry_pN_per_um": 1.0e-10 if conservative else 5.0,
            "max_abs_entry_pN_per_um": 1.0e3,
        },
        "families": families,
    }


def _native(*, lambda_max_resolved: bool = True) -> dict[str, object]:
    return {"projected_PKP": {"at_sigma_2theta": {"lambda_max_resolved": lambda_max_resolved}}}


def test_a_run_with_no_certification_stage_refuses_instead_of_scoring_itself() -> None:
    """Lanczos cannot check symmetry; run against an asymmetric operator it reports no warning at all."""
    with pytest.raises(ValueError, match="unchecked premise"):
        driver.form_verdict(None, _native())


def test_the_verdict_passes_only_when_all_three_conjuncts_hold() -> None:
    passed, residual, signal = driver.form_verdict(_certification(), _native())
    assert passed is True
    assert residual == pytest.approx(1.0e-10)
    assert signal == pytest.approx(1.0e3)
    assert residual / signal < 0.01          # and so the void ceiling admits the verdict


def test_a_single_asymmetric_family_fails_the_gate_even_when_the_assembly_looks_clean() -> None:
    """One family's asymmetry can cancel another's in the sum; that cancellation is configuration-specific.

    Scoring only the assembled operator would let a defect that happens to cancel at THIS build be
    carried to a population where it does not.
    """
    passed, _, _ = driver.form_verdict(
        _certification(conservative=True, family_conservative=False), _native())
    assert passed is False


def test_an_unresolved_large_end_fails_the_gate() -> None:
    """An unconverged Ritz value is a property of the iteration count, not of the operator."""
    passed, _, _ = driver.form_verdict(_certification(), _native(lambda_max_resolved=False))
    assert passed is False


def test_an_unexercised_family_cannot_decide_the_gate_in_either_direction() -> None:
    """Its tangent is exactly zero, so it carries no evidence — neither a pass nor a fail."""
    passed, _, _ = driver.form_verdict(
        _certification(unexercised_family_asymmetric=True), _native())
    assert passed is True


def test_the_certification_alone_is_scorable_without_a_native_stage() -> None:
    """The reverse of the refusal above: symmetry stands on its own; a spectrum does not."""
    passed, _, _ = driver.form_verdict(_certification(), None)
    assert passed is True
    assert driver.form_verdict(_certification(conservative=False), None)[0] is False


# ── The ablation itself, against a synthetic operator: the one piece of novel logic here ──────────────
_N_NODES = 4
_N_DOF = 3 * _N_NODES


def _symmetric(seed: int, scale: float) -> np.ndarray:
    root = np.random.default_rng(seed).standard_normal((_N_DOF, _N_DOF))
    return scale * (root + root.T) / 2.0


def _synthetic_operator(cell: SimpleNamespace, solver: SimpleNamespace, *, round_off: bool = False):
    """Return ``(matvec_factory, blocks)`` for a cell whose families are the toggled attributes.

    Mirrors the real operator's structure exactly where it matters: each family contributes a matrix
    only when its gating attribute is live, the contributions are ADDED, and one family contributes a
    zero matrix while its gate is live — which is what the WCA steric kernel does on an overlap-free
    build and is the case a launch-dimension census gets wrong.

    Args:
        cell: The cell whose attributes gate the families.
        solver: The solver carrying the membrane edge topology.
        round_off: Quantize each response to one ulp of the CURRENTLY assembled operator's largest
            entry.  That is what a float64 accumulation actually does, and it is the mechanism that
            makes a leave-one-OUT difference carry the whole assembly's noise instead of the family's.
    """
    # magnitudes are the MEASURED ones from the real cortex build (certify_smoke, 2026-07-28): the
    # crosslink family is ~5 decades above the membrane edge springs, and that spread is the whole
    # reason the differencing construction failed.
    blocks = {
        "bending": _symmetric(1, 4.2),
        "crosslink": _symmetric(2, 7.3e6),
        "membrane_area_edges": _symmetric(3, 3.1e1),
        "steric_wca": np.zeros((_N_DOF, _N_DOF)),          # launches, contributes nothing
        "myosin": _symmetric(5, 2.0),
        "nucleus_linc": _symmetric(6, 9.9e1),
        "membrane_erm": _symmetric(7, 4.5e3),
    }

    def assembled() -> np.ndarray:
        total = np.zeros((_N_DOF, _N_DOF))
        if cell.n_tri:
            total += blocks["bending"]
        if cell.n_xl:
            total += blocks["crosslink"]
        if solver.mem_edges_d is not None:
            total += blocks["membrane_area_edges"]
        if cell.steric is not None:
            total += blocks["steric_wca"]
        if cell.myosin is not None:
            total += blocks["myosin"]
        if cell.nucleus is not None:
            total += blocks["nucleus_linc"]
        if cell.membrane is not None:
            total += blocks["membrane_erm"]
        return total

    def factory(projected: bool):
        # read the attributes NOW, i.e. after whatever toggle is in force
        matrix = assembled()
        if projected:                                     # a crude stand-in projector: drop node 0
            keep = np.ones(_N_DOF)
            keep[:3] = 0.0
            matrix = (keep[:, None] * matrix) * keep[None, :]

        largest = float(np.max(np.abs(matrix))) if matrix.size else 0.0
        ulp = np.spacing(largest) if (round_off and largest > 0.0) else 0.0

        def matvec(vector: np.ndarray) -> np.ndarray:
            response = (matrix @ np.asarray(vector, float).reshape(-1))
            if ulp > 0.0:
                response = np.round(response / ulp) * ulp
            return response.reshape(_N_NODES, 3)

        return matvec

    return factory, blocks


def test_the_ablation_isolates_each_family_and_reproduces_the_assembly() -> None:
    """``K_on − K_off`` per family, summing back to ``K``: this is what licenses a per-family gradient."""
    cell, solver = _fake_cell(), _fake_solver(n=_N_NODES)
    factory, blocks = _synthetic_operator(cell, solver)
    census = driver._stiffness_launch_census(cell, solver)

    record = driver.certify_families(
        cell, solver, launch_census=census, matvec_factory=factory)

    assert record["multilinearity"]["relative_residual"] == pytest.approx(0.0, abs=1e-12)
    for name in ("bending", "crosslink", "myosin", "membrane_erm"):
        measured = record["families"][name]["frobenius_norm_pN_per_um"]
        assert measured == pytest.approx(float(np.linalg.norm(blocks[name], "fro")), rel=1e-10)
    # and every attribute is back where it started, or the LAST family's difference would be wrong
    assert cell.steric is not None and cell.myosin is not None
    assert cell.n_tri == 30 and cell.n_xl == 40 and solver.mem_edges_d is not None


def test_a_family_that_launches_but_contributes_nothing_is_reported_uncertified() -> None:
    """The exact failure a launch-dimension census hides, and the reason the ablation exists."""
    cell, solver = _fake_cell(), _fake_solver(n=_N_NODES)
    factory, _ = _synthetic_operator(cell, solver)
    record = driver.certify_families(
        cell, solver, launch_census=driver._stiffness_launch_census(cell, solver),
        matvec_factory=factory)

    steric = record["families"]["steric_wca"]
    assert steric["exercised"] is False
    assert steric["symmetry"] is None
    assert steric["launch_dimension"] > 0            # it DID launch — that is the point
    assert "steric_wca" in record["families_unexercised"]
    assert "NOT certified" in steric["why_not_exercised"]


def test_an_asymmetric_family_is_caught_even_when_it_is_a_small_part_of_the_assembly() -> None:
    """A family's asymmetry has to be scored against ITS OWN floor, not diluted by the whole operator."""
    cell, solver = _fake_cell(), _fake_solver(n=_N_NODES)
    factory, blocks = _synthetic_operator(cell, solver)
    # make the smallest family non-conservative: an assembly-level check would barely register it
    skew = np.zeros((_N_DOF, _N_DOF))
    skew[0, 1], skew[1, 0] = 0.4, -0.4
    blocks["membrane_area_edges"] = blocks["membrane_area_edges"] + skew

    record = driver.certify_families(
        cell, solver, launch_census=driver._stiffness_launch_census(cell, solver),
        matvec_factory=factory)

    assert record["families"]["membrane_area_edges"]["symmetry"]["conservative"] is False
    assert record["families"]["crosslink"]["symmetry"]["conservative"] is True
    assert driver.form_verdict(record, None)[0] is False


def test_a_small_family_is_not_buried_under_a_large_one_s_round_off() -> None:
    """REGRESSION, measured 2026-07-28: isolating by ``K_on − K_off`` made every family fail.

    The difference of two matrices each accumulated at the scale of the WHOLE assembly carries the
    whole assembly's round-off, so the smallest family's own tangent sat under the largest family's
    noise while being scored against a floor derived from its own much smaller magnitude.  On the real
    cortex build that made all seven families report the identical asymmetry 9.3132e-10 — one half-ulp
    of ``max|K| = 7.30e6`` — and four of them "fail".  Assembling each family ALONE removes the
    subtraction entirely, which is what this asserts.
    """
    cell, solver = _fake_cell(), _fake_solver(n=_N_NODES)
    factory, blocks = _synthetic_operator(cell, solver, round_off=True)
    record = driver.certify_families(
        cell, solver, launch_census=driver._stiffness_launch_census(cell, solver),
        matvec_factory=factory)

    # the smallest exercised family is ~5 decades below the largest; under differencing its asymmetry
    # would be pinned at the assembly's ulp, which sits far ABOVE this family's own round-off floor
    smallest = record["families"]["membrane_area_edges"]["symmetry"]
    assembly_ulp = float(np.spacing(record["assembled_symmetry"]["max_abs_entry_pN_per_um"]))
    assert smallest["round_off_floor_pN_per_um"] < assembly_ulp
    assert smallest["max_abs_asymmetry_pN_per_um"] <= smallest["round_off_floor_pN_per_um"]
    assert smallest["conservative"] is True

    for name in ("bending", "crosslink", "nucleus_linc", "membrane_erm"):
        assert record["families"][name]["symmetry"]["conservative"] is True, name
    assert driver.form_verdict(record, None)[0] is True


# ── The console summary: a pure function BECAUSE its absence cost a native run ────────────────────────
def _stability() -> dict[str, object]:
    return driver._explicit_step_verdict(_fake_cell(dt_mu=5.0e-8, kmax=2.0e6), 1.0e7)


def _projected(**extra: object) -> dict[str, object]:
    block = {
        "at_sigma_2theta": {
            "lambda_max_pN_per_um": 8.33e6, "lambda_max_resolved": True,
            "lambda_min_pN_per_um": 52.1, "lambda_min_resolved": False,
        },
        "shift_invariance": {"both_resolved": False},
    }
    block.update(extra)
    return block


def test_the_native_summary_asks_only_for_keys_the_verdict_actually_produces() -> None:
    """REGRESSION, and it cost a full-native run on 2026-07-28.

    The summary was three inline prints; one asked ``stability`` for ``fire_max_step_stable``, a key
    ``_explicit_step_verdict`` deliberately does not produce (FIRE is damped MD; the forward-Euler
    limit does not adjudicate it).  Renaming it to ``fire_max_verdict`` updated the tests and the
    figure and missed the print, so a ``KeyError`` fired AFTER the expensive Lanczos pass had already
    produced its answer — the physics finished and the run was lost at a format string, with no
    ``record.json`` written.  ``main()`` is not under test; this function is.
    """
    lines = driver.native_summary_lines(_stability(), _projected())
    assert any("λ_max(PKP)" in line for line in lines)
    assert any("dt_mu·λ_max" in line for line in lines)
    assert any("not adjudicated here" in line for line in lines)


def test_the_summary_raises_in_a_unit_test_rather_than_on_the_device() -> None:
    """The failure mode has to stay LOUD and CHEAP; silently skipping a line would hide a rename."""
    broken = dict(_stability())
    del broken["stability_product_at_fire_max"]
    with pytest.raises(KeyError):
        driver.native_summary_lines(broken, _projected())


def test_the_summary_reports_the_sweep_only_when_one_was_run() -> None:
    assert not any("sweep" in line for line in
                   driver.native_summary_lines(_stability(), _projected()))
    with_sweep = _projected(small_end_sweep={
        "by_dimension": [{"n_iterations": 120, "lambda_min_pN_per_um": 52.1,
                          "bound_over_value": 25.2, "relative_disagreement": 0.048}],
        "verdict": "OUT_OF_REACH"})
    lines = driver.native_summary_lines(_stability(), with_sweep)
    assert any("m= 120" in line for line in lines)
    assert any("OUT_OF_REACH" in line for line in lines)


# ── The small-end sweep: turning one ambiguous point into a verdict ───────────────────────────────────
def _entry(m: int, bound: float, disagreement: float, *, resolved: bool = False) -> dict[str, object]:
    return {"n_iterations": m, "bound_over_value": bound,
            "relative_disagreement": disagreement, "resolved": resolved}


def test_a_sweep_of_one_dimension_is_refused_because_it_measures_nothing() -> None:
    with pytest.raises(ValueError, match="at least two"):
        driver.classify_small_end_sweep([_entry(120, 25.0, 0.05)])


def test_both_witnesses_falling_is_CONVERGING() -> None:
    """Then the remedy is a larger Krylov dimension, and the sweep says so."""
    result = driver.classify_small_end_sweep(
        [_entry(120, 25.0, 0.05), _entry(240, 8.0, 0.02), _entry(480, 2.0, 0.004)])
    assert result["verdict"] == "CONVERGING"
    assert result["krylov_span"] == pytest.approx(4.0)


def test_a_static_bound_is_OUT_OF_REACH_and_says_do_not_re_queue() -> None:
    """The permanent negative result: the point of running the sweep at all."""
    result = driver.classify_small_end_sweep(
        [_entry(120, 25.0, 0.048), _entry(240, 24.6, 0.047), _entry(480, 24.9, 0.049)])
    assert result["verdict"] == "OUT_OF_REACH"
    assert "DO NOT RE-QUEUE" in result["reason"]


def test_one_witness_improving_is_not_enough() -> None:
    """A falling bound with a static shift disagreement is the start vector settling, not the operator."""
    result = driver.classify_small_end_sweep(
        [_entry(120, 25.0, 0.048), _entry(240, 12.0, 0.048), _entry(480, 6.0, 0.049)])
    assert result["verdict"] == "OUT_OF_REACH"


def test_an_actually_resolved_end_outranks_the_trajectory_test() -> None:
    """If Parlett's bound got under 1% of the value, the answer is the operator's, however it got there."""
    result = driver.classify_small_end_sweep(
        [_entry(120, 25.0, 0.05), _entry(240, 0.004, 1e-5, resolved=True)])
    assert result["verdict"] == "RESOLVED"


def test_the_iterative_path_is_scored_against_the_dense_eigensolve_of_the_same_operator() -> None:
    """Lanczos is used where nothing can check it, so it is checked where something can."""
    cell, solver = _fake_cell(), _fake_solver(n=_N_NODES)
    factory, _ = _synthetic_operator(cell, solver)
    record = driver.certify_families(
        cell, solver, launch_census=driver._stiffness_launch_census(cell, solver),
        matvec_factory=factory)

    check = record["instrument_check_against_ground_truth"]
    assert check["relative_difference"] == pytest.approx(0.0, abs=1e-10)
    assert check["gershgorin_over_lambda_max"] >= 1.0
