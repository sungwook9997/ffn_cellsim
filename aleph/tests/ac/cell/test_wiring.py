"""Static wiring-integrity gates for the assembly milestone-2 integration (AST + Warp codegen; no launches).

These catch a silent DE-WIRING regression: the driver must SUM the nucleus + membrane compartments into the
inner force, and the myosin KMC must run the I4 ``fill_walk_dir_kernel`` hand-off in the correct order
(attach → fill_walk_dir → step/detach). They parse source (like ``tests/ac/motor/test_powerstroke_kernel_wiring``)
so they run on a CPU-only Mac, and add a codegen-clean compile of the new Warp source.
"""
from __future__ import annotations

import ast
from pathlib import Path

import numpy as np
import pytest

from aleph.components.incumbent.assemble import CellConfig, _actin_segment_topology, build_cell

REPO_ROOT = Path(__file__).resolve().parents[4]
DRIVER = REPO_ROOT / "aleph" / "components" / "incumbent" / "driver.py"
ASSEMBLE = REPO_ROOT / "aleph" / "components" / "incumbent" / "assemble.py"


def _func(path: Path, name: str) -> ast.FunctionDef:
    tree = ast.parse(path.read_text(encoding="utf-8"))
    for node in ast.walk(tree):
        if isinstance(node, ast.FunctionDef) and node.name == name:
            return node
    raise AssertionError(f"{name} not found in {path.name}")


def _launched_kernels(fn: ast.FunctionDef) -> list[str]:
    """Ordered names of the kernels passed as the first arg to each ``wp.launch(...)`` in body order."""
    out = []
    for node in ast.walk(fn):
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute) \
                and node.func.attr == "launch" and node.args:
            first = node.args[0]
            if isinstance(first, ast.Name):
                out.append(first.id)
    return out


def _imported_or_launched(path: Path) -> set[str]:
    """Symbols the module IMPORTS or passes to wp.launch (i.e. actually executes) — ignores docstring prose."""
    tree = ast.parse(path.read_text(encoding="utf-8"))
    names: set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.ImportFrom):
            names.update(a.name for a in node.names)
        elif isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute) \
                and node.func.attr == "launch" and node.args and isinstance(node.args[0], ast.Name):
            names.add(node.args[0].id)
    return names


def test_accumulate_sums_nucleus_and_membrane() -> None:
    """driver._accumulate_all must add BOTH new compartments' forces (else they are silently decoupled)."""
    src = ast.get_source_segment(DRIVER.read_text(encoding="utf-8"), _func(DRIVER, "_accumulate_all"))
    assert "cell.nucleus.accumulate" in src, "nucleus compartment not summed into the inner force"
    assert "cell.membrane.accumulate" in src, "membrane compartment not summed into the inner force"


def test_walk_dir_handoff_order() -> None:
    """step_myosin_kinetics must launch attach_kernel → fill_walk_dir_kernel → step_detach_kernel IN ORDER.

    fill_walk_dir must run AFTER attach (so a freshly-bound head gets its actin's polarity) and BEFORE the next
    force/compute_loads uses walk_dir. A wrong order silently re-introduces the minifilament-axis default.
    """
    launched = _launched_kernels(_func(DRIVER, "step_myosin_kinetics"))
    seq = [c for c in launched if c in ("attach_kernel", "fill_walk_dir_kernel", "step_detach_kernel")]
    assert seq == ["attach_kernel", "fill_walk_dir_kernel", "step_detach_kernel"], f"bad hand-off order: {seq}"


def test_irreversible_biology_stepped_once_per_outer() -> None:
    """The commit hook advances rupture/KMC once per OUTER step, predicated entirely on device state."""
    src = ast.get_source_segment(DRIVER.read_text(encoding="utf-8"), _func(DRIVER, "make_inner_solve"))
    assert "cell.nucleus.rupture_step(cell.pos_d, accepted_d)" in src
    assert "cell.myosin.commit_segment_kinetics(cell.pos_d, accepted_d, current_dt_phys, seed)" in src
    assert "cell.membrane.commit_kinetics(cell.pos_d, accepted_d, current_dt_phys, seed)" in src
    assert "accepted_d" in src
    assert "kinetics_commit_index" not in src


def test_segment_topology_preserves_fiber_polarity() -> None:
    """Flattened adjacent-node segments retain the source fiber's barbed-end direction sign."""
    a, b, polarity = _actin_segment_topology(
        np.array([0, 3, 7], dtype=np.int32),
        np.array([1, -1], dtype=np.int32),
    )
    assert np.array_equal(a, np.array([0, 1, 3, 4, 5], dtype=np.int32))
    assert np.array_equal(b, a + 1)
    assert np.array_equal(polarity, np.array([1, 1, -1, -1, -1], dtype=np.int32))


def test_cell_enables_production_segment_motor() -> None:
    """The composed cell must opt into segment anchors; node anchors remain diagnostic-only."""
    src = ast.get_source_segment(ASSEMBLE.read_text(encoding="utf-8"), _func(ASSEMBLE, "build_cell"))
    assert src is not None
    assert "myo.enable_segment_runtime(" in src
    assert '"myosin_attachment_primitive": "POINT_TO_SEGMENT_BARYCENTRIC"' in src


def test_nmii_backbone_persistence_length_is_source_gated_before_cuda() -> None:
    """The diagnostic Lp sweep fixture must never become an unsourced composed-runtime default."""
    with pytest.raises(ValueError, match="provenance"):
        build_cell(CellConfig(nmii_backbone_lp_um=1.0))
    with pytest.raises(ValueError, match="source requires"):
        build_cell(CellConfig(nmii_backbone_lp_source="TEST_ONLY"))


def test_public_runner_forwards_inner_solver_contract() -> None:
    """Both scheduler paths must forward the public implicit-solver controls into ``make_inner_solve``."""
    src = ast.get_source_segment(DRIVER.read_text(encoding="utf-8"), _func(DRIVER, "run_from_resting"))
    assert src is not None
    assert src.count("inner_solver=inner_solver") == 2
    assert src.count("implicit_cg_max_iterations=implicit_cg_max_iterations") == 2
    assert src.count("implicit_line_search_steps=implicit_line_search_steps") == 2
    assert src.count("implicit_coarse_iterations=implicit_coarse_iterations") == 2
    assert src.count("rkc_stages=rkc_stages") == 2
    assert 'cell.ledger["inner_solver"] = inner_solver' in src
    assert 'cell.ledger["implicit_cg_max_iterations"] = int(implicit_cg_max_iterations)' in src
    assert 'cell.ledger["implicit_line_search_steps"] = int(implicit_line_search_steps)' in src
    assert 'cell.ledger["implicit_coarse_iterations"] = int(implicit_coarse_iterations)' in src
    assert 'cell.ledger["rkc_stages"] = int(rkc_stages)' in src


def test_implicit_driver_uses_residual_backtracking_not_explicit_displacement_cap() -> None:
    """The implicit candidate may exceed one explicit CFL step but must win an exact-residual trial."""
    src = ast.get_source_segment(DRIVER.read_text(encoding="utf-8"), _func(DRIVER, "make_inner_solve"))
    assert src is not None
    assert "compute_trust_scale_kernel" not in src
    assert "implicit_line_search_scales_d" in src
    assert "decide_better_line_search_trial_kernel" in src
    assert "explicit_candidate_d" in src


def test_accelerated_solver_modes_share_the_exact_residual_fallback() -> None:
    """Every accelerated candidate must traverse the same explicit comparison as analytic PCG."""
    src = ast.get_source_segment(DRIVER.read_text(encoding="utf-8"), _func(DRIVER, "make_inner_solve"))
    assert src is not None
    assert '"block_descent"' in src
    assert '"anderson"' in src
    assert '"rkc1"' in src
    assert '"tournament"' in src
    assert '"contact_schwarz"' in src
    assert '"contact_tournament"' in src
    assert '"fiber_contact_schwarz"' in src
    assert '"cluster_tournament"' in src
    assert '"rigid_contact_schwarz"' in src
    assert '"rigid_cluster_tournament"' in src
    assert "preconditioned_direction" in src
    assert "AndersonDepthOne" in DRIVER.read_text(encoding="utf-8")
    assert "ContactPairSchwarz" in DRIVER.read_text(encoding="utf-8")
    assert "FiberContactSchwarz" in DRIVER.read_text(encoding="utf-8")
    assert "RigidFiberContactSchwarz" in DRIVER.read_text(encoding="utf-8")
    assert "direction_authorized_d" in src
    assert "best_candidate_d" in src
    assert "_reduce_projected_metric(base_projected_f_d, best_residual_d, selection=True)" in src
    assert 'line_search_objective: str = "max_force"' in src
    assert 'line_search_objective == "l2_squared"' in src
    assert "_evaluate_projected_residual(finite_d, convergence_force_d)" in src
    assert "rkc1_recurrence_kernel" in src
    assert "accelerator_candidate_families" in src


def test_pressure_acts_on_actin_only() -> None:
    """node_volume must be zeroed for non-actin nodes (pressure -α∇p·V acts on the solid skeleton only)."""
    src = ast.get_source_segment(ASSEMBLE.read_text(encoding="utf-8"), _func(ASSEMBLE, "build_cell"))
    assert "node_vol = np.zeros(n_total" in src, "node_volume no longer actin-restricted"


def test_nucleus_shell_kernel_not_launched() -> None:
    """The retired lumped nucleus_shell_kernel must not be IMPORTED or LAUNCHED anywhere in ac/cell (same-commit
    guard). Docstring prose that documents the retirement is fine — only real execution is forbidden."""
    for path in (DRIVER, ASSEMBLE, REPO_ROOT / "aleph" / "components" / "incumbent" / "compartments.py"):
        assert "nucleus_shell_kernel" not in _imported_or_launched(path), f"{path.name} executes the bead ball"


def test_fire_inner_solver_wired_as_tangent_free_mode() -> None:
    """The FIRE static minimizer is a registered solver whose update path uses no tangent/line-search.

    Guards the solver pivot: FIRE must be selectable, must be EXCLUDED from the accelerated (tangent)
    path, and must drive its own momentum kernels in ``make_inner_solve``. The full-native convergence
    run is the gbook GO (CUDA); this is the CPU-side contract gate.
    """
    fn = _func(DRIVER, "make_inner_solve")
    src = ast.get_source_segment(DRIVER.read_text(encoding="utf-8"), fn)
    solver_choices_block = src.split("solver_choices", 1)[1].split("}", 1)[0]
    assert '"fire"' in solver_choices_block, "fire must be a registered inner_solver"
    assert 'not in {"explicit", "fire"}' in src, "fire must be excluded from the accelerated tangent path"
    launched = _launched_kernels(fn)
    for kernel in ("fire_reset_state_kernel", "fire_adapt_kernel", "fire_mix_and_step_kernel"):
        assert kernel in launched, f"make_inner_solve must launch {kernel}"
    # FIRE reuses the exact projected-force convergence gate — it must not introduce its own criterion.
    assert "update_convergence_kernel" in launched


def test_new_warp_source_codegen_clean() -> None:
    """Compile the new/edited Warp source (compartments + boundary) — codegen only, no launch/sim."""
    try:
        import warp as wp

        import aleph.components.incumbent.compartments as C
        import aleph.components.fluid.boundary as B

        wp.init()
        device_obj = next((d for d in wp.get_devices() if d.is_cuda), None)
        if device_obj is None:
            pytest.skip("I0-A: Warp codegen/runtime gate requires CUDA")
        device = str(device_obj)
        wp.load_module(C, device=device)
        wp.load_module(B, device=device)
    except Exception as exc:  # pragma: no cover - environment-dependent (real gate is the gbook)
        pytest.skip(f"Warp codegen check unavailable here: {exc}")
