"""AST-only guards for the Codex NG-2/3/6 foundation hardening (no Warp launches)."""

from __future__ import annotations

import ast
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[4]
CELL = ROOT / "aleph" / "components" / "incumbent"
FLUID = ROOT / "aleph" / "components" / "fluid"


def _function(path: Path, name: str) -> ast.FunctionDef:
    tree = ast.parse(path.read_text(encoding="utf-8"))
    for node in ast.walk(tree):
        if isinstance(node, ast.FunctionDef) and node.name == name:
            return node
    raise AssertionError(f"{name} not found in {path}")


def _source(path: Path, node: ast.AST) -> str:
    return ast.get_source_segment(path.read_text(encoding="utf-8"), node) or ""


def _assert_no_d2h(path: Path, function: str) -> None:
    src = _source(path, _function(path, function))
    assert ".numpy(" not in src, f"authoritative D2H read in hot function {path.name}:{function}"
    assert "pressure_to_host(" not in src
    assert "total_content(" not in src


def test_physical_and_inner_hot_paths_have_no_d2h_readback() -> None:
    _assert_no_d2h(FLUID / "domain.py", "remap")
    _assert_no_d2h(FLUID / "scheduler.py", "outer_step")
    _assert_no_d2h(FLUID / "scheduler.py", "run")
    _assert_no_d2h(CELL / "driver.py", "make_inner_solve")


def test_live_mask_refresh_has_no_host_geometry_download() -> None:
    path = CELL / "live_mesh_domain.py"
    _assert_no_d2h(path, "_refresh")
    src = path.read_text(encoding="utf-8")
    assert "self.mesh.refit()" in src
    assert "support_winding_number=True" in src


def test_pressure_traction_is_in_canonical_force_assembly() -> None:
    path = CELL / "driver.py"
    src = _source(path, _function(path, "_accumulate_all"))
    assert "cell.membrane_pressure.accumulate(pos, f)" in src
    trace_src = (FLUID / "surface_trace.py").read_text(encoding="utf-8")
    assert "masked_affine_pressure_trace" in trace_src
    assert "wp.inverse(normal)" in trace_src


def test_live_membrane_not_static_sphere_when_compartment_exists() -> None:
    path = CELL / "assemble.py"
    src = _source(path, _function(path, "build_cell"))
    assert "LiveMeshMembraneMaskProvider" in src
    assert "fluid_membrane_provider = mem_provider or StaticSphereMembraneProvider" in src


#: The canonical runtime surface. `ac/` split into these three on 2026-08-09; naming the old path
#: would have scanned an empty directory and passed.
CANONICAL_ROOTS = (
    ROOT / "aleph" / "engine",
    ROOT / "aleph" / "components",
    ROOT / "aleph" / "laws",
)


def test_the_canonical_roots_are_not_empty() -> None:
    """A static contract over a directory that does not exist is not a contract."""
    for root in CANONICAL_ROOTS:
        assert root.is_dir() and any(root.rglob("*.py")), f"{root} would be scanned vacuously"


def _without_docstrings(source: str) -> str:
    """Source with docstrings **and `#` comments** blanked, so a rule about code is not tripped by prose.

    Widening this scan to `laws/` on 2026-08-09 produced exactly one hit, and it was
    ``device="cuda:0"`` inside a docstring in `network_warp.py` — an example of the thing the rule
    forbids, written down to explain the rule. A text scan cannot tell an example from an
    instruction; the parser can.
    """
    try:
        tree = ast.parse(source)
    except SyntaxError:
        return source
    spans: list[tuple[int, int]] = []
    for node in ast.walk(tree):
        if not isinstance(node, (ast.Module, ast.ClassDef, ast.FunctionDef, ast.AsyncFunctionDef)):
            continue
        body = getattr(node, "body", None)
        if not body:
            continue
        first = body[0]
        if (
            isinstance(first, ast.Expr)
            and isinstance(first.value, ast.Constant)
            and isinstance(first.value.value, str)
            and first.lineno is not None
            and first.end_lineno is not None
        ):
            spans.append((first.lineno, first.end_lineno))
    lines = source.splitlines(keepends=True)
    for start, end in spans:
        for index in range(start - 1, min(end, len(lines))):
            lines[index] = "\n"

    # `#` comments too. Stripping docstrings alone still left three prose mentions of
    # ``device="cuda:0"`` in `gamma_floor.py`, written to explain the very rule being checked.
    # String *literals* are deliberately kept: `device: str = "cuda:0"` is the defect.
    out: list[str] = []
    for line in lines:
        quote = None
        cut = None
        for index, ch in enumerate(line):
            if quote:
                if ch == quote and (index == 0 or line[index - 1] != "\\"):
                    quote = None
            elif ch in "\"'":
                quote = ch
            elif ch == "#":
                cut = index
                break
        out.append(line if cut is None else line[:cut] + "\n")
    return "".join(out)


#: Real hard-coded ordinals still in the tree, found when this scan widened to `laws/` on
#: 2026-08-09. `gamma_floor_sweep.py` takes `device: str = "cuda:0"` as a function default — the
#: charter forbids exactly this ("No hard-coded device IDs"). **Ratcheted, not fixed**: changing a
#: default changes what the driver runs on, and this module's numbers are already blocked by
#: `STATE.md` (c) 2. It may only shrink.
CUDA_ORDINAL_DEBT = 1


def test_canonical_path_has_no_hard_coded_cuda_ordinal() -> None:
    offenders = [
        f"{path.relative_to(ROOT)}"
        for root in CANONICAL_ROOTS
        for path in root.rglob("*.py")
        if re.search(r"cuda:\d+", _without_docstrings(path.read_text(encoding="utf-8")))
    ]
    assert len(offenders) <= CUDA_ORDINAL_DEBT, (
        f"hard-coded device ordinals rose from {CUDA_ORDINAL_DEBT} to {len(offenders)}: {offenders}"
    )


def test_the_docstring_stripper_only_strips_docstrings() -> None:
    """Vacuity: a stripper that blanked everything would make the rule above unfailable."""
    sample = '''"""a docstring mentioning cuda:0"""\nDEVICE = "cuda:0"\n'''
    stripped = _without_docstrings(sample)
    assert "cuda:0" in stripped, "the real assignment must survive"
    assert stripped.count("cuda:0") == 1, "the docstring mention must not"


def test_ac_tests_do_not_launch_production_kernels_on_cpu() -> None:
    for path in (ROOT / "aleph" / "tests" / "ac").rglob("*.py"):
        src = path.read_text(encoding="utf-8")
        forbidden = (
            r"\bDEV(?:ICE)?\s*=\s*['\"]cpu['\"]",
            r"\bdevice\s*=\s*['\"]cpu['\"]",
            r"wp\.get_device\(\s*['\"]cpu['\"]",
            r"wp\.ScopedDevice\(\s*['\"]cpu['\"]",
        )
        assert not any(re.search(pattern, src) for pattern in forbidden), f"Warp-CPU launch in {path}"


def test_runnable_native_gates_require_and_pass_resolved_cuda() -> None:
    paths = (
        ROOT / "aleph/components/motor/native_gates/ng0_motor_parity.py",
        ROOT / "aleph/components/nucleus/native_gates/ng0_nucleus_parity.py",
        ROOT / "aleph/components/solid/native_gates/ng0_steric_parity.py",
    )
    for path in paths:
        src = path.read_text(encoding="utf-8")
        assert "if not dev.is_cuda:" in src, f"native gate can fall back to Warp CPU: {path}"
        assert "device = str(dev)" in src
        assert "device=DEVICE" not in src


def test_nucleus_live_mask_has_no_per_remap_host_download() -> None:
    path = CELL / "compartments.py"
    src = _source(path, _function(path, "classify_nucleus"))
    assert ".numpy(" not in src
    assert "self._live.classify_nucleus(grid, mask)" in src


def test_inner_convergence_checks_constraints_and_finite_state() -> None:
    src = (CELL / "driver.py").read_text(encoding="utf-8")
    assert "max_constraint_error_kernel" in src
    assert "finite_vec3_kernel" in src
    assert "project_constraint_forces_kernel" in src
    assert "inputs=[pos, dt_mu_d, projected_f_d, active_d]" in src
    convergence = (CELL / "inner_mechanics.py").read_text(encoding="utf-8")
    assert "P = I - J^T (J J^T)^-1 J" in convergence
    assert "max_constraint_error[0] <= tolerance_um" in convergence
    assert "max_projected_force[0] * dt_mu[0] <= tolerance_um" in convergence
    assert "if finite[0] == 0" in convergence
    assert "active[0] = 0" in convergence


def test_membrane_source_does_not_delete_interior_source() -> None:
    src = _source(FLUID / "scheduler.py", _function(FLUID / "scheduler.py", "outer_step"))
    assert "g.s_water.zero_()" not in src
    assert "self.membrane_bc.apply" in src
    field_src = (FLUID / "field_grid.py").read_text(encoding="utf-8")
    assert "self.s_membrane" in field_src and "self.s_total" in field_src


def test_production_membrane_flux_uses_live_triangle_area_not_voxel_faces() -> None:
    boundary_src = (FLUID / "boundary.py").read_text(encoding="utf-8")
    assert "live_triangle_membrane_flux_kernel" in boundary_src
    assert "q_face = influx * wp.float64(0.5) * area2" in boundary_src
    assert "q_face * w / w_sum * inv_cell_volume" in boundary_src
    assemble_src = _source(CELL / "assemble.py", _function(CELL / "assemble.py", "build_cell"))
    assert "faces=membrane.faces_d" in assemble_src
    assert "pos=pos_d" in assemble_src


def test_rejected_outer_step_restores_authoritative_device_fields() -> None:
    src = _source(FLUID / "scheduler.py", _function(FLUID / "scheduler.py", "outer_step"))
    for field in ("g.p", "g.p_new", "g.mask", "g.div_vs", "g.s_membrane", "g.s_total"):
        assert field in src
    assert "_conditional_restore_f64_3d_kernel" in src
    assert "_conditional_restore_i32_3d_kernel" in src
    assert "_advance_time_if_accepted_kernel" in src
    assert "rollback(self._outer_accepted_d)" in src
    assert "_finite_f64_3d_kernel" in src
    assert "dt_phys must be finite and positive" in src
    assert src.index("_finite_f64_3d_kernel") < src.index("self.inner_solve(dt_phys, self._fluid_finite_d)")


def test_erm_force_is_pure_and_rupture_is_commit_only() -> None:
    erm_path = CELL / "erm_tether.py"
    force_src = _source(erm_path, _function(erm_path, "erm_tether_force_kernel"))
    bell_force_src = _source(erm_path, _function(erm_path, "erm_bell_force_kernel"))
    assert "bound[k] = 0" not in force_src
    assert "bound[k] = 0" not in bell_force_src
    assert "rest[k] =" not in bell_force_src
    assert "tension <= wp.float64(0.0)" in force_src
    assert "tension > f_rupt" in force_src
    assert "f_rupt" not in bell_force_src
    compartments = (CELL / "compartments.py").read_text(encoding="utf-8")
    assert "erm_tether_force_kernel" in compartments
    assert "erm_bell_force_kernel" in compartments
    assert "erm_bell_kmc_kernel" in compartments
    assert "wp.launch(erm_tether_kernel" not in compartments
    driver = _source(CELL / "driver.py", _function(CELL / "driver.py", "make_inner_solve"))
    assert "cell.membrane.commit_kinetics(cell.pos_d, accepted_d, current_dt_phys, seed)" in driver
    scheduler = _source(FLUID / "scheduler.py", _function(FLUID / "scheduler.py", "outer_step"))
    assert "commit_irreversible(self._outer_accepted_d)" in scheduler
    gates = (CELL / "native_gates" / "ng2_ng3_ng6.py").read_text(encoding="utf-8")
    assert '"erm_commit_only_rupture"' in gates
    assert 'native_ledger.get("erm_density_mcf7_production") is True' in gates
    assert 'native_ledger.get("erm_kinetics_mode") == "BELL_SLIP_ON_OFF"' in gates
    assert 'native_ledger.get("erm_rebind_rest_policy") == "formation_length"' in gates
    assert 'native_ledger.get("preload_capacity_force_basis") == "SOURCE_GROUNDED_SINGLE_ERM_FORCE"' in gates


def test_state_dump_cannot_label_rejected_candidate_stable() -> None:
    dump = _source(CELL / "dump_state.py", _function(CELL / "dump_state.py", "dump_state"))
    for required in ("inner_converged", "outer_accepted", "outer_rolled_back", "residual_committed"):
        assert required in dump
    assert "np.isfinite(f_total_post_all)" in dump


def test_ng2_coupled_ledger_is_executable_while_ng3_still_requires_native_evidence() -> None:
    gates = (CELL / "native_gates" / "ng2_ng3_ng6.py").read_text(encoding="utf-8")
    assert '"device_terzaghi_dirichlet"' in gates
    assert '"device_green_impulse"' in gates
    assert '"native_config_pass"' in gates
    assert '"native_population_pass"' in gates
    assert '"native_diagnostics_pass"' in gates
    assert '"canonical_coupled_outer_step"' in gates
    assert '"chunk_partition_trajectory"' in gates
    assert '"membrane_surface_flux_integral"' in gates
    assert '"solid_dilatation_integral"' in gates
    assert '"moving_face_delta"' in gates
    assert '"solver_conservation_relative_error"' in gates
    assert '"physical_conservation_relative_error"' in gates
    assert '"time_index_contract"' in gates
    assert "fsi_coupling.update" in gates
    assert '"rigid_translation_zero_divergence"' in gates
    assert '"affine_dilation_integral_identity"' in gates
    assert '"spread_partition_conservation"' in gates
    assert '"spatial_affine_pressure_bulk_plus_surface_net_force"' in gates
    assert '"fsi_transfer_adjoint_work_identity"' in gates
    assert "_node_work_kernel" in gates
    assert "_grid_pressure_work_kernel" in gates
    assert "n_filaments=70686" in gates
    assert '"FAIL_COUPLED"' in gates


def test_ng5_uses_one_masked_adjoint_and_the_live_discrete_volume() -> None:
    """The native NG-5 pass must come from matched operators/geometry, not a post-hoc force correction."""
    gates = (CELL / "native_gates" / "ng2_ng3_ng6.py").read_text(encoding="utf-8")
    biot = (FLUID / "biot_substrate.py").read_text(encoding="utf-8")
    fsi = (CELL / "fsi_coupling.py").read_text(encoding="utf-8")
    assemble = (CELL / "assemble.py").read_text(encoding="utf-8")
    assert "normal = wp.mat44d()" in biot and "coefficients = wp.mul(wp.inverse(normal), rhs)" in biot
    assert "weight_sum" in biot and "mask[ii, jj, kk] == _FLUID" in biot
    assert "weight_sum" in fsi and "mask[ii, jj, kk] != _FLUID" in fsi
    assert "g.mask" in biot and "g.mask" in fsi
    assert "signed_volume(mem_verts, membrane.mesh.faces)" in assemble
    assert '"bulk_oracle_relative_error"' in gates
    assert '"net_force_relative_error"' in gates
    assert '"nonnegative_dissipation_pass"' in gates


def test_scheduler_ledger_integrates_and_predicates_every_coupled_source() -> None:
    src = (FLUID / "scheduler.py").read_text(encoding="utf-8")
    for required in (
        "_accumulate_grid_conservation_kernel",
        "_accumulate_scalar_rate_kernel",
        "_finalize_conservation_kernel",
        "membrane_surface_flux_integral_d",
        "membrane_grid_source_integral_d",
        "interior_source_integral_d",
        "solid_dilatation_integral_d",
        "membrane_deposition_error_d",
        "solver_conservation_error_d",
        "physical_conservation_error_d",
    ):
        assert required in src
    outer = _source(FLUID / "scheduler.py", _function(FLUID / "scheduler.py", "outer_step"))
    assert outer.index("_accumulate_grid_conservation_kernel") < outer.index("sub.step(dt_sub)")
    assert outer.index("self.domain.remap()") < outer.index("_finalize_conservation_kernel")
    assert "_zero_if_rejected_kernel" in outer


def test_moving_remap_reads_immutable_old_pressure() -> None:
    src = _source(FLUID / "domain.py", _function(FLUID / "domain.py", "_remap_transport_kernel"))
    assert "p_old" in src
    assert "s += p_old" in src


def test_scheduler_wires_solid_dilatation_through_inner_solve() -> None:
    src = _source(CELL / "driver.py", _function(CELL / "driver.py", "make_inner_solve"))
    assert "SolidDilatationCoupling" in src
    assert "update_from_displacement" in src


def test_retry_chunks_are_device_predicated_and_do_not_change_physical_time() -> None:
    driver = _source(CELL / "driver.py", _function(CELL / "driver.py", "make_inner_solve"))
    mechanics = (CELL / "inner_mechanics.py").read_text(encoding="utf-8")
    assert "max_inner_retries" in driver
    assert "begin_inner_attempt_kernel" in driver
    assert "attempt * n_inner + local_it + 1" in driver
    assert "dt_phys" not in _source(
        CELL / "inner_mechanics.py", _function(CELL / "inner_mechanics.py", "conditional_axpy_kernel"))
    assert "if active[0] != 0" in mechanics
    assert "attempts[0] += 1" in mechanics


def test_multistep_scheduler_uses_device_fail_stop_without_acceptance_readback() -> None:
    scheduler_path = FLUID / "scheduler.py"
    outer = _source(scheduler_path, _function(scheduler_path, "outer_step"))
    run = _source(scheduler_path, _function(scheduler_path, "run"))
    assert "_attempted" not in outer
    assert "_deactivate_run_if_rejected_kernel" in outer
    assert "wp.copy(self._attempt_enabled_d, self._run_active_d)" in outer
    assert "for _ in range(n_phys)" in run
    assert ".numpy(" not in outer and ".numpy(" not in run
