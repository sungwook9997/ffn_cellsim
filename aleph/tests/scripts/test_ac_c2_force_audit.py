"""Structural gate for rejected-candidate C-2 force-family instrumentation."""

from pathlib import Path


def test_candidate_force_audit_reads_the_device_capture_without_accepting_it() -> None:
    script = Path(__file__).resolve().parents[2] / "scripts" / "ac_c2_force_audit.py"
    source = script.read_text(encoding="utf-8")
    assert 'parser.add_argument("--candidate-solver"' in source
    assert "capture_candidate=True" in source
    assert "implicit_coarse_modes=args.coarse_modes" in source
    assert "implicit_coarse_iterations=args.coarse_iterations" in source
    assert "line_search_objective=args.line_search_objective" in source
    assert "audited_position_d = inner_solve.candidate_pos_d" in source
    assert '"authoritative_geometry_rolled_back": True' in source
    assert "audit_cell = copy.copy(cell)" in source
    assert "node_pos=audited_position_d" in source
    assert "_accumulate_all(audit_cell, audited_position_d" in source
    assert "_project(audit_cell, audited_position_d, work_d)" in source
