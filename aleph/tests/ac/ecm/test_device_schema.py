"""Host-visible configuration and source contracts for ECM device topology."""

from __future__ import annotations

import ast
from pathlib import Path

import pytest

from aleph.components.ecm.mikado_topology import MikadoInitConfig

ECM_ROOT = Path(__file__).resolve().parents[3] / "components" / "ecm"


def _config(**changes: object) -> MikadoInitConfig:
    values: dict[str, object] = {
        "box_lo_um": (-5.0, -4.0, -3.0),
        "box_hi_um": (5.0, 4.0, 3.0),
        "n_fibers": 12,
        "fiber_length_um": 5.0,
        "target_segment_um": 0.6,
        "crosslink_capture_um": 0.2,
        "pin_faces": ("z_lo", "x_hi"),
        "pin_margin_um": 0.5,
        "rng_seed": 17,
        "max_refinement_level": 2,
        "persistent_id_base": 1_000,
    }
    values.update(changes)
    return MikadoInitConfig(**values)  # type: ignore[arg-type]


def test_capacity_derivation_never_coarsens_above_target_segment_length() -> None:
    cfg = _config(fiber_length_um=5.0, target_segment_um=0.6)
    assert cfg.nodes_per_fiber == 10
    assert cfg.max_segment_length_um == pytest.approx(5.0 / 9.0)
    assert cfg.max_segment_length_um <= cfg.target_segment_um
    assert cfg.segment_slots_per_fiber == 36
    assert cfg.node_slots_per_fiber == 37
    assert cfg.bend_slots_per_fiber == 35


def test_hashgrid_shape_is_derived_from_geometry_and_exact_query_bound() -> None:
    cfg = _config()
    radius = cfg.max_segment_length_um + cfg.crosslink_capture_um
    assert cfg.grid_dimensions == tuple(
        max(1, __import__("math").ceil(extent / radius)) for extent in (10.0, 8.0, 6.0)
    )


@pytest.mark.parametrize(
    ("change", "message"),
    [
        ({"pin_faces": ()}, "free network"),
        ({"pin_faces": ("none",)}, "unknown far-field"),
        ({"n_fibers": 0}, "positive integer"),
        ({"n_fibers": 1.5}, "positive integer"),
        ({"target_segment_um": 0.0}, "finite and positive"),
        ({"box_hi_um": (-5.0, 4.0, 3.0)}, "must exceed"),
        ({"box_hi_um": (5.0, 4.0)}, "three box coordinates"),
        ({"rng_seed": 2**31}, "int32"),
        ({"rng_seed": 1.5}, "int32"),
        ({"max_refinement_level": -1}, "nonnegative integer"),
        ({"persistent_id_base": 1.5}, "must be an integer"),
    ],
)
def test_config_rejects_nonphysical_or_unaddressable_initialization(
    change: dict[str, object],
    message: str,
) -> None:
    with pytest.raises(ValueError, match=message):
        _config(**change)


def test_runtime_package_has_no_host_topology_or_coupled_core_dependency() -> None:
    forbidden_roots = {
        "numpy",
        "scipy",
        "aleph.laws.ecm_library",
        "aleph.laws.ecm_mikado",
        "aleph.components.incumbent",
        "aleph.engine.ecm_world",
    }
    violations: list[str] = []
    for path in ECM_ROOT.glob("*.py"):
        tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        for node in ast.walk(tree):
            names: list[str] = []
            if isinstance(node, ast.Import):
                names = [alias.name for alias in node.names]
            elif isinstance(node, ast.ImportFrom) and node.module:
                names = [node.module]
            for name in names:
                if any(name == root or name.startswith(f"{root}.") for root in forbidden_roots):
                    violations.append(f"{path.name}:{node.lineno}: {name}")
    assert not violations, "forbidden ECM topology dependencies:\n" + "\n".join(violations)


def test_topology_source_carries_material_point_and_generation_contract() -> None:
    source = (ECM_ROOT / "mikado_topology.py").read_text(encoding="utf-8")
    for token in (
        "out_segment_a",
        "out_segment_b",
        "out_u_a",
        "out_u_b",
        "out_generation_a",
        "out_generation_b",
        "_closest_segment_points",
    ):
        assert token in source
    assert "cKDTree" not in source
    assert "30--100 Pa versus 5--100 Pa" in source


def test_remodel_runtime_is_commit_predicated_and_has_no_device_to_host_read() -> None:
    source = (ECM_ROOT / "remodel_runtime.py").read_text(encoding="utf-8")
    for token in (
        "prepare_candidate",
        "validate_remap_acknowledgements",
        "candidate_valid_d",
        "topology_dirty_d",
        "refinement_root_segment_id_d",
        "refinement_path_d",
        "wp.utils.array_scan",
    ):
        assert token in source
    for forbidden in (".numpy(", ".to_list(", "wp.to_torch("):
        assert forbidden not in source
