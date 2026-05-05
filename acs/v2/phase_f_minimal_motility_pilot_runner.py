"""Phase F minimal cell motility pilot runner.

Sister with :mod:`acs.v2.phase_e_v2_pilot_runner`. This runner repeatedly
calls the locked :func:`acs.v2.dynamics.phase_f_minimal_motility.step_phase_f_minimal_motility`
on a synthetic single-cell + multi-FA fixture. Per step it:

- updates the active contour via FA-traction Newton's 3rd law reaction
  (locked Phase F minimal pilot at commit ``bf03e83`` /
  ``docs/v2_phase_f_minimal_cell_motility_pilot_locked.md``)
- runs the Phase E v2 wrapper for ECM-side diagnostics in the loop
  (HB#4 multipliers strictly diagnostic-only)
- writes one HDF5 frame snapshot of the cell + ECM cluster
- records per-step CSV diagnostics (centroid, attached/non-attached
  vertex displacement, contour perimeter / area, multiplier summary)
- persists JSON run metadata with explicit honest-scope notes

Display-only artifacts: this runner is plumbing over locked physics.
It introduces no new dynamics, no new measurement metric, and no
PI-data fitting path. The runner outputs are visualizations of an
existing locked composition, NOT validation evidence.

This pilot shows FA-traction-driven contour motion with Phase E v2
ECM-side diagnostics running in the loop. It does not yet implement
ECM->cell motility feedback, because HB#4 rate multipliers are
diagnostic-only here.
"""

from __future__ import annotations

import csv
import json
import os
import time
from dataclasses import asdict, dataclass
from typing import Optional

import numpy as np

from acs.v2.active_contour import (
    ActiveContourParameters,
    ActiveContourState,
    regular_polygon_vertices,
)
from acs.v2.cell_cluster import CellClusterState
from acs.v2.dynamics.phase_f_minimal_motility import (
    PhaseFStepResult,
    step_phase_f_minimal_motility,
)
from acs.v2.ecm_substrate import ECMSubstrateState
from acs.v2.focal_adhesion import FocalAdhesionState
from acs.v2.measurement_boundary import MeasurementBoundary
from acs.v2.output.frame_dump import write_frame
from acs.v2.single_cell import SingleCellState


@dataclass(frozen=True, slots=True)
class PhaseFPilotConfig:
    """Configuration for the synthetic Phase F minimal motility pilot runner.

    Defaults produce a polygon with ~visible motion under a single FA
    pulling along ``-x`` (so the cell vertex receives a Newton 3rd
    reaction along ``+x``). All values are smoke-run inputs, NOT
    biological defaults.
    """

    n_steps: int = 200
    dt_cell_s: float = 1e-3
    frame_interval: int = 20
    grid_n: int = 8
    spacing_um: float = 0.5
    n_vertices: int = 16
    radius_um: float = 1.0
    cell_center_um: tuple[float, float] = (2.0, 2.0)
    lambda_c_nN: float = 0.05
    sigma_c_nN_per_um: float = 0.05
    k_a_nN_per_um: float = 0.0
    target_area_um2_factor: float = 1.0
    xi_line_nN_s_per_um2: float = 1.0
    effective_height_um: float = 1.0
    k_active: float = 1.0
    fa_traction_nN: float = 0.5
    fa_count: int = 1

    def validate(self) -> None:
        _require_positive_int(self.n_steps, "n_steps")
        _require_positive_int(self.frame_interval, "frame_interval")
        _require_positive_int(self.grid_n, "grid_n")
        _require_positive_int(self.n_vertices, "n_vertices")
        _require_positive_int(self.fa_count, "fa_count")
        if self.fa_count > self.n_vertices:
            raise ValueError(
                f"fa_count={self.fa_count} cannot exceed n_vertices={self.n_vertices}"
            )
        for name, value in (
            ("dt_cell_s", self.dt_cell_s),
            ("spacing_um", self.spacing_um),
            ("radius_um", self.radius_um),
            ("lambda_c_nN", self.lambda_c_nN),
            ("sigma_c_nN_per_um", self.sigma_c_nN_per_um),
            ("k_a_nN_per_um", self.k_a_nN_per_um),
            ("target_area_um2_factor", self.target_area_um2_factor),
            ("xi_line_nN_s_per_um2", self.xi_line_nN_s_per_um2),
            ("effective_height_um", self.effective_height_um),
            ("k_active", self.k_active),
            ("fa_traction_nN", self.fa_traction_nN),
        ):
            if isinstance(value, bool) or not np.isfinite(value):
                raise ValueError(f"{name} must be finite numeric, got {value!r}")
        if self.dt_cell_s <= 0.0:
            raise ValueError(f"dt_cell_s must be positive, got {self.dt_cell_s!r}")
        if self.spacing_um <= 0.0:
            raise ValueError(f"spacing_um must be positive, got {self.spacing_um!r}")
        if self.radius_um <= 0.0:
            raise ValueError(f"radius_um must be positive, got {self.radius_um!r}")
        if self.k_active <= 0.0:
            raise ValueError(f"k_active must be positive, got {self.k_active!r}")
        if self.xi_line_nN_s_per_um2 <= 0.0:
            raise ValueError(
                f"xi_line_nN_s_per_um2 must be positive, "
                f"got {self.xi_line_nN_s_per_um2!r}"
            )
        if self.effective_height_um <= 0.0:
            raise ValueError(
                f"effective_height_um must be positive, "
                f"got {self.effective_height_um!r}"
            )
        if self.target_area_um2_factor <= 0.0:
            raise ValueError(
                f"target_area_um2_factor must be positive, "
                f"got {self.target_area_um2_factor!r}"
            )


@dataclass(frozen=True, slots=True)
class PhaseFPilotStepDiagnostics:
    """One row of runner diagnostics after a Phase F step."""

    step_index: int
    time_s: float
    centroid_x_um: float
    centroid_y_um: float
    attached_vertex_displacement_um: float
    max_non_attached_vertex_displacement_um: float
    perimeter_um: float
    area_um2: float
    v_active_um2: float
    multiplier_min: float
    multiplier_max: float
    multiplier_mean: float


@dataclass(frozen=True, slots=True)
class PhaseFPilotRunResult:
    """Artifacts and final state from one Phase F pilot run."""

    output_dir: str
    config: PhaseFPilotConfig
    diagnostics: tuple[PhaseFPilotStepDiagnostics, ...]
    frame_paths: tuple[str, ...]
    csv_path: str
    metadata_path: str
    final_contour: ActiveContourState
    final_phase_e_v2: Optional["object"]
    wall_clock_s: float


def _require_positive_int(value: int, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise TypeError(f"{name} must be int, got {type(value).__name__}")
    if value <= 0:
        raise ValueError(f"{name} must be positive, got {value!r}")
    return value


def make_phase_f_pilot_fixture(
    config: PhaseFPilotConfig,
) -> tuple[
    ActiveContourState,
    ECMSubstrateState,
    tuple[FocalAdhesionState, ...],
    dict[str, int],
]:
    """Build the synthetic Phase F minimal pilot fixture.

    Returns (contour_state, ecm, adhesions, fa_to_vertex_index). The
    contour is a regular polygon centered at ``cell_center_um``. FAs
    are placed at the leftmost vertex (single-FA default) or evenly
    distributed for multi-FA configs. Each FA pulls along ``-x``
    (cell-on-substrate convention: substrate pulls cell along ``+x``
    via Newton's 3rd law -> visible +x cell motion).
    """

    config.validate()

    # P1 active contour state
    target_area = float(np.pi * config.radius_um ** 2 * config.target_area_um2_factor)
    params = ActiveContourParameters(
        k_a_nN_per_um=config.k_a_nN_per_um,
        target_area_um2=target_area,
        dt_cell_s=config.dt_cell_s,
        n_vertices=config.n_vertices,
        lambda_c_nN=config.lambda_c_nN,
        sigma_c_nN_per_um=config.sigma_c_nN_per_um,
        xi_line_nN_s_per_um2=config.xi_line_nN_s_per_um2,
        effective_height_um=config.effective_height_um,
    ).validate()
    vertices = regular_polygon_vertices(
        config.n_vertices, config.radius_um, center=config.cell_center_um
    )
    contour = ActiveContourState(
        cell_id="phase-f-cell",
        params=params,
        vertices_xy_um=vertices,
    )

    # ECM with isotropic 0.5*I orientation IC (HB#4-active automatic-neutral)
    orientation = np.zeros(
        (config.grid_n, config.grid_n, 2, 2), dtype=np.float64
    )
    orientation[..., 0, 0] = 0.5
    orientation[..., 1, 1] = 0.5
    ecm = ECMSubstrateState(
        origin_um_xy=(0.0, 0.0),
        spacing_um=config.spacing_um,
        stiffness_kpa=np.zeros((config.grid_n, config.grid_n), dtype=np.float64),
        ligand_density=np.zeros((config.grid_n, config.grid_n), dtype=np.float64),
        fiber_density=np.zeros((config.grid_n, config.grid_n), dtype=np.float64),
        orientation_tensor=orientation,
        accumulated_traction_nNs_per_um2=np.zeros(
            (config.grid_n, config.grid_n), dtype=np.float64
        ),
        source="synthetic",
    )

    # FAs at leftmost vertex(es), pulling -x. Newton 3rd reaction -> +x cell motion.
    # Sort vertices by x-coordinate; pick the leftmost fa_count.
    x_order = np.argsort(vertices[:, 0])
    fa_vertex_indices = [int(x_order[i]) for i in range(config.fa_count)]
    adhesions: list[FocalAdhesionState] = []
    fa_to_vertex_index: dict[str, int] = {}
    for k, v_idx in enumerate(fa_vertex_indices):
        fa_id = f"phase-f-fa-{k}"
        fa = FocalAdhesionState(
            adhesion_id=fa_id,
            cell_id="phase-f-cell",
            position_um_xy=(
                float(vertices[v_idx, 0]),
                float(vertices[v_idx, 1]),
            ),
            age_s=0.0,
            maturity=1.0,
            bound_fraction=1.0,
            state="mature",
            traction_force_nN_xy=(-config.fa_traction_nN, 0.0),
        )
        adhesions.append(fa)
        fa_to_vertex_index[fa_id] = v_idx

    return contour, ecm, tuple(adhesions), fa_to_vertex_index


def _cluster_from_contour(
    contour: ActiveContourState,
    adhesions: tuple[FocalAdhesionState, ...],
    ecm: ECMSubstrateState,
    *,
    time_s: float,
) -> CellClusterState:
    """Build a CellClusterState snapshot from current contour + FAs."""
    boundary = MeasurementBoundary.from_array(
        np.asarray(contour.vertices_xy_um, dtype=np.float64),
        coordinate_convention="world_um_y_up",
        source_modality="phase_f_minimal_motility_pilot_synthetic",
        object_id=contour.cell_id,
    )
    cell = SingleCellState(
        cell_id=contour.cell_id,
        time_s=float(time_s),
        measurement_boundary=boundary,
        adhesions=list(adhesions),
        cell_state="alive",
    )
    cell.validate()
    return CellClusterState(cells={cell.cell_id: cell}, ecm=ecm)


def _write_frame(
    output_dir: str,
    contour: ActiveContourState,
    adhesions: tuple[FocalAdhesionState, ...],
    ecm: ECMSubstrateState,
    *,
    step_index: int,
    time_s: float,
) -> str:
    cluster = _cluster_from_contour(contour, adhesions, ecm, time_s=time_s)
    path = os.path.join(output_dir, f"frame_{step_index:06d}.h5")
    write_frame(path, cluster, time_s=float(time_s))
    return path


def _compute_diagnostics(
    step_index: int,
    time_s: float,
    contour_before: ActiveContourState,
    result: PhaseFStepResult,
    fa_to_vertex_index: dict[str, int],
) -> PhaseFPilotStepDiagnostics:
    contour_after = result.contour
    centroid = contour_after.vertices_xy_um.mean(axis=0)
    displacement = (
        contour_after.vertices_xy_um - contour_before.vertices_xy_um
    )
    displacement_magnitude = np.linalg.norm(displacement, axis=-1)
    attached = sorted(set(fa_to_vertex_index.values()))
    non_attached = [
        i for i in range(contour_after.vertices_xy_um.shape[0]) if i not in attached
    ]
    attached_max = (
        float(displacement_magnitude[attached].max()) if attached else 0.0
    )
    non_attached_max = (
        float(displacement_magnitude[non_attached].max())
        if non_attached
        else 0.0
    )
    multipliers = np.asarray(
        result.phase_e_v2.ecm_to_fa_bias.multipliers_per_fa, dtype=np.float64
    )
    return PhaseFPilotStepDiagnostics(
        step_index=int(step_index),
        time_s=float(time_s),
        centroid_x_um=float(centroid[0]),
        centroid_y_um=float(centroid[1]),
        attached_vertex_displacement_um=attached_max,
        max_non_attached_vertex_displacement_um=non_attached_max,
        perimeter_um=float(contour_after.perimeter_um()),
        area_um2=float(contour_after.area_um2()),
        v_active_um2=float(result.phase_e_v2.lyapunov_metric.v_active_um2),
        multiplier_min=float(multipliers.min()) if multipliers.size else 1.0,
        multiplier_max=float(multipliers.max()) if multipliers.size else 1.0,
        multiplier_mean=float(multipliers.mean()) if multipliers.size else 1.0,
    )


def _write_diagnostics_csv(
    output_dir: str, diagnostics: tuple[PhaseFPilotStepDiagnostics, ...]
) -> str:
    path = os.path.join(output_dir, "diagnostics.csv")
    fieldnames = [
        "step_index",
        "time_s",
        "centroid_x_um",
        "centroid_y_um",
        "attached_vertex_displacement_um",
        "max_non_attached_vertex_displacement_um",
        "perimeter_um",
        "area_um2",
        "v_active_um2",
        "multiplier_min",
        "multiplier_max",
        "multiplier_mean",
    ]
    with open(path, "w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        for row in diagnostics:
            writer.writerow(asdict(row))
    return path


def _write_metadata_json(
    output_dir: str,
    *,
    config: PhaseFPilotConfig,
    diagnostics: tuple[PhaseFPilotStepDiagnostics, ...],
    frame_paths: tuple[str, ...],
    git_commit_hash: str,
    wall_clock_s: float,
) -> str:
    final = diagnostics[-1] if diagnostics else None
    payload = {
        "runner": "phase_f_minimal_motility_pilot",
        "tier": "B",
        "contract": (
            "runner/plumbing over locked step_phase_f_minimal_motility "
            "(commit bf03e83)"
        ),
        "git_commit_hash": git_commit_hash or "unrecorded",
        "config": asdict(config),
        "n_diagnostic_rows": len(diagnostics),
        "n_frames": len(frame_paths),
        "frame_paths": [os.path.basename(p) for p in frame_paths],
        "csv_path": "diagnostics.csv",
        "wall_clock_s": float(wall_clock_s),
        "final": asdict(final) if final is not None else None,
        "visual_smoke_only": False,
        "not_mechanistic": True,
        "physics_layer": "FA-vertex Newton 3rd law (locked)",
        "ecm_to_cell_feedback": False,
        "notes": [
            "synthetic single-cell + multi-FA fixture",
            "FA traction drives attached cell vertex via Newton 3rd law",
            "P1 active contour cortex+area dynamics on all vertices",
            "Phase E v2 wrapper runs ECM-side diagnostics in the loop",
            "HB#4 rate multipliers are diagnostic-only (NOT used for force scaling)",
            "no ECM->cell motility feedback (separate Phase F proper A-tier cycle)",
            "no FA dynamics (assembly/disassembly not modeled)",
            "no PI experimental data use (Hard Rule 1)",
            (
                "this pilot shows FA-traction-driven contour motion with "
                "Phase E v2 ECM-side diagnostics running in the loop; it "
                "does not yet implement ECM->cell motility feedback, "
                "because HB#4 rate multipliers are diagnostic-only here"
            ),
        ],
    }
    path = os.path.join(output_dir, "metadata.json")
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(payload, fh, indent=2)
    return path


def run_phase_f_minimal_motility_pilot(
    output_dir: str,
    config: PhaseFPilotConfig | None = None,
    *,
    git_commit_hash: str = "",
) -> PhaseFPilotRunResult:
    """Run the synthetic Phase F minimal cell motility pilot."""

    cfg = config or PhaseFPilotConfig()
    cfg.validate()
    os.makedirs(output_dir, exist_ok=True)

    contour, ecm, adhesions, fa_to_vertex_index = make_phase_f_pilot_fixture(cfg)
    diagnostics: list[PhaseFPilotStepDiagnostics] = []
    frame_paths: list[str] = []
    final_phase_e_v2 = None
    start = time.perf_counter()

    # Initial frame at t=0
    frame_paths.append(
        _write_frame(
            output_dir, contour, adhesions, ecm, step_index=0, time_s=0.0
        )
    )

    for step_index in range(1, cfg.n_steps + 1):
        contour_before = contour
        result = step_phase_f_minimal_motility(
            contour,
            adhesions,
            ecm,
            fa_to_vertex_index,
            k_active=cfg.k_active,
        )
        contour = result.contour
        adhesions = result.adhesions
        ecm = result.updated_ecm
        final_phase_e_v2 = result.phase_e_v2
        time_s = float(step_index) * float(cfg.dt_cell_s)
        diagnostics.append(
            _compute_diagnostics(
                step_index, time_s, contour_before, result, fa_to_vertex_index
            )
        )
        if step_index % cfg.frame_interval == 0 or step_index == cfg.n_steps:
            frame_paths.append(
                _write_frame(
                    output_dir,
                    contour,
                    adhesions,
                    ecm,
                    step_index=step_index,
                    time_s=time_s,
                )
            )

    wall_clock_s = time.perf_counter() - start
    diagnostics_tuple = tuple(diagnostics)
    frame_paths_tuple = tuple(frame_paths)
    csv_path = _write_diagnostics_csv(output_dir, diagnostics_tuple)
    metadata_path = _write_metadata_json(
        output_dir,
        config=cfg,
        diagnostics=diagnostics_tuple,
        frame_paths=frame_paths_tuple,
        git_commit_hash=git_commit_hash,
        wall_clock_s=wall_clock_s,
    )
    return PhaseFPilotRunResult(
        output_dir=output_dir,
        config=cfg,
        diagnostics=diagnostics_tuple,
        frame_paths=frame_paths_tuple,
        csv_path=csv_path,
        metadata_path=metadata_path,
        final_contour=contour,
        final_phase_e_v2=final_phase_e_v2,
        wall_clock_s=wall_clock_s,
    )
