"""Production-like Phase E v2 pilot runner.

This module is B-tier runner/plumbing. It repeatedly calls the locked
``step_closed_loop_phase_e_v2`` wrapper on a synthetic one-cell / two-FA
fixture, writes HDF5 frame dumps, records per-step CSV diagnostics, and
persists JSON metadata. It introduces no new dynamics, no new measurement
metric, and no PI-data fitting path.

The synthetic fixture parameters are smoke-run inputs only. They are not model
defaults and must not be interpreted as biological rates or fitted constants.
"""

from __future__ import annotations

import csv
import json
import os
import time
from dataclasses import asdict, dataclass
from typing import Literal, Optional

import numpy as np

from acs.v2.cell_cluster import CellClusterState
from acs.v2.dynamics.closed_loop_phase_e import (
    PhaseEStepResult,
    step_closed_loop_phase_e_v2,
)
from acs.v2.dynamics.ecm_constitutive_response import (
    K_ORIENT_PER_S,
    TRACTION_REF_NN_PER_UM2,
)
from acs.v2.ecm_substrate import ECMSubstrateState
from acs.v2.focal_adhesion import FocalAdhesionState
from acs.v2.measurement_boundary import MeasurementBoundary
from acs.v2.output.frame_dump import write_frame
from acs.v2.single_cell import SingleCellState


@dataclass(frozen=True, slots=True)
class PhaseEV2PilotConfig:
    """Configuration for the synthetic Phase E v2 pilot runner.

    The optional ``cell_motility`` field enables an explicit
    *visualization-only* smoke mode. When set to ``"translation_smoke"``,
    the runner switches to a single-FA fixture (so net traction is
    nonzero) and translates the cell boundary by
    ``smoke_displacement_um_per_step`` along the FA traction direction
    each step. This is **NOT** a biological motility model and the field
    name uses ``smoke_displacement`` (not ``motility_speed``) to avoid
    suggesting a physical velocity. See Codex ``id=1906``/``id=1912``
    A-modified guardrails: smoke evidence only, no mechanistic claim.

    Default: ``cell_motility="off"`` preserves the previous Phase E v2
    pilot runner behavior bit-for-bit.
    """

    n_steps: int = 4
    dt_s: float = 0.0
    frame_interval: int = 1
    grid_n: int = 4
    spacing_um: float = 1.0
    k_active: float = 1.0
    traction_ref_nN_per_um2: float = TRACTION_REF_NN_PER_UM2
    k_orient_per_s: float = K_ORIENT_PER_S
    cell_motility: Literal["off", "translation_smoke"] = "off"
    smoke_displacement_um_per_step: float = 0.0

    def validate(self) -> None:
        _require_non_negative_int(self.n_steps, "n_steps")
        _require_positive_int(self.frame_interval, "frame_interval")
        _require_positive_int(self.grid_n, "grid_n")
        for name, value in (
            ("dt_s", self.dt_s),
            ("spacing_um", self.spacing_um),
            ("k_active", self.k_active),
            ("traction_ref_nN_per_um2", self.traction_ref_nN_per_um2),
            ("k_orient_per_s", self.k_orient_per_s),
            ("smoke_displacement_um_per_step", self.smoke_displacement_um_per_step),
        ):
            if isinstance(value, bool) or not np.isfinite(value):
                raise ValueError(f"{name} must be finite numeric, got {value!r}")
        if self.dt_s < 0.0:
            raise ValueError(f"dt_s must be non-negative, got {self.dt_s!r}")
        if self.spacing_um <= 0.0:
            raise ValueError(f"spacing_um must be positive, got {self.spacing_um!r}")
        if self.k_active <= 0.0:
            raise ValueError(f"k_active must be positive, got {self.k_active!r}")
        if self.traction_ref_nN_per_um2 <= 0.0:
            raise ValueError(
                "traction_ref_nN_per_um2 must be positive, "
                f"got {self.traction_ref_nN_per_um2!r}"
            )
        if self.k_orient_per_s <= 0.0:
            raise ValueError(f"k_orient_per_s must be positive, got {self.k_orient_per_s!r}")
        if self.cell_motility not in ("off", "translation_smoke"):
            raise ValueError(
                f"cell_motility must be 'off' or 'translation_smoke', "
                f"got {self.cell_motility!r}"
            )
        if self.smoke_displacement_um_per_step < 0.0:
            raise ValueError(
                f"smoke_displacement_um_per_step must be non-negative, "
                f"got {self.smoke_displacement_um_per_step!r}"
            )
        if self.cell_motility == "translation_smoke" and self.smoke_displacement_um_per_step == 0.0:
            raise ValueError(
                "cell_motility='translation_smoke' requires "
                "smoke_displacement_um_per_step > 0 to produce visible motion"
            )


@dataclass(frozen=True, slots=True)
class PhaseEV2PilotStepDiagnostics:
    """One row of runner diagnostics after a Phase E v2 step."""

    step_index: int
    time_s: float
    v_active_um2: float
    multiplier_min: float
    multiplier_max: float
    multiplier_mean: float
    max_traction_norm_nN_per_um2: float
    orientation_tensor_abs_max: float


@dataclass(frozen=True, slots=True)
class PhaseEV2PilotRunResult:
    """Artifacts and final state from one pilot run."""

    output_dir: str
    config: PhaseEV2PilotConfig
    diagnostics: tuple[PhaseEV2PilotStepDiagnostics, ...]
    frame_paths: tuple[str, ...]
    csv_path: str
    metadata_path: str
    final_ecm: ECMSubstrateState
    final_phase_e_result: Optional[PhaseEStepResult]
    wall_clock_s: float


def _require_non_negative_int(value: int, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise TypeError(f"{name} must be int, got {type(value).__name__}")
    if value < 0:
        raise ValueError(f"{name} must be non-negative, got {value!r}")
    return value


def _require_positive_int(value: int, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise TypeError(f"{name} must be int, got {type(value).__name__}")
    if value <= 0:
        raise ValueError(f"{name} must be positive, got {value!r}")
    return value


def make_phase_e_v2_pilot_fixture(
    *, grid_n: int = 4, spacing_um: float = 1.0, motility_smoke: bool = False
) -> tuple[ECMSubstrateState, SingleCellState, tuple[FocalAdhesionState, ...]]:
    """Create the synthetic one-cell/two-FA fixture used by the runner.

    The ECM orientation tensor is traceless anisotropic everywhere:
    ``[[1, 0], [0, -1]]``. Two FAs share one cell-centered position and
    pull along perpendicular directions, matching the Phase E v2
    anti-collapse fixture.

    When ``motility_smoke=True``, switches to a **single-FA** variant
    pulling along ``(1, 0)`` only, so net traction is nonzero and the
    visualization-only smoke mode (Codex ``id=1906`` A-modified) shows a
    clear cell centroid translation. This variant is for *visual smoke
    only* and is **NOT** a biological motility setup.
    """

    _require_positive_int(grid_n, "grid_n")
    if not np.isfinite(spacing_um) or spacing_um <= 0.0:
        raise ValueError(f"spacing_um must be positive finite, got {spacing_um!r}")

    orientation = np.zeros((grid_n, grid_n, 2, 2), dtype=np.float64)
    orientation[..., 0, 0] = 1.0
    orientation[..., 1, 1] = -1.0
    ecm = ECMSubstrateState(
        origin_um_xy=(0.0, 0.0),
        spacing_um=float(spacing_um),
        stiffness_kpa=np.zeros((grid_n, grid_n), dtype=np.float64),
        ligand_density=np.zeros((grid_n, grid_n), dtype=np.float64),
        fiber_density=np.zeros((grid_n, grid_n), dtype=np.float64),
        orientation_tensor=orientation,
        accumulated_traction_nNs_per_um2=np.zeros((grid_n, grid_n), dtype=np.float64),
        source="synthetic",
    )
    ecm.validate()

    center_idx = min(max(grid_n // 2, 0), grid_n - 1)
    center = (center_idx + 0.5) * float(spacing_um)
    half = 0.35 * float(spacing_um)
    vertices = np.array(
        [
            [center - half, center - half],
            [center + half, center - half],
            [center + half, center + half],
            [center - half, center + half],
        ],
        dtype=np.float64,
    )
    boundary = MeasurementBoundary.from_array(
        vertices,
        coordinate_convention="world_um_y_up",
        source_modality="phase_e_v2_pilot_synthetic",
        object_id="pilot-cell",
    )
    if motility_smoke:
        # Single-FA variant pulling along (1, 0) so net traction is
        # nonzero — visualization-only smoke fixture (Codex id=1906).
        fas = (
            FocalAdhesionState(
                adhesion_id="pilot-fa-x",
                cell_id="pilot-cell",
                position_um_xy=(center, center),
                age_s=0.0,
                maturity=1.0,
                bound_fraction=1.0,
                state="mature",
                traction_force_nN_xy=(1.0, 0.0),
            ),
        )
    else:
        # Default: two perpendicular FAs at the same position (anti-collapse
        # fixture; net traction is zero, so the smoke mode would do nothing).
        fas = (
            FocalAdhesionState(
                adhesion_id="pilot-fa-x",
                cell_id="pilot-cell",
                position_um_xy=(center, center),
                age_s=0.0,
                maturity=1.0,
                bound_fraction=1.0,
                state="mature",
                traction_force_nN_xy=(1.0, 0.0),
            ),
            FocalAdhesionState(
                adhesion_id="pilot-fa-y",
                cell_id="pilot-cell",
                position_um_xy=(center, center),
                age_s=0.0,
                maturity=1.0,
                bound_fraction=1.0,
                state="mature",
                traction_force_nN_xy=(0.0, 1.0),
            ),
        )
    cell = SingleCellState(
        cell_id="pilot-cell",
        time_s=0.0,
        measurement_boundary=boundary,
        adhesions=list(fas),
        cell_state="alive",
    )
    cell.validate()
    return ecm, cell, fas


def _cell_at_time(template: SingleCellState, time_s: float) -> SingleCellState:
    return SingleCellState(
        cell_id=template.cell_id,
        time_s=float(time_s),
        measurement_boundary=template.measurement_boundary,
        height_um=template.height_um,
        polarity_xy=template.polarity_xy,
        protrusions=list(template.protrusions),
        adhesions=list(template.adhesions),
        cell_state=template.cell_state,
        cell_age_s=template.cell_age_s + float(time_s),
        cell_cycle_phase=template.cell_cycle_phase,
        division_count=template.division_count,
        parent_cell_id=template.parent_cell_id,
        mechanosignal_yap_taz=template.mechanosignal_yap_taz,
        neighbor_cell_ids=template.neighbor_cell_ids,
    )


def _apply_motility_smoke(
    cell: SingleCellState,
    adhesions: tuple[FocalAdhesionState, ...],
    *,
    smoke_displacement_um: float,
) -> SingleCellState:
    """Translate the cell boundary along the net traction direction.

    **Visual smoke only**: this is a display-time translation, NOT a
    biological motility model. The displacement parameter is named
    ``smoke_displacement`` (not ``motility_speed``) per Codex
    ``id=1906`` A-modified guardrail. The function does NOT update FA
    positions, ECM, or any physics state — it only shifts the
    measurement boundary so a viewer can see the cell move.
    """

    if not adhesions:
        return cell
    net = np.zeros(2, dtype=np.float64)
    for fa in adhesions:
        net += np.asarray(fa.traction_force_nN_xy, dtype=np.float64)
    norm = float(np.linalg.norm(net))
    if norm <= 0.0:
        # Net traction zero → no visible direction; skip translation
        # (smoke fixture should provide a single-FA setup with nonzero net).
        return cell
    direction = net / norm
    delta = direction * float(smoke_displacement_um)

    boundary = cell.measurement_boundary
    new_vertices = np.asarray(boundary.vertices_xy_um, dtype=np.float64) + delta
    new_boundary = MeasurementBoundary.from_array(
        new_vertices,
        coordinate_convention=boundary.coordinate_convention,
        source_modality=boundary.source_modality,
        object_id=boundary.object_id,
    )
    return SingleCellState(
        cell_id=cell.cell_id,
        time_s=cell.time_s,
        measurement_boundary=new_boundary,
        height_um=cell.height_um,
        polarity_xy=cell.polarity_xy,
        protrusions=list(cell.protrusions),
        adhesions=list(cell.adhesions),
        cell_state=cell.cell_state,
        cell_age_s=cell.cell_age_s,
        cell_cycle_phase=cell.cell_cycle_phase,
        division_count=cell.division_count,
        parent_cell_id=cell.parent_cell_id,
        mechanosignal_yap_taz=cell.mechanosignal_yap_taz,
        neighbor_cell_ids=cell.neighbor_cell_ids,
    )


def _write_cluster_frame(
    output_dir: str,
    ecm: ECMSubstrateState,
    cell_template: SingleCellState,
    *,
    step_index: int,
    time_s: float,
) -> str:
    cell = _cell_at_time(cell_template, time_s)
    cluster = CellClusterState(cells={cell.cell_id: cell}, ecm=ecm)
    path = os.path.join(output_dir, f"frame_{step_index:06d}.h5")
    write_frame(path, cluster, time_s=float(time_s))
    return path


def _diagnostics_from_result(
    step_index: int, time_s: float, result: PhaseEStepResult
) -> PhaseEV2PilotStepDiagnostics:
    multipliers = np.asarray(result.ecm_to_fa_bias.multipliers_per_fa, dtype=np.float64)
    traction = np.asarray(result.traction_density_xy, dtype=np.float64)
    traction_norm = np.linalg.norm(traction, axis=-1)
    orientation = np.asarray(result.updated_ecm.orientation_tensor, dtype=np.float64)
    return PhaseEV2PilotStepDiagnostics(
        step_index=int(step_index),
        time_s=float(time_s),
        v_active_um2=float(result.lyapunov_metric.v_active_um2),
        multiplier_min=float(multipliers.min()) if multipliers.size else 1.0,
        multiplier_max=float(multipliers.max()) if multipliers.size else 1.0,
        multiplier_mean=float(multipliers.mean()) if multipliers.size else 1.0,
        max_traction_norm_nN_per_um2=(
            float(traction_norm.max()) if traction_norm.size else 0.0
        ),
        orientation_tensor_abs_max=float(np.abs(orientation).max()) if orientation.size else 0.0,
    )


def _write_diagnostics_csv(
    output_dir: str, diagnostics: tuple[PhaseEV2PilotStepDiagnostics, ...]
) -> str:
    path = os.path.join(output_dir, "diagnostics.csv")
    fieldnames = [
        "step_index",
        "time_s",
        "v_active_um2",
        "multiplier_min",
        "multiplier_max",
        "multiplier_mean",
        "max_traction_norm_nN_per_um2",
        "orientation_tensor_abs_max",
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
    config: PhaseEV2PilotConfig,
    diagnostics: tuple[PhaseEV2PilotStepDiagnostics, ...],
    frame_paths: tuple[str, ...],
    git_commit_hash: str,
    wall_clock_s: float,
) -> str:
    final = diagnostics[-1] if diagnostics else None
    motility_notes: list[str] = []
    if config.cell_motility == "translation_smoke":
        motility_notes = [
            "cell_motility=translation_smoke (visual smoke only)",
            "smoke_displacement_um_per_step is a display parameter, NOT biological speed",
            "no FA-to-cell force coupling; boundary translated post-step for viewer",
        ]
    payload = {
        "runner": "phase_e_v2_pilot",
        "tier": "B",
        "contract": "runner/plumbing over locked step_closed_loop_phase_e_v2",
        "git_commit_hash": git_commit_hash or "unrecorded",
        "config": asdict(config),
        "n_diagnostic_rows": len(diagnostics),
        "n_frames": len(frame_paths),
        "frame_paths": [os.path.basename(p) for p in frame_paths],
        "csv_path": "diagnostics.csv",
        "wall_clock_s": float(wall_clock_s),
        "final": asdict(final) if final is not None else None,
        "visual_smoke_only": config.cell_motility == "translation_smoke",
        "not_mechanistic": True,
        "cell_motility_mode": config.cell_motility,
        "notes": [
            "synthetic one-cell/two-FA fixture (default) or single-FA (motility smoke)",
            "no FA dynamics",
            "no new physics law",
            "no new measurement semantics",
            "no PI experimental data use",
            *motility_notes,
        ],
    }
    path = os.path.join(output_dir, "metadata.json")
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(payload, fh, indent=2)
    return path


def run_phase_e_v2_pilot(
    output_dir: str,
    config: PhaseEV2PilotConfig | None = None,
    *,
    git_commit_hash: str = "",
) -> PhaseEV2PilotRunResult:
    """Run the synthetic Phase E v2 pilot and persist artifacts."""

    cfg = config or PhaseEV2PilotConfig()
    cfg.validate()
    os.makedirs(output_dir, exist_ok=True)

    motility_smoke = cfg.cell_motility == "translation_smoke"
    ecm, cell, adhesions = make_phase_e_v2_pilot_fixture(
        grid_n=cfg.grid_n,
        spacing_um=cfg.spacing_um,
        motility_smoke=motility_smoke,
    )
    diagnostics: list[PhaseEV2PilotStepDiagnostics] = []
    frame_paths: list[str] = []
    final_result: Optional[PhaseEStepResult] = None
    start = time.perf_counter()

    frame_paths.append(
        _write_cluster_frame(output_dir, ecm, cell, step_index=0, time_s=0.0)
    )

    for step_index in range(1, cfg.n_steps + 1):
        result = step_closed_loop_phase_e_v2(
            adhesions,
            ecm,
            cfg.dt_s,
            k_active=cfg.k_active,
            traction_ref_nN_per_um2=cfg.traction_ref_nN_per_um2,
            k_orient_per_s=cfg.k_orient_per_s,
        )
        ecm = result.updated_ecm
        final_result = result
        time_s = float(step_index) * float(cfg.dt_s)
        diagnostics.append(_diagnostics_from_result(step_index, time_s, result))
        if motility_smoke:
            cell = _apply_motility_smoke(
                cell,
                adhesions,
                smoke_displacement_um=cfg.smoke_displacement_um_per_step,
            )
        if step_index % cfg.frame_interval == 0 or step_index == cfg.n_steps:
            frame_paths.append(
                _write_cluster_frame(
                    output_dir, ecm, cell, step_index=step_index, time_s=time_s
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
    return PhaseEV2PilotRunResult(
        output_dir=output_dir,
        config=cfg,
        diagnostics=diagnostics_tuple,
        frame_paths=frame_paths_tuple,
        csv_path=csv_path,
        metadata_path=metadata_path,
        final_ecm=ecm,
        final_phase_e_result=final_result,
        wall_clock_s=wall_clock_s,
    )
