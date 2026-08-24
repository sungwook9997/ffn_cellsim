"""Validated bridge from a READY AFM preflight manifest to the CUDA loading driver.

This module chooses no physical value.  It checks that one seed/speed loading path is fully present in the
preflight manifest, requires a zero-depth baseline followed by strictly increasing indentation, and derives
the initial spherical-indenter centre from the actual accepted membrane mesh.  The centre derivation is the
lowest axial position that leaves every discrete membrane vertex at or beyond the declared centre-contact
distance; therefore depth zero is structurally force-free without a fitted clearance.
"""

from __future__ import annotations

import json
import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np

__all__ = ["ReadyAFMLoadingPath", "derive_force_free_initial_centre", "load_ready_afm_path"]


@dataclass(frozen=True, slots=True)
class ReadyAFMLoadingPath:
    """One monotone speed/seed path selected from a fail-closed preflight manifest."""

    seed: int
    speed_um_s: float
    depths_um: tuple[float, ...]
    indenter_radius_um: float
    surface_contact_gap_um: float
    contact_stiffness_pn_per_um: float
    exterior_viscosity_pa_s: float
    pi0_pa: float
    result_class: str
    apparatus_source: str
    pi0_source: str
    contact_source: str


def _positive(value: object, label: str) -> float:
    number = float(value)
    if not math.isfinite(number) or number <= 0.0:
        raise ValueError(f"{label} must be finite and positive")
    return number


def _source(value: object, label: str) -> str:
    source = str(value).strip()
    if not source or source == "None":
        raise ValueError(f"{label} provenance/source is required")
    return source


def load_ready_afm_path(path: Path, *, seed: int, speed_um_s: float) -> ReadyAFMLoadingPath:
    """Load exactly one monotone loading path; a BLOCKED manifest can never reach CUDA."""
    payload = json.loads(path.read_text(encoding="utf-8"))
    if payload.get("schema") != "afm-sweep-preflight@1":
        raise ValueError("unsupported AFM preflight schema")
    if payload.get("status") != "READY_FOR_CUDA_DRIVER" or payload.get("blockers"):
        raise RuntimeError("AFM CUDA launch requires a READY preflight manifest with no blockers")
    result_class = str(payload.get("result_class_requested", ""))
    if result_class != "mechanism_demo_not_quantitative":
        raise RuntimeError(
            "quantitative AFM launch remains blocked until a post-induction, equilibrated, "
            "dt-converged active-cortex state handoff closes STATE (c) 3/17"
        )
    protocol = payload.get("protocol")
    pi0 = payload.get("pi0")
    contact = payload.get("contact")
    if not isinstance(protocol, dict) or not isinstance(pi0, dict) or not isinstance(contact, dict):
        raise RuntimeError("READY AFM manifest must carry protocol, Pi_0 and contact provenance blocks")
    if (
        protocol.get("geometry") != "suspended_round"
        or protocol.get("holder") != "exterior_stokes_six_rigid_modes"
    ):
        raise RuntimeError("AFM loading requires the declared suspended-round/exterior-Stokes protocol")
    selected_speed = _positive(speed_um_s, "loading speed [um/s]")
    selected_seed = int(seed)
    if selected_seed < 0:
        raise ValueError("seed must be nonnegative")
    rows = [
        row for row in payload.get("points", [])
        if int(row["seed"]) == selected_seed and float(row["speed_um_s"]) == selected_speed
    ]
    depths = tuple(sorted({float(row["depth_um"]) for row in rows}))
    if not depths or depths[0] != 0.0 or any(
        not math.isfinite(depth) or depth < 0.0 for depth in depths
    ):
        raise ValueError("an AFM loading path must contain a finite zero-depth baseline")
    if len(depths) < 2:
        raise ValueError("an AFM loading path needs zero plus at least one positive depth")
    expected = {
        (selected_speed, depth, selected_seed)
        for depth in (float(value) for value in payload["axes"]["depth_um"])
    }
    observed = {
        (float(row["speed_um_s"]), float(row["depth_um"]), int(row["seed"]))
        for row in rows
    }
    if observed != expected:
        raise ValueError("selected AFM path is incomplete or contains duplicate/mismatched points")
    return ReadyAFMLoadingPath(
        seed=selected_seed,
        speed_um_s=selected_speed,
        depths_um=depths,
        indenter_radius_um=_positive(protocol["indenter_radius_um"], "indenter radius [um]"),
        surface_contact_gap_um=_positive(contact["contact_distance_um"], "surface contact gap [um]"),
        contact_stiffness_pn_per_um=_positive(
            contact["stiffness_pn_per_um"], "contact stiffness [pN/um]"
        ),
        exterior_viscosity_pa_s=_positive(
            protocol["exterior_viscosity_pa_s"], "exterior viscosity [Pa*s]"
        ),
        pi0_pa=_positive(pi0["value_pa"], "Pi_0 [Pa]"),
        result_class=result_class,
        apparatus_source=_source(protocol["apparatus_source"], "apparatus"),
        pi0_source=_source(pi0["source"], "Pi_0"),
        contact_source=_source(contact["provenance"], "contact calibration"),
    )


def derive_force_free_initial_centre(
    membrane_position_um: np.ndarray,
    *,
    centre_contact_distance_um: float,
) -> tuple[float, float, float]:
    """Derive the lowest axial sphere centre whose discrete gaps are all nonnegative."""
    position = np.asarray(membrane_position_um, dtype=np.float64)
    rest = _positive(centre_contact_distance_um, "centre contact distance [um]")
    if position.ndim != 2 or position.shape[1] != 3 or position.shape[0] < 1:
        raise ValueError("membrane position must have shape (N, 3) with N>=1")
    if not np.isfinite(position).all():
        raise ValueError("membrane position must be finite")
    centre_xy = position[:, :2].mean(axis=0)
    radial_sq = np.sum((position[:, :2] - centre_xy) ** 2, axis=1)
    eligible = radial_sq <= rest * rest
    if not np.any(eligible):
        raise ValueError("no membrane vertex lies within the indenter's axial contact footprint")
    clearance = np.sqrt(np.maximum(0.0, rest * rest - radial_sq[eligible]))
    centre_z = float(np.max(position[eligible, 2] + clearance))
    centre = np.array([centre_xy[0], centre_xy[1], centre_z], dtype=np.float64)
    distances = np.linalg.norm(position - centre, axis=1)
    roundoff = 64.0 * np.finfo(np.float64).eps * max(1.0, rest, abs(centre_z))
    if float(distances.min()) < rest - roundoff:
        raise RuntimeError("derived indenter centre overlaps the accepted membrane at depth zero")
    return tuple(float(value) for value in centre)
