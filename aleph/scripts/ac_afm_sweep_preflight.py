#!/usr/bin/env python3
"""Fail-closed manifest builder for a suspended-round cortical-tension AFM sweep.

This host-side script launches no physics.  It validates the committed three-seed C-2 ensemble and
enumerates the physical speed/depth/seed grid only when every PI-owned apparatus and material input is
present with provenance.  A blocked invocation writes a blocker manifest with an empty point list and exits
2, so a scheduler wrapper cannot accidentally interpret incomplete inputs as a runnable experiment.
"""

from __future__ import annotations

import argparse
import json
from dataclasses import asdict
from pathlib import Path
from typing import Any

from aleph.engine.afm_cortical_tension import (
    AFMContactCalibration,
    AFMSweepDeclaration,
    C2Acceptance,
    CellGeometry,
    Holder,
    OsmoticSetpoint,
    Pi0Evidence,
    SphericalIndenterProtocol,
)


def load_c2_ensemble(path: Path) -> C2Acceptance:
    """Validate and reduce an ``afm-c2-accepted-ensemble@1`` record."""
    record = json.loads(path.read_text(encoding="utf-8"))
    if record.get("schema") != "afm-c2-accepted-ensemble@1":
        raise ValueError("unsupported C-2 ensemble schema")
    config = record["configuration"]
    rows = record["records"]
    declared_seeds = tuple(int(seed) for seed in record["accepted_seeds"])
    row_seeds = tuple(int(row["seed"]) for row in rows)
    if record.get("status") != "ACCEPTED_3_SEEDS" or row_seeds != declared_seeds:
        raise ValueError("C-2 ensemble status/seeds do not describe the stored records")
    if not rows:
        raise ValueError("C-2 ensemble contains no records")
    if any(
        float(row["max_projected_force_pn"])
        > float(row["projected_force_tolerance_pn"])
        for row in rows
    ):
        raise ValueError("C-2 ensemble contains a force-rejected record")
    return C2Acceptance(
        build_commit=str(record["declared_build_commit"]),
        inner_converged=True,
        outer_accepted=True,
        committed_time_s=float(record["committed_time_s_per_seed"]),
        max_projected_force_pn=max(float(row["max_projected_force_pn"]) for row in rows),
        projected_force_tolerance_pn=min(
            float(row["projected_force_tolerance_pn"]) for row in rows
        ),
        membrane_subdivisions=int(config["membrane_subdivisions"]),
        full_native_population=bool(config["full_native_population"]),
        full_compartments=bool(config["full_compartments"]),
        accepted_seeds=declared_seeds,
    )


def _all_or_none(values: tuple[object | None, ...], name: str) -> bool:
    present = tuple(value is not None for value in values)
    if any(present) and not all(present):
        raise ValueError(f"{name} fields must be supplied together")
    return all(present)


def build_manifest(args: argparse.Namespace) -> dict[str, Any]:
    """Build a blocked or runnable manifest without launching CUDA work."""
    protocol = SphericalIndenterProtocol(
        target_protocol=args.target_protocol,
        target_cell_line=args.target_cell_line,
        geometry=CellGeometry.SUSPENDED_ROUND,
        holder=Holder.EXTERIOR_STOKES,
        indenter_radius_um=args.indenter_radius_um,
        exterior_viscosity_pa_s=args.exterior_viscosity_pa_s,
        apparatus_source=args.apparatus_source,
    )
    pi0 = None
    if _all_or_none((args.pi0_pa, args.pi0_evidence, args.pi0_source), "Pi_0"):
        pi0 = OsmoticSetpoint(
            value_pa=args.pi0_pa,
            evidence=Pi0Evidence(args.pi0_evidence),
            source=args.pi0_source,
        )
    contact = None
    if _all_or_none(
        (args.contact_distance_um, args.contact_stiffness_pn_per_um, args.contact_source),
        "contact calibration",
    ):
        contact = AFMContactCalibration(
            contact_distance_um=args.contact_distance_um,
            stiffness_pn_per_um=args.contact_stiffness_pn_per_um,
            provenance=args.contact_source,
        )
    c2 = load_c2_ensemble(args.c2_ensemble)
    declaration = AFMSweepDeclaration(
        protocol=protocol,
        pi0=pi0,
        contact=contact,
        c2=c2,
        speeds_um_s=tuple(args.speed_um_s),
        depths_um=tuple(args.depth_um),
        seeds=tuple(args.seed),
    )
    require_quantitative = not args.allow_mechanism_demo
    blockers = declaration.blockers(require_quantitative=require_quantitative)
    points = () if blockers else declaration.points(require_quantitative=require_quantitative)
    return {
        "schema": "afm-sweep-preflight@1",
        "status": "BLOCKED" if blockers else "READY_FOR_CUDA_DRIVER",
        "result_class_requested": "mechanism_demo_not_quantitative"
        if args.allow_mechanism_demo
        else "quantitative",
        "reaction_channel": declaration.reaction_channel,
        "protocol": asdict(protocol),
        "pi0": None if pi0 is None else asdict(pi0),
        "contact": None if contact is None else asdict(contact),
        "c2": asdict(c2),
        "axes": {
            "speed_um_s": list(declaration.speeds_um_s),
            "depth_um": list(declaration.depths_um),
            "seed": list(declaration.seeds),
        },
        "blockers": list(blockers),
        "points": [asdict(point) for point in points],
        "scope": "Preflight only; READY_FOR_CUDA_DRIVER is not an executed force-indentation result.",
    }


def parser() -> argparse.ArgumentParser:
    """Return the command-line parser without assigning physical defaults."""
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--c2-ensemble", type=Path, required=True)
    result.add_argument("--target-protocol", required=True)
    result.add_argument("--target-cell-line", required=True)
    result.add_argument("--indenter-radius-um", type=float, required=True)
    result.add_argument("--exterior-viscosity-pa-s", type=float, required=True)
    result.add_argument("--apparatus-source", required=True)
    result.add_argument("--pi0-pa", type=float)
    result.add_argument("--pi0-evidence", choices=[value.value for value in Pi0Evidence])
    result.add_argument("--pi0-source")
    result.add_argument("--contact-distance-um", type=float)
    result.add_argument("--contact-stiffness-pn-per-um", type=float)
    result.add_argument("--contact-source")
    result.add_argument("--speed-um-s", action="append", type=float, required=True)
    result.add_argument("--depth-um", action="append", type=float, required=True)
    result.add_argument("--seed", action="append", type=int, required=True)
    result.add_argument("--allow-mechanism-demo", action="store_true")
    result.add_argument("--output", type=Path, required=True)
    return result


def main() -> None:
    args = parser().parse_args()
    manifest = build_manifest(args)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    if manifest["status"] == "BLOCKED":
        raise SystemExit(2)


if __name__ == "__main__":
    main()
