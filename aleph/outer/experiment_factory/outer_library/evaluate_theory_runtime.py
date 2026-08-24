#!/usr/bin/env python3
"""Verify equation residuals, domain refusals, and exact ambiguity outputs."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from theory_runtime import (
    evaluate_expert,
    persistent_random_walk_msd_2d,
    standard_linear_solid_relaxation,
    tether_ambiguity_set,
    tether_apparent_tension,
    tether_force,
)


SCHEMA = "outer.biophysics.theory_runtime_evaluation.v1"


def evaluate(registry_path: Path) -> dict[str, object]:
    registry = json.loads(registry_path.read_text(encoding="utf-8"))
    if registry["schema"] != "outer.biophysics.theory_registry.v1":
        raise ValueError("wrong theory registry schema")
    if registry["aleph_dependency"]:
        raise ValueError("independent theory registry acquired an Aleph dependency")
    theory_ids = [item["theory_id"] for item in registry["theories"]]
    if len(theory_ids) != len(set(theory_ids)):
        raise ValueError("duplicate theory ID")
    implemented = [item for item in registry["theories"] if item["status"] == "implemented"]
    if any(not item["implementation"] for item in implemented):
        raise ValueError("implemented theory lacks runtime binding")

    rng = np.random.default_rng(20260806)
    kappa = np.exp(rng.uniform(np.log(0.01), np.log(1.0), 2_000))
    sigma = np.exp(rng.uniform(np.log(1e-4), np.log(100.0), 2_000))
    force = np.asarray([
        tether_force(a, b)["tether_force_pn"] for a, b in zip(kappa, sigma)
    ])
    recovered = np.asarray([
        tether_apparent_tension(f, a)["apparent_tension_pn_per_um"]
        for f, a in zip(force, kappa)
    ])
    tether_relative_residual = np.abs(recovered - sigma) / sigma

    ambiguity = tether_ambiguity_set(21.4841, points=257)
    ambiguity_residual = np.max(np.abs(
        ambiguity["sigma_bilayer_pn_per_um"]
        + ambiguity["attachment_energy_density_pn_per_um"]
        - 21.4841
    ))

    time = np.linspace(0.0, 50.0, 501)
    sls = standard_linear_solid_relaxation(time, 7.0, 2.0, 8.5)["response"]
    sls_expected = 2.0 + 5.0 * np.exp(-time / 8.5)
    sls_residual = float(np.max(np.abs(sls - sls_expected)))

    lag = np.geomspace(1e-6, 1e5, 2_000)
    speed, tau = 2.3, 17.0
    msd = persistent_random_walk_msd_2d(lag, speed, tau)[
        "mean_squared_displacement_um2"
    ]
    short_ballistic_ratio = float(msd[0] / (speed * speed * lag[0] ** 2))
    long_diffusive_ratio = float(msd[-1] / (2.0 * speed * speed * tau * lag[-1]))

    refusal_cases = {
        "negative_tension": evaluate_expert(
            "membrane_tether_quasistatic",
            {"bending_rigidity_pn_um": 0.08, "apparent_tension_pn_per_um": -1.0},
        ),
        "negative_time": evaluate_expert(
            "standard_linear_solid_single_mode",
            {
                "time_s": np.asarray([-1.0, 0.0]),
                "initial_response": 1.0,
                "long_time_response": 0.0,
                "relaxation_time_s": 1.0,
            },
        ),
        "unimplemented_expert": evaluate_expert("molecular_clutch_family", {}),
    }
    if any(not case.get("refused") for case in refusal_cases.values()):
        raise ValueError("theory domain refusal failed")

    # A real measured force can constrain sigma_app only conditional on kappa;
    # this example intentionally returns an ambiguity set, not an Aleph axis.
    measured_force = 18.6974
    assumed_kappa = 0.082
    apparent = tether_apparent_tension(measured_force, assumed_kappa)
    measured_ambiguity = tether_ambiguity_set(
        apparent["apparent_tension_pn_per_um"], points=101
    )

    checks = {
        "registry_has_no_Aleph_dependency": registry["aleph_dependency"] is False,
        "tether_roundtrip_relative_residual_below_1e-12": bool(tether_relative_residual.max() < 1e-12),
        "tether_decomposition_ambiguity_is_exact": bool(ambiguity_residual < 1e-12),
        "single_mode_relaxation_equation_residual_below_1e-12": bool(sls_residual < 1e-12),
        "persistent_walk_short_time_is_ballistic": bool(abs(short_ballistic_ratio - 1.0) < 1e-5),
        "persistent_walk_long_time_is_diffusive": bool(abs(long_diffusive_ratio - 1.0) < 1e-3),
        "invalid_domains_return_typed_refusal": True,
    }
    return {
        "schema": SCHEMA,
        "implemented_expert_count": len(implemented),
        "planned_expert_count": len(registry["theories"]) - len(implemented),
        "checks": checks,
        "all_checks_passed": all(checks.values()),
        "equation_residuals": {
            "tether_roundtrip_max_relative": float(tether_relative_residual.max()),
            "tether_ambiguity_sum_max_absolute": float(ambiguity_residual),
            "single_mode_relaxation_max_absolute": sls_residual,
            "persistent_walk_short_ballistic_ratio": short_ballistic_ratio,
            "persistent_walk_long_diffusive_ratio": long_diffusive_ratio,
        },
        "typed_refusal_examples": refusal_cases,
        "measured_tether_conditional_inference_example": {
            "measured_force_pn": measured_force,
            "assumed_not_measured_bending_rigidity_pn_um": assumed_kappa,
            "apparent_tension_pn_per_um": apparent["apparent_tension_pn_per_um"],
            "ambiguity_relation": measured_ambiguity["relation"],
            "ambiguity_point_count": measured_ambiguity["point_count"],
            "may_choose_bilayer_vs_attachment_decomposition": False,
        },
        "synthetic_or_equation_pass_grants_real_evidence_authority": False,
        "aleph_dependency": False,
        "may_emit_Aleph_parameter": False,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--registry", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    report = evaluate(args.registry)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"all_checks_passed": report["all_checks_passed"], "equation_residuals": report["equation_residuals"]}, sort_keys=True))


if __name__ == "__main__":
    main()
