#!/usr/bin/env python3
"""Lab-invariant, within-study mechanistic contrast evaluation.

Absolute calcium traces are acquisition-domain dependent.  This evaluator uses
predeclared within-study control/treatment contrasts and keeps the independent
lab untouched until the final confirmatory calculation.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from evaluate_piezo1_transfer import _curves


def _hedges_g(treated: np.ndarray, control: np.ndarray) -> float:
    nt, nc = len(treated), len(control)
    pooled = np.sqrt(
        ((nt - 1) * treated.var(ddof=1) + (nc - 1) * control.var(ddof=1))
        / (nt + nc - 2)
    )
    if pooled == 0:
        return 0.0
    correction = 1 - 3 / (4 * (nt + nc) - 9)
    return float(correction * (treated.mean() - control.mean()) / pooled)


def _contrast(
    treated: np.ndarray,
    control: np.ndarray,
    rng: np.random.Generator,
    draws: int,
) -> tuple[dict, np.ndarray]:
    effects = np.empty(draws, dtype=float)
    for index in range(draws):
        effects[index] = (
            treated[rng.integers(len(treated), size=len(treated))].mean()
            - control[rng.integers(len(control), size=len(control))].mean()
        )
    result = {
        "treated_n": len(treated),
        "control_n": len(control),
        "mean_difference": float(treated.mean() - control.mean()),
        "hedges_g": _hedges_g(treated, control),
        "bootstrap_ci95": [float(value) for value in np.quantile(effects, [0.025, 0.975])],
        "bootstrap_probability_negative": float(np.mean(effects < 0)),
    }
    return result, effects


def evaluate(
    discovery_manifest: Path,
    external_manifest: Path,
    draws: int = 20_000,
) -> dict:
    discovery_x, discovery_y, _, _ = _curves(discovery_manifest, "hela")
    external_x, external_y, external_conditions, _ = _curves(
        external_manifest, "external"
    )
    # curve_features index 40 is the normalized trapezoidal AUC.
    discovery_auc = discovery_x[:, 40]
    external_auc = external_x[:, 40]
    rng = np.random.default_rng(20260805)
    discovery, _ = _contrast(
        discovery_auc[discovery_y == 1], discovery_auc[discovery_y == 0], rng, draws
    )

    # Predeclared from the author README: LooPINS and CaPINS are the two
    # PIEZO1-targeting nanoswitches. BINPs and no-NP arms are controls, not
    # confirmatory active-actuator strata.
    confirmatory_prefixes = ("LooPINS", "CaPINS")
    confirmatory: dict[str, dict] = {}
    bootstrap_effects = []
    for prefix in confirmatory_prefixes:
        selected = np.asarray([
            condition.startswith(prefix) for condition in external_conditions
        ])
        contrast, effects = _contrast(
            external_auc[selected & (external_y == 1)],
            external_auc[selected & (external_y == 0)],
            rng,
            draws,
        )
        confirmatory[prefix] = contrast
        bootstrap_effects.append(effects)
    equal_stratum_effects = np.stack(bootstrap_effects).mean(axis=0)
    pooled_ci = [float(value) for value in np.quantile(equal_stratum_effects, [0.025, 0.975])]
    negative_groups = sum(value["mean_difference"] < 0 for value in confirmatory.values())
    passed = negative_groups == len(confirmatory) and pooled_ci[1] < 0
    return {
        "task": "within_study_Piezo1_inhibition_calcium_AUC_contrast",
        "representation": "within-study treated-minus-matched-control effect",
        "discovery": {
            "dataset": "zenodo-18495219 HeLa",
            "contrast": "GsMTx4 minus untreated control",
            **discovery,
        },
        "external_confirmatory": {
            "dataset": "zenodo-19890985 EA.hy926",
            "lab_holdout": True,
            "predeclared_active_strata": list(confirmatory_prefixes),
            "strata": confirmatory,
            "direction_consistency": f"{negative_groups}/{len(confirmatory)}",
            "equal_stratum_mean_difference": float(equal_stratum_effects.mean()),
            "equal_stratum_bootstrap_ci95": pooled_ci,
            "bootstrap_probability_negative": float(np.mean(equal_stratum_effects < 0)),
        },
        "mechanism_holdout_passed": passed,
        "cell_state_classifier_holdout_passed": False,
        "observable_constraint_eligible": passed,
        "aleph_parameter_proposal_eligible": False,
        "aleph_authority": "none",
        "limitations": [
            "confirmatory endpoint is a condition-level contrast, not per-cell classification",
            "biological replicate identity is not available per EA.hy926 cell",
            "the two active nanoswitch strata are heterogeneous and are weighted equally",
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--discovery-manifest", type=Path, required=True)
    parser.add_argument("--external-manifest", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--bootstrap-draws", type=int, default=20_000)
    args = parser.parse_args()
    result = evaluate(
        args.discovery_manifest, args.external_manifest, args.bootstrap_draws
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(result["external_confirmatory"], sort_keys=True))


if __name__ == "__main__":
    main()
