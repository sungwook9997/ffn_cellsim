#!/usr/bin/env python3
"""Exercise the pre-sweep promotion contract on membrane-tether inference."""

from __future__ import annotations

import json
from pathlib import Path

from external_inference import infer
from latent_promotion import evaluate


def main() -> None:
    inference = infer({
        "modality": "AFM_membrane_tether_force",
        "value": 18.697398,
        "unit": "pN",
        "standard_error": 0.5,
        "bending_rigidity_pn_um": 0.082,
        "provenance": {"dataset": "elife-72381", "role": "worked_example"},
    })
    result = evaluate(
        inference,
        external_dataset_holdout_passed=False,
        external_lab_holdout_passed=False,
        uncertainty_calibrated=False,
        operator_mapping=None,
        forward_sensitivity_passed=False,
    )
    output = Path(__file__).resolve().parent / "results" / "latent_promotion_report.json"
    output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()
