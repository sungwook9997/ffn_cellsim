#!/usr/bin/env python
r"""Controls for the interface-residual gate — does it fail when the interface is broken?

WHY A SCRIPT AND NOT THE TEST.  `tests/ac/engine/test_coupled_solve.py` carries these same checks, but
the GPU host's environment has no ``pytest``, so on the only machine where the gate can actually run
the test cannot. This driver runs the identical assertions with no framework, and writes a record.

WHAT IT PROVES.  `measure_interface_residual` scores a two-body cut by reducing each sub-body's force
array into one of the ledger's two independently sourced channels and requiring cancellation within
the D8-derived tolerance. A gate that always passes and a gate that cannot fail are the same artifact,
so the useful checks are the ones that must FAIL:

  PAIR        equal and opposite -> BALANCED         (the gate can pass)
  ONE-SIDED   one node's reaction never scattered    -> IMBALANCED, residual == the injected force
  SIGN-FLIP   scattered with the wrong sign          -> IMBALANCED
  UNBOUND     no runtime bound at all                -> every cut UNBOUND, none reported balanced

The ONE-SIDED case is quantitative on purpose: the residual must equal the force that was withheld, to
round-off. A gate that merely says "not balanced" could be failing for any reason; one that returns
the injected magnitude is measuring what it claims to measure.

Sanity Gate: dimensional — forces vec3d [pN], residual [pN]. Boundary — UNBOUND is asserted to be
distinct from BALANCED. Conservation — the two channels come from two never-merged arrays. Sign sense
— the sign-flip case is the second failure mode and is checked separately from the one-sided case.
CFL/precision — no integration; one force state is read. Measurement protocol — host readback only
after each cut's device reduction completes.

Runtime: CUDA only. Raises on a non-CUDA device.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import warp as wp

from aleph.engine.contracts import reference_cell_architecture
from aleph.engine.coupled_solve import measure_interface_residual

#: The withheld force in the one-sided case [pN]. Not physical and not a tolerance — an injected
#: signal whose EXACT recovery is the assertion, the same way the sf_motor lane injects a stall force.
INJECTED_PN = 0.5


def main() -> None:
    ap = argparse.ArgumentParser(description="Interface-residual gate controls (CUDA).")
    ap.add_argument("--out", default="")
    args = ap.parse_args()

    wp.init()
    device = wp.get_device()
    if not device.is_cuda:
        raise SystemExit("I0-A: the balance gate requires CUDA; there is no CPU simulation path.")
    dev = str(device)

    arch = reference_cell_architecture()
    cut = next(iter(arch.attributable_cuts()))

    class _Body:
        def __init__(self, arr: np.ndarray) -> None:
            self.force_d = wp.array(np.ascontiguousarray(arr), dtype=wp.vec3d, device=dev)

    def score(force_a: np.ndarray, force_b: np.ndarray):
        bodies = {cut[0]: _Body(force_a), cut[1]: _Body(force_b)}
        report = measure_interface_residual(
            world=None, architecture=arch, device=dev,
            component_runtime=lambda name: bodies.get(name),
        )
        return next(c for c in report.cuts if c.components == cut), report

    f = np.array([[1.0, -2.0, 0.5], [0.25, 0.0, -3.0]])
    results: dict[str, dict] = {}

    pair, pair_report = score(f, -f)
    results["pair"] = {"status": pair.status, "residual_pn": pair.residual_pn,
                       "n_imbalanced": len(pair_report.imbalanced)}

    withheld = -f.copy()
    withheld[0, 0] += INJECTED_PN
    one_sided, _ = score(f, withheld)
    results["one_sided"] = {"status": one_sided.status, "residual_pn": one_sided.residual_pn,
                            "injected_pn": INJECTED_PN}

    flipped, _ = score(f, f)
    results["sign_flip"] = {"status": flipped.status, "residual_pn": flipped.residual_pn}

    unbound = measure_interface_residual(
        world=None, architecture=arch, device=dev, component_runtime=lambda _n: None,
    )
    results["unbound"] = {"n_cuts": len(unbound.cuts), "n_measured": len(unbound.measured),
                          "n_imbalanced": len(unbound.imbalanced),
                          "all_unbound": all(c.status == "UNBOUND" for c in unbound.cuts)}

    print(f"[control] device={dev}  attributable cut = {cut}  connector = {arch.attributable_cuts()[cut]}")
    for name, row in results.items():
        print(f"[control] {name:11s} {row}")

    failures: list[str] = []
    if results["pair"]["status"] != "BALANCED":
        failures.append("a genuine Newton pair did not balance — the gate cannot pass")
    if results["one_sided"]["status"] != "IMBALANCED":
        failures.append("a one-sided scatter was NOT caught — the gate cannot fail")
    elif abs((results["one_sided"]["residual_pn"] or 0.0) - INJECTED_PN) > 1e-9 * INJECTED_PN:
        failures.append(
            f"the one-sided residual {results['one_sided']['residual_pn']!r} pN does not recover the "
            f"injected {INJECTED_PN} pN — the gate fails for some other reason than the injection"
        )
    if results["sign_flip"]["status"] != "IMBALANCED":
        failures.append("a sign-flipped scatter was NOT caught")
    if not results["unbound"]["all_unbound"] or results["unbound"]["n_measured"]:
        failures.append("an unbound world did not report every cut UNBOUND")

    record = {
        "schema": "ac.engine.observe/run-record@2",
        "run_label": "interface_residual_gate_controls",
        "kind": "diagnostic",
        "device": dev,
        "evidence": "CUDA_UNIT",
        "quantitative_claim_status": "BLOCKED",
        "evidence_basis": (
            "controls on the two-body force-balance gate: a genuine pair balances, a one-sided "
            "scatter is caught and returns the injected magnitude, a sign flip is caught, and an "
            "unbound world reports UNBOUND rather than balanced. BLOCKED because no physics ran — "
            "these are injected forces on a synthetic pair, not a cell."
        ),
        "cut": {"components": list(cut), "connector": arch.attributable_cuts()[cut]},
        "results": results,
        "verdict": "PASS" if not failures else "FAIL",
        "failures": failures,
    }
    if args.out:
        out = Path(args.out)
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(json.dumps(record, indent=2), encoding="utf-8")
        print(f"[control] wrote {out}")

    if failures:
        for line in failures:
            print(f"[control] FAIL — {line}")
        raise SystemExit("[control] the interface-residual gate did not pass its own controls.")
    print("[control] PASS — the gate passes on a true pair and FAILS on both broken interfaces.")


if __name__ == "__main__":
    main()
