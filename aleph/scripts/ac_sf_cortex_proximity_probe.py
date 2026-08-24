#!/usr/bin/env python
r"""Can ``sf_cortex_transient`` be wired at all? Measure the SF↔cortex separation first.

`sf_cortex_transient` is the one unwired power port whose every constant is already SOURCED:
α-actinin stiffness 4.6e5 pN/µm (Ferrer 2008, PI-approved 2026-06-30), capture radius 0.06 µm,
``k_on`` 10/s, ``k_off0`` 0.066/s. So its pairing topology is DERIVED — SF nodes within the sourced
capture radius of a cortex node — rather than chosen, which is what makes it the next wiring target.

**But derived-from-geometry is only a topology if the geometry supplies pairs.** The SF arc population
and the cortex are built by different builders in different frames; nothing has ever checked that they
are within 60 nm of each other. If they are not, generating specs would produce zero joints and the
"wiring" would be vacuous — the same trap the interface residual just hit at rest, where a BALANCED
verdict meant only that there was nothing to balance.

So this measures the separation before anything is built on the assumption. It reports the minimum
SF→cortex distance, how many SF nodes have a cortex node inside the sourced capture radius, and the
distance distribution — enough to say whether the edge is wireable, needs a different pairing rule, or
is geometrically absent in this configuration.

NOTHING IS CHOSEN HERE. The only length compared against is `ALPHA_ACTININ.capture_radius_um`, read
from the sourced hand-parameter record. Host-side geometry at BUILD time, which is where topology
generation belongs (`build_sf_mechanics_topology` and the Mikado builder both work this way); no
per-step readback and no physical-time loop exists in this script.

Sanity Gate: dimensional — every distance µm. Boundary — "no pair within the radius" is reported as
that, never as a zero-joint success. Conservation — read-only, allocates no force. Sign sense — n/a.
CFL/precision — no integration. Protocol — one build, one host distance query, no device loop.

Runtime: CUDA for the cell build; the distance query is host NumPy on build-time positions.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import warp as wp


def main() -> None:
    ap = argparse.ArgumentParser(description="SF↔cortex separation vs the sourced capture radius.")
    ap.add_argument("--k-axial", type=float, required=True, help="SF axial stiffness [pN/um] (PI-GAP)")
    ap.add_argument("--cortex-filaments", type=int, default=70686)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--commit", type=str, default="")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    wp.init()
    if not wp.get_device().is_cuda:
        raise SystemExit("I0-A: requires CUDA for the cell build.")
    dev = str(wp.get_device())

    from aleph.components.incumbent.assemble import (
        NMII_BACKBONE_LP_DIAGNOSTIC_UM, CellConfig, build_cell,
    )
    from aleph.engine.sf_mechanics import ALPHA_ACTININ
    from aleph.engine.sf_population import build_sf_arc_population

    capture_um = float(ALPHA_ACTININ.capture_radius_um)

    cell = build_cell(CellConfig(
        n_filaments=int(args.cortex_filaments), with_myosin=True, overlap_free_cortex=True,
        membrane_subdivisions=6, nucleus_subdivisions=3, erm_radial_pairing=True,
        resting_bound_myosin_fraction=None, nmii_straddle_placement=True,
        nmii_backbone_lp_um=NMII_BACKBONE_LP_DIAGNOSTIC_UM,
        nmii_backbone_lp_source="PROXIMITY_PROBE_FIXTURE_NOT_PRODUCTION_PI_GAP",
        seed=int(args.seed),
    ))
    sf_pop = build_sf_arc_population(n_ventral=8, n_dorsal=4, n_arc=4, n_cap=4, n_per_fiber=9)
    sf_pop.assert_partitioned()

    sf = np.ascontiguousarray(sf_pop.pos, dtype=np.float64).reshape(-1, 3)
    cortex = np.ascontiguousarray(cell.pos_d.numpy(), dtype=np.float64).reshape(-1, 3)[: cell.n_actin]

    try:
        from scipy.spatial import cKDTree
        tree = cKDTree(cortex)
        nearest, _ = tree.query(sf, k=1)
        within = tree.query_ball_point(sf, r=capture_um)
        n_pairs = int(sum(len(hits) for hits in within))
    except ImportError:  # brute force is affordable at this SF size and keeps the probe runnable
        d = np.linalg.norm(sf[:, None, :] - cortex[None, :, :], axis=-1)
        nearest = d.min(axis=1)
        n_pairs = int((d <= capture_um).sum())

    sf_with_partner = int((nearest <= capture_um).sum())
    record = {
        "schema": "ac.engine.observe/run-record@2",
        "run_label": "sf_cortex_proximity_probe",
        "kind": "diagnostic",
        "build": {"commit": args.commit or None, "source": "declared"},
        "device": dev,
        "evidence": "NATIVE",
        "quantitative_claim_status": "BLOCKED",
        "evidence_basis": (
            "build-time separation between the sf_arc population and the cortex actin nodes, compared "
            "against the SOURCED alpha-actinin capture radius. BLOCKED: a geometric readiness check "
            "for wiring, not a physics result — nothing was solved and no force was computed."
        ),
        "capture_radius_um": capture_um,
        "capture_radius_provenance": "SOURCED — ALPHA_ACTININ.capture_radius_um (Ferrer 2008 line)",
        "census": {"cortex_filaments": int(args.cortex_filaments), "n_cortex_actin_nodes": int(cortex.shape[0]),
                   "n_sf_nodes": int(sf.shape[0])},
        "separation_um": {
            "min": float(nearest.min()), "median": float(np.median(nearest)),
            "max": float(nearest.max()), "mean": float(nearest.mean()),
        },
        "pairs_within_capture": {"n_sf_nodes_with_a_partner": sf_with_partner,
                                 "n_sf_nodes": int(sf.shape[0]), "n_candidate_pairs": n_pairs},
        "wireable": bool(sf_with_partner > 0),
    }
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out).write_text(json.dumps(record, indent=2), encoding="utf-8")

    s = record["separation_um"]
    print(f"\n[proximity] capture radius {capture_um:.4f} um (SOURCED)  "
          f"SF nodes {sf.shape[0]}  cortex actin nodes {cortex.shape[0]}")
    print(f"[proximity] SF->cortex nearest distance  min {s['min']:.6f}  median {s['median']:.6f}  "
          f"max {s['max']:.6f} um")
    print(f"[proximity] SF nodes with a partner inside the capture radius: "
          f"{sf_with_partner} / {sf.shape[0]}   candidate pairs: {n_pairs}")
    if sf_with_partner == 0:
        print(f"[proximity] NOT WIREABLE AS DECLARED — the nearest SF node is {s['min']:.4f} um from any "
              f"cortex node, {s['min'] / capture_um:.1f}x the sourced capture radius. Generating specs "
              f"here would emit ZERO joints, and a connector with no joints is not a wired connector. "
              f"The populations are built in frames that do not meet; that is a PLACEMENT question and "
              f"belongs to whoever owns the SF inventory, not to a capture radius that is already sourced.")
    else:
        print("[proximity] WIREABLE — the pairing topology is derivable from the sourced radius.")
    print(f"[proximity] wrote {args.out}")


if __name__ == "__main__":
    main()
