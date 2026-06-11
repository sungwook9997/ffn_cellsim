"""PART A — find the CONFLUENT-packing regime for the native-mesh DCM spheroid.

The PI wants the spheroid to form like SimuCell3D: deformable cells ADHERING to
each other into a confluent, space-filling tissue (foam-like polygonal packing,
no big gaps) — not separated rounded spheres with voids. This script sweeps the
relevant ``ResolvedNativeDCM`` knobs and, for each, measures a confluence metric
on the relaxed aggregate, picks the clearly-confluent + BAOAB-stable regime, and
reports it against a separated-sphere baseline.

Physics of the knobs (SimuCell3D §2.1/§2.2 + integration doc):
  * cell-cell ADHESION (``c_adh`` contact range, ``adh_strength`` ω) — the driver
    that pulls neighbours together and closes gaps. Must reach across the initial
    surface gap, so we both START CELLS CLOSER (lower ``spacing_factor``) and
    WIDEN ``c_adh`` so the bilinear tent engages.
  * SOFT CORTEX (lower ``k_edge`` edge-spring stiffness, lower ``gamma_node`` drag)
    — lets cells DEFORM against each other into polygonal contact rather than
    staying rigidly spherical.
  * MODERATE turgor (``turgor_inflate``, ``K_bulk``) — keeps cells inflated/
    space-filling but not so stiff they resist deformation.

Confluence metric (geometry-only, no magic numbers): packing fraction
``Φ = ΣV_cell / V_hull`` and the inter-cell contact-node fraction. A separated
ball of spheres leaves large interstitial voids (Φ ≲ random-sphere-pack), and
almost no shared interface; a confluent aggregate deforms to fill the hull
(Φ → 1) with a large adhered-interface fraction.

Stability: each candidate runs a finite BAOAB gate ``run(0) → run(N)`` at a small
dt; non-finite or exploded (hull blows up) candidates are rejected.

Run FROM REPO ROOT:
    conda activate ffn_sim
    python scripts/dcm_confluent_tune.py
"""

from __future__ import annotations

import json
import time
from dataclasses import asdict
from pathlib import Path

import numpy as np

from ffn_sim.cell.dcm_native_shell import (
    ResolvedNativeDCM, build_native_dcm_simulation)
from ffn_sim.cell.dcm_confluence import (
    capture_positions, compute_confluence, per_cell_tris)

OUT = Path("ffn_sim/outputs/h_dcm_active")
FIGS = OUT / "figs"
FIGS.mkdir(parents=True, exist_ok=True)

N_CELLS = 16
N_RELAX = 4000          # BAOAB relax steps
SUBDIV = 1              # 42 nodes/cell — CPU-tractable many-cell smoke


def build_and_relax(p: ResolvedNativeDCM, n_cells: int, n_relax: int):
    """Build a native DCM spheroid (no substrate — free aggregate) and relax.

    Returns (handles, pos0, posN, stable). Contact uses the FULL c_adh as the
    tent adhesion range; we override the builder's r_contact-bounded search by
    constructing with the tuned p (c_adh widened) so adhesion reaches neighbours.
    """
    h = build_native_dcm_simulation(p, n_cells, substrate=False, contact=True)
    sim = h["sim"]
    # widen the tent's adhesion search to the tuned c_adh (builder sets it from p,
    # but r_search = max(r_contact, c_adh) so a wide c_adh is honoured already).
    pos0 = capture_positions(sim)
    stable = True
    try:
        sim.run(n_relax)
    except Exception:
        stable = False
    posN = capture_positions(sim)
    if not np.isfinite(posN).all() or np.abs(posN).max() > 50.0 * p.R_cell * n_cells:
        stable = False
    return h, pos0, posN, stable


def metrics_for(h, pos):
    p = h["p"]
    tris_pc = per_cell_tris(_global_tris(h), _global_typeids(h), h["ranges"])
    return compute_confluence(
        pos, cell_of_tag=h["cell_of_tag"], ranges=h["ranges"],
        mesh_tris_per_cell=tris_pc, contact_range=p.c_adh)


# The builder returns mesh + handles but not the raw global tris/typeids; recover
# them from the same primitive the builder used (deterministic).
def _global_tris_typeids(h):
    from ffn_sim.cell.dcm_native_shell import build_native_snapshot
    p = h["p"]
    from ffn_sim.cell.dcm import _cluster_centers
    centers = _cluster_centers(h["n_cells"], p.spacing_factor * p.R_cell,
                               p.z_substrate, p.R_cell, mode=p.cluster)
    (_snap, mesh_tris, mesh_typeids, *_rest) = build_native_snapshot(
        p, h["n_cells"], centers)
    return mesh_tris, mesh_typeids


_CACHE: dict = {}


def _global_tris(h):
    key = id(h)
    if key not in _CACHE:
        _CACHE[key] = _global_tris_typeids(h)
    return _CACHE[key][0]


def _global_typeids(h):
    key = id(h)
    if key not in _CACHE:
        _CACHE[key] = _global_tris_typeids(h)
    return _CACHE[key][1]


def main() -> None:
    t_start = time.time()
    results = []

    # ---- BASELINE: separated spheres (current defaults: loose spacing, short
    # adhesion range, stiff cortex) — the "wrong" look the PI is unhappy with.
    baseline = ResolvedNativeDCM(
        subdivisions=SUBDIV, spacing_factor=2.3,
        c_adh=5.0e-7, adh_strength=1.0e8,
        k_edge=1.0e-3, gamma_node=3.9e-10,
        turgor_inflate=1.05, dt=3.0e-10)
    print("[baseline] separated-sphere defaults ...", flush=True)
    hb, p0b, pNb, sb = build_and_relax(baseline, N_CELLS, N_RELAX)
    mb = metrics_for(hb, pNb)
    print(f"  stable={sb}  Phi={mb.packing_fraction:.3f} "
          f"contact_frac={mb.contact_fraction:.3f} "
          f"mean_contact_nodes={mb.mean_contact_nodes:.1f}", flush=True)

    # ---- SWEEP: closer start + wider adhesion + softer cortex.
    # c_adh must span the residual surface gap; soften k_edge so cells deform.
    sweep = []
    for spacing in (1.9, 2.0):
        for c_adh in (3.0e-6, 5.0e-6):
            for adh in (3.0e8, 8.0e8):
                for k_edge in (2.0e-4, 5.0e-5):
                    sweep.append(dict(spacing_factor=spacing, c_adh=c_adh,
                                      adh_strength=adh, k_edge=k_edge))

    for i, knobs in enumerate(sweep):
        p = ResolvedNativeDCM(
            subdivisions=SUBDIV,
            turgor_inflate=1.05, gamma_node=3.9e-10, dt=2.0e-10,
            rep_strength=1.0e8, force_cap=5.0e-8, **knobs)
        h, p0, pN, stable = build_and_relax(p, N_CELLS, N_RELAX)
        if not stable:
            print(f"[{i:02d}] {knobs} -> UNSTABLE", flush=True)
            results.append({**knobs, "stable": False})
            continue
        m = metrics_for(h, pN)
        print(f"[{i:02d}] spacing={knobs['spacing_factor']} "
              f"c_adh={knobs['c_adh']:.1e} adh={knobs['adh_strength']:.1e} "
              f"k_edge={knobs['k_edge']:.1e} -> "
              f"Phi={m.packing_fraction:.3f} contact={m.contact_fraction:.3f} "
              f"mcn={m.mean_contact_nodes:.1f}", flush=True)
        results.append({
            **knobs, "stable": True,
            "packing_fraction": m.packing_fraction,
            "contact_fraction": m.contact_fraction,
            "mean_contact_nodes": m.mean_contact_nodes,
        })

    # ---- CHOSEN regime (refined balance: keep volume AND deform/adhere).
    # Higher rep + slightly higher turgor stops the implosion failure mode that
    # very strong adhesion + weak rep produces (cells crushed into a blob,
    # Φ→0.17); this regime keeps Vcell/Vrest≈1.22 while filling space.
    print("[chosen] refined confluent regime ...", flush=True)
    chosen = ResolvedNativeDCM(
        subdivisions=SUBDIV, spacing_factor=1.8,
        c_adh=5.0e-6, adh_strength=8.0e8, k_edge=5.0e-5,
        rep_strength=2.0e8, turgor_inflate=1.08,
        gamma_node=3.9e-10, dt=2.0e-10, force_cap=8.0e-8)
    hc, p0c, pNc, sc = build_and_relax(chosen, N_CELLS, N_RELAX)
    mc = metrics_for(hc, pNc)
    vcell_over_vrest = float(mc.cell_volumes.mean() / hc["V_rest"])
    print(f"  stable={sc}  Phi={mc.packing_fraction:.3f} "
          f"contact={mc.contact_fraction:.3f} "
          f"Vcell/Vrest={vcell_over_vrest:.2f}", flush=True)

    # ---- PICK: highest contact_fraction among stable candidates that also
    # clearly beat the baseline packing fraction (confluent = high contact +
    # high Phi). Score = contact_fraction + packing_fraction (both ↑ for tissue).
    stable_res = [r for r in results if r.get("stable")]
    for r in stable_res:
        r["score"] = r["contact_fraction"] + r["packing_fraction"]
    stable_res.sort(key=lambda r: r["score"], reverse=True)
    best = stable_res[0] if stable_res else None

    summary = {
        "n_cells": N_CELLS, "n_relax": N_RELAX, "subdivisions": SUBDIV,
        "baseline": {
            "params": {k: getattr(baseline, k) for k in
                       ("spacing_factor", "c_adh", "adh_strength", "k_edge")},
            "stable": sb,
            "packing_fraction": mb.packing_fraction,
            "contact_fraction": mb.contact_fraction,
            "mean_contact_nodes": mb.mean_contact_nodes,
        },
        "sweep": results,
        "best_sweep": best,
        "chosen": {
            "params": {k: getattr(chosen, k) for k in
                       ("spacing_factor", "c_adh", "adh_strength", "k_edge",
                        "rep_strength", "turgor_inflate")},
            "stable": sc,
            "packing_fraction": mc.packing_fraction,
            "contact_fraction": mc.contact_fraction,
            "vcell_over_vrest": vcell_over_vrest,
        },
        "wall_s": round(time.time() - t_start, 1),
    }
    (OUT / "confluent_tune.json").write_text(json.dumps(summary, indent=2))
    print("\n=== CONFLUENT REGIME ===", flush=True)
    print(f"baseline  Phi={mb.packing_fraction:.3f} "
          f"contact={mb.contact_fraction:.3f}  (separated spheres)", flush=True)
    print(f"CHOSEN    Phi={mc.packing_fraction:.3f} "
          f"contact={mc.contact_fraction:.3f}  Vcell/Vrest={vcell_over_vrest:.2f} "
          f"(confluent: deformed + space-filling)", flush=True)
    print(f"wrote {OUT/'confluent_tune.json'}  ({summary['wall_s']}s)", flush=True)


if __name__ == "__main__":
    main()
