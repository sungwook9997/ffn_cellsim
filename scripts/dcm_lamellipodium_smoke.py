"""Smoke test for the per-rim-cell mechanistic lamellipodium (CPU).

Validates :mod:`ffn_sim.cell.dcm_lamellipodium` — the REAL bead-resolved
lamellipodium that replaces the coarse ``ActiveRimTraction`` body force.

Stages (per the H.X smoke convention — finite gate FIRST):
  1. Build a small spheroid (~8-14 cells) + per-rim-cell lamellipodium.
  2. Finite gate: run(0) then run(2000) at small dt; assert all finite.
  3. Spread run: run N batches; measure
       - A/A₀ (basal footprint convex-hull area)
       - rim-cell outward centroid displacement
       - actin beads advanced (real protrusion) + clutches engaged
       - finite/stable
  4. Passive control (no lamellipodium updaters) for comparison.

Run from repo root:
    PYTHONPATH=. python scripts/dcm_lamellipodium_smoke.py
"""

from __future__ import annotations

import numpy as np

import hoomd

from ffn_sim.cell.dcm_lamellipodium import (
    ResolvedLamellipodiumSpheroid,
    build_lamellipodium_spheroid,
    footprint_area,
    rim_centroid_radius,
)
from ffn_sim.cell.dcm_native_shell import ResolvedNativeDCM


def _all_finite(sim) -> bool:
    snap = sim.state.get_snapshot()
    pos = np.asarray(snap.particles.position)
    return bool(np.all(np.isfinite(pos)))


def main() -> None:
    n_cells = 9

    base = ResolvedNativeDCM(
        subdivisions=1,        # coarse shell (42 nodes/cell) — keep CPU smoke light
        dt=3.0e-10,            # small dt for the finite gate (drop to 1e-10 if needed)
        spacing_factor=2.1,    # compact cluster, no t=0 overlap
    )
    p = ResolvedLamellipodiumSpheroid(
        base=base,
        n_seed_per_cell=6,
        pool_per_cell=200,
        batch_steps=100,
    )

    print(f"=== build (n_cells={n_cells}) ===")
    h = build_lamellipodium_spheroid(p, n_cells, device=hoomd.device.CPU(notice_level=0))
    sim = h["sim"]
    rim = h["rim_cells"]
    print(f"rim cells: {rim.tolist()}  (of {n_cells})")
    print(f"actin seed beads: {h['n_seed']}   ligand beads: {h['n_lig']}   "
          f"mem nodes: {h['n_mem']}")

    cell_of_node = h["cell_of_tag"]
    ranges = h["ranges"]
    mem_tid = h["mem_typeid"]
    z_basal = h["z_basal"]
    band = h["mean_edge"] * 1.5

    A0 = footprint_area(sim, cell_of_node=cell_of_node, z_basal=z_basal,
                        band=band, mem_typeid=mem_tid)
    r0 = rim_centroid_radius(sim, rim_cells=rim, ranges=ranges, mem_typeid=mem_tid)
    print(f"A0 = {A0:.3e} m^2   rim centroid radius r0 = {r0:.3e} m")

    # ---- finite gate FIRST ----
    print("\n=== finite gate ===")
    sim.run(0)
    assert _all_finite(sim), "non-finite after run(0)"
    sim.run(2000)
    ok = _all_finite(sim)
    print(f"run(2000) finite: {ok}")
    assert ok, "non-finite after run(2000) — BAOAB guard would have raised"

    # ---- spread run ----
    print("\n=== spread run ===")
    n_batches = 120
    steps_per = p.batch_steps * 5
    for i in range(n_batches):
        sim.run(steps_per)
        if not _all_finite(sim):
            print(f"  blew up at batch {i}")
            break
    A = footprint_area(sim, cell_of_node=cell_of_node, z_basal=z_basal,
                       band=band, mem_typeid=mem_tid)
    r = rim_centroid_radius(sim, rim_cells=rim, ranges=ranges, mem_typeid=mem_tid)
    n_adv = h["nucleator"].n_promoted
    n_grip = h["clutch"].n_grips
    n_teth = h["tether"].n_tethered
    snap = sim.state.get_snapshot()
    n_actin = int(np.count_nonzero(
        np.asarray(snap.particles.typeid) == h["actin_typeid"]))

    print(f"actin beads advanced (protrusion events): {n_adv}")
    print(f"clutches engaged (grip bonds): {n_grip}")
    print(f"membrane nodes tethered (last eval): {n_teth}")
    print(f"total actin_lamel beads now: {n_actin} (seed was {h['n_seed']})")
    print(f"A/A0 = {A / A0:.4f}   (A={A:.3e}, A0={A0:.3e})")
    print(f"rim centroid radius: r0={r0:.3e} -> r={r:.3e}  "
          f"(Δ={r - r0:+.3e} m, {1e6 * (r - r0):+.3f} µm outward)")
    print(f"finite/stable: {_all_finite(sim)}")

    # ---- passive control (no lamellipodium) ----
    print("\n=== passive control (no lamellipodium updaters) ===")
    hp = build_lamellipodium_spheroid(p, n_cells,
                                      device=hoomd.device.CPU(notice_level=0))
    simp = hp["sim"]
    # strip the nucleation + clutch updaters (keep BAOAB), so no protrusion drives.
    keep = []
    for u in list(simp.operations.updaters):
        act = getattr(u, "_action", None) or getattr(u, "action", None)
        cn = type(act).__name__ if act is not None else ""
        if cn in ("LeadingEdgeNucleationUpdater", "FrontClutchRatchetUpdater"):
            continue
        keep.append(u)
    simp.operations.updaters.clear()
    for u in keep:
        simp.operations.updaters.append(u)
    A0p = footprint_area(simp, cell_of_node=hp["cell_of_tag"], z_basal=hp["z_basal"],
                         band=band, mem_typeid=hp["mem_typeid"])
    r0p = rim_centroid_radius(simp, rim_cells=hp["rim_cells"], ranges=hp["ranges"],
                              mem_typeid=hp["mem_typeid"])
    simp.run(0)
    for i in range(n_batches):
        simp.run(steps_per)
        if not _all_finite(simp):
            break
    Ap = footprint_area(simp, cell_of_node=hp["cell_of_tag"], z_basal=hp["z_basal"],
                        band=band, mem_typeid=hp["mem_typeid"])
    rp = rim_centroid_radius(simp, rim_cells=hp["rim_cells"], ranges=hp["ranges"],
                             mem_typeid=hp["mem_typeid"])
    print(f"passive A/A0 = {Ap / A0p:.4f}   rim Δr = {1e6 * (rp - r0p):+.3f} µm")

    dr_active = 1e6 * (r - r0)
    dr_passive = 1e6 * (rp - r0p)
    print("\n=== SUMMARY ===")
    print(f"lamellipodium: A/A0={A / A0:.4f}  rimΔr={dr_active:+.3f}µm  "
          f"adv={n_adv} grip={n_grip}")
    print(f"passive:       A/A0={Ap / A0p:.4f}  rimΔr={dr_passive:+.3f}µm")
    print(f"\nMIGRATION DIFFERENTIAL (active - passive rim outward Δr): "
          f"{dr_active - dr_passive:+.3f} µm  "
          f"({100 * (dr_active - dr_passive) / max(dr_passive, 1e-9):+.1f}% over passive)")
    print("Mechanism check — real protrusion (adv), real FA grip (grip), real "
          "traction (rim cells migrate outward MORE than passive turgor wetting):",
          "PASS" if (n_adv > 0 and n_grip > 0 and dr_active > dr_passive) else "WEAK")


if __name__ == "__main__":
    main()
