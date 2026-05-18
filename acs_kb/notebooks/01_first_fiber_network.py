"""Phase 1 Unit 1.1 — first 2D fiber network walkthrough (P1-corrected).

Generates three figures + JSON summary in ``acs_kb/outputs/phase1/``:

  - 01_visual_iso.png        — isotropic Mikado network with emergent XLs
  - 02_emergent_z.png        — measured ⟨z⟩ vs theory across 10 seeds
  - 03_order_check.png       — measured S vs target sweep
"""

from __future__ import annotations

import json
import time
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from acs_kb.common.derived_params import load_config
from acs_kb.common.sanity_gate import gate_phase1_ecm
from acs_kb.ecm.cross_links import (
    compute_xl_energy_and_forces,
    generate_cross_links,
    measure_coordination,
    measure_xl_per_fiber,
)
from acs_kb.ecm.fiber_mechanics import compute_energy, compute_forces
from acs_kb.ecm.fiber_network import (
    generate_2d_fiber_network,
    measure_nematic_order,
)
from acs_kb.ecm.visualization import plot_2d_static

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs" / "phase1"
OUT.mkdir(parents=True, exist_ok=True)
CFG_PATH = ROOT / "configs" / "phase1_unit1.yaml"


def main() -> None:
    cfg = load_config(CFG_PATH)
    ecm = cfg["ecm"]; d = ecm["derived"]

    # 1. Generate
    t0 = time.perf_counter()
    net = generate_2d_fiber_network(
        L_box=ecm["L_box"], n_fibers=d["n_fibers"], L_fiber=ecm["L_fiber"],
        beads_per_fiber=ecm["beads_per_fiber"], S_order=0.0, seed=ecm["seed"],
    )
    t_gen = time.perf_counter() - t0

    t0 = time.perf_counter()
    links = generate_cross_links(net, stiffness=ecm["xl_stiffness"])
    t_xl = time.perf_counter() - t0

    n_beads_total = net.bead_positions.shape[0] * net.bead_positions.shape[1]
    z = measure_coordination(links, n_beads_total, ecm["beads_per_fiber"])
    xl_per_f = measure_xl_per_fiber(links, d["n_fibers"])
    measured_ell_c = ecm["L_fiber"] / max(xl_per_f, 1e-9)

    # 2. Mechanics (backbone + XL)
    t0 = time.perf_counter()
    E_total = compute_energy(
        net.bead_positions, net.rest_length,
        ecm["stretching_modulus"], ecm["bending_modulus"], net.box_size,
        cross_links=links,
    )
    F_total = compute_forces(
        net.bead_positions, net.rest_length,
        ecm["stretching_modulus"], ecm["bending_modulus"], net.box_size,
        cross_links=links,
    )
    E_xl, F_xl = compute_xl_energy_and_forces(net.bead_positions, links, net.box_size)
    t_phys = time.perf_counter() - t0

    # 3. Sanity Gate (P1-corrected)
    gate = gate_phase1_ecm(
        biological_mesh_target=ecm["biological_mesh"],
        measured_segment_length=measured_ell_c,
        expected_segment_length=d["mikado_ell_c_predicted"],
        measured_z=z,
        expected_z=d["expected_total_z"],
        ku13_reference_z=ecm["target_coordination_ref"],
        n_fibers=d["n_fibers"], n_beads_total=n_beads_total, n_links=len(links),
        acceptance=ecm["acceptance"], demo_mode=ecm["demo_mode"],
    )
    gate.raise_if_failed()

    # Figure 1: visual
    plot_2d_static(
        net, cross_links=links,
        save_path=OUT / "01_visual_iso.png",
        title=(f"Phase 1 Unit 1.1 — Mikado, N_f={d['n_fibers']}, "
               f"ℓ_c={measured_ell_c*1e6:.2f} μm, ⟨z⟩={z:.3f}"),
    )

    # Figure 2: emergent ⟨z⟩ across seeds vs theory and KU-1.3
    zs, ells = [], []
    for s in range(10):
        n = generate_2d_fiber_network(
            L_box=ecm["L_box"], n_fibers=d["n_fibers"], L_fiber=ecm["L_fiber"],
            beads_per_fiber=ecm["beads_per_fiber"], S_order=0.0, seed=100 + s,
        )
        lk = generate_cross_links(n, stiffness=ecm["xl_stiffness"])
        zs.append(measure_coordination(lk, n_beads_total, ecm["beads_per_fiber"]))
        x_pf = measure_xl_per_fiber(lk, d["n_fibers"])
        ells.append(ecm["L_fiber"] / max(x_pf, 1e-9))
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.errorbar(range(10), zs, fmt="o", color="steelblue", label="emergent ⟨z⟩")
    ax.axhline(d["expected_total_z"], color="darkorange", ls="--",
               label=f"Mikado prediction = {d['expected_total_z']:.3f}")
    ax.axhline(ecm["target_coordination_ref"], color="red", ls=":",
               label=f"KU-1.3 biological ref = {ecm['target_coordination_ref']:.2f}")
    z_lo, z_hi = ecm["acceptance"]["z_range"]
    ax.axhspan(z_lo, z_hi, color="green", alpha=0.1, label="sub-isostatic band")
    ax.set_xlabel("seed index"); ax.set_ylabel("emergent ⟨z⟩")
    ax.set_title("Emergent coordination — Mikado theory vs measurement vs biology")
    ax.legend(fontsize=8, loc="lower right")
    fig.tight_layout(); fig.savefig(OUT / "02_emergent_z.png", dpi=150); plt.close(fig)

    # Figure 3: order parameter sweep
    targets = [0.0, 0.1, 0.2, 0.3, 0.5, 0.7]
    measured_S = []
    for st in targets:
        n = generate_2d_fiber_network(
            L_box=ecm["L_box"], n_fibers=d["n_fibers"], L_fiber=ecm["L_fiber"],
            beads_per_fiber=ecm["beads_per_fiber"], S_order=st, seed=ecm["seed"],
        )
        measured_S.append(measure_nematic_order(n.fiber_orientations))
    sigma = 1.0 / np.sqrt(2.0 * d["n_fibers"])
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot([0, 1], [0, 1], ls="--", color="gray", label="y=x")
    ax.errorbar(targets, measured_S, yerr=sigma, fmt="o", capsize=4,
                color="darkorange", label=f"measured (±1/√(2N)={sigma:.3f})")
    ax.set_xlim(-0.05, 0.8); ax.set_ylim(-0.05, 0.8)
    ax.set_xlabel("target S_order"); ax.set_ylabel("measured S")
    ax.set_title("Nematic order check (KU-1.9)")
    ax.legend(fontsize=8)
    fig.tight_layout(); fig.savefig(OUT / "03_order_check.png", dpi=150); plt.close(fig)

    summary = {
        "n_fibers": d["n_fibers"],
        "n_beads_total": n_beads_total,
        "n_links": len(links),
        "xl_per_fiber_measured": xl_per_f,
        "xl_per_fiber_predicted": d["expected_xl_per_fiber"],
        "ell_c_measured_m": measured_ell_c,
        "ell_c_predicted_m": d["mikado_ell_c_predicted"],
        "biological_mesh_target_m": ecm["biological_mesh"],
        "z_measured": z,
        "z_predicted": d["expected_total_z"],
        "z_KU13_reference": ecm["target_coordination_ref"],
        "z_10seeds_mean": float(np.mean(zs)),
        "z_10seeds_std": float(np.std(zs)),
        "S_measured_at_targets": dict(zip(targets, measured_S)),
        "energy_total_J": E_total,
        "energy_xl_J": E_xl,
        "max_force_N": float(np.abs(F_total).max()),
        "max_force_xl_N": float(np.abs(F_xl).max()),
        "timings_s": {
            "fiber_generation": t_gen,
            "cross_link_generation": t_xl,
            "energy_plus_forces": t_phys,
            "total_phase1_step": t_gen + t_xl + t_phys,
        },
        "sanity_gate_passed": gate.passed,
        "sanity_gate_report": gate.summary(),
    }
    with open(OUT / "01_summary.json", "w") as f:
        json.dump(summary, f, indent=2, default=float)
    print(json.dumps(summary, indent=2, default=float))


if __name__ == "__main__":
    main()
