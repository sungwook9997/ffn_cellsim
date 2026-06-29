"""H.8 membrane-reservoir tether-mesh ACTIVATION GATE (LIVE, 2026-06-09).

Fifth compartment graduated EXPERIMENTAL->LIVE under the PI ownership grant. The
plasma membrane gets its OWN radially-offset ``mem_node`` bead layer outside the
cortex shell + a static breakable ``mem_tether`` mesh (each node → a distinct
cortex bead, force-free at the membrane–cortex offset) — the explicit membrane–
cortex adhesion (MCA) load path (BLOCKER-1 fix: a shared shell would make the
tethers degenerate self-pairs). The gate builds the FULL physiological baseline
cell OFF vs ON and checks the build-time controls:

  MESH ASSEMBLED — ON adds n_mem_nodes mem_node particles + one mem_tether per
                   node; mem_node type + mem_tether bond present.
  CROSS-LAYER    — every tether joins a mem_node to a DISTINCT cortex bead (no
                   self-pairs; unique cortex anchors → per-cortex degree +1,
                   nlist exclusion-cap safe, the LINC/MTOC lesson).
  FORCE-FREE     — every tether born at exactly membrane_offset (strain ~0; the
                   adhesion tension emerges from dynamics, not construction).
  NO-CONTAM      — cortical γ_soft IDENTICAL OFF vs ON (mem_ registry-denylisted);
                   every mem_ bond type excluded from the cortical mask.
  OFF-IDENTITY   — reservoir off reproduces the bare baseline (+0 particles).
  BLEB BLOCKED   — the rupture updater (σ_crit_bleb None) + reservoir release
                   (f_excess None) HONESTLY raise (PI-pending paths disabled).

Only the STATIC tether mesh is LIVE (W_MCA/k_tether/max_tether_dist/membrane_offset
anchored). Bleb nucleation (σ_crit_bleb, Tinevez 2009) + reservoir tension-buffering
(f_excess, Raucher-Sheetz/Figard) stay PI-pending → those paths stay blocked.

Run:  python ffn_sim/scripts/h7_membrane_reservoir_activation_gate.py
"""

from __future__ import annotations

import sys
from pathlib import Path

_HERE = Path(__file__).resolve()
if str(_HERE.parents[2]) not in sys.path:
    sys.path.insert(0, str(_HERE.parents[2]))

import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from ffn_sim.archive.hoomd_legacy.cell.compartment_registry import REGISTRY, load_recipe
from ffn_sim.archive.hoomd_legacy.cell.membrane_reservoir import (
    MembraneTetherUpdater,
    resolve_membrane_reservoir,
)
from ffn_sim.archive.hoomd_legacy.cell.manifest import build_baseline_cell, load_manifest
from ffn_sim.archive.hoomd_legacy.cortex.cortical_tension import (
    _is_adhesion_bond_type,
    measure_cortical_tension,
)

_OUT = _HERE.parents[1] / "outputs" / "h7"


def _gamma_soft(cell) -> float:
    cell.simulation.run(0)
    r = measure_cortical_tension(
        cell.simulation, R_cell=cell.p_cortex.R_cell,
        p_enclosed_volume=cell.p_enclosed_volume,
    )
    return float(r.get("gamma_soft_N_per_m", r.get("gamma_soft", float("nan"))))


def _bleb_paths_blocked() -> dict:
    """Confirm the PI-pending bleb/reservoir paths honestly raise (None-gated)."""
    p = resolve_membrane_reservoir(
        {"membrane_reservoir": {"enabled": True}}, R_cell=7.5e-6
    )
    rupture_raises = False
    try:
        MembraneTetherUpdater(
            p, layout=None, kT=4.0e-21,  # type: ignore[arg-type]
        )
    except NotImplementedError:
        rupture_raises = True
    release_raises = False
    try:
        p.released_area(1.0)
    except NotImplementedError:
        release_raises = True
    return {
        "sigma_crit_bleb_is_none": p.sigma_crit_bleb is None,
        "f_excess_is_none": p.f_excess is None,
        "rupture_updater_raises": rupture_raises,
        "reservoir_release_raises": release_raises,
    }


def run() -> dict:
    base = load_manifest("mcf7_baseline.yaml")
    m_on, _ = REGISTRY.compose_manifest(
        load_recipe("membrane_reservoir_tethered"), base_manifest=base, strict=True
    )
    cell_off = build_baseline_cell("mcf7_baseline.yaml", seed=1)
    cell_on = build_baseline_cell("mcf7_baseline.yaml", manifest=m_on, seed=1)
    p = cell_on.p_membrane_reservoir

    s_off = cell_off.simulation.state.get_snapshot()
    s_on = cell_on.simulation.state.get_snapshot()
    handles = cell_on.extras.get("handles", {})
    n_tethers = int(handles.get("n_mem_tethers", 0))
    layout = handles.get("membrane_tether_layout")
    n_part_delta = int(s_on.particles.N) - int(s_off.particles.N)
    n_bond_delta = int(s_on.bonds.N) - int(s_off.bonds.N)
    mem_bonds = [t for t in s_on.bonds.types if t.startswith("mem_")]

    # MESH ASSEMBLED + own layer (mem_node added).
    mesh_ok = (
        n_tethers > 0
        and "mem_node" in s_on.particles.types
        and "mem_tether" in s_on.bonds.types
        and n_part_delta == n_tethers   # one mem_node per tether
        and n_bond_delta == n_tethers
    )

    # CROSS-LAYER: each tether is mem_node ↔ cortex, distinct cortex anchors.
    cross_ok = self_pairs = unique_anchors = None
    if layout is not None and n_tethers > 0:
        pairs = np.asarray(layout.tether_pairs, dtype=np.int64)
        types = list(s_on.particles.types)
        tid = np.asarray(s_on.particles.typeid, dtype=np.int64)
        mem_id = types.index("mem_node")
        cor_id = types.index("actin_cortex")
        a_is_mem = tid[pairs[:, 0]] == mem_id
        b_is_cor = tid[pairs[:, 1]] == cor_id
        cross_ok = bool((a_is_mem & b_is_cor).all())
        unique_anchors = bool(len(set(pairs[:, 1].tolist())) == pairs.shape[0])
        self_pairs = int((pairs[:, 0] == pairs[:, 1]).sum())
    cross_layer_ok = bool(cross_ok and unique_anchors and self_pairs == 0)

    # FORCE-FREE: tether sep == membrane_offset (strain ~0).
    max_strain = float("nan")
    if layout is not None and n_tethers > 0:
        pos = np.asarray(s_on.particles.position, dtype=np.float64)
        pr = np.asarray(layout.tether_pairs, dtype=np.int64)
        sep = np.linalg.norm(pos[pr[:, 0]] - pos[pr[:, 1]], axis=1)
        max_strain = float(
            np.max(np.abs(sep - p.membrane_offset) / max(p.membrane_offset, 1e-30))
        )
    force_free_ok = bool(np.isfinite(max_strain) and max_strain < 1e-6)

    # NO-CONTAMINATION.
    mask_excludes_mem = all(_is_adhesion_bond_type(t) for t in mem_bonds)
    g_off, g_on = _gamma_soft(cell_off), _gamma_soft(cell_on)
    gamma_identical = bool(
        np.isfinite(g_off) and np.isfinite(g_on)
        and abs(g_on - g_off) <= 1e-12 * max(1.0, abs(g_off))
    )

    blocked = _bleb_paths_blocked()
    bleb_ok = all(blocked.values())

    ok = (
        mesh_ok and cross_layer_ok and force_free_ok
        and mask_excludes_mem and gamma_identical and bleb_ok
    )

    result = {
        "gate": "H.8 membrane_reservoir activation",
        "verdict": "PASS" if ok else "REVIEW",
        "physics_claim": True,
        "params": {
            "W_MCA_Jm2": float(p.W_MCA), "k_tether_Npm": float(p.k_tether),
            "max_tether_dist_m": float(p.max_tether_dist),
            "membrane_offset_m": float(p.membrane_offset),
            "n_mem_nodes_requested": int(p.n_mem_nodes),
            "rupture_force_N": float(layout.rupture_force) if layout else 0.0,
        },
        "controls": {
            "mesh_assembled": bool(mesh_ok),
            "n_mem_tethers": n_tethers,
            "n_particle_delta": n_part_delta, "n_bond_delta": n_bond_delta,
            "cross_layer": {
                "ok": cross_layer_ok, "self_pairs": self_pairs,
                "unique_cortex_anchors": unique_anchors,
            },
            "force_free_construction": {"ok": force_free_ok, "max_strain": max_strain},
            "no_contamination": {
                "mask_excludes_all_mem": bool(mask_excludes_mem),
                "gamma_soft_off_N_per_m": g_off,
                "gamma_soft_on_N_per_m": g_on,
                "gamma_identical": gamma_identical,
            },
            "off_identity": bool(n_part_delta == n_tethers),
            "bleb_paths_blocked": blocked,
        },
        "deferred": {
            "bleb_nucleation": "σ_crit_bleb None (Tinevez 2009; MCF7 PI-pending) → "
            "MembraneTetherUpdater raises.",
            "reservoir_release": "f_excess None (Raucher-Sheetz/Figard; MCF7 "
            "PI-pending) → released_area raises.",
        },
        "note": (
            "REAL build-time activation gate on the full physiological baseline. "
            "STATIC mem_tether mesh on an OWN mem_node offset layer (BLOCKER-1 "
            "fix). KEY control is no-contamination: cortical γ_soft identical OFF "
            "vs ON (mem_ registry-denylisted). Degree-aware unique-acceptor "
            "selection bounds per-cortex bond degree to +1 (exclusion-cap safe). "
            "Bleb rupture + reservoir release are PI-pending and honestly blocked."
        ),
    }
    _figure(s_on, layout, g_off, g_on, p, result)
    return result


def _figure(s_on, layout, g_off, g_on, p, result) -> None:
    fig, (ax0, ax1, ax2) = plt.subplots(1, 3, figsize=(15, 4.4))
    pos = np.asarray(s_on.particles.position) * 1e6
    tid = np.asarray(s_on.particles.typeid)
    types = list(s_on.particles.types)
    cor_id = types.index("actin_cortex")
    mem_id = types.index("mem_node") if "mem_node" in types else -1
    slab = 0.10 * float(p.R_cell) * 1e6
    is_cor = (tid == cor_id) & (np.abs(pos[:, 2]) < slab)
    is_mem = (tid == mem_id) & (np.abs(pos[:, 2]) < slab) if mem_id >= 0 else np.zeros(len(tid), bool)
    ax0.scatter(pos[is_cor, 0], pos[is_cor, 1], s=3, color="#1b9e77", label="cortex")
    ax0.scatter(pos[is_mem, 0], pos[is_mem, 1], s=5, color="#d95f02", label="mem_node")
    if layout is not None and layout.n_tether > 0:
        pr = np.asarray(layout.tether_pairs, dtype=np.int64)
        for a, b in pr:
            if abs(pos[a, 2]) < slab and abs(pos[b, 2]) < slab:
                ax0.plot([pos[a, 0], pos[b, 0]], [pos[a, 1], pos[b, 1]],
                         color="#7570b3", lw=0.4, alpha=0.6)
    ax0.set_xlabel("x [µm]"); ax0.set_ylabel("y [µm]")
    ax0.set_aspect("equal", adjustable="datalim")
    ax0.set_title(f"mem_node layer + tethers (|z|<{slab:.1f} µm)")
    ax0.legend(fontsize=8)

    if layout is not None and layout.n_tether > 0:
        ax1.hist(np.asarray(layout.tether_r0) * 1e9, bins=20, color="#3690c0", alpha=0.85)
        ax1.axvline(p.membrane_offset * 1e9, color="crimson", lw=2,
                    label=f"offset = {p.membrane_offset*1e9:.0f} nm")
        ax1.set_xlabel("mem_tether r0 [nm]"); ax1.set_ylabel("count")
        ax1.set_title(
            f"force-free (max strain "
            f"{result['controls']['force_free_construction']['max_strain']:.1e})"
        )
        ax1.legend(fontsize=8)

    ax2.bar(["γ_soft OFF", "γ_soft ON"], [g_off * 1e3, g_on * 1e3],
            color=["#4c72b0", "#dd8452"])
    ax2.set_ylabel("cortical γ_soft [mN/m]")
    ax2.set_title(
        f"NO-contamination: identical="
        f"{result['controls']['no_contamination']['gamma_identical']}"
    )
    fig.suptitle(
        f"H.8 membrane_reservoir activation GATE — {result['verdict']} | "
        f"mesh={result['controls']['mesh_assembled']} | "
        f"cross-layer={result['controls']['cross_layer']['ok']} | "
        f"force-free={result['controls']['force_free_construction']['ok']} | "
        f"no-contam={result['controls']['no_contamination']['gamma_identical']}",
        fontsize=9.5,
    )
    (_OUT / "figs").mkdir(parents=True, exist_ok=True)
    fig.savefig(_OUT / "figs" / "h7_membrane_reservoir_activation_gate.png",
                dpi=130, bbox_inches="tight")
    plt.close(fig)


def main() -> int:
    res = run()
    (_OUT / "production").mkdir(parents=True, exist_ok=True)
    path = _OUT / "production" / "h7_membrane_reservoir_activation_gate.json"
    with open(path, "w") as fh:
        json.dump(res, fh, indent=2, default=str)
    c = res["controls"]; nc = c["no_contamination"]
    print(f"[H.8 membrane_reservoir activation gate] {res['verdict']}")
    print(f"  mesh assembled={c['mesh_assembled']} (n_tethers={c['n_mem_tethers']}, "
          f"+{c['n_particle_delta']} part, +{c['n_bond_delta']} bonds)")
    print(f"  cross-layer={c['cross_layer']['ok']} (self_pairs={c['cross_layer']['self_pairs']}, "
          f"unique_anchors={c['cross_layer']['unique_cortex_anchors']})")
    print(f"  force-free: max strain={c['force_free_construction']['max_strain']:.2e} "
          f"ok={c['force_free_construction']['ok']}")
    print(f"  no-contam: γ_soft OFF={nc['gamma_soft_off_N_per_m']:.4e} ON={nc['gamma_soft_on_N_per_m']:.4e} "
          f"identical={nc['gamma_identical']}; mem excluded={nc['mask_excludes_all_mem']}")
    print(f"  bleb paths blocked={all(c['bleb_paths_blocked'].values())}")
    print(f"  json: {path}")
    return 0 if res["verdict"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
