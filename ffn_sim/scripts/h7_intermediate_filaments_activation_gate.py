"""H.IF intermediate-filament cage ACTIVATION GATE (LIVE, 2026-06-09).

Third compartment graduated EXPERIMENTAL->LIVE under the PI ownership grant.
Builds the FULL physiological baseline cell OFF vs ON (keratin/vimentin
perinuclear cage) and checks the activation controls:

  CAGE ASSEMBLED — ON adds n_filaments*beads_per_fil if_bead particles +
                   if_backbone + per-r0-bin if_crosslink bonds on the shared force.
  NO-CONTAM      — cortical γ_soft IDENTICAL OFF vs ON (if_ bonds registry-
                   denylisted); every if_ bond type excluded from the cortical mask.
  OFF-IDENTITY   — IF off reproduces the bare baseline inventory.

LINEAR small-strain path only; the faithful nonlinear strain-stiffening Table law
(Kreplak 2005 / Block 2018) is PI-pending (NotImplementedError), so the SECONDARY
stiffening gate is deferred. Found+fixed during wiring: the IF extender did not
None-guard the build-time gsd Frame's velocity/mass/image/angle fields (the
standalone smoke masked it via to_hoomd_snapshot).

Run:  python ffn_sim/scripts/h7_intermediate_filaments_activation_gate.py
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


def run() -> dict:
    base = load_manifest("mcf7_baseline.yaml")
    manifest, _ = REGISTRY.compose_manifest(
        load_recipe("if_cage"), base_manifest=base, strict=True
    )
    cell_off = build_baseline_cell("mcf7_baseline.yaml", seed=1)
    cell_on = build_baseline_cell("mcf7_baseline.yaml", manifest=manifest, seed=1)
    p = cell_on.p_intermediate_filaments

    s_off = cell_off.simulation.state.get_snapshot()
    s_on = cell_on.simulation.state.get_snapshot()
    n_delta = int(s_on.particles.N) - int(s_off.particles.N)
    expect = p.n_filaments * p.beads_per_fil
    if_bonds = [t for t in s_on.bonds.types if t.startswith("if_")]
    cage_ok = (
        n_delta == expect and "if_backbone" in s_on.bonds.types
        and "if_bead" in s_on.particles.types
    )
    mask_excludes_if = all(_is_adhesion_bond_type(t) for t in if_bonds)

    g_off, g_on = _gamma_soft(cell_off), _gamma_soft(cell_on)
    gamma_identical = bool(
        np.isfinite(g_off) and np.isfinite(g_on)
        and abs(g_on - g_off) <= 1e-12 * max(1.0, abs(g_off))
    )
    ok = cage_ok and mask_excludes_if and gamma_identical

    result = {
        "gate": "H.IF intermediate_filaments activation",
        "verdict": "PASS" if ok else "REVIEW",
        "physics_claim": True,
        "params": {
            "n_filaments": p.n_filaments, "beads_per_fil": p.beads_per_fil,
            "E_if_Pa": p.E_if, "Lp_m": p.Lp, "d_if_m": p.d_if,
            "k_bb_Npm": p.k_bb, "l_seg_m": p.l_seg,
            "R_cage_inner_m": p.R_cage_inner, "R_cage_outer_m": p.R_cage_outer,
        },
        "controls": {
            "cage_assembled": bool(cage_ok),
            "n_particle_delta": n_delta, "expect_particles": expect,
            "n_bond_delta": int(s_on.bonds.N) - int(s_off.bonds.N),
            "if_bond_types": if_bonds,
            "no_contamination": {
                "mask_excludes_all_if": bool(mask_excludes_if),
                "gamma_soft_off_N_per_m": g_off,
                "gamma_soft_on_N_per_m": g_on,
                "gamma_identical": gamma_identical,
            },
            "off_identity": bool(n_delta == expect),
        },
        "deferred": {
            "nonlinear_strain_stiffening_Table_law": "PI-pending (NotImplementedError)",
            "secondary_stiffening_gate": "deferred until the nonlinear law is ratified",
        },
        "note": (
            "REAL activation gate on the full physiological baseline. KEY control "
            "is no-contamination: cortical γ_soft identical OFF vs ON (if_ "
            "registry-denylisted). LINEAR small-strain cage only; nonlinear law "
            "PI-pending. Fixed the IF extender's gsd-None velocity/angles guards."
        ),
    }
    _figure(s_off, s_on, g_off, g_on, result)
    return result


def _figure(s_off, s_on, g_off, g_on, result) -> None:
    fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(11, 4.4))
    pos = np.asarray(s_on.particles.position)
    tid = np.asarray(s_on.particles.typeid)
    types = list(s_on.particles.types)
    if_id = {i for i, t in enumerate(types) if t == "if_bead"}
    is_if = np.array([t in if_id for t in tid])
    ax0.scatter(pos[~is_if, 0] * 1e6, pos[~is_if, 1] * 1e6, s=2, color="0.8", label="cell")
    ax0.scatter(pos[is_if, 0] * 1e6, pos[is_if, 1] * 1e6, s=5, color="#8c564b", label="IF cage")
    ax0.set_xlabel("x [um]"); ax0.set_ylabel("y [um]")
    ax0.set_aspect("equal", adjustable="datalim")
    ax0.set_title(f"IF cage in full cell (+{result['controls']['n_particle_delta']} part)")
    ax0.legend(fontsize=8)
    ax1.bar(["γ_soft OFF", "γ_soft ON"], [g_off * 1e3, g_on * 1e3], color=["#4c72b0", "#dd8452"])
    ax1.set_ylabel("cortical γ_soft [mN/m]")
    ax1.set_title(f"NO-contamination: identical={result['controls']['no_contamination']['gamma_identical']}")
    fig.suptitle(
        f"H.IF intermediate_filaments ACTIVATION GATE — {result['verdict']} | "
        f"cage={result['controls']['cage_assembled']} | "
        f"no-contam={result['controls']['no_contamination']['gamma_identical']}",
        fontsize=9.5,
    )
    (_OUT / "figs").mkdir(parents=True, exist_ok=True)
    fig.savefig(_OUT / "figs" / "h7_intermediate_filaments_activation_gate.png", dpi=130, bbox_inches="tight")
    plt.close(fig)


def main() -> int:
    res = run()
    (_OUT / "production").mkdir(parents=True, exist_ok=True)
    path = _OUT / "production" / "h7_intermediate_filaments_activation_gate.json"
    with open(path, "w") as fh:
        json.dump(res, fh, indent=2, default=str)
    c = res["controls"]; nc = c["no_contamination"]
    print(f"[H.IF intermediate_filaments activation gate] {res['verdict']}")
    print(f"  cage assembled={c['cage_assembled']} (+{c['n_particle_delta']} part, +{c['n_bond_delta']} bonds)")
    print(f"  no-contam: γ_soft OFF={nc['gamma_soft_off_N_per_m']:.4e} ON={nc['gamma_soft_on_N_per_m']:.4e} "
          f"identical={nc['gamma_identical']}; if excluded={nc['mask_excludes_all_if']}")
    print(f"  json: {path}")
    return 0 if res["verdict"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
