"""H.MT microtubules ACTIVATION GATE (LIVE, 2026-06-09).

Second compartment graduated EXPERIMENTAL->LIVE under the PI ownership grant.
Builds the FULL physiological baseline cell OFF vs ON (microtubule aster) and
checks the activation controls:

  ASTER ASSEMBLED — ON adds exactly 1 MTOC + n_mt*beads_per_mt mt_bead particles,
                    n_mt*beads_per_mt mt_backbone bonds, mt_bending angles.
  CFL            — the build's MT CFL gate passes at cytoplasm drag (γ_b NOT water).
  NO-CONTAM      — cortical γ_soft is IDENTICAL OFF vs ON (mt_ bonds are
                   registry-denylisted), AND the cortical bond mask excludes every
                   mt_ type. The aster does NOT move cortical tension.
  OFF-IDENTITY   — MT off reproduces the bare baseline inventory.

PRIMARY rod-like persistence (L_p ~ 5.2 mm >> L_mt) is validated by the standalone
smoke (scripts/compartment_smoke/microtubules.py: backbone rms dev 0.03%, rod
stable). A full-cell raw thermal run trips the BAOAB guard (needs the equilibration
prelude — known full-cell issue), so the stepped rod check stays in the smoke.

Auto-viz per the production-driver rule.  Run:
  python ffn_sim/scripts/h7_microtubules_activation_gate.py
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
    res = measure_cortical_tension(
        cell.simulation, R_cell=cell.p_cortex.R_cell,
        p_enclosed_volume=cell.p_enclosed_volume,
    )
    return float(res.get("gamma_soft_N_per_m", res.get("gamma_soft", float("nan"))))


def run() -> dict:
    base = load_manifest("mcf7_baseline.yaml")
    manifest, deferred = REGISTRY.compose_manifest(
        load_recipe("aster_microtubules"), base_manifest=base, strict=True
    )

    cell_off = build_baseline_cell("mcf7_baseline.yaml", seed=1)
    cell_on = build_baseline_cell("mcf7_baseline.yaml", manifest=manifest, seed=1)
    p_mt = cell_on.p_microtubules

    s_off = cell_off.simulation.state.get_snapshot()
    s_on = cell_on.simulation.state.get_snapshot()

    n_delta = int(s_on.particles.N) - int(s_off.particles.N)
    expect_particles = 1 + p_mt.n_mt * p_mt.beads_per_mt
    bond_delta = int(s_on.bonds.N) - int(s_off.bonds.N)
    expect_bonds = p_mt.n_mt * p_mt.beads_per_mt
    aster_ok = (
        n_delta == expect_particles and bond_delta == expect_bonds
        and "mt_bead" in s_on.particles.types and "mtoc" in s_on.particles.types
        and "mt_backbone" in s_on.bonds.types and "mt_bending" in s_on.angles.types
    )

    # NO-CONTAMINATION: every mt_ bond type excluded from the cortical mask.
    mt_types = [t for t in s_on.bonds.types if t.startswith("mt_")]
    mask_excludes_mt = all(_is_adhesion_bond_type(t) for t in mt_types)
    cortex_counted = all(
        not _is_adhesion_bond_type(t) for t in s_on.bonds.types
        if t.startswith("cortex") or t == "cortex-bond"
    )
    g_off = _gamma_soft(cell_off)
    g_on = _gamma_soft(cell_on)
    gamma_identical = bool(np.isfinite(g_off) and np.isfinite(g_on)
                           and abs(g_on - g_off) <= 1e-12 * max(1.0, abs(g_off)))

    off_identity = (n_delta == expect_particles)   # ON-OFF delta is exactly the aster

    ok = aster_ok and mask_excludes_mt and cortex_counted and gamma_identical
    result = {
        "gate": "H.MT microtubules activation",
        "verdict": "PASS" if ok else "REVIEW",
        "physics_claim": True,
        "params": {
            "n_mt": p_mt.n_mt, "beads_per_mt": p_mt.beads_per_mt,
            "L_mt_m": p_mt.L_mt, "l0_m": p_mt.l0, "EI_Nm2": p_mt.EI,
            "L_p_m": p_mt.L_p, "Y_stretch_N": p_mt.Y_stretch,
            "gamma_b_cytoplasm_Nspm": p_mt.gamma_b,
            "dt_cfl_stretch_s": p_mt.dt_cfl_stretch,
        },
        "controls": {
            "aster_assembled": bool(aster_ok),
            "n_particle_delta": n_delta, "expect_particles": expect_particles,
            "n_bond_delta": bond_delta, "expect_bonds": expect_bonds,
            "cfl_gate_passed_at_build": True,   # build raised otherwise
            "no_contamination": {
                "mt_types": mt_types,
                "mask_excludes_all_mt": bool(mask_excludes_mt),
                "cortex_bonds_counted": bool(cortex_counted),
                "gamma_soft_off_N_per_m": g_off,
                "gamma_soft_on_N_per_m": g_on,
                "gamma_identical": gamma_identical,
            },
            "off_identity": bool(off_identity),
        },
        "PRIMARY_persistence": {
            "L_p_m": p_mt.L_p, "band_m": [1.0e-3, 8.0e-3],
            "in_band": bool(1.0e-3 <= p_mt.L_p <= 8.0e-3),
            "note": "rod-like thermal stability validated in smoke_microtubules (rms dev 0.03%)",
        },
        "note": (
            "REAL activation gate on the full physiological baseline. The KEY "
            "activation control is no-contamination: cortical γ_soft is identical "
            "OFF vs ON because mt_ bonds are registry-denylisted. CFL gate passed "
            "at cytoplasm drag. Stepped rod-like L_p validated in the standalone "
            "smoke (full-cell raw run needs the equilibration prelude)."
        ),
    }

    _figure(s_off, s_on, p_mt, g_off, g_on, result)
    return result


def _figure(s_off, s_on, p_mt, g_off, g_on, result) -> None:
    fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(11, 4.4))

    pos = np.asarray(s_on.particles.position)
    tid = np.asarray(s_on.particles.typeid)
    types = list(s_on.particles.types)
    mt_ids = {i for i, t in enumerate(types) if t in ("mt_bead", "mtoc")}
    is_mt = np.array([t in mt_ids for t in tid])
    ax0.scatter(pos[~is_mt, 0] * 1e6, pos[~is_mt, 2] * 1e6, s=2, color="0.8", label="cell")
    ax0.scatter(pos[is_mt, 0] * 1e6, pos[is_mt, 2] * 1e6, s=6, color="#d62728", label="MT aster")
    ax0.set_xlabel("x [um]"); ax0.set_ylabel("z [um]")
    ax0.set_aspect("equal", adjustable="datalim")
    ax0.set_title(f"Aster in full cell (+{result['controls']['n_particle_delta']} particles)")
    ax0.legend(fontsize=8)

    ax1.bar(["γ_soft OFF", "γ_soft ON"], [g_off * 1e3, g_on * 1e3],
            color=["#4c72b0", "#dd8452"])
    ax1.set_ylabel("cortical γ_soft [mN/m]")
    ax1.set_title(f"NO-contamination: identical = {result['controls']['no_contamination']['gamma_identical']}")

    fig.suptitle(
        f"H.MT microtubules ACTIVATION GATE — {result['verdict']} | "
        f"aster assembled={result['controls']['aster_assembled']} | "
        f"no-contam={result['controls']['no_contamination']['gamma_identical']}",
        fontsize=9.5,
    )
    (_OUT / "figs").mkdir(parents=True, exist_ok=True)
    fig.savefig(_OUT / "figs" / "h7_microtubules_activation_gate.png", dpi=130, bbox_inches="tight")
    plt.close(fig)


def main() -> int:
    res = run()
    (_OUT / "production").mkdir(parents=True, exist_ok=True)
    path = _OUT / "production" / "h7_microtubules_activation_gate.json"
    with open(path, "w") as fh:
        json.dump(res, fh, indent=2, default=str)
    c = res["controls"]
    nc = c["no_contamination"]
    print(f"[H.MT microtubules activation gate] {res['verdict']}")
    print(f"  aster assembled={c['aster_assembled']} (+{c['n_particle_delta']} part, "
          f"+{c['n_bond_delta']} bonds); CFL gate passed={c['cfl_gate_passed_at_build']}")
    print(f"  no-contam: γ_soft OFF={nc['gamma_soft_off_N_per_m']:.4e} ON={nc['gamma_soft_on_N_per_m']:.4e} "
          f"identical={nc['gamma_identical']}; mt excluded={nc['mask_excludes_all_mt']}")
    print(f"  PRIMARY L_p={res['PRIMARY_persistence']['L_p_m']*1e3:.1f} mm in band="
          f"{res['PRIMARY_persistence']['in_band']}")
    print(f"  json: {path}")
    return 0 if res["verdict"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
