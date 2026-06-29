"""GATE-J — explicit E-cadherin trans-dimer junction ACTIVATION (LIVE 2026-06-09).

Seventh compartment graduated EXPERIMENTAL->LIVE under the PI ownership grant, and
the first MULTICELL one. The explicit E-cadherin trans-dimer adherens junction is
intrinsically TWO-cell, so its activation needs the new two-cell assembler
:func:`ffn_sim.archive.hoomd_legacy.cell.doublet.build_cell_doublet`: two cortex shells in one box,
facing across an interface, with cadherins seeded on each cell's facing cap and the
trans-dimers SEEDED pre-bound (the physiological engaged-junction baseline — a bare
junction would need ~1e3 binder batches to engage, invalid as a baseline; the
catch-bond binder then MAINTAINS them). Build-time GATE-J controls:

  DOUBLET ASSEMBLED — two cortex cells (2·n_cortex) + 2·n_cad cadherins +
                      cadherin_trans + cadherin_anchor bond types.
  TRANS ENGAGED     — n_trans seeded trans-dimers > 0, EVERY one A↔B (no intra-cell
                      cadherin bond) — a genuine two-cell junction.
  FORCE-FREE        — each seeded trans-dimer born at exactly r0_trans (strain ~0);
                      cohesion emerges from dynamics, not construction pre-stress.
  NO-CONTAM         — every cadherin_ bond type (cadherin_trans / cadherin_anchor)
                      is excluded from the cortical-tension mask (registry-denylisted)
                      → the junction never enters cortical γ_soft.
  BINDER ATTACHED   — the catch-slip CadherinTransJunctionUpdater is attached, ready
                      to maintain the junction (break/rebind via Rakshit k_off(F)).

The DYNAMIC catch-slip maintenance run (the bound-fraction φ ≈ k_on/(k_on+k_off)
equilibrium + the Iturri ~6.5 nN ensemble de-adhesion observable) needs an
equilibrated run (the raw run trips the BAOAB guard — equilibration prelude) and is
DEFERRED, exactly as the other compartments' active observables.

Run:  python ffn_sim/scripts/h7_cadherin_junction_activation_gate.py
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

from ffn_sim.archive.hoomd_legacy.cell.doublet import build_cell_doublet
from ffn_sim.archive.hoomd_legacy.cortex.cortical_tension import (
    _is_adhesion_bond_type,
    measure_cortical_tension,
)

_OUT = _HERE.parents[1] / "outputs" / "h7"


def run() -> dict:
    d = build_cell_doublet(
        "mcf7_baseline.yaml", n_cad_per_cell=60, seed=1, run_binder_batches=0,
    )
    sim = d.simulation
    s = sim.state.get_snapshot()
    types = list(s.particles.types)
    tid = np.asarray(s.particles.typeid, dtype=np.int64)
    pos = np.asarray(s.particles.position, dtype=np.float64)
    bn = list(s.bonds.types)
    bg = np.asarray(s.bonds.group, dtype=np.int64)
    bt = np.asarray(s.bonds.typeid, dtype=np.int64)
    p = d.p_cadherin

    n_cad = int((tid == types.index("cadherin")).sum()) if "cadherin" in types else 0
    cad_bonds = [t for t in bn if t.startswith("cadherin_")]

    # DOUBLET ASSEMBLED.
    doublet_ok = (
        d.n_cortex_per_cell > 0
        and "cadherin" in types
        and n_cad == 2 * int(p.n_cad_per_cell)
        and "cadherin_trans" in bn and "cadherin_anchor" in bn
    )

    # TRANS ENGAGED: seeded dimers, every one A↔B.
    trans = bg[bt == bn.index("cadherin_trans")] if "cadherin_trans" in bn else np.empty((0, 2), int)
    a_set = set(d.cell_a_cadherin_tags.tolist())
    b_set = set(d.cell_b_cadherin_tags.tolist())
    cross = sum(
        1 for u, v in trans
        if ((u in a_set and v in b_set) or (u in b_set and v in a_set))
    )
    intra = int(trans.shape[0]) - cross
    engaged_ok = bool(trans.shape[0] > 0 and cross == trans.shape[0] and intra == 0)

    # FORCE-FREE: trans-dimer construction strain.
    max_strain = float("nan")
    if trans.shape[0] > 0:
        L = np.linalg.norm(pos[trans[:, 0]] - pos[trans[:, 1]], axis=1)
        max_strain = float(np.max(np.abs(L - p.r0_trans) / p.r0_trans))
    force_free_ok = bool(np.isfinite(max_strain) and max_strain < 1e-6)

    # NO-CONTAM: every cadherin_ bond type is excluded by the cortical mask.
    mask_excludes_cad = all(_is_adhesion_bond_type(t) for t in cad_bonds)
    sim.run(0)
    g = measure_cortical_tension(sim, R_cell=d.p_cortex.R_cell, p_enclosed_volume=None)
    g_soft = float(g.get("gamma_soft_N_per_m", g.get("gamma_soft", float("nan"))))
    no_contam_ok = bool(mask_excludes_cad and np.isfinite(g_soft))

    binder_ok = d.cadherin_binder is not None

    ok = doublet_ok and engaged_ok and force_free_ok and no_contam_ok and binder_ok

    result = {
        "gate": "GATE-J cadherin_junction activation (two-cell doublet)",
        "verdict": "PASS" if ok else "REVIEW",
        "physics_claim": True,
        "params": {
            "n_cad_per_cell": int(p.n_cad_per_cell),
            "k_trans_Npm": float(p.k_trans), "r0_trans_m": float(p.r0_trans),
            "r_bind_m": float(p.r_bind), "k_on_s": float(p.k_on),
            "batch_steps": int(p.batch_steps),
            "n_cad_scale_bridge": float(p.n_cad),
        },
        "controls": {
            "doublet_assembled": bool(doublet_ok),
            "n_cortex_per_cell": int(d.n_cortex_per_cell),
            "n_cadherin_total": n_cad,
            "trans_engaged": {
                "ok": engaged_ok, "n_trans_dimers": int(trans.shape[0]),
                "cross_AB": int(cross), "intra_cell": int(intra),
            },
            "force_free_construction": {"ok": force_free_ok, "max_strain": max_strain},
            "no_contamination": {
                "ok": no_contam_ok,
                "mask_excludes_all_cadherin": bool(mask_excludes_cad),
                "cadherin_bond_types": cad_bonds,
                "gamma_soft_N_per_m": g_soft,
            },
            "binder_attached": bool(binder_ok),
        },
        "deferred": {
            "dynamic_catch_slip_maintenance": "bound-fraction φ ≈ k_on/(k_on+k_off) "
            "equilibrium + the Iturri ~6.5 nN ensemble de-adhesion observable — "
            "needs an equilibrated run (raw run trips the BAOAB guard).",
            "iturri_SourceEvidence": "Iturri 2020 de-adhesion anchor unregistered in "
            "the Notion SoT — register before citing in a deliverable (PI queue).",
        },
        "note": (
            "REAL build-time GATE-J on the new two-cell doublet assembler. The "
            "engaged junction is SEEDED pre-bound (physiological baseline; a bare "
            "junction needs ~1e3 binder batches to engage). KEY controls: every "
            "trans-dimer is A↔B (genuine two-cell junction), force-free at r0_trans, "
            "and cadherin_ bonds are excluded from cortical γ (registry-denylisted). "
            "Dynamic catch-slip maintenance + the Iturri de-adhesion observable are "
            "deferred (equilibrated run)."
        ),
    }
    _figure(s, d, trans, g_soft, result)
    return result


def _figure(s, d, trans, g_soft, result) -> None:
    fig, (ax0, ax1, ax2) = plt.subplots(1, 3, figsize=(15, 4.4))
    pos = np.asarray(s.particles.position) * 1e6
    tid = np.asarray(s.particles.typeid)
    types = list(s.particles.types)
    cor = tid == types.index("actin_cortex")
    cad = tid == types.index("cadherin") if "cadherin" in types else np.zeros(len(tid), bool)
    n_a = d.n_cortex_per_cell
    # colour cortex by cell (A vs B) using tag order.
    cor_idx = np.flatnonzero(cor)
    a_cor = cor_idx[cor_idx < n_a]; b_cor = cor_idx[cor_idx >= n_a]
    ax0.scatter(pos[a_cor, 0], pos[a_cor, 2], s=1, color="#9ecae1", label="cell A cortex")
    ax0.scatter(pos[b_cor, 0], pos[b_cor, 2], s=1, color="#fdae6b", label="cell B cortex")
    ax0.scatter(pos[cad, 0], pos[cad, 2], s=8, color="#d62728", label="cadherin")
    ax0.set_xlabel("x [µm]"); ax0.set_ylabel("z [µm]")
    ax0.set_aspect("equal", adjustable="datalim")
    ax0.set_title("Two-cell doublet (x = interface axis)")
    ax0.legend(fontsize=7)

    te = result["controls"]["trans_engaged"]
    if trans.shape[0] > 0:
        L = np.linalg.norm(
            np.asarray(s.particles.position)[trans[:, 0]]
            - np.asarray(s.particles.position)[trans[:, 1]], axis=1)
        ax1.hist(L * 1e9, bins=15, color="#d62728", alpha=0.85)
        ax1.axvline(d.p_cadherin.r0_trans * 1e9, color="navy", lw=2,
                    label=f"r0_trans={d.p_cadherin.r0_trans*1e9:.0f} nm")
        ax1.set_xlabel("cadherin_trans length [nm]"); ax1.set_ylabel("count")
        ax1.set_title(f"{te['n_trans_dimers']} trans-dimers ({te['cross_AB']} A↔B, "
                      f"{te['intra_cell']} intra) force-free")
        ax1.legend(fontsize=8)

    ax2.bar(["A↔B", "intra-cell"], [te["cross_AB"], te["intra_cell"]],
            color=["#2ca02c", "#d62728"])
    ax2.set_ylabel("trans-dimer count")
    ax2.set_title(f"genuine 2-cell junction (intra=0); γ_soft excl. cadherin_\n"
                  f"no-contam={result['controls']['no_contamination']['ok']}")
    fig.suptitle(
        f"GATE-J cadherin_junction activation — {result['verdict']} | "
        f"doublet={result['controls']['doublet_assembled']} | "
        f"engaged={te['ok']} | force-free={result['controls']['force_free_construction']['ok']} | "
        f"no-contam={result['controls']['no_contamination']['ok']}",
        fontsize=9.5,
    )
    (_OUT / "figs").mkdir(parents=True, exist_ok=True)
    fig.savefig(_OUT / "figs" / "h7_cadherin_junction_activation_gate.png",
                dpi=130, bbox_inches="tight")
    plt.close(fig)


def main() -> int:
    res = run()
    (_OUT / "production").mkdir(parents=True, exist_ok=True)
    path = _OUT / "production" / "h7_cadherin_junction_activation_gate.json"
    with open(path, "w") as fh:
        json.dump(res, fh, indent=2, default=str)
    c = res["controls"]; te = c["trans_engaged"]; nc = c["no_contamination"]
    print(f"[GATE-J cadherin_junction activation] {res['verdict']}")
    print(f"  doublet assembled={c['doublet_assembled']} "
          f"(2×{c['n_cortex_per_cell']} cortex, {c['n_cadherin_total']} cadherins)")
    print(f"  trans engaged={te['ok']} (n={te['n_trans_dimers']}, A↔B={te['cross_AB']}, intra={te['intra_cell']})")
    print(f"  force-free: max strain={c['force_free_construction']['max_strain']:.2e} "
          f"ok={c['force_free_construction']['ok']}")
    print(f"  no-contam={nc['ok']} (cadherin_ excluded={nc['mask_excludes_all_cadherin']}, "
          f"γ_soft={nc['gamma_soft_N_per_m']:.4e})")
    print(f"  binder attached={c['binder_attached']}")
    print(f"  json: {path}")
    return 0 if res["verdict"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
