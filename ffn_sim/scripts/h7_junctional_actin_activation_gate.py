"""H.junctional_actin activation gate — α-catenin/vinculin cadherin↔cortex clutch
(STUB→LIVE 2026-06-09).

Eighth (LAST) compartment graduated under the PI ownership grant. The junctional-
actin belt is the THIRD leg of cell–cell adhesion: cadherin trans-dimers carry the
cell-to-cell load (cadherin_junction), and this belt carries that load from the
cadherin tail into each cell's cortical actin via an explicit α-catenin/vinculin
CATCH clutch (Buckley 2014). Built on the two-cell doublet
(:func:`ffn_sim.cell.doublet.build_cell_doublet` ``with_junctional_actin=True``):
one ``junc_actin`` coupling head per interface cadherin that has a SAME-cell cortex
bead within the catch reach, a permanent ``junc_actin_anchor`` (head↔cadherin) and a
force-free per-r0-binned ``junc_actin_couple_b{i}`` (head↔cortex). Build-time controls:

  BELT ASSEMBLED — junc_actin heads > 0 + junc_actin_anchor + junc_actin_couple_b{i}.
  SAME-CELL      — every coupling joins a head to a cortex bead of the head's OWN
                   cell (cadherin-A heads → cell-A cortex; never cross-cell).
  FORCE-FREE     — anchor bonds at anchor_r0 (strain ~0); coupling bonds per-r0-bin.
  NO-CONTAM      — every junc_actin bond type excluded from the cortical mask
                   (registry-denylisted) → the belt never enters cortical γ.
  CATCH LAW      — catch_off_rate(F) is biphasic (falls to a minimum at F* then
                   rises) — the α-catenin force-stabilised catch signature.

The catch-set constants are PI-CANDIDATES (Buckley x_catch/x_slip SOLID; k_catch0/
k_slip0 ORDER; k_couple/k_anchor/k_on/max_couple_dist/anchor_r0 DERIVED H.3 transfers;
module defaults None → un-anchored build raises). The DYNAMIC catch-slip maintenance
(the JunctionalActinCouplingUpdater as a live Action, the engaged-fraction φ
equilibrium) needs an equilibrated run — DEFERRED. ⚠ Fidelity caveat: the single-
particle cadherin is the ectodomain tip at the interface (~0.5 µm from cortex), so
only tips within the α-catenin reach couple → a SPARSE belt; a faithful dense belt
needs a cadherin-tail particle near the cortex (follow-on; PI queue).

Run:  python ffn_sim/scripts/h7_junctional_actin_activation_gate.py
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

from ffn_sim.cell.doublet import build_cell_doublet
from ffn_sim.cortex.cortical_tension import _is_adhesion_bond_type
from ffn_sim.junction.junctional_actin import catch_off_rate

_OUT = _HERE.parents[1] / "outputs" / "h7"


def run() -> dict:
    d = build_cell_doublet(
        "mcf7_baseline.yaml", n_cad_per_cell=60, seed=1,
        with_junctional_actin=True, run_binder_batches=0,
    )
    sim = d.simulation
    s = sim.state.get_snapshot()
    types = list(s.particles.types)
    tid = np.asarray(s.particles.typeid, dtype=np.int64)
    pos = np.asarray(s.particles.position, dtype=np.float64)
    bn = list(s.bonds.types)
    bg = np.asarray(s.bonds.group, dtype=np.int64)
    bt = np.asarray(s.bonds.typeid, dtype=np.int64)
    p = d.p_junctional_actin
    n_a = d.n_cortex_per_cell

    n_head = int((tid == types.index("junc_actin")).sum()) if "junc_actin" in types else 0
    couple_names = [b for b in bn if b.startswith("junc_actin_couple")]
    junc_bonds = [b for b in bn if b.startswith("junc_actin")]

    # BELT ASSEMBLED.
    belt_ok = (
        n_head > 0
        and "junc_actin" in types
        and "junc_actin_anchor" in bn
        and len(couple_names) > 0
    )

    # FORCE-FREE: anchor bonds at anchor_r0.
    anc = bg[bt == bn.index("junc_actin_anchor")] if "junc_actin_anchor" in bn else np.empty((0, 2), int)
    anc_strain = float("nan")
    if anc.shape[0] > 0:
        La = np.linalg.norm(pos[anc[:, 0]] - pos[anc[:, 1]], axis=1)
        anc_strain = float(np.max(np.abs(La - p.anchor_r0) / p.anchor_r0))
    # coupling bonds force-free vs their per-bin r0.
    from ffn_sim.junction.junctional_actin import junc_actin_couple_bin_rest_lengths
    bin_r0 = junc_actin_couple_bin_rest_lengths(p.n_bins, float(p.max_couple_dist))
    couple_strain = 0.0
    couple_ids = {bn.index(nm): i for i, nm in enumerate(
        [f"junc_actin_couple_b{i}" for i in range(p.n_bins)]) if nm in bn}
    n_couple = 0
    for tid_b, bin_i in couple_ids.items():
        m = bt == tid_b
        if not m.any():
            continue
        seg = bg[m]
        n_couple += int(seg.shape[0])
        Lc = np.linalg.norm(pos[seg[:, 0]] - pos[seg[:, 1]], axis=1)
        couple_strain = max(couple_strain, float(np.max(np.abs(Lc - bin_r0[bin_i]) / max(bin_r0[bin_i], 1e-30))))
    # per-r0-bin coupling: residual ≤ half a bin width (force-free to bin resolution).
    bin_w = float(p.max_couple_dist) / p.n_bins
    force_free_ok = bool(
        np.isfinite(anc_strain) and anc_strain < 1e-6
        and couple_strain <= (0.5 * bin_w / bin_r0.min()) + 1e-9
    )

    # SAME-CELL coupling: head's coupling cortex bead in the same cell as its
    # parent cadherin. Use the junc_actin handle (built per cell).
    ja = d.extras["handles"].get("junc_actin")
    same_cell_ok = True
    cross_cell = 0
    if ja is not None:
        heads = np.asarray(ja["head_tags"], dtype=np.int64)
        cortex_g = np.asarray(ja["cortex_global"], dtype=np.int64)
        cad = np.asarray(ja["coupled_cadherin"], dtype=np.int64)
        a_set = set(d.cell_a_cadherin_tags.tolist())
        for cad_t, cor_t in zip(cad, cortex_g):
            cad_is_a = cad_t in a_set
            cor_is_a = cor_t < n_a
            if cad_is_a != cor_is_a:
                cross_cell += 1
        same_cell_ok = cross_cell == 0

    # NO-CONTAM.
    mask_excludes_junc = all(_is_adhesion_bond_type(t) for t in junc_bonds)

    # CATCH LAW: biphasic k_off(F) with a minimum at F* (catch signature).
    Fgrid = np.linspace(0.0, 60.0e-12, 200)
    koff = np.array([catch_off_rate(p, float(f)) for f in Fgrid])
    imin = int(np.argmin(koff))
    catch_ok = bool(0 < imin < len(Fgrid) - 1 and koff[imin] < koff[0] and koff[-1] > koff[imin])
    F_star = float(Fgrid[imin])

    ok = belt_ok and force_free_ok and same_cell_ok and mask_excludes_junc and catch_ok

    result = {
        "gate": "H.junctional_actin activation (cadherin↔cortex catch clutch)",
        "verdict": "PASS" if ok else "REVIEW",
        "physics_claim": True,
        "params": {
            "k_couple_Npm": float(p.k_couple), "k_anchor_Npm": float(p.k_anchor),
            "anchor_r0_m": float(p.anchor_r0), "max_couple_dist_m": float(p.max_couple_dist),
            "k_catch0_s": float(p.k_catch0), "x_catch_m": float(p.x_catch),
            "k_slip0_s": float(p.k_slip0), "x_slip_m": float(p.x_slip),
            "n_bins": int(p.n_bins),
        },
        "controls": {
            "belt_assembled": bool(belt_ok),
            "n_junc_actin_heads": n_head, "n_coupling_bonds": n_couple,
            "n_interface_cadherins": int(2 * d.p_cadherin.n_cad_per_cell),
            "same_cell_coupling": {"ok": same_cell_ok, "cross_cell": cross_cell},
            "force_free_construction": {
                "ok": force_free_ok, "anchor_max_strain": anc_strain,
                "couple_max_strain": couple_strain,
            },
            "no_contamination": {
                "ok": bool(mask_excludes_junc),
                "mask_excludes_all_junc_actin": bool(mask_excludes_junc),
                "junc_actin_bond_types": junc_bonds[:6],
            },
            "catch_signature": {
                "ok": catch_ok, "F_star_pN": F_star * 1e12,
                "koff0_per_s": float(koff[0]), "koff_min_per_s": float(koff[imin]),
            },
        },
        "deferred": {
            "dynamic_catch_slip_maintenance": "JunctionalActinCouplingUpdater as a "
            "live hoomd Action (engaged-fraction φ equilibrium, break/rebind) — "
            "needs an equilibrated run (BAOAB guard).",
            "dense_belt_fidelity": "single-particle cadherin = ectodomain tip at the "
            "interface (~0.5 µm from cortex) → only tips within the α-catenin reach "
            "couple (sparse belt). A faithful dense belt needs a cadherin-tail "
            "particle near the cortex (follow-on; PI queue).",
        },
        "note": (
            "REAL build-time activation gate on the two-cell doublet with the "
            "junctional-actin belt. The reserved STUB build path is now IMPLEMENTED "
            "(heads + anchor + per-r0-bin coupling). Catch-set constants are "
            "PI-candidates (module default None → un-anchored build raises). KEY "
            "controls: same-cell coupling (no cross-cell), force-free, junc_actin_ "
            "excluded from cortical γ, and the α-catenin catch signature (biphasic "
            "k_off with F*). Dynamic maintenance deferred."
        ),
    }
    _figure(s, d, Fgrid, koff, F_star, result)
    return result


def _figure(s, d, Fgrid, koff, F_star, result) -> None:
    fig, (ax0, ax1, ax2) = plt.subplots(1, 3, figsize=(15, 4.4))
    pos = np.asarray(s.particles.position) * 1e6
    tid = np.asarray(s.particles.typeid)
    types = list(s.particles.types)
    n_a = d.n_cortex_per_cell
    cor = np.flatnonzero(tid == types.index("actin_cortex"))
    a_cor = cor[cor < n_a]; b_cor = cor[cor >= n_a]
    ax0.scatter(pos[a_cor, 0], pos[a_cor, 2], s=1, color="#cfe2f3", label="cell A cortex")
    ax0.scatter(pos[b_cor, 0], pos[b_cor, 2], s=1, color="#fce5cd", label="cell B cortex")
    if "cadherin" in types:
        cad = tid == types.index("cadherin")
        ax0.scatter(pos[cad, 0], pos[cad, 2], s=8, color="#1f77b4", label="cadherin")
    if "junc_actin" in types:
        ja = tid == types.index("junc_actin")
        ax0.scatter(pos[ja, 0], pos[ja, 2], s=12, color="#2ca02c", label="junc_actin head")
    ax0.set_xlabel("x [µm]"); ax0.set_ylabel("z [µm]")
    ax0.set_aspect("equal", adjustable="datalim")
    ax0.set_title(f"{result['controls']['n_junc_actin_heads']} α-catenin coupling heads")
    ax0.legend(fontsize=7)

    c = result["controls"]
    ax1.bar(["heads", "couplings", "cross-cell"],
            [c["n_junc_actin_heads"], c["n_coupling_bonds"], c["same_cell_coupling"]["cross_cell"]],
            color=["#2ca02c", "#1f77b4", "#d62728"])
    ax1.set_ylabel("count")
    ax1.set_title(f"belt (same-cell={c['same_cell_coupling']['ok']}, "
                  f"no-contam={c['no_contamination']['ok']})")

    ax2.plot(Fgrid * 1e12, koff, color="#9467bd", lw=2)
    ax2.axvline(F_star * 1e12, color="crimson", ls="--", lw=1.5, label=f"F*={F_star*1e12:.1f} pN")
    ax2.set_xlabel("force [pN]"); ax2.set_ylabel("k_off(F) [1/s]")
    ax2.set_title("α-catenin CATCH signature (biphasic)")
    ax2.legend(fontsize=8)
    fig.suptitle(
        f"H.junctional_actin activation GATE — {result['verdict']} | "
        f"belt={c['belt_assembled']} | same-cell={c['same_cell_coupling']['ok']} | "
        f"force-free={c['force_free_construction']['ok']} | catch={c['catch_signature']['ok']}",
        fontsize=9.5,
    )
    (_OUT / "figs").mkdir(parents=True, exist_ok=True)
    fig.savefig(_OUT / "figs" / "h7_junctional_actin_activation_gate.png",
                dpi=130, bbox_inches="tight")
    plt.close(fig)


def main() -> int:
    res = run()
    (_OUT / "production").mkdir(parents=True, exist_ok=True)
    path = _OUT / "production" / "h7_junctional_actin_activation_gate.json"
    with open(path, "w") as fh:
        json.dump(res, fh, indent=2, default=str)
    c = res["controls"]
    print(f"[H.junctional_actin activation gate] {res['verdict']}")
    print(f"  belt assembled={c['belt_assembled']} (heads={c['n_junc_actin_heads']}, "
          f"couplings={c['n_coupling_bonds']} of {c['n_interface_cadherins']} cadherins)")
    print(f"  same-cell={c['same_cell_coupling']['ok']} (cross_cell={c['same_cell_coupling']['cross_cell']})")
    print(f"  force-free: anchor strain={c['force_free_construction']['anchor_max_strain']:.2e} "
          f"couple strain={c['force_free_construction']['couple_max_strain']:.3f} ok={c['force_free_construction']['ok']}")
    print(f"  no-contam={c['no_contamination']['ok']} (junc_actin excluded)")
    print(f"  catch signature={c['catch_signature']['ok']} (F*={c['catch_signature']['F_star_pN']:.1f} pN)")
    print(f"  json: {path}")
    return 0 if res["verdict"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
