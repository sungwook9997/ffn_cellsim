"""A1: emergent A/A0(R0) across Bare/Pre/Lam4 via ligand→ACTIVE-traction (the faithful driver).

L2.6 mapped the ligand condition onto the COARSE *passive* substrate-adhesion depth and the
three curves barely separated (a cohesive MCF7 caps regardless). This A1 driver maps each
condition onto the cell's **active edge-traction** through the per-species integrin-clutch
kinetics (``spheroid.ligand_traction.resolve_ligand_traction``: col-I vs laminin-111 occupancy
φ·F_s, density axis flagged) — the validated spreading knob — and reports the three emergent
A/A0 = a + b/R + c/R² curves. The passive substrate wall is held COMMON across conditions so
the ONLY thing that differs is the ligand-set active traction (isolates the A1 mechanism).

The relative traction ordering is set by MEASURED kinetics, never guessed to match the poster:
col-I clutch (Bare/Pre) is stronger than the laminin clutch (Lam4, 0.61×); Bare<Pre is the
flagged density axis (pV4D4 col-I adsorption density is a literature gap → collaborator).
PI A/A0 is overlay-only (A2); A1 reports what the mechanism produces (the single-cell↔
collective laminin split is expected and stated, REPORT §A1).

Usage: python -m ffn_sim.scripts.layer2_ligand_traction_sweep
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import yaml

from ffn_sim.archive.hoomd_legacy.spheroid.cadherin_bonds import resolve_cadherin
from ffn_sim.archive.hoomd_legacy.spheroid.ligand_traction import resolve_ligand_traction
from ffn_sim.archive.hoomd_legacy.spheroid.observables import effective_radius
from ffn_sim.archive.hoomd_legacy.spheroid.params import resolve_layer2, resolve_proliferation
from ffn_sim.archive.hoomd_legacy.spheroid.proliferation import run_growth_pooled
from ffn_sim.archive.hoomd_legacy.spheroid.substrate import resolve_substrate
from ffn_sim.validation.oracles.spheroid.aa0_law import aa0_model, fit_aa0

_ROOT = Path(__file__).resolve().parents[1]
_CFG = _ROOT / "configs" / "layer2_cbm.yaml"
_FIG_DIR = _ROOT / "outputs" / "layer2" / "figs"
_CKPT = _ROOT / "outputs" / "layer2" / "ligand_traction_sweep.checkpoint.json"

_CONDITIONS = ("Bare", "Pre", "Lam4")


def main(argv: list[str] | None = None) -> int:
    cfg = yaml.safe_load(_CFG.read_text())
    resolved = resolve_layer2(cfg)
    prolif = resolve_proliferation(cfg, resolved)
    cad = resolve_cadherin(resolved)
    # passive substrate held COMMON (reference adhesion) across conditions — isolates the
    # active-traction mechanism (L2.6 already showed the passive axis is nearly flat).
    sub = resolve_substrate(resolved, adhesion_ratio=1.0)

    # resolve the three ligand→traction anchors (mechanism, not fit)
    ligs = {c: resolve_ligand_traction(c) for c in _CONDITIONS}
    print("[A1] ligand → active edge-traction (clutch kinetics):")
    for c in _CONDITIONS:
        r = ligs[c]
        print(f"   {c:5s} {r.ligand:12s} φ={r.engaged_fraction:.3f} F_s={r.F_s*1e12:4.1f}pN "
              f"strength={r.clutch_strength:.3f} density={r.density_factor:.2f} → "
              f"f_traction={r.f_traction*1e9:.2f} nN{'  (proxy)' if r.proxy else ''}")

    n_list = [120, 250, 400, 600]  # 4 R0 → non-degenerate 3-param fit (1 dof); first-pass CPU
    n_seeds = 2
    total_time = 2.0 * prolif.cycle_time_mean

    ckpt = json.loads(_CKPT.read_text()) if _CKPT.exists() else {}
    for c in _CONDITIONS:
        f_tr = ligs[c].f_traction
        for n in n_list:
            key = f"{c}:{n}"
            if key in ckpt:
                continue
            r0r, aar = [], []
            n_ejected = 0
            for s in range(n_seeds):
                res = run_growth_pooled(
                    resolved, prolif, n_cells_init=n, total_time=total_time,
                    epoch_steps=1200, settle_steps=1000, seed=3000 + s, max_cells=3000,
                    cohesion="catch", cad=cad, substrate=sub, f_traction=f_tr,
                )
                if res.get("ejected"):
                    n_ejected += 1  # surfaced (no silent truncation)
                r0r.append(effective_radius(res["a0_core"]))
                aar.append(float(res["area_core_over_a0"][-1]))
            ckpt[key] = {"R0": float(np.mean(r0r)), "aa": float(np.mean(aar)),
                         "aa_sd": float(np.std(aar)), "f_traction": f_tr, "ejected": n_ejected}
            _CKPT.write_text(json.dumps(ckpt, indent=2))
            ej = f"  ⚠EJECTED {n_ejected}/{n_seeds}" if n_ejected else ""
            print(f"  {c:5s}(f={f_tr*1e9:.2f}nN) N0={n:4d} R0={ckpt[key]['R0']*1e6:5.1f}µm "
                  f"A/A0={ckpt[key]['aa']:.2f}±{ckpt[key]['aa_sd']:.2f}{ej}", flush=True)

    # ---- per-condition fits ----
    fits = {}
    print("\n[A1] emergent a/b/c per condition (mechanism → law; NOT fit to PI):")
    for c in _CONDITIONS:
        R0 = np.array([ckpt[f"{c}:{n}"]["R0"] for n in n_list])
        AA = np.array([ckpt[f"{c}:{n}"]["aa"] for n in n_list])
        fit = fit_aa0(R0, AA)
        fits[c] = {"R0": R0, "AA": AA, "fit": fit, "f": ligs[c].f_traction}
        print(f"   {c:5s}: A/A0 = {fit['a']:.3f} + ({fit['b']*1e6:.2f} µm)/R + "
              f"({fit['c']*1e12:.1f} µm²)/R²   r²={fit['r_squared']:.3f}   "
              f"(f_traction={ligs[c].f_traction*1e9:.2f} nN)")

    # spread separation at a common mid R0 (does the active mechanism separate the conditions?)
    R_mid = float(np.median([fits[c]["R0"] for c in _CONDITIONS]))
    aa_mid = {c: float(aa0_model(np.array([R_mid]), **{k: fits[c]["fit"][k] for k in "abc"})[0])
              for c in _CONDITIONS}
    sep = max(aa_mid.values()) - min(aa_mid.values())
    print(f"\n[A1] A/A0 at R0≈{R_mid*1e6:.0f} µm: " +
          "  ".join(f"{c}={aa_mid[c]:.2f}" for c in _CONDITIONS) +
          f"   → separation Δ={sep:.3f} "
          f"({'SEPARATED' if sep > 0.10 else 'still mild'} vs L2.6 passive)")

    try:
        _make_figure(fits, ligs, R_mid, sep)
    except Exception as exc:  # noqa: BLE001 — viz best-effort
        print(f"[viz] skipped figure: {exc}")
    return 0


def _make_figure(fits, ligs, R_mid, sep) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    _FIG_DIR.mkdir(parents=True, exist_ok=True)
    colors = {"Bare": "tab:gray", "Pre": "tab:blue", "Lam4": "tab:green"}
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(12.8, 5.2))

    allR = np.concatenate([fits[c]["R0"] for c in _CONDITIONS])
    Rg = np.linspace(allR.min() * 0.92, allR.max() * 1.05, 200)
    for c in _CONDITIONS:
        d = fits[c]
        fit = d["fit"]
        proxy = "  (proxy)" if ligs[c].proxy else ""
        ax.plot(d["R0"] * 1e6, d["AA"], "o", ms=8, color=colors[c], zorder=3)
        ax.plot(Rg * 1e6, aa0_model(Rg, fit["a"], fit["b"], fit["c"]), "-", color=colors[c],
                label=f"{c} (f={d['f']*1e9:.2f} nN, {ligs[c].ligand}{proxy})  r²={fit['r_squared']:.2f}")
    ax.axhline(1.0, color="black", ls="--", lw=0.8, label="A/A₀ = 1 (no spreading)")
    ax.set_xlabel("initial effective radius R₀ (µm)")
    ax.set_ylabel("A / A₀  (connected-core spread ratio)")
    ax.set_title("A1 — emergent MCF7 spreading law across ligand conditions\n"
                 "ligand → ACTIVE edge-traction (integrin-clutch kinetics); passive wall common",
                 fontsize=9)
    ax.legend(fontsize=8)

    # right: the mechanism — resolved f_traction per condition (the anchored ordering)
    xs = np.arange(len(_CONDITIONS))
    fvals = [ligs[c].f_traction * 1e9 for c in _CONDITIONS]
    bars = ax2.bar(xs, fvals, color=[colors[c] for c in _CONDITIONS], alpha=0.85)
    ax2.axhline(3.0, color="crimson", ls=":", lw=1.2, label="B1 stable ceiling (3 nN)")
    for x, c in zip(xs, _CONDITIONS):
        ax2.text(x, fvals[x] + 0.05, f"{ligs[c].ligand}\nstrength {ligs[c].clutch_strength:.2f}\n"
                 f"dens {ligs[c].density_factor:.2f}", ha="center", va="bottom", fontsize=7)
    ax2.set_xticks(xs); ax2.set_xticklabels(_CONDITIONS)
    ax2.set_ylabel("resolved edge-traction f (nN)")
    ax2.set_ylim(0, 3.4)
    ax2.set_title("mechanism: ligand → f_traction = T_ref·density·(φ·F_s)\n"
                  "col-I clutch > laminin (0.61×, measured); Bare<Pre density (flagged)",
                  fontsize=9)
    ax2.legend(fontsize=8)

    fig.tight_layout()
    out = _FIG_DIR / "fig_layer2_ligand_traction_conditions.png"
    fig.savefig(out, dpi=130); plt.close(fig)
    print(f"[viz] wrote {out}")


if __name__ == "__main__":
    sys.exit(main())
