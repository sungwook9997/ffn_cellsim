"""L2.6 ligand-condition axis: A/A0(R0) under Bare / Pre / Lam4 (the experiment's variable).

The PI experiment compares three substrate ligand conditions (Bare untreated pV4D4, Pre =
adsorbed col-I, Lam4 = + soluble laminin-111), each giving a different A/A0 = a + b/R + c/R²
curve. In the Layer-2 model the ligand condition sets the cell-substrate adhesion strength
(``substrate.resolve_substrate(adhesion_ratio=…)``, the z=0 wall depth = the integrin-ECM
adhesion the ligand presents). This driver sweeps R0 for each condition and overlays the three
emergent curves — the platform analog of the experiment's condition comparison.

HONEST SCOPE: this maps the ligand condition to the COARSE substrate-adhesion knob
(adhesion_ratio); MCF7 is cohesive (low-invasion) so passive adhesion changes the footprint
only mildly (it forms a 3D cap, not a monolayer — L2.6 finding). The full ligand dependence
(per-species integrin catch-slip kinetics in ``bridge/ligand_species.py`` driving ACTIVE edge
traction) is the faithful upgrade; this is the first-pass adhesion-axis result. PI A/A0 is
overlay-only (never fit).

Usage: python -m ffn_sim.scripts.layer2_ligand_sweep
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import yaml

from ffn_sim.spheroid.cadherin_bonds import resolve_cadherin
from ffn_sim.spheroid.observables import effective_radius
from ffn_sim.spheroid.params import resolve_layer2, resolve_proliferation
from ffn_sim.spheroid.proliferation import run_growth_pooled
from ffn_sim.spheroid.substrate import resolve_substrate
from ffn_sim.validation.oracles.spheroid.aa0_law import aa0_model, fit_aa0

_ROOT = Path(__file__).resolve().parents[1]
_CFG = _ROOT / "configs" / "layer2_cbm.yaml"
_FIG_DIR = _ROOT / "outputs" / "layer2" / "figs"
_CKPT = _ROOT / "outputs" / "layer2" / "ligand_sweep.checkpoint.json"

# Ligand condition → substrate-adhesion ratio (× cell-cell cohesion D_e). Coarse first-pass
# mapping (Bare weakest, Lam4 strongest); flagged modeling axis (see module docstring).
_LIGANDS = {"Bare": 0.5, "Pre": 1.0, "Lam4": 2.0}


def main(argv: list[str] | None = None) -> int:
    cfg = yaml.safe_load(_CFG.read_text())
    resolved = resolve_layer2(cfg)
    prolif = resolve_proliferation(cfg, resolved)
    cad = resolve_cadherin(resolved)
    n_list = [120, 250, 450]
    n_seeds = 2
    total_time = 2.0 * prolif.cycle_time_mean

    ckpt = json.loads(_CKPT.read_text()) if _CKPT.exists() else {}
    for lig, ratio in _LIGANDS.items():
        sub = resolve_substrate(resolved, adhesion_ratio=ratio)
        for n in n_list:
            key = f"{lig}:{n}"
            if key in ckpt:
                continue
            r0r, aar = [], []
            n_ejected = 0
            for s in range(n_seeds):
                res = run_growth_pooled(
                    resolved, prolif, n_cells_init=n, total_time=total_time,
                    epoch_steps=1200, settle_steps=1000, seed=2000 + s, max_cells=3000,
                    cohesion="catch", cad=cad, substrate=sub,
                )
                if res.get("ejected"):
                    n_ejected += 1  # surfaced below (no silent truncation)
                r0r.append(effective_radius(res["a0_core"]))
                aar.append(float(res["area_core_over_a0"][-1]))
            ckpt[key] = {"R0": float(np.mean(r0r)), "aa": float(np.mean(aar)),
                         "aa_sd": float(np.std(aar)), "ejected": n_ejected}
            _CKPT.write_text(json.dumps(ckpt, indent=2))
            ej = f"  ⚠EJECTED {n_ejected}/{n_seeds}" if n_ejected else ""
            print(f"  {lig:5s}(adh×{ratio}) N0={n:4d} R0={ckpt[key]['R0']*1e6:5.1f}µm "
                  f"A/A0={ckpt[key]['aa']:.2f}±{ckpt[key]['aa_sd']:.2f}{ej}", flush=True)

    fits = {}
    for lig in _LIGANDS:
        R0 = np.array([ckpt[f"{lig}:{n}"]["R0"] for n in n_list])
        aa = np.array([ckpt[f"{lig}:{n}"]["aa"] for n in n_list])
        f = fit_aa0(R0, aa)
        fits[lig] = (R0, aa, f)
        print(f"[{lig}] A/A0 = {f['a']:.2f} + ({f['b']*1e6:.1f}µm)/R + ({f['c']*1e12:.0f}µm²)/R²"
              f"  r²={f['r_squared']:.3f}")

    try:
        _figure(fits)
    except Exception as exc:  # noqa: BLE001
        print(f"[viz] skipped: {exc}")
    return 0


def _figure(fits) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    _FIG_DIR.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(7.5, 5.5))
    colors = {"Bare": "gray", "Pre": "steelblue", "Lam4": "crimson"}
    for lig, (R0, aa, f) in fits.items():
        Rg = np.linspace(R0.min() * 0.95, R0.max() * 1.05, 100)
        ax.plot(R0 * 1e6, aa, "o", color=colors[lig], ms=8)
        ax.plot(Rg * 1e6, aa0_model(Rg, f["a"], f["b"], f["c"]), "-", color=colors[lig],
                label=f"{lig}: a+b/R+c/R² (r²={f['r_squared']:.2f})")
    ax.axhline(1.0, color="black", ls="--", lw=0.8)
    ax.set_xlabel("initial effective radius R₀ (µm)")
    ax.set_ylabel("A / A₀ (core, after 2 doublings)")
    ax.set_title("L2.6 — emergent A/A₀(R₀) across ligand conditions (substrate adhesion axis)\n"
                 "Bare / Pre / Lam4 = weak / mid / strong cell-substrate adhesion", fontsize=9)
    ax.legend(fontsize=8)
    fig.tight_layout()
    out = _FIG_DIR / "fig_layer2_ligand_conditions.png"
    fig.savefig(out, dpi=130); plt.close(fig)
    print(f"[viz] wrote {out}")


if __name__ == "__main__":
    sys.exit(main())
