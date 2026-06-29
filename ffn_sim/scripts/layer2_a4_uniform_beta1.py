"""A4: does Lam4's UNIFORM β1 (collective mechanism) flip the emergent ordering to Lam4-highest?

A2 confirmed the single-cell↔collective laminin SPLIT: the PI *collective* ranking is
Lam4 > Pre > Bare, but the A1 *single-cell* active-traction puts Lam4 ≈ Bare (laminin is the
weaker single-cell clutch). A2's implication: Lam4's enhancement is NOT single-cell traction —
it is a COLLECTIVE mechanism rooted in laminin's measured "uniform β1" IF pattern (vs Bare
"diffuse" / Pre "peripheral"; LIGAND_PRESENTATION_MECHANISM.md).

A4 tests that mechanism with ONE faithful change: the β1 distribution sets the traction
SCREENING LENGTH (`ligand_traction.resolve_ligand_traction(uniform_beta1=…)` → Lp). Bare/Pre =
edge-localised (Lp = 11 µm, peripheral β1, the A1 baseline); Lam4 = UNIFORM (Lp ≫ spheroid,
"uniform β1") so every basal cell — not just the rim — transmits traction. The MAGNITUDE stays
the A1 clutch-kinetics value (Lam4 1.53 nN); only the DISTRIBUTION differs. Everything else
(catch cohesion, common substrate, proliferation) is identical to A1 — isolating the mechanism.

Hypothesis: uniform traction engages the interior (∝ N ∝ volume) so it lifts LARGE spheroids
relatively more → the small-size penalty is removed / the curve flattens ("scale-independent",
the documented Lam4 phenotype) and Lam4 rises above Pre → the emergent ordering flips to
Lam4-highest, matching the PI collective ranking that single-cell traction could not. Reported
honestly whether it matches, overshoots, or fails (overlay-only; PI never fit).

Usage: python -m ffn_sim.scripts.layer2_a4_uniform_beta1
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
_A1_CKPT = _ROOT / "outputs" / "layer2" / "ligand_traction_sweep.checkpoint.json"
_A4_CKPT = _ROOT / "outputs" / "layer2" / "a4_uniform_beta1.checkpoint.json"

_CONDITIONS = ("Bare", "Pre", "Lam4")
_COLORS = {"Bare": "tab:gray", "Pre": "tab:blue", "Lam4": "tab:green"}
_N_LIST = [120, 250, 400, 600]
_N_SEEDS = 2


def _run_lam4_uniform() -> dict:
    """Run only the Lam4-UNIFORM-β1 sweep (Bare/Pre are unchanged from A1 → reuse A1 ckpt)."""
    cfg = yaml.safe_load(_CFG.read_text())
    resolved = resolve_layer2(cfg)
    prolif = resolve_proliferation(cfg, resolved)
    cad = resolve_cadherin(resolved)
    sub = resolve_substrate(resolved, adhesion_ratio=1.0)
    lam = resolve_ligand_traction("Lam4")  # uniform β1 → Lp = LP_UNIFORM, f = 1.53 nN
    print(f"[A4] Lam4 UNIFORM β1: f={lam.f_traction*1e9:.2f} nN, Lp={lam.Lp*1e6:.0f} µm "
          f"(uniform={lam.uniform_beta1}); Bare/Pre = A1 edge (reused)")
    ck = json.loads(_A4_CKPT.read_text()) if _A4_CKPT.exists() else {}
    for n in _N_LIST:
        key = str(n)
        if key in ck:
            continue
        r0r, aar = [], []
        n_ej = 0
        for s in range(_N_SEEDS):
            res = run_growth_pooled(
                resolved, prolif, n_cells_init=n, total_time=2.0 * prolif.cycle_time_mean,
                epoch_steps=1200, settle_steps=1000, seed=3000 + s, max_cells=3000,
                cohesion="catch", cad=cad, substrate=sub,
                f_traction=lam.f_traction, Lp=lam.Lp,
            )
            if res.get("ejected"):
                n_ej += 1
            r0r.append(effective_radius(res["a0_core"]))
            aar.append(float(res["area_core_over_a0"][-1]))
        ck[key] = {"R0": float(np.mean(r0r)), "aa": float(np.mean(aar)),
                   "aa_sd": float(np.std(aar)), "ejected": n_ej}
        _A4_CKPT.write_text(json.dumps(ck, indent=2))
        ej = f"  ⚠EJECTED {n_ej}/{_N_SEEDS}" if n_ej else ""
        print(f"  Lam4-uniform N0={n:4d} R0={ck[key]['R0']*1e6:5.1f}µm "
              f"A/A0={ck[key]['aa']:.2f}±{ck[key]['aa_sd']:.2f}{ej}", flush=True)
    return ck


def _load_a1(cond: str) -> tuple[np.ndarray, np.ndarray]:
    ck = json.loads(_A1_CKPT.read_text())
    keys = sorted((k for k in ck if k.startswith(f"{cond}:")), key=lambda k: int(k.split(":")[1]))
    return (np.array([ck[k]["R0"] for k in keys]), np.array([ck[k]["aa"] for k in keys]))


def main(argv: list[str] | None = None) -> int:
    if not _A1_CKPT.exists():
        print("[A4] need the A1 checkpoint first (run layer2_ligand_traction_sweep).")
        return 1
    a4 = _run_lam4_uniform()

    # assemble the three A4 conditions: Bare/Pre = A1 edge, Lam4 = A4 uniform
    R0_lam = np.array([a4[str(n)]["R0"] for n in _N_LIST])
    AA_lam = np.array([a4[str(n)]["aa"] for n in _N_LIST])
    curves = {
        "Bare": _load_a1("Bare"),
        "Pre": _load_a1("Pre"),
        "Lam4": (R0_lam, AA_lam),
    }
    lam_edge = _load_a1("Lam4")  # the A1 edge-localised Lam4, for the A/B

    print("\n[A4] emergent A/A0(R0) — Bare/Pre edge (A1), Lam4 UNIFORM (A4):")
    fits = {}
    for c in _CONDITIONS:
        R0, AA = curves[c]
        fit = fit_aa0(R0, AA)
        fits[c] = fit
        print(f"   {c:5s}: a={fit['a']:.2f} b={fit['b']*1e6:+.0f}µm c={fit['c']*1e12:+.0f}µm² "
              f"r²={fit['r_squared']:.3f}  (A/A0: {AA.min():.2f}–{AA.max():.2f})")

    # mid-R0 ordering (A4) vs A1 edge vs PI collective
    R_mid = float(np.median([curves[c][0] for c in _CONDITIONS]))
    aa_mid = {c: float(aa0_model(np.array([R_mid]), fits[c]["a"], fits[c]["b"], fits[c]["c"])[0])
              for c in _CONDITIONS}
    a4_rank = sorted(_CONDITIONS, key=lambda c: aa_mid[c], reverse=True)
    # size-dependence change for Lam4 (edge → uniform): does the small-size penalty lift?
    lam_edge_slope = np.corrcoef(lam_edge[0], lam_edge[1])[0, 1]
    lam_unif_slope = np.corrcoef(R0_lam, AA_lam)[0, 1]
    print(f"\n[A4] === verdict ===")
    print(f"  ordering at R0≈{R_mid*1e6:.0f}µm: " + "  ".join(f"{c}={aa_mid[c]:.2f}" for c in _CONDITIONS)
          + f"  → {' > '.join(a4_rank)}")
    print(f"  A1 edge ordering = Pre > Lam4 ≳ Bare  |  PI collective = Lam4 > Pre > Bare")
    flipped = a4_rank[0] == "Lam4"
    print(f"  → Lam4 now highest? {flipped}  ({'MATCHES PI collective ordering' if flipped else 'no flip'})")
    print(f"  Lam4 corr(R0,A/A0): edge={lam_edge_slope:+.2f} → uniform={lam_unif_slope:+.2f}  "
          f"(uniform engages interior ∝ volume → small-size penalty {'lifted/flattened' if lam_unif_slope>lam_edge_slope else 'unchanged'})")

    try:
        _make_figure(curves, lam_edge, fits, R_mid, aa_mid, a4_rank)
    except Exception as exc:  # noqa: BLE001
        print(f"[viz] skipped figure: {exc}")
    return 0


def _make_figure(curves, lam_edge, fits, R_mid, aa_mid, a4_rank) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    _FIG_DIR.mkdir(parents=True, exist_ok=True)
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(13.0, 5.3))

    allR = np.concatenate([curves[c][0] for c in _CONDITIONS])
    Rg = np.linspace(allR.min() * 0.92, allR.max() * 1.05, 200)
    for c in _CONDITIONS:
        R0, AA = curves[c]
        fit = fits[c]
        lbl = f"{c} ({'UNIFORM β1' if c == 'Lam4' else 'edge β1'})"
        ax.plot(R0 * 1e6, AA, "o", ms=8, color=_COLORS[c], zorder=3)
        ax.plot(Rg * 1e6, aa0_model(Rg, fit["a"], fit["b"], fit["c"]), "-", color=_COLORS[c], label=lbl)
    # the A1 edge Lam4 (the A/B reference) as faint dashed
    ax.plot(lam_edge[0] * 1e6, lam_edge[1], "x", ms=8, color="darkgreen", alpha=0.6)
    fe = fit_aa0(*lam_edge)
    ax.plot(Rg * 1e6, aa0_model(Rg, fe["a"], fe["b"], fe["c"]), "--", color="darkgreen",
            alpha=0.6, lw=1.2, label="Lam4 EDGE β1 (A1, for A/B)")
    ax.axhline(1.0, color="black", ls="--", lw=0.7)
    ax.set_xlabel("initial effective radius R₀ (µm)")
    ax.set_ylabel("A / A₀  (connected-core spread ratio)")
    ax.set_title("A4 — Lam4 UNIFORM β1 (collective) vs A1 edge β1\n"
                 "uniform traction engages the interior → lifts large spheroids / Lam4 rises",
                 fontsize=9)
    ax.legend(fontsize=8)

    # right: the ordering flip — mid-R0 A/A0, A1-edge vs A4, annotated against the PI ranking
    xs = np.arange(len(_CONDITIONS))
    a1_mid = {}
    for c in _CONDITIONS:
        src = lam_edge if c == "Lam4" else curves[c]
        f = fit_aa0(*src)
        a1_mid[c] = float(aa0_model(np.array([R_mid]), f["a"], f["b"], f["c"])[0])
    w = 0.36
    ax2.bar(xs - w / 2, [a1_mid[c] for c in _CONDITIONS], w, color=[_COLORS[c] for c in _CONDITIONS],
            alpha=0.4, hatch="//", label="A1 edge β1 (single-cell)")
    ax2.bar(xs + w / 2, [aa_mid[c] for c in _CONDITIONS], w, color=[_COLORS[c] for c in _CONDITIONS],
            alpha=0.9, label="A4 (Lam4 uniform β1)")
    ax2.set_xticks(xs); ax2.set_xticklabels(_CONDITIONS)
    ax2.set_ylabel(f"A/A₀ at R₀≈{R_mid*1e6:.0f} µm")
    ax2.set_title(f"ordering: A1 edge = Pre>Lam4≳Bare  →  A4 = {' > '.join(a4_rank)}\n"
                  f"PI collective = Lam4>Pre>Bare  ({'MATCH' if a4_rank[0]=='Lam4' else 'no flip'})",
                  fontsize=9)
    ax2.legend(fontsize=8)

    fig.tight_layout()
    out = _FIG_DIR / "fig_layer2_a4_uniform_beta1.png"
    fig.savefig(out, dpi=130); plt.close(fig)
    print(f"[viz] wrote {out}")


if __name__ == "__main__":
    sys.exit(main())
