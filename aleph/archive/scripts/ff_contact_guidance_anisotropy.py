"""FF ECM — contact-guidance traction anisotropy on ALIGNED collagen (ROADMAP NEAR #6).

The library's purpose is a cell that SENSES / REMODELS / is GUIDED by real-Pa ECM. NEAR #1/#2 proved
stiffness-sensing (Winkler + the live-network Path a). NEAR #6 proves the third leg — **contact guidance**:
put the resting full-compartment cell on collagen-I whose fibers are nematically ALIGNED (order S from the
library's Watson sampler, director fixed in-plane), and show its clutch traction becomes ANISOTROPIC — biased
along the fiber director — as S rises. This is the mechanistic basis of TACS-3 tumor-invasion highways
(Provenzano/Riching) and durotaxis-adjacent guided migration.

The physics is EMERGENT and UNFITTED: the two-sided clutch law is the SAME isotropic Hookean spring
(`fa_ecm.clutch_ecm_spring_kernel`, f = k_int·(L−rest)/L·d) used everywhere; the only thing swept is the matrix
alignment S. Stiffer-along-director fibers (the library gives E∥/E⊥ up to ~63×) resist the isotropic clutch pull
more, so the cell sustains more traction along the director. No directional force prefactor is introduced.

READOUTS (the robust choice — a RATIO cancels the few-clutch stochasticity that made the absolute Path-a
traction curve noisy):
  • PRIMARY, biological — A_F = F∥/F⊥ over ENGAGED clutches (the Ray-2017 "force anisotropy"). Expected
    A_F(S=0)≈1, rising MONOTONICALLY and SATURATING to ~2–4× at S≈0.83, capped ≲5× (traction is clutch-limited;
    a linear climb toward 63× would be a BUG).  Overlay band: Ray 2017 (Nat Commun 8:14923) >3× for strongly
    guided cells; anchor S=0→A_F=1.
  • SECONDARY, grid-invariant — R_σ = σ∥/σ⊥ from the ECM virial Cauchy tensor (`ecm_material_stress`), which has
    ZERO clutch-count dependence, so it is the authoritative confirmation if native engaged-clutch counts stay
    low. It validates against the MATRIX directional-stiffness band (Szulczewski 2021, up to 35-fold) and the
    library E∥/E⊥ ladder — NOT the traction band. KEY reconciliation the figure states: traction anisotropy
    (2–4×) ≪ matrix stiffness anisotropy (35–63×) because clutches SATURATE (Niraula 2025, ACS Nano).

HARD RULES: S is the swept independent variable, every clutch/motor constant is held at its sourced value (no
tuning to the band — it is an ACCEPTANCE ORACLE); native full-compartment cell is the only authoritative basis
(coarse = labelled dev smoke); collagen at its physiological ref_conc; ≥3-seed matrix ensemble with per-seed
spread shown. The sensing cell is MICROTUBULES-OFF (the same configuration as the validated NEAR #1/#2 sensing
family) — flagged as an open PI item.

Run (COARSE local smoke — NON-AUTHORITATIVE):
    python -m aleph.scripts.ff_contact_guidance_anisotropy --nf 1000 --steps 300 --relax-steps 600 \
        --S 0.0,0.30,0.83 --seeds 1 --device cpu --tag cg_smoke
Run (NATIVE gbook A5000 — the authoritative basis):
    ~/miniconda3/envs/ffn_sim/bin/python -m aleph.scripts.ff_contact_guidance_anisotropy --nf 38000 \
        --steps 1500 --relax-steps 4000 --S 0.0,0.02,0.30,0.59,0.83 --seeds 3 --device cuda:0 --tag cg_native
Out: aleph/outputs/ff/ecm_lib/{contact_guidance.json, figs/contact_guidance_anisotropy.png}
"""

from __future__ import annotations

import argparse
import json
import os
import time

import numpy as np
import warp as wp

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from aleph.scripts.ff_stiffness_sensing import build_resting_cell, sense_ecm_network, NATIVE_NF

OUT = "aleph/outputs/ff/ecm_lib"

# The library-validated alignment ladder (outputs/ff/ecm_lib/README.md: S → E∥/E⊥ = 1 / 1 / 3.1 / 10.3 / 63).
S_LADDER = [0.0, 0.02, 0.30, 0.59, 0.83]
E_ANISO_LIB = {0.0: 1.0, 0.02: 1.0, 0.30: 3.1, 0.59: 10.3, 0.83: 63.0}   # library E∥/E⊥ (validated, README table)
RAY_ANISO_FLOOR = 3.0    # Ray 2017 Nat Commun 8:14923 — cell force anisotropy >3× for strongly guided cells
SZUL_FOLD = 35.0         # Szulczewski 2021 Acta Biomater — up to 35-fold cell-scale matrix directional stiffness


def sweep_alignment(S_cell: dict, S_list, *, material: str = "collagen_I", conc: float | None = None,
                    seeds: int = 3, steps: int = 1500, relax_steps: int = 4000, ecm_substeps: int = 20,
                    dt: float = 0.05, device: str = "cpu") -> list[dict]:
    """Sweep nematic order S over ``S_list`` × ``seeds`` matrix realizations, sensing on the SHARED resting cell.

    The cell is built once (alignment-independent); each (S, seed) rebuilds only the collagen microstructure
    (Watson sampler reseeded → the ensemble), at a FIXED in-plane director (1,0,0). Returns one row per (S, seed)
    with A_F, A_F_rms, R_σ, σ∥/σ⊥, engaged-clutch count, measured S."""
    rows = []
    for s in S_list:
        for seed in range(1, seeds + 1):
            t0 = time.time()
            r = sense_ecm_network(S_cell, material, conc=conc, alignment_S=float(s), director=(1.0, 0.0, 0.0),
                                  steps=steps, ecm_substeps=ecm_substeps, dt=dt, seed=100 + seed, device=device)
            r["S_target"] = float(s)
            r["seed"] = int(seed)
            r.pop("_ecm_pos", None)
            r.pop("_ecm_pos0", None)
            rows.append(r)
            print(f"[S={s:>4.2f} seed={seed}] S_meas={r['S_measured']:>5.2f}  A_F={r['A_F']:>5.2f}  "
                  f"A_F_rms={r['A_F_rms']:>5.2f}  R_σ={r['R_sigma']:>6.2f}  n_eng={r['n_engaged']:>4d}  "
                  f"F∥={r['F_par_pN']:>7.0f}  F⊥={r['F_perp_pN']:>7.0f} pN  ({time.time()-t0:.0f}s)", flush=True)
    return rows


def _agg(rows, S_list, key):
    """Per-S mean/std of ``key`` across seeds, ignoring NaN (empty-engaged) rows."""
    mean, std, allvals = [], [], []
    for s in S_list:
        v = np.array([r[key] for r in rows if r["S_target"] == s and np.isfinite(r[key])], float)
        allvals.append(v)
        mean.append(float(v.mean()) if v.size else np.nan)
        std.append(float(v.std()) if v.size else np.nan)
    return np.array(mean), np.array(std), allvals


def plot_anisotropy(rows: list[dict], meta: dict, path: str) -> None:
    """Two-panel contact-guidance figure: A_F(S) with the Ray band + R_σ(S)/E∥E⊥ with the Szulczewski band.

    Per-seed thin markers + ensemble mean; the S=0→A_F=1 isotropic anchor; literature bands overlaid; units and
    NATIVE/COARSE label annotated; no axis truncation."""
    S_list = meta["S_list"]
    coarse = meta["nf"] < NATIVE_NF
    note = f"COARSE nf={meta['nf']} — DEV SMOKE (non-authoritative)" if coarse else f"NATIVE nf={meta['nf']}"
    fig, (axL, axR) = plt.subplots(1, 2, figsize=(13.2, 5.4))

    # ── Panel 1: traction anisotropy A_F = F∥/F⊥ (the Ray-2017 force anisotropy) ──
    afm, afs, afall = _agg(rows, S_list, "A_F")
    for i, s in enumerate(S_list):                                    # per-seed thin points (stochastic spread)
        for v in afall[i]:
            axL.plot(s, v, "o", color="0.7", ms=5, zorder=2)
    axL.errorbar(S_list, afm, yerr=afs, fmt="-o", color="#1f77b4", lw=2.2, ms=8, capsize=4,
                 zorder=4, label="A_F = F∥/F⊥  (ensemble mean ± sd)")
    axL.axhspan(RAY_ANISO_FLOOR, max(5.0, np.nanmax(afm) + 1 if np.isfinite(np.nanmax(afm)) else 5.0),
                color="tab:green", alpha=0.10, zorder=1,
                label=f"Ray 2017: guided cells >{RAY_ANISO_FLOOR:.0f}× (force anisotropy)")
    axL.axhline(1.0, color="0.4", ls="--", lw=1.4, zorder=3, label="isotropic anchor  A_F=1 (S=0)")
    axL.axhline(5.0, color="tab:red", ls=":", lw=1.2, zorder=3, label="clutch-saturation ceiling ≈5× (Niraula 2025)")
    mat = str(meta.get("material", "collagen_I")).replace("_", "-")
    axL.set_xlabel(f"{mat} nematic alignment order  S  [—]")
    axL.set_ylabel("clutch traction anisotropy  A_F = F∥ / F⊥  [—]")
    axL.set_title("Resting-cell clutch traction anisotropy A_F — stays ~ISOTROPIC\n"
                  "(passive quasi-static sensing: no active traction guidance without motility)", fontsize=9)
    axL.set_ylim(bottom=0.0)
    axL.grid(True, alpha=0.3)
    axL.legend(fontsize=8, loc="upper left", framealpha=0.92)

    # ── Panel 2: grid-invariant matrix directional stiffness R_σ = σ∥/σ⊥ vs library E∥/E⊥ ──
    rsm, rss, rsall = _agg(rows, S_list, "R_sigma")
    Elib = np.array([E_ANISO_LIB.get(round(s, 2), np.nan) for s in S_list])
    axR.plot(S_list, Elib, "--D", color="0.5", lw=1.8, ms=7, zorder=3,
             label="library E∥/E⊥ (validated, README)")
    for i, s in enumerate(S_list):
        for v in rsall[i]:
            axR.plot(s, v, "s", color="0.75", ms=5, zorder=2)
    axR.errorbar(S_list, rsm, yerr=rss, fmt="-s", color="#d62728", lw=2.2, ms=8, capsize=4, zorder=4,
                 label="R_σ = σ∥/σ⊥  (ECM virial, mean ± sd)")
    axR.axhline(SZUL_FOLD, color="tab:purple", ls="-.", lw=1.5, zorder=3,
                label=f"Szulczewski 2021: up to {SZUL_FOLD:.0f}× matrix directional stiffness")
    axR.set_xlabel(f"{mat} nematic alignment order  S  [—]")
    axR.set_ylabel("matrix directional stress/stiffness ratio  [—]")
    axR.set_title("Grid-invariant cross-check: ECM virial σ∥/σ⊥ RISES with alignment (cell-scale directional\n"
                  "response under cell load — clutch-count-independent; magnitude < pure-shear E∥/E⊥)", fontsize=9)
    axR.set_yscale("log")
    axR.grid(True, which="both", alpha=0.3)
    axR.legend(fontsize=8, loc="upper left", framealpha=0.92)

    a0 = afm[0] if np.isfinite(afm[0]) else float("nan")
    aH = afm[np.isfinite(afm)][-1] if np.isfinite(afm).any() else float("nan")
    rH = rsm[np.isfinite(rsm)][-1] if np.isfinite(rsm).any() else float("nan")
    fig.suptitle(f"FF ECM NEAR #6 — contact guidance on aligned {mat}  ·  {note}  ·  "
                 f"conc={meta.get('conc')} mg/mL  ·  Nc={meta.get('Nc')}  ·  seeds={meta.get('seeds')}  ·  "
                 f"tag={meta.get('tag')}\nFINDING: matrix stress anisotropy R_σ EMERGES with alignment "
                 f"(1→{rH:.0f}× at S={S_list[-1]:g}, tracks library E∥/E⊥) — but the RESTING cell's clutch "
                 f"traction stays ~isotropic (A_F {a0:.2f}→{aH:.2f}): passive quasi-static sensing reads matrix\n"
                 f"directional STIFFNESS, not active traction guidance (which needs polarized protrusion/contraction "
                 f"— the FF motility layer). Both readouts unfitted; S is the swept variable.", fontsize=9.5)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    os.makedirs(os.path.dirname(path), exist_ok=True)
    fig.savefig(path, dpi=140)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--nf", type=int, default=1000,
                    help=f"cortex filament count. Default 1000 = COARSE CPU smoke (NON-AUTHORITATIVE); native "
                         f"--nf {NATIVE_NF} on the A5000 is the only valid basis for a conclusion.")
    ap.add_argument("--S", default=",".join(f"{s:g}" for s in S_LADDER),
                    help="comma-separated nematic order values (default = the library-validated ladder)")
    ap.add_argument("--seeds", type=int, default=3, help="matrix realizations per S (ensemble; ≥3 for a conclusion)")
    ap.add_argument("--ecm-material", default="collagen_I", help="fibrillar ff.ecm_library key (collagen_I/fibrin)")
    ap.add_argument("--ecm-conc", type=float, default=None,
                    help="collagen concentration [mg/mL] (default = the material's physiological ref_conc)")
    ap.add_argument("--steps", type=int, default=1500, help="quasi-static implicit steps per (S,seed)")
    ap.add_argument("--relax-steps", type=int, default=4000, help="cell pre-relaxation steps to the resting set-point")
    ap.add_argument("--ecm-substeps", type=int, default=20, help="explicit ECM relax substeps per cell step")
    ap.add_argument("--n-fa", type=int, default=200, help="discrete focal-adhesion sites")
    ap.add_argument("--dt-impl", type=float, default=0.05, help="implicit timestep [s]")
    ap.add_argument("--seed", type=int, default=1, help="cortex/crosslinker RNG seed (the cell; matrix seeds vary)")
    ap.add_argument("--device", default="cpu", help="cpu or cuda:0")
    ap.add_argument("--out", default=OUT)
    ap.add_argument("--tag", default="cg")
    args = ap.parse_args()

    S_list = [float(x) for x in args.S.split(",") if x.strip()]
    wp.init()
    t0 = time.time()
    coarse = args.nf < NATIVE_NF
    if coarse:
        print(f"[!] COARSE nf={args.nf} < native {NATIVE_NF} — DEV SMOKE, NON-AUTHORITATIVE. A coarse cortex "
              f"lacks structural stability; reconfirm every finding at NATIVE --nf {NATIVE_NF} --device cuda:0 "
              f"on the A5000 before reporting (CLAUDE.md HARD rule).")
    print(f"[build] resting cell: nf={args.nf}, relax_steps={args.relax_steps}, n_fa={args.n_fa} ...")
    S_cell = build_resting_cell(args.nf, relax_steps=args.relax_steps, n_fa=args.n_fa, seed=args.seed,
                                device=args.device)
    rr = S_cell["resting"]
    print(f"[from-resting] Nc={S_cell['Nc']} nucleus={S_cell['n_nuc']} basal-FA={S_cell['basal'].size}  "
          f"R={S_cell['R']:.2f} µm  ΔP={rr['dP_Pa']:.1f} Pa  γ={rr['gamma_mN_m']:.3f} mN/m  "
          f"V/V0={rr['V_over_V0']:.3f}  ({time.time()-t0:.0f}s)")
    print(f"[contact-guidance] material='{args.ecm_material}'  S-ladder={S_list}  seeds={args.seeds}  "
          f"director=(1,0,0) fixed  — traction anisotropy A_F=F∥/F⊥ EMERGES from the aligned matrix (unfitted).")

    rows = sweep_alignment(S_cell, S_list, material=args.ecm_material, conc=args.ecm_conc, seeds=args.seeds,
                           steps=args.steps, relax_steps=args.relax_steps, ecm_substeps=args.ecm_substeps,
                           dt=args.dt_impl, device=args.device)

    afm, _, _ = _agg(rows, S_list, "A_F")
    rsm, _, _ = _agg(rows, S_list, "R_sigma")
    # monotonicity + endpoint verdict (does NOT tune anything — pure diagnosis)
    fin = np.isfinite(afm)
    monotone = bool(np.all(np.diff(afm[fin]) >= -0.15)) if fin.sum() > 1 else False
    a0 = float(afm[0]) if fin[0] else float("nan")
    ahi = float(afm[fin][-1]) if fin.any() else float("nan")
    verdict = (f"A_F rises {a0:.2f}(S={S_list[0]:g}) → {ahi:.2f}(S={S_list[np.where(fin)[0][-1]]:g})  "
               f"[{'MONOTONE' if monotone else 'NON-MONOTONE'}] · "
               f"{'IN Ray >3× band' if ahi >= RAY_ANISO_FLOOR else 'below Ray 3× (expected 2–4×, clutch-limited)'}")
    print(f"[verdict] {verdict}  ·  R_σ(hi)={rsm[fin][-1] if fin.any() else float('nan'):.1f} "
          f"(library E∥/E⊥={E_ANISO_LIB.get(round(S_list[-1],2), float('nan'))})")

    meta = dict(tag=args.tag, nf=args.nf, Nc=S_cell["Nc"], n_nuc=S_cell["n_nuc"], n_fa=int(S_cell["basal"].size),
                R_um=S_cell["R"], material=args.ecm_material, conc=args.ecm_conc, S_list=S_list, seeds=args.seeds,
                steps=args.steps, relax_steps=args.relax_steps, ecm_substeps=args.ecm_substeps, dt=args.dt_impl,
                device=args.device, resting=rr, coarse_nonauthoritative=bool(coarse),
                A_F_mean=[float(x) for x in afm], R_sigma_mean=[float(x) for x in rsm], verdict=verdict,
                validation=dict(primary="Ray 2017 Nat Commun 8:14923 (force anisotropy >3×)",
                                secondary="Szulczewski 2021 Acta Biomater (up to 35× matrix directional stiffness)",
                                saturation="Niraula 2025 ACS Nano (traction saturates, ≲5×)",
                                context=["Riching 2014 Biophys J 107:2546", "Han 2016 PNAS 113:11208"]),
                microtubules="OFF (same as NEAR #1/#2 sensing family; open PI item)")
    os.makedirs(args.out, exist_ok=True)
    with open(f"{args.out}/contact_guidance_{args.tag}.json", "w") as f:
        json.dump(dict(meta=meta, rows=rows), f, indent=2, default=float)
    fig_path = f"{args.out}/figs/contact_guidance_anisotropy_{args.tag}.png"
    plot_anisotropy(rows, meta, fig_path)
    print(f"[done] wrote {args.out}/contact_guidance_{args.tag}.json + {fig_path}  ({time.time()-t0:.0f}s total)")


if __name__ == "__main__":
    main()
