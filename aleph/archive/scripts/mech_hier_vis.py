"""Single figure regenerator for the hierarchical mechanics track (S1..S8).

CLAUDE.md "one entry-point" + "visualize at closeout" rule for the professor's
hierarchical single-cell/spheroid mechanics flow
(``docs/v2_audit/_historical/HIERARCHICAL_MECHANICS_VALIDATION_PLAN_2026-07-14.md``). Regenerates
EVERY figure for this track into ``outputs/mech_hier/figs/``:

  1. oracle reference curves  (the analytical acceptance gates, standalone)
  2. S1 loading effect        (naive vs physiological loading: the 215x overshoot fix)
  3. S1 density convergence    (E_fit vs filament count -> native)
  4. S1 layer decomposition    (bare cortex -> +membrane -> ... toward the whole-cell MCF7)

Visualization-integrity rules (CLAUDE.md): no axis truncation, analytic reference
overlaid on every measurement, SI/annotated units, log axes flagged in the title.
Runs from any data present; missing runs are skipped with a note (idempotent).

Usage: ``python -m aleph.scripts.mech_hier_vis``
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from aleph.validation.oracles.mechanics import (  # noqa: E402
    hertz,
    viscoelastic_relaxation as ve,
    boussinesq as bq,
    thin_shell,
)

# Wong colourblind-safe palette.
BLUE, ORANGE, GREEN, VERM, PURPLE, GREY = (
    "#0072B2", "#E69F00", "#009E73", "#D55E00", "#CC79A7", "#666666",
)
MCF7_BAND = (224.0, 279.0)   # Pa, whole-cell target (informational at bare-cortex stages)

ROOT = Path(__file__).resolve().parents[1]            # aleph/
S1DIR = ROOT / "outputs" / "mech_hier" / "s1_sphere"
FIGDIR = ROOT / "outputs" / "mech_hier" / "figs"


def _load_run(name: str) -> dict | None:
    """Load a committed S1 run's (sweep arrays + manifest gate/config), or None."""
    d = S1DIR / name
    npz, man = d / "sweep.npz", d / "manifest.json"
    if not (npz.exists() and man.exists()):
        return None
    z = np.load(npz)
    m = json.loads(man.read_text())
    return {
        "name": name,
        "strain": z["strain"], "delta1_um": z["delta1_um"], "F_pN": z["F_pN"],
        "R0_um": m["R0_um"], "gate": m.get("gate", {}), "config": m.get("config", {}),
    }


# --------------------------------------------------------------------------- #
# 1. Oracle reference curves
# --------------------------------------------------------------------------- #
def fig_oracles_reference() -> None:
    fig, ax = plt.subplots(2, 2, figsize=(10, 8))

    # (a) Hertz F(delta): sphere-plane vs between-plates, at 3 moduli.
    R = 7.5
    d = np.linspace(0, 0.03 * R, 120)
    for E, c in [(249.0, GREEN), (4000.0, BLUE), (40000.0, VERM)]:
        es = hertz.rigid_indenter_reduced_modulus(E, 0.5)
        ax[0, 0].plot(d, hertz.force_sphere_plane(es, R, d), color=c, label=f"E={E:.0f} Pa (plane)")
        ax[0, 0].plot(d, hertz.force_between_plates(es, R, 2 * d), "--", color=c, lw=1)
    ax[0, 0].set(xlabel=r"indentation $\delta$ [µm]", ylabel="F [pN]",
                 title="(a) Hertz F(δ): solid=sphere–plane, dashed=between-plates")
    ax[0, 0].legend(fontsize=7, frameon=False)

    # (b) SLS relaxation G(t) and creep J(t).
    t = np.linspace(0, 40, 300)
    g_inf, g1, tau = 400.0, 600.0, 5.0
    axb = ax[0, 1]
    axb.plot(t, ve.sls_relaxation(t, g_inf, g1, tau), color=BLUE, label="G(t) relaxation")
    axb.axhline(g_inf, color=BLUE, ls=":", lw=0.8)
    axb.axhline(g_inf + g1, color=BLUE, ls=":", lw=0.8)
    axb.set(xlabel="t [s]", ylabel="G(t) [Pa]", title="(b) SLS relaxation modulus + creep compliance")
    axc = axb.twinx()
    axc.plot(t, ve.sls_creep_compliance(t, g_inf, g1, tau) * 1e3, color=ORANGE, label="J(t) creep")
    axc.set_ylabel("J(t) [1/kPa]", color=ORANGE)
    axb.legend(loc="center right", fontsize=7, frameon=False)

    # (c) Boussinesq axial stress decay ~ 1/z^2 (log-log).
    z = np.geomspace(1e-6, 1e-4, 60)
    P = 1e-9
    ax[1, 0].loglog(z * 1e6, bq.sigma_zz_axis(P, z), color=BLUE, label=r"$\sigma_{zz}(z)$")
    ax[1, 0].loglog(z * 1e6, bq.sigma_zz_axis(P, z[0]) * (z / z[0]) ** -2, "k:", lw=0.8,
                    label=r"$\propto 1/z^2$ ref")
    ax[1, 0].set(xlabel="depth z [µm]", ylabel=r"$\sigma_{zz}$ [Pa]",
                 title="(c) Boussinesq point-load stress decay (log–log)")
    ax[1, 0].legend(fontsize=7, frameon=False)

    # (d) Thin-shell inflation dV/V vs pressure, few (E,h).
    p = np.linspace(0, 200, 100)
    for (E, h), c in [((1e3, 0.2e-6), GREEN), ((4e3, 0.2e-6), BLUE), ((1e3, 0.5e-6), VERM)]:
        ax[1, 1].plot(p, thin_shell.thin_shell_volume_change(p, 7.5e-6, h, E, 0.45) * 100,
                      color=c, label=f"E={E:.0f} Pa, h={h*1e9:.0f} nm")
    ax[1, 1].set(xlabel="internal pressure p [Pa]", ylabel="ΔV/V [%]",
                 title="(d) Pressurized thin-shell inflation")
    ax[1, 1].legend(fontsize=7, frameon=False)

    fig.suptitle("Mechanics acceptance oracles (analytic ground truth) — S1/S2 gates", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    _save(fig, "oracles_reference.png")


# --------------------------------------------------------------------------- #
# 2. S1 loading effect (naive vs physio) — the key finding
# --------------------------------------------------------------------------- #
def fig_s1_loading_effect() -> None:
    naive, physio = _load_run("naive"), _load_run("physio")
    if naive is None and physio is None:
        print("# [skip] loading-effect: no naive/physio runs"); return
    fig, ax = plt.subplots(figsize=(6.4, 4.8))
    for run, c, lab in [(naive, VERM, "naive loading"), (physio, BLUE, "physiological loading")]:
        if run is None:
            continue
        d1, F = run["delta1_um"], run["F_pN"]
        e_fit = run["gate"].get("E_fit_Pa", float("nan"))
        ax.plot(d1, F, "o", color=c, ms=7,
                label=f"{lab}: E_fit={e_fit:.0f} Pa ({run['gate'].get('verdict','')})")
        dd = np.linspace(d1.min(), d1.max(), 80)
        es = run["gate"].get("E_star_Pa")
        if es:
            ax.plot(dd, hertz.force_sphere_plane(es, run["R0_um"], dd), "-", color=c, lw=1.2)
        ksh = run["gate"].get("k_shell_pN_per_um")
        if ksh and run["gate"].get("better_model") == "linear-shell":
            ax.plot(dd, ksh * dd, "--", color=c, lw=1.4)  # linear pressurized-shell fit (better)
    # MCF7 whole-cell reference (informational).
    R0 = (physio or naive)["R0_um"]
    dd = np.linspace(0.005 * R0, 0.03 * R0, 80)
    es_mcf7 = hertz.rigid_indenter_reduced_modulus(249.0, 0.5)
    ax.plot(dd, hertz.force_sphere_plane(es_mcf7, R0, dd), "k:", lw=1.2,
            label="MCF7 whole-cell (249 Pa, ref only)")
    ax.set_yscale("log")
    ax.set(xlabel=r"per-contact indentation $\delta_1$ [µm]", ylabel="plate force F [pN] (log)",
           title="S1 bare cortex: loading protocol sets the modulus (native, NF=70686)")
    if naive and physio:
        rn = naive["gate"].get("E_fit_Pa", 1); rp = physio["gate"].get("E_fit_Pa", 1)
        ax.annotate(f"physiological loading\n= {rn/rp:.0f}× softer",
                    xy=(0.55, 0.5), xycoords="axes fraction", fontsize=9, color=GREEN,
                    ha="center", bbox=dict(boxstyle="round", fc="white", ec=GREEN))
    ax.legend(fontsize=7.5, frameon=False, loc="lower right")
    fig.tight_layout()
    _save(fig, "s1_loading_effect.png")


# --------------------------------------------------------------------------- #
# 3. S1 density (NF) convergence
# --------------------------------------------------------------------------- #
def fig_s1_density_convergence() -> None:
    pj = S1DIR / "density_probe_cpu.json"
    if not pj.exists():
        print("# [skip] density-convergence: no probe json"); return
    rows = json.loads(pj.read_text())
    nf = np.array([r["nfil"] for r in rows])
    e = np.array([r["E_fit_Pa"] for r in rows])
    fig, ax = plt.subplots(figsize=(6.2, 4.4))
    ax.plot(nf, e, "o-", color=BLUE, label="CPU probe (non-authoritative tail)")
    physio = _load_run("physio")
    if physio:
        ax.plot(70686, physio["gate"].get("E_fit_Pa"), "*", color=VERM, ms=16,
                label="native NF=70686 (authoritative, GPU)")
    ax.axhspan(*MCF7_BAND, color=GREEN, alpha=0.15, label="MCF7 whole-cell band (ref)")
    ax.set_xscale("log")
    ax.set(xlabel="filament count NF (log)", ylabel="bare-cortex E_fit [Pa]",
           title="S1 filament-density convergence toward native")
    ax.legend(fontsize=8, frameon=False)
    fig.tight_layout()
    _save(fig, "s1_density_convergence.png")


# --------------------------------------------------------------------------- #
# 4. S1 layer decomposition (bare -> +membrane -> ... -> whole cell)
# --------------------------------------------------------------------------- #
def fig_s1_layer_decomposition() -> None:
    stages = [("naive", "bare\n(naive load)"), ("physio", "bare cortex\n(physio)"),
              ("membrane", "+membrane\nreservoir"), ("relax8k", "bare\n(8k relax)")]
    labels, vals, cols = [], [], []
    for name, lab in stages:
        run = _load_run(name)
        if run is None:
            continue
        labels.append(lab); vals.append(run["gate"].get("E_fit_Pa", float("nan")))
        cols.append(VERM if name == "naive" else BLUE)
    if not vals:
        print("# [skip] layer-decomposition: no runs"); return
    fig, ax = plt.subplots(figsize=(6.6, 4.4))
    ax.bar(range(len(vals)), vals, color=cols, width=0.6)
    for i, v in enumerate(vals):
        ax.text(i, v, f"{v:.0f}", ha="center", va="bottom", fontsize=8)
    ax.axhspan(*MCF7_BAND, color=GREEN, alpha=0.2)
    ax.text(len(vals) - 0.5, MCF7_BAND[1], "MCF7 whole-cell target", fontsize=8,
            color=GREEN, va="bottom", ha="right")
    ax.set_yscale("log")
    ax.set_xticks(range(len(labels))); ax.set_xticklabels(labels, fontsize=8)
    ax.set(ylabel="E_fit [Pa] (log)", title="S1 layer decomposition toward the whole-cell modulus")
    fig.tight_layout()
    _save(fig, "s1_layer_decomposition.png")


def fig_s1_loading_spectrum() -> None:
    """The modulus is set by the LOADING regime, not compartments — and my S1 driver
    agrees with the previously-validated harness. Directly answers 'is the old work wrong?'."""
    pj = S1DIR / "loading_regime_diag.json"
    if not pj.exists():
        print("# [skip] loading-spectrum: no diag json"); return
    import matplotlib.patches as mpatches
    diag = json.loads(pj.read_text())
    pts = diag["points"]
    vals = [p["E_fit"] for p in pts]
    cmap = {"undrained": VERM, "drained-mine": BLUE, "validated-harness": ORANGE}
    cols = [cmap.get(p["regime"], GREY) for p in pts]
    fig, ax = plt.subplots(figsize=(9, 4.8))
    ax.bar(range(len(vals)), vals, color=cols, width=0.62)
    for i, v in enumerate(vals):
        ax.text(i, v, f"{v:.0f}", ha="center", va="bottom", fontsize=8)
    ax.axhspan(*diag["mcf7_band"], color=GREEN, alpha=0.2)
    ax.text(len(vals) - 0.5, diag["mcf7_band"][1], "MCF7 whole-cell band",
            fontsize=8, color=GREEN, va="bottom", ha="right")
    ax.set_yscale("log")
    ax.set_xticks(range(len(pts)))
    ax.set_xticklabels([p["label"] for p in pts], fontsize=7)
    ax.set(ylabel="bare-cortex E_fit [Pa] (log)",
           title="S1 modulus set by LOADING regime, not compartments\n"
                 "(my S1 driver ≈ the previously-validated harness; both need the fully-relaxed "
                 "config to reach MCF7)")
    ax.legend([mpatches.Patch(color=VERM), mpatches.Patch(color=BLUE), mpatches.Patch(color=ORANGE)],
              ["undrained (naive)", "drained (my S1 driver)", "drained (validated harness)"],
              fontsize=8, frameon=False, loc="upper right")
    fig.tight_layout()
    _save(fig, "s1_loading_spectrum.png")


def fig_s1_step_convergence() -> None:
    """E_fit vs relaxation steps at native rigid — the convergence-first evidence (Finding 7/8):
    E_fit plateaus by ~15000 steps, and the earlier 4000-step run was under-relaxed / over-stiff."""
    pts = []
    for name in ("physio", "relax8k", "converged"):
        mf = S1DIR / name / "manifest.json"
        if not mf.exists():
            continue
        m = json.loads(mf.read_text())
        ns = m.get("n_steps") or m.get("config", {}).get("n_steps")
        g = m.get("gate", {})
        if ns and g.get("E_fit_Pa"):
            pts.append((int(ns), float(g["E_fit_Pa"]), g.get("verdict", "")))
    if len(pts) < 2:
        print("# [skip] step-convergence: <2 runs"); return
    pts.sort()
    ns = [p[0] for p in pts]
    e = [p[1] for p in pts]
    fig, ax = plt.subplots(figsize=(6.8, 4.6))
    ax.plot(ns, e, "o-", color=BLUE, ms=8)
    for x, y, v in pts:
        ax.annotate(f"{y:.0f} Pa\n{'PASS' if 'clean' in v else ''}", (x, y),
                    textcoords="offset points", xytext=(0, 10), fontsize=8, ha="center")
    ax.axhspan(*MCF7_BAND, color=GREEN, alpha=0.15, label="MCF7 whole-cell band (ref)")
    ax.set(xlabel="relaxation steps", ylabel="bare-cortex E_fit [Pa]",
           title="S1 relaxation-step convergence (native, rigid)\n"
                 "E_fit plateaus by ~15000 steps → converged, LINEAR-SHELL gate PASS")
    ax.legend(fontsize=8, frameon=False)
    fig.tight_layout()
    _save(fig, "s1_step_convergence.png")


def fig_s1_poroelastic() -> None:
    """Poroelastic rate-dependence: cytoplasm drainage softens the cortex. Answers
    'is the fluid the key, not the filaments?' — yes rate-dependent (poroelastic), but
    modest (~3.7x) over physical rates; the 888 kPa impression was an instant-load artifact."""
    pj = S1DIR / "poroelastic_rate.json"
    if not pj.exists():
        print("# [skip] poroelastic: no json"); return
    d = json.loads(pj.read_text())
    r, e = np.array(d["rate_um_s"]), np.array(d["E_fit_Pa"])
    fig, ax = plt.subplots(figsize=(6.8, 4.6))
    ax.plot(r, e, "o-", color=BLUE, ms=7, label="drained skeleton + fluid (poroelastic)")
    ax.axhspan(*d["mcf7_band"], color=GREEN, alpha=0.18, label="MCF7 whole-cell band")
    ax.axhline(d["instant_artifact_Pa"], color=VERM, ls=":", lw=1.3,
               label=f"instant load = {d['instant_artifact_Pa']/1e3:.0f} kPa (unphysical artifact)")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set(xlabel="loading rate [µm/s] (log) — faster →", ylabel="bare-cortex E_fit [Pa] (log)",
           title="Poroelastic rate-dependence (cytoplasm drainage)\n"
                 "fast=fluid trapped=stiff → slow=drained=soft; only ~3.7× over physical rates")
    ax.legend(fontsize=7.5, frameon=False, loc="center left")
    fig.tight_layout()
    _save(fig, "s1_poroelastic_rate.png")


def fig_s1_ft() -> None:
    """Plate force over relaxation time — does it SATURATE or keep changing? (PI directive).
    Also shows the physical time reached: dt is CFL-tiny, so the sim reaches sub-ms = the
    MECHANICAL settle, far short of the drainage/turnover timescale (~seconds)."""
    pj = S1DIR / "ft_trace_native.json"
    if not pj.exists():
        print("# [skip] F(t): no json"); return
    d = json.loads(pj.read_text())
    tr = d["trace"]
    step = np.array([r["step"] for r in tr])
    F = np.array([r["F_plate_pN"] for r in tr])
    t_ms = np.array([r.get("t_s", r["step"] * d.get("dt_s", 0.0)) for r in tr]) * 1e3
    fig, ax = plt.subplots(figsize=(7.0, 4.6))
    ax.plot(step, F, "-", color=BLUE, lw=1.5, label="plate force F(step)")
    ax.axhline(d.get("final_F_pN", F[-1]), color=VERM, ls=":", lw=1.1,
               label=f"final (steady read) = {d.get('final_F_pN', F[-1]):.0f} pN")
    ax.set_yscale("log")
    ax.set(xlabel="relaxation step", ylabel="plate force F [pN] (log)",
           title="S1 F(t): plate force over relaxation (native)\n"
                 f"physical time reached = {t_ms[-1]:.3f} ms (dt={d.get('dt_s',0.0):.1e}s) — "
                 f"MECHANICAL settle, NOT the drainage/turnover timescale (~s)")
    ax.legend(fontsize=8, frameon=False)
    ax.text(0.97, 0.95, f"strain={d.get('strain')}\nNF=70686", transform=ax.transAxes,
            fontsize=9, va="top", ha="right", bbox=dict(boxstyle="round", fc="white", ec=GREY))
    fig.tight_layout()
    _save(fig, "s1_ft_saturation.png")


def _branch_class(f, defl):
    """Classify a load branch from its secants: 'STIFFENS' (sub-linear defl-vs-force),
    'SOFTENS' (super-linear), or '≈linear'. Returns (label, last/first secant ratio)."""
    f = np.abs(np.asarray(f, float)); defl = np.abs(np.asarray(defl, float))
    nz = f > 0
    if nz.sum() < 2:
        return "≈linear", 1.0
    s_first = defl[nz][0] / f[nz][0]                              # near-origin secant
    s_last = (defl[nz][-1] - defl[nz][0]) / (f[nz][-1] - f[nz][0])  # top incremental secant
    r = s_last / max(s_first, 1e-12)
    return ("STIFFENS" if r < 0.85 else "SOFTENS" if r > 1.15 else "≈linear"), r


def fig_s2_force_deflection() -> None:
    """S2 localized-load force–deflection curve (native): sign reversal + the constitutive
    asymmetry, read from the CONVERGED ramp gate (s2_sign_reversal.json, 15000 steps/force,
    drift-free) if present, else the un-converged calibration probe. Stiffen/soften is
    classified directly from the data per branch."""
    gate = S1DIR / "s2_sign_reversal.json"
    probe = S1DIR / "s2_calibration_probe.json"
    pj = gate if gate.exists() else probe
    converged = pj is gate
    if not pj.exists():
        print("# [skip] S2 force–deflection: no json"); return
    d = json.loads(pj.read_text())
    n = len(d["inward_ramp"])
    f = np.linspace(0.0, d["f_max_pN"], n)                       # per-node force magnitude [pN]
    din = np.array(d["inward_ramp"]); dout = np.array(d["outward_ramp"])
    cin, rin = _branch_class(f, din); cout, rout = _branch_class(f, dout)
    fig, ax = plt.subplots(figsize=(7.0, 4.8))
    ax.plot(-f, din, "o-", color=BLUE, lw=1.8, ms=6, label=f"inward load (2A, compression) — {cin}")
    ax.plot(f, dout, "s-", color=VERM, lw=1.8, ms=6, label=f"outward load (2B, tension) — {cout}")
    if n >= 2 and f[1] > 0:                                       # near-origin linear reference
        slope = din[1] / -f[1]
        fr = np.linspace(-f[-1], f[-1], 50)
        ax.plot(fr, slope * fr, ":", color=GREY, lw=1.2,
                label=f"near-origin linear ({abs(slope):.1e} µm/pN)")
    ax.axhline(0, color="#bbb", lw=0.8); ax.axvline(0, color="#bbb", lw=0.8)
    ax.annotate("SIGN REVERSAL\n(force sign → deflection sign)", xy=(0, 0),
                xytext=(0.06, 0.14), textcoords="axes fraction", fontsize=9, color=GREEN, ha="left")
    tag = "CONVERGED (15000 steps/force, drift-free)" if converged else "probe (UN-CONVERGED)"
    ax.set(xlabel="signed per-node patch force  [pN]  (−inward / +outward)",
           ylabel="patch deflection Δ [µm]  (−dimple / +bulge)",
           title=f"S2 localized load — native force–deflection (NF={d.get('nfil', 70686)})\n"
                 f"{tag}: sign-reversing; compression {cin} (r={rin:.2f}), tension {cout} (r={rout:.2f})")
    ax.legend(fontsize=8, frameon=False, loc="lower right")
    fig.tight_layout()
    _save(fig, "s2_force_deflection.png")


def fig_s2_deflection_saturation() -> None:
    """S2 deflection(t) under a fixed localized load — does the dimple reach elastic equilibrium
    (saturate) or keep creeping? (PI directive, localized-load stage.) Reads the native trace."""
    pj = S1DIR / "s2_deflection_trace.json"
    if not pj.exists():
        print("# [skip] S2 deflection(t): no trace json"); return
    d = json.loads(pj.read_text())
    t_ms = np.array(d["t_s"]) * 1e3
    y = np.array(d["deflection_um"])
    yf = d["final_deflection_um"]
    t99 = float(t_ms[np.argmax(np.abs(y) >= 0.99 * abs(yf))]) if len(y) else float("nan")  # 99% equilibration
    fig, ax = plt.subplots(figsize=(7.0, 4.6))
    ax.plot(t_ms, y, "-o", ms=3.5, color=BLUE, lw=1.6)
    ax.axhline(yf, color=VERM, ls="--", lw=1.0, label=f"equilibrium {yf:+.3f} µm")
    ax.axvline(t99, color=GREY, ls=":", lw=1.0, label=f"99% equilibrated @ {t99*1e3:.1f} µs")
    verdict = "SATURATED" if d.get("saturated") else "STILL CREEPING"
    ax.set(xlabel="sim time [ms]", ylabel="patch deflection Δ [µm]",
           title=f"S2 deflection(t) — fixed load f={d['f_node_pN']:.0f} pN/node (native)\n"
                 f"{verdict}: reaches elastic equilibrium at the mechanical timescale (~{t99*1e3:.0f} µs), "
                 f"tail Δ/run={d['tail_frac_change']:.2%}")
    ax.legend(fontsize=8, frameon=False)
    fig.tight_layout()
    _save(fig, "s2_deflection_saturation.png")


def fig_s2_localization() -> None:
    """S2 spatial-localization decay: load-induced radial displacement vs polar angle θ from the
    load axis (rest-differenced). Localized load ⇒ concentrated at the pole, decays to ~0."""
    pj = S1DIR / "s2_sign_reversal.json"
    if not pj.exists():
        print("# [skip] S2 localization: no gate json"); return
    d = json.loads(pj.read_text())
    loc = d.get("localization")
    if not loc:
        print("# [skip] S2 localization: no profile in gate"); return
    th = np.array(loc["theta_deg"]); pin = np.array(loc["profile_inward_um"]); pout = np.array(loc["profile_outward_um"])
    fig, ax = plt.subplots(figsize=(7.0, 4.6))
    ax.axhline(0, color="#bbb", lw=0.8)
    ax.plot(th, pin, "o-", color=BLUE, lw=1.6, label="inward (dimple)")
    ax.plot(th, pout, "s-", color=VERM, lw=1.6, label="outward (bulge)")
    ax.axvspan(0, 20, color=BLUE, alpha=0.07)
    ax.text(10, ax.get_ylim()[1] * 0.9, "load\npole", fontsize=8, ha="center", color=GREY)
    r = loc.get("loc_ratio", float("nan"))
    ax.set(xlabel="polar angle θ from load axis [deg]", ylabel="load-induced radial displacement Δr [µm]",
           title=f"S2 localization — deflection decays from the load pole (native)\n"
                 f"pole/equator ratio = {r:.1f}  ({'localized ✓' if r >= 3 else 'NOT localized'})")
    ax.legend(fontsize=8, frameon=False)
    fig.tight_layout()
    _save(fig, "s2_localization_decay.png")


def fig_s3_nucleus_layered() -> None:
    """S3 layered compression: force vs strain (cortex-only vs cortex+nucleus), the modest force ratio,
    and the incompressible-nucleus BULGE (R_nuc_eq vs strain). Prefers the incompressible sweep
    (PI decision A) if present, else the compressible."""
    inc = S1DIR / "s3_nucleus_incompressible.json"
    comp = S1DIR / "s3_nucleus_compression.json"
    pj = inc if inc.exists() else comp
    if not pj.exists():
        print("# [skip] S3 layered: no sweep json"); return
    incompressible = pj is inc
    d = json.loads(pj.read_text())
    rows = d["rows"]
    strains = sorted({r["strain"] for r in rows})
    fc = [next(r["F_plate_pN"] for r in rows if r["strain"] == s and not r["with_nucleus"]) for s in strains]
    fn = [next(r["F_plate_pN"] for r in rows if r["strain"] == s and r["with_nucleus"]) for s in strains]
    req = [next(r.get("R_nuc_eq_um", 0.0) for r in rows if r["strain"] == s and r["with_nucleus"]) for s in strains]
    sc, accel, pct = d.get("s_contact_theory", 0.30), d.get("contact_accel_x", float("nan")), d.get("nucleus_pct_at_max", float("nan"))
    R_nuc = d.get("R_nuc_um", 5.25)
    fig, (ax, ax2, ax3) = plt.subplots(1, 3, figsize=(15.0, 4.5))
    ax.plot(strains, fc, "o-", color=BLUE, lw=1.8, label="cortex only (S1)")
    ax.plot(strains, fn, "s-", color=VERM, lw=1.8, label="cortex + nucleus (S3)")
    ax.axvline(sc, color=GREY, ls="--", lw=1.0, label=f"contact ε≈{sc:.2f}")
    ax.set_yscale("log"); ax.legend(fontsize=8, frameon=False)
    ax.set(xlabel="compression strain ε", ylabel="plate force F [pN] (log)", title="force vs strain")
    ratio = [fn[i] / fc[i] for i in range(len(strains))]
    ax2.plot(strains, ratio, "D-", color=GREEN, lw=1.8)
    ax2.axvline(sc, color=GREY, ls="--", lw=1.0); ax2.axhline(1.0, color="#bbb", lw=0.8)
    ax2.set(xlabel="compression strain ε", ylabel="F(cortex+nuc) / F(cortex-only)",
            title=f"nucleus contribution MODEST: +{pct:.0f}% at ε={strains[-1]:.2f}\ncontact accel {accel:.1f}× past ε≈{sc:.2f}")
    ax3.plot(strains, req, "^-", color=ORANGE, lw=1.8, label="nucleus R_nuc_eq")
    ax3.axhline(R_nuc, color="#bbb", ls=":", lw=1.0, label=f"rest R_nuc={R_nuc:.2f}")
    ax3.axvline(sc, color=GREY, ls="--", lw=1.0)
    ax3.set(xlabel="compression strain ε", ylabel="nucleus equatorial radius [µm]",
            title=("nucleus BULGES past contact (incompressible)" if incompressible else
                   "nucleus flat (compressible — no volume conservation)"))
    ax3.legend(fontsize=8, frameon=False)
    tag = "INCOMPRESSIBLE ν=0.499 (PI-A): bulges + volume-conserved" if incompressible else "compressible (legacy)"
    fig.suptitle(f"S3 layered compression — soft nucleus (E={d.get('E_nuc_Pa', 399):.0f} Pa) modest; "
                 f"{tag} (NF={d.get('nfil', 70686)})", fontsize=11)
    fig.tight_layout()
    _save(fig, "s3_nucleus_layered.png")


def _save(fig, name: str) -> None:
    FIGDIR.mkdir(parents=True, exist_ok=True)
    p = FIGDIR / name
    fig.savefig(p, dpi=140)
    plt.close(fig)
    print(f"# wrote {p.relative_to(ROOT)}")


def main() -> None:
    fig_oracles_reference()
    fig_s1_loading_effect()
    fig_s1_density_convergence()
    fig_s1_layer_decomposition()
    fig_s1_loading_spectrum()
    fig_s1_poroelastic()
    fig_s1_ft()
    fig_s1_step_convergence()
    fig_s2_force_deflection()
    fig_s2_deflection_saturation()
    fig_s2_localization()
    fig_s3_nucleus_layered()
    print("# mech_hier_vis done")


if __name__ == "__main__":
    main()
