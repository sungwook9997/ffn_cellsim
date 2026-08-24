"""S.1 driver — homogeneous (visco)elastic sphere under distributed (parallel-plate) load.

Stage 1 of the advisor's hierarchical mechanics flow
(``docs/v2_audit/_historical/HIERARCHICAL_MECHANICS_VALIDATION_PLAN_2026-07-14.md``,
contract ``docs/briefs/S1_homogeneous_sphere.md``).

The S1 object is a SINGLE-COMPARTMENT FF cortex shell (nucleus OFF, microtubules OFF,
no membrane-extras) = an effective homogeneous viscoelastic pressurized shell. This
driver builds it, runs a parallel-plate compression strain sweep, records
force-displacement, and gates ``F(delta)`` FORWARD against the committed Hertz oracle
(``validation/oracles/mechanics/hertz.py``) — the advisor's "match analytical/benchmark
+ consistent parameter response" gate, written before the run.

Engine units: pN, um, s. Because ``1 Pa == 1 pN/um^2`` exactly, the Hertz oracle fed
``E*`` in Pa and ``R, delta`` in um returns force directly in pN (no conversion).

Model notes (grounded in the FF-build audit + FF_AFM_REALDATA_VALIDATION_2026-07-02):
- Build FRESH and let the compression driver relax it. ``--from-resting`` is a
  substrate-ADHERED checkpoint (crawl), NOT valid for a FREE compressed sphere.
- Regulated turgor (``pressure_setpoint=40 Pa`` + ``K_drained=300 Pa``) is the validated
  liquid-drop closure that keeps the apparent tension confinement-independent.
- Small-strain window (delta/R <= 0.03, i.e. strain <= ~3%) is where the Hertz benchmark
  is valid; the soft plate is validated there. Larger strain / rigid-plate is for the
  full F-delta curve (see ``--max-strain`` / ``--rigid-plate``).
- CPU / small-NF runs are DEV SMOKE tests (labelled non-authoritative); the authoritative
  S1 result is the native GPU run (``--mode native``). Never draw an S1 conclusion from a
  smoke run (CLAUDE.md native+full HARD rule).
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np

from aleph.laws.gamma_floor import (
    CortexParams,
    NMIIA_MINIFIL_STALL_PN,
    build_crosslinked_cortex,
)
from aleph.laws.network_warp import simulate_whole_cell_compression_on_device
from aleph.laws.compartments import resolve_membrane
from aleph.validation.oracles.mechanics import hertz

NU = 0.5                          # incompressible cell (ff_hertz_validation convention)
E_MCF7_ZBIRAL = 249.0             # Pa, whole-cell colloidal bead (KB-6.1.1, Zbiral 2023)
E_MCF7_BAND = (224.0, 279.0)      # Zbiral breast panel band
TURGOR_DP0 = 40.0                 # Pa, resting osmotic pressure (physiological baseline)
K_DRAINED = 300.0                 # Pa, Moeendarbary drained bulk modulus
LP_CENTRAL = 1.0e-7               # um/(s*Pa), COS-7 + airway Pf water permeability
ETA_CYTO = 65.9                   # Pa*s, MCF7 cytoplasm bulk viscosity (Dessard 2024)
XL_KOFF = 0.066                   # 1/s, alpha-actinin off-rate (Ferrer 2008)
SMALL_STRAIN_MAX = 0.03           # Hertz validity: delta1/R0 <= 0.03

# Loading protocols passed to the compression driver. The default was NON-physiological
# (instant strain, no drainage/turnover/bulk-viscosity, soft plate) and reproduced the
# documented AFM overshoot (native ~21900x too stiff, project-ff-afm-overshoot-was-artifact).
# PHYSIO is the validated biphasic + turnover + bulk-viscosity + rigid-plate protocol that
# brought the native cortex to 0.73x MCF7 (IN BAND) — mandated by the physiological-baseline
# HARD rule (loading must be at real viscosity/rate, not instant/convenient).
LOAD_NAIVE = {"pressure_setpoint": TURGOR_DP0, "K_drained_Pa": K_DRAINED, "rigid_plate": False}
LOAD_PHYSIO = {
    "pressure_setpoint": None,       # None enables biphasic drainage in the driver
    "Lp_um_s_Pa": LP_CENTRAL,        # water efflux (Kedem-Katchalsky drainage)
    "K_drained_Pa": K_DRAINED,       # Terzaghi drained-solid bulk modulus
    "load_time_s": 300.0,            # long physical load -> fully drained equilibrium
    "xl_koff_per_s": XL_KOFF,        # alpha-actinin turnover -> emergent stress relaxation
    "eta_bulk_Pa_s": ETA_CYTO,       # real cytoplasm viscosity (not water)
    "rigid_plate": True,             # hard-clamp reaction (avoids soft-penalty NF-scaling)
}


def build_bare_cortex_shell(n_filaments: int, seed: int = 1, R_um: float = 7.5):
    """Build a single-compartment cortex shell (cortex + turgor only), free (unadhered).

    Args:
        n_filaments: Cortex filament count (``--cortex-fil``). Crosslinks = n_filaments,
            myosin = n_filaments // 10.
        seed: RNG seed for the random cortex network.
        R_um: Rest cell radius (um).

    Returns:
        A ``CrosslinkedCortex`` with ``R0_mean`` (self-consistent turgor volume reference)
        set to the mean node radius, ready for the compression driver.
    """
    cx = build_crosslinked_cortex(
        CortexParams(R_um=R_um),
        n_filaments=n_filaments,
        n_xl=n_filaments,
        n_myo=n_filaments // 10,
        rng=np.random.default_rng(seed),
    )
    cx.R0_mean = float(np.linalg.norm(cx.net.pos - cx.net.pos.mean(0), axis=1).mean())
    return cx


def _cortex_stress_p95(cx, pos_all, m) -> float:
    """95th-percentile per-node von-Mises stress [Pa] on the compressed cortex block.

    Best-effort: returns NaN if the cortex crosslink/myosin index arrays are not present
    on this build (so a smoke run still yields F-delta even if the stress carrier changes).
    """
    try:
        Nc = int(m["Nc"])
        pcx = np.asarray(pos_all)[:Nc]
        R0 = float(cx.R0_mean)
        V0 = (4.0 / 3.0) * np.pi * R0**3
        from aleph.laws.ff_virial_stress import cortex_node_stress

        svm, _ = cortex_node_stress(
            pcx,
            np.full(Nc, V0 / Nc),
            xl_ij=np.stack([cx.xl_i, cx.xl_j], 1).astype(np.int64),
            k_xl=cx.xl_k,
            r0_xl=cx.xl_rest,
            myo_ij=np.stack([cx.myo_i, cx.myo_j], 1).astype(np.int64),
            f_myo=float(NMIIA_MINIFIL_STALL_PN),
            dP=float(m["dP_turgor_Pa"]),
        )
        return float(np.percentile(svm, 95))
    except Exception as exc:  # noqa: BLE001 - stress is a best-effort auxiliary output
        print(f"# [warn] stress field skipped: {type(exc).__name__}: {exc}", flush=True)
        return float("nan")


def s1_compression_sweep(
    strains,
    *,
    n_filaments: int = 2000,
    n_steps: int = 1600,
    device: str = "cpu",
    load: dict | None = None,
    membrane=None,
    seed: int = 1,
    with_stress: bool = True,
    quiet: bool = False,
) -> dict:
    """Run a parallel-plate compression strain sweep on the bare cortex shell.

    One FRESH build + sim per strain point (the driver returns only the final-state
    plate force). Per-contact indentation ``delta1 = R0 * strain``; total plate approach
    ``D = 2 * delta1``.

    Args:
        strains: Iterable of nominal compressive strains (dimensionless).
        n_filaments: Cortex filament count.
        n_steps: Quasi-static relaxation step budget per point.
        device: ``"cpu"`` (dev smoke) or ``"cuda:0"`` (native GPU).
        load: Loading-protocol kwargs splatted into the compression driver
            (see ``LOAD_PHYSIO`` / ``LOAD_NAIVE``): pressure_setpoint, Lp_um_s_Pa,
            K_drained_Pa, load_time_s, xl_koff_per_s, eta_bulk_Pa_s, rigid_plate, ...
        membrane: Optional ``ResolvedMembrane`` (area reservoir) — a higher layer than the
            bare cortex; ``None`` = pure single-compartment S1 shell.
        seed: RNG seed.
        with_stress: Compute the per-node von-Mises stress p95 at each strain.
        quiet: Suppress per-point printing.

    Returns:
        Dict with ``R0_um``, ``rows`` (list of per-strain dicts), and the config echo.
    """
    load = dict(load if load is not None else LOAD_NAIVE)
    R0 = build_bare_cortex_shell(n_filaments, seed=seed).R0_mean
    if not quiet:
        print(
            f"# S1 sweep | N_fil={n_filaments} R0={R0:.3f}um device={device} "
            f"membrane={'on' if membrane is not None else 'off'} load={load}",
            flush=True,
        )
        print("strain  delta1[um]  F[pN]     dP[Pa]  V/V0   svm95[Pa]", flush=True)
    rows = []
    for s in strains:
        cx = build_bare_cortex_shell(n_filaments, seed=seed)
        pos_all, m = simulate_whole_cell_compression_on_device(
            cx,
            NMIIA_MINIFIL_STALL_PN,
            strain=float(s),
            nucleus=None,
            membrane=membrane,
            microtubule=None,
            n_steps=n_steps,
            turgor_every=(50 if device.startswith("cuda") else 20),
            device=device,
            **load,
        )
        delta1 = R0 * float(s)
        F = float(m["F_plate_pN"])
        svm95 = _cortex_stress_p95(cx, pos_all, m) if with_stress else float("nan")
        rows.append(
            {
                "strain": float(s),
                "delta1_um": delta1,
                "delta_total_um": 2.0 * delta1,
                "F_pN": F,
                "dP_Pa": float(m["dP_turgor_Pa"]),
                "V_over_V0": float(m.get("V_over_V0", float("nan"))),
                "gamma_mN_m": float(m.get("gamma_apparent_mN_m", float("nan"))),
                "svm_p95_Pa": svm95,
                "delta1_over_R": delta1 / R0,
            }
        )
        if not quiet:
            print(
                f"{s*100:5.1f}%  {delta1:8.3f}  {F:9.1f}  {m['dP_turgor_Pa']:6.0f} "
                f"{m.get('V_over_V0', float('nan')):5.3f}  {svm95:9.1f}",
                flush=True,
            )
    return {
        "R0_um": R0,
        "rows": rows,
        "config": {
            "n_filaments": n_filaments,
            "n_steps": n_steps,
            "device": device,
            "membrane": membrane is not None,
            "load": load,
            "seed": seed,
        },
    }


def gate_s1_hertz(
    sweep: dict, *, nu: float = NU, band=E_MCF7_BAND, small_strain_max: float = SMALL_STRAIN_MAX
) -> dict:
    """Gate the swept ``F(delta)`` FORWARD against the committed Hertz oracle.

    Fits the reduced modulus from the small-strain points (``delta1/R0 <= small_strain_max``)
    with :func:`hertz.fit_reduced_modulus`, converts to an apparent Young's modulus, and
    reports the fit quality, the pointwise-E flatness (constitutive linearity), the
    small-strain validity, and whether ``E_fit`` lands in the measured MCF7 band.

    Args:
        sweep: Output of :func:`s1_compression_sweep`.
        nu: Poisson ratio (0.5 incompressible).
        band: Acceptance band for the apparent Young's modulus [Pa].
        small_strain_max: Hertz small-strain validity threshold on ``delta1/R0``.

    Returns:
        Dict with ``E_star_Pa``, ``E_fit_Pa``, ``r_squared``, ``flatness``,
        ``n_points_used``, ``all_small_strain_valid``, ``forward_max_resid_frac``,
        and ``verdict`` ("IN BAND" / "TOO STIFF" / "TOO SOFT" / "NO VALID POINTS").
    """
    R0 = sweep["R0_um"]
    rows = [r for r in sweep["rows"] if r["strain"] > 0.0]
    small = [r for r in rows if r["delta1_over_R"] <= small_strain_max and r["F_pN"] > 0.0]
    if len(small) < 1:
        return {"verdict": "NO VALID POINTS", "n_points_used": 0}

    d1 = np.array([r["delta1_um"] for r in small])
    F = np.array([r["F_pN"] for r in small])
    fit = hertz.fit_reduced_modulus(d1, F, R0)
    e_star = float(fit["E_star"])
    e_fit = hertz.youngs_from_reduced(e_star, nu)

    # Forward check: does the oracle prediction reproduce the measured points?
    F_pred = hertz.force_sphere_plane(e_star, R0, d1)
    resid_frac = float(np.max(np.abs(F_pred - F) / np.maximum(np.abs(F), 1e-12)))

    valid = hertz.valid_small_strain(d1, R0, small_strain_max)
    r2 = float(fit["r_squared"])
    flat = float(fit["E_star_flatness"])

    # S1 ANALYTICAL gate: does the bare body follow a clean CONSTITUTIVE law? Fit BOTH
    # the solid-Hertz law (F ~ delta^1.5) and the pressurized-thin-shell / Reissner-tension
    # law (F ~ delta, linear) — the better fit identifies the constitutive CLASS. A
    # turgor-pressurized cortex, once drained, indents ~linearly (shell); an undrained
    # water-filled body responds ~Hertzian (solid). This is NOT "does E_fit match the
    # MCF7 WHOLE-cell band" (a sub-component is not the whole cell; the MCF7 ratio is an
    # informational overlay only, and S1 forbids a biology/measured-modulus claim).
    xs = d1
    k_shell = float(np.dot(xs, F) / np.dot(xs, xs))            # linear stiffness [pN/um]
    shell_pred = k_shell * xs
    ss_tot = float(np.sum((F - F.mean()) ** 2))
    shell_r2 = 1.0 - float(np.sum((F - shell_pred) ** 2)) / ss_tot if ss_tot > 0 else 1.0
    shell_resid = float(np.max(np.abs(shell_pred - F) / np.maximum(np.abs(F), 1e-12)))

    hertz_clean = r2 >= 0.99 and flat <= 0.10 and resid_frac <= 0.10
    shell_clean = shell_r2 >= 0.99 and shell_resid <= 0.10
    better = "hertz-solid" if r2 >= shell_r2 else "linear-shell"
    clean = (hertz_clean if better == "hertz-solid" else shell_clean) and bool(valid["all_valid"])
    verdict = f"{better.upper()}" + (" (clean)" if clean else " (needs rigor)")

    return {
        "E_star_Pa": e_star,
        "E_fit_Pa": e_fit,
        "r_squared": r2,                       # solid-Hertz fit
        "flatness": flat,
        "forward_max_resid_frac": resid_frac,
        "k_shell_pN_per_um": k_shell,          # linear pressurized-shell fit
        "shell_r_squared": shell_r2,
        "shell_resid_frac": shell_resid,
        "better_model": better,
        "model_clean": clean,
        "verdict": verdict,
        "n_points_used": len(small),
        "all_small_strain_valid": bool(valid["all_valid"]),
        # informational overlay only (NOT the S1 pass/fail):
        "E_fit_over_MCF7": e_fit / E_MCF7_ZBIRAL,
        "mcf7_band_Pa": list(band),
    }


# --------------------------------------------------------------------------- #
# Output + best-effort figure
# --------------------------------------------------------------------------- #
def _write_outputs(outdir: Path, sweep: dict, gate: dict) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    (outdir / "manifest.json").write_text(
        json.dumps({"config": sweep["config"], "R0_um": sweep["R0_um"], "gate": gate}, indent=2)
    )
    rows = sweep["rows"]
    np.savez(
        outdir / "sweep.npz",
        strain=np.array([r["strain"] for r in rows]),
        delta1_um=np.array([r["delta1_um"] for r in rows]),
        F_pN=np.array([r["F_pN"] for r in rows]),
        dP_Pa=np.array([r["dP_Pa"] for r in rows]),
        svm_p95_Pa=np.array([r["svm_p95_Pa"] for r in rows]),
    )
    _write_report(outdir, sweep, gate)
    _maybe_plot(outdir, sweep, gate)


def _write_report(outdir: Path, sweep: dict, gate: dict) -> None:
    cfg = sweep["config"]
    auth = "NATIVE (authoritative)" if cfg["device"].startswith("cuda") and cfg["n_filaments"] >= 38000 else (
        "DEV SMOKE (non-authoritative)"
    )
    lines = [
        "# S.1 — homogeneous sphere, distributed load — REPORT",
        "",
        f"**Run class:** {auth}  ·  N_fil={cfg['n_filaments']}, device={cfg['device']}, "
        f"steps={cfg['n_steps']}, membrane={cfg['membrane']}",
        f"**Loading:** {cfg['load']}",
        f"**R0:** {sweep['R0_um']:.3f} um",
        "",
        "## S1 analytical gate (constitutive self-consistency: Hertz-solid vs pressurized-shell)",
        "",
        f"- **{gate.get('verdict')}** — better-fitting law = {gate.get('better_model')}",
        f"- solid-Hertz (F~δ^1.5): r^2 = {gate.get('r_squared', float('nan')):.4f}, "
        f"pointwise-E flatness = {gate.get('flatness', float('nan')):.3f}, "
        f"residual = {gate.get('forward_max_resid_frac', float('nan')):.1%}",
        f"- linear-shell (F~δ, Reissner/tension): r^2 = {gate.get('shell_r_squared', float('nan')):.4f}, "
        f"residual = {gate.get('shell_resid_frac', float('nan')):.1%}, "
        f"k = {gate.get('k_shell_pN_per_um', float('nan')):.0f} pN/µm",
        f"- Hertz E_fit = {gate.get('E_fit_Pa', float('nan')):.1f} Pa "
        f"(E* = {gate.get('E_star_Pa', float('nan')):.1f} Pa), small-strain valid = "
        f"{gate.get('all_small_strain_valid')}, n_points = {gate.get('n_points_used')}",
        f"- **INFORMATIONAL (not S1 pass/fail):** E_fit x{gate.get('E_fit_over_MCF7', float('nan')):.1f} "
        f"the MCF7 WHOLE-cell band {gate.get('mcf7_band_Pa')} Pa (bare cortex is a sub-component; MCF7 is a "
        f"downstream layer-assembled target).",
        "",
        "## Force–displacement",
        "",
        "| strain | delta1 [um] | F [pN] | dP [Pa] | V/V0 | svm_p95 [Pa] |",
        "|---|---|---|---|---|---|",
    ]
    for r in sweep["rows"]:
        lines.append(
            f"| {r['strain']*100:.1f}% | {r['delta1_um']:.3f} | {r['F_pN']:.1f} | "
            f"{r['dP_Pa']:.0f} | {r['V_over_V0']:.3f} | {r['svm_p95_Pa']:.1f} |"
        )
    lines += [
        "",
        "## Figures",
        "",
        "- `figs/s1_force_displacement.png` — F vs delta1 with the fitted Hertz oracle overlaid.",
        "",
    ]
    (outdir / "REPORT.md").write_text("\n".join(lines))


def _maybe_plot(outdir: Path, sweep: dict, gate: dict) -> None:
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        rows = [r for r in sweep["rows"] if r["strain"] > 0.0]
        d1 = np.array([r["delta1_um"] for r in rows])
        F = np.array([r["F_pN"] for r in rows])
        figs = outdir / "figs"
        figs.mkdir(exist_ok=True)
        fig, ax = plt.subplots(figsize=(5.2, 4.0))
        ax.plot(d1, F, "o", ms=6, label="FF cortex shell (measured)")
        if "E_star_Pa" in gate:
            dd = np.linspace(d1.min(), d1.max(), 100)
            ax.plot(
                dd,
                hertz.force_sphere_plane(gate["E_star_Pa"], sweep["R0_um"], dd),
                "-",
                lw=1.6,
                label=f"Hertz oracle (E_fit={gate['E_fit_Pa']:.0f} Pa)",
            )
        ax.set_xlabel(r"per-contact indentation $\delta_1$ [µm]")
        ax.set_ylabel("plate force $F$ [pN]")
        ax.set_title(f"S1 force–displacement · {gate.get('verdict','')}")
        ax.legend(frameon=False, fontsize=8)
        fig.tight_layout()
        fig.savefig(figs / "s1_force_displacement.png", dpi=140)
        plt.close(fig)
    except Exception as exc:  # noqa: BLE001 - viz is best-effort
        print(f"# [warn] figure skipped: {type(exc).__name__}: {exc}", flush=True)


def main() -> None:
    ap = argparse.ArgumentParser(description="S.1 homogeneous-sphere distributed-load driver (FF).")
    ap.add_argument("--mode", choices=["smoke", "native"], default="smoke",
                    help="smoke=CPU dev (non-authoritative); native=GPU authoritative.")
    ap.add_argument("--nfil", type=int, default=None, help="cortex filament count (override).")
    ap.add_argument("--steps", type=int, default=None)
    ap.add_argument("--device", default=None)
    ap.add_argument("--max-strain", type=float, default=0.03,
                    help="max strain in the sweep (0.03 = Hertz small-strain window).")
    ap.add_argument("--min-strain", type=float, default=0.005,
                    help="min (nonzero) strain; Hertz needs delta>0, F(0) is a pre-contact baseline "
                         "(matches the validated ff_hertz_validation window).")
    ap.add_argument("--n-strain", type=int, default=6)
    ap.add_argument("--physio", action="store_true",
                    help="validated physiological loading (biphasic drainage + turnover + bulk-eta "
                         "+ rigid plate); default is the naive instant/soft protocol.")
    ap.add_argument("--f-excess", type=float, default=0.0,
                    help="membrane area reservoir (0 = bare cortex S1; >0 adds the membrane layer).")
    ap.add_argument("--no-stress", action="store_true")
    ap.add_argument("--out", default="aleph/outputs/mech_hier/s1_sphere")
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()

    if a.mode == "native":
        nfil = a.nfil if a.nfil is not None else 70686   # density-anchored native (AFM-validation NF)
        steps = a.steps if a.steps is not None else 4000
        device = a.device if a.device is not None else "cuda:0"
    else:
        nfil = a.nfil if a.nfil is not None else 1500
        steps = a.steps if a.steps is not None else 1200
        device = a.device if a.device is not None else "cpu"

    load = dict(LOAD_PHYSIO if a.physio else LOAD_NAIVE)
    membrane = resolve_membrane(f_excess=a.f_excess) if a.f_excess > 0.0 else None
    strains = np.linspace(a.min_strain, a.max_strain, a.n_strain)
    t0 = time.time()
    sweep = s1_compression_sweep(
        strains,
        n_filaments=nfil,
        n_steps=steps,
        device=device,
        load=load,
        membrane=membrane,
        seed=a.seed,
        with_stress=not a.no_stress,
    )
    gate = gate_s1_hertz(sweep)
    print(
        f"# GATE {gate.get('verdict')} · E_fit={gate.get('E_fit_Pa', float('nan')):.1f} Pa "
        f"(x{gate.get('E_fit_over_MCF7', float('nan')):.1f} MCF7) | "
        f"Hertz r2={gate.get('r_squared', float('nan')):.3f} shell r2={gate.get('shell_r_squared', float('nan')):.3f} "
        f"({time.time()-t0:.0f}s)",
        flush=True,
    )
    _write_outputs(Path(a.out), sweep, gate)
    print(f"# wrote {a.out}/REPORT.md", flush=True)


if __name__ == "__main__":
    main()
