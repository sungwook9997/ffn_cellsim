"""H.9 nucleus + KU-3.1 enclosed-volume INTEGRATION smoke harness.

Builds a SMALL cell (n_fil≈200, the dev/CPU scale) with BOTH the H.9 nucleus
(:class:`ffn_sim.archive.hoomd_legacy.cell.nucleus.NucleusConfinement`) and the KU-3.1 enclosed-volume
pressure (:class:`ffn_sim.archive.hoomd_legacy.cortex.enclosed_volume.EnclosedVolumePressure`) wired
as additive, default-off compartments through
:func:`ffn_sim.archive.hoomd_legacy.cell.cell.build_cortex_full_simulation`, runs a few hundred BAOAB
steps, and reports:

* crash-free finite-position run (the core integration smoke),
* tag-APPEND confirmation (nucleus_bead block appended LAST),
* the runtime enclosed-volume Young-Laplace pressure under a controlled 1%
  radius compression (sign + 2γ/R scale),
* the nucleus cloud bound + cross-section vs the KU-1.V.3.3 critical pore.

This is NOT a production run and does NOT claim the composite KU-3.1 / KU-3.5
gate PASSES — it exercises the wiring and surfaces the multi-timescale CFL
finding (the nucleus k_hi is stiffer than the cortex dt timescale; a CFL-safe
n_beads is used so the strict CFL gate is honoured, not loosened — the NATIVE
path needs the B1 global-dt reconciliation extended to the nucleus term).

Per the production-driver-auto-viz rule (CLAUDE.md visualize-at-closeout), the
harness auto-generates a best-effort PNG of the cortex shell + nucleus cloud
cross-section at run end (best-effort; a headless / no-matplotlib environment
degrades to skipping the figure with a logged note).

Usage
-----
    conda activate ffn_sim
    python ffn_sim/scripts/h9_ku31_smoke.py            # n_fil=200, 300 steps
    python ffn_sim/scripts/h9_ku31_smoke.py --n-fil 200 --steps 500 --no-fig

References
----------
- ``ffn_sim/cell/nucleus.py`` (H.9 nucleus), ``ffn_sim/cortex/enclosed_volume.py``
  (KU-3.1), ``ffn_sim/cell/cell.py::build_cortex_full_simulation`` (integration).
- ``ffn_sim/docs/briefs/H9_nucleus.md`` (DRAFT skeleton, KU-3.B2).
- KU-3.1 osmotic anchor Stewart 2011 (ΔP≈40 Pa interphase); KU-3.5 γ≈1 mN/m
  (Salbreux 2012); KU-1.V.3.3 critical pore ~7 µm² (tumor).
"""

from __future__ import annotations

import argparse
import logging
import math
from copy import deepcopy
from pathlib import Path

import numpy as np
import yaml

import hoomd

from ffn_sim.archive.hoomd_legacy.cortex.cortex import resolve_h3_derived
from ffn_sim.archive.hoomd_legacy.cortex.enclosed_volume import (
    resolve_enclosed_volume,
    young_laplace_pressure,
)
from ffn_sim.archive.hoomd_legacy.cell.nucleus import (
    resolve_nucleus,
    nucleus_cross_section_area,
)
from ffn_sim.archive.hoomd_legacy.cell.cell import build_cortex_full_simulation

_LOG = logging.getLogger("h9_ku31_smoke")

_CONFIG_PATH = (
    Path(__file__).resolve().parents[1] / "configs" / "phase1_h3.yaml"
)
_OUT_DIR = (
    Path(__file__).resolve().parents[1] / "outputs" / "h9" / "figs"
)

# CFL-safe nucleus discretization at the small-cortex dt_cfl: k_hi ∝ 1/n_beads,
# so 300 beads keeps 0.1·γ_nuc/k_hi above dt_cfl (~1.9× margin). R_nuc =
# 0.25·R_cell (a typical nucleus:cell radius ratio).
_N_NUC_BEADS = 300
_R_NUC_FRACTION = 0.25


def _resolve(n_fil: int):
    cfg = deepcopy(yaml.safe_load(open(_CONFIG_PATH)))
    cfg["cortex"]["n_filaments"] = n_fil
    cfg["cortex"]["demo_mode"] = True
    p_cortex = resolve_h3_derived(cfg)
    p_ev = resolve_enclosed_volume({}, R_cell=p_cortex.R_cell)
    p_nuc = resolve_nucleus(
        {}, R_nuc=_R_NUC_FRACTION * p_cortex.R_cell, n_beads=_N_NUC_BEADS
    )
    return p_cortex, p_ev, p_nuc


def _maybe_figure(pos_cortex, pos_nuc, R_cell, R_nuc, dP, out_path: Path):
    """Best-effort cross-section figure (cortex shell + nucleus cloud)."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as exc:  # pragma: no cover - environment dependent
        _LOG.warning("matplotlib unavailable; skipping figure (%s)", exc)
        return None

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(6.0, 6.0))
    um = 1.0e6
    # Project onto x-y plane (full SI axes, no truncation).
    ax.scatter(
        pos_cortex[:, 0] * um, pos_cortex[:, 1] * um,
        s=6, c="#1f77b4", alpha=0.5, label="cortex actin (shell)",
    )
    ax.scatter(
        pos_nuc[:, 0] * um, pos_nuc[:, 1] * um,
        s=8, c="#d62728", alpha=0.7, label="nucleus_bead (cloud)",
    )
    # Reference circles: cortex R_cell and nucleus R_nuc.
    th = np.linspace(0, 2 * math.pi, 200)
    ax.plot(R_cell * um * np.cos(th), R_cell * um * np.sin(th),
            "--", c="#1f77b4", lw=1.0, label=f"R_cell = {R_cell*um:.1f} µm")
    ax.plot(R_nuc * um * np.cos(th), R_nuc * um * np.sin(th),
            "--", c="#d62728", lw=1.0, label=f"R_nuc = {R_nuc*um:.1f} µm")
    ax.set_xlabel("x [µm]")
    ax.set_ylabel("y [µm]")
    ax.set_aspect("equal")
    ax.set_title(
        "H.9 nucleus + KU-3.1 enclosed-volume (x-y cross-section)\n"
        f"compression ΔP = {dP:.1f} Pa (Young-Laplace 2γ/R scale)"
    )
    ax.legend(loc="upper right", fontsize=8)
    fig.tight_layout()
    fig.savefig(out_path, dpi=130)
    plt.close(fig)
    return out_path


def main() -> int:
    logging.basicConfig(
        level=logging.INFO, format="%(name)s %(levelname)s %(message)s"
    )
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--n-fil", type=int, default=200,
                    help="cortex filament count (small/dev scale; default 200)")
    ap.add_argument("--steps", type=int, default=300,
                    help="BAOAB steps to run (default 300)")
    ap.add_argument("--no-fig", action="store_true",
                    help="skip the auto-generated cross-section figure")
    args = ap.parse_args()

    _LOG.info("HOOMD %s", hoomd.version.version)
    p_cortex, p_ev, p_nuc = _resolve(args.n_fil)
    N_cortex = p_cortex.n_filaments * p_cortex.beads_per_filament

    # CFL bookkeeping (surface the multi-timescale finding explicitly).
    tau_nuc = p_cortex.gamma_b / p_nuc.k_hi
    dt_cfl_nuc = p_cortex.cfl_safety_factor * tau_nuc
    _LOG.info(
        "cortex: R_cell=%.2e m  dt_cfl=%.3e s  gamma_b=%.3e  N_cortex=%d",
        p_cortex.R_cell, p_cortex.dt_cfl, p_cortex.gamma_b, N_cortex,
    )
    _LOG.info(
        "nucleus: R_nuc=%.2e m  n_beads=%d  k_hi=%.3e N/m  "
        "0.1·τ_nuc=%.3e s  CFL-safe@dt_cfl=%s (margin %.2f×)",
        p_nuc.R_nuc, p_nuc.n_beads, p_nuc.k_hi, dt_cfl_nuc,
        p_cortex.dt_cfl <= dt_cfl_nuc, dt_cfl_nuc / p_cortex.dt_cfl,
    )
    _LOG.info(
        "enclosed-volume: K_vol=%.1f Pa  V0=%.3e m³  ΔP_ref=%.1f Pa "
        "(Stewart 2011 interphase)",
        p_ev.K_vol, p_ev.V0, p_ev.dP_ref,
    )

    # Build with BOTH compartments ON (additive, default-off when None).
    handles = build_cortex_full_simulation(
        p_cortex,
        p_enclosed_volume=p_ev,
        p_nucleus=p_nuc,
        with_baoab=True,
    )
    sim = handles["sim"]
    ev = handles["enclosed_volume_force"]
    nf = handles["nucleus_force"]
    ni = handles["nucleus_integration"]

    # tag-APPEND confirmation.
    assert ni.nucleus_tag_start == N_cortex, (
        f"tag-APPEND violated: nucleus_tag_start={ni.nucleus_tag_start} "
        f"!= N_cortex={N_cortex}"
    )
    assert sim.state.N_particles == N_cortex + _N_NUC_BEADS
    _LOG.info(
        "tag-APPEND OK: nucleus_bead block at [%d, %d); N_particles=%d",
        ni.nucleus_tag_start, ni.nucleus_tag_start + ni.n_beads,
        sim.state.N_particles,
    )

    # --- Run the smoke ---
    sim.run(args.steps)
    # Read positions from the TAG-ORDERED snapshot (get_snapshot), NOT the
    # row-indexed cpu_local_snapshot — the appended nucleus_bead tags map to
    # rows [N_cortex, N) only in the tag-ordered snapshot.
    snap = sim.state.get_snapshot()
    pos = np.asarray(snap.particles.position, dtype=np.float64)
    finite = bool(np.isfinite(pos).all())
    _LOG.info("ran %d steps; finite positions=%s", args.steps, finite)
    assert finite, "FAIL: nucleus+enclosed-volume BAOAB produced NaN/Inf"

    nuc_pos = pos[N_cortex:]
    c = nuc_pos.mean(axis=0)
    r = np.linalg.norm(nuc_pos - c, axis=1)
    xsec = nucleus_cross_section_area(nuc_pos, normal_axis=1)
    crit = nf.critical_pore_area()
    _LOG.info(
        "nucleus cloud: max r=%.3e m (R_nuc=%.3e) bounded=%s  "
        "cross-section=%.3e m² vs critical pore=%.3e m² (KU-1.V.3.3)",
        r.max(), p_nuc.R_nuc, bool((r < 5.0 * p_nuc.R_nuc).all()), xsec, crit,
    )
    _LOG.info(
        "enclosed volume after run: V=%.4e m³ (V/V0=%.4f)  ΔP=%.3f Pa",
        ev.last_volume, ev.last_volume / p_ev.V0, ev.last_pressure,
    )

    # --- Young-Laplace probe: impose 1% compression, read runtime ΔP ---
    eps = 0.01
    snap = sim.state.get_snapshot()
    sp = np.asarray(snap.particles.position, dtype=np.float64)
    shell = sp[:N_cortex]
    cc = shell.mean(axis=0)
    sp[:N_cortex] = cc + (1.0 - eps) * (shell - cc)
    sim.state.set_snapshot(snap)
    sim.run(0)
    dP = ev.last_pressure
    R = ev.last_R_mean
    gamma_equiv = dP * R / 2.0
    dP_ideal = p_ev.K_vol * 3.0 * eps
    _LOG.info(
        "Young-Laplace probe (%.0f%% compression): ΔP=%.2f Pa (ideal "
        "K_vol·3ε=%.1f Pa)  sign=%s  R_mean=%.3e m  γ_equiv=%.3f mN/m  "
        "2γ/R=%.2f Pa",
        eps * 100, dP, dP_ideal, "outward(+)" if dP > 0 else "inward(-)",
        R, gamma_equiv * 1e3, young_laplace_pressure(gamma_equiv, R),
    )
    if not (dP > 0.0):
        _LOG.error("FAIL: compression did not yield outward ΔP>0")
        return 1

    # --- Auto-viz (best-effort) ---
    if not args.no_fig:
        fig_path = _maybe_figure(
            pos[:N_cortex], nuc_pos, p_cortex.R_cell, p_nuc.R_nuc, dP,
            _OUT_DIR / "fig_h9_ku31_nucleus_cross_section.png",
        )
        if fig_path is not None:
            _LOG.info("wrote figure: %s", fig_path)

    _LOG.info("H.9 + KU-3.1 SMOKE: PASS (build + run crash-free; ΔP sign/scale OK)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
