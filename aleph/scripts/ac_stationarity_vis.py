#!/usr/bin/env python
r"""Figures for a cortex-motor stationarity run — the lane's own viz entry point.

WHY THIS EXISTS.  The charter's visualization rule says every run writes its figures beside its record,
from the lane's own script, and commits them with the data.  `ac_gate_b_cortex_motor_native.py` wrote
none: the run that produced this project's first settled native measurement landed as a JSON trajectory
and nothing else, which is the shape the rule exists to forbid — a text artefact hides the physics, and
the visual check is what catches a wrong sign or an outlier that a band passes.

WHAT THE THREE FIGURES ARE FOR, and each is a question the number alone cannot answer:

1. **Trajectory + verdict.** The series with its discarded transient shaded, the analysed window marked,
   and the steady-state mean drawn with its ±sem band.  A verdict of STATIONARY is a claim about a
   PICTURE — that the trend is smaller than the fluctuation — and this is that picture.
2. **The source/network split.** Whether the crossbridge stress term or the loaded network carries the
   tension is the whole reading of a blebbistatin sweep, and the two differ by more than two decades, so
   they are drawn on separate axes rather than one plot where the smaller is a flat line.
3. **Autocorrelation.** τ_int is the number the error bar rests on, and it was measured AT the sampling
   floor — the ACF plot is where a reader sees that the series is already uncorrelated at one sample
   interval, i.e. that τ is an upper bound set by the protocol rather than an estimate of the physics.

Figure integrity follows the charter: no axis truncation, units on every axis, and the fluctuation band
drawn rather than implied.

Usage:
    python aleph/scripts/ac_stationarity_vis.py --record outputs/ac/stationarity/<run>/record.json

Sanity Gate:
    * dimensional: time axes are SECONDS of physical time, taken from the record's ``t_s`` column, never
      a step index — the whole point of the run is that the step count is not the measurement.
    * boundary: a record whose verdict is not STATIONARY has no mean and no sem; the figure then draws
      the series and the verdict, and omits the band rather than inventing one.
    * conservation/invariant: γ_total is drawn against γ_source + γ_network on the same figure, so a
      split that does not sum is visible rather than asserted.
    * numerical: the ACF is the same FFT estimator ``observe/stationarity.py`` uses, so the plot and the
      recorded τ cannot disagree.
    * sign-sense: the drift is annotated SIGNED, because a tension falling to its steady state and one
      rising to it are different physics.
    * measurement-protocol: the contract (min_windows, drift_sigma, driving τ) is printed on the figure,
      so a verdict is never seen apart from the criterion that produced it.

Runs on the dev machine: pure host-side plotting of an artifact a run already wrote.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

#: Observable → axis label, in the order the figures present them.
LABELS = {
    "gamma_total_pn_per_um": r"$\gamma_{\rm total}$  [pN/µm]",
    "gamma_network_pn_per_um": r"$\gamma_{\rm network}$  [pN/µm]",
    "gamma_source_pn_per_um": r"$\gamma_{\rm source}$ (crossbridge)  [pN/µm]",
    "bound_fraction": "bound head fraction  [-]",
}


def _series(record: dict, key: str) -> tuple[np.ndarray, np.ndarray]:
    """Return ``(t_s, values)`` for one observable."""
    rows = record["measurements"]["trajectory"]
    return (np.array([r["t_s"] for r in rows], dtype=float),
            np.array([r[key] for r in rows], dtype=float))


def _contract_text(record: dict, key: str) -> str:
    """Return the one-line pre-declared contract, for printing next to the verdict."""
    st = record["measurements"]["stationarity"][key]
    c = st.get("contract", {})
    return (f"contract: window ≥ {c.get('min_windows', '?')}× driving τ "
            f"({c.get('driving_correlation_time_s', float('nan')):.3g} s), "
            f"|drift| < {c.get('drift_sigma', '?')}σ")


def trajectory_figure(record: dict, out: Path) -> Path:
    """Draw γ_total with its transient, analysed window, mean and ±sem band."""
    key = "gamma_total_pn_per_um"
    t, y = _series(record, key)
    st = record["measurements"]["stationarity"][key]
    fig, ax = plt.subplots(figsize=(9.5, 4.6))

    ax.plot(t, y, lw=0.9, color="#1f4e79", label="per-sample measurement")
    t_eq = float(st.get("equilibration_time_s", 0.0))
    if t_eq > 0:
        ax.axvspan(t[0], t_eq, color="#bbbbbb", alpha=0.35, lw=0,
                   label=f"discarded transient ({t_eq:.2f} s, chosen by max $n_{{\\rm eff}}$)")
    mean, sem = st.get("mean"), st.get("sem")
    if mean is not None:
        ax.axhline(mean, color="#c00000", lw=1.4,
                   label=f"steady-state mean {mean:.5g} ± {sem:.2g} (sem)")
        ax.axhspan(mean - sem, mean + sem, color="#c00000", alpha=0.20, lw=0)
    ax.set_xlabel("physical time  [s]")
    ax.set_ylabel(LABELS[key])
    cfg = record["config"]
    ax.set_title(
        f"{st['verdict']} — {record['census']['cortex_filaments']:,} cortical filaments, "
        f"$k_{{xb}}$ = {cfg['k_xb']:g} pN/µm, {record['measurements']['lifetimes_covered']:.1f} "
        f"bound-head lifetimes\n{_contract_text(record, key)}   ·   "
        f"drift {st.get('drift_per_s', float('nan')):+.2g} /s "
        f"= {st.get('drift_over_std', float('nan')):.2f}σ over the window   ·   "
        f"$\\tau_{{\\rm int}}$ {st.get('tau_int_s', float('nan')):.3g} s, "
        f"$n_{{\\rm eff}}$ {st.get('n_eff', 0):.0f}",
        fontsize=8.5)
    ax.legend(fontsize=8, loc="lower right")
    ax.margins(x=0.01)
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    plt.close(fig)
    return out


def split_figure(record: dict, out: Path) -> Path:
    """Draw the source/network split on separate axes, plus the closure check."""
    t, total = _series(record, "gamma_total_pn_per_um")
    _, net = _series(record, "gamma_network_pn_per_um")
    _, src = _series(record, "gamma_source_pn_per_um")
    fig, axes = plt.subplots(3, 1, figsize=(9.5, 7.4), sharex=True)

    for ax, (key, y, colour) in zip(axes, (
            ("gamma_network_pn_per_um", net, "#1f4e79"),
            ("gamma_source_pn_per_um", src, "#c00000"),
            (None, total - (net + src), "#444444"))):
        ax.plot(t, y, lw=0.9, color=colour)
        ax.set_ylabel(LABELS[key] if key else r"$\gamma_{\rm total}-(\gamma_{\rm net}+\gamma_{\rm src})$"
                                              "\n[pN/µm]", fontsize=8.5)
        if key:
            st = record["measurements"]["stationarity"][key]
            ax.set_title(f"{key} — {st['verdict']}"
                         + (f", mean {st['mean']:.5g} ± {st['sem']:.2g}" if st.get("mean") is not None
                            else ", no mean returned"), fontsize=8.5, loc="left")
        else:
            ax.set_title("closure: the split must sum to the total", fontsize=8.5, loc="left")
    axes[-1].set_xlabel("physical time  [s]")
    frac = float(np.mean(src)) / float(np.mean(total)) if np.mean(total) else float("nan")
    fig.suptitle(f"crossbridge stress term is {100 * frac:.2f}% of the reported tension — "
                 f"myosin loads the NETWORK, and that is where its effect lives", fontsize=9)
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    plt.close(fig)
    return out


def acf_figure(record: dict, out: Path) -> Path:
    """Draw the normalised autocorrelation of the analysed window against the sampling floor."""
    key = "gamma_total_pn_per_um"
    t, y = _series(record, key)
    st = record["measurements"]["stationarity"][key]
    dt = float(record["measurements"]["sample_dt_s"])
    start = int(round(float(st.get("equilibration_time_s", 0.0)) / dt))
    tail = y[start:]
    x = tail - tail.mean()
    size = 1 << (2 * x.size - 1).bit_length()
    spec = np.fft.rfft(x, size)
    acf = np.fft.irfft(spec * np.conjugate(spec), size)[: x.size].real
    acf /= acf[0]
    lags = np.arange(acf.size) * dt

    fig, ax = plt.subplots(figsize=(9.5, 4.2))
    n_show = max(8, min(acf.size, int(60)))
    ax.axhline(0.0, color="#888888", lw=0.8)
    ax.plot(lags[:n_show], acf[:n_show], marker="o", ms=3, lw=1.0, color="#1f4e79")
    tau = float(st.get("tau_int_s", float("nan")))
    ax.axvline(tau, color="#c00000", lw=1.3, ls="--",
               label=fr"measured $\tau_{{\rm int}}$ = {tau:.4g} s")
    ax.axvline(dt / 2.0, color="#888800", lw=1.2, ls=":",
               label=fr"estimator floor $\Delta t/2$ = {dt / 2:.4g} s (uncorrelated series)")
    drive = float(record["measurements"]["bound_head_lifetime_s"])
    ax.set_xlabel(f"lag  [s]      (driving bound-head lifetime = {drive:.3g} s, off this axis)")
    ax.set_ylabel(r"normalised autocorrelation  $\rho_k$  [-]")
    ax.set_title(
        r"$\tau_{\rm int}$ sits AT the sampling floor: the series is already uncorrelated at one sample "
        "interval,\nso this is an UPPER BOUND set by the protocol, not an estimate of the physics — and "
        f"it is {drive / tau:.0f}× shorter than the process driving it", fontsize=8.5)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    plt.close(fig)
    return out


def main() -> None:
    """Render the three figures next to a run record."""
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--record", required=True, help="run-record@2 JSON written by the cortex-motor driver")
    ap.add_argument("--outdir", default=None, help="default: <record dir>/figs")
    args = ap.parse_args()

    path = Path(args.record)
    record = json.loads(path.read_text())
    outdir = Path(args.outdir) if args.outdir else path.parent / "figs"
    outdir.mkdir(parents=True, exist_ok=True)

    made = [
        trajectory_figure(record, outdir / "stationarity_trajectory.png"),
        split_figure(record, outdir / "stationarity_source_network_split.png"),
        acf_figure(record, outdir / "stationarity_autocorrelation.png"),
    ]
    for figure in made:
        print(f"[vis] {figure}")


if __name__ == "__main__":
    main()
