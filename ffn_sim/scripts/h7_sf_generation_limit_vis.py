"""B5(ii) — SF generation-limit FIGURE (waterfall + cortical-γ unification).

Visualises the SF tension generation-limit synthesis
(docs/v2_audit/SF_GENERATION_LIMIT_SYNTHESIS_2026-06-09.md): a labeled waterfall of
the raw single-SF force budget up to the Kumar 2006 band, and a side panel showing the
SAME ½·n·f·ℓ generation-limit structure for the cortical active-γ floor (the unification).

Integrity (CLAUDE.md): no axis truncation hidden (log axis NOTED on the y-label), the
Kumar reference band overlaid on the measurement, SI units annotated, every multiplier
labeled with its physical cause. Reads the live JSON the budget script emits.

Run:  python ffn_sim/scripts/h7_sf_generation_limit_vis.py
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

_HERE = Path(__file__).resolve()
if str(_HERE.parents[2]) not in sys.path:
    sys.path.insert(0, str(_HERE.parents[2]))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

_OUT = _HERE.parents[1] / "outputs" / "h7"
_JSON = _OUT / "production" / "h7_basal_sf_force_budget.json"


def _load() -> dict:
    if not _JSON.exists():
        # generate it if missing
        import subprocess
        subprocess.run([sys.executable, str(_HERE.parent / "h7_basal_sf_force_budget.py")],
                       check=True)
    return json.loads(_JSON.read_text())


def make_figure(rep: dict):
    b = rep["budget"]
    d = rep["decomposition"]
    raw_pN = b["T_cable_raw_N"] * 1e12
    fid = d["per_minifilament_fidelity_factor"]
    cnt_lo = d["cross_sectional_nmii_count_to_band_lo"]
    cnt_hi = d["cross_sectional_nmii_count_to_band_centre"]
    band_lo_pN = b["kumar_band_N"][0] * 1e12
    band_hi_pN = b["kumar_band_N"][1] * 1e12

    # waterfall stages (in pN): raw → ×fidelity → ×count(floor→centre spans the band)
    after_fid = raw_pN * fid
    after_cnt_lo = after_fid * cnt_lo
    after_cnt_hi = after_fid * cnt_hi

    fig, (axL, axR) = plt.subplots(1, 2, figsize=(13.5, 6.4))

    # ---- LEFT: SF waterfall ----
    stages = ["raw budget\n(brief-literal\nminifilament)",
              "× ~%.0f\nper-minifilament\nfidelity (→lit 56 pN;\n2.8×heads × 4×stall,\nparams)" % fid,
              "× ~%.0f–%.0f\ncross-sectional\nNMII count / SF\n(Route-B density datum)" % (cnt_lo, cnt_hi)]
    xs = [0, 1, 2]
    vals_lo = [raw_pN, after_fid, after_cnt_lo]
    vals_hi = [raw_pN, after_fid, after_cnt_hi]
    axL.plot(xs, vals_hi, "o-", color="#c0392b", lw=2, ms=9, label="to band centre (20 nN)")
    axL.plot(xs, vals_lo, "o--", color="#e67e22", lw=2, ms=7, label="to band floor (10 nN)")
    for x, v in zip(xs, vals_hi):
        axL.annotate(f"{v:.0f} pN" if v < 1e3 else f"{v/1e3:.1f} nN",
                     (x, v), textcoords="offset points", xytext=(8, 8), fontsize=9)
    # Kumar reference band
    axL.axhspan(band_lo_pN, band_hi_pN, color="#2ecc71", alpha=0.25,
                label="Kumar 2006 single-SF band\n(10–30 nN)")
    axL.set_yscale("log")
    axL.set_ylabel("single-SF tension  [pN]   (LOG axis)")
    axL.set_xticks(xs)
    axL.set_xticklabels(stages, fontsize=8.5)
    axL.set_title("SF tension is GENERATION-limited (N_filaments-INDEPENDENT)\n"
                  f"raw ≈ {raw_pN:.1f} pN → ~{b['gap_to_band_lo']:.0f}× under Kumar floor",
                  fontsize=10.5)
    axL.legend(loc="upper left", fontsize=8.5)
    axL.grid(True, which="both", alpha=0.25)
    axL.text(0.02, 0.02,
             "UPPER bound on a COHERENT bundle.\n§11 SF-2c: random mixed-polarity → −61 pN\n"
             "(slackening) → needs SARCOMERIC organisation.\n"
             f"Engagement realism (φ0.99→11%): ~{d['engagement_realism_penalty_factor']:.0f}× worse.",
             transform=axL.transAxes, fontsize=8, va="bottom", ha="left",
             bbox=dict(boxstyle="round", fc="#fdf6e3", ec="#b58900", alpha=0.9))

    # ---- RIGHT: unification table-as-text ----
    axR.axis("off")
    axR.set_title("ONE ½·n·f·ℓ generation-limit  —  cortex ∥ SF", fontsize=11)
    lines = [
        ("",                       "cortical active-γ",        "ventral SF"),
        ("budget",                 "½·n₂D·f·ℓ  (N/m)",          "N_cross·f  (N)"),
        ("dipole f",               "56 pN",                     "56 pN (lit ref)"),
        ("density",                "areal n₂D [µm⁻²]",          "cross-sect. N [/SF]"),
        ("target",                 "0.39–0.41 mN/m (Hosseini21)", "10–30 nN (Kumar06)"),
        ("gap @ proxy",            "~40–80× under",             "~180× under*"),
        ("needs",                  "n₂D≈16–47 µm⁻²",            "N_cross≈180–360"),
        ("MCF7 datum",             "NONE (HeLa proxy)",         "NONE (session i)"),
        ("right lever",            "actin OVERLAP",             "sarcomeric POLARITY"),
        ("",                       "(Chugh/Truong Quang)",      "(Hotulainen-Lapp.)"),
    ]
    y = 0.92
    for i, (a, c, s) in enumerate(lines):
        weight = "bold" if i == 0 else "normal"
        color = "#1f3a93" if i == 0 else "black"
        axR.text(0.01, y, a, fontsize=9, fontweight="bold", transform=axR.transAxes)
        axR.text(0.34, y, c, fontsize=8.7, fontweight=weight, color=color, transform=axR.transAxes)
        axR.text(0.70, y, s, fontsize=8.7, fontweight=weight, color=color, transform=axR.transAxes)
        y -= 0.092
    axR.text(0.01, 0.02,
             "* session (i) REFUTE: the cross-section count datum is\n"
             "MISSING + GEOMETRICALLY IMPOSSIBLE (band needs ~590–1760,\n"
             "a 50–250 nm cross-section holds O(5–15)). Reframe (Kassianidou/\n"
             "Kumar 2017): Kumar 10–30 nN is mostly network/prestress;\n"
             "single-fiber ACTIVE ≈6 nN — still ~7× under. Both floors:\n"
             "generation-bound, same missing motor-density datum.",
             transform=axR.transAxes, fontsize=7.3, va="bottom",
             bbox=dict(boxstyle="round", fc="#eef5ff", ec="#1f3a93", alpha=0.9))

    fig.suptitle("H.SF (ii) — SF tension generation-limit, unified with the cortical γ-floor "
                 "(2026-06-09)", fontsize=12, y=0.99)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    return fig


def main() -> int:
    rep = _load()
    fig = make_figure(rep)
    (_OUT / "figs").mkdir(parents=True, exist_ok=True)
    out = _OUT / "figs" / "h7_sf_generation_limit.png"
    fig.savefig(out, dpi=130, bbox_inches="tight")
    print(f"[SF generation-limit fig] {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
