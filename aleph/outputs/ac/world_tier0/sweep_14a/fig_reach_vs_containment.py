#!/usr/bin/env python
"""The two constraints on the NMII shell radius, on one axis — PI queue item 14.

Item 14 exists because reach and containment were never plotted against the same variable: the
membrane assert lives in ``world/geometry.py`` and the shell excursion guard in
``world/build/nmii.py``, each correct in its own scope. This renders the declared radius ladder
against BOTH, so the overlap is read rather than argued.

⚠ Neither axis is a magnitude claim. The separation is a geometric property of the BUILT cell (no
force evaluated); the 0.210 µm capture radius drawn as a reference band is the ``hand_kmc`` NMIIA
proxy and is a PI-GAP — it is drawn so the curve can be READ against it, never as a threshold this
run passed. ``n_heads_per_side=30`` is UNRATIFIED (PI queue 2) and ``support=706.86 µm²`` is PI
queue 19's wrong area, held FIXED across every point so the ladder has one cause.

Run:  python fig_reach_vs_containment.py
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = Path(__file__).parent
CAPTURE_PROXY_UM = 0.210        # hand_kmc NMIIA capture radius — PI-GAP, drawn as reference only
REACH_POINT_UM = 0.20           # the ladder rung nearest the proxy; a declared axis point

rows = []
for f in sorted(HERE.glob("r_*.json")):
    d = json.loads(f.read_text())
    n = d["rows"][0]["nearest_cortex_node_um"]
    st = {round(r["reach_um"], 3): r.get("n_stationed") for r in d["rows"]}
    rows.append((d["shells_um"]["nmii_radius_used"], n["min"], n["median"], n["p95"],
                 d["containment"]["max_nmii_node_radius_um"],
                 d["containment"]["n_nmii_nodes_outside_membrane"],
                 d["containment"]["guard_predicted_shell_excursion_um"],
                 st.get(REACH_POINT_UM), d["nmii"]["n_minifilaments"],
                 d["shells_um"]["r_cell"], d["shells_um"]["nmii_span"][1]))
rows.sort()
r, lo, med, hi, maxr, nout, guard, nst, nmf, r_cell, span_hi = (
    list(c) for c in zip(*rows, strict=True))

fig, ax = plt.subplots(1, 3, figsize=(15.5, 4.6))

ax[0].axhspan(0.0, CAPTURE_PROXY_UM, color="tab:green", alpha=0.12,
              label=f"within capture proxy {CAPTURE_PROXY_UM} µm (PI-GAP)")
ax[0].fill_between(r, lo, hi, alpha=0.25, color="tab:blue", label="p05–p95 across minifilaments")
ax[0].plot(r, med, "o-", color="tab:blue", label="median")
ax[0].plot(r, lo, ".--", color="tab:blue", lw=1, label="min")
ax[0].set_ylabel("distance to nearest cortex node  [µm]")
ax[0].set_title("REACH — can a motor touch the actin?")
ax[0].set_ylim(0.0, max(hi) * 1.08)

# ⚠ THE MIDDLE PANEL IS IN NANOMETRES, ON PURPOSE, AND THAT IS THE SECOND REDRAW. Draft 1 plotted
# the measured and nominal outer edges as absolute radii: over a 0.45 µm span a 3.3 nm difference is
# one line width, so the panel read "they agree" — the opposite of its finding. Draft 2 plotted the
# margin to the membrane, and the same 3.3 nm vanished again under a 0.45 µm range. Both were found
# by LOOKING at the render. The overshoot is the finding, so the overshoot gets the axis; where the
# build leaves the membrane is a SHADED REGION on x instead of a second squashed curve.
over_nm = [(m - sh) * 1e3 for m, sh in zip(maxr, span_hi, strict=True)]
margin_nm = [(rc - m) * 1e3 for rc, m in zip(r_cell, maxr, strict=True)]
bad = [x for x, k in zip(r, nout, strict=True) if k]
if bad:
    ax[1].axvspan(min(bad) - 0.025, max(r) + 0.02, color="tab:red", alpha=0.10,
                  label="assert_inside_membrane REFUSES this build")
# A CURVE, not a line: the excursion is offset^2/2r + L^2/8r, so it depends on the radius it is
# drawn at. Plotted flat it would assert a constant the formula does not have.
ax[1].plot(r, [g * 1e3 for g in guard], ":", color="tab:green", lw=1.4,
           label="excursion build/nmii.py predicts (its own formula)")
ax[1].plot(r, over_nm, "o-", color="tab:orange",
           label="MEASURED overshoot past the shell edge (r + t/2)")
for x, y, k, mg in zip(r, over_nm, nout, margin_nm, strict=True):
    if k:
        ax[1].annotate(f"{k}/32,708 nodes OUT\n{-mg:.1f} nm past the membrane", (x, y),
                       textcoords="offset points", xytext=(-160, 40), fontsize=8, color="tab:red")
ax[1].axhline(0.0, color="0.3", lw=1)
ax[1].set_ylim(min(0.0, min(over_nm) * 1.6), max(max(over_nm), max(guard) * 1e3) * 1.45)
ax[1].set_ylabel("radial overshoot past the nominal shell  [nm]")
ax[1].set_title("CONTAINMENT — AFTER the inset: the guard now promises something true")
ax[2].axhline(nmf[0], color="0.5", lw=1, ls=":", label=f"all {nmf[0]} minifilaments")
ax[2].plot(r, nst, "o-", color="tab:purple")
ax[2].set_ylabel(f"minifilaments with a legal station @ reach {REACH_POINT_UM} µm")
ax[2].set_title("STATIONS — does the population pair?")
ax[2].set_ylim(0, nmf[0] * 1.08)

for a in ax:
    a.set_xlabel("NMII shell radius passed to build_nmii  [µm]")
    a.grid(alpha=0.3)
    a.legend(fontsize=8, loc="best")
fig.suptitle("PI queue 14 (a) — the draw is INSET by the excursion. Overshoot is now NEGATIVE at every "
             "radius and 7.40 contains. "
             "Full native, 442 minifilaments, H=30 (UNRATIFIED). No magnitude claimed.", fontsize=10)
fig.tight_layout()
out = HERE / "figs" / "reach_vs_containment.png"
fig.savefig(out, dpi=150)
print(f"wrote {out}")
