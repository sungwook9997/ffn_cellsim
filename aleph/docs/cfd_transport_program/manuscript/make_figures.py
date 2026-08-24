#!/usr/bin/env python
"""Publication-quality figures for the CFD-transport / field-actuation manuscript.

All figures are ORIGINAL (schematics + physics plots) generated from the numbers
established in the session's physics analyses; no external/copyrighted images.
Run: python make_figures.py  -> writes figs/fig01..fig13 .png (200 dpi).
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Circle, Rectangle, Ellipse, Wedge
from matplotlib.lines import Line2D

HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(HERE, "figs")
os.makedirs(FIG, exist_ok=True)

# ---- consistent style ----
plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 10, "axes.titlesize": 11,
    "axes.titleweight": "bold", "axes.labelsize": 10, "figure.dpi": 200,
    "savefig.dpi": 200, "savefig.bbox": "tight", "axes.spines.top": False,
    "axes.spines.right": False, "axes.linewidth": 0.9, "mathtext.default": "regular",
})
C = {  # colorblind-safe
    "blue": "#1f6feb", "orange": "#e8710a", "green": "#188038", "red": "#c5221f",
    "purple": "#8e24aa", "teal": "#0f9d9d", "gray": "#5f6b7a", "ltblue": "#d6e6fb",
    "ltorange": "#fde3cf", "ltgreen": "#d9 efd9".replace(" ", ""), "ltgray": "#eef1f5",
}
def save(fig, name):
    fig.savefig(os.path.join(FIG, name), facecolor="white")
    plt.close(fig); print("wrote", name)

def arrow(ax, xy1, xy2, color="k", lw=2, style="-|>", ms=12, alpha=1.0, ls="-"):
    ax.add_patch(FancyArrowPatch(xy1, xy2, arrowstyle=style, mutation_scale=ms,
                 color=color, lw=lw, alpha=alpha, linestyle=ls, zorder=5))

def box(ax, x, y, w, h, text, fc, ec="#33445a", fs=9, tc="#12202e", lw=1.2):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.02,rounding_size=0.03",
                 fc=fc, ec=ec, lw=lw, zorder=3))
    ax.text(x + w/2, y + h/2, text, ha="center", va="center", fontsize=fs, color=tc, zorder=4)

# =====================================================================
# FIG 1 — Overdamped regime: v = F/gamma + the four force sources
# =====================================================================
def fig01():
    fig, (ax, axr) = plt.subplots(1, 2, figsize=(9.2, 4.2), gridspec_kw={"width_ratios":[1.35,1]})
    ax.set_xlim(0, 10); ax.set_ylim(0, 10); ax.axis("off")
    # cell
    ax.add_patch(Circle((5, 5), 3.6, fc=C["ltgray"], ec=C["gray"], lw=1.5))
    ax.text(5, 8.9, "Overdamped cell interior (Re ≈ 10$^{-10}$)", ha="center", fontsize=10, weight="bold")
    # a filament
    fx = np.linspace(3.2, 6.8, 20); fy = 5 + 0.25*np.sin((fx-3.2)*1.5)
    ax.plot(fx, fy, color=C["blue"], lw=3, solid_capstyle="round")
    ax.text(6.9, 5.2, "F-actin", color=C["blue"], fontsize=9)
    # force arrows onto the filament
    arrow(ax, (2.2, 5.0), (3.25, 5.0), color=C["green"], lw=2.4)      # polymerization push
    ax.text(1.2, 5.5, "polymerization\nratchet (push)", color=C["green"], fontsize=8, ha="center")
    arrow(ax, (7.9, 4.6), (6.75, 4.85), color=C["red"], lw=2.4)        # myosin pull
    ax.text(8.7, 4.2, "myosin-II\n(pull/slide)", color=C["red"], fontsize=8, ha="center")
    arrow(ax, (5.0, 7.3), (5.0, 5.6), color=C["purple"], lw=2.0, style="-|>")
    ax.text(5.0, 7.7, "crosslink network (transmit)", color=C["purple"], fontsize=8, ha="center")
    for dx in (-0.5, 0.5):
        arrow(ax, (5+dx, 3.0), (5+dx, 4.4), color=C["gray"], lw=1.2, style="-|>", ms=8, alpha=0.7)
    ax.text(5.0, 2.6, "thermal (Brownian)", color=C["gray"], fontsize=8, ha="center")
    # right: equation panel
    axr.axis("off"); axr.set_xlim(0,1); axr.set_ylim(0,1)
    axr.text(0.5, 0.9, "No inertia — force sets velocity", ha="center", weight="bold", fontsize=10.5)
    axr.text(0.5, 0.72, r"$\gamma\,\dot{x} = F_{\rm total}$", ha="center", fontsize=20, color=C["blue"])
    axr.text(0.5, 0.55, r"$\Rightarrow\ v = F/\gamma$", ha="center", fontsize=17)
    axr.text(0.5, 0.40, "moves only while a force is applied;\nstops instantly when it ceases",
             ha="center", fontsize=9, color="#33445a")
    axr.text(0.5, 0.20, r"$\gamma = 6\pi\eta R / N_c$", ha="center", fontsize=13)
    axr.text(0.5, 0.09, r"$\eta_{\rm MCF7}\approx 65.9\ \mathrm{Pa\cdot s}$  (no free parameter)",
             ha="center", fontsize=8.5, color="#33445a")
    axr.add_patch(FancyBboxPatch((0.04,0.03),0.92,0.94, boxstyle="round,pad=0.01",
                  fc="#f7fafd", ec=C["blue"], lw=1.2))
    fig.suptitle("Figure 1  |  How filaments actually move inside a cell", x=0.02, ha="left", fontsize=11.5)
    save(fig, "fig01_overdamped_forces.png")

# =====================================================================
# FIG 2 — Retrograde flow + molecular clutch (schematic + biphasic curve)
# =====================================================================
def fig02():
    fig, (ax, axr) = plt.subplots(1, 2, figsize=(9.4, 4.0), gridspec_kw={"width_ratios":[1.3,1]})
    ax.set_xlim(0, 10); ax.set_ylim(0, 6.5); ax.axis("off")
    # membrane leading edge
    ax.plot([1, 9], [5.6, 5.6], color="#111", lw=2.5)
    ax.text(9.1, 5.6, "membrane", fontsize=8, va="center")
    # actin network flowing rearward
    for i, x0 in enumerate(np.linspace(2.0, 8.0, 7)):
        ax.plot([x0, x0-1.1], [5.3, 3.6], color=C["blue"], lw=2, alpha=0.8)
    arrow(ax, (7.2, 4.4), (3.4, 4.4), color=C["blue"], lw=2.6)
    ax.text(5.3, 4.75, "retrograde F-actin flow  (10–140 nm/s)", color=C["blue"], fontsize=9, ha="center")
    # polymerization at front
    arrow(ax, (1.6, 5.35), (2.6, 5.35), color=C["green"], lw=2.2, ms=10)
    ax.text(1.7, 5.9, "polymerization\npush (front)", color=C["green"], fontsize=8)
    # myosin at rear
    ax.add_patch(Ellipse((7.7, 3.4), 0.9, 0.35, fc=C["red"], ec="k", lw=0.6))
    ax.text(8.6, 3.3, "myosin-II\n(rear)", color=C["red"], fontsize=8)
    # clutch
    for x0 in (3.2, 4.6, 6.0):
        ax.plot([x0, x0], [3.6, 2.6], color=C["orange"], lw=2)
        ax.add_patch(Circle((x0, 2.5), 0.14, fc=C["orange"], ec="k", lw=0.5))
    ax.plot([2.5, 7.0], [2.3, 2.3], color="#5a4632", lw=4)
    ax.text(4.75, 1.9, "ECM / substrate", fontsize=8, ha="center", color="#5a4632")
    ax.text(4.6, 3.0, "clutch (integrin–talin–vinculin)", color=C["orange"], fontsize=8, ha="center")
    ax.text(5, 6.25, "Molecular-clutch traction transmission", ha="center", fontsize=10, weight="bold")
    # right: biphasic traction vs flow (Gardel/Chan-Odde)
    v = np.linspace(0, 60, 300)
    thr = 9.0
    trac = np.where(v < thr, v/thr, np.exp(-(v-thr)/22))  # rise then fall (biphasic)
    axr.plot(v, trac, color=C["red"], lw=2.4)
    axr.axvline(thr, ls="--", color=C["gray"], lw=1)
    axr.text(thr+1, 0.15, "threshold\n~8–10 nm/s", fontsize=8, color=C["gray"])
    axr.set_xlabel("retrograde flow speed (nm/s)"); axr.set_ylabel("traction (norm.)")
    axr.set_title("Biphasic clutch  (Gardel 2008; Chan–Odde 2008)", fontsize=9.5)
    axr.set_ylim(0, 1.15); axr.annotate("grip", (thr*0.5, 0.55), fontsize=8, color=C["green"])
    axr.annotate("slip", (35, 0.45), fontsize=8, color=C["blue"])
    fig.suptitle("Figure 2  |  Retrograde flow + molecular clutch — the real migration engine",
                 x=0.02, ha="left", fontsize=11.5)
    save(fig, "fig02_retrograde_clutch.png")

# =====================================================================
# FIG 3 — Stress-fiber three-route assembly
# =====================================================================
def fig03():
    fig, axes = plt.subplots(1, 3, figsize=(10.5, 3.7))
    titles = ["Dorsal SF", "Transverse arc", "Ventral SF"]
    subs = ["formin (mDia1) de novo\nat the focal adhesion",
            "condensation of lamellipodial\nfilaments by myosin + α-actinin",
            "fusion of 2 dorsal SF + arc\n+ sarcomeric maturation"]
    for ax, t, s in zip(axes, titles, subs):
        ax.set_xlim(0, 10); ax.set_ylim(0, 10); ax.axis("off")
        ax.set_title(t, fontsize=10.5)
        ax.text(5, 0.4, s, ha="center", fontsize=8, color="#33445a")
    # dorsal
    ax = axes[0]
    ax.add_patch(Ellipse((5, 2.0), 2.6, 0.8, fc=C["orange"], ec="k", lw=0.6))
    ax.text(5, 2.0, "FA", ha="center", va="center", fontsize=8, color="w", weight="bold")
    for dx in (-0.15, 0, 0.15):
        ax.plot([5+dx, 5+dx*3], [2.4, 8.5], color=C["blue"], lw=2)
    ax.add_patch(Circle((5, 8.6), 0.22, fc=C["green"], ec="k"))
    ax.text(6.1, 8.5, "formin\n(barbed tip)", fontsize=7.5, color=C["green"])
    ax.annotate("", (5.1, 8.2), (5.05, 6.5), arrowprops=dict(arrowstyle="-|>", color=C["blue"], lw=2))
    ax.text(3.0, 6.0, "grows\nfrom FA", fontsize=8, color=C["blue"])
    # arc
    ax = axes[1]
    for y0, a in zip(np.linspace(7.5, 5.5, 4), [0.4,0.6,0.8,1.0]):
        xs = np.linspace(2, 8, 30); ys = y0 + 0.3*np.random.default_rng(int(y0*7)).standard_normal(30)*0.0
        ax.plot(xs, np.full_like(xs, y0), color=C["blue"], lw=1.4, alpha=a)
    arrow(ax, (5, 5.2), (5, 3.6), color=C["red"], lw=2.2)
    ax.text(5.4, 4.3, "myosin\ncondenses", fontsize=8, color=C["red"])
    ax.plot(np.linspace(2,8,50), 3.0+0.0*np.linspace(2,8,50), color=C["blue"], lw=3.2)
    ax.text(5, 2.4, "condensed arc", ha="center", fontsize=8, color=C["blue"])
    # ventral
    ax = axes[2]
    for x0 in (2.5, 7.5):
        ax.add_patch(Ellipse((x0, 2.0), 1.6, 0.7, fc=C["orange"], ec="k", lw=0.6))
        ax.text(x0, 2.0, "FA", ha="center", va="center", fontsize=7.5, color="w", weight="bold")
    ax.plot([2.5, 7.5], [2.3, 2.3], color=C["blue"], lw=3.4)
    # sarcomeric alpha-actinin + myosin bands
    for x0 in np.linspace(3.1, 6.9, 5):
        ax.plot([x0, x0], [2.1, 2.5], color=C["purple"], lw=2)   # z-body
    for x0 in np.linspace(3.6, 6.4, 4):
        ax.add_patch(Ellipse((x0, 2.3), 0.35, 0.16, fc=C["red"], ec="k", lw=0.4))
    ax.text(5, 3.3, "α-actinin Z-bodies (bars) + NMIIA bands (ovals)", ha="center", fontsize=7.5)
    fig.suptitle("Figure 3  |  Stress fibers assemble by three routes — not by filaments 'escaping' the cortex",
                 x=0.02, ha="left", fontsize=11.5)
    save(fig, "fig03_sf_three_routes.png")

# =====================================================================
# FIG 4 — Intracellular transport regime map (Peclet)
# =====================================================================
def fig04():
    fig, ax = plt.subplots(figsize=(7.6, 4.6))
    v = np.logspace(-3, 3, 400)   # um/s
    D = 5.0  # um^2/s
    for L, ls, lab in [(1.0, "-", "L = 1 µm (near FA)"), (10.0, "--", "L = 10 µm (cell)")]:
        Pe = v * L / D
        ax.loglog(v, Pe, ls, color=C["blue"], lw=2, label=lab)
    ax.axhline(1, color=C["red"], lw=1.5)
    ax.text(1.2e-3, 1.5, "Pe = 1  (advection = diffusion)", color=C["red"], fontsize=8.5)
    ax.axhspan(1e-4, 1, color=C["ltblue"], alpha=0.5)
    ax.axhspan(1, 1e5, color=C["ltorange"], alpha=0.4)
    ax.text(3e2, 3e-3, "DIFFUSION-dominated", color=C["blue"], fontsize=9, ha="right", weight="bold")
    ax.text(3e2, 2e4, "ADVECTION matters", color=C["orange"], fontsize=9, ha="right", weight="bold")
    # markers
    pts = [(0.03, "slow mesenchymal\nanimal cell", C["gray"]),
           (0.5, "v$_{thr}$ = D/L", C["red"]),
           (75, "plant streaming\n(myosin-XI)", C["green"])]
    for vp, lab, col in pts:
        ax.axvline(vp, color=col, ls=":", lw=1.1, alpha=0.8)
        ax.text(vp, 4e4, lab, rotation=90, fontsize=7.5, color=col, va="top", ha="right")
    ax.set_xlabel("flow speed v (µm/s)"); ax.set_ylabel("Péclet number  Pe = vL/D")
    ax.set_ylim(1e-4, 1e5); ax.legend(loc="lower right", fontsize=8.5, frameon=False)
    ax.set_title("Figure 4  |  Advection beats diffusion only above v ≈ 0.5 µm/s (D$_{G\\text{-actin}}$≈5 µm²/s)",
                 fontsize=10.5, loc="left")
    save(fig, "fig04_peclet_regime.png")

# =====================================================================
# FIG 5 — Three-gate framework for field actuation
# =====================================================================
def fig05():
    fig, ax = plt.subplots(figsize=(9.6, 4.4)); ax.set_xlim(0, 12); ax.set_ylim(0, 6); ax.axis("off")
    box(ax, 0.2, 2.6, 1.9, 1.0, "Idea:\nfield drives\ninternal flow", C["ltgray"], fs=8.5)
    gates = [
        (2.7, "Gate 1\nPéclet\nv ≳ 0.5 µm/s?", C["ltblue"], "needs fast\ncoherent flow", C["blue"]),
        (5.4, "Gate 2\nMembrane\nshielding", C["ltorange"], "DC attenuated\n~10$^3$–10$^4$×", C["red"]),
        (8.1, "Gate 3\nRate-limiting\nat FA?", C["ltgreen"], "formin/reaction\nlimited (ε≲0.2)", C["red"]),
    ]
    xprev = 2.1
    for x, t, fc, note, nc in gates:
        box(ax, x, 2.6, 1.9, 1.0, t, fc, fs=8.5)
        arrow(ax, (xprev, 3.1), (x, 3.1), lw=1.8)
        ax.text(x+0.95, 2.2, note, ha="center", fontsize=7.5, color=nc)
        xprev = x + 1.9
    arrow(ax, (10.0, 3.1), (10.6, 3.1), lw=1.8)
    box(ax, 10.6, 3.3, 1.3, 1.4, "(d2) roll\nby field\n= DEAD", "#f7d7d5", ec=C["red"], tc=C["red"], fs=8.5)
    box(ax, 10.6, 0.9, 1.3, 1.6, "survivor:\nmotor-driven\nRAD transport\n(modelable)", "#d7efd9", ec=C["green"], tc=C["green"], fs=8)
    arrow(ax, (9.05, 2.55), (11.2, 2.55), lw=1.6, color=C["green"], ls=(0,(4,2)))
    ax.text(5.9, 5.4, "Every external-field 'drive the interior directly' idea must pass three quantitative gates",
            ha="center", fontsize=9.5, weight="bold")
    fig.suptitle("Figure 5  |  A three-gate physical test for field actuation of the cytoskeleton",
                 x=0.02, ha="left", fontsize=11.5)
    save(fig, "fig05_three_gates.png")

# =====================================================================
# FIG 6 — Reynolds number: turbulence forbidden
# =====================================================================
def fig06():
    fig, ax = plt.subplots(figsize=(7.8, 4.4))
    labels = ["MCF7 cytoplasm\nη=65.9, v=1µm/s", "water η,\nv=1µm/s",
              "water η, v=1mm/s\n(1000× fantasy)", "turbulence\nonset"]
    vals = [1.5e-10, 1.0e-5, 1.0e-2, 2.0e3]
    cols = [C["blue"], C["teal"], C["orange"], C["red"]]
    y = np.arange(len(vals))[::-1]
    ax.barh(y, vals, color=cols, height=0.55, log=True, zorder=3)
    for yi, v, l in zip(y, vals, labels):
        ax.text(v*1.6, yi, f"Re ≈ {v:.0e}".replace("e-0","e-").replace("e+0","e"),
                va="center", fontsize=8.5)
    ax.set_yticks(y); ax.set_yticklabels(labels, fontsize=8.5)
    ax.set_xlabel("Reynolds number  (log)"); ax.set_xlim(1e-11, 1e6)
    ax.axvline(2e3, color=C["red"], ls="--", lw=1.2)
    ax.text(2e3, 3.35, "Re ≳ 2000 needed", color=C["red"], fontsize=8.5, ha="center")
    ax.annotate("", (1.5e-10, 0.0), (2e3, 0.0),
                arrowprops=dict(arrowstyle="<->", color=C["gray"], lw=1.2))
    ax.text(3e-4, 0.28, "~10–13 orders short → turbulence physically forbidden",
            fontsize=8.5, color=C["gray"], ha="center")
    ax.set_title("Figure 6  |  Cytosol is deep-Stokes: internal turbulence cannot exist at cell scale",
                 fontsize=10.5, loc="left")
    save(fig, "fig06_reynolds.png")

# =====================================================================
# FIG 7 — Membrane dielectric shielding vs frequency (Schwan)
# =====================================================================
def fig07():
    fig, ax = plt.subplots(figsize=(7.8, 4.4))
    f = np.logspace(0, 9, 400)   # Hz
    fc = 1.6e6                    # Maxwell-Wagner crossover
    ratio = 2e-4 + (1 - 2e-4) / np.sqrt(1 + (fc/np.maximum(f,1))**2)
    ax.loglog(f, ratio, color=C["orange"], lw=2.4)
    ax.axhline(2e-4, color=C["gray"], ls=":", lw=1)
    ax.text(2, 3e-4, "DC interior field\nE$_{int}$/E$_{ext}$ ≈ 2×10$^{-4}$", fontsize=8, color=C["gray"])
    ax.axvline(fc, color=C["red"], ls="--", lw=1.2)
    ax.text(fc*1.2, 1e-3, "β-dispersion\ncrossover ~1.6 MHz", fontsize=8, color=C["red"])
    ax.axhspan(1e-4, 2e-4*3, color=C["ltorange"], alpha=0.4)
    ax.text(3, 1.2e-4, "cell interior shielded (DC / low-f)", fontsize=8.5, color=C["orange"])
    ax.text(3e7, 0.5, "AC penetrates\n→ but zero net drift", fontsize=8.5, color="#33445a", ha="center")
    ax.set_xlabel("applied field frequency (Hz)")
    ax.set_ylabel("interior / external field  E$_{int}$/E$_{ext}$")
    ax.set_ylim(1e-4, 2); ax.set_title("Figure 7  |  The membrane dielectrically shields the interior from external DC fields",
                                       fontsize=10.5, loc="left")
    save(fig, "fig07_membrane_shielding.png")

# =====================================================================
# FIG 8 — Electro-orientation + Debye screening collapse
# =====================================================================
def fig08():
    fig, ax = plt.subplots(figsize=(7.8, 4.4))
    I = np.logspace(-1, 2.6, 300)   # mM ionic strength
    debye = 0.304 / np.sqrt(I/1000.0)  # nm (approx, water 25C)  lambda_D = 0.304/sqrt(I[M])
    ax2 = ax.twinx()
    ax.loglog(I, debye, color=C["blue"], lw=2.4)
    ax.set_ylabel("Debye length λ$_D$ (nm)", color=C["blue"])
    ax.tick_params(axis="y", colors=C["blue"])
    # orientation "torque efficiency" ~ collapses when debye << filament scale
    eff = 1/(1+ (0.8/debye)**2)   # arbitrary monotone collapse with screening
    ax2.plot(I, eff, color=C["red"], lw=2.2, ls="--")
    ax2.set_ylabel("electro-orientation torque (rel.)", color=C["red"])
    ax2.tick_params(axis="y", colors=C["red"]); ax2.set_ylim(-0.03, 1.05)
    ax.axvspan(140, 160, color=C["gray"], alpha=0.25)
    ax.text(150, 30, "physiological\n150 mM", ha="center", fontsize=8, color="#33445a")
    ax.text(0.15, 2.5, "low salt:\nfilaments align\n(500–1900 V/cm)", fontsize=8, color=C["green"])
    ax.text(80, 0.5, "λ$_D$≈0.8 nm\n→ torque collapses", fontsize=8, color=C["red"], ha="center")
    ax.set_xlabel("ionic strength (mM)")
    ax.set_title("Figure 8  |  F-actin electro-orients only at low salt; physiological screening kills it",
                 fontsize=10.3, loc="left")
    save(fig, "fig08_electro_orientation.png")

# =====================================================================
# FIG 9 — Galvanotaxis: signaling-mediated, not direct force
# =====================================================================
def fig09():
    fig, ax = plt.subplots(figsize=(9.4, 4.0)); ax.set_xlim(0,12); ax.set_ylim(0,6); ax.axis("off")
    ax.text(3.0, 5.6, "External DC field (0.1–6 V/cm)", ha="center", fontsize=9.5, weight="bold")
    arrow(ax, (0.6, 3.3), (5.2, 3.3), color=C["purple"], lw=3, ms=16)
    ax.text(2.9, 3.65, "E", color=C["purple"], fontsize=12, style="italic")
    # membrane with receptors redistributing
    ax.add_patch(Circle((7.4, 3.0), 2.2, fc=C["ltgray"], ec="k", lw=1.5))
    for ang, s in [(160,0.5),(180,0.6),(200,0.5),(150,0.4),(210,0.4)]:
        x = 7.4 + 2.2*np.cos(np.radians(ang)); y = 3.0 + 2.2*np.sin(np.radians(ang))
        ax.add_patch(Circle((x,y), 0.13*(1+s), fc=C["red"], ec="k", lw=0.4))
    for ang in (10,-10,30,-30):
        x = 7.4 + 2.2*np.cos(np.radians(ang)); y = 3.0 + 2.2*np.sin(np.radians(ang))
        ax.add_patch(Circle((x,y), 0.10, fc=C["red"], ec="k", lw=0.3, alpha=0.5))
    ax.text(5.0, 1.0, "charged membrane receptors\nelectrophorese → asymmetry", fontsize=8, color=C["red"], ha="center")
    ax.text(9.9, 4.3, "PI3Kγ / PTEN\npolarized signaling", fontsize=8, color=C["green"], ha="center")
    arrow(ax, (7.4, 3.0), (9.2, 3.0), color=C["green"], lw=2)
    ax.text(9.6, 3.0, "biased actin\npolymerization", fontsize=8, color=C["blue"])
    # crossed-out direct force
    box(ax, 0.4, 0.2, 3.2, 0.9, "direct electrophoretic\nforce on cytoskeleton", "#f7d7d5", ec=C["red"], fs=8)
    ax.plot([0.5,3.5],[0.25,1.05], color=C["red"], lw=2.5); ax.plot([0.5,3.5],[1.05,0.25], color=C["red"], lw=2.5)
    ax.text(2.0, -0.15, "interior shielded → NOT this", fontsize=7.5, color=C["red"], ha="center")
    ax.text(6.0, 5.9, "Galvanotaxis is a SIGNAL the cell's own engine responds to", ha="center", fontsize=9.5, weight="bold")
    fig.suptitle("Figure 9  |  Galvanotaxis is signaling-mediated — the field steers, it does not push",
                 x=0.02, ha="left", fontsize=11.5)
    save(fig, "fig09_galvanotaxis.png")

# =====================================================================
# FIG 10 — Reaction- vs transport-limited SF growth
# =====================================================================
def fig10():
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(9.6, 4.0))
    # panel A: FA sink draws only ~6% of diffusion ceiling
    a1.bar([0,1], [3800, 232], color=[C["blue"], C["orange"]], width=0.6, log=True, zorder=3)
    a1.set_xticks([0,1]); a1.set_xticklabels(["diffusion\nceiling I$_{max}$", "one barbed-end\nconsumption J"])
    a1.set_ylabel("monomer flux (subunits/s)")
    a1.text(0, 4600, "3,800/s", ha="center", fontsize=9)
    a1.text(1, 300, "232/s\n(~6%)", ha="center", fontsize=9)
    a1.set_title("FA is a weak sink", fontsize=10)
    a1.set_ylim(50, 1e4)
    # panel B: flux elasticity vs regime
    regimes = ["resting FA", "global pool\ndepletion", "fast\nlamellipodium\n(Pe≈3)"]
    eps = [0.15, 0.9, 0.7]
    cols = [C["red"], C["green"], C["green"]]
    a2.bar(range(3), eps, color=cols, width=0.6, zorder=3)
    a2.axhline(1.0, ls=":", color=C["gray"], lw=1); a2.text(2.4, 1.02, "proportional", fontsize=7.5, color=C["gray"])
    a2.set_xticks(range(3)); a2.set_xticklabels(regimes, fontsize=8.5)
    a2.set_ylabel("flux elasticity  d ln(growth)/d ln(flux)")
    a2.set_title("Flux is a weak lever at a resting FA", fontsize=10)
    a2.set_ylim(0, 1.15)
    a2.text(0, 0.18, "≲0.2\n(buffered)", ha="center", fontsize=8, color=C["red"])
    fig.suptitle("Figure 10  |  SF growth at a resting FA is reaction/formin-limited, not transport-limited",
                 x=0.02, ha="left", fontsize=11.5)
    fig.tight_layout(rect=(0,0,1,0.93))
    save(fig, "fig10_reaction_limited.png")

# =====================================================================
# FIG 11 — RAD monomer-transport engine architecture
# =====================================================================
def fig11():
    fig, ax = plt.subplots(figsize=(9.8, 4.6)); ax.set_xlim(0,12); ax.set_ylim(0,7); ax.axis("off")
    box(ax, 0.3, 4.6, 2.3, 1.2, "Biot pressure\ngrid  ∂p/∂t=c$_v$∇²p\n[EXISTS]", C["ltgreen"], fs=8.2)
    box(ax, 3.2, 4.6, 2.3, 1.2, "Darcy velocity\nu = −(k/µ)∇p\n[derivable now]", C["ltgreen"], fs=8.2)
    box(ax, 0.3, 2.3, 2.3, 1.2, "G-actin field\nc(x,t)  RAD\n[clone Biot]", C["ltblue"], fs=8.2)
    box(ax, 3.2, 2.3, 2.3, 1.2, "advection\n−u·∇c\n[new, days]", C["ltblue"], fs=8.2)
    box(ax, 6.2, 3.4, 2.4, 1.2, "FA barbed-end\nsink (formin)\n[clone source kernel]", C["ltorange"], fs=8.2)
    box(ax, 9.1, 3.4, 2.5, 1.2, "growth v$_0$=f(local c)\n→ emergent SF\n[un-scalarize]", C["ltorange"], fs=8.2)
    arrow(ax, (2.6, 5.2), (3.2, 5.2)); arrow(ax, (4.35, 4.55), (4.35, 3.55), color=C["blue"])
    arrow(ax, (2.6, 2.9), (3.2, 2.9)); arrow(ax, (5.5, 3.0), (6.2, 3.6), color=C["orange"])
    arrow(ax, (5.5, 5.0), (6.4, 4.6), color=C["green"])
    arrow(ax, (8.6, 4.0), (9.1, 4.0))
    ax.text(6.0, 6.4, "∂c/∂t = D∇²c − u·∇c + src − sink$_{FA}$",
            ha="center", fontsize=13, color=C["blue"])
    ax.text(6.0, 0.9, "~70% of the plumbing already ships (green); Darcy-u build ≈ 1–2 weeks.  "
                      "Independent u–p velocity DOF + electro-osmotic constants = PI-gated.",
            ha="center", fontsize=8, color="#33445a")
    leg = [Line2D([0],[0], marker="s", color="w", markerfacecolor=C["ltgreen"], markersize=12, label="exists"),
           Line2D([0],[0], marker="s", color="w", markerfacecolor=C["ltblue"], markersize=12, label="new (reuse)"),
           Line2D([0],[0], marker="s", color="w", markerfacecolor=C["ltorange"], markersize=12, label="coupling")]
    ax.legend(handles=leg, loc="lower right", fontsize=8, frameon=False)
    fig.suptitle("Figure 11  |  The salvageable core: a reaction-advection-diffusion G-actin field in the CFD engine",
                 x=0.02, ha="left", fontsize=11.2)
    save(fig, "fig11_rad_engine.png")

# =====================================================================
# FIG 12 — TFM discriminator (traction kinetics)
# =====================================================================
def fig12():
    fig, ax = plt.subplots(figsize=(7.8, 4.4))
    t = np.linspace(-2, 12, 400)
    field_on = (t >= 0) & (t <= 8)
    # signaling: gradual rise over minutes, relax after off
    sig = np.zeros_like(t)
    sig[t>=0] = 1 - np.exp(-(t[t>=0])/2.2)
    sig[t>8] = sig[np.argmin(np.abs(t-8))]*np.exp(-(t[t>8]-8)/2.2)
    # direct force: instantaneous step
    direct = np.where(field_on, 1.0, 0.0)
    ax.axvspan(0, 8, color="#eef4fb", zorder=0)
    ax.text(4, 1.12, "field ON", ha="center", fontsize=9, color=C["gray"])
    ax.plot(t, sig, color=C["green"], lw=2.6, label="signaling-mediated\n(gradual, min; drug-sensitive)")
    ax.plot(t, direct, color=C["red"], lw=2.2, ls="--", label="direct force\n(instantaneous; drug-insensitive)")
    ax.set_xlabel("time (min)"); ax.set_ylabel("traction-polarity reorientation (norm.)")
    ax.set_ylim(-0.05, 1.25); ax.legend(loc="center right", fontsize=8, frameon=False)
    ax.set_title("Figure 12  |  TFM traction on/off kinetics discriminate signaling from direct force",
                 fontsize=10.3, loc="left")
    ax.annotate("the discriminator", (1.2, 0.45), (3.2, 0.25),
                arrowprops=dict(arrowstyle="->", color=C["gray"]), fontsize=8.5, color=C["gray"])
    save(fig, "fig12_tfm_discriminator.png")

# =====================================================================
# FIG 13 — Unified physics thread
# =====================================================================
def fig13():
    fig, ax = plt.subplots(figsize=(9.8, 4.6)); ax.set_xlim(0,12); ax.set_ylim(0,7); ax.axis("off")
    ideas = [("(d2) turbulent cytosol\n→ roll the cell", "Re≈10$^{-10}$ + scallop\n+ Biot ≠ Navier-Stokes", 5.4),
             ("field → electrophorese\nmonomers", "membrane shielding\n10$^3$–10$^4$×; Debye collapse", 3.5),
             ("field flux → control\nSF growth", "reaction-limited at FA\n(ε ≲ 0.2)", 1.6)]
    for t, wall, y in ideas:
        box(ax, 0.3, y-0.55, 2.7, 1.1, t, C["ltgray"], fs=8.2)
        arrow(ax, (3.05, y), (4.5, y), color=C["red"], lw=1.8)
        box(ax, 4.5, y-0.55, 3.2, 1.1, wall, "#f7d7d5", ec=C["red"], tc=C["red"], fs=7.8)
        ax.text(6.1, y-0.85, "✗ WALL", ha="center", fontsize=7.5, color=C["red"])
    # converge to survivor
    box(ax, 8.6, 2.3, 3.1, 1.9, "SURVIVOR\nmotor-driven advective\ntransport of a\nmass-conserving\nG-actin field", "#d7efd9",
        ec=C["green"], tc=C["green"], fs=8.5)
    for y in (5.4, 3.5, 1.6):
        arrow(ax, (7.75, y), (8.55, 3.2), color=C["green"], lw=1.3, ls=(0,(3,2)), alpha=0.7)
    ax.text(6.0, 6.6, "Three field-actuation ideas, one recurring physical wall — and the modelable path that survives",
            ha="center", fontsize=9.5, weight="bold")
    ax.text(10.15, 1.9, "'directed advection is a MOTOR problem, not a FIELD problem'", ha="center",
            fontsize=7.6, style="italic", color=C["green"])
    fig.suptitle("Figure 13  |  The unified physics thread of the field-actuation analysis",
                 x=0.02, ha="left", fontsize=11.5)
    save(fig, "fig13_unified_thread.png")

for fn in [fig01,fig02,fig03,fig04,fig05,fig06,fig07,fig08,fig09,fig10,fig11,fig12,fig13]:
    fn()
print("ALL FIGURES DONE ->", FIG)
