"""M6 follow-up — adhesion-density (n_fa) sweep: what more FA clutches actually buy.

After the M6 contractility sweep showed migration is adhesion-limited (traction≈0 at n_fa=100), this sweep raises the
adhesion count (n_fa = 100 → 500 → 1000) at baseline contractility to see what higher adhesion does. Reads the native
jsons (emt_contract1 = n_fa 100, emt_nfa500, emt_nfa1000). Robust to missing files.

The honest story it tells: more static adhesion buys **traction and matrix REMODELLING**, but NOT migration — the cell
grips harder and pulls the collagen in place rather than translocating (migration needs the dynamic front-grip/rear-
release walking cycle, not more static clutches). Reported as measured.
"""
import json
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

FIGDIR = "aleph/outputs/ff/figs"
POINTS = [(100, "emt_contract1"), (500, "emt_nfa500"), (1000, "emt_nfa1000")]
rows = []
for nfa, tag in POINTS:
    p = os.path.join(FIGDIR, f"{tag}.json")
    if not os.path.exists(p):
        rows.append({"nfa": nfa, "pending": True}); continue
    d = json.load(open(p)); co = d.get("clutch_on", {}); em = d.get("ecm_metrics", {})
    rows.append({"nfa": nfa, "pending": False, "v": co.get("v_crawl_nm_s", np.nan),
                 "traction": co.get("traction_nN", np.nan), "dens": em.get("dens_nm", np.nan),
                 "recruit": em.get("recruit_nm", np.nan), "bound": co.get("bound_frac", np.nan)})
have = [r for r in rows if not r["pending"]]
pending = [r["nfa"] for r in rows if r["pending"]]

fig, ax = plt.subplots(1, 2, figsize=(14, 5.2))
nn = [r["nfa"] for r in have]

# (1) traction + matrix remodelling RISE with adhesion
a = ax[0]
if have:
    a.plot(nn, [r["traction"] for r in have], "s-", color="#1f77b4", lw=2.2, ms=9, label="traction [nN]")
    a2 = a.twinx()
    a2.plot(nn, [r["dens"] for r in have], "o-", color="#d62728", lw=2.2, ms=9, label="densification [nm]")
    a2.plot(nn, [r["recruit"] for r in have], "^--", color="#ff7f0e", lw=1.8, ms=8, label="recruit [nm]")
    a2.set_ylabel("matrix remodel [nm]", color="#d62728"); a2.tick_params(axis="y", colors="#d62728")
    a.legend(loc="upper left", fontsize=8); a2.legend(loc="lower right", fontsize=8)
a.set_xlabel("FA adhesion count  n_fa"); a.set_ylabel("traction [nN]", color="#1f77b4")
a.set_title("(1) More adhesion → more traction + matrix REMODELLING", fontsize=10.5); a.grid(True, alpha=0.25)

# (2) ...but migration does NOT rise — it stays ~0 (vs the KB band)
b = ax[1]
b.axhspan(10, 30, color="#2ca02c", alpha=0.13, zorder=0)
b.text(np.mean(nn) if nn else 500, 19, "KB physiological 10–30 nm/s", fontsize=8.5, color="#2ca02c", ha="center")
if have:
    b.plot(nn, [max(r["v"], 1e-3) for r in have], "D-", color="#9467bd", lw=2.4, ms=10, label="v∥ (measured)")
    for r in have:
        b.annotate(f"{r['v']:.2f}", (r["nfa"], max(r["v"], 1e-3)), textcoords="offset points", xytext=(0, 9), fontsize=8, ha="center")
b.set_yscale("log"); b.set_ylim(1e-3, 40); b.set_xlabel("FA adhesion count  n_fa"); b.set_ylabel("crawl speed v∥ [nm/s]")
b.set_title("(2) ...but MIGRATION stays ≈0 (needs walking cycle, not static clutches)", fontsize=10.5)
b.legend(fontsize=8, loc="lower left"); b.grid(True, which="both", alpha=0.25)

state = "n_fa " + "/".join(str(r["nfa"]) for r in have) if have else "none"
if pending: state += f"  (pending {pending})"
fig.suptitle(f"FF M6 adhesion-density sweep — adhesion buys TRACTION + REMODELLING, not MIGRATION  ·  {state}", fontsize=11)
fig.tight_layout(rect=[0, 0, 1, 0.94])
out = os.path.join(FIGDIR, "ff_m6_adhesion.png")
fig.savefig(out, dpi=130); print(f"wrote {out}  ({len(have)}/{len(POINTS)}, pending={pending})")
