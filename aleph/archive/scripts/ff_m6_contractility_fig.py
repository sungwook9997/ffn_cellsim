"""M6 KB-3.14 EMT contractility sweep — migration response vs the KB physiological band.

PI 2026-07-11/12: the KB-grounded attempt at physiological migration (raise myosin into the KB-3.14 2–10× EMT band on
the motile cell-type). This reads the native sweep results (emt_contract{1,5,10}.json, written by run_emt_watcher.sh
when the A5000 frees) and plots v_crawl(contractility) + traction/coherence, with the KB physiological band overlaid
(viz-integrity). Robust to missing files — plots whatever has landed, labels the rest "pending".

Honest by construction: whatever the curve does — rises toward physiological, saturates, or destabilises — is reported
as measured. This is NOT tuned to a target (contractility is set to KB-3.14 band values, not swept to hit a speed).
"""
import json
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

FIGDIR = "aleph/outputs/ff/figs"
MULTS = [1, 5, 10]                       # KB-3.14 band: 1× baseline, 5× mid, 10× max
rows = []
for m in MULTS:
    p = os.path.join(FIGDIR, f"emt_contract{m}.json")
    if not os.path.exists(p):
        rows.append({"mult": m, "pending": True}); continue
    d = json.load(open(p))
    co = d.get("clutch_on", {}); em = d.get("ecm_metrics", {})
    rows.append({"mult": m, "pending": False,
                 "v": co.get("v_crawl_nm_s", np.nan), "disp": co.get("disp_along_um", np.nan) * 1e3,
                 "traction": co.get("traction_nN", np.nan), "bound": co.get("bound_frac", np.nan),
                 "coh": em.get("coh", np.nan), "dens": em.get("dens_nm", np.nan)})

have = [r for r in rows if not r["pending"]]
pending = [r["mult"] for r in rows if r["pending"]]

fig, ax = plt.subplots(1, 2, figsize=(14, 5.2))

# (1) crawl speed vs contractility, KB physiological band overlaid
a = ax[0]
a.axhspan(10, 30, color="#2ca02c", alpha=0.13, zorder=0)
a.text(5.5, 19, "KB physiological\nmesenchymal 10–30 nm/s", fontsize=8.5, color="#2ca02c", va="center")
if have:
    mm = [r["mult"] for r in have]; vv = [r["v"] for r in have]
    a.plot(mm, vv, "o-", color="#d62728", lw=2.4, ms=11, label="FF native v∥ (measured)")
    for r in have:
        a.annotate(f"{r['v']:.2f}", (r["mult"], r["v"]), textcoords="offset points", xytext=(0, 9), fontsize=8, ha="center")
a.set_yscale("log"); a.set_xlabel("EMT myosin contractility multiplier (KB-3.14 band)"); a.set_ylabel("crawl speed v∥ [nm/s]")
a.set_xticks(MULTS); a.set_ylim(0.01, 40)
ttl = "(1) Migration vs EMT contractility (KB-3.14)"
if pending:
    ttl += f"  ·  pending: {pending}×"
a.set_title(ttl, fontsize=10.5); a.legend(fontsize=8, loc="lower right"); a.grid(True, which="both", alpha=0.25)

# (2) the underlying levers — traction, coherence, translocation vs contractility
b = ax[1]
if have:
    mm = [r["mult"] for r in have]
    b.plot(mm, [r["traction"] for r in have], "s-", color="#1f77b4", label="traction [nN]")
    b.plot(mm, [r["disp"] for r in have], "^-", color="#ff7f0e", label="COM disp∥ [nm]")
    b2 = b.twinx()
    b2.plot(mm, [r["coh"] for r in have], "D--", color="#2ca02c", label="remodel coherence")
    b2.set_ylabel("remodel coherence", color="#2ca02c"); b2.tick_params(axis="y", colors="#2ca02c"); b2.set_ylim(0, 1)
    b.legend(fontsize=8, loc="upper left")
b.set_xlabel("EMT contractility multiplier"); b.set_ylabel("traction [nN] · disp∥ [nm]"); b.set_xticks(MULTS)
b.set_title("(2) Traction / translocation / coherence vs contractility", fontsize=10.5); b.grid(True, alpha=0.25)

state = "landed: " + ", ".join(f"{r['mult']}×(v={r['v']:.2f})" for r in have) if have else "NO RESULTS YET (watcher waiting for a free A5000)"
fig.suptitle(f"FF M6 — KB-3.14 EMT contractility migration sweep vs physiological band  ·  {state}", fontsize=11)
fig.tight_layout(rect=[0, 0, 1, 0.94])
out = os.path.join(FIGDIR, "ff_m6_contractility.png")
fig.savefig(out, dpi=130); print(f"wrote {out}  ({len(have)}/{len(MULTS)} results, pending={pending})")
