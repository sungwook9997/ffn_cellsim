"""Build the paper-ready Layer-2 synthesis as a Word (.docx) report.

Consolidates the Layer-2 multicellular-CBM line: the PI spreading-law FORM reproduced from a fully
mechanistic model (zero calibration), the magnitude gap bounded as the center-based structural limit
on six independent axes, the aggregate-wetting explanation, and the (c) ligand-condition FORM-level
separation. Reads the committed result files (production fit, ligand_prod.jsonl, PI medians) at build
time and embeds the key figures. Output: outputs/layer2/Layer2_Report.docx.

Usage:  python -m aleph.scripts.layer2_report_docx
"""

from __future__ import annotations

import json
from collections import defaultdict
from pathlib import Path

from docx import Document
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.shared import Inches, Pt, RGBColor

_ROOT = Path(__file__).resolve().parents[1]
_OUT = _ROOT / "outputs" / "layer2"
_FIGS = _OUT / "figs"
_DOCX = _OUT / "Layer2_Report.docx"


def _pi_medians():
    p = _OUT / "pi_overlay_summary.json"
    if not p.exists():
        return {}
    d = json.loads(p.read_text())
    return {c: d[c]["pi_AA_med"] for c in ("Bare", "Pre", "Lam4")
            if c in d and isinstance(d[c], dict) and "pi_AA_med" in d[c]}


def _ligand_table():
    """(c) ligand-production results, aggregated by condition+N0 from ligand_prod.jsonl."""
    p = _OUT / "ligand_prod" / "ligand_prod.jsonl"
    if not p.exists():
        return None
    rows = [json.loads(l) for l in p.read_text().splitlines() if l.strip().startswith("{")]
    if not rows:
        return None
    by = defaultdict(lambda: defaultdict(list))
    for r in rows:
        by[r["condition"]][r["n"]].append(r)
    import numpy as np
    out = {}
    for c in ("Bare", "Pre", "Lam4"):
        if c not in by:
            continue
        ns = sorted(by[c])
        out[c] = [(int(np.mean([x["R0_um"] for x in by[c][n]])),
                   float(np.mean([x["aa0_core"] for x in by[c][n]])),
                   float(np.std([x["aa0_core"] for x in by[c][n]]))) for n in ns]
    return out


def _h(doc, text, level):
    doc.add_heading(text, level=level)


def _fig(doc, name, caption):
    path = _FIGS / name
    if path.exists():
        doc.add_picture(str(path), width=Inches(6.2))
        doc.paragraphs[-1].alignment = WD_ALIGN_PARAGRAPH.CENTER
        cap = doc.add_paragraph()
        run = cap.add_run(caption)
        run.italic = True; run.font.size = Pt(9)
        cap.alignment = WD_ALIGN_PARAGRAPH.CENTER
    else:
        doc.add_paragraph(f"[figure {name} not found]")


def build() -> Path:
    doc = Document()
    # base style
    doc.styles["Normal"].font.name = "Calibri"
    doc.styles["Normal"].font.size = Pt(10.5)

    title = doc.add_heading(
        "Mechanistic reproduction of the MCF7 spreading law and a six-axis bounding "
        "of the magnitude gap in a center-based multicellular spheroid model (Layer-2)", level=0)
    sub = doc.add_paragraph()
    r = sub.add_run("ffn_cellsim Layer-2 line · branch layer2/spheroid-cbm · 2026-06-04 · "
                    "all PI experimental A/A₀ comparisons are overlay-only (never fitted)")
    r.italic = True; r.font.size = Pt(9); r.font.color.rgb = RGBColor(0x55, 0x55, 0x55)

    # ---- Abstract ----
    _h(doc, "Abstract", 1)
    doc.add_paragraph(
        "The PI MCF7-spheroid-on-substrate platform exhibits a genuinely novel spreading law, "
        "A/A₀ = a + b/R + c/R² (initial-radius dependence of the projected spread-area ratio), with no "
        "published closed-form analog. We reproduce this law’s FORM from a fully mechanistic, fine-grained "
        "center-based model (CBM) — explicit per-cell overdamped dynamics, a force-dependent E-cadherin "
        "catch-bond cohesion (Rakshit 2012), contact-inhibited proliferation, and an adhesive substrate — "
        "with ZERO calibration to the PI data, across the entire PI experimental R₀ range (104–394 µm, "
        "r² = 0.998, 35/35 GPU runs). The law’s MAGNITUDE, however, is ~5–9× below the PI medians. We show "
        "this magnitude gap is not a parameter to tune but the center-based one-particle structural limit, "
        "eliminating every CBM-expressible mechanism on six independent axes (scale, observable definition, "
        "cohesion, lateral collective coordination, passive wetting, active wetting). The aggregate-wetting "
        "framework (Douezan–Brochard-Wyart) explains it: the PI magnitude is the ACTIVE complete-wetting "
        "precursor-monolayer regime, an intrinsically shape-resolved process a point-cell model cannot host. "
        "Finally, the three PI ligand conditions (Bare/Pre/Lam4) separate and order correctly at the FORM "
        "level from mechanistic integrin-clutch kinetics. The magnitude is thus a precisely-located scope "
        "boundary that belongs to the fine-grained single-cell line.")

    # ---- 1. Model & methods ----
    _h(doc, "1. Model and methods", 1)
    doc.add_paragraph(
        "The CBM runs on the same HOOMD-blue + Leimkuhler–Matthews BAOAB-limit stack as the fine-grained "
        "single-cell line. Each cell is one overdamped particle, r(t+dt) = r(t) + (F/γ)·dt. The four cell "
        "scales are all measured/derived-anchored, never fitted to the PI spreading data:")
    for b in [
        "Cell diameter 15 µm (Wagner 2011, measured).",
        "Cell–cell cohesion force-anchored to the measured MCF7–MCF7 de-adhesion 6.5 nN (Iturri 2020), "
        "realized as the faithful Rakshit-2012 E-cadherin catch bond (N_cad = 223 = 6.5 nN / f₀); the "
        "force-dependent lifetime is adiabatically eliminated to an effective force-dependent cohesion "
        "(catch lifetime ~0.1 s ≪ CBM step).",
        "Per-cell migration drag γ = 0.30 N·s/m = clutch-ensemble friction (KU-2.18 Bangasser 2013), not "
        "water-Stokes.",
        "Contact-inhibited proliferation by the Drasdo–Höhme free-space division rule (MCF7 30 h cycle).",
    ]:
        doc.add_paragraph(b, style="List Bullet")
    doc.add_paragraph(
        "Production runs on an RTX A5000 GPU via the B2 GPU-main port (the per-step host sync removed; "
        "throughput ~N-flat). Spread area is reported both as a fragmentation-robust connected-core area "
        "and as a raw union-of-disks footprint (the image-segmentation analog of the PI readout).")

    # ---- 2. Headline ----
    _h(doc, "2. The PI law FORM is reproduced with zero calibration (headline)", 1)
    doc.add_paragraph(
        "Across the full PI experimental R₀ range, the fully mechanistic model reproduces the spreading-law "
        "FORM with no calibration:")
    p = doc.add_paragraph()
    rr = p.add_run("A/A₀ = 1.000 + (47.39 µm)/R + (97.46 µm²)/R²    r² = 0.998    (R₀ 104–394 µm, 35/35 GPU runs)")
    rr.bold = True
    doc.add_paragraph(
        "The signs match the PI law (b > 0 the 1/R surface/volume term; the asymptote a = 1.000 the "
        "cohesion-locked no-spread limit). The raw union-of-disks footprint equals the connected core to "
        "within 0.97–0.99 at every R₀, so the experiment’s own area definition does not close the gap — the "
        "magnitude difference is genuine physics, not an observable-definition artifact.")
    _fig(doc, "fig_layer2_prod_rlaw.png",
         "Figure 1. (A) Production R₀-law (35/35, real A5000): the mechanistic model (blue) reproduces the "
         "PI law FORM across the full PI R₀ range; PI medians (markers) sit ~5–9× higher. (B) raw/core ≈ "
         "0.97–0.99 — the gap is not an observable artifact.")

    # ---- 3. Six-axis elimination ----
    _h(doc, "3. The magnitude gap is the center-based structural limit (six-axis elimination)", 1)
    doc.add_paragraph(
        "The platform under-spreads ~5–9× in magnitude relative to the PI medians (7.2/7.5/10.0). We "
        "eliminated every CBM-expressible cause on six independent axes:")
    _fig(doc, "fig_layer2_morphology.png",
         "Figure 2. The spreading process rendered (Lam4, R₀≈151 µm, t = 0 → 61 h): the actual spatial "
         "reality behind A/A₀. TOP — the projected footprint grows only ~1.0→1.5×. BOTTOM — the aggregate "
         "stays a 3D CAP on the dish; it does NOT melt into a flat monolayer. This is the structural limit, "
         "made visible (an animated version is fig_layer2_morphology_anim.gif).")
    axes = [
        ("Scale (§B2)", "the native-N GPU sweep reaches the PI R₀ range; the gap persists at matched R₀ — "
                        "not a small-N/statistics artifact."),
        ("Observable (§B2)", "raw union-of-disks footprint ≈ connected core (0.97–0.99) — not the area "
                             "definition."),
        ("Cohesion (§C–§D)", "the lamellipodium→CBM motility bridge localizes the gap to the ACTIVE side "
                             "(cohesion is anchored two independent ways); the protrusion-scale force ejects "
                             "rim cells, and the Kadzik-grounded turnover-remodeled (ductile) cohesion removes "
                             "the ejection (flow not fracture) yet A/A₀ stays ~1.3 — the gap survives every "
                             "cohesion variant (catch/morse, brittle/yield)."),
        ("Lateral coordination (§E)", "an SPP-plithotaxis polarity field (Smeets-2016 MCF10A anchors: "
                                      "persistence D_r, free-edge CIL; correlations emergent from cohesion, no "
                                      "imposed Vicsek) keeps A/A₀ flat (~1.3) across the anchored force range — "
                                      "the persistent + CIL coordination is rim-limited at native N."),
        ("Passive wetting (§F)", "raising cell–substrate adhesion 16× does not flatten the cohesive cap "
                                 "(kinetic trap) — A/A₀ flat ~1.5."),
        ("Active wetting (§F)", "active motility + high adhesion lifts the raw footprint only modestly "
                                "(1.72→1.94, a partial precursor-film) while the core stays ~1.5 — no "
                                "complete-wetting monolayer."),
    ]
    for name, desc in axes:
        para = doc.add_paragraph(style="List Bullet")
        run = para.add_run(f"{name}: "); run.bold = True
        para.add_run(desc)
    doc.add_paragraph(
        "The aggregate-wetting framework (Douezan & Brochard-Wyart 2011; Beaune 2014; Gonzalez-Rodriguez "
        "2012) makes the limit precise. A spheroid on an adhesive substrate wets per the spreading "
        "coefficient S = W_cs − 2γ (cell–substrate adhesion vs cell–cell cohesion): S < 0 → a rounded 3D "
        "cap (our state); S > 0 → complete wetting as a motility-driven precursor MONOLAYER film (~1-cell "
        "height). Breast-carcinoma MCF7 spreading is ACTIVE (footprint 3–4× in 24 h on collagen-I; "
        "EpCAM-knockdown drives a flat coherent monolayer by 48 h — Aslemarz/Gupta 2024), pulled by motile "
        "edge cells, with cell-contact contractility (cohesion) as the dominant lever. A geometric "
        "flattening estimate A/A₀ ≈ 4R/3h gives ≈7–13 for R≈150 µm and h≈15–30 µm (1–2 cell layers) — the "
        "PI 7–10. The center-based point-cell model cannot represent this 3D-pile→2D-monolayer transition; "
        "it is intrinsically shape-resolved/subcellular.")
    _fig(doc, "fig_layer2_motility_bridge.png",
         "Figure 2. The active-traction force ladder (§C): the protrusion anchor (9.4 nN) exceeds cohesion "
         "(6.5 nN) → detachment, not spreading; the ×5.9 protrusion/whole-cell ratio ≈ the observed gap, "
         "localized to the active side (cohesion anchored two ways).")
    _fig(doc, "fig_layer2_daxis_yield.png",
         "Figure 3. Cohesion eliminated (§D): turnover-remodeled (ductile) cohesion removes the "
         "protrusion-force ejection (flow not fracture, Kadzik 2026) yet A/A₀ stays ~1.3–1.5 ≪ PI.")
    _fig(doc, "fig_layer2_plithotaxis.png",
         "Figure 4. Lateral coordination eliminated (§E): SPP-plithotaxis (persistent + free-edge CIL, "
         "MCF10A anchors) is flat at ~1.3 across the full anchored force range, overlapping the radial-crawl "
         "null — rim-limited at native N.")
    _fig(doc, "fig_layer2_wetting.png",
         "Figure 5. Wetting axis eliminated (§F): passive adhesion (no motility) cannot flatten the cap "
         "(kinetic trap); active motility + adhesion lifts the raw footprint only modestly (partial "
         "precursor-film) while the core stays ~1.5 ≪ PI 7–10.")

    # ---- 4. Ligand axis (c) ----
    _h(doc, "4. The ligand conditions separate and order correctly at the FORM level (c)", 1)
    doc.add_paragraph(
        "The three PI ligand conditions are driven mechanistically: each maps to an active edge-traction via "
        "per-species integrin-clutch kinetics (collagen-I clutch for Bare/Pre, the weaker laminin clutch for "
        "Lam4 at 0.61×; Bare<Pre is the flagged pV4D4 collagen-density axis), with Lam4 additionally given "
        "the A4′ partial-β1-uniformity (a larger traction screening length Lp≈40 µm engaging the interior — "
        "the collective resolution of the single-cell↔collective laminin split). Cohesion, substrate, and "
        "proliferation are held common so the only difference is the ligand-set traction. The magnitude is "
        "overlay-only (the known structural-limit scope boundary); the deliverable is the FORM-level "
        "separation and ordering of the three emergent A/A₀(R₀) curves.")
    lig = _ligand_table()
    pim = _pi_medians()
    if lig:
        t = doc.add_table(rows=1, cols=4); t.style = "Light Grid Accent 1"
        hdr = t.rows[0].cells
        hdr[0].text = "Condition"; hdr[1].text = "R₀ (µm) → A/A₀ core (mean±sd)"
        hdr[2].text = "f_traction"; hdr[3].text = "PI median (overlay)"
        fmap = {"Bare": "1.50 nN, Lp 11 µm", "Pre": "2.50 nN, Lp 11 µm", "Lam4": "1.53 nN, Lp 40 µm (A4′)"}
        for c in ("Bare", "Pre", "Lam4"):
            if c not in lig:
                continue
            row = t.add_row().cells
            row[0].text = c
            row[1].text = "; ".join(f"{r0}:{aa:.2f}±{sd:.2f}" for r0, aa, sd in lig[c])
            row[2].text = fmap.get(c, "")
            row[3].text = f"{pim.get(c, float('nan')):.1f}" if c in pim else "—"
        doc.add_paragraph(
            "The conditions separate (mechanistic clutch ordering) and the A4′ Lam4 partial-uniformity lifts "
            "Lam4’s large-size A/A₀ toward the PI Lam4>Pre>Bare collective ordering; the absolute magnitude "
            "remains the structural-limit scope boundary (PI medians overlaid).")
    else:
        doc.add_paragraph("[ligand-production results pending — table/figure inserted on run completion]")
    _fig(doc, "fig_layer2_ligand_production.png",
         "Figure 6. (c) Ligand-condition production: emergent A/A₀(R₀) for Bare/Pre/Lam4 from mechanistic "
         "clutch traction (+A4′ Lam4 partial-uniformity); PI medians overlay-only. FORM-level "
         "separation/ordering; magnitude = structural limit.")

    # ---- 5. Discussion / verdict ----
    _h(doc, "5. Discussion and verdict", 1)
    doc.add_paragraph(
        "The Layer-2 CBM delivers a clean two-part result. (i) POSITIVE: the PI spreading law’s FORM is an "
        "emergent property of a fully mechanistic model with zero calibration, across the entire PI R₀ range "
        "— a strong validation of the underlying single-cell mechanics carried into the multicellular "
        "regime. (ii) BOUNDING: the law’s MAGNITUDE (~5–9× under the PI medians) is the center-based "
        "one-particle structural limit, eliminated on six independent axes and explained precisely by the "
        "aggregate-wetting framework — the PI magnitude is the active complete-wetting precursor-monolayer "
        "regime, which a point-cell model cannot represent. The magnitude therefore belongs to the "
        "fine-grained single-cell line (a shape-resolved lamellipodium spreading on, and reacting traction "
        "against, the substrate); a CBM ‘wetting extension’ would require a lumped 3D→2D state-switch, which "
        "the project’s fine-grained-mechanistic principle disfavours. Throughout, the discipline was "
        "literature-first and overlay-only: every parameter is a measured/derived anchor, and the PI A/A₀ "
        "was never fitted.")

    # ---- 6. Anchors & integrity ----
    _h(doc, "6. Anchors and citation integrity", 1)
    doc.add_paragraph(
        "New literature anchors were sourced by adversarially-verified deep-research passes and recorded as "
        "SourceEvidence candidates (not auto-registered): the collective-migration cluster (Smeets 2016 "
        "MCF10A CIL-SPP; Petitjean 2010; Garcia 2015; Tambe 2011; Trepat 2009; Poujade 2007; Bi/Manning "
        "2016) and the aggregate-wetting cluster (Douezan 2011; Beaune 2014; Gonzalez-Rodriguez 2012; Nagle "
        "2022 breast-epithelial tissue surface tension 21–45 mN/m; Aslemarz/Gupta 2024 MCF7). Citation "
        "integrity was enforced (e.g. the tug-of-war buildup quote attributed to Tambe 2011, not Trepat "
        "2009; Ryan 2001 PMID corrected; no MCF7-specific W_cs exists — MCF10A/MDCK/S180 are flagged "
        "proxies). Open item: the §D2 single-cell→aggregate surface-tension bridge gave γ≈0.57 mN/m versus "
        "the Nagle breast-epithelial 21–45 mN/m (~40–80×) — to be reconciled by the surface-tension-bridge "
        "owner.")

    # ---- Appendix ----
    _h(doc, "Appendix — artifacts", 1)
    doc.add_paragraph(
        "Code: spheroid/{cbm,cadherin_bonds,proliferation,substrate,spreading,motility_bridge,"
        "substrate_crawl,plithotaxis,ligand_traction}.py; integrator/baoab_device.py (GPU port). Drivers: "
        "scripts/layer2_{prod_rlaw,gpu_scaleup,decisive_traction,daxis,plithotaxis,wetting,ligand_production}"
        ".py (+ auto-viz). Full technical record: outputs/layer2/REPORT.md (§HEADLINE, §B2, §C, §D, §E, §F). "
        "Commits on branch layer2/spheroid-cbm.")

    doc.save(str(_DOCX))
    return _DOCX


if __name__ == "__main__":
    out = build()
    print(f"wrote {out}  ({out.stat().st_size} bytes)")
