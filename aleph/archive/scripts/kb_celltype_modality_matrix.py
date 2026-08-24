"""Build the cell-type x measurement-method matrix from the KB, so a validation target can be chosen.

WHY THIS EXISTS.  The ratified library is K=2 (MCF7, MDA-MB-231), and on 2026-07-29 a session read that
as "the data we have is two cell lines" and ranked virtual measurement protocols on a single pair's
separability.  The KB holds **22 distinct cell types/lines** carrying numeric claims.  Collapsing to two
is not a smaller analysis, it is a different and weaker one: a model can be fitted to two lines, while a
panel that spans four tissue origins, a four-stage melanoma progression, and a non-adherent
(cortex-only, no FA/clutch) lineage constrains the SAME mechanistic parameters from directions that
cannot all be satisfied by tuning.

WHAT IT IS AND IS NOT.  This is an INDEX over `knowledge_claim`, extracted by pattern from free-text
`value_range_si` / `title` fields.  It does not re-derive anything and it is not a data table: every row
carries its `kb_id` so the claim itself remains the source, and the `source_audit` verdict of the
citation is carried beside it because `CLAUDE.md` requires that check before a source is cited in a
deliverable.  A cell appearing here means "the KB says something numeric about this pair", never "this
value is usable".

Read-only.  Writes nothing but its own report.
"""

from __future__ import annotations

import argparse
import re
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
KB = ROOT / "aleph" / "outputs" / "tag_kb" / "kb.duckdb"

#: Cell line / type -> the pattern that identifies it in KB free text, with its lineage.  Lineage is the
#: axis that makes the panel a constraint rather than a list: a parameter law must transfer across it.
CELL_TYPES: dict[str, tuple[str, str]] = {
    "MCF10A":       (r"MCF-?10A", "breast · normal epithelial"),
    "MCF7":         (r"MCF-?7", "breast · non-invasive carcinoma"),
    "MDA-MB-231":   (r"MDA(?:-MB)?-?231", "breast · metastatic"),
    "IOSE364":      (r"IOSE", "ovary · normal surface epithelial"),
    "SKOV3":        (r"SKOV-?3", "ovary · carcinoma"),
    "WM35":         (r"WM35", "melanoma · RGP (stage 1/4)"),
    "ME10538":      (r"ME10538", "melanoma · VGP (stage 2/4)"),
    "A375P":        (r"A375P", "melanoma · metastatic (stage 3/4)"),
    "A375M":        (r"A375M", "melanoma · metastatic+ (stage 4/4)"),
    "HL60":         (r"HL-?60", "blood · AML — NON-ADHERENT, cortex-only"),
    "Jurkat":       (r"Jurkat", "blood · ALL — NON-ADHERENT, cortex-only"),
    "K562":         (r"K-?562", "blood · CML — NON-ADHERENT, cortex-only"),
    "neutrophil":   (r"neutrophil", "blood · primary — NON-ADHERENT"),
    "PC-3":         (r"PC-?3", "prostate · metastatic"),
    "PANC-1":       (r"PANC-?1", "pancreas · carcinoma"),
    "HeLa":         (r"HeLa", "cervix · carcinoma"),
    "MDCK":         (r"MDCK", "kidney · epithelial (canine)"),
    "fibroblast":   (r"fibroblast", "connective · primary"),
    "hMSC":         (r"hMSC|mesenchymal stem", "stem · mesenchymal"),
    "RBC":          (r"\bRBC\b|erythrocyte", "blood · anucleate ARCHETYPE"),
    "blastomere":   (r"blastomere", "embryo · ARCHETYPE"),
    "C. elegans":   (r"C\.? ?elegans", "model organism"),
}

#: Measurement method -> pattern.  The point of the project's cross-validation is that these are NOT
#: interchangeable: the KB records a 10-100x spread between them on the same quantity (`KB-1.14`), and
#: `KB-6.1.6` classifies three "cortical tension" numbers as DISTINCT QUANTITIES by protocol.
METHODS: dict[str, str] = {
    "AFM (sharp-tip)":      r"sharp-?tip",
    "AFM (colloidal/bead)": r"colloidal|10\s*um bead|spherical AFM",
    "AFM (parallel-plate)": r"parallel-?plate",
    "AFM (unspecified)":    r"\bAFM\b",
    "micropipette":         r"micropipette|aspiration",
    "MRS (nanowire)":       r"\bMRS\b|magnetic nanowire",
    "TFM (traction)":       r"\bTFM\b|traction force microscop",
    "SICM":                 r"\bSICM\b",
    "microwell":            r"microwell",
    "tether pulling":       r"tether",
    "cone-plate rheology":  r"cone-?plate|rheolog",
    "optical/stretcher":    r"optical stretcher|optical tweez",
    "ektacytometry":        r"ektacytometry",
}


def main() -> int:
    """Emit the matrix; return 2 if the KB cannot be opened, else 0."""
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=Path, default=None, help="write the markdown report here")
    args = ap.parse_args()

    try:
        import duckdb
    except ImportError:
        print("duckdb is not installed in this interpreter")
        return 2
    if not KB.is_file():
        print(f"KB not found at {KB}")
        return 2

    con = duckdb.connect(str(KB), read_only=True)
    claims = con.execute(
        "select kb_id, title, value_range_si, status, coalesce(citations_text,'') "
        "from knowledge_claim where value_range_si is not null"
    ).fetchall()
    audit = {k: v for k, v in con.execute(
        "select citation_key, verdict from source_audit").fetchall()}

    # cell -> method -> [(kb_id, status, observable)]
    grid: dict[str, dict[str, list[tuple[str, str, str]]]] = defaultdict(lambda: defaultdict(list))
    per_cell_claims: dict[str, set[str]] = defaultdict(set)
    method_only: dict[str, set[str]] = defaultdict(set)

    for kb_id, title, value, status, cites in claims:
        blob = f"{title} {value} {cites}"
        cells = [n for n, (pat, _) in CELL_TYPES.items() if re.search(pat, blob, re.I)]
        meths = [m for m, pat in METHODS.items() if re.search(pat, blob, re.I)]
        # "AFM (unspecified)" only when no more specific AFM variant matched
        if len([m for m in meths if m.startswith("AFM")]) > 1:
            meths = [m for m in meths if m != "AFM (unspecified)"]
        for cell in cells:
            per_cell_claims[cell].add(kb_id)
            for m in meths or ["(method not stated)"]:
                grid[cell][m].append((kb_id, status, title[:58]))
        if not cells:
            for m in meths:
                method_only[m].add(kb_id)

    lines: list[str] = []
    w = lines.append
    w("# Cell-type x measurement-method matrix — extracted from `kb.duckdb`\n")
    w(f"{len(claims)} numeric `knowledge_claim` rows scanned · "
      f"**{len(per_cell_claims)} cell types** matched · {len(METHODS)} methods probed\n")
    w("A cell means the KB says something NUMERIC about that (cell, method) pair. It does not mean the "
      "value is usable: check the claim's own status and its citation's `source_audit` verdict.\n")

    order = sorted(per_cell_claims, key=lambda c: (-len(per_cell_claims[c]), c))
    w("\n## Coverage per cell type\n")
    w("| cell type | lineage | claims | methods named |")
    w("|---|---|---:|---|")
    for cell in order:
        ms = [m for m in grid[cell] if m != "(method not stated)"]
        w(f"| **{cell}** | {CELL_TYPES[cell][1]} | {len(per_cell_claims[cell])} | "
          f"{', '.join(sorted(ms)) if ms else '—'} |")

    w("\n## The matrix\n")
    used = [m for m in METHODS if any(m in grid[c] for c in order)] + ["(method not stated)"]
    w("| cell type | " + " | ".join(used) + " |")
    w("|---" * (len(used) + 1) + "|")
    for cell in order:
        row = [f"**{cell}**"]
        for m in used:
            hits = grid[cell].get(m, [])
            row.append(" ".join(f"`{k}`" for k, _, _ in hits) if hits else "·")
        w("| " + " | ".join(row) + " |")

    w("\n## Methods that appear with NO cell type attached\n")
    w("These are generic/model claims — a method the KB discusses without binding it to a line.\n")
    for m, ids in sorted(method_only.items()):
        w(f"- **{m}** — {', '.join(sorted(ids)[:8])}{' …' if len(ids) > 8 else ''}")

    w("\n## Citation audit of the sources behind these claims\n")
    verdicts: dict[str, int] = defaultdict(int)
    for v in audit.values():
        verdicts[v] += 1
    w(f"`source_audit` over {len(audit)} citation keys: "
      + ", ".join(f"**{k}** {n}" for k, n in sorted(verdicts.items(), key=lambda kv: -kv[1])))
    w("\n`CLAUDE.md`: confirm a source's verdict is `OK` before citing it in a deliverable. "
      "`STATE.md` (c) 11 records that `CHECK` is ~95% a year-parsing regex artifact rather than a "
      "citation-quality statement, but that is a reason to look, not a reason to skip the look.\n")

    text = "\n".join(lines)
    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(text, encoding="utf-8")
        print(f"wrote {args.out}")
    else:
        print(text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
