#!/usr/bin/env python3
"""#2 — add CODE + TEST + DOC nodes to the Obsidian vault, linked to the KB
claims they reference. Run AFTER notion_to_obsidian.py (it appends to the vault).

Scans ffn_sim/ for .py modules, tests, and docs/*.md. For each file:
 - creates a note  <Code> <relpath>.md  (type: Code / Test / Doc)
 - parses KU-x.y / KB-x.y tokens in the file -> links `implements:: [[KB-x.y]]`
   (KU-5.x in code = lamellipodium = H5-LP-x where applicable; see ALIAS map)
This makes the graph span literature -> claim -> contract -> CODE -> test.

Pointer, not mirror: the note stores the repo path + the KB refs found, not the
code body. git stays the source of truth.
"""
from __future__ import annotations
import pathlib, re

ROOT = pathlib.Path(__file__).resolve().parents[2]   # .../ffn_cellsim/ffn_sim
VAULT = pathlib.Path(__file__).parent / "vault"
REPO = ROOT.parent                                    # .../ffn_cellsim

KU = re.compile(r'KU-(\d+)\.(\d+[A-Za-z]?\.?\d*)')    # KU-3.5, KU-1.V.3.1, KU-3.B2.1
KB_TOK = re.compile(r'KB-[\dA-Z.]+')

# code's KU-5.x means lamellipodium (H5-LP), NOT collective KB-5.x — map by file area
def ku_to_kb(unit, rest, path):
    tok = f"KU-{unit}.{rest}"
    if unit == "5" and ("lamell" in path or "h5" in path or "lp" in path.lower()):
        # code lamellipodium KU-5.1/2/3 -> H5-LP-1/2/3
        m = re.match(r'(\d+)', rest)
        if m:
            return f"H5-LP-{m.group(1)}"
    return "KB-" + f"{unit}.{rest}"

SCAN = [
    ("Code", ROOT, "*.py", ("scripts", "outputs", "__pycache__", "tests")),
    ("Test", ROOT / "tests", "*.py", ("__pycache__",)),
    ("Doc",  ROOT / "docs", "*.md", ()),
]

def existing_kb():
    return {p.stem for p in VAULT.glob("*.md")}

def main():
    if not VAULT.exists():
        raise SystemExit("run notion_to_obsidian.py first")
    kb_files = existing_kb()
    n_made = n_link = 0
    for ntype, base, pat, skip in SCAN:
        if not base.exists():
            continue
        for f in base.rglob(pat):
            rel = f.relative_to(REPO).as_posix()
            if any(s in rel for s in skip):
                continue
            try:
                txt = f.read_text(encoding="utf-8", errors="ignore")
            except Exception:
                continue
            refs = set()
            for u, rest in KU.findall(txt):
                tgt = ku_to_kb(u, rest.rstrip("."), rel)
                if tgt in kb_files:
                    refs.add(tgt)
            for tok in KB_TOK.findall(txt):
                if tok in kb_files:
                    refs.add(tok)
            if not refs:
                continue  # only add code nodes that actually touch the knowledge graph
            stem = ntype[0] + "_" + f.stem      # C_myosin, T_test_cortex, D_H3_cortex
            # disambiguate
            base_stem = stem
            i = 2
            while (VAULT / f"{stem}.md").exists():
                stem = f"{base_stem}_{i}"; i += 1
            links = " ".join(f"[[{r}]]" for r in sorted(refs))
            out = ["---", f"type: {ntype}", f"path: {rel}",
                   f"kb_refs: {len(refs)}", "---", "",
                   f"# {rel}", "",
                   f"\U0001f4c1 `{rel}`  ({len(refs)} KB refs)", "",
                   f"implements:: {links}"]
            (VAULT / f"{stem}.md").write_text("\n".join(out) + "\n", encoding="utf-8")
            n_made += 1; n_link += len(refs)
    print(f"added {n_made} code/test/doc nodes, {n_link} links to KB claims")

    # extend graph.json color groups with the 3 new types
    import json
    gj = VAULT / ".obsidian" / "graph.json"
    g = json.loads(gj.read_text())
    new = {"Code": 0x00bcd4, "Test": 0xcddc39, "Doc": 0x8d6e63}
    have = {grp["query"] for grp in g["colorGroups"]}
    for t, rgb in new.items():
        q = f'["type":"{t}"]'
        if q not in have:
            g["colorGroups"].append({"query": q, "color": {"a": 1, "rgb": rgb}})
    gj.write_text(json.dumps(g, indent=2), encoding="utf-8")
    print("graph.json color groups extended: Code=cyan Test=lime Doc=brown")

if __name__ == "__main__":
    main()
