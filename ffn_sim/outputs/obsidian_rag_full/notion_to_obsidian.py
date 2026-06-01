#!/usr/bin/env python3
"""Export the ffn_cellsim Notion Contract-Graph (8 DBs) to an Obsidian vault.

Notion = source of truth. This produces a READ-ONLY mirror: ~240 markdown notes
(one per row) + YAML frontmatter + typed [[wikilinks]] for every relation, plus
a pre-baked `.obsidian/graph.json` so the graph opens already color-coded by
node type. Re-run any time Notion changes — the vault is regenerated from scratch.

SETUP (PI, one-time):
  1. https://www.notion.so/my-integrations -> New integration -> copy the
     "Internal Integration Secret" (ntn_... / secret_...).
  2. Open the "Contract Graph" Notion page -> ... -> Connections -> add the
     integration (the 8 child DBs inherit access).
  3. Put the token where this script reads it:
        echo 'NOTION_TOKEN=<paste>' > .notion_token       (next to this script)
     or export NOTION_TOKEN=<paste> in the shell.

RUN:
  conda activate ffn_sim
  python notion_to_obsidian.py
  open -a Obsidian ./vault

The 8 data-source IDs are the live collection IDs of the Contract Graph DBs.
"""
from __future__ import annotations
import json, os, pathlib, re, sys, time
import requests
try:
    from ku_dependencies import DEPS  # KU-x.y -> [KU-x.y, ...] dependency edges
except Exception:
    DEPS = {}
try:
    from ku_dependencies import DEPS_V2  # KB-x -> [KB-x, ...] for KU v2 layer
except Exception:
    DEPS_V2 = {}
try:
    from ku_dependencies import DEPS_PIV  # KB-PIV-x -> [KB-x, ...] PI-exp validation layer
    DEPS_V2 = {**DEPS_V2, **DEPS_PIV}     # merge: both are KB-keyed
except Exception:
    pass

HERE = pathlib.Path(__file__).parent
VAULT = HERE / "vault"
NOTION_VERSION = "2022-06-28"
API = "https://api.notion.com/v1"

# --- 8 DATABASE IDs (not collection ids) — the /databases/{id}/query endpoint --
DATA_SOURCES = {
    "SourceEvidence":  "bffea918b1c44581aa482f4f1ffc3d9d",
    "KnowledgeClaim":  "5b2d50409458479b8f774407c6e99719",
    "ModelContract":   "fc987bb7e4d64212928608b889177929",
    "Parameter":       "fd40391fcd034f6cac47a845634b3d66",
    "ValidationGate":  "277e1d5165e5403d959d5bf15ca4a317",
    "CodeMapping":     "e6df88ae27334a2280596941803eab2c",
    "RunResult":       "aa1a80911f6649edae956dab2bf96921",
    "DecisionLedger":  "793d56b19963444fbb8793ae6efe0191",
}

# graph.json color groups by node type (Obsidian color = 0xRRGGBB int)
TYPE_COLORS = {
    "SourceEvidence": 0x9e9e9e,  # gray  — papers (leaves)
    "KnowledgeClaim": 0x4caf50,  # green — literature claims
    "ModelContract":  0x2196f3,  # blue  — adopted approximations
    "Parameter":      0xffc107,  # amber — concrete values
    "ValidationGate": 0xf44336,  # red   — pass/fail gates
    "CodeMapping":    0x795548,  # brown — git pointers
    "RunResult":      0x607d8b,  # slate — run pointers
    "DecisionLedger": 0x9c27b0,  # purple— decisions
}


def get_token() -> str:
    tok = os.environ.get("NOTION_TOKEN", "").strip()
    if not tok:
        tf = HERE / ".notion_token"
        if tf.exists():
            txt = tf.read_text().strip()
            tok = txt.split("=", 1)[1].strip() if txt.startswith("NOTION_TOKEN=") else txt
    if not tok:
        sys.exit("NOTION_TOKEN not found (env or .notion_token). See setup in docstring.")
    return tok


def headers(tok):
    return {"Authorization": f"Bearer {tok}",
            "Notion-Version": NOTION_VERSION,
            "Content-Type": "application/json"}


def query_all(db_id, tok):
    """Return all rows of a database (handles pagination). Notion-Version 2022-06-28
    + /databases/{id}/query is the combination verified to work for this workspace."""
    rows, cursor = [], None
    url = f"{API}/databases/{db_id}/query"
    while True:
        body = {"page_size": 100}
        if cursor:
            body["start_cursor"] = cursor
        r = requests.post(url, headers=headers(tok), json=body, timeout=30)
        r.raise_for_status()
        data = r.json()
        rows.extend(data["results"])
        if not data.get("has_more"):
            break
        cursor = data["next_cursor"]
        time.sleep(0.2)
    return rows


def plain(prop):
    """Flatten one Notion property value to a string or list of related page-ids."""
    t = prop["type"]
    v = prop[t]
    if t == "title" or t == "rich_text":
        return "".join(x["plain_text"] for x in v)
    if t == "select":
        return v["name"] if v else ""
    if t == "multi_select":
        return ", ".join(x["name"] for x in v)
    if t == "url":
        return v or ""
    if t == "number":
        return "" if v is None else str(v)
    if t == "date":
        return (v or {}).get("start", "") if v else ""
    if t == "unique_id":
        return f"{v.get('prefix','')}{v.get('number','')}" if v else ""
    if t == "relation":
        return [x["id"].replace("-", "") for x in v]  # list of target page ids
    if t == "formula":
        return str(v.get(v["type"], ""))
    return ""


def slugify(s):
    s = re.sub(r"[\\/:*?\"<>|#^\[\]]", "", s).strip()
    return re.sub(r"\s+", " ", s)[:120] or "untitled"


def main():
    tok = get_token()
    if VAULT.exists():
        import shutil
        shutil.rmtree(VAULT)
    VAULT.mkdir(parents=True)

    # which property holds the stable ID per DB (so filenames = IDs, not prose)
    ID_PROP = {"KnowledgeClaim": "KB ID", "ModelContract": "MC ID",
               "Parameter": "PARAM ID", "ValidationGate": "VG ID",
               "CodeMapping": "CM ID", "RunResult": "RUN ID",
               "DecisionLedger": "DEC ID"}  # SourceEvidence: title IS the key

    # pass 1: pull every row, build id -> (filename, claim-title, type)
    raw = {}        # ds_name -> list of page dicts
    id2name = {}    # page-id(no dashes) -> note filename
    id2title = {}   # page-id -> human claim/title for the H1
    id2type = {}
    for ds, dsid in DATA_SOURCES.items():
        print(f"reading {ds} ...", flush=True)
        rows = query_all(dsid, tok)
        raw[ds] = rows
        for pg in rows:
            pid = pg["id"].replace("-", "")
            props = pg["properties"]
            title = next((plain(p) for p in props.values() if p["type"] == "title"), pid)
            # filename = explicit ID prop if present & non-empty, else title
            fname = title
            idp = ID_PROP.get(ds)
            if idp and idp in props:
                idv = plain(props[idp]).strip()
                if idv:
                    fname = idv
            id2name[pid] = slugify(fname)
            id2title[pid] = title
            id2type[pid] = ds
        print(f"  {len(rows)} rows")

    # disambiguate duplicate filenames
    seen = {}
    for pid, name in list(id2name.items()):
        if name in seen and seen[name] != pid:
            id2name[pid] = f"{name} ({id2type[pid][:3]}-{pid[:4]})"
        seen.setdefault(name, pid)

    # lookup for free-text cross-refs: every note filename, + KU-alias -> filename
    # (so "KB-3.15" or "KU-3.15" mentioned in any Notes/value becomes a wikilink)
    all_files = {id2name[p] for p in id2name}
    token2file = {f: f for f in all_files}            # KB-3.5 -> KB-3.5
    for pid, ds in id2type.items():
        if ds == "KnowledgeClaim":
            # map the KU alias to this claim's filename
            pg = next(r for r in raw[ds] if r["id"].replace("-", "") == pid)
            al = next((plain(p) for k, p in pg["properties"].items() if k == "Aliases"), "")
            for tok_ in re.findall(r'KU-[\dA-Z.]+', al):
                token2file.setdefault(tok_, id2name[pid])
    ref_re = re.compile(r'\b(KB-[\dA-Z.]+|KU-[\dA-Z.]+|MC-[\w.-]+|VG-[\w.-]+)\b')

    # pass 2: write notes
    n_notes = n_edges = 0
    for ds, rows in raw.items():
        for pg in rows:
            pid = pg["id"].replace("-", "")
            fm = {"type": ds, "notion_id": pid}
            body_props, link_lines = [], []
            for pname, prop in pg["properties"].items():
                val = plain(prop)
                if prop["type"] == "relation":
                    if val:
                        targets = " ".join(f"[[{id2name.get(t, t)}]]" for t in val if t in id2name)
                        if targets:
                            link_lines.append(f"{pname.lower().replace(' ','-')}:: {targets}")
                            n_edges += sum(1 for t in val if t in id2name)
                elif prop["type"] == "title":
                    continue
                elif val:
                    key = pname.lower().replace(" ", "_").replace("(", "").replace(")", "")
                    # never clobber reserved frontmatter keys (node type, id)
                    if key in ("type", "notion_id"):
                        key = "kind" if key == "type" else "src_id"
                    # short scalar -> frontmatter; long text -> body
                    if len(str(val)) <= 80 and "\n" not in str(val):
                        fm[key] = str(val).replace(":", " -")
                    else:
                        body_props.append(f"**{pname}:** {val}")
            fname = id2name[pid]
            h1 = id2title[pid]  # human claim sentence as the visible heading
            # KnowledgeClaim: add KU->KU dependency wikilinks (claim<->claim web)
            if ds == "KnowledgeClaim":
                my_alias = fm.get("aliases", "")
                ku = my_alias.split()[0] if my_alias else ""  # e.g. "KU-3.5"
                deps = [("KB-" + d.split("-", 1)[1]) for d in DEPS.get(ku, [])]  # KU-3.4->KB-3.4
                deps += DEPS_V2.get(fname, [])  # KU v2 layer: already KB-keyed
                dep_links = []
                for target in dict.fromkeys(deps):
                    if target != fname and ((VAULT / f"{target}.md").exists() or target in all_files):
                        dep_links.append(f"[[{target}]]")
                if dep_links:
                    link_lines.append("depends-on:: " + " ".join(dep_links))
                    n_edges += len(dep_links)
            # free-text cross-refs: scan all text values for KB/KU/MC/VG tokens,
            # link to their note (rescues orphan VGs whose Notes cite "KB-3.15", etc.)
            already = set(re.findall(r'\[\[([^\]]+)\]\]', " ".join(link_lines)))
            text_blob = " ".join(str(plain(p)) for p in pg["properties"].values()
                                  if p["type"] in ("rich_text", "title"))
            rel = []
            for tok_ in dict.fromkeys(ref_re.findall(text_blob)):  # ordered-unique
                tgt = token2file.get(tok_)
                if tgt and tgt != fname and tgt not in already and f"[[{tgt}]]" not in rel:
                    rel.append(f"[[{tgt}]]")
            if rel:
                link_lines.append("relates:: " + " ".join(rel))
                n_edges += len(rel)
            # external links: Notion page + source DOI/URL (clickable in Obsidian)
            ext = [f"🔗 [Open in Notion](https://www.notion.so/{pid})"]
            doi = ""
            for pn, pr in pg["properties"].items():
                if pr["type"] == "url" and pr["url"]:
                    u = pr["url"]
                    label = "DOI / source" if "doi" in u.lower() else pn
                    ext.append(f"📄 [{label}]({u})")
                    doi = u
            if doi:
                fm["doi"] = doi
            out = ["---"] + [f"{k}: {v}" for k, v in fm.items()] + ["---", "", f"# {h1}", ""]
            out += ext + [""]
            out += body_props + [""] + link_lines
            (VAULT / f"{fname}.md").write_text("\n".join(out) + "\n", encoding="utf-8")
            n_notes += 1

    # pass 3: .obsidian/graph.json with color groups by type
    obs = VAULT / ".obsidian"
    obs.mkdir(exist_ok=True)
    # Obsidian color-group query = its search syntax; frontmatter match is ["key":"value"]
    groups = [{"query": f'["type":"{t}"]', "color": {"a": 1, "rgb": rgb}}
              for t, rgb in TYPE_COLORS.items()]
    graph = {
        "collapse-filter": True, "search": "", "showTags": False,
        "showAttachments": False, "hideUnresolved": False, "showOrphans": True,
        "collapse-color-groups": False, "colorGroups": groups,
        "collapse-display": False, "showArrow": True, "textFadeMultiplier": 0,
        "nodeSizeMultiplier": 1.1, "lineSizeMultiplier": 1,
        "collapse-forces": False, "centerStrength": 0.3, "repelStrength": 12,
        "linkStrength": 1, "linkDistance": 90, "scale": 1,
    }
    (obs / "graph.json").write_text(json.dumps(graph, indent=2), encoding="utf-8")
    # app.json: enable frontmatter, show in graph
    (obs / "app.json").write_text(json.dumps(
        {"propertyAttachment": True, "showFrontmatter": True}, indent=2), encoding="utf-8")

    print(f"\nwrote {n_notes} notes, {n_edges} relation edges to {VAULT}")
    print("color groups baked into .obsidian/graph.json")
    print(f"\nopen -a Obsidian \"{VAULT}\"")


if __name__ == "__main__":
    main()
