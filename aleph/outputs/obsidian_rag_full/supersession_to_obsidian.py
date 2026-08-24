#!/usr/bin/env python3
"""#4 — surface supersession / authoritative-chain edges in the Obsidian vault.

Run LAST in refresh.sh (after notion_to_obsidian.py wipes+rebuilds the vault,
add_code_nodes.py adds the D_<doc> nodes, and devlogs_to_obsidian.py adds the
day-logs). This step reads the SAME explicit supersession source as the TAG
layer — the `kb_record:` front-matter in `ffn_sim/docs/**/*.md`, parsed by
`tag_kb/supersession.py` (single source of truth, no local inference) — and:

 - annotates each record's vault note with `supersedes::` / `superseded-by::`
   wikilinks + a `record_status` / `authoritative_as_of` frontmatter field, and
   `concerns::` links to the KB claim / VG gate it anchors;
 - creates a small `R_<id>` stub note for any record that has no doc note (e.g.
   a fileless commit-message claim cited only as a supersession target);
 - adds graph.json color groups so authoritative (gold) vs superseded (red)
   records are visually distinct.

The point: in the graph, the latest authoritative record visibly points back at
the findings it overturned, so a reader never mistakes a superseded hypothesis
for the current conclusion.
"""
from __future__ import annotations
import json
import pathlib
import re
import sys

from notion_to_obsidian import VAULT, slugify

# reuse the TAG layer's explicit-source parser (DRY single source of truth)
TAG_KB = pathlib.Path(__file__).resolve().parent.parent / "tag_kb"
sys.path.insert(0, str(TAG_KB))
from supersession import collect_records, resolve_authoritative  # noqa: E402

AUTH_COLOR = 0xffd700   # gold  — current authoritative record
SUPERSEDED_COLOR = 0xb71c1c  # dark red — superseded / overturned record


def node_for(record_id: str) -> tuple[str, bool]:
    """(vault note stem, exists_as_doc_node) for a record id.

    A doc record's note is D_<stem> (add_code_nodes.py naming); fall back to a
    disambiguated D_<stem>_N; otherwise it's a fileless record → R_<slug>."""
    d = f"D_{record_id}"
    if (VAULT / f"{d}.md").exists():
        return d, True
    cand = sorted(VAULT.glob(f"D_{record_id}_*.md"))
    if cand:
        return cand[0].stem, True
    return f"R_{slugify(record_id).replace(' ', '_')}", False


def _inject(stem: str, status: str, as_of: str | None,
            supersedes: list[str], superseded_by: list[str],
            anchors: list[str]) -> None:
    """Add record_status frontmatter + supersession wikilinks to an existing note."""
    path = VAULT / f"{stem}.md"
    text = path.read_text(encoding="utf-8")
    if "record_status:" in text:
        return  # already annotated this refresh — keep idempotent on standalone re-runs
    lines = text.splitlines()
    # locate the frontmatter close (2nd '---')
    fence = [i for i, ln in enumerate(lines) if ln.strip() == "---"]
    if len(fence) >= 2:
        ins = fence[1]
        fm = [f"record_status: {status}"]
        if as_of:
            fm.append(f"authoritative_as_of: {as_of}")
        lines[ins:ins] = fm
    tail: list[str] = [""]
    if supersedes:
        tail.append("supersedes:: " + " ".join(f"[[{s}]]" for s in supersedes))
    if superseded_by:
        tail.append("superseded-by:: " + " ".join(f"[[{s}]]" for s in superseded_by))
    if as_of:
        tail.append(f"authoritative-as-of:: {as_of}")
    if anchors:
        tail.append("concerns:: " + " ".join(f"[[{a}]]" for a in anchors))
    path.write_text("\n".join(lines + tail) + "\n", encoding="utf-8")


def _make_stub(stem: str, record_id: str, topic: str, status: str,
               as_of: str | None, conclusion: str | None,
               supersedes: list[str], superseded_by: list[str],
               anchors: list[str]) -> None:
    fm = ["---", "type: Record", f"record_id: {record_id}", f"topic: {topic}",
          f"record_status: {status}"]
    if as_of:
        fm.append(f"authoritative_as_of: {as_of}")
    fm.append("---")
    body = ["", f"# {record_id}", "",
            conclusion or "(supersession target — no standalone doc note)", ""]
    if supersedes:
        body.append("supersedes:: " + " ".join(f"[[{s}]]" for s in supersedes))
    if superseded_by:
        body.append("superseded-by:: " + " ".join(f"[[{s}]]" for s in superseded_by))
    if anchors:
        body.append("concerns:: " + " ".join(f"[[{a}]]" for a in anchors))
    (VAULT / f"{stem}.md").write_text("\n".join(fm + body) + "\n", encoding="utf-8")


def main():
    if not VAULT.exists():
        raise SystemExit("run notion_to_obsidian.py first (vault missing)")

    records, edges = collect_records()      # offline: docs front-matter + stubs
    if not records:
        print("no kb_record supersession metadata found — nothing to mirror")
        return

    by_id = {r["record_id"]: r for r in records}
    # per-record supersedes / superseded_by from the directional edges
    sup_of: dict[str, list[str]] = {}
    supby_of: dict[str, list[str]] = {}
    for e in edges:
        sup_of.setdefault(e["newer_id"], []).append(e["older_id"])
        supby_of.setdefault(e["older_id"], []).append(e["newer_id"])

    stems = {rid: node_for(rid) for rid in by_id}
    n_ann = n_stub = n_edge = 0
    annotated_topics = set()

    for rid, r in by_id.items():
        stem, is_doc = stems[rid]
        # anchors: link the record to its KB claim / VG gate node when present
        anchors = []
        for a in (r.get("claim"), r.get("gate")):
            if a and (VAULT / f"{a}.md").exists():
                anchors.append(a)
        sup = [stems[o][0] for o in sup_of.get(rid, []) if o in stems]
        supby = [stems[n][0] for n in supby_of.get(rid, []) if n in stems]
        n_edge += len(sup) + len(supby)
        if is_doc:
            _inject(stem, r["status"], r["authoritative_as_of"], sup, supby, anchors)
            n_ann += 1
        else:
            _make_stub(stem, rid, r["topic"], r["status"],
                       r["authoritative_as_of"], r.get("conclusion"),
                       sup, supby, anchors)
            n_stub += 1
        annotated_topics.add(r["topic"])

    # graph.json color groups: authoritative=gold, superseded=red
    gj = VAULT / ".obsidian" / "graph.json"
    if gj.exists():
        g = json.loads(gj.read_text())
        have = {grp["query"] for grp in g["colorGroups"]}
        for q, rgb in (('["record_status":"authoritative"]', AUTH_COLOR),
                       ('["record_status":"superseded"]', SUPERSEDED_COLOR)):
            if q not in have:
                g["colorGroups"].append({"query": q, "color": {"a": 1, "rgb": rgb}})
        gj.write_text(json.dumps(g, indent=2), encoding="utf-8")

    print(f"supersession mirror: annotated {n_ann} doc note(s), "
          f"{n_stub} stub record node(s), {n_edge} supersession wikilink(s) "
          f"across {len(annotated_topics)} topic(s)")
    for t in sorted(annotated_topics):
        auth = resolve_authoritative(records, t)
        if auth:
            print(f"  [{t}] authoritative -> {node_for(auth['record_id'])[0]}")
    print("graph.json color groups: authoritative=gold, superseded=red")


if __name__ == "__main__":
    main()
