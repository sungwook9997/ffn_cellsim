#!/usr/bin/env python3
"""Build a browsable HTML gallery of all ffn_cellsim output figures + movies.

Lives at outputs/deck_assets/GALLERY.html. Lets the PI eyeball every test figure,
read its filename/path, and drag the real file from Finder into Claude Design.
"""
import os
from urllib.parse import quote

OUT = "/Users/sw1/ffn_cellsim/ffn_sim/outputs"
DA = os.path.join(OUT, "deck_assets")
THUMBS = os.path.join(DA, "_thumbs")

def flat(rel):
    return rel.replace("/", "__").replace(" ", "__")

# collect assets relative to OUT, skipping deck_assets
pngs, mp4s, gifs = [], [], []
for root, dirs, files in os.walk(OUT):
    if os.path.abspath(root).startswith(os.path.abspath(DA)):
        dirs[:] = []
        continue
    for f in files:
        rel = os.path.relpath(os.path.join(root, f), OUT)
        ext = f.lower().rsplit(".", 1)[-1] if "." in f else ""
        if ext == "png":
            pngs.append(rel)
        elif ext == "mp4":
            mp4s.append(rel)
        elif ext == "gif":
            gifs.append(rel)

# group everything by directory
groups = {}
def add(rel, kind):
    d = os.path.dirname(rel)
    groups.setdefault(d, []).append((rel, kind))

for r in pngs: add(r, "png")
for r in mp4s: add(r, "mp4")
for r in gifs: add(r, "gif")

def thumb_src(rel, kind):
    if kind == "png":
        t = os.path.join(THUMBS, flat(rel))
        return "_thumbs/" + quote(flat(rel)) if os.path.exists(t) else "../" + quote(rel)
    if kind == "mp4":
        name = flat(rel)[:-4] + "__MOV.png"
        t = os.path.join(THUMBS, name)
        return "_thumbs/" + quote(name) if os.path.exists(t) else ""
    if kind == "gif":
        return "../" + quote(rel)  # gifs are small, reference original
    return ""

def card(rel, kind):
    fname = os.path.basename(rel)
    href = "../" + quote(rel)
    ts = thumb_src(rel, kind)
    abspath = os.path.join(OUT, rel)
    badge = {"mp4": '<span class="badge mov">MOVIE</span>',
             "gif": '<span class="badge gif">GIF</span>'}.get(kind, "")
    img = (f'<img loading="lazy" src="{ts}" alt="{fname}">' if ts
           else '<div class="noimg">no preview</div>')
    return (f'<div class="card">'
            f'<a href="{href}" target="_blank">{img}{badge}</a>'
            f'<div class="fn">{fname}</div>'
            f'<div class="pp">{abspath}</div>'
            f'</div>')

# order: put the headline dirs first
PRIORITY = ["warp_decohesion/figs", "h_dcm_two_stage/figs", "layer2/figs",
            "h_dcm_gpu_lod/figs", "h3/figs", "h7/figs"]
def sort_key(d):
    return (PRIORITY.index(d) if d in PRIORITY else len(PRIORITY) + 1, d)

total = len(pngs) + len(mp4s) + len(gifs)
parts = [f"""<!doctype html><html><head><meta charset="utf-8">
<title>ffn_cellsim figure gallery</title>
<style>
:root{{--bg:#0E1116;--card:#161b22;--teal:#34D2C8;--coral:#FF6B5C;--dim:#8b949e}}
body{{background:var(--bg);color:#e6edf3;font-family:Inter,system-ui,Arial,sans-serif;margin:0;padding:24px}}
h1{{font-size:20px;margin:0 0 4px}} .sub{{color:var(--dim);font-size:13px;margin-bottom:18px}}
h2{{font-size:15px;color:var(--teal);border-bottom:1px solid #21262d;padding:18px 0 6px;margin:24px 0 10px;position:sticky;top:0;background:var(--bg)}}
.count{{color:var(--dim);font-weight:400;font-size:12px}}
.grid{{display:grid;grid-template-columns:repeat(auto-fill,minmax(220px,1fr));gap:14px}}
.card{{background:var(--card);border:1px solid #21262d;border-radius:8px;overflow:hidden}}
.card a{{display:block;position:relative;line-height:0}}
.card img{{width:100%;height:150px;object-fit:contain;background:#0b0e12}}
.noimg{{height:150px;display:flex;align-items:center;justify-content:center;color:var(--dim);font-size:12px}}
.badge{{position:absolute;top:6px;right:6px;font-size:10px;font-weight:700;padding:2px 6px;border-radius:4px;color:#0b0e12}}
.badge.mov{{background:var(--coral)}} .badge.gif{{background:var(--teal)}}
.fn{{font:600 12px/1.3 'JetBrains Mono',ui-monospace,monospace;padding:8px 8px 2px;word-break:break-all}}
.pp{{font:11px/1.3 ui-monospace,monospace;color:var(--dim);padding:0 8px 8px;word-break:break-all}}
.toc{{font-size:12px;color:var(--dim);margin-bottom:8px}} .toc a{{color:var(--teal);text-decoration:none;margin-right:12px}}
</style></head><body>
<h1>ffn_cellsim — figure &amp; movie gallery</h1>
<div class="sub">{total} assets ({len(pngs)} png · {len(mp4s)} movies · {len(gifs)} gif). Click a tile to open full-res. To use in Claude Design: find the file in Finder at the gray path and drag it in. Movie tiles show a poster frame — the matching GIF/poster is staged in <b>deck_assets/curated/</b>.</div>
<div class="toc">"""]

ordered = sorted(groups.keys(), key=sort_key)
for d in ordered:
    anchor = d.replace("/", "_") or "root"
    parts.append(f'<a href="#{anchor}">{d or "(root)"} ({len(groups[d])})</a>')
parts.append("</div>")

for d in ordered:
    anchor = d.replace("/", "_") or "root"
    items = sorted(groups[d], key=lambda x: x[0])
    parts.append(f'<h2 id="{anchor}">{d or "(root)"} <span class="count">— {len(items)}</span></h2><div class="grid">')
    parts.extend(card(rel, kind) for rel, kind in items)
    parts.append("</div>")

parts.append("</body></html>")

with open(os.path.join(DA, "GALLERY.html"), "w") as fh:
    fh.write("\n".join(parts))
print(f"  GALLERY.html: {total} assets, {len(ordered)} folders")
