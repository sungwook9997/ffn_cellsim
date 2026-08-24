"""Stage-render emit point for the staged ac/engine build (P0 resting / P1 SF, and any future slice).

This is the viz-side entry point the physics stages call to EMIT a render when a slice reaches CUDA_UNIT / GO.
It does NOT touch physics: the stage (on the gbook A5000) writes a v2 state dump with
``scripts/dump_state.py`` and then invokes this — as a subprocess or an import — with the dump path + a stage
label. Here (dev Mac, numpy only) it renders the ONE self-contained interactive WebGL HTML, browser-checks a
curated set of scenes for that stage (proving a real render — not a blank), extracts flagship turntable GIFs for
the gh-pages gallery, and writes a gallery card. Rendering is decoupled from CUDA so it also runs on the Mac
against a dump synced from gbook.

Stage recipes (curated scene filters — everything else in the dump is still in the viewer dropdown):

* **P0 — resting cell + residual heatmap:** the composed whole cell (cortex coloured by the composed per-node
  |F| — that IS the residual/force heatmap) + the steric-interpenetration hotspots. Answers "is the resting
  baseline force field where we expect, and where does it blow up?"
* **P1 — SF isolation + FA/cortex connector:** the ``sf_arc`` component ALONE + the ``transient-actin`` and
  ``MOTOR`` connectors between ``sf_arc`` and ``cortex`` (the FA/cortex load path), so the SF slice is readable
  before it is allowed a connected/production claim.

    # gbook stage, after emitting a dump (CUDA):
    python -m aleph.scripts.dump_state --incumbent --out .../p0_resting_v2.npz         # or --engine (P1)
    # then (here / synced Mac):
    python -m aleph.scripts.ac_stage_render --stage P0 --npz .../p0_resting_v2.npz
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from aleph.scripts.ac_cell_assembled_viz import build as build_html
from aleph.scripts.ac_viz_common import build_scenes, load_dump

_ROOT = Path(__file__).resolve().parents[1]
GALLERY = _ROOT / "outputs" / "ac" / "stage_gallery"

# curated scene-name SUBSTRINGS per stage (matched case-insensitively against the dump's actual scene keys)
STAGE_SCENES = {
    "P0": ["composed cell", "steric interpenetration", "◦ cortex only", "compose"],
    "P1": ["◦ sf_arc only", "transient-actin · sf_arc", "MOTOR · sf_arc", "compose"],
    "ALL": ["composed cell", "◦ ", "⇄ ", "▸ compose"],
}


def _match_scenes(all_scenes: list[str], filters: list[str]) -> list[str]:
    out = []
    for s in all_scenes:
        if any(f.lower() in s.lower() for f in filters):
            out.append(s)
    return out or all_scenes[:1]


def render_stage(npz: Path, stage: str, *, out_dir: Path | None = None, make_gifs: bool = True,
                 budget_ms: int = 20000) -> dict:
    """Render + browser-check + (optional) GIF a stage dump; write a gallery card. Returns the card dict."""
    out_dir = out_dir or (GALLERY / stage.lower())
    out_dir.mkdir(parents=True, exist_ok=True)
    dump = load_dump(npz)
    scenes, _ = build_scenes(dump)
    scene_names = list(scenes.keys())
    filters = STAGE_SCENES.get(stage.upper(), STAGE_SCENES["ALL"])
    chosen = _match_scenes(scene_names, filters)

    html = out_dir / f"ac_stage_{stage.lower()}.html"
    _, counts = build_html(Path(npz), html)

    # browser-check every curated scene → a real-render proof (blank/JS-error detection is in browser_check)
    from aleph.scripts.browser_check import check
    shots, failed = [], []
    for s in chosen:
        tag = "".join(ch if ch.isalnum() else "_" for ch in s)[:40]
        png = out_dir / f"shot_{tag}.png"
        rc = check(html, png, scene=s, budget_ms=budget_ms, timeout_s=150.0)
        (shots if rc == 0 else failed).append(s)

    gifs = []
    if make_gifs:
        try:
            from aleph.scripts.flagship_gif import extract
            gpaths = extract(html, out_dir / "gifs", prefix=f"ac_{stage.lower()}", mode="turntable",
                             scenes_filter=[f for f in filters if f not in ("compose", "▸ compose")],
                             n_frames=36, width=720, duration_ms=60, ready_timeout=180.0,
                             window=(1000, 1000), still=False, zoom=1.15, hide_ui=True)
            gifs = [str(p.relative_to(out_dir)) for p in gpaths]
        except Exception as e:  # GIFs are a gallery nicety, never block the stage render
            print(f"[ac_stage_render] GIF extraction skipped: {e}")

    stage_meta = dump.meta.get("stage", {})
    card = {
        "stage": stage, "stage_name": stage_meta.get("name", stage),
        "evidence": stage_meta.get("evidence", ""), "npz": str(npz), "html": str(html),
        "components": [{"key": c.key, "nodes": c.n_nodes, "actors": c.n_actors} for c in dump.components],
        "connectors": [{"family": cn.family, "a": cn.a_component, "b": cn.b_component, "n": cn.n}
                       for cn in dump.connectors],
        "scenes_total": len(scene_names), "scenes_checked": shots, "scenes_failed": failed,
        "gifs": gifs, "counts": counts,
    }
    (out_dir / "card.json").write_text(json.dumps(card, indent=2, default=float))
    _write_gallery_index()
    status = "OK" if not failed else f"{len(failed)} scene(s) FAILED render"
    print(f"[ac_stage_render] {stage}: {len(shots)}/{len(chosen)} scenes rendered ({status}); "
          f"{len(gifs)} GIFs; card → {out_dir/'card.json'}")
    return card


def _write_gallery_index() -> None:
    """Assemble a single gh-pages-ready index over every stage card written so far (data-driven)."""
    cards = []
    for cj in sorted(GALLERY.glob("*/card.json")):
        try:
            cards.append((cj.parent.name, json.loads(cj.read_text())))
        except Exception:
            pass
    rows = []
    for sub, c in cards:
        comps = ", ".join(f"{x['key']}({x['actors']:,})" for x in c["components"])
        cons = ", ".join(f"{x['family']}:{x['n']}" for x in c["connectors"]) or "—"
        gif_imgs = "".join(
            f'<img src="{sub}/{g}" loading="lazy" style="max-width:320px;margin:4px;border-radius:6px">'
            for g in c.get("gifs", []))
        rows.append(
            f'<div class="card"><h2>{c["stage"]} · {c["stage_name"]} '
            f'<span class="ev">{c.get("evidence","")}</span></h2>'
            f'<p class="meta">components: {comps}</p><p class="meta">connectors: {cons}</p>'
            f'<p class="meta">{c["scenes_total"]} scenes · {len(c.get("scenes_checked",[]))} render-verified'
            f'{" · <b>"+str(len(c["scenes_failed"]))+" FAILED</b>" if c.get("scenes_failed") else ""}</p>'
            f'<div class="gifs">{gif_imgs}</div>'
            f'<p><a href="{sub}/{Path(c["html"]).name}">open interactive viewer →</a></p></div>')
    html = (
        '<!doctype html><meta charset="utf-8"><title>AC engine — staged-build gallery</title>'
        '<style>body{background:#0e1117;color:#c9d1d9;font-family:system-ui;margin:0;padding:24px}'
        'h1{font-weight:650}.card{background:#161b22;border:1px solid #30363d;border-radius:10px;'
        'padding:14px 18px;margin:14px 0}.meta{color:#8b949e;font-size:13px;margin:3px 0}'
        '.ev{color:#37d67a;font-size:12px;border:1px solid #37d67a55;border-radius:5px;padding:1px 6px}'
        'a{color:#7fd4ff}.gifs{display:flex;flex-wrap:wrap}</style>'
        '<h1>AC engine — staged component/connector build</h1>'
        '<p class="meta">Each card = one staged slice, rendered per-component (isolation) + per-connector '
        '(load path) + cumulative composition. Data-driven from the dump manifest — new components/connectors '
        'appear automatically.</p>' + "".join(rows))
    (GALLERY / "index.html").write_text(html)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--stage", required=True, help="P0 | P1 | ALL | <label>")
    ap.add_argument("--npz", required=True, type=Path, help="v2 state dump for this stage")
    ap.add_argument("--out-dir", type=Path, default=None)
    ap.add_argument("--no-gifs", action="store_true", help="skip flagship GIF extraction")
    ap.add_argument("--budget-ms", type=int, default=20000)
    a = ap.parse_args()
    render_stage(a.npz, a.stage, out_dir=a.out_dir, make_gifs=not a.no_gifs, budget_ms=a.budget_ms)


if __name__ == "__main__":
    main()
