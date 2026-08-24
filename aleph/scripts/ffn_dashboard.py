#!/usr/bin/env python
r"""One page answering "what is going on right now" — lease, lanes, evidence, runs, doc debt.

    python aleph/scripts/ffn_dashboard.py            # -> aleph/outputs/dashboard.html
    python aleph/scripts/ffn_dashboard.py --print    # same figures to stdout

Every number here already existed somewhere; none of it was ever in one place, so the questions a
session actually asks at minute zero — *is the GPU free? whose lane is this file in? has this been
measured already? which of these documents pre-date the reframe?* — each cost a search. That is
also why duplicated work kept happening: the cheapest way to find out whether something had been
done was to do it.

READS ONLY, and says so when it cannot read. A missing lease file means the GPU is free; a missing
run index means nothing has run through the enforced path yet; neither is an error. What it will not
do is infer: an absent input is rendered as absent, never as zero, because "no record" and "measured
zero" are the distinction this project has spent the most time repairing.

Sanity Gate:
    * boundary: every input is optional; each panel degrades to "not available" independently.
    * conservation/invariant: tier-(a) rows are counted from the append-only YAML, not scraped from
      the rendered table, so the dashboard cannot disagree with the renderer.
    * measurement-protocol: the run index shows each run's POPULATION string beside it — the whole
      point of the 2026-07-28 retraction was that scale was invisible at a glance.
    * numerical / dimensional / sign-sense: not applicable — no quantity is computed.
"""

from __future__ import annotations

import argparse
import html
import importlib.util
import json
import subprocess
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from aleph.coordination.gpu_lease import RunIndex, read_lease  # noqa: E402
from aleph.coordination.ownership import load_ownership  # noqa: E402

REPO = Path(__file__).resolve().parents[2]


def _load_contracts():
    """Load ``ac/engine/contracts.py`` WITHOUT importing the engine package.

    Importing ``aleph.engine`` runs its eager ``__init__``, which pulls in Warp — and a
    read-only status page has no business initialising a GPU runtime. ``contracts.py`` is pure
    dataclasses and enums, so it loads standalone.

    Returns:
        The module, or ``None`` if it cannot be loaded.
    """
    path = REPO / "aleph" / "engine" / "contracts.py"
    try:
        spec = importlib.util.spec_from_file_location("_ffn_contracts", path)
        module = importlib.util.module_from_spec(spec)
        # Register BEFORE executing: `dataclasses` resolves annotations through
        # `sys.modules[cls.__module__]`, so a module absent from the table raises while building the
        # first dataclass. This is the documented way to load a module standalone, not a workaround.
        sys.modules[spec.name] = module
        spec.loader.exec_module(module)
        return module
    except Exception:
        return None


def collect() -> dict:
    """Gather every panel's data, each independently optional."""
    data: dict = {"generated_at": time.time(), "repo": str(REPO)}

    lease = read_lease()
    data["lease"] = (
        {"held": True, "holder": lease.holder, "host": lease.host, "reason": lease.reason,
         "minutes_left": round(lease.seconds_remaining() / 60.0, 1)}
        if lease is not None and lease.is_live() else {"held": False}
    )

    index = RunIndex()
    data["runs"] = {
        "summary": index.summary(),
        "recent": [
            {"when": time.strftime("%m-%d %H:%M", time.localtime(r.started_at)),
             "build": r.build_commit[:8], "holder": r.holder,
             "population": r.population or "(not declared)", "label": r.run_label}
            for r in list(index)[-15:]
        ],
    }

    try:
        own = load_ownership(REPO / "ownership.yaml")
        data["lanes"] = {name: len(globs) for name, globs in own.sessions.items()}
        data["shared_files"] = sorted(own.shared)
    except Exception:
        data["lanes"], data["shared_files"] = None, []

    contracts = _load_contracts()
    if contracts is not None:
        data["ladder"] = [r.value for r in contracts.EVIDENCE_ORDER]
        data["claims"] = [c.value for c in contracts.QuantitativeClaim]
    else:
        data["ladder"], data["claims"] = None, None

    rows_path = REPO / "aleph" / "docs" / "state_rows.yaml"
    if rows_path.exists():
        import yaml
        rows = (yaml.safe_load(rows_path.read_text()) or {}).get("rows") or []
        data["tier_a"] = {
            "n": len(rows),
            "populations": [r.get("population", "")[:60] for r in rows],
        }
    else:
        data["tier_a"] = None

    audit = REPO / "aleph" / "docs" / "v2_audit"
    if audit.exists():
        docs = sorted(audit.glob("*.md"))
        pre = 0
        for d in docs:
            got = subprocess.run(["git", "-C", str(REPO), "log", "-1",
                                  "--format=%ad", "--date=format:%Y-%m-%d", "--", str(d)],
                                 capture_output=True, text=True, check=False).stdout.strip()
            if got and got < "2026-07-25":
                pre += 1
        data["docs"] = {"root_md": len(docs), "pre_reframe": pre}
    else:
        data["docs"] = None
    return data


def _panel(title: str, body: str, note: str = "") -> str:
    tail = f'<p class="note">{html.escape(note)}</p>' if note else ""
    return f'<section><h2>{html.escape(title)}</h2>{body}{tail}</section>'


def to_html(d: dict) -> str:
    """Render the collected data as one self-contained page."""
    esc = html.escape
    lease = d["lease"]
    if lease["held"]:
        lease_body = (f'<p class="big held">HELD by {esc(lease["holder"])}</p>'
                      f'<p>{esc(lease["reason"] or "(reason unstated)")} — '
                      f'{lease["minutes_left"]} min left, on {esc(lease["host"])}</p>')
    else:
        lease_body = '<p class="big free">FREE</p><p>nothing holds the device</p>'

    runs = d["runs"]["summary"]
    rows = "".join(
        f'<tr><td>{esc(r["when"])}</td><td><code>{esc(r["build"])}</code></td>'
        f'<td>{esc(r["holder"])}</td><td>{esc(r["population"])}</td><td>{esc(r["label"])}</td></tr>'
        for r in d["runs"]["recent"])
    runs_body = (
        f'<p>{runs["n_runs"]} run(s) · {runs["n_distinct_configs"]} distinct config(s) · '
        f'{runs["n_builds"]} build(s)</p>'
        + (f'<table><tr><th>when</th><th>build</th><th>session</th><th>population</th>'
           f'<th>driver</th></tr>{rows}</table>'
           if rows else '<p class="note">nothing has run through the enforced path yet</p>'))

    if d["lanes"]:
        lanes_body = "".join(f'<li><b>{esc(k)}</b> — {v} glob(s)</li>' for k, v in d["lanes"].items())
        lanes_body = f"<ul>{lanes_body}</ul>"
        lanes_body += ('<p class="note">shared, owned by nobody: '
                       + ", ".join(f"<code>{esc(s)}</code>" for s in d["shared_files"]) + "</p>")
    else:
        lanes_body = '<p class="note">ownership declaration not readable</p>'

    if d["ladder"]:
        ladder_body = ('<p>' + " → ".join(f'<code>{esc(r)}</code>' for r in d["ladder"]) + '</p>'
                       + '<p>× ' + " / ".join(f'<code>{esc(c)}</code>' for c in d["claims"]) + '</p>')
    else:
        ladder_body = '<p class="note">contracts module not readable</p>'

    if d["tier_a"]:
        pops = "".join(f"<li>{esc(p)}</li>" for p in d["tier_a"]["populations"])
        tier_body = (f'<p class="big">{d["tier_a"]["n"]}</p><p>quotable result(s), each with its '
                     f'population:</p><ul class="small">{pops}</ul>')
    else:
        tier_body = '<p class="note">state_rows.yaml not present</p>'

    if d["docs"]:
        pct = 100.0 * d["docs"]["pre_reframe"] / max(1, d["docs"]["root_md"])
        docs_body = (f'<p class="big">{d["docs"]["pre_reframe"]} / {d["docs"]["root_md"]}</p>'
                     f'<p>v2_audit root documents that pre-date the 2026-07-25 reframe '
                     f'({pct:.0f}%) — written for a different objective and not yet marked.</p>')
    else:
        docs_body = '<p class="note">v2_audit not present</p>'

    stamp = time.strftime("%Y-%m-%d %H:%M", time.localtime(d["generated_at"]))
    return f"""<!doctype html><meta charset="utf-8"><title>ffn_cellsim — now</title>
<style>
 :root {{ color-scheme: light dark; }}
 body {{ font: 14px/1.55 -apple-system, system-ui, sans-serif; margin: 0 auto; padding: 2rem;
        max-width: 1100px; }}
 h1 {{ font-size: 1.4rem; margin: 0 0 .2rem; }} h2 {{ font-size: .95rem; margin: 0 0 .6rem;
        text-transform: uppercase; letter-spacing: .06em; opacity: .65; }}
 section {{ border: 1px solid color-mix(in srgb, currentColor 18%, transparent); border-radius: 10px;
        padding: 1rem 1.2rem; margin: .8rem 0; }}
 .grid {{ display: grid; gap: .8rem; grid-template-columns: repeat(auto-fit, minmax(300px, 1fr)); }}
 .big {{ font-size: 1.8rem; font-weight: 650; margin: .1rem 0; }}
 .free {{ color: #2e7d32; }} .held {{ color: #c1502e; }}
 .note {{ opacity: .62; font-size: .85rem; }} .small li {{ font-size: .82rem; opacity: .8; }}
 table {{ border-collapse: collapse; width: 100%; font-size: .82rem; }}
 th, td {{ text-align: left; padding: .25rem .5rem; border-bottom: 1px solid
        color-mix(in srgb, currentColor 12%, transparent); }}
 th {{ opacity: .6; font-weight: 600; }} code {{ font-size: .85em; }}
 ul {{ margin: .3rem 0; padding-left: 1.1rem; }}
</style>
<h1>ffn_cellsim — what is going on right now</h1>
<p class="note">generated {stamp} · read-only · an absent input renders as absent, never as zero</p>
<div class="grid">
{_panel("GPU lease", lease_body, "one A5000; the lease self-expires so a dead session cannot wedge it")}
{_panel("Tier-(a) quotable results", tier_body, "rendered from the append-only state_rows.yaml")}
</div>
{_panel("Runs that actually executed", runs_body,
        "the enforced path records (build, config); a repeat of both is refused unless forced")}
<div class="grid">
{_panel("Session lanes", lanes_body, "isolation by declared ownership — one tree, no sibling checkouts")}
{_panel("Evidence axes", ladder_body, "structural rung × whether a magnitude may be quoted (D1-B)")}
</div>
{_panel("Documentation debt", docs_body,
        "a document written before the objective changed reads as current unless it is marked")}
"""


def main() -> int:
    """Collect and render."""
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--out", default="aleph/outputs/dashboard.html")
    ap.add_argument("--print", dest="show", action="store_true", help="figures to stdout instead")
    args = ap.parse_args()

    data = collect()
    if args.show:
        print(json.dumps({k: v for k, v in data.items() if k != "runs"}, indent=2, default=str))
        print(json.dumps(data["runs"]["summary"], indent=2))
        return 0
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(to_html(data))
    print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
