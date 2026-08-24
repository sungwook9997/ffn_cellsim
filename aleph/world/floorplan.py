"""The floorplan: what is claimed, by whom, and on whose authority — as a page you can watch.

WHY THIS AND NOT A 3D VIEW.  Ten million nodes cannot be drawn individually at interactive rates, and
the thing worth watching while a cell is being assembled is not shape anyway — it is WHAT IS WIRED. That
payload is kilobytes, it is host-side, it touches no device array, and it is therefore safe to emit
between accepted steps at any cadence. The spatial views come later and are bounded by the viewport
rather than by the population; this one is bounded by the number of populations, which is a dozen.

THE VIEW THIS IS MODELLED ON.  A chip floorplan, and specifically the screen a designer actually keeps
open: not the layout alone but the layout with the design-rule violations drawn ON it. The equivalent
here is that every bond family is shown together with the SOURCE CLASS of the count that produced it,
so a population built on an unsourced density is visible as such rather than buried in a JSON field.
Under the PI's 2026-08-15 ruling that every physiological value is physiological FOR SOME PARTICULAR
CELL, "which parts of this cell rest on numbers nobody has sourced" is a first-class question, and this
is the view that answers it at a glance.

WHAT IT DELIBERATELY DOES NOT DO.  It reads. It computes no physics, launches no kernel, mutates no
arena, and takes no device. It also does not sample: every number on the page is a total over the whole
population, because a viewer that quietly showed a subset would misreport exactly the thing it exists
to report.

Sanity Gate (recorded per CLAUDE.md, before first execution):
  * dimensional — counts are dimensionless; bytes are bytes; a density carries the unit its basis
    declares and is printed with it, never bare.
  * boundary — an arena with nothing claimed renders an empty die rather than raising, because "nothing
    is wired yet" is the state this view exists to show at the start of a build.
  * conservation/invariant — the page's per-population totals are read from the arena's own census, so
    they cannot disagree with it; the utilisation bar is live/capacity and cannot exceed 1 by
    construction, since a claim past capacity raises.
  * CFL/precision — not applicable; nothing is integrated.
  * sign sense — not applicable; nothing is accumulated.
  * measurement protocol — one host-side read of host-side bookkeeping. No readback, no device, and
    therefore nothing that could sit inside a physical-time loop even by accident.

Runtime: pure host. Emits a self-contained HTML page with no external asset of any kind.
"""

from __future__ import annotations

import html
import json
from collections.abc import Sequence

from aleph.world.arena import Kind, WorldArena
from aleph.world.bond import BondFamily

__all__ = ["render_floorplan", "write_floorplan"]

#: How each source class is drawn. A count nobody sourced must not look like one that was sourced —
#: this mapping is the whole "violations on the layout" idea, and it is why the class travels with the
#: family rather than living in a separate report.
_CLASS_STYLE: dict[str, tuple[str, str]] = {
    "SOURCED": ("ok", "traceable to a literature anchor"),
    "DERIVED": ("ok", "computed from a sourced anchor or from geometry"),
    "CONVENIENCE": ("warn", "a numerical guard, explicitly NOT a physiological magnitude"),
    "UNRATIFIED_PROXY": ("bad", "sourced, but for a different cell / state / assay than this run"),
    "PI_GAP": ("bad", "the physiological value is unknown and the code says so"),
}

_CSS = """
:root{--bg:#0B1020;--panel:#151C2F;--line:#2C3654;--ink:#E6ECFA;--dim:#8494B0;
--ok:#5BC8AF;--warn:#E8C547;--bad:#FF7A59;--accent:#7C9BFF}
@media(prefers-color-scheme:light){:root{--bg:#F7F6F2;--panel:#FFF;--line:#DAD5C9;--ink:#181B22;
--dim:#6C7486;--ok:#0E8F76;--warn:#9A7407;--bad:#C2452A;--accent:#3A5BD9}}
*{box-sizing:border-box}
body{margin:0;background:var(--bg);color:var(--ink);font:14px/1.55 ui-sans-serif,system-ui,sans-serif}
.wrap{max-width:1000px;margin:0 auto;padding:32px 24px 64px;display:flex;flex-direction:column;gap:26px}
h1{font-size:24px;margin:0;letter-spacing:-.02em}
h2{font-size:13px;margin:0;letter-spacing:.11em;text-transform:uppercase;color:var(--dim);font-weight:650}
.sub{color:var(--dim);font-size:13px}
.panel{background:var(--panel);border:1px solid var(--line);border-radius:9px;padding:18px 20px;
display:flex;flex-direction:column;gap:13px}
.mono{font-family:ui-monospace,SFMono-Regular,Menlo,monospace;font-variant-numeric:tabular-nums}
.bar{height:22px;background:var(--line);border-radius:4px;overflow:hidden;position:relative}
.bar>i{display:block;height:100%;background:var(--accent);border-radius:4px 0 0 4px}
.row{display:grid;grid-template-columns:96px 1fr 190px;gap:12px;align-items:center}
.die{display:flex;gap:3px;height:64px;border:1px solid var(--line);border-radius:6px;overflow:hidden}
.die>div{display:flex;align-items:center;justify-content:center;font-size:11px;overflow:hidden;
white-space:nowrap;color:#0B1020;font-weight:650}
table{border-collapse:collapse;width:100%;font-size:13px}
th,td{padding:8px 10px;text-align:left;border-bottom:1px solid var(--line);vertical-align:top}
th{font-size:11px;letter-spacing:.08em;text-transform:uppercase;color:var(--dim);white-space:nowrap}
tr:last-child td{border-bottom:0}
.tag{display:inline-block;padding:2px 8px;border-radius:999px;font-size:11px;font-weight:650;
border:1px solid currentColor}
.ok{color:var(--ok)}.warn{color:var(--warn)}.bad{color:var(--bad)}
.empty{color:var(--dim);font-style:italic}
.note{border-left:3px solid var(--bad);padding:10px 14px;background:#FF7A5915;border-radius:0 6px 6px 0}
"""

#: Deterministic population colours. Fixed and cycled rather than hashed, so the same build renders the
#: same way twice — a die map whose colours moved between refreshes would be unreadable as a progress view.
_HUES = ["#5BC8AF", "#7C9BFF", "#E8C547", "#FF7A59", "#B98CFF", "#59C2E8", "#E8899B", "#9BD45B"]


def _esc(x: object) -> str:
    return html.escape(str(x), quote=True)


def render_floorplan(
    arena: WorldArena,
    families: Sequence[BondFamily] = (),
    *,
    title: str = "Arena floorplan",
    standalone: bool = True,
) -> str:
    """Render the arena's census and its bond families as HTML.

    Args:
        arena: the world to read. Not mutated, and no device array is touched.
        families: the bond families built into it, in build order. Their counts' source classes are
            what the violation panel is drawn from.
        title: page title.
        standalone: ``True`` emits a complete document, which is what a local viewer needs. ``False``
            emits the style and content only, for a host that supplies its own document skeleton — the
            page is identical either way, and neither form references an external asset.

    Returns:
        The HTML, with no external asset of any kind.
    """
    c = arena.census()
    live: dict[str, int] = c["live"]          # type: ignore[assignment]
    cap: dict[str, int] = c["capacity"]       # type: ignore[assignment]
    pops: dict[str, dict[str, int]] = c["populations"]  # type: ignore[assignment]

    # ── the die: populations by NODE count, proportional, in claim order ──────────────────────────
    node_pops = [(name, kinds.get("node", 0)) for name, kinds in pops.items() if kinds.get("node", 0)]
    total_nodes = sum(n for _, n in node_pops)
    # Proportional among the LIVE populations, NOT against capacity. A die drawn against a reservation
    # is a sliver at low utilisation and shows nothing — and utilisation is already its own panel, so
    # drawing it twice costs the one view that has to stay readable at every stage of a build.
    if node_pops:
        die = "".join(
            f'<div style="flex:{n};background:{_HUES[i % len(_HUES)]}" '
            f'title="{_esc(name)}: {n:,} nodes ({n / total_nodes:.1%} of live)">'
            f'{_esc(name) if n / total_nodes > 0.06 else ""}</div>'
            for i, (name, n) in enumerate(node_pops)
        )
    else:
        die = '<div style="flex:1;background:var(--line)"></div>'

    # ── utilisation, one row per reserved kind ───────────────────────────────────────────────────
    rows = []
    for kind in Kind:
        k = str(kind)
        if k not in cap:
            continue
        n, total = live.get(k, 0), cap[k]
        pct = (n / total * 100.0) if total else 0.0
        rows.append(
            f'<div class="row"><span class="mono">{_esc(k)}</span>'
            f'<span class="bar"><i style="width:{pct:.3f}%"></i></span>'
            f'<span class="mono sub">{n:,} / {total:,} &nbsp;({pct:.2f}%)</span></div>'
        )

    # ── bond families, with the source class of the count that produced each ─────────────────────
    if families:
        frows = []
        for f in families:
            row = f.provenance_row()
            cls = str(row["source_class"])
            style, meaning = _CLASS_STYLE.get(cls, ("warn", "unrecognised source class"))
            pairs = ", ".join(
                f"{a}–{b}: {n:,}" for (a, b), n in sorted(f.component_pairs(arena).items())
            ) or "<span class='empty'>none</span>"
            density = (f"{row['density']:g} per {row['support_unit']}"
                       if row["basis"] != "explicit" else "explicit")
            frows.append(
                f"<tr><td><strong>{_esc(row['family'])}</strong><br>"
                f"<span class='sub mono'>{_esc(row['chemistry_card'])}</span></td>"
                f"<td class='mono'>{row['n_bonds']:,}</td>"
                f"<td class='mono'>{_esc(density)}<br>"
                f"<span class='sub'>&times; {row['support']:g} {_esc(row['support_unit'])}</span></td>"
                f"<td>{pairs}</td>"
                f"<td><span class='tag {style}'>{_esc(cls)}</span><br>"
                f"<span class='sub'>{_esc(row['scope'])}</span><br>"
                f"<span class='sub'>{_esc(row['provenance'])}</span></td></tr>"
            )
        fam_html = (
            "<table><thead><tr><th>family</th><th>bonds</th><th>count basis</th>"
            "<th>joins (derived)</th><th>authority</th></tr></thead><tbody>"
            + "".join(frows) + "</tbody></table>"
        )
    else:
        fam_html = "<p class='empty'>No bond family built yet — the die is allocated but nothing is wired.</p>"

    # ── the violation panel: which parts of this cell rest on numbers nobody sourced ──────────────
    unsourced = [f for f in families
                 if _CLASS_STYLE.get(str(f.count.source_class), ("", ""))[0] in ("warn", "bad")]
    if unsourced:
        items = "".join(
            f"<li><strong>{_esc(f.name)}</strong> — {f.n_bonds:,} bonds on a "
            f"<span class='{_CLASS_STYLE[str(f.count.source_class)][0]}'>"
            f"{_esc(f.count.source_class)}</span> count: {_esc(f.count.provenance)}</li>"
            for f in unsourced
        )
        note = (f"<div class='note'><strong>{len(unsourced)} of {len(families)} families rest on a count "
                f"that is not a sourced measurement for this cell.</strong><ul>{items}</ul></div>")
    elif families:
        note = "<div class='panel'><span class='ok'>Every family's count is SOURCED or DERIVED.</span></div>"
    else:
        note = ""

    payload = _esc(json.dumps({"census": c, "families": [f.provenance_row() for f in families]}))
    head = (f"""<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>{_esc(title)}</title><style>{_CSS}</style></head><body>"""
            if standalone else f"<title>{_esc(title)}</title><style>{_CSS}</style>")
    tail = "</div></body></html>" if standalone else "</div>"
    return f"""{head}<div class="wrap">
<div>
  <h1>{_esc(title)}</h1>
  <p class="sub mono">{total_nodes:,} nodes live of {cap.get('node', 0):,} reserved &middot;
  {len(node_pops)} node population(s) &middot; {len(families)} bond famil(ies) &middot;
  {c['n_claims']} claim(s) &middot; device {_esc(c['device'] or 'none — bookkeeping only')}</p>
</div>
<div class="panel"><h2>Die — populations by node count, proportional among live</h2>
  <div class="die">{die}</div>
  <p class="sub">Reservation is shown in <em>Capacity</em> below; this bar is the live prefix only,
  so it stays readable whether the arena is 2% claimed or 100%.</p></div>
<div class="panel"><h2>Capacity</h2>{''.join(rows)}</div>
<div class="panel"><h2>Bond families</h2>{fam_html}</div>
{note}
<script type="application/json" id="floorplan-data">{payload}</script>
{tail}"""


def write_floorplan(arena: WorldArena, path: str, families: Sequence[BondFamily] = (), **kw) -> str:
    """Render and write the floorplan to ``path``; returns the path written."""
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(render_floorplan(arena, families, **kw))
    return path


def _demo() -> None:
    """Self-check: a build with one sourced and one unsourced family renders both, and says so."""
    import numpy as np

    from aleph.world.bond import BondCount, SourceClass, build_bond_family
    from aleph.world.strand import build_strand

    arena = WorldArena(capacity={Kind.NODE: 2_000, Kind.SEGMENT: 2_000,
                                 Kind.ANGLE3: 2_000, Kind.BOND: 2_000})

    empty = render_floorplan(arena)
    assert "nothing is wired" in empty, "an empty die must say so rather than raise"

    strands = [build_strand(arena, "cortex", start=(0, 0.05 * i, 0), direction=(1, 0, 0),
                            contour_um=1.0, seg_um=0.1) for i in range(8)]
    build_strand(arena, "membrane", start=(0, 0, 1.0), direction=(1, 0, 0),
                 contour_um=1.0, seg_um=0.05)

    pairs = np.stack([
        np.concatenate([s.nodes.lo + np.arange(s.n_nodes) for s in strands[:-1]]),
        np.concatenate([s.nodes.lo + np.arange(s.n_nodes) for s in strands[1:]]),
    ], axis=1)

    fams = [
        build_bond_family(
            arena, "alpha_actinin", chemistry_card="alpha_actinin_ferrer2008",
            count=BondCount(basis="explicit", value=40.0, scope="MCF7 interphase adherent 37C",
                            source_class=SourceClass.SOURCED,
                            provenance="Ferrer 2008 PNAS AFM, 455 pN/nm; PI-approved 2026-06-30"),
            support=1.0, pairs=pairs, rest_um=0.05, stiffness_pn_per_um=4.6e5),
        build_bond_family(
            arena, "erm", chemistry_card="ezrin",
            count=BondCount(basis="areal", value=58.0, scope="none — subdiv-6 mesh, not a cell",
                            source_class=SourceClass.UNRATIFIED_PROXY,
                            provenance="subdiv-6 icosphere vertex density; NOT a physiological datum"),
            support=0.2, pairs=pairs, rest_um=0.02, stiffness_pn_per_um=4.6e3),
    ]

    page = render_floorplan(arena, fams, title="Arena floorplan — self-check")
    for must in ("cortex", "membrane", "alpha_actinin", "erm",
                 "UNRATIFIED_PROXY", "Ferrer 2008", "1 of 2 families"):
        assert must in page, f"missing from the page: {must!r}"
    assert "SOURCED" in page and "NOT a physiological datum" in page
    assert page.startswith("<!doctype html>") and page.rstrip().endswith("</html>")
    assert "http://" not in page and "https://" not in page, "the page must be self-contained"

    frag = render_floorplan(arena, fams, title="fragment", standalone=False)
    assert not frag.startswith("<!doctype") and "<body>" not in frag and "</html>" not in frag
    assert "cortex" in frag and "UNRATIFIED_PROXY" in frag, "the fragment is the same page, unwrapped"
    print(f"floorplan self-check OK — {len(page):,} B, "
          f"{arena.n_live(Kind.NODE):,} nodes, {len(fams)} families, 1 flagged")


if __name__ == "__main__":
    _demo()
