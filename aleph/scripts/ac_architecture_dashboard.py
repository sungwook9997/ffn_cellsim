#!/usr/bin/env python3
"""One-page map of the `ac/engine` cell: the declared component/connector graph vs what has actually run.

Every other viewer in this repo shows ONE compartment or ONE experiment. The thing that was
missing is the project itself: 14 components / 36 connectors, which of them own state, which
connectors actually carry force, and — the part that is easy to lose — **how far each one has
really climbed**. That last axis is the reason this script exists, and the reason it does not
invent anything:

  * the GRAPH is read from code (:func:`aleph.engine.contracts.reference_cell_architecture`),
    so it cannot drift from the contract;
  * the EVIDENCE is read from the state dumps on disk, each of which declares its own
    ``stage.evidence`` (e.g. ``CONNECTED`` for the GATE-B cortex-motor run, ``SCHEMA-PROOF``
    for the synthetic topology graph). A component is credited with the *highest* evidence
    level of any dump it actually appears in, and the dump is named next to it.

Nothing here is hand-authored status. If a component has only ever appeared in a
``SCHEMA-PROOF`` dump, this page says so — that is the honest signal ("declared, not yet run"),
and it is exactly what a per-component roadmap prose paragraph tends to blur.

Evidence vocabulary lives in :mod:`aleph.engine.contracts`, not here: the 8-state ladder
(`CONTRACTED → SEAMED → KERNEL_BOUND → CUDA_UNIT → CONNECTED → NATIVE → OPTIMISED → PRODUCTION`,
PI 2026-07-22) plus two dump-only markers that sit BELOW the ladder because they prove wiring
rather than physics: ``SCHEMA-PROOF`` (topology only) and ``CENSUS-WIRED`` (census, no physics gate).

**Two axes, not one** (PI decision D1 option B, 2026-07-28). The rung above is STRUCTURAL — how far
something has been wired and executed. Whether its MAGNITUDES may be quoted is a second, orthogonal
field, ``stage.quantitative_claim_status`` (``BLOCKED`` / ``OPEN`` / ``CONFIRMED``), and this page
shows both because a high rung with a blocked magnitude is the expected success case under the
2026-07-25 reframe, not a shortfall. Dumps predating the vocabulary show ``BLOCKED``, which is
accurate: none of them validated a magnitude.

Usage:
    python -m aleph.scripts.ac_architecture_dashboard --out outputs/ac/figs/architecture.html
"""
from __future__ import annotations

import argparse
import html
import json
import math
from pathlib import Path

import numpy as np

# The rung vocabulary and its order are the ENGINE's, not this script's (PI D1-B, 2026-07-28). They
# used to be a literal tuple right here, which is how drivers ended up hand-typing rungs as free
# strings; only the colours below are this script's business.
from aleph.engine.contracts import (  # noqa: E402 - vocabulary import, kept beside its use
    EVIDENCE_ORDER,
    LADDER_FLOOR,
    QuantitativeClaim,
    rung_rank as evidence_rank,
)

EVIDENCE_COLOR = {
    "SCHEMA-PROOF": "#5b6478", "CENSUS-WIRED": "#6d7fa3",
    "CONTRACTED": "#7d6b4f", "SEAMED": "#a07d3c", "KERNEL_BOUND": "#c08a2e",
    "CUDA_UNIT": "#b8a022", "CONNECTED": "#3f9d5a", "NATIVE": "#2bb673",
    "OPTIMISED": "#1f9ec4", "PRODUCTION": "#7b5cd6",
}
ROLE_COLOR = {
    "surface_body": "#e06c9f", "filament_network": "#4aa3ff", "fluid_volume": "#39c0c8",
    "deformable_body": "#c084fc", "actuator": "#ff8a5c", "environment": "#8a9bb5",
    "adhesion": "#ffd166", "boundary": "#7f8ea3",
}


def scan_dumps(outputs: Path) -> dict[str, dict]:
    """Scan v2 state dumps and return per-component observed evidence.

    Args:
        outputs: The ``aleph/outputs`` root to walk.

    Returns:
        ``{component_key: {"evidence": str, "quantitative": str, "dump": str, "stage": str,
        "seen_in": [names]}}`` where ``evidence`` is the highest STRUCTURAL rung of any dump the
        component appears in, and ``quantitative`` is that same dump's claim status.

    Note:
        The two axes are reported from the same dump rather than maxed independently: a rung and a
        claim status only mean anything as a pair, so taking the best of each across different runs
        would manufacture a run that never happened.  Dumps written before the two-axis vocabulary
        carry no claim status, and are reported as ``BLOCKED`` — the conservative default, matching
        the fact that no such artifact ever validated a magnitude.
    """
    best: dict[str, dict] = {}
    for path in sorted(outputs.rglob("*.npz")):
        try:
            data = np.load(path, allow_pickle=True)
            if "manifest_json" not in data.files:
                continue
            manifest = json.loads(str(data["manifest_json"]))
        except Exception:
            continue
        stage = manifest.get("stage") or {}
        ev = stage.get("evidence")
        if ev is None:
            continue
        declared = "quantitative_claim_status" in stage
        claim = str(stage.get("quantitative_claim_status") or QuantitativeClaim.BLOCKED.value)
        for comp in manifest.get("components", []):
            key = comp.get("key")
            if not key:
                continue
            slot = best.setdefault(
                key,
                {"evidence": None, "quantitative": "", "rung_verified": False,
                 "dump": "", "stage": "", "seen_in": []},
            )
            slot["seen_in"].append(path.name)
            if evidence_rank(ev) > evidence_rank(slot["evidence"]):
                slot.update(
                    evidence=str(ev), quantitative=claim, rung_verified=declared,
                    dump=path.name, stage=str(stage.get("name", "")),
                )
    return best


def collect() -> dict:
    """Assemble the dashboard payload: declared graph (code) + observed evidence (disk)."""
    from aleph.engine.contracts import reference_cell_architecture

    arch = reference_cell_architecture()
    outputs = Path(__file__).resolve().parents[1] / "outputs"
    observed = scan_dumps(outputs)

    comps = []
    for c in arch.components:
        obs = observed.get(c.name, {})
        comps.append({
            "name": c.name,
            "role": getattr(c.role, "value", str(c.role)),
            "representation": c.representation,
            "solver": c.solver,
            "owns_geometry": bool(c.owns_geometry),
            "dynamically_evolving": bool(c.dynamically_evolving),
            "has_events": bool(c.has_events),
            "evidence": obs.get("evidence"),
            "quantitative": obs.get("quantitative", ""),
            "rung_verified": bool(obs.get("rung_verified", False)),
            "evidence_dump": obs.get("dump", ""),
            "evidence_stage": obs.get("stage", ""),
            "n_dumps": len(obs.get("seen_in", [])),
        })
    conns = []
    for k in arch.connectors:
        conns.append({
            "name": k.name,
            "family": getattr(k.family, "value", str(k.family)),
            "a": k.component_a,
            "b": k.component_b,
            "kinetics": bool(k.kinetics),
            "commit_on_accept": bool(k.commit_on_accept),
            "bidirectional": bool(k.bidirectional),
            "adjoint": bool(k.adjoint_transfer_required),
            "scope": getattr(k.scope, "value", str(k.scope)),
            "chemistry": k.chemistry_card or "",
            "role_a": k.endpoint_role_a,
            "role_b": k.endpoint_role_b,
        })
    return {"components": comps, "connectors": conns}


def _layout(n: int, radius: float = 250.0) -> list[tuple[float, float]]:
    """Evenly place ``n`` nodes on a circle centred at the origin."""
    return [(radius * math.cos(2 * math.pi * i / n - math.pi / 2),
             radius * math.sin(2 * math.pi * i / n - math.pi / 2)) for i in range(n)]


def render(payload: dict, out: Path, title: str) -> Path:
    """Write the self-contained dashboard HTML (inline SVG + JS; no external assets)."""
    comps = payload["components"]
    pos = {c["name"]: p for c, p in zip(comps, _layout(len(comps)))}
    counts: dict[str, int] = {}
    for c in comps:
        key = c["evidence"] or "—"
        counts[key] = counts.get(key, 0) + 1
    on_ladder = sum(1 for c in comps if evidence_rank(c["evidence"]) >= LADDER_FLOOR)

    doc = {
        "components": comps, "connectors": payload["connectors"],
        "pos": {k: [round(v[0], 1), round(v[1], 1)] for k, v in pos.items()},
        "evidence_color": EVIDENCE_COLOR, "role_color": ROLE_COLOR,
        "order": [str(rung) for rung in EVIDENCE_ORDER], "ladder_floor": LADDER_FLOOR,
        "summary": {"n_components": len(comps), "n_connectors": len(payload["connectors"]),
                    "on_ladder": on_ladder, "counts": counts},
        "title": title,
    }
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(_HTML.replace("/*__PAYLOAD__*/", json.dumps(doc)), encoding="utf-8")
    return out


_HTML = r"""<!doctype html><html><head><meta charset="utf-8"><title>ac/engine architecture</title>
<style>
 :root{--bg:#0e1117;--panel:#161b22;--line:#30363d;--fg:#c9d1d9;--dim:#8b949e}
 *{box-sizing:border-box} body{margin:0;background:var(--bg);color:var(--fg);
   font-family:system-ui,-apple-system,Segoe UI,sans-serif;font-size:13px}
 header{padding:14px 18px;border-bottom:1px solid var(--line)}
 h1{margin:0 0 4px;font-size:16px;font-weight:650}
 .sub{color:var(--dim);font-size:12px;line-height:1.5}
 .wrap{display:flex;gap:14px;padding:14px;align-items:flex-start;flex-wrap:wrap}
 .card{background:var(--panel);border:1px solid var(--line);border-radius:10px;padding:12px}
 #graph{flex:1 1 640px;min-width:340px}
 svg{width:100%;height:auto;max-height:74vh;display:block}   /* keep the whole ring on screen */
 #side{flex:0 1 380px;min-width:300px;max-height:78vh;overflow-y:auto}
 .node{cursor:pointer} .node text{font-size:11px;fill:var(--fg);pointer-events:none}
 .edge{stroke-opacity:.5;cursor:pointer} .edge.dim{stroke-opacity:.06}
 .node.dim{opacity:.18}
 .chip{display:inline-block;padding:1px 7px;border-radius:999px;font-size:11px;
   border:1px solid var(--line);margin:2px 3px 2px 0;white-space:nowrap}
 .row{display:flex;justify-content:space-between;gap:10px;padding:4px 0;border-bottom:1px solid #ffffff10}
 .k{color:var(--dim)} .v{text-align:right}
 .bar{display:flex;height:9px;border-radius:5px;overflow:hidden;margin:8px 0 4px}
 .bar i{display:block;height:100%}
 .legend{display:flex;flex-wrap:wrap;gap:6px;margin-top:6px}
 button{background:#21262d;color:var(--fg);border:1px solid var(--line);border-radius:6px;
   padding:4px 9px;cursor:pointer;font-size:12px} button:hover{background:#2d333b}
 .warn{color:#e3b341}
</style></head><body>
<header>
  <h1 id="ttl"></h1>
  <div class="sub" id="sub"></div>
</header>
<div class="wrap">
  <div class="card" id="graph">
    <div style="margin-bottom:8px">
      <button onclick="setFilter('all')">all</button>
      <button onclick="setFilter('kinetic')">kinetic connectors</button>
      <button onclick="setFilter('ladder')">on the ladder</button>
      <button onclick="setFilter('unrun')">declared, not yet run</button>
      <span class="sub" id="fnote" style="margin-left:8px"></span>
    </div>
    <svg id="svg" viewBox="-310 -300 620 610" preserveAspectRatio="xMidYMid meet"></svg>
    <div class="legend" id="leg"></div>
  </div>
  <div class="card" id="side"><div class="sub">Click a component or connector.</div></div>
</div>
<script>
const P = /*__PAYLOAD__*/;
const svg=document.getElementById('svg'), side=document.getElementById('side');
const rank=n=>P.order.indexOf(String(n));
const evColor=n=>P.evidence_color[n]||'#3a4250';
document.getElementById('ttl').textContent=P.title;
const s=P.summary;
document.getElementById('sub').innerHTML=
  `${s.n_components} components · ${s.n_connectors} connectors — graph read from `+
  `<code>contracts.reference_cell_architecture()</code>; evidence read from each state dump's own `+
  `<code>stage.evidence</code>. <b>${s.on_ladder}/${s.n_components}</b> components have appeared in a dump `+
  `at or above <code>CONTRACTED</code>; the rest are declared but have only ever appeared in wiring-level dumps.`;

// ---- edges first so nodes draw on top
const NS='http://www.w3.org/2000/svg';
function el(t,a){const e=document.createElementNS(NS,t);for(const k in a)e.setAttribute(k,a[k]);return e;}
const famColors={};
let hue=0; for(const c of P.connectors){ if(!(c.family in famColors)){ famColors[c.family]=`hsl(${hue},58%,60%)`; hue=(hue+47)%360; } }

const edgeEls=[], nodeEls=[];
for(const c of P.connectors){
  const a=P.pos[c.a], b=P.pos[c.b]; if(!a||!b) continue;
  const bow=(a[0]*b[1]-a[1]*b[0])>0?0.18:-0.18;         // curve so opposite pairs don't overlap
  const mx=(a[0]+b[0])/2 - (b[1]-a[1])*bow, my=(a[1]+b[1])/2 + (b[0]-a[0])*bow;
  const p=el('path',{d:`M${a[0]},${a[1]} Q${mx},${my} ${b[0]},${b[1]}`,fill:'none',
    stroke:famColors[c.family],'stroke-width':c.kinetics?2.1:1.2,
    'stroke-dasharray':c.kinetics?'':'4,3',class:'edge'});
  p.addEventListener('click',e=>{e.stopPropagation();showConn(c);});
  p.__c=c; svg.appendChild(p); edgeEls.push(p);
}
for(const c of P.components){
  const [x,y]=P.pos[c.name];
  const g=el('g',{class:'node',transform:`translate(${x},${y})`});
  const r=13+Math.min(9,c.n_dumps*1.6);
  g.appendChild(el('circle',{r:r,fill:evColor(c.evidence),
    stroke:P.role_color[c.role]||'#8a9bb5','stroke-width':3}));
  const t=el('text',{y:r+13,'text-anchor':'middle'}); t.textContent=c.name; g.appendChild(t);
  g.addEventListener('click',e=>{e.stopPropagation();showComp(c);});
  g.__c=c; svg.appendChild(g); nodeEls.push(g);
}
// legend
document.getElementById('leg').innerHTML =
  P.order.filter(o=>P.components.some(c=>c.evidence===o)||rank(o)>=P.ladder_floor)
   .map(o=>`<span class="chip" style="border-color:${evColor(o)};color:${evColor(o)}">${o}</span>`).join('')
  +`<span class="chip" style="color:var(--dim)">solid edge = kinetic · dashed = static</span>`;

function esc(t){return String(t==null?'':t).replace(/[&<>]/g,m=>({'&':'&amp;','<':'&lt;','>':'&gt;'}[m]));}
function row(k,v){return `<div class="row"><span class="k">${k}</span><span class="v">${v}</span></div>`;}
function showComp(c){
  const onLadder=rank(c.evidence)>=P.ladder_floor;
  side.innerHTML=`<div style="font-weight:650;font-size:14px;margin-bottom:2px">${esc(c.name)}</div>
   <div class="sub" style="margin-bottom:8px">${esc(c.role)}</div>
   <span class="chip" style="border-color:${evColor(c.evidence)};color:${evColor(c.evidence)}">
     ${esc(c.evidence||'no dump')}</span>
   ${c.quantitative?`<span class="chip${c.quantitative==='CONFIRMED'?'':' warn'}">`+
     `magnitude ${esc(c.quantitative)}</span>`:''}
   ${onLadder?'':'<span class="chip warn">wiring-level only</span>'}
   ${c.evidence&&!c.rung_verified?'<span class="chip warn">rung self-declared '+
     '(dump predates the two-axis vocabulary)</span>':''}
   ${row('representation',esc(c.representation))}
   ${row('solver',esc(c.solver))}
   ${row('owns geometry',c.owns_geometry?'yes':'no')}
   ${row('dynamically evolving',c.dynamically_evolving?'yes':'no')}
   ${row('has events',c.has_events?'yes':'no')}
   ${row('appears in dumps',c.n_dumps)}
   ${c.evidence_dump?row('best evidence from',esc(c.evidence_dump)):''}
   ${c.evidence_stage?`<div class="sub" style="margin-top:6px">${esc(c.evidence_stage)}</div>`:''}
   <div style="margin-top:10px;font-weight:600">connectors</div>`+
   P.connectors.filter(k=>k.a===c.name||k.b===c.name).map(k=>
     `<div class="row"><span class="k">${esc(k.family)}</span>
      <span class="v">${esc(k.a===c.name?k.b:k.a)}${k.kinetics?' · kinetic':''}</span></div>`).join('');
}
function showConn(k){
  side.innerHTML=`<div style="font-weight:650;font-size:14px;margin-bottom:2px">${esc(k.name)}</div>
   <div class="sub" style="margin-bottom:8px">${esc(k.family)} · ${esc(k.scope)}</div>
   ${row('path',`${esc(k.a)} ⟷ ${esc(k.b)}`)}
   ${row('endpoint A',esc(k.role_a))}
   ${row('endpoint B',esc(k.role_b))}
   ${row('kinetic',k.kinetics?'yes':'no')}
   ${row('commit on accept',k.commit_on_accept?'yes':'no')}
   ${row('bidirectional',k.bidirectional?'yes':'no')}
   ${row('adjoint transfer',k.adjoint?'required':'—')}
   ${k.chemistry?row('chemistry card',esc(k.chemistry)):''}`;
}
function setFilter(mode){
  let note='';
  for(const g of nodeEls){
    const c=g.__c; let on=true;
    if(mode==='ladder') on=rank(c.evidence)>=P.ladder_floor;
    if(mode==='unrun')  on=rank(c.evidence)<P.ladder_floor;
    g.classList.toggle('dim',!on);
  }
  for(const p of edgeEls){
    const k=p.__c; let on=true;
    if(mode==='kinetic') on=k.kinetics;
    if(mode==='ladder'||mode==='unrun'){
      const ca=P.components.find(c=>c.name===k.a), cb=P.components.find(c=>c.name===k.b);
      const ok=c=>c&&(mode==='ladder'?rank(c.evidence)>=P.ladder_floor:rank(c.evidence)<P.ladder_floor);
      on=ok(ca)&&ok(cb);
    }
    p.classList.toggle('dim',!on);
  }
  if(mode==='kinetic') note=P.connectors.filter(k=>k.kinetics).length+' of '+P.connectors.length+' connectors are kinetic (commit only on an accepted step)';
  if(mode==='unrun')   note='these have only appeared in wiring-level dumps (SCHEMA-PROOF / CENSUS-WIRED)';
  if(mode==='ladder')  note='appeared in a dump at or above CONTRACTED';
  document.getElementById('fnote').textContent=note;
}
setFilter('all');
</script></body></html>
"""


def main() -> None:
    """CLI entry point."""
    ap = argparse.ArgumentParser(description=__doc__.split("\n", 1)[0])
    ap.add_argument("--out", type=Path,
                    default=Path("aleph/outputs/ac/figs/architecture_dashboard.html"))
    ap.add_argument("--title", default="ac/engine — declared cell graph vs what has actually run")
    args = ap.parse_args()

    payload = collect()
    out = render(payload, args.out, args.title)
    comps = payload["components"]
    on_ladder = sum(1 for c in comps if evidence_rank(c["evidence"]) >= LADDER_FLOOR)
    print(f"[dashboard] {len(comps)} components / {len(payload['connectors'])} connectors")
    print(f"[dashboard] on the ladder (>= CONTRACTED): {on_ladder}/{len(comps)}")
    for c in sorted(comps, key=lambda c: -evidence_rank(c["evidence"])):
        claim = f"magnitude {c['quantitative']}" if c["quantitative"] else ""
        mark = "" if c["rung_verified"] or not c["evidence"] else "  [rung self-declared]"
        print(
            f"    {c['name']:22s} {str(c['evidence'] or '—'):14s} {claim:22s}"
            f"{c['evidence_dump']}{mark}"
        )
    print(f"[dashboard] wrote {out}")


if __name__ == "__main__":
    main()
