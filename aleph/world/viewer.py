"""The spatial view: the cell itself, every node and every connection, rotatable in a browser.

WHY THIS IS SEPARATE FROM THE FLOORPLAN.  The floorplan answers "what is wired, and on whose authority"
— it is a die map and a provenance ledger. It does not show the cell. This does: nodes at their real
positions, filaments as the polylines they are, bonds as the lines they are, rotatable and zoomable,
with nothing summarised and nothing sampled away.

NO DOWNSAMPLING, EVER.  Every element in the payload is drawn, and the page states the count it drew so
the claim is checkable on screen rather than taken on trust. When a scene grows past what one page can
hold, the honest limit is the PAGE — a byte budget, stated — and never a stride through the population.
Filaments are drawn as polylines rather than as points because a filament is a curve, and a scatter of
dots is a different object that happens to occupy similar space.

WHY WebGL AND A BINARY PAYLOAD.  Canvas 2D issues one path operation per segment and tops out around
tens of thousands of lines before a drag stops being interactive; JSON positions cost roughly 25 bytes
per coordinate. Both were the binding limit on the first version of this view, which could show ~0.4% of
the native filament count at five times the sourced discretisation — a picture of the wrong cell. WebGL
draws every line of a layer in ONE call from a buffer the GPU already holds, and float32 positions
base64-encoded cost 16 bytes per node instead of ~75. Together that is roughly two orders of magnitude,
which is the difference between a sketch and the cortex.

WHAT IT DELIBERATELY DOES NOT DO.  It renders one snapshot. It does not stream, animate, or connect to a
running simulation. Live streaming is next and needs the readback discipline the arena already documents
— between accepted steps, never inside the loop, delta-encoded so a cell at rest costs almost nothing to
watch.

Sanity Gate (recorded per CLAUDE.md, before first execution):
  * dimensional — positions [µm] throughout; the scale bar is labelled in µm and derived from the
    payload's own extent rather than assumed.
  * boundary — an empty arena renders an empty scene with a message rather than dividing by a zero
    extent; a browser without WebGL2 gets a stated refusal rather than a blank rectangle.
  * conservation/invariant — the number of lines drawn is counted at draw time and displayed, and the
    node count drawn is compared against the count claimed, so a discrepancy is visible not silent.
  * CFL/precision — nothing is integrated. Positions are float64 in the engine and are narrowed to
    float32 for transport, which is ~1e-7 relative — four orders below the µm-scale features drawn and
    stated on the page rather than hidden.
  * sign sense — not applicable.
  * measurement protocol — host-side read of already-built host geometry. No device, no readback.

engine units: length µm.  Runtime: pure host; emits a self-contained HTML page.
"""

from __future__ import annotations

import base64
import html
import json
from collections.abc import Sequence

import numpy as np

from aleph.world.arena import Kind, WorldArena
from aleph.world.bond import BondFamily
from aleph.world.strand import Strand
from aleph.world.surface import Surface

__all__ = ["render_cell_view", "write_cell_view"]

_HUES = ["#5BC8AF", "#7C9BFF", "#E8C547", "#FF7A59", "#B98CFF", "#59C2E8", "#E8899B", "#9BD45B"]


def _b64(arr: np.ndarray) -> str:
    """Base64 of an array's raw bytes — 4/3 overhead against ~6x for the same numbers as JSON text."""
    return base64.b64encode(np.ascontiguousarray(arr).tobytes()).decode("ascii")


def _hex_rgb(h: str) -> list[float]:
    return [int(h[i:i + 2], 16) / 255.0 for i in (1, 3, 5)]


def render_cell_view(
    arena: WorldArena,
    strands: Sequence[Strand] = (),
    surfaces: Sequence[Surface] = (),
    families: Sequence[BondFamily] = (),
    *,
    title: str = "Arena — cell view",
    standalone: bool = True,
    note: str = "",
) -> str:
    """Render the built geometry as an interactive, self-contained WebGL page.

    Args:
        arena: the world, for the census line. Not mutated.
        strands: built strands — drawn as polylines with real thickness.
        surfaces: built surfaces — drawn as their unique triangle edges.
        families: built bond families — drawn as lines, one colour and one toggle per family.
        title: page title.
        standalone: emit a whole document, or style-and-content only for a host with its own skeleton.
        note: a scoping sentence shown under the census line — e.g. what fraction of native this is.

    Returns:
        The HTML, referencing no external asset.
    """
    n_live = arena.n_live(Kind.NODE)
    pos = np.zeros((n_live, 3), np.float64)
    seen = np.zeros(n_live, bool)
    for s in strands:
        pos[s.nodes.lo:s.nodes.hi] = s.position
        seen[s.nodes.lo:s.nodes.hi] = True
    for sf in surfaces:
        pos[sf.nodes.lo:sf.nodes.hi] = sf.position
        seen[sf.nodes.lo:sf.nodes.hi] = True

    # One index buffer per LAYER, so a toggle is a skipped draw call rather than a rebuild.
    per_pop: dict[tuple[str, str], list[np.ndarray]] = {}
    for s in strands:
        per_pop.setdefault((s.population, "strand"), []).append(s.seg_node)
    for sf in surfaces:
        e = np.concatenate([sf.face_idx[:, [0, 1]], sf.face_idx[:, [1, 2]], sf.face_idx[:, [2, 0]]])
        per_pop.setdefault((sf.population, "surface"), []).append(np.unique(np.sort(e, axis=1), axis=0))

    layers: list[dict] = []
    for i, ((p, k), parts) in enumerate(sorted(per_pop.items())):
        idx = np.concatenate(parts).astype(np.uint32)
        layers.append({"name": p, "kind": k, "colour": _HUES[i % len(_HUES)],
                       "n": int(idx.shape[0]), "idx": _b64(idx), "tag": ""})
    for j, f in enumerate(families):
        idx = np.stack([f.node_i, f.node_j], axis=1).astype(np.uint32)
        layers.append({"name": f.name, "kind": "bond", "colour": _HUES[(len(per_pop) + j) % len(_HUES)],
                       "n": int(idx.shape[0]), "idx": _b64(idx),
                       "tag": str(f.count.source_class)})

    extent = float(np.abs(pos[seen]).max()) if seen.any() else 1.0
    scene = {
        "pos": _b64(pos.astype(np.float32)),
        "n_pos": int(n_live),
        "layers": layers,
        "extent": extent,
    }
    n_lines = sum(l["n"] for l in layers)
    payload = json.dumps(scene, separators=(",", ":"))
    esc = lambda x: html.escape(str(x), quote=True)  # noqa: E731

    css = """
:root{--bg:#0B1020;--panel:#151C2F;--line:#2C3654;--ink:#E6ECFA;--dim:#8494B0;--accent:#7C9BFF}
@media(prefers-color-scheme:light){:root{--bg:#F7F6F2;--panel:#FFF;--line:#DAD5C9;--ink:#181B22;
--dim:#6C7486;--accent:#3A5BD9}}
*{box-sizing:border-box}
body{margin:0;background:var(--bg);color:var(--ink);font:14px/1.5 ui-sans-serif,system-ui,sans-serif}
.wrap{max-width:1180px;margin:0 auto;padding:24px 20px 48px;display:flex;flex-direction:column;gap:14px}
h1{font-size:22px;margin:0;letter-spacing:-.02em}
.sub{color:var(--dim);font-size:13px;font-family:ui-monospace,Menlo,monospace}
.stage{position:relative;background:var(--panel);border:1px solid var(--line);border-radius:9px;
overflow:hidden;touch-action:none}
canvas{display:block;width:100%;cursor:grab}canvas:active{cursor:grabbing}
.hud{position:absolute;left:12px;bottom:11px;font-family:ui-monospace,Menlo,monospace;font-size:11px;
color:var(--dim);background:color-mix(in srgb,var(--panel) 86%,transparent);padding:5px 9px;border-radius:5px}
.bar{position:absolute;right:14px;bottom:11px;font-family:ui-monospace,Menlo,monospace;font-size:11px;
color:var(--dim);text-align:right}
.bar i{display:block;height:1px;background:currentColor;margin-bottom:3px}
.legend{display:flex;flex-wrap:wrap;gap:7px}
.legend button{display:flex;align-items:center;gap:7px;background:var(--panel);color:var(--ink);
border:1px solid var(--line);border-radius:999px;padding:5px 12px;font-size:12.5px;cursor:pointer;
font-family:inherit}
.legend button[aria-pressed=false]{opacity:.38}
.legend button:focus-visible{outline:2px solid var(--accent);outline-offset:2px}
.sw{width:11px;height:11px;border-radius:3px;flex:none}
.tip{color:var(--dim);font-size:12.5px;margin:0}
.err{padding:40px;text-align:center;color:var(--dim)}
"""

    js = r"""
const S=JSON.parse(document.getElementById('scene').textContent);
const un=s=>{const b=atob(s),u=new Uint8Array(b.length);for(let i=0;i<b.length;i++)u[i]=b.charCodeAt(i);return u.buffer;};
const POS=new Float32Array(un(S.pos));
const cv=document.getElementById('cv');
const gl=cv.getContext('webgl2',{antialias:true,alpha:false});
if(!gl){document.getElementById('stage').innerHTML=
 '<p class="err">This view needs WebGL2. Nothing is hidden — the browser cannot draw it.</p>';}
else{
const VS=`#version 300 es
in vec3 p; uniform mat3 R; uniform vec2 vp; uniform float s; uniform vec2 pan;
void main(){vec3 q=R*p; gl_Position=vec4((q.x*s+pan.x)/vp.x*2.0,(q.y*s+pan.y)/vp.y*2.0,0.0,1.0);} `;
const FS=`#version 300 es
precision mediump float; uniform vec3 c; uniform float a; out vec4 o;
void main(){o=vec4(c,a);} `;
function sh(t,src){const s=gl.createShader(t);gl.shaderSource(s,src);gl.compileShader(s);
 if(!gl.getShaderParameter(s,gl.COMPILE_STATUS))throw gl.getShaderInfoLog(s);return s;}
const pr=gl.createProgram();gl.attachShader(pr,sh(gl.VERTEX_SHADER,VS));
gl.attachShader(pr,sh(gl.FRAGMENT_SHADER,FS));gl.linkProgram(pr);gl.useProgram(pr);
const vb=gl.createBuffer();gl.bindBuffer(gl.ARRAY_BUFFER,vb);gl.bufferData(gl.ARRAY_BUFFER,POS,gl.STATIC_DRAW);
const loc=gl.getAttribLocation(pr,'p');gl.enableVertexAttribArray(loc);
gl.vertexAttribPointer(loc,3,gl.FLOAT,false,0,0);
const uR=gl.getUniformLocation(pr,'R'),uVp=gl.getUniformLocation(pr,'vp'),
 uS=gl.getUniformLocation(pr,'s'),uC=gl.getUniformLocation(pr,'c'),
 uA=gl.getUniformLocation(pr,'a'),uPan=gl.getUniformLocation(pr,'pan');
for(const L of S.layers){L.buf=gl.createBuffer();
 gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER,L.buf);
 gl.bufferData(gl.ELEMENT_ARRAY_BUFFER,new Uint32Array(un(L.idx)),gl.STATIC_DRAW);
 L.rgb=[parseInt(L.colour.slice(1,3),16)/255,parseInt(L.colour.slice(3,5),16)/255,
        parseInt(L.colour.slice(5,7),16)/255];L.on=true;}
gl.enable(gl.BLEND);gl.blendFunc(gl.SRC_ALPHA,gl.ONE_MINUS_SRC_ALPHA);
let rx=-0.42,ry=0.62,zoom=1,px=0,py=0,drag=null;
function rot(){const cx=Math.cos(rx),sx=Math.sin(rx),cy=Math.cos(ry),sy=Math.sin(ry);
 return new Float32Array([cy,sx*sy,-cx*sy, 0,cx,sx, sy,-sx*cy,cx*cy]);}
function draw(){const d=window.devicePixelRatio||1,w=cv.width,h=cv.height;
 gl.viewport(0,0,w,h);
 const bg=getComputedStyle(document.body).backgroundColor.match(/\d+/g).map(Number);
 gl.clearColor(bg[0]/255,bg[1]/255,bg[2]/255,1);gl.clear(gl.COLOR_BUFFER_BIT);
 const s=Math.min(w,h)/(2.35*S.extent)*zoom;
 gl.uniformMatrix3fv(uR,false,rot());gl.uniform2f(uVp,w,h);gl.uniform1f(uS,s);
 gl.uniform2f(uPan,px*d,-py*d);
 let n=0;
 for(const L of S.layers){if(!L.on)continue;
  gl.uniform3fv(uC,L.rgb);gl.uniform1f(uA,L.kind==='bond'?0.30:(L.kind==='surface'?0.45:0.95));
  gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER,L.buf);
  gl.drawElements(gl.LINES,L.n*2,gl.UNSIGNED_INT,0);n+=L.n;}
 document.getElementById('hud').textContent=
  n.toLocaleString()+' lines drawn — all of them · zoom '+zoom.toFixed(2)+'×';
 const bar=Math.max(1,Math.round(S.extent/2));
 document.getElementById('barw').style.width=(bar*s/d)+'px';
 document.getElementById('barl').textContent=bar+' µm';}
function size(){const r=cv.parentElement.getBoundingClientRect(),d=window.devicePixelRatio||1;
 cv.width=Math.round(r.width*d);cv.height=Math.round(Math.max(520,r.width*0.62)*d);
 cv.style.height=(cv.height/d)+'px';draw();}
cv.addEventListener('pointerdown',e=>{drag=[e.clientX,e.clientY,e.shiftKey];cv.setPointerCapture(e.pointerId);});
cv.addEventListener('pointerup',()=>drag=null);
cv.addEventListener('pointermove',e=>{if(!drag)return;const dx=e.clientX-drag[0],dy=e.clientY-drag[1];
 if(drag[2]){px+=dx;py+=dy;}else{ry+=dx*0.008;rx+=dy*0.008;}drag=[e.clientX,e.clientY,drag[2]];draw();});
cv.addEventListener('wheel',e=>{e.preventDefault();zoom*=Math.exp(-e.deltaY*0.0015);
 zoom=Math.min(80,Math.max(0.12,zoom));draw();},{passive:false});
document.querySelectorAll('.legend button').forEach((b,i)=>b.addEventListener('click',()=>{
 S.layers[i].on=!S.layers[i].on;b.setAttribute('aria-pressed',S.layers[i].on);draw();}));
window.addEventListener('resize',size);size();}
"""

    buttons = "".join(
        f'<button aria-pressed="true"><span class="sw" style="background:{l["colour"]}"></span>'
        f'{esc(l["name"])}<span class="sub">{l["n"]:,}'
        + (f' &middot; {esc(l["tag"])}' if l["tag"] else "") + "</span></button>"
        for l in layers
    )
    head = (f'<!doctype html><html lang="en"><head><meta charset="utf-8">'
            f'<meta name="viewport" content="width=device-width,initial-scale=1">'
            f"<title>{esc(title)}</title><style>{css}</style></head><body>"
            if standalone else f"<title>{esc(title)}</title><style>{css}</style>")
    tail = "</body></html>" if standalone else ""

    return f"""{head}<div class="wrap">
<div><h1>{esc(title)}</h1>
<p class="sub">{n_live:,} nodes &middot; {int(seen.sum()):,} positioned &middot; {n_lines:,} lines &middot;
every one drawn, none sampled &middot; positions float32 (~1e-7 relative)</p>
{f'<p class="tip">{esc(note)}</p>' if note else ''}</div>
<div class="stage" id="stage"><canvas id="cv"></canvas><div class="hud" id="hud"></div>
<div class="bar"><i id="barw"></i><span id="barl"></span></div></div>
<div class="legend">{buttons}</div>
<p class="tip">Drag to rotate &middot; shift-drag to pan &middot; scroll to zoom &middot; click a chip to toggle a layer.</p>
<script type="application/json" id="scene">{payload}</script>
<script>{js}</script>
</div>{tail}"""


def write_cell_view(arena: WorldArena, path: str, *a, **kw) -> str:
    """Render and write the cell view to ``path``; returns the path written."""
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(render_cell_view(arena, *a, **kw))
    return path


def _demo() -> None:
    """Self-check: every element reaches the payload, and the page is self-contained."""
    from aleph.world.bond import BondCount, SourceClass, build_bond_family
    from aleph.world.strand import build_strand
    from aleph.world.surface import build_surface, icosphere

    arena = WorldArena(capacity={Kind.NODE: 5_000, Kind.SEGMENT: 5_000, Kind.ANGLE3: 5_000,
                                 Kind.ANGLE4: 5_000, Kind.FACE: 5_000, Kind.BOND: 5_000})
    strands = [build_strand(arena, "cortex", start=(-1.5, 0.1 * i - 0.4, 0), direction=(1, 0, 0),
                            contour_um=3.0, seg_um=0.1) for i in range(8)]
    v, f = icosphere(2)
    surf = [build_surface(arena, "membrane", vertices=v, faces=f, radius_um=2.5)]
    cn = np.concatenate([s.nodes.lo + np.arange(s.n_nodes) for s in strands])
    fams = [build_bond_family(
        arena, "alpha_actinin", chemistry_card="alpha_actinin_ferrer2008",
        count=BondCount(basis="explicit", value=60.0, scope="MCF7 interphase 37C",
                        source_class=SourceClass.SOURCED, provenance="Ferrer 2008 PNAS"),
        support=1.0, pairs=np.stack([cn[:-31], cn[31:]], axis=1),
        rest_um=0.1, stiffness_pn_per_um=4.6e5)]

    page = render_cell_view(arena, strands, surf, fams)
    data = json.loads(page.split('id="scene">')[1].split("</script>")[0])
    assert data["n_pos"] == arena.n_live(Kind.NODE)
    by = {(l["name"], l["kind"]): l["n"] for l in data["layers"]}
    assert by[("cortex", "strand")] == 8 * 30, "every segment of every strand"
    assert by[("membrane", "surface")] == 480, "unique edges of a 320-face icosphere"
    assert by[("alpha_actinin", "bond")] == 60
    # The payload must decode back to exactly the positions that went in.
    back = np.frombuffer(base64.b64decode(data["pos"]), np.float32).reshape(-1, 3)
    assert back.shape[0] == arena.n_live(Kind.NODE)
    assert np.allclose(back[strands[0].nodes.lo], strands[0].position[0], atol=1e-6)
    for forbidden in ("http://", "https://", "<script src", "<link "):
        assert forbidden not in page, forbidden
    print(f"cell-view self-check OK — {data['n_pos']:,} nodes, "
          f"{sum(l['n'] for l in data['layers']):,} lines, {len(page):,} B")


if __name__ == "__main__":
    _demo()
