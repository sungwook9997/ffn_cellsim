"""Self-contained interactive three.js HTML viewer for FF filament networks + compartments.

The FF-engine counterpart of ``dcm_mesh_viewer_html`` (PI 2026-07-02: "dcm 엔진처럼 인터렉티브하게"). FF
primitives are FILAMENTS (line segments) and point clouds (nucleus beads), not cell meshes — so this renders
``THREE.LineSegments`` + ``THREE.Points`` with the same UX: ROTATE / zoom / pan (OrbitControls) + a per-LAYER
CHECKBOX panel (PI 2026-07-23) that toggles each compartment / connector's visibility INDEPENDENTLY — many
visible at once, grouped compartments-vs-connectors — with a PRESET dropdown (the composed cell is the default,
all layers on; the isolation / cut-away / connector presets are convenience layer-sets that just re-seed the
checkboxes). Point size is a fixed sensible default (no UI slider). Optional rigid AFM plates are drawn as
translucent planes. One double-clickable offline HTML (three.js r128 from CDN; all geometry base64-embedded).

Use programmatically::

    from aleph.scripts.ff_viewer_html import build_viewer
    build_viewer(scenes={"isotropic": [layer, ...], ...}, out="viewer.html", title="FF cortex")

where each ``layer`` is ``{"name", "kind": "lines"|"points", "verts": np.ndarray, "color": "#rrggbb",
"size": float}``; ``lines`` verts are segment endpoint PAIRS shaped (Nseg, 2, 3) (or flat 6·Nseg); ``points``
verts are (Npts, 3). Pass ``plates=[z0, z1]`` per scene via a ``"plates"`` layer with kind ``"plates"``.

ANIMATION (PI 2026-07-02: "dcm 엔진처럼" — the DCM viewer plays frames): a ``lines`` / ``points`` layer may
carry ``"frames"``: a list of vert arrays (same shape as ``verts``), one per timestep. When any layer in the
vault has frames, a ▶ play button + frame slider appear and the geometry positions swap per frame — e.g. a
filopodium bundle visibly ELONGATING over the polymerization steps. ``verts`` (frame 0) is still required for
the bounding box.
"""
from __future__ import annotations

import base64
import json

import numpy as np


def _b64(arr: np.ndarray, dtype) -> str:
    return base64.b64encode(np.ascontiguousarray(arr, dtype=dtype).tobytes()).decode("ascii")


def _quantize_pos(arr: np.ndarray) -> tuple[np.ndarray, list[float]]:
    """Quantize an (N,3) position buffer to int16 within its own bounding box.

    Halves the payload versus float32 (gzip does not help here: float32 coordinate mantissas
    are high-entropy, measured 1.1x, whereas this is an exact 2.0x). The error is bounded by
    span/65534 per axis — for the 15 um native cell bbox that is ~2.3e-4 um, i.e. ~15x finer
    than one actin filament radius, so it is invisible at any zoom the viewer allows. This is
    a PRECISION choice in the transport encoding, never a downsample: every filament, node and
    segment is still present (PI 2026-07-07).

    Args:
        arr: (N,3) world-space positions.

    Returns:
        ``(q, meta)`` where ``q`` is the (N,3) int16 buffer and ``meta`` is
        ``[ox, oy, oz, sx, sy, sz]``: the per-axis origin and span needed to invert it.
    """
    a = np.ascontiguousarray(arr, dtype=np.float32).reshape(-1, 3)
    if a.size == 0:                      # empty layer (a component present but with no drawn nodes):
        return a.astype(np.int16), [0.0, 0.0, 0.0, 1.0, 1.0, 1.0]   # identity meta, empty buffer
    lo = a.min(axis=0)
    span = np.maximum(a.max(axis=0) - lo, 1e-9)
    q = np.rint((a - lo) / span * 65534.0 - 32767.0).astype(np.int16)
    return q, [float(lo[0]), float(lo[1]), float(lo[2]),
               float(span[0]), float(span[1]), float(span[2])]


def build_viewer(scenes: dict, out: str, title: str = "FF viewer", cbars: dict | None = None) -> str:
    """Write a self-contained interactive HTML viewer; returns the out path.

    Geometry buffers are POOLED + DEDUPLICATED: an identical vertex/colour/face buffer that recurs across scenes
    (e.g. the cortex line-vertex buffer, embedded once per scene it appears in — composed, isolation, cumulative,
    connector, milestone) is base64-encoded ONCE and every layer references it by integer index into ``blobs``.
    This is NOT a downsample — every filament/node is kept at full native resolution; only the repeated bytes are
    collapsed. Without it the native 2.9 M-node cortex appears ~6× and inflates the inline payload past what the
    browser can load (750 MB → blank viewer, empty preset dropdown); pooling brings that to a single copy.
    """
    all_pts = []
    payload_scenes = {}
    blobs: list[str] = []
    blob_index: dict[str, int] = {}
    quant: dict[str, list[float]] = {}   # blob index (as str) -> [ox,oy,oz,sx,sy,sz] for int16 position blobs

    def _ref(arr, dtype) -> int:
        """Base64-encode ``arr`` once and return its index in the shared ``blobs`` pool (dedup by exact bytes)."""
        b = _b64(arr, dtype)
        i = blob_index.get(b)
        if i is None:
            i = len(blobs)
            blob_index[b] = i
            blobs.append(b)
        return i

    def _ref_pos(arr) -> int:
        """Pool an (N,3) POSITION buffer as int16 + its dequantization meta (see :func:`_quantize_pos`)."""
        q, meta = _quantize_pos(arr)
        b = _b64(q, np.int16)
        key = "q:" + b                    # namespaced so a quantized blob never aliases a raw one
        i = blob_index.get(key)
        if i is None:
            i = len(blobs)
            blob_index[key] = i
            blobs.append(b)
            quant[str(i)] = meta
        return i

    for sname, layers in scenes.items():
        pl = []
        for L in layers:
            kind = L["kind"]
            if kind == "plates":
                pl.append({"name": L.get("name", "plates"), "kind": "plates",
                           "z": [float(z) for z in L["verts"]], "half_xy": float(L.get("half_xy", 12.0)),
                           "color": L.get("color", "#888888")})
                continue
            v = np.asarray(L["verts"], dtype=np.float32).reshape(-1, 3)
            all_pts.append(v)
            entry = {"name": L["name"], "kind": kind, "color": L.get("color", "#4aa3ff"),
                     "size": float(L.get("size", 2.0)), "n": int(v.shape[0]),
                     "b64": _ref_pos(v),
                     "opacity": float(L.get("opacity", 0.75 if kind == "lines" else 1.0)),
                     "on_top": bool(L.get("on_top", False)),   # on_top: depthTest off → draws through occluders
                     "clip": bool(L.get("clip", False)),       # clip: this layer is cut by the CUT plane
                     # per-layer checkbox panel (PI 2026-07-23): group = which panel section (compartment /
                     # connector / reference); swatch = the categorical colour chip even when color is white.
                     "group": L.get("group", "compartment"),
                     "swatch": L.get("swatch", L.get("color", "#4aa3ff"))}
            if "pt_size" in L:                                    # optional per-layer POINT size (else global PT_SIZE)
                entry["pt_size"] = float(L["pt_size"])            # small dots → a dense overlay doesn't overdraw the field
            if kind == "mesh":                                    # triangle surface: verts + faces + opacity
                entry["faces_b64"] = _ref(np.asarray(L["faces"], np.uint32).reshape(-1, 3), np.uint32)
                entry["opacity"] = float(L.get("opacity", 1.0))
            if L.get("frames") is not None:                       # animated layer: swap positions per timestep
                fr = [np.asarray(f, dtype=np.float32).reshape(-1, 3) for f in L["frames"]]
                entry["frames"] = [_ref_pos(f) for f in fr]
                all_pts.append(fr[-1])                            # last frame bounds the box (elongated tip)
            if L.get("color_frames") is not None:                 # FEM field: per-vertex per-frame colors (N,3) uint8
                entry["color_frames"] = [_ref(np.asarray(cf, np.uint8).reshape(-1, 3), np.uint8) for cf in L["color_frames"]]
                entry["cbar"] = L.get("cbar", "")                 # legend string for the field range
            pl.append(entry)
        payload_scenes[sname] = pl
    allp = np.concatenate(all_pts, axis=0) if all_pts else np.zeros((1, 3), np.float32)
    lo = allp.min(0).tolist(); hi = allp.max(0).tolist()
    n_frames = max([len(L.get("frames", [])) for ls in payload_scenes.values() for L in ls] + [0])
    payload = {"title": title, "scenes": payload_scenes, "n_frames": int(n_frames),
               "lo": lo, "hi": hi, "scene_names": list(scenes.keys()), "cbars": cbars or {},
               "blobs": blobs,   # shared base64 buffer pool; layer .b64/.faces_b64/.frames/.color_frames are indices
               "quant": quant}   # blob-index -> [ox,oy,oz,sx,sy,sz]; present iff that blob is an int16 position buffer
    html = _HTML.replace("/*__PAYLOAD__*/", json.dumps(payload))
    with open(out, "w") as f:
        f.write(html)
    return out


_HTML = r"""<!doctype html><html><head><meta charset="utf-8"><title>FF viewer</title>
<style>
 body{margin:0;overflow:hidden;background:#0e1117;font-family:system-ui,sans-serif;color:#c9d1d9}
 #ui{position:absolute;top:10px;left:10px;background:#161b22cc;padding:10px 12px;border-radius:8px;font-size:13px}
 #ui .row{margin:5px 0;display:flex;align-items:center;gap:8px}
 #ui label{min-width:52px;color:#8b949e} select{background:#0d1117;color:#c9d1d9;border:1px solid #30363d;border-radius:5px;padding:2px 6px}
 .hint{color:#6e7681;font-size:11px} #title{font-weight:600;margin-bottom:6px}
 #leg{position:absolute;bottom:10px;left:10px;background:#161b22ee;padding:8px 10px;border-radius:8px;font-size:12px;max-height:66vh;overflow-y:auto;min-width:190px;max-width:290px}
 #leg .panelhdr{display:flex;align-items:center;gap:8px;margin-bottom:4px;font-weight:600}
 #leg .panelhdr .lnk{color:#58a6ff;cursor:pointer;font-weight:400;font-size:11px}
 #leg .panelhdr .lnk:hover{text-decoration:underline}
 #leg .ghdr{display:flex;align-items:center;gap:6px;margin:7px 0 2px;color:#8b949e;font-weight:700;text-transform:uppercase;font-size:10px;letter-spacing:.05em}
 #leg .lrow{display:flex;align-items:center;gap:6px;margin:2px 0 2px 4px;cursor:pointer;line-height:1.35}
 #leg .lname{white-space:nowrap;overflow:hidden;text-overflow:ellipsis;max-width:220px}
 #leg input[type=checkbox]{margin:0;cursor:pointer;flex:none}
 .sw{display:inline-block;width:11px;height:11px;border-radius:2px;vertical-align:middle;flex:none}
 #cbar{position:absolute;bottom:16px;right:16px;background:#161b22dd;padding:8px 12px 8px 10px;border-radius:8px;font-size:11px;display:none}
 /* scale bar: the viewer is in physical um, so a zoom-reactive ruler keeps every screenshot
    self-describing (visualization rule: always annotate units). */
 #scale{position:absolute;bottom:16px;left:50%;transform:translateX(-50%);text-align:center;
   font-size:11px;color:#c9d1d9;text-shadow:0 1px 3px #000;pointer-events:none}
 #scalebar{height:3px;background:#c9d1d9;border-radius:2px;margin:0 auto 3px}
 #cbartitle{margin-bottom:6px;max-width:180px;color:#c9d1d9;font-weight:600}
 #cbarbar{width:16px;height:170px;border-radius:3px;display:inline-block;vertical-align:top;border:1px solid #30363d}
 #cbarticks{display:inline-block;height:170px;vertical-align:top;position:relative;margin-left:6px;width:64px}
 #cbarticks span{position:absolute;left:0;color:#c9d1d9;white-space:nowrap}
 #cbarticks span:before{content:"– ";color:#6e7681}
</style></head><body>
<div id="ui">
  <div id="title"></div>
  <div class="row"><label>preset</label><select id="scene"></select></div>
  <div class="row"><label>cut</label><input type="checkbox" id="cut"><select id="cutax"><option>x</option><option>y</option><option>z</option></select><input id="cutpos" type="range" min="-1" max="1" value="0" step="0.02" style="width:90px"><span class="hint">reveal nucleus/aster inside</span></div>
  <div class="row"><label>thickness</label><input id="lw" type="range" min="0.4" max="4" value="2" step="0.1" style="width:90px"><span class="hint" id="lwv">2.0×</span><span class="hint">filament line width</span></div>
  <div class="row" id="animrow" style="display:none"><label>frame</label>
    <button id="play">▶</button><input id="frame" type="range" min="0" max="0" value="0" step="1" style="width:120px">
    <span id="fnum" class="hint"></span></div>
  <div class="row hint">drag = rotate · scroll = zoom · right-drag = pan</div>
</div>
<div id="leg"></div>
<div id="cbar"><div id="cbartitle"></div><div id="cbarbar"></div><div id="cbarticks"></div></div>
<div id="scale"><div id="scalebar"></div><span id="scaletxt"></span></div>
<script src="https://cdnjs.cloudflare.com/ajax/libs/three.js/r128/three.min.js"></script>
<script src="https://cdn.jsdelivr.net/npm/three@0.128.0/examples/js/controls/OrbitControls.js"></script>
<!-- fat lines: WebGL ignores LineBasicMaterial.linewidth (every line is 1 px), which turns a dense
     filament network into a field of dots. LineSegments2 draws each segment as an instanced quad so
     the width is real and the network reads as filaments. -->
<script src="https://cdn.jsdelivr.net/npm/three@0.128.0/examples/js/lines/LineSegmentsGeometry.js"></script>
<script src="https://cdn.jsdelivr.net/npm/three@0.128.0/examples/js/lines/LineMaterial.js"></script>
<script src="https://cdn.jsdelivr.net/npm/three@0.128.0/examples/js/lines/LineSegments2.js"></script>
<script>
const P = /*__PAYLOAD__*/;
const B = P.blobs || [];                             // shared base64 buffer pool; layer buffer fields are indices into it
const PT_SIZE = 2.5;                                 // fixed sensible point size (PI 2026-07-23: dropped the UI slider)
document.getElementById('title').textContent = P.title;
function dec(b64){const s=atob(b64);const u=new Uint8Array(s.length);
  for(let i=0;i<s.length;i++)u[i]=s.charCodeAt(i);return new Float32Array(u.buffer);}
function dec32(b64){const s=atob(b64);const u=new Uint8Array(s.length);
  for(let i=0;i<s.length;i++)u[i]=s.charCodeAt(i);return new Uint32Array(u.buffer);}
function decU8(b64){const s=atob(b64);const u=new Uint8Array(s.length);
  for(let i=0;i<s.length;i++)u[i]=s.charCodeAt(i);return u;}   // (N*3,) uint8 per-vertex colors (FEM field)
const QM = P.quant || {};                            // blob index -> [ox,oy,oz,sx,sy,sz] for int16 position blobs
// POSITIONS: decode blob `i` back to float32 xyz. int16-quantized blobs (half the bytes of float32,
// error ~span/65534 ≈ 2e-4 um — far below one filament radius) carry a dequantization entry in QM;
// a blob without one is a legacy raw float32 buffer and is returned as-is.
function decPos(i){
  const u=decU8(B[i]); const m=QM[i];
  if(!m) return new Float32Array(u.buffer);
  const q=new Int16Array(u.buffer); const out=new Float32Array(q.length);
  const ox=m[0],oy=m[1],oz=m[2],sx=m[3],sy=m[4],sz=m[5];
  for(let k=0;k<q.length;k+=3){
    out[k]  =ox+(q[k]  +32767)/65534*sx;
    out[k+1]=oy+(q[k+1]+32767)/65534*sy;
    out[k+2]=oz+(q[k+2]+32767)/65534*sz;
  }
  return out;
}

const renderer=new THREE.WebGLRenderer({antialias:true});
renderer.setSize(innerWidth,innerHeight); renderer.setPixelRatio(devicePixelRatio);
renderer.localClippingEnabled=true;                 // CUT view: reveal interior compartments (nucleus/aster)
document.body.appendChild(renderer.domElement);
const clipPlane=new THREE.Plane(new THREE.Vector3(-1,0,0),0);   // normal set by the cut slider
let clipMats=[];                                    // materials of clip-tagged layers (cortex/outline)
let lineMats=[];                                    // LineMaterial instances (need .resolution on resize + width slider)
// 2.0 is the calibrated default: at 1.0 a dense native network still reads as a dot field (the
// thing fat lines exist to fix), and by ~3.5 the filaments merge into an opaque shell that hides
// the nucleus behind them. Purely a display control — no geometry or data depends on it.
let LW_SCALE=2.0;                                   // user line-width multiplier (thickness slider)
function applyClip(){ const on=document.getElementById('cut').checked;
  for(const m of clipMats){ m.clippingPlanes = on ? [clipPlane] : []; m.needsUpdate=true; } }  // needsUpdate: recompile shader for clipping
function updateClip(){ const ax=document.getElementById('cutax').value;
  // normal points toward the camera-near +axis half, so THAT half is clipped away → interior (nucleus/aster) revealed
  const nrm = ax==='x'? new THREE.Vector3(-1,0,0) : ax==='y'? new THREE.Vector3(0,-1,0) : new THREE.Vector3(0,0,-1);
  const t=parseFloat(document.getElementById('cutpos').value);  // -1..1: slide the cut plane across the cell
  const pt=new THREE.Vector3(cx,cy,cz);
  if(ax==='x') pt.x=cx+t*R*0.5; else if(ax==='y') pt.y=cy+t*R*0.5; else pt.z=cz+t*R*0.5;
  clipPlane.normal.copy(nrm); clipPlane.constant=-pt.dot(nrm);   // plane through pt with normal nrm
  applyClip(); }
const scene=new THREE.Scene(); scene.background=new THREE.Color(0x0e1117);
const cam=new THREE.PerspectiveCamera(45,innerWidth/innerHeight,0.05,5000);
const lo=P.lo, hi=P.hi;
const cx=(lo[0]+hi[0])/2, cy=(lo[1]+hi[1])/2, cz=(lo[2]+hi[2])/2;
const R=Math.max(hi[0]-lo[0],hi[1]-lo[1],hi[2]-lo[2])||20;
cam.position.set(cx+R*1.7,cy+R*1.1,cz+R*1.7);
const ctr=new THREE.OrbitControls(cam,renderer.domElement); ctr.target.set(cx,cy,cz); ctr.update();
scene.add(new THREE.AmbientLight(0xffffff,0.7));
const dl=new THREE.DirectionalLight(0xffffff,0.6); dl.position.set(1,1.5,2); scene.add(dl);

let current=[];   // three objects for the current scene
let animated=[];  // {obj, frames:[b64]} for layers that play over time
function clearScene(){ for(const o of current){ scene.remove(o); if(o.geometry)o.geometry.dispose(); if(o.material)o.material.dispose(); } current=[]; animated=[]; clipMats=[]; }
function _fmt(v){ const a=Math.abs(v); return (a>=1000||(a>0&&a<0.01))? v.toExponential(1) : (a>=10? v.toFixed(0) : v.toPrecision(3)); }
function updateCbar(name){   // paper-style colorbar: turbo gradient + value ticks for the active FEM field scene
  const cb=(P.cbars||{})[name]; const el=document.getElementById('cbar');
  if(!cb){ el.style.display='none'; return; }
  el.style.display='block';
  document.getElementById('cbartitle').textContent = cb.label + (cb.unit? '  ['+cb.unit+']' : '');
  document.getElementById('cbarbar').style.background = 'linear-gradient(to top,'+cb.grad.join(',')+')';
  const t=document.getElementById('cbarticks'); t.innerHTML=''; const H=170;
  for(let f=0; f<=5; f++){ const v=cb.lo+(cb.hi-cb.lo)*f/5; const s=document.createElement('span');
    s.style.bottom=(f/5*(H-2))+'px'; s.textContent=_fmt(v); t.appendChild(s); }
}

function setFrame(i){
  for(const a of animated){
    if(a.frames && i<a.frames.length){
      const v=decPos(a.frames[i]);
      if(a.obj.userData.fat){ a.obj.geometry.setPositions(v); }   // LineSegments2 keeps instanceStart/End, not `position`
      else { a.obj.geometry.setAttribute('position',new THREE.BufferAttribute(v,3)); }
      a.obj.geometry.attributes.position.needsUpdate=true; a.obj.geometry.computeBoundingSphere();
      if(a.obj.type==='Mesh') a.obj.geometry.computeVertexNormals();
    }
    if(a.color_frames && i<a.color_frames.length){   // FEM field: swap per-vertex colors this frame
      if(a.obj.userData.fat){ const c=decU8(B[a.color_frames[i]]); const cf=new Float32Array(c.length);
        for(let k=0;k<c.length;k++) cf[k]=c[k]/255; a.obj.geometry.setColors(cf); }
      else { a.obj.geometry.setAttribute('color', new THREE.BufferAttribute(decU8(B[a.color_frames[i]]),3,true));
      a.obj.geometry.attributes.color.needsUpdate=true; }
    }
  }
  document.getElementById('frame').value=i;
  document.getElementById('fnum').textContent=`${i}/${P.n_frames-1}`;
}

function showScene(name){
  clearScene();
  const layers=P.scenes[name]; const rows=[];   // rows: {name,color,swatch,group,objs} → one per-layer checkbox
  for(const L of layers){
    if(L.kind==='plates'){
      // plates are PERPENDICULAR to the compression axis (data-z); PlaneGeometry's default normal is +z,
      // so place the plane in the xy-plane at z=z_plate (NO rotation). (Earlier bug: rotated to xz + put at
      // y=z_plate → plates faced the wrong axis / floated off the cell.)
      const pobjs=[];
      for(const z of L.z){
        const g=new THREE.PlaneGeometry(L.half_xy*2,L.half_xy*2);
        const m=new THREE.MeshBasicMaterial({color:L.color,transparent:true,opacity:0.2,side:THREE.DoubleSide});
        const mesh=new THREE.Mesh(g,m); mesh.position.set(cx,cy,z);
        scene.add(mesh); current.push(mesh); pobjs.push(mesh);
      }
      rows.push({name:L.name||'rigid plate (⊥ compression axis)', color:L.color||'#888888',
                 group:L.group||'reference', objs:pobjs}); continue;
    }
    const v=decPos(L.b64); const geo=new THREE.BufferGeometry();
    geo.setAttribute('position',new THREE.BufferAttribute(v,3));
    let o;
    if(L.kind==='lines'){
      const op=(L.opacity===undefined)?0.75:L.opacity;
      const hasField=!!L.color_frames;                                   // per-segment FIELD colour (ECM tension/reorientation/displacement)
      // FAT LINES: a 1-px LineBasicMaterial makes a dense filament network read as a dot field.
      // LineSegments2 gives each segment a real screen-space width.
      const lg=new THREE.LineSegmentsGeometry(); lg.setPositions(v);
      if(hasField){ const c=decU8(B[L.color_frames[0]]); const cf=new Float32Array(c.length);
        for(let i=0;i<c.length;i++) cf[i]=c[i]/255; lg.setColors(cf); }
      const mat=new THREE.LineMaterial({color:new THREE.Color(L.color),vertexColors:hasField,
        transparent:true,opacity:op,depthTest:!L.on_top,
        linewidth:(L.size||1.5)*LW_SCALE, worldUnits:false, dashed:false});
      mat.resolution.set(innerWidth,innerHeight);
      lineMats.push(mat);
      o=new THREE.LineSegments2(lg,mat); o.computeLineDistances(); o.userData.fat=true;
      if(L.on_top){o.renderOrder=999;} scene.add(o); current.push(o);
    }else if(L.kind==='mesh'){
      const idx=dec32(B[L.faces_b64]); geo.setIndex(new THREE.BufferAttribute(idx,1)); geo.computeVertexNormals();
      const op=(L.opacity===undefined)?1.0:L.opacity;
      const hasField=!!L.color_frames;
      if(hasField){ geo.setAttribute('color', new THREE.BufferAttribute(decU8(B[L.color_frames[0]]),3,true)); }  // FEM field per-vertex color
      const mat=new THREE.MeshStandardMaterial({color:L.color,vertexColors:hasField,transparent:op<1.0,opacity:op,
        side:THREE.DoubleSide,roughness:0.55,metalness:0.0,flatShading:false,
        depthWrite:(op>=0.5)});   // faint envelopes (cortex hull) don't write depth → don't occlude the nucleus inside
      o=new THREE.Mesh(geo,mat); scene.add(o); current.push(o);
    }else{
      const hasFieldP=!!L.color_frames;                                   // per-point FIELD colour (e.g. CFD pore pressure)
      if(hasFieldP){ geo.setAttribute('color', new THREE.BufferAttribute(decU8(B[L.color_frames[0]]),3,true)); }
      const mat=new THREE.PointsMaterial({color:L.color,size:(L.pt_size||PT_SIZE),sizeAttenuation:true,vertexColors:hasFieldP});
      o=new THREE.Points(geo,mat); o.userData.pts=true; scene.add(o); current.push(o);
    }
    if(L.clip){ clipMats.push(o.material); }   // cortex/outline get cut; nucleus/aster stay whole
    // animate only when there is >1 frame to play (a single static color_frames snapshot is NOT animated)
    if((L.frames&&L.frames.length>1) || (L.color_frames&&L.color_frames.length>1)){ animated.push({obj:o, frames:L.frames, color_frames:L.color_frames}); }
    rows.push({name:L.name, color:L.color, swatch:L.swatch, group:L.group||'compartment', objs:[o]});
  }
  applyClip();
  updateCbar(name);
  const hasAnim = animated.length>0;
  document.getElementById('animrow').style.display = hasAnim ? 'flex' : 'none';
  if(hasAnim){ const nf=Math.max(...animated.map(a=>a.frames.length)); document.getElementById('frame').max=nf-1; setFrame(0); }
  buildPanel(rows);
}

// Per-layer visibility panel (PI 2026-07-23): one checkbox per compartment/connector, MANY visible at once,
// grouped compartments-vs-connectors. Replaces the passive legend + the old mutually-exclusive scene picker as
// the primary control — the preset dropdown just re-seeds these on scene rebuild. Data-driven: any new layer in
// the dump auto-gets a row. Toggling flips that layer's THREE object(s) .visible live (no geometry rebuild).
let panelRows=[];   // current scene's rows, each carrying a live ._cb checkbox
const _GROUPS=[['compartment','compartments'],['connector','connectors'],['reference','reference']];
function buildPanel(rows){
  panelRows=rows;
  const leg=document.getElementById('leg'); leg.innerHTML='';
  const hdr=document.createElement('div'); hdr.className='panelhdr';
  const t=document.createElement('span'); t.textContent='layers'; hdr.appendChild(t);
  const aOn=document.createElement('span'); aOn.className='lnk'; aOn.textContent='all'; aOn.onclick=()=>setAll(true);
  const aOff=document.createElement('span'); aOff.className='lnk'; aOff.textContent='none'; aOff.onclick=()=>setAll(false);
  hdr.appendChild(aOn); hdr.appendChild(aOff); leg.appendChild(hdr);
  for(const [gk,gl] of _GROUPS){
    const grows=rows.filter(r=>(r.group||'compartment')===gk);
    if(!grows.length) continue;
    const gh=document.createElement('div'); gh.className='ghdr';
    const gcb=document.createElement('input'); gcb.type='checkbox'; gcb.checked=true;   // group master → toggle whole group
    const gs=document.createElement('span'); gs.textContent=gl+' ('+grows.length+')';
    gh.appendChild(gcb); gh.appendChild(gs); leg.appendChild(gh);
    for(const r of grows){
      const row=document.createElement('label'); row.className='lrow';
      const cb=document.createElement('input'); cb.type='checkbox'; cb.checked=(r.on!==false);
      const sw=document.createElement('span'); sw.className='sw'; sw.style.background=(r.swatch||r.color||'#888');
      const nm=document.createElement('span'); nm.className='lname'; nm.textContent=r.name; nm.title=r.name;
      row.appendChild(cb); row.appendChild(sw); row.appendChild(nm); leg.appendChild(row);
      r._cb=cb; for(const o of r.objs) o.visible=cb.checked;
      cb.onchange=()=>{ for(const o of r.objs) o.visible=cb.checked;
        gcb.checked=grows.every(x=>x._cb.checked); requestRender(); };
    }
    gcb.onchange=()=>{ for(const r of grows){ r._cb.checked=gcb.checked; for(const o of r.objs) o.visible=gcb.checked; }
      requestRender(); };
    gcb.checked=grows.every(x=>x._cb.checked);
  }
}
function setAll(v){ for(const r of panelRows){ if(r._cb) r._cb.checked=v; for(const o of r.objs) o.visible=v; }
  document.querySelectorAll('#leg .ghdr input').forEach(g=>{ g.checked=v; }); requestRender(); }
const sel=document.getElementById('scene');
for(const n of P.scene_names){ const o=document.createElement('option'); o.value=n; o.textContent=n; sel.appendChild(o); }
sel.onchange=()=>showScene(sel.value);
document.getElementById('cut').onchange=updateClip;
document.getElementById('cutax').onchange=updateClip;
document.getElementById('cutpos').oninput=updateClip;

let playing=false, tacc=0;
document.getElementById('play').onclick=()=>{ playing=!playing; document.getElementById('play').textContent=playing?'⏸':'▶'; };
document.getElementById('frame').oninput=e=>{ playing=false; document.getElementById('play').textContent='▶'; setFrame(parseInt(e.target.value)); };
const _hp=new URLSearchParams(location.hash.slice(1));                 // #scene=NAME&frame=N|last auto-selects (for scripted per-scene/frame screenshots)
const _hs=_hp.get('scene');
const _init=(_hs&&P.scene_names.indexOf(_hs)>=0)?_hs:P.scene_names[0];
sel.value=_init; showScene(_init); updateClip();
const _hf=_hp.get('frame');
if(_hf!==null){ const _nf=parseInt(document.getElementById('frame').max)||0;
  const _fi=(_hf==='last')?_nf:Math.max(0,Math.min(_nf,parseInt(_hf)||0)); setFrame(_fi); }

// event-driven rendering: NO rAF when idle. Was an unconditional 60fps renderer.render →
// idle tabs pegged CPU/GPU forever, catastrophic under software-WebGL (2026-07-08).
// The loop reschedules only while playback runs or a repaint is pending, then stops →
// a still viewer costs zero CPU and headless capture can settle to idle.
let dirty=true, running=false, last=performance.now();
function wake(){ if(!running){ running=true; last=performance.now(); requestAnimationFrame(loop); } }
function requestRender(){ dirty=true; wake(); }
// zoom-reactive scale bar: convert the on-screen bar length to world um at the orbit target's
// depth, then snap it to a 1/2/5 x 10^n figure so the number stays readable.
function updateScale(){
  const dist=cam.position.distanceTo(ctr.target);
  const worldPerPx=2*dist*Math.tan(cam.fov*Math.PI/360)/innerHeight;
  if(!isFinite(worldPerPx)||worldPerPx<=0) return;
  const want=worldPerPx*130;                       // aim for a ~130 px bar
  const pow=Math.pow(10,Math.floor(Math.log10(want))), m=want/pow;
  const nice=(m<1.5?1:m<3.5?2:m<7.5?5:10)*pow;
  document.getElementById('scalebar').style.width=(nice/worldPerPx).toFixed(1)+'px';
  document.getElementById('scaletxt').textContent=
    (nice<0.1?nice.toFixed(3):nice<1?nice.toFixed(2):nice<10?nice.toFixed(1):nice.toFixed(0))+' µm';
}
ctr.addEventListener('change', ()=>{updateScale(); requestRender();});
updateScale();
for(const ev of ['input','change','click']) document.addEventListener(ev, requestRender);
// thickness slider: fat-line width is screen-space, so this is a pure display control —
// it changes no geometry and no data, only how readable the network is at a given zoom.
document.getElementById('lw').addEventListener('input',e=>{
  const f=parseFloat(e.target.value); document.getElementById('lwv').textContent=f.toFixed(1)+'×';
  const r=f/LW_SCALE; LW_SCALE=f;
  for(const m of lineMats){ m.linewidth*=r; m.needsUpdate=true; }
  requestRender();
});
addEventListener('resize',()=>{cam.aspect=innerWidth/innerHeight;cam.updateProjectionMatrix();renderer.setSize(innerWidth,innerHeight);
  for(const m of lineMats){ m.resolution.set(innerWidth,innerHeight); }   // fat lines are screen-space: width depends on it
  requestRender();});
function loop(){ ctr.update();
  const now=performance.now(); const dt=(now-last)/1000; last=now;
  if(playing && animated.length){ tacc+=dt; if(tacc>0.06){ tacc=0;
    let i=(parseInt(document.getElementById('frame').value)+1); const nf=parseInt(document.getElementById('frame').max);
    if(i>nf) i=0; setFrame(i); dirty=true; } }
  if(dirty){ dirty=false; renderer.render(scene,cam); }
  if(playing || dirty){ requestAnimationFrame(loop); } else { running=false; }}
requestRender();
</script></body></html>"""
