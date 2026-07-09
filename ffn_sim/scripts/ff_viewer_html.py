"""Self-contained interactive three.js HTML viewer for FF filament networks + compartments.

The FF-engine counterpart of ``dcm_mesh_viewer_html`` (PI 2026-07-02: "dcm 엔진처럼 인터렉티브하게"). FF
primitives are FILAMENTS (line segments) and point clouds (nucleus beads), not cell meshes — so this renders
``THREE.LineSegments`` + ``THREE.Points`` with the same UX: ROTATE / zoom / pan (OrbitControls) + a SCENE
dropdown to switch between views (e.g. the three filament arrangements, or the whole-cell compartments).
Optional rigid AFM plates are drawn as translucent planes. One double-clickable offline HTML (three.js r128
from CDN; all geometry base64-embedded).

Use programmatically::

    from ffn_sim.scripts.ff_viewer_html import build_viewer
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


def build_viewer(scenes: dict, out: str, title: str = "FF viewer", cbars: dict | None = None) -> str:
    """Write a self-contained interactive HTML viewer; returns the out path."""
    all_pts = []
    payload_scenes = {}
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
                     "b64": _b64(v, np.float32),
                     "opacity": float(L.get("opacity", 0.75 if kind == "lines" else 1.0)),
                     "on_top": bool(L.get("on_top", False)),   # on_top: depthTest off → draws through occluders
                     "clip": bool(L.get("clip", False))}       # clip: this layer is cut by the CUT plane
            if kind == "mesh":                                    # triangle surface: verts + faces + opacity
                entry["faces_b64"] = _b64(np.asarray(L["faces"], np.uint32).reshape(-1, 3), np.uint32)
                entry["opacity"] = float(L.get("opacity", 1.0))
            if L.get("frames") is not None:                       # animated layer: swap positions per timestep
                fr = [np.asarray(f, dtype=np.float32).reshape(-1, 3) for f in L["frames"]]
                entry["frames"] = [_b64(f, np.float32) for f in fr]
                all_pts.append(fr[-1])                            # last frame bounds the box (elongated tip)
            if L.get("color_frames") is not None:                 # FEM field: per-vertex per-frame colors (N,3) uint8
                entry["color_frames"] = [_b64(np.asarray(cf, np.uint8).reshape(-1, 3), np.uint8) for cf in L["color_frames"]]
                entry["cbar"] = L.get("cbar", "")                 # legend string for the field range
            pl.append(entry)
        payload_scenes[sname] = pl
    allp = np.concatenate(all_pts, axis=0) if all_pts else np.zeros((1, 3), np.float32)
    lo = allp.min(0).tolist(); hi = allp.max(0).tolist()
    n_frames = max([len(L.get("frames", [])) for ls in payload_scenes.values() for L in ls] + [0])
    payload = {"title": title, "scenes": payload_scenes, "n_frames": int(n_frames),
               "lo": lo, "hi": hi, "scene_names": list(scenes.keys()), "cbars": cbars or {}}
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
 #leg{position:absolute;bottom:10px;left:10px;background:#161b22cc;padding:8px 10px;border-radius:8px;font-size:12px}
 .sw{display:inline-block;width:11px;height:11px;border-radius:2px;margin-right:5px;vertical-align:middle}
 #cbar{position:absolute;bottom:16px;right:16px;background:#161b22dd;padding:8px 12px 8px 10px;border-radius:8px;font-size:11px;display:none}
 #cbartitle{margin-bottom:6px;max-width:180px;color:#c9d1d9;font-weight:600}
 #cbarbar{width:16px;height:170px;border-radius:3px;display:inline-block;vertical-align:top;border:1px solid #30363d}
 #cbarticks{display:inline-block;height:170px;vertical-align:top;position:relative;margin-left:6px;width:64px}
 #cbarticks span{position:absolute;left:0;color:#c9d1d9;white-space:nowrap}
 #cbarticks span:before{content:"– ";color:#6e7681}
</style></head><body>
<div id="ui">
  <div id="title"></div>
  <div class="row"><label>scene</label><select id="scene"></select></div>
  <div class="row"><label>points</label><input id="psize" type="range" min="1" max="8" value="3" step="0.5"></div>
  <div class="row"><label>cut</label><input type="checkbox" id="cut"><select id="cutax"><option>x</option><option>y</option><option>z</option></select><input id="cutpos" type="range" min="-1" max="1" value="0" step="0.02" style="width:90px"><span class="hint">reveal nucleus/aster inside</span></div>
  <div class="row" id="animrow" style="display:none"><label>frame</label>
    <button id="play">▶</button><input id="frame" type="range" min="0" max="0" value="0" step="1" style="width:120px">
    <span id="fnum" class="hint"></span></div>
  <div class="row hint">drag = rotate · scroll = zoom · right-drag = pan</div>
</div>
<div id="leg"></div>
<div id="cbar"><div id="cbartitle"></div><div id="cbarbar"></div><div id="cbarticks"></div></div>
<script src="https://cdnjs.cloudflare.com/ajax/libs/three.js/r128/three.min.js"></script>
<script src="https://cdn.jsdelivr.net/npm/three@0.128.0/examples/js/controls/OrbitControls.js"></script>
<script>
const P = /*__PAYLOAD__*/;
document.getElementById('title').textContent = P.title;
function dec(b64){const s=atob(b64);const u=new Uint8Array(s.length);
  for(let i=0;i<s.length;i++)u[i]=s.charCodeAt(i);return new Float32Array(u.buffer);}
function dec32(b64){const s=atob(b64);const u=new Uint8Array(s.length);
  for(let i=0;i<s.length;i++)u[i]=s.charCodeAt(i);return new Uint32Array(u.buffer);}
function decU8(b64){const s=atob(b64);const u=new Uint8Array(s.length);
  for(let i=0;i<s.length;i++)u[i]=s.charCodeAt(i);return u;}   // (N*3,) uint8 per-vertex colors (FEM field)

const renderer=new THREE.WebGLRenderer({antialias:true});
renderer.setSize(innerWidth,innerHeight); renderer.setPixelRatio(devicePixelRatio);
renderer.localClippingEnabled=true;                 // CUT view: reveal interior compartments (nucleus/aster)
document.body.appendChild(renderer.domElement);
const clipPlane=new THREE.Plane(new THREE.Vector3(-1,0,0),0);   // normal set by the cut slider
let clipMats=[];                                    // materials of clip-tagged layers (cortex/outline)
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
      const v=dec(a.frames[i]); a.obj.geometry.setAttribute('position',new THREE.BufferAttribute(v,3));
      a.obj.geometry.attributes.position.needsUpdate=true; a.obj.geometry.computeBoundingSphere();
      if(a.obj.type==='Mesh') a.obj.geometry.computeVertexNormals();
    }
    if(a.color_frames && i<a.color_frames.length){   // FEM field: swap per-vertex colors this frame
      a.obj.geometry.setAttribute('color', new THREE.BufferAttribute(decU8(a.color_frames[i]),3,true));
      a.obj.geometry.attributes.color.needsUpdate=true;
    }
  }
  document.getElementById('frame').value=i;
  document.getElementById('fnum').textContent=`${i}/${P.n_frames-1}`;
}

function showScene(name){
  clearScene();
  const layers=P.scenes[name]; const leg=[];
  for(const L of layers){
    if(L.kind==='plates'){
      // plates are PERPENDICULAR to the compression axis (data-z); PlaneGeometry's default normal is +z,
      // so place the plane in the xy-plane at z=z_plate (NO rotation). (Earlier bug: rotated to xz + put at
      // y=z_plate → plates faced the wrong axis / floated off the cell.)
      for(const z of L.z){
        const g=new THREE.PlaneGeometry(L.half_xy*2,L.half_xy*2);
        const m=new THREE.MeshBasicMaterial({color:L.color,transparent:true,opacity:0.2,side:THREE.DoubleSide});
        const mesh=new THREE.Mesh(g,m); mesh.position.set(cx,cy,z);
        scene.add(mesh); current.push(mesh);
      }
      leg.push(['#888888','rigid plate (⊥ compression axis)']); continue;
    }
    const v=dec(L.b64); const geo=new THREE.BufferGeometry();
    geo.setAttribute('position',new THREE.BufferAttribute(v,3));
    let o;
    if(L.kind==='lines'){
      const op=(L.opacity===undefined)?0.75:L.opacity;
      const hasField=!!L.color_frames;                                   // per-segment FIELD colour (ECM tension/reorientation/displacement)
      if(hasField){ geo.setAttribute('color', new THREE.BufferAttribute(decU8(L.color_frames[0]),3,true)); }
      const mat=new THREE.LineBasicMaterial({color:L.color,vertexColors:hasField,transparent:true,opacity:op,
        depthTest:!L.on_top});
      o=new THREE.LineSegments(geo,mat); if(L.on_top){o.renderOrder=999;} scene.add(o); current.push(o);
    }else if(L.kind==='mesh'){
      const idx=dec32(L.faces_b64); geo.setIndex(new THREE.BufferAttribute(idx,1)); geo.computeVertexNormals();
      const op=(L.opacity===undefined)?1.0:L.opacity;
      const hasField=!!L.color_frames;
      if(hasField){ geo.setAttribute('color', new THREE.BufferAttribute(decU8(L.color_frames[0]),3,true)); }  // FEM field per-vertex color
      const mat=new THREE.MeshStandardMaterial({color:L.color,vertexColors:hasField,transparent:op<1.0,opacity:op,
        side:THREE.DoubleSide,roughness:0.55,metalness:0.0,flatShading:false,
        depthWrite:(op>=0.5)});   // faint envelopes (cortex hull) don't write depth → don't occlude the nucleus inside
      o=new THREE.Mesh(geo,mat); scene.add(o); current.push(o);
    }else{
      const mat=new THREE.PointsMaterial({color:L.color,size:parseFloat(document.getElementById('psize').value),sizeAttenuation:true});
      o=new THREE.Points(geo,mat); o.userData.pts=true; scene.add(o); current.push(o);
    }
    if(L.clip){ clipMats.push(o.material); }   // cortex/outline get cut; nucleus/aster stay whole
    if(L.frames || L.color_frames){ animated.push({obj:o, frames:L.frames, color_frames:L.color_frames}); }
    leg.push([L.color,L.name]);
  }
  applyClip();
  updateCbar(name);
  const hasAnim = animated.length>0;
  document.getElementById('animrow').style.display = hasAnim ? 'flex' : 'none';
  if(hasAnim){ const nf=Math.max(...animated.map(a=>a.frames.length)); document.getElementById('frame').max=nf-1; setFrame(0); }
  document.getElementById('leg').innerHTML=leg.map(([c,n])=>`<div><span class="sw" style="background:${c}"></span>${n}</div>`).join('');
}
const sel=document.getElementById('scene');
for(const n of P.scene_names){ const o=document.createElement('option'); o.value=n; o.textContent=n; sel.appendChild(o); }
sel.onchange=()=>showScene(sel.value);
document.getElementById('psize').oninput=e=>{ for(const o of current) if(o.userData.pts) o.material.size=parseFloat(e.target.value); };
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
ctr.addEventListener('change', requestRender);
for(const ev of ['input','change','click']) document.addEventListener(ev, requestRender);
addEventListener('resize',()=>{cam.aspect=innerWidth/innerHeight;cam.updateProjectionMatrix();renderer.setSize(innerWidth,innerHeight);requestRender();});
function loop(){ ctr.update();
  const now=performance.now(); const dt=(now-last)/1000; last=now;
  if(playing && animated.length){ tacc+=dt; if(tacc>0.06){ tacc=0;
    let i=(parseInt(document.getElementById('frame').value)+1); const nf=parseInt(document.getElementById('frame').max);
    if(i>nf) i=0; setFrame(i); dirty=true; } }
  if(dirty){ dirty=false; renderer.render(scene,cam); }
  if(playing || dirty){ requestAnimationFrame(loop); } else { running=false; }}
requestRender();
</script></body></html>"""
