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


def build_viewer(scenes: dict, out: str, title: str = "FF viewer") -> str:
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
                     "b64": _b64(v, np.float32)}
            if L.get("frames") is not None:                       # animated layer: swap positions per timestep
                fr = [np.asarray(f, dtype=np.float32).reshape(-1, 3) for f in L["frames"]]
                entry["frames"] = [_b64(f, np.float32) for f in fr]
                all_pts.append(fr[-1])                            # last frame bounds the box (elongated tip)
            pl.append(entry)
        payload_scenes[sname] = pl
    allp = np.concatenate(all_pts, axis=0) if all_pts else np.zeros((1, 3), np.float32)
    lo = allp.min(0).tolist(); hi = allp.max(0).tolist()
    n_frames = max([len(L.get("frames", [])) for ls in payload_scenes.values() for L in ls] + [0])
    payload = {"title": title, "scenes": payload_scenes, "n_frames": int(n_frames),
               "lo": lo, "hi": hi, "scene_names": list(scenes.keys())}
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
</style></head><body>
<div id="ui">
  <div id="title"></div>
  <div class="row"><label>scene</label><select id="scene"></select></div>
  <div class="row"><label>points</label><input id="psize" type="range" min="1" max="8" value="3" step="0.5"></div>
  <div class="row" id="animrow" style="display:none"><label>frame</label>
    <button id="play">▶</button><input id="frame" type="range" min="0" max="0" value="0" step="1" style="width:120px">
    <span id="fnum" class="hint"></span></div>
  <div class="row hint">drag = rotate · scroll = zoom · right-drag = pan</div>
</div>
<div id="leg"></div>
<script src="https://cdnjs.cloudflare.com/ajax/libs/three.js/r128/three.min.js"></script>
<script src="https://cdn.jsdelivr.net/npm/three@0.128.0/examples/js/controls/OrbitControls.js"></script>
<script>
const P = /*__PAYLOAD__*/;
document.getElementById('title').textContent = P.title;
function dec(b64){const s=atob(b64);const u=new Uint8Array(s.length);
  for(let i=0;i<s.length;i++)u[i]=s.charCodeAt(i);return new Float32Array(u.buffer);}

const renderer=new THREE.WebGLRenderer({antialias:true});
renderer.setSize(innerWidth,innerHeight); renderer.setPixelRatio(devicePixelRatio);
document.body.appendChild(renderer.domElement);
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
function clearScene(){ for(const o of current){ scene.remove(o); if(o.geometry)o.geometry.dispose(); if(o.material)o.material.dispose(); } current=[]; animated=[]; }

function setFrame(i){
  for(const a of animated){
    if(i>=a.frames.length) continue;
    const v=dec(a.frames[i]); a.obj.geometry.setAttribute('position',new THREE.BufferAttribute(v,3));
    a.obj.geometry.attributes.position.needsUpdate=true; a.obj.geometry.computeBoundingSphere();
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
      const mat=new THREE.LineBasicMaterial({color:L.color,transparent:true,opacity:0.75});
      o=new THREE.LineSegments(geo,mat); scene.add(o); current.push(o);
    }else{
      const mat=new THREE.PointsMaterial({color:L.color,size:parseFloat(document.getElementById('psize').value),sizeAttenuation:true});
      o=new THREE.Points(geo,mat); o.userData.pts=true; scene.add(o); current.push(o);
    }
    if(L.frames){ animated.push({obj:o, frames:L.frames}); }
    leg.push([L.color,L.name]);
  }
  const hasAnim = animated.length>0;
  document.getElementById('animrow').style.display = hasAnim ? 'flex' : 'none';
  if(hasAnim){ const nf=Math.max(...animated.map(a=>a.frames.length)); document.getElementById('frame').max=nf-1; setFrame(0); }
  document.getElementById('leg').innerHTML=leg.map(([c,n])=>`<div><span class="sw" style="background:${c}"></span>${n}</div>`).join('');
}
const sel=document.getElementById('scene');
for(const n of P.scene_names){ const o=document.createElement('option'); o.value=n; o.textContent=n; sel.appendChild(o); }
sel.onchange=()=>showScene(sel.value);
document.getElementById('psize').oninput=e=>{ for(const o of current) if(o.userData.pts) o.material.size=parseFloat(e.target.value); };

let playing=false, tacc=0;
document.getElementById('play').onclick=()=>{ playing=!playing; document.getElementById('play').textContent=playing?'⏸':'▶'; };
document.getElementById('frame').oninput=e=>{ playing=false; document.getElementById('play').textContent='▶'; setFrame(parseInt(e.target.value)); };
showScene(P.scene_names[0]);

addEventListener('resize',()=>{cam.aspect=innerWidth/innerHeight;cam.updateProjectionMatrix();renderer.setSize(innerWidth,innerHeight);});
let last=performance.now();
(function loop(){requestAnimationFrame(loop);ctr.update();
  const now=performance.now(); const dt=(now-last)/1000; last=now;
  if(playing && animated.length){ tacc+=dt; if(tacc>0.06){ tacc=0;
    let i=(parseInt(document.getElementById('frame').value)+1); const nf=parseInt(document.getElementById('frame').max);
    if(i>nf) i=0; setFrame(i); } }
  renderer.render(scene,cam);})();
</script></body></html>"""
