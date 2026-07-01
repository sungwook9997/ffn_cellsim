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
            pl.append({"name": L["name"], "kind": kind, "color": L.get("color", "#4aa3ff"),
                       "size": float(L.get("size", 2.0)), "n": int(v.shape[0]),
                       "b64": _b64(v, np.float32)})
        payload_scenes[sname] = pl
    allp = np.concatenate(all_pts, axis=0) if all_pts else np.zeros((1, 3), np.float32)
    lo = allp.min(0).tolist(); hi = allp.max(0).tolist()
    payload = {"title": title, "scenes": payload_scenes,
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
function clearScene(){ for(const o of current){ scene.remove(o); if(o.geometry)o.geometry.dispose(); if(o.material)o.material.dispose(); } current=[]; }

function showScene(name){
  clearScene();
  const layers=P.scenes[name]; const leg=[];
  for(const L of layers){
    if(L.kind==='plates'){
      for(const z of L.z){
        const g=new THREE.PlaneGeometry(L.half_xy*2,L.half_xy*2);
        const m=new THREE.MeshBasicMaterial({color:L.color,transparent:true,opacity:0.18,side:THREE.DoubleSide});
        const mesh=new THREE.Mesh(g,m); mesh.rotation.x=Math.PI/2; mesh.position.set(cx,z,cz);
        scene.add(mesh); current.push(mesh);
      }
      leg.push(['#888888','rigid plate']); continue;
    }
    const v=dec(L.b64); const geo=new THREE.BufferGeometry();
    geo.setAttribute('position',new THREE.BufferAttribute(v,3));
    if(L.kind==='lines'){
      const mat=new THREE.LineBasicMaterial({color:L.color,transparent:true,opacity:0.75});
      const o=new THREE.LineSegments(geo,mat); scene.add(o); current.push(o);
    }else{
      const mat=new THREE.PointsMaterial({color:L.color,size:parseFloat(document.getElementById('psize').value),sizeAttenuation:true});
      o=new THREE.Points(geo,mat); o.userData.pts=true; scene.add(o); current.push(o);
    }
    leg.push([L.color,L.name]);
  }
  document.getElementById('leg').innerHTML=leg.map(([c,n])=>`<div><span class="sw" style="background:${c}"></span>${n}</div>`).join('');
}
const sel=document.getElementById('scene');
for(const n of P.scene_names){ const o=document.createElement('option'); o.value=n; o.textContent=n; sel.appendChild(o); }
sel.onchange=()=>showScene(sel.value);
document.getElementById('psize').oninput=e=>{ for(const o of current) if(o.userData.pts) o.material.size=parseFloat(e.target.value); };
showScene(P.scene_names[0]);

addEventListener('resize',()=>{cam.aspect=innerWidth/innerHeight;cam.updateProjectionMatrix();renderer.setSize(innerWidth,innerHeight);});
(function loop(){requestAnimationFrame(loop);ctr.update();renderer.render(scene,cam);})();
</script></body></html>"""
