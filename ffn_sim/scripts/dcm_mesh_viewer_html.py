"""Self-contained interactive three.js HTML viewer for a DCM frames npz.

ONE double-clickable HTML that does all three things at once (PI 2026-06-30):
  * ROTATE / zoom / pan (OrbitControls)
  * PLAY all frames (play/pause + frame slider + fps) — animation
  * FREE CROSS-SECTION — a clipping plane you move with an axis selector + position slider,
    in real time, while rotating and playing.

Size is kept down for large N by (1) storing the face topology ONCE (constant across frames),
(2) quantising per-frame node positions to uint16 (half of float32), (3) per-vertex cell colours
as uint8 — all base64-embedded so the file opens offline by double-click (three.js from CDN).

    python -m ffn_sim.scripts.dcm_mesh_viewer_html --npz frames.npz --out viewer.html [--title ...]
"""
from __future__ import annotations

import argparse
import base64
import json

import numpy as np


def _b64(arr: np.ndarray) -> str:
    return base64.b64encode(np.ascontiguousarray(arr).tobytes()).decode("ascii")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--npz", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--title", default="DCM mesh viewer")
    ap.add_argument("--fps", type=int, default=6)
    args = ap.parse_args()

    d = np.load(args.npz, allow_pickle=True)
    frames = np.asarray(d["frames"], dtype=np.float64) * 1e6   # → µm
    faces = np.asarray(d["faces"], dtype=np.int64)
    cof = np.asarray(d["cof"], dtype=np.int64)
    F, N, _ = frames.shape
    M = faces.shape[0]

    # only render faces whose nodes are live (cof>=0); drop dormant
    live_face = cof[faces[:, 0]] >= 0
    faces = faces[live_face]
    M = faces.shape[0]

    # quantise positions to uint16 over the global bbox (per-axis scale/offset for JS dequant)
    lo = frames.reshape(-1, 3).min(0)
    hi = frames.reshape(-1, 3).max(0)
    span = np.maximum(hi - lo, 1e-9)
    q = np.clip(((frames - lo) / span) * 65535.0, 0, 65535).astype(np.uint16)  # (F,N,3)

    # per-vertex colour by cell (distinct vivid hue), uint8 RGB
    import colorsys
    cells = np.unique(cof[cof >= 0])
    cidx = {int(c): i for i, c in enumerate(cells)}
    palette = np.zeros((len(cells), 3), np.uint8)
    for i in range(len(cells)):
        r, g, b = colorsys.hsv_to_rgb((i * 7 % max(len(cells), 1)) / max(len(cells), 1), 0.62, 0.95)
        palette[i] = (int(r * 255), int(g * 255), int(b * 255))
    colors = np.zeros((N, 3), np.uint8)
    for n in range(N):
        c = int(cof[n])
        colors[n] = palette[cidx[c]] if c >= 0 else (120, 120, 120)

    meta = {
        "F": int(F), "N": int(N), "M": int(M),
        "lo": lo.tolist(), "span": span.tolist(),
        "fps": int(args.fps), "title": args.title,
    }
    payload = {
        "meta": meta,
        "faces_b64": _b64(faces.astype(np.uint32)),     # (M,3) uint32, constant
        "colors_b64": _b64(colors),                     # (N,3) uint8, constant
        "q_b64": _b64(q),                               # (F,N,3) uint16
    }
    html = _HTML.replace("/*__PAYLOAD__*/", json.dumps(payload))
    with open(args.out, "w") as f:
        f.write(html)
    import os
    mb = os.path.getsize(args.out) / 1e6
    print(f"wrote {args.out}  ({F} frames, {N} nodes, {M} faces, {len(cells)} cells, {mb:.1f} MB)")


_HTML = r"""<!DOCTYPE html><html><head><meta charset="utf-8">
<title>DCM viewer</title>
<style>
  body{margin:0;background:#0e1117;color:#ddd;font:13px sans-serif;overflow:hidden}
  #ui{position:absolute;top:8px;left:8px;z-index:10;background:rgba(20,24,33,.85);
      padding:10px 12px;border-radius:8px;max-width:340px}
  #ui label{display:inline-block;min-width:64px}
  #ui input[type=range]{width:170px;vertical-align:middle}
  #title{font-weight:bold;margin-bottom:6px}
  button{background:#2b6cb0;color:#fff;border:0;border-radius:5px;padding:4px 10px;cursor:pointer}
  .row{margin:5px 0}
</style></head><body>
<div id="ui">
  <div id="title"></div>
  <div class="row"><button id="play">▶ play</button>
     <span id="fnum"></span></div>
  <div class="row"><label>frame</label><input id="frame" type="range" min="0" value="0"></div>
  <div class="row"><label>fps</label><input id="fps" type="range" min="1" max="30" value="6"><span id="fpsv"></span></div>
  <hr style="border-color:#333">
  <div class="row"><label><input id="clipon" type="checkbox"> cut</label>
     <select id="axis"><option value="0">x</option><option value="1">y</option>
        <option value="2">z</option></select>
     <label style="min-width:auto"> flip</label><input id="flip" type="checkbox"></div>
  <div class="row"><label>cut pos</label><input id="clip" type="range" min="0" max="1000" value="500"></div>
  <div class="row" style="color:#888">drag=rotate · scroll=zoom · right-drag=pan</div>
</div>
<script src="https://cdnjs.cloudflare.com/ajax/libs/three.js/r128/three.min.js"></script>
<script src="https://cdn.jsdelivr.net/npm/three@0.128.0/examples/js/controls/OrbitControls.js"></script>
<script>
const P = /*__PAYLOAD__*/;
const m = P.meta;
function dec(b64,Ctor){const s=atob(b64);const u=new Uint8Array(s.length);
  for(let i=0;i<s.length;i++)u[i]=s.charCodeAt(i);return new Ctor(u.buffer);}
const faces = dec(P.faces_b64, Uint32Array);          // M*3
const colU  = dec(P.colors_b64, Uint8Array);          // N*3
const Q     = dec(P.q_b64, Uint16Array);              // F*N*3
const N=m.N, F=m.F, lo=m.lo, span=m.span;
document.getElementById('title').textContent = m.title;

// dequantised positions for frame f -> Float32Array(N*3) in µm
function framePos(f){
  const out=new Float32Array(N*3); const base=f*N*3;
  for(let i=0;i<N;i++){
    out[i*3]  = Q[base+i*3]  /65535*span[0]+lo[0];
    out[i*3+1]= Q[base+i*3+1]/65535*span[1]+lo[1];
    out[i*3+2]= Q[base+i*3+2]/65535*span[2]+lo[2];
  }
  return out;
}
// vertex colors (0..1)
const colF=new Float32Array(N*3);
for(let i=0;i<N*3;i++) colF[i]=colU[i]/255;

const renderer=new THREE.WebGLRenderer({antialias:true});
renderer.setSize(innerWidth,innerHeight); renderer.localClippingEnabled=true;
document.body.appendChild(renderer.domElement);
const scene=new THREE.Scene(); scene.background=new THREE.Color(0x0e1117);
const cam=new THREE.PerspectiveCamera(45,innerWidth/innerHeight,0.1,5000);
const cx=lo[0]+span[0]/2, cy=lo[1]+span[1]/2, cz=lo[2]+span[2]/2;
const R=Math.max(span[0],span[1],span[2]);
cam.position.set(cx+R*1.6,cy+R*1.1,cz+R*1.6);
const ctr=new THREE.OrbitControls(cam,renderer.domElement); ctr.target.set(cx,cy,cz); ctr.update();
scene.add(new THREE.AmbientLight(0xffffff,0.55));
const dl=new THREE.DirectionalLight(0xffffff,0.85); dl.position.set(1,1.5,2); scene.add(dl);

const geo=new THREE.BufferGeometry();
geo.setIndex(new THREE.BufferAttribute(faces,1));
geo.setAttribute('position',new THREE.BufferAttribute(framePos(0),3));
geo.setAttribute('color',new THREE.BufferAttribute(colF,3));
geo.computeVertexNormals();
const clipPlane=new THREE.Plane(new THREE.Vector3(1,0,0),0);
const mat=new THREE.MeshLambertMaterial({vertexColors:true,side:THREE.DoubleSide,
  clippingPlanes:[], clipShadows:true});
const mesh=new THREE.Mesh(geo,mat); scene.add(mesh);

// ---- UI ----
const elFrame=document.getElementById('frame'); elFrame.max=F-1;
const elPlay=document.getElementById('play'), elFnum=document.getElementById('fnum');
const elFps=document.getElementById('fps'), elFpsv=document.getElementById('fpsv');
const elClipOn=document.getElementById('clipon'), elClip=document.getElementById('clip');
const elAxis=document.getElementById('axis'), elFlip=document.getElementById('flip');
let cur=0, playing=false, fps=m.fps, last=0; elFps.value=fps; elFpsv.textContent=fps;
function setFrame(f){ cur=((f%F)+F)%F; geo.attributes.position.array.set(framePos(cur));
  geo.attributes.position.needsUpdate=true; geo.computeVertexNormals();
  elFrame.value=cur; elFnum.textContent='frame '+cur+'/'+(F-1); }
function updClip(){
  if(!elClipOn.checked){ mat.clippingPlanes=[]; mat.needsUpdate=true; return; }
  const ax=+elAxis.value, sgn=elFlip.checked?-1:1;
  const nrm=new THREE.Vector3(0,0,0); nrm.setComponent(ax,sgn);
  const t=elClip.value/1000; const val=lo[ax]+t*span[ax];
  clipPlane.normal.copy(nrm); clipPlane.constant=-sgn*val;
  mat.clippingPlanes=[clipPlane]; mat.needsUpdate=true;
}
elFrame.oninput=()=>{playing=false;elPlay.textContent='▶ play';setFrame(+elFrame.value);};
elPlay.onclick=()=>{playing=!playing;elPlay.textContent=playing?'❚❚ pause':'▶ play';};
elFps.oninput=()=>{fps=+elFps.value;elFpsv.textContent=fps;};
[elClipOn,elClip,elAxis,elFlip].forEach(e=>e.oninput=e.onchange=updClip);
setFrame(0); updClip();
addEventListener('resize',()=>{cam.aspect=innerWidth/innerHeight;cam.updateProjectionMatrix();
  renderer.setSize(innerWidth,innerHeight);});
function loop(t){requestAnimationFrame(loop); ctr.update();
  if(playing && t-last>1000/fps){last=t; setFrame(cur+1);}
  renderer.render(scene,cam);}
requestAnimationFrame(loop);
</script></body></html>"""


if __name__ == "__main__":
    main()
