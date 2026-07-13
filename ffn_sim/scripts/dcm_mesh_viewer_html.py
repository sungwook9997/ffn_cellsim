"""Self-contained interactive three.js HTML viewer for a DCM frames npz.

ONE double-clickable HTML that does everything at once (PI 2026-06-30):
  * ROTATE / zoom / pan (OrbitControls)
  * PLAY all frames (play/pause + frame slider + fps) — animation
  * FREE CROSS-SECTION with a clip plane (axis selector + position slider + flip),
    in real time, while rotating and playing.

The headline fix (PI: "단면이 제대로 안 보인다"): a plain clip plane cuts each closed
cell shell HOLLOW (three.js r128 clips per-fragment with no back-fill, so you see into
empty shells). This viewer offers THREE cross-section modes, all of which show a SOLID,
clean cross-section:

  * peel  (DEFAULT) — hide every cell whose centroid is beyond the plane; only WHOLE
    closed cells are drawn, so every drawn cell is a solid shell (never sliced → never
    hollow). Per-cell centroids are recomputed per frame from `cof`.
  * slab  — like peel but keep only cells whose centroid lies within ±thickness of the
    plane (a one-cell-thick cross-sectional layer). Adds a thickness slider.
  * cut   — true geometric clip WITH a stencil-buffer cap, so each cut closed cell shows
    a SOLID filled cross-section instead of a hollow shell. Uses a single whole-scene
    stencil pass (back-incr / front-decr parity) + one neutral cap plane (see JS comment
    for why a single pass is correct and per-cell capping is avoided).

Size is kept down for large N by (1) storing the face topology ONCE (constant across
frames), (2) quantising per-frame node positions to uint16, (3) per-vertex cell colours
as uint8, (4) a per-face compact cell-index array (for O(M) peel/slab) — all base64
embedded so the file opens offline by double-click (three.js r128 from CDN).

    python -m ffn_sim.scripts.dcm_mesh_viewer_html --npz frames.npz --out viewer.html \
        [--title ...] [--fps N]
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
    ap.add_argument("--r-nuc-factor", type=float, default=0.0, dest="r_nuc_factor",
                    help="if >0, render a NUCLEUS sphere per cell at its per-frame centroid, radius = "
                         "factor*cell_radius (PROPORTIONAL / N-C-ratio look). Cortex/membrane shells become "
                         "translucent so the nucleus (and cytoplasm interior) are visible inside. WARNING: "
                         "proportional drawing makes a real cell-SIZE gradient look like nuclear compression — "
                         "prefer --r-nuc-abs to match the sim (which uses ONE fixed R_nuc for every cell).")
    ap.add_argument("--r-nuc-abs", type=float, default=0.0, dest="r_nuc_abs",
                    help="if >0, render every nucleus at this FIXED absolute radius [µm] (e.g. 1.88 = the sim's "
                         "R_nuc = 0.25*7.5µm). Faithful to the sim (uniform R_nuc for all cells); overrides "
                         "--r-nuc-factor for sizing. Use this for honest morphology; the proportional mode is a "
                         "display convention only.")
    ap.add_argument("--frame-stride", type=int, default=1, dest="frame_stride",
                    help="keep every Nth trajectory frame (default 1 = all). Use >1 for large-N runs to "
                         "bound the embedded HTML size (e.g. N=2000 has ~320k nodes → stride 3 keeps ~8 frames).")
    ap.add_argument("--max-frames", type=int, default=0, dest="max_frames",
                    help="if >0, cap the number of kept frames to this many (evenly spaced across the "
                         "trajectory, always including the last). Applied after --frame-stride.")
    args = ap.parse_args()

    d = np.load(args.npz, allow_pickle=True)
    frames = np.asarray(d["frames"], dtype=np.float64) * 1e6   # → µm
    faces = np.asarray(d["faces"], dtype=np.int64)
    cof = np.asarray(d["cof"], dtype=np.int64)
    F0 = frames.shape[0]
    # frame subsampling (bounds embedded HTML size for large-N runs). `keep` indexes the ORIGINAL
    # trajectory; every per-frame array (frames + mech overlays + bondpairs) is indexed through it so
    # the animation stays self-consistent. Always keep the last frame (the settled/endpoint state).
    keep = np.arange(0, F0, max(1, int(args.frame_stride)))
    if args.max_frames and len(keep) > args.max_frames:
        keep = np.unique(np.linspace(0, F0 - 1, int(args.max_frames)).round().astype(int))
    if keep[-1] != F0 - 1:
        keep = np.append(keep, F0 - 1)
    frames = frames[keep]
    F, N, _ = frames.shape

    # only render faces whose nodes are live (cof>=0); drop dormant
    live_face = cof[faces[:, 0]] >= 0
    faces = faces[live_face]
    M = faces.shape[0]

    # quantise positions to uint16 over the global bbox (per-axis scale/offset for JS dequant)
    lo = frames.reshape(-1, 3).min(0)
    hi = frames.reshape(-1, 3).max(0)
    span = np.maximum(hi - lo, 1e-9)
    q = np.clip(((frames - lo) / span) * 65535.0, 0, 65535).astype(np.uint16)  # (F,N,3)

    # compact cell indices: live cells → 0..(C-1). per-NODE compact id and per-FACE compact id.
    cells = np.unique(cof[cof >= 0])
    C = int(len(cells))
    remap = {int(c): i for i, c in enumerate(cells)}
    node_cell = np.full(N, -1, np.int32)
    for n in range(N):
        c = int(cof[n])
        if c >= 0:
            node_cell[n] = remap[c]
    # a face's three nodes all share one cell (cells don't share nodes); use vertex 0.
    face_cell = node_cell[faces[:, 0]].astype(np.int32)  # (M,) compact cell id per face

    # per-vertex vivid colour by cell, uint8 RGB
    import colorsys
    palette = np.zeros((C, 3), np.uint8)
    for i in range(C):
        # spread hues with a step coprime-ish to C for vivid neighbour contrast
        h = (i * 7 % max(C, 1)) / max(C, 1)
        r, g, b = colorsys.hsv_to_rgb(h, 0.66, 0.97)
        palette[i] = (int(r * 255), int(g * 255), int(b * 255))
    colors = np.zeros((N, 3), np.uint8)
    for n in range(N):
        ci = int(node_cell[n])
        colors[n] = palette[ci] if ci >= 0 else (120, 120, 120)

    # Per-cell FIXED nucleus radius [µm]: the tightest cell inradius over the WHOLE trajectory,
    # clamped to R_nuc. A stiff nucleus is ~constant — the old PER-FRAME clamp made the drawn
    # nucleus GROW as cells relaxed. Bake ONE radius per cell so it stays constant across the
    # animation AND never pokes through the (compressed) membrane in any frame.
    rnuc_abs = float(args.r_nuc_abs)
    rnuc_cell = np.full(C, rnuc_abs, np.float32)
    if rnuc_abs > 0.0:
        cell_nodes = [[] for _ in range(C)]
        for n in range(N):
            ci = int(node_cell[n])
            if ci >= 0:
                cell_nodes[ci].append(n)
        for i in range(C):
            nb = cell_nodes[i]
            if len(nb) < 4:
                continue
            nb = np.asarray(nb)
            rmin = 1.0e30
            for fi in range(F):
                p = frames[fi, nb]
                dr = float(np.linalg.norm(p - p.mean(0), axis=1).min())
                if dr < rmin:
                    rmin = dr
            rnuc_cell[i] = min(rnuc_abs, rmin * 0.99)

    # Per-NODE overlays the viewer colours by (REAL sim mechanics when the run saved them):
    #   stress   = virial/Cauchy von-Mises σ_vm per node [Pa] (svm) → the REAL stress heatmap
    #   pressure = hydrostatic pressure per node [Pa] (spress, turgor-dominated)
    #   residual = |net residual force| per node [N] (fmag; ∝ velocity, NOT a stress — kept, honest name)
    #   junction = cadherin bond count per node (nbond)
    # plus the actual bond node-pairs per frame (bondpairs) → draw the JUNCTIONS as line segments.
    def _norm_pctl(A, hi=99.0):
        A = A.astype(np.float32)
        top = float(np.percentile(A, hi)) if A.size else 1.0
        return (np.clip(A / top, 0.0, 1.0) if top > 0 else A), top

    have_mech = "fmag" in d.files and "nbond" in d.files
    have_stress = "svm" in d.files                          # real virial stress (new runs); else fall back to fmag
    if have_mech:
        if have_stress:
            stress, stress_pk = _norm_pctl(np.asarray(d["svm"])[keep])     # von-Mises σ_vm [Pa]
        else:
            stress, stress_pk = _norm_pctl(np.asarray(d["fmag"])[keep])    # legacy: residual as "stress"
        residual, resid_pk = _norm_pctl(np.asarray(d["fmag"])[keep])       # |net residual force| [N]
        junctionN, _ = _norm_pctl(np.asarray(d["nbond"])[keep])
        if "spress" in d.files:
            pressure, press_pk = _norm_pctl(np.abs(np.asarray(d["spress"]))[keep])
        else:
            pressure, press_pk = np.zeros((F, N), np.float32), 0.0
    else:                                                      # legacy npz (no mech): flat overlays
        stress = residual = junctionN = pressure = np.zeros((F, N), np.float32)
        stress_pk = resid_pk = press_pk = 0.0
    # bond pairs → concatenated index array + per-frame offsets (for LineSegments)
    bond_idx_list, bond_off = [], [0]
    if "bondpairs" in d.files:
        bp_all = d["bondpairs"]
        for fi in keep:
            bp = np.asarray(bp_all[fi], dtype=np.int32).reshape(-1, 2)
            bond_idx_list.append(bp)
            bond_off.append(bond_off[-1] + bp.shape[0])
        bond_idx = (np.concatenate(bond_idx_list, axis=0) if bond_idx_list else
                    np.zeros((0, 2), np.int32))
    else:
        bond_idx = np.zeros((0, 2), np.int32)
        bond_off = [0] * (F + 1)
    bond_off = np.asarray(bond_off, np.int32)

    meta = {
        "F": int(F), "N": int(N), "M": int(M), "C": C,
        "lo": lo.tolist(), "span": span.tolist(),
        "fps": int(args.fps), "title": args.title,
        "rNuc": float(args.r_nuc_factor),
        "rNucAbs": float(args.r_nuc_abs),
    }
    payload = {
        "meta": meta,
        "faces_b64": _b64(faces.astype(np.uint32)),      # (M,3) uint32, constant
        "facecell_b64": _b64(face_cell),                 # (M,)  int32  compact cell per face
        "nodecell_b64": _b64(node_cell),                 # (N,)  int32  compact cell per node
        "colors_b64": _b64(colors),                      # (N,3) uint8, constant
        "palette_b64": _b64(palette),                    # (C,3) uint8 per-cell colour (nucleus tint)
        "rnuc_b64": _b64(rnuc_cell),                      # (C,) float32 per-cell FIXED nucleus radius [µm]
        "stress_b64": _b64(stress),                      # (F,N) float32 per-node virial von-Mises STRESS [0,1]
        "residual_b64": _b64(residual),                  # (F,N) float32 per-node RESIDUAL force [0,1]
        "pressure_b64": _b64(pressure),                  # (F,N) float32 per-node hydrostatic pressure |·| [0,1]
        "junctionN_b64": _b64(junctionN),                # (F,N) float32 per-node cadherin bond count [0,1]
        "bondidx_b64": _b64(bond_idx.astype(np.int32)),  # (Btot,2) int32 bonded node pairs, all frames
        "bondoff_b64": _b64(bond_off.astype(np.int32)),  # (F+1,) int32 per-frame slice offsets
        "hasmech": 1 if have_mech else 0,
        "hasstress": 1 if have_stress else 0,            # 1 = real virial σ_vm; 0 = fmag fallback in the stress slot
        "stress_pk": float(stress_pk), "resid_pk": float(resid_pk), "press_pk": float(press_pk),
        "q_b64": _b64(q),                                # (F,N,3) uint16
    }
    html = _HTML.replace("/*__PAYLOAD__*/", json.dumps(payload))
    with open(args.out, "w") as f:
        f.write(html)
    import os
    mb = os.path.getsize(args.out) / 1e6
    print(f"wrote {args.out}  ({F} frames, {N} nodes, {M} faces, {C} cells, {mb:.1f} MB)")


_HTML = r"""<!DOCTYPE html><html><head><meta charset="utf-8">
<title>DCM viewer</title>
<style>
  body{margin:0;background:#0e1117;color:#ddd;font:13px sans-serif;overflow:hidden}
  #ui{position:absolute;top:8px;left:8px;z-index:10;background:rgba(20,24,33,.88);
      padding:10px 12px;border-radius:8px;max-width:340px}
  #ui label{display:inline-block;min-width:64px}
  #ui input[type=range]{width:160px;vertical-align:middle}
  #title{font-weight:bold;margin-bottom:6px}
  button{background:#2b6cb0;color:#fff;border:0;border-radius:5px;padding:4px 10px;cursor:pointer}
  .row{margin:5px 0}
  select{background:#1a1f2b;color:#ddd;border:1px solid #333;border-radius:4px;padding:2px}
  .hint{color:#888}
</style></head><body>
<div id="ui">
  <div id="title"></div>
  <div class="row"><button id="play">▶ play</button>
     <span id="fnum"></span></div>
  <div class="row"><label>frame</label><input id="frame" type="range" min="0" value="0"></div>
  <div class="row"><label>fps</label><input id="fps" type="range" min="1" max="30" value="6"><span id="fpsv"></span></div>
  <div class="row" id="nucrow"><label><input id="nucon" type="checkbox"> nucleus</label>
     <span class="hint" style="margin-left:8px">off → cells opaque (see stress + contacts)</span></div>
  <div class="row"><label>colour</label>
     <select id="cmode">
        <option value="stress" selected>stress (virial σ_vm)</option>
        <option value="pressure">pressure (hydrostatic)</option>
        <option value="residual">residual force</option>
        <option value="junction">junction density</option>
        <option value="cell">cell id</option>
     </select></div>
  <div class="row" id="junrow"><label><input id="junon" type="checkbox" checked> show cadherin bonds</label></div>
  <div class="row hint" id="cmodehint"></div>
  <hr style="border-color:#333">
  <div class="row"><label><input id="clipon" type="checkbox"> section</label>
     <select id="mode">
        <option value="peel" selected>peel (solid)</option>
        <option value="slab">slab (layer)</option>
        <option value="cut">cut (capped)</option>
     </select></div>
  <div class="row" id="axisrow" style="display:none"><label>axis</label>
     <select id="axis"><option value="0">x</option><option value="1">y</option>
        <option value="2">z</option></select>
     <label style="min-width:auto"> flip</label><input id="flip" type="checkbox"></div>
  <div class="row" id="cutrow" style="display:none"><label>cut pos</label><input id="clip" type="range" min="0" max="1000" value="500"></div>
  <div class="row" id="slabrow" style="display:none"><label>thickness</label>
     <input id="thick" type="range" min="1" max="500" value="120"><span id="thickv"></span></div>
  <div class="row hint">drag=rotate · scroll=zoom · right-drag=pan</div>
  <div class="row hint" id="modehint"></div>
</div>
<script src="https://cdnjs.cloudflare.com/ajax/libs/three.js/r128/three.min.js"></script>
<script src="https://cdn.jsdelivr.net/npm/three@0.128.0/examples/js/controls/OrbitControls.js"></script>
<script>
const P = /*__PAYLOAD__*/;
const m = P.meta;
function dec(b64,Ctor){const s=atob(b64);const u=new Uint8Array(s.length);
  for(let i=0;i<s.length;i++)u[i]=s.charCodeAt(i);return new Ctor(u.buffer);}
const facesAll = dec(P.faces_b64, Uint32Array);       // M*3  (constant topology)
const faceCell = dec(P.facecell_b64, Int32Array);     // M    compact cell id per face
const nodeCell = dec(P.nodecell_b64, Int32Array);     // N    compact cell id per node
const colU     = dec(P.colors_b64, Uint8Array);       // N*3
const Q        = dec(P.q_b64, Uint16Array);           // F*N*3
const N=m.N, F=m.F, M=m.M, C=m.C, lo=m.lo, span=m.span;
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

const renderer=new THREE.WebGLRenderer({antialias:true,stencil:true});
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

// ---------- geometry ----------
const geo=new THREE.BufferGeometry();
// peel/slab rewrite the index in place; allocate the full-size index buffer up front
// (Uint32) and drive how many tris are drawn with setDrawRange.
const idxBuf = new Uint32Array(facesAll.length);
idxBuf.set(facesAll);
geo.setIndex(new THREE.BufferAttribute(idxBuf,1));
geo.setAttribute('position',new THREE.BufferAttribute(framePos(0),3));
geo.setAttribute('color',new THREE.BufferAttribute(colF,3));

// ---- per-NODE overlays: cell-id | stress (virial σ_vm) | pressure | residual force | junction ----
const stressArr   = dec(P.stress_b64, Float32Array);     // F*N in [0,1] (virial von Mises, or fmag fallback)
const residualArr = P.residual_b64 ? dec(P.residual_b64, Float32Array) : stressArr;
const pressureArr = P.pressure_b64 ? dec(P.pressure_b64, Float32Array) : stressArr;
const junctionArr = dec(P.junctionN_b64, Float32Array);  // F*N in [0,1]
function heat(t){ t=Math.max(0,Math.min(1,t));           // blue→cyan→green→yellow→red
  return [Math.max(0,Math.min(1,1.5-Math.abs(4*t-3))),
          Math.max(0,Math.min(1,1.5-Math.abs(4*t-2))),
          Math.max(0,Math.min(1,1.5-Math.abs(4*t-1)))]; }
function overlayArr(mode){ return mode==='pressure'?pressureArr : mode==='residual'?residualArr
  : mode==='junction'?junctionArr : stressArr; }
function applyColors(f){
  const mode=document.getElementById('cmode').value;
  if(mode==='cell'){ for(let i=0;i<N*3;i++) colF[i]=colU[i]/255; }
  else {
    const arr=overlayArr(mode), base=f*N;
    for(let n=0;n<N;n++){
      if(nodeCell[n]<0){ colF[n*3]=colF[n*3+1]=colF[n*3+2]=0.47; continue; }
      const g=heat(arr[base+n]); colF[n*3]=g[0]; colF[n*3+1]=g[1]; colF[n*3+2]=g[2]; }
  }
  geo.attributes.color.array.set(colF); geo.attributes.color.needsUpdate=true;
}

// ---- cadherin JUNCTION lines: draw each active bond as a segment between its two nodes ----
const bondIdx = dec(P.bondidx_b64, Int32Array);          // (Btot*2)
const bondOff = dec(P.bondoff_b64, Int32Array);          // (F+1)
let jLine=null;
if(P.hasmech && bondOff.length > 1){
  const maxSeg = (bondOff[F]||0);
  const jpos=new Float32Array(Math.max(maxSeg,1)*2*3);
  const jgeo=new THREE.BufferGeometry();
  jgeo.setAttribute('position', new THREE.BufferAttribute(jpos,3));
  const jmat=new THREE.LineBasicMaterial({color:0xffe64d});   // bright yellow bonds
  jLine=new THREE.LineSegments(jgeo, jmat); jLine.visible=false; scene.add(jLine);
}
function applyJunctions(f){
  if(!jLine) return;
  const on=document.getElementById('junon').checked;
  jLine.visible=on;
  if(!on) return;
  const a=bondOff[f], b=bondOff[f+1], jpos=jLine.geometry.attributes.position.array;
  let w=0;
  for(let k=a;k<b;k++){
    const i=bondIdx[k*2], j=bondIdx[k*2+1];
    jpos[w++]=curPos[i*3]; jpos[w++]=curPos[i*3+1]; jpos[w++]=curPos[i*3+2];
    jpos[w++]=curPos[j*3]; jpos[w++]=curPos[j*3+1]; jpos[w++]=curPos[j*3+2];
  }
  jLine.geometry.setDrawRange(0, (b-a)*2);
  jLine.geometry.attributes.position.needsUpdate=true;
}
geo.computeVertexNormals();

// ---------- clip plane (shared by all clip-using materials) ----------
const clipPlane=new THREE.Plane(new THREE.Vector3(1,0,0),0);
const clipArr=[];   // empty when section OFF or in peel/slab; [clipPlane] only in cut mode

// (1) visible cell surface — the normal look. Reads clipArr (only populated in cut mode).
const mat=new THREE.MeshLambertMaterial({vertexColors:true,side:THREE.DoubleSide,
  clippingPlanes:clipArr, clipShadows:true});
const nucOn = (m.rNuc>0 || m.rNucAbs>0);
if(nucOn){ mat.transparent=true; mat.opacity=0.42; mat.depthWrite=false; } // translucent membrane/cortex → the distinct blue nucleus reads clearly inside while the cell boundary stays visible
const mesh=new THREE.Mesh(geo,mat); scene.add(mesh);

// ---------- NUCLEUS compartment: one instanced sphere per cell at its centroid, radius R_nuc ----------
const palette = nucOn ? dec(P.palette_b64, Uint8Array) : null;   // (C,3) uint8
const rnucFixed = nucOn ? dec(P.rnuc_b64, Float32Array) : null;  // (C,) per-cell FIXED nucleus radius [µm] — constant across frames
let nucMesh=null;
if(nucOn){
  const nucMat=new THREE.MeshLambertMaterial({clippingPlanes:clipArr, clipShadows:true});
  nucMesh=new THREE.InstancedMesh(new THREE.IcosahedronGeometry(1,2), nucMat, C);
  nucMesh.instanceColor=new THREE.InstancedBufferAttribute(new Float32Array(C*3),3);
  scene.add(nucMesh);
}

// ---------- STENCIL CAPPING (single whole-scene pass) for "cut" mode ----------
// We use ONE whole-scene stencil pass + ONE neutral cap plane, NOT 400 per-cell groups.
// The back-incr / front-decr parity is additive across disjoint closed manifolds: with
// GL_INCR_WRAP on back faces and GL_DECR_WRAP on front faces (depth test passing, i.e.
// in front of the clip plane), the stencil buffer ends up nonzero exactly on pixels that
// are interior to SOME closed solid at the cut plane, and zero elsewhere. So a single pass
// over the whole merged geometry caps every cell correctly — no need to isolate cells.
// Per-cell cap colour is intentionally NOT done: the stencil mask is a scalar with no cell
// identity, so colouring caps per-cell would force 400 per-cell passes (the slow path).
// A neutral grey cap reads cleanly as "interior". 8-bit stencil INCR_WRAP wraps at 256 —
// safe here (you'd need 256 back-faces stacked at one pixel; ~400 disjoint cells never do).
const matBack=new THREE.MeshBasicMaterial();
matBack.depthWrite=false; matBack.colorWrite=false;
matBack.side=THREE.BackSide;                 // render back faces of the solids
matBack.clippingPlanes=clipArr;
matBack.stencilWrite=true;
matBack.stencilFunc=THREE.AlwaysStencilFunc;
matBack.stencilFail=THREE.IncrementWrapStencilOp;
matBack.stencilZFail=THREE.IncrementWrapStencilOp;
matBack.stencilZPass=THREE.IncrementWrapStencilOp;
const meshBack=new THREE.Mesh(geo,matBack); meshBack.renderOrder=1; scene.add(meshBack);

const matFront=new THREE.MeshBasicMaterial();
matFront.depthWrite=false; matFront.colorWrite=false;
matFront.side=THREE.FrontSide;               // render front faces of the solids
matFront.clippingPlanes=clipArr;
matFront.stencilWrite=true;
matFront.stencilFunc=THREE.AlwaysStencilFunc;
matFront.stencilFail=THREE.DecrementWrapStencilOp;
matFront.stencilZFail=THREE.DecrementWrapStencilOp;
matFront.stencilZPass=THREE.DecrementWrapStencilOp;
const meshFront=new THREE.Mesh(geo,matFront); meshFront.renderOrder=2; scene.add(meshFront);

// the cap: a big quad living ON the clip plane, drawn only where stencil != 0, then the
// stencil is cleared (zero) so it doesn't leak into later frames.
const capGeo=new THREE.PlaneGeometry(1,1);
const capMat=new THREE.MeshLambertMaterial({color:0x9099a8,side:THREE.DoubleSide});
capMat.stencilWrite=true;
capMat.stencilRef=0;
capMat.stencilFunc=THREE.NotEqualStencilFunc;        // draw where stencil != 0 (interior)
capMat.stencilFail=THREE.ReplaceStencilOp;           // and reset to ref(0) as we go
capMat.stencilZFail=THREE.ReplaceStencilOp;
capMat.stencilZPass=THREE.ReplaceStencilOp;
const capMesh=new THREE.Mesh(capGeo,capMat); capMesh.renderOrder=3; scene.add(capMesh);
// the cap quad spans the whole bbox; size it generously and reorient onto the plane.
const capSize=R*3;
function orientCap(ax,sgn,val){
  // place + orient the cap quad on the clip plane
  capMesh.scale.set(capSize,capSize,1);
  const px=lo[0]+span[0]/2, py=lo[1]+span[1]/2, pz=lo[2]+span[2]/2;
  if(ax===0){ capMesh.position.set(val,py,pz); capMesh.rotation.set(0,Math.PI/2,0); }
  else if(ax===1){ capMesh.position.set(px,val,pz); capMesh.rotation.set(Math.PI/2,0,0); }
  else { capMesh.position.set(px,py,val); capMesh.rotation.set(0,0,0); }
}

// ---------- per-cell centroids (recomputed when frame OR cut moves) ----------
const cent=new Float32Array(C*3);
const cnt=new Float32Array(C);
function computeCentroids(pos){
  cent.fill(0); cnt.fill(0);
  for(let n=0;n<N;n++){
    const ci=nodeCell[n]; if(ci<0) continue;
    cent[ci*3]  += pos[n*3];
    cent[ci*3+1]+= pos[n*3+1];
    cent[ci*3+2]+= pos[n*3+2];
    cnt[ci]++;
  }
  for(let c=0;c<C;c++){ const k=cnt[c]||1; cent[c*3]/=k; cent[c*3+1]/=k; cent[c*3+2]/=k; }
}
// per-cell radius (mean node distance to centroid). Nucleus size is EITHER a fixed absolute
// radius rNucAbs [µm] (faithful to the sim, which uses ONE global R_nuc for every cell) OR the
// legacy per-cell rNuc·cellR proportional render (N/C-ratio look — do NOT read it as the sim's
// nucleus; proportional drawing makes a real cell-SIZE gradient look like nuclear compression).
const cellR=new Float32Array(C);
const cellRmin=new Float32Array(C);   // per-cell inradius (min node→centroid) — the nucleus can't exceed it
function computeCellR(pos){
  cellR.fill(0); cellRmin.fill(1e30); const k=new Float32Array(C);
  for(let n=0;n<N;n++){ const ci=nodeCell[n]; if(ci<0) continue;
    const dx=pos[n*3]-cent[ci*3], dy=pos[n*3+1]-cent[ci*3+1], dz=pos[n*3+2]-cent[ci*3+2];
    const dr=Math.sqrt(dx*dx+dy*dy+dz*dz);
    cellR[ci]+=dr; if(dr<cellRmin[ci]) cellRmin[ci]=dr; k[ci]++; }
  for(let c=0;c<C;c++){ cellR[c]/=(k[c]||1); if(cellRmin[c]>1e29) cellRmin[c]=cellR[c]; }
}
const _nm=new THREE.Matrix4(), _nc=new THREE.Color();
function updateNucleus(){
  if(!nucMesh) return;
  computeCentroids(curPos); computeCellR(curPos);
  const on=elClipOn.checked, mode=elMode.value;
  for(let c=0;c<C;c++){
    // The nucleus is a stiff ~constant body. Use the BAKED per-cell fixed radius (tightest inradius
    // over the whole trajectory) so it never grows/pulses across frames and never pokes through the
    // (compressed) membrane. Proportional legacy path unchanged.
    let rn = (m.rNucAbs>0.0) ? rnucFixed[c] : m.rNuc*cellR[c];
    if(on && (mode==='peel'||mode==='slab') && !keep[c]) rn=0;  // hide nuclei of peeled/hidden cells
    _nm.makeScale(rn,rn,rn); _nm.setPosition(cent[c*3],cent[c*3+1],cent[c*3+2]);
    nucMesh.setMatrixAt(c,_nm);
    _nc.setRGB(0.16,0.34,0.88);  // DAPI-like blue nucleus — a distinct compartment colour, NOT a shade of the cell membrane, so nucleus vs cell boundary read apart
    nucMesh.setColorAt(c,_nc);
  }
  nucMesh.instanceMatrix.needsUpdate=true;
  if(nucMesh.instanceColor) nucMesh.instanceColor.needsUpdate=true;
}
// per-cell keep flag for peel/slab
const keep=new Uint8Array(C);
function rebuildIndex(){
  // walk faces, write the visible subset into idxBuf, set draw range. O(M).
  let w=0;
  for(let f=0;f<M;f++){
    if(keep[faceCell[f]]){
      idxBuf[w]   = facesAll[f*3];
      idxBuf[w+1] = facesAll[f*3+1];
      idxBuf[w+2] = facesAll[f*3+2];
      w+=3;
    }
  }
  geo.index.needsUpdate=true;
  geo.setDrawRange(0,w);
}
function showAllIndex(){
  idxBuf.set(facesAll);
  geo.index.needsUpdate=true;
  geo.setDrawRange(0,facesAll.length);
}

// ---------- UI ----------
const elFrame=document.getElementById('frame'); elFrame.max=F-1;
const elPlay=document.getElementById('play'), elFnum=document.getElementById('fnum');
const elFps=document.getElementById('fps'), elFpsv=document.getElementById('fpsv');
const elClipOn=document.getElementById('clipon'), elClip=document.getElementById('clip');
const elAxis=document.getElementById('axis'), elFlip=document.getElementById('flip');
const elMode=document.getElementById('mode'), elModeHint=document.getElementById('modehint');
const elThick=document.getElementById('thick'), elThickv=document.getElementById('thickv');
const elSlabRow=document.getElementById('slabrow');
let cur=0, playing=false, fps=m.fps, last=0; elFps.value=fps; elFpsv.textContent=fps;
let curPos=framePos(0);

const HINTS={
  peel:'peel: whole near-side cells only — solid shells, never hollow.',
  slab:'slab: only cells within ±thickness of the plane — a one-cell layer.',
  cut :'cut: true clip with stencil cap — solid grey interior face.',
};

function setFrame(f){
  cur=((f%F)+F)%F;
  curPos=framePos(cur);
  geo.attributes.position.array.set(curPos);
  geo.attributes.position.needsUpdate=true;
  geo.computeVertexNormals();
  elFrame.value=cur; elFnum.textContent='frame '+cur+'/'+(F-1);
  applySection();   // peel/slab depend on the (moved) positions
  updateNucleus();  // nucleus spheres track the (moved) centroids
  applyColors(cur); // recolour by the current overlay (cell / stress / junction)
  applyJunctions(cur); // redraw the active cadherin bonds for this frame
}

function applySection(){
  const on=elClipOn.checked;
  const mode=elMode.value;
  const ax=+elAxis.value, sgn=elFlip.checked?-1:1;
  const t=elClip.value/1000; const val=lo[ax]+t*span[ax];
  // toggle which meshes are active
  const cutOn = on && mode==='cut';
  meshBack.visible=meshFront.visible=capMesh.visible=cutOn;

  if(!on){
    // section off: show everything, no clip, no cap
    clipArr.length=0; mat.needsUpdate=true;
    showAllIndex();
    return;
  }

  if(mode==='cut'){
    // geometric clip on the surface material + stencil cap pass
    clipArr.length=0;
    clipArr.push(clipPlane);
    const nrm=new THREE.Vector3(0,0,0); nrm.setComponent(ax,sgn);
    clipPlane.normal.copy(nrm); clipPlane.constant=-sgn*val;
    mat.needsUpdate=true;
    showAllIndex();
    orientCap(ax,sgn,val);
  } else {
    // peel / slab: no clip plane; select whole cells by centroid, rebuild index
    clipArr.length=0; mat.needsUpdate=true;
    computeCentroids(curPos);
    if(mode==='peel'){
      // keep cells on the near side of the plane (centroid*sgn <= val*sgn)
      for(let c=0;c<C;c++) keep[c] = (cent[c*3+ax]*sgn <= val*sgn) ? 1 : 0;
    } else { // slab
      const half=(elThick.value/1000)*span[ax]*0.5;
      for(let c=0;c<C;c++) keep[c] = (Math.abs(cent[c*3+ax]-val) <= half) ? 1 : 0;
    }
    rebuildIndex();
  }
}

function updMode(){
  const on=elClipOn.checked;
  document.getElementById('axisrow').style.display = on ? '' : 'none';   // section detail only when ON
  document.getElementById('cutrow').style.display  = on ? '' : 'none';
  elSlabRow.style.display = (on && elMode.value==='slab') ? '' : 'none';
  elModeHint.textContent = on ? HINTS[elMode.value] : '';
  applySection(); updateNucleus();
}

elFrame.oninput=()=>{playing=false;elPlay.textContent='▶ play';setFrame(+elFrame.value);};
elPlay.onclick=()=>{playing=!playing;elPlay.textContent=playing?'❚❚ pause':'▶ play';};
elFps.oninput=()=>{fps=+elFps.value;elFpsv.textContent=fps;};
const secRefresh=()=>{applySection();updateNucleus();};
elThick.oninput=()=>{elThickv.textContent=(elThick.value/10).toFixed(0)+'%';secRefresh();};
[elClipOn,elClip,elAxis,elFlip].forEach(e=>{e.oninput=secRefresh;e.onchange=secRefresh;});
elClipOn.addEventListener('change',updMode);
elMode.onchange=updMode;

// ---- nucleus on/off toggle: nucleus ON → translucent membrane (see the nucleus inside);
//      nucleus OFF → OPAQUE cells so the cell-cell CONTACT regions read clearly ----
const elNucOn=document.getElementById('nucon');
function applyNuc(){
  if(nucMesh) nucMesh.visible = (nucOn && elNucOn.checked);
  if(nucOn){
    mat.transparent = elNucOn.checked;
    mat.opacity     = elNucOn.checked ? 0.42 : 1.0;
    mat.depthWrite  = !elNucOn.checked;
    mat.needsUpdate = true;
  }
}
if(!nucOn){ const nr=document.getElementById('nucrow'); if(nr) nr.style.display='none'; }
elNucOn.onchange=applyNuc;
applyNuc();

// ---- colour-overlay selector (cell id / junction contact / deformation load) ----
const elCMode=document.getElementById('cmode'), elCHint=document.getElementById('cmodehint');
const _sp=(P.hasstress?('virial von-Mises σ_vm, peak '+(P.stress_pk||0).toPrecision(2)+' Pa')
                       :'⚠ residual-force fallback (this npz has no virial stress) — re-run for σ_vm');
const CHINTS={cell:'each cell a distinct colour',
  stress:'blue→red = REAL '+_sp+' — concentrates at load-bearing junctions',
  pressure:'blue→red = |hydrostatic pressure| per node (turgor-dominated), peak '+(P.press_pk||0).toPrecision(2)+' Pa',
  residual:'blue→red = |net RESIDUAL force| per node (∝ velocity, →0 at convergence — NOT a stress)',
  junction:'blue→red = number of active cadherin bonds per node'};
elCMode.onchange=()=>{ elCHint.textContent=CHINTS[elCMode.value]||''; applyColors(cur); };
elCHint.textContent=CHINTS[elCMode.value]||'';
const elJunOn=document.getElementById('junon');
if(!P.hasmech){ const jr=document.getElementById('junrow'); if(jr) jr.style.display='none';
  ['stress','pressure','residual','junction'].forEach(v=>{
    const o=elCMode.querySelector('option[value='+v+']'); if(o) o.remove(); }); }
if(elJunOn) elJunOn.onchange=()=>applyJunctions(cur);

elThickv.textContent=(elThick.value/10).toFixed(0)+'%';
elFpsv.textContent=fps;

// DEFAULT: section OFF — show the WHOLE spheroid surface first, coloured by the FEM stress
// field with cadherin bonds overlaid (the clearest "deformation/stress + junctions" view the
// user asked for). Peel is preselected so turning section ON gives a clean SOLID cross-section
// (peel never slices a shell → no hollow interior); cut/slab available from the dropdown.
elClipOn.checked=false;
elMode.value='peel';
// deep-link overrides so a specific cross-section is shareable / headless-renderable:
//   ?section=cut|peel|slab|off  &axis=0|1|2  &pos=0..1  &flip=0|1  &thick=0..1  &frame=N|last
let _f0=F-1;   // open on the SETTLED last frame (full stress + all cadherin bonds); scrub back to watch assembly
try{
  const q=new URLSearchParams(location.search);
  if(q.has('section')){const s=q.get('section');
    if(s==='off'){elClipOn.checked=false;}
    else{elClipOn.checked=true; if(['peel','slab','cut'].includes(s)) elMode.value=s;}}
  if(q.has('axis')){const a=+q.get('axis'); if(a>=0&&a<=2) elAxis.value=String(a);}
  if(q.has('pos')){elClip.value=String(Math.round(Math.max(0,Math.min(1,+q.get('pos')))*1000));}
  if(q.has('flip')){elFlip.checked=(q.get('flip')==='1'||q.get('flip')==='true');}
  if(q.has('thick')){elThick.value=String(Math.round(Math.max(0,Math.min(1,+q.get('thick')))*1000));
    elThickv.textContent=(elThick.value/10).toFixed(0)+'%';}
  const fr=q.get('frame');
  if(fr==='last') _f0=F-1; else if(fr!==null && !isNaN(+fr)) _f0=+fr;
}catch(e){}
setFrame(_f0);
updMode();

// event-driven rendering: NO rAF when idle. Was an unconditional 60fps renderer.render →
// idle tabs pegged CPU/GPU forever, catastrophic under software-WebGL (2026-07-08).
// The loop reschedules itself only while playback runs or a repaint is pending, then
// stops — a still viewer costs zero CPU and headless capture can settle to idle.
let dirty=true, running=false;
function wake(){ if(!running){ running=true; requestAnimationFrame(loop); } }
function requestRender(){ dirty=true; wake(); }
ctr.addEventListener('change', requestRender);
for(const ev of ['input','change','click']) document.addEventListener(ev, requestRender);
addEventListener('resize',()=>{cam.aspect=innerWidth/innerHeight;cam.updateProjectionMatrix();
  renderer.setSize(innerWidth,innerHeight); requestRender();});
function loop(t){ ctr.update();
  if(playing && t-last>1000/fps){last=t; setFrame(cur+1); dirty=true;}
  if(dirty){ dirty=false; renderer.render(scene,cam); }
  if(playing || dirty){ requestAnimationFrame(loop); } else { running=false; }}
requestRender();
</script></body></html>"""


if __name__ == "__main__":
    main()
