"""Interactive 'grab & stretch' cell playground — self-contained HTML (Three.js + in-browser PBD).

PI 2026-07-07: a TOY, not a physics unit. The goal is to *play* with the FF fiber-network cell in a
browser — grab it with the mouse, stretch/squish it, watch it deform and spring back — purely to feel
how the cell moves, NOT to be fine-grained (that carve-out is explicit; this runs a fast real-time
approximate solver, never the Warp kernels). It therefore lives outside the mechanistic hard-rule
regime that governs the runtime engine.

Compartments (PI 2026-07-07 "멤브레인이랑 핵 등 컴파트먼트 보이게만"): a translucent MEMBRANE surface
and a blue NUCLEUS are drawn for a cell-like look. They are VISUAL ONLY — not extra physics bodies:
each is an icosphere whose vertices are SKINNED (weighted to their nearest cortex nodes at rest), so it
follows/deforms with the grabbed cortex but adds no forces, no volume of its own, no collision. The
membrane rides just outside the fiber shell; the nucleus is a ~0.55 R scaled-down copy of the same
deformed shell. Toggle either off in the UI.

Density (PI 2026-07-07 "필라멘트가 너무 적은거 아닌가"): the page embeds several DENSITY TIERS
(filament count) selectable at runtime. The real-time grab solver has a hard ceiling — a V8 benchmark
of this exact hot-loop puts ~2–3k filaments (~15–23k nodes) at 5–9 ms/frame (60 fps-comfortable) and
~5k filaments (~40k nodes) at ~15 ms. Native cortex is ≈38k–70k filaments — out of real-time reach.

What is real vs. approximated
-----------------------------
REAL (reused verbatim from the engine): the cell GEOMETRY + topology is a genuine FF cortex from
:func:`ff.cortex_assembly.build_cortex_network` (great-circle actin arcs on the MCF7 R=7.5 µm sphere,
KU-3.17); segments carry true rest lengths; crosslinks are the crosslinker layer that makes the loose
fiber set one cohesive elastic shell.

APPROXIMATED (real-time interactivity, NOT literature-anchored — do not read physics off this):
Position-Based Dynamics (Verlet + distance constraints, overdamped); a soft GLOBAL volume-proxy
turgor (Σ rᵢ³ → V₀); stiffness/turgor sliders in arbitrary units; membrane/nucleus are kinematic skins.
Absent entirely: bending rigidity, thermal noise, excluded volume, real η, active myosin, adhesion.

Output: one double-clickable offline HTML (three.js r128 from CDN; geometry base64-embedded).
"""
from __future__ import annotations

import argparse
import base64
import json

import numpy as np

from ffn_sim.ff.cortex_assembly import CortexParams, build_cortex_network

SKIN_K = 6  # cortex nodes each compartment vertex is skinned to


def _b64(arr: np.ndarray, dtype) -> str:
    """Base64-encode a numpy array as a contiguous little-endian buffer of ``dtype``."""
    return base64.b64encode(np.ascontiguousarray(arr, dtype=dtype).tobytes()).decode("ascii")


def _icosphere(subdiv: int) -> tuple[np.ndarray, np.ndarray]:
    """Unit-sphere icosphere (subdivided icosahedron). Returns (verts (V,3), faces (F,3) uint32)."""
    t = (1.0 + 5.0 ** 0.5) / 2.0
    verts = [(-1, t, 0), (1, t, 0), (-1, -t, 0), (1, -t, 0), (0, -1, t), (0, 1, t),
             (0, -1, -t), (0, 1, -t), (t, 0, -1), (t, 0, 1), (-t, 0, -1), (-t, 0, 1)]
    verts = [tuple(np.asarray(v) / np.linalg.norm(v)) for v in verts]
    faces = [(0, 11, 5), (0, 5, 1), (0, 1, 7), (0, 7, 10), (0, 10, 11), (1, 5, 9), (5, 11, 4),
             (11, 10, 2), (10, 7, 6), (7, 1, 8), (3, 9, 4), (3, 4, 2), (3, 2, 6), (3, 6, 8),
             (3, 8, 9), (4, 9, 5), (2, 4, 11), (6, 2, 10), (8, 6, 7), (9, 8, 1)]
    for _ in range(subdiv):
        cache: dict[tuple[int, int], int] = {}
        new_faces = []

        def mid(a: int, b: int) -> int:
            key = (a, b) if a < b else (b, a)
            if key in cache:
                return cache[key]
            m = np.asarray(verts[a]) + np.asarray(verts[b])
            m = m / np.linalg.norm(m)
            verts.append(tuple(m))
            cache[key] = len(verts) - 1
            return cache[key]

        for a, b, c in faces:
            ab, bc, ca = mid(a, b), mid(b, c), mid(c, a)
            new_faces += [(a, ab, ca), (b, bc, ab), (c, ca, bc), (ab, bc, ca)]
        faces = new_faces
    return np.asarray(verts, dtype=np.float32), np.asarray(faces, dtype=np.uint32)


def _skin(dirs: np.ndarray, cortex_pos: np.ndarray, R: float, k: int) -> tuple[np.ndarray, np.ndarray]:
    """Skin each compartment vertex (unit direction) to its ``k`` nearest cortex rest nodes.

    Query point is ``R·dir`` (the rest shell in that direction); weights are inverse-distance,
    row-normalised. Returns ``(idx (V,k) uint32, w (V,k) float32)`` — the browser rebuilds each
    vertex every frame as the weighted average of those (moving) cortex nodes, so the compartment
    surface follows the deformed shell without being a physics body.
    """
    pts = dirs * R
    try:
        from scipy.spatial import cKDTree  # type: ignore

        d, idx = cKDTree(cortex_pos).query(pts, k=k)
    except Exception:  # pragma: no cover - brute-force fallback
        idx = np.zeros((len(pts), k), np.int64)
        d = np.zeros((len(pts), k), np.float64)
        for i, p in enumerate(pts):
            dd = np.linalg.norm(cortex_pos - p, axis=1)
            nn = np.argpartition(dd, k)[:k]
            nn = nn[np.argsort(dd[nn])]
            idx[i] = nn
            d[i] = dd[nn]
    d = np.atleast_2d(d)
    idx = np.atleast_2d(idx)
    w = 1.0 / (d + 1e-6)
    w = w / w.sum(axis=1, keepdims=True)
    return idx.astype(np.uint32), w.astype(np.float32)


def build_tier(*, n_filaments: int, cutoff: float, k_per_node: int, R_um: float, beads: int,
               seg_um: float, length_dist: str, seed: int,
               mem_dirs: np.ndarray, nuc_dirs: np.ndarray) -> dict:
    """Build one coarse cortex + crosslinks tier + membrane/nucleus skins; pack browser arrays."""
    params = CortexParams(R_um=R_um, beads_per_filament=beads, seg_um=seg_um)
    rng = np.random.default_rng(seed)
    net, meta = build_cortex_network(params, rng=rng, n_filaments=n_filaments,
                                     orientation="isotropic", length_dist=length_dist)
    pos = net.pos.astype(np.float32)
    fid = np.zeros(net.n_nodes, dtype=np.int64)
    for f in range(len(net.fiber_offsets) - 1):
        fid[net.fiber_offsets[f]:net.fiber_offsets[f + 1]] = f

    seg = net.segments.astype(np.int64)
    seg_rest = net.seg_rest.astype(np.float32)
    max_links = int(1.8 * net.n_nodes)
    xl = _crosslinks(net.pos, fid, cutoff=cutoff, k_per_node=k_per_node, max_links=max_links, rng=rng)
    xl_rest = (np.linalg.norm(net.pos[xl[:, 1]] - net.pos[xl[:, 0]], axis=1).astype(np.float32)
               if len(xl) else np.zeros(0, np.float32))

    si = np.concatenate([seg[:, 0], xl[:, 0]]).astype(np.uint32) if len(xl) else seg[:, 0].astype(np.uint32)
    sj = np.concatenate([seg[:, 1], xl[:, 1]]).astype(np.uint32) if len(xl) else seg[:, 1].astype(np.uint32)
    srest = np.concatenate([seg_rest, xl_rest]).astype(np.float32)
    skind = np.concatenate([np.zeros(len(seg), np.uint8), np.ones(len(xl), np.uint8)])

    mem_idx, mem_w = _skin(mem_dirs, net.pos, R_um, SKIN_K)
    nuc_idx, nuc_w = _skin(nuc_dirs, net.pos, R_um, SKIN_K)

    return {
        "name": f"{n_filaments:,} filaments",
        "n_filaments": int(n_filaments), "n_nodes": int(pos.shape[0]),
        "n_segments": int(len(seg)), "n_crosslinks": int(len(xl)), "n_springs": int(len(srest)),
        "pos_b64": _b64(pos, np.float32),
        "si_b64": _b64(si, np.uint32), "sj_b64": _b64(sj, np.uint32),
        "srest_b64": _b64(srest, np.float32), "skind_b64": _b64(skind, np.uint8),
        "mem_idx_b64": _b64(mem_idx, np.uint32), "mem_w_b64": _b64(mem_w, np.float32),
        "nuc_idx_b64": _b64(nuc_idx, np.uint32), "nuc_w_b64": _b64(nuc_w, np.float32),
    }


def _crosslinks(pos: np.ndarray, fiber_id: np.ndarray, *, cutoff: float, k_per_node: int,
                max_links: int, rng: np.random.Generator) -> np.ndarray:
    """Proximity crosslinks between nodes on DIFFERENT fibers (greedy shortest-first, degree-capped)."""
    n = pos.shape[0]
    pairs: list[tuple[float, int, int]] = []
    try:
        from scipy.spatial import cKDTree  # type: ignore

        tree = cKDTree(pos)
        for i, j in tree.query_pairs(r=cutoff):
            if fiber_id[i] != fiber_id[j]:
                pairs.append((float(np.linalg.norm(pos[i] - pos[j])), int(i), int(j)))
    except Exception:  # pragma: no cover
        grid: dict[tuple[int, int, int], list[int]] = {}
        keys = np.floor(pos / cutoff).astype(int)
        for idx, k in enumerate(map(tuple, keys)):
            grid.setdefault(k, []).append(idx)
        seen: set[tuple[int, int]] = set()
        for idx in range(n):
            kx, ky, kz = keys[idx]
            for dx in (-1, 0, 1):
                for dy in (-1, 0, 1):
                    for dz in (-1, 0, 1):
                        for j in grid.get((kx + dx, ky + dy, kz + dz), ()):
                            if j == idx or fiber_id[idx] == fiber_id[j]:
                                continue
                            a, b = (idx, j) if idx < j else (j, idx)
                            if (a, b) in seen:
                                continue
                            d = float(np.linalg.norm(pos[a] - pos[b]))
                            if d <= cutoff:
                                seen.add((a, b))
                                pairs.append((d, a, b))
    pairs.sort(key=lambda t: t[0])
    deg = np.zeros(n, dtype=np.int32)
    out: list[tuple[int, int]] = []
    seen2: set[tuple[int, int]] = set()
    for d, i, j in pairs:
        a, b = (i, j) if i < j else (j, i)
        if (a, b) in seen2 or deg[a] >= k_per_node or deg[b] >= k_per_node:
            continue
        seen2.add((a, b))
        out.append((a, b))
        deg[a] += 1
        deg[b] += 1
        if len(out) >= max_links:
            break
    return np.asarray(out, dtype=np.int64) if out else np.zeros((0, 2), np.int64)


def build_payload(*, tiers: list[int], default_tier: int, cutoff: float, k_per_node: int,
                  R_um: float, beads: int, seg_um: float, length_dist: str, seed: int) -> dict:
    """Build every density tier + the shared membrane/nucleus icosphere topology."""
    mem_v, mem_f = _icosphere(4)     # 2562 verts — smooth membrane
    nuc_v, nuc_f = _icosphere(2)     # 162 verts — nucleus blob
    tier_dicts = [build_tier(n_filaments=n, cutoff=cutoff, k_per_node=k_per_node, R_um=R_um,
                             beads=beads, seg_um=seg_um, length_dist=length_dist, seed=seed,
                             mem_dirs=mem_v, nuc_dirs=nuc_v)
                  for n in tiers]
    return {
        "R_um": float(R_um), "skin_k": SKIN_K,
        "default_tier": int(np.clip(default_tier, 0, len(tier_dicts) - 1)),
        "membrane": {"verts_b64": _b64(mem_v, np.float32), "faces_b64": _b64(mem_f, np.uint32),
                     "nv": int(len(mem_v)), "inflate": 1.05},
        "nucleus": {"verts_b64": _b64(nuc_v, np.float32), "faces_b64": _b64(nuc_f, np.uint32),
                    "nv": int(len(nuc_v)), "scale": 0.55},
        "tiers": tier_dicts,
    }


def build_html(payload: dict, out: str, title: str) -> str:
    """Write the self-contained playground HTML; returns the out path."""
    html = _HTML.replace("/*__PAYLOAD__*/", json.dumps({**payload, "title": title}))
    with open(out, "w") as f:
        f.write(html)
    return out


_HTML = r"""<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1"><title>FF cell playground</title>
<style>
 :root{--bg:#0b0e14;--panel:#141a24f2;--line:#232c3b;--fg:#e6edf3;--mut:#9aa7b8;--acc:#5ac8fa;--grn:#76e0a0}
 body{margin:0;overflow:hidden;background:var(--bg);font:14px/1.5 -apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,sans-serif;color:var(--fg)}
 #ui{position:absolute;top:12px;left:12px;background:var(--panel);padding:12px 14px;border-radius:10px;
     border:1px solid var(--line);width:250px}
 #ui h1{font-size:15px;margin:0 0 2px} #ui .sub{color:var(--mut);font-size:11.5px;margin-bottom:10px}
 #ui .row{margin:9px 0;display:flex;align-items:center;gap:8px}
 #ui label{min-width:74px;color:var(--mut);font-size:12px}
 input[type=range]{flex:1;accent-color:var(--acc)}
 select{flex:1;background:#0d1117;color:var(--fg);border:1px solid var(--line);border-radius:5px;padding:3px 6px;font-size:12px}
 .val{min-width:30px;text-align:right;color:var(--fg);font-variant-numeric:tabular-nums;font-size:12px}
 button{background:#1c2635;color:var(--fg);border:1px solid var(--line);border-radius:6px;padding:5px 9px;cursor:pointer;font-size:12px}
 button:hover{border-color:var(--acc)} button.on{background:var(--acc);color:#04121e;border-color:var(--acc);font-weight:600}
 .toggles{display:flex;flex-wrap:wrap;gap:6px;margin-top:4px}
 #hint{position:absolute;bottom:12px;left:12px;background:var(--panel);border:1px solid var(--line);
       border-radius:10px;padding:9px 12px;color:var(--mut);font-size:12px;max-width:300px}
 #hint b{color:var(--fg)} .k{color:var(--grn)}
 #leg{position:absolute;bottom:12px;right:12px;background:var(--panel);border:1px solid var(--line);
      border-radius:10px;padding:9px 12px;font-size:11.5px}
 #leg .bar{height:9px;width:150px;border-radius:5px;margin:5px 0 3px;
   background:linear-gradient(90deg,#4a86ff,#8fb6ff,#ffcc66,#ff8a3d,#ff3b3b)}
 #leg .lab{display:flex;justify-content:space-between;color:var(--mut)}
 #stat{position:absolute;top:12px;right:12px;background:var(--panel);border:1px solid var(--line);
       border-radius:10px;padding:8px 11px;font-size:11.5px;color:var(--mut);text-align:right;font-variant-numeric:tabular-nums}
 #stat b{color:var(--fg)}
</style></head><body>
<div id="ui">
  <h1>FF cell — grab &amp; stretch</h1>
  <div class="sub">real fiber cortex · membrane + nucleus (visual) · toy solver</div>
  <div class="row"><label>density</label><select id="tier"></select></div>
  <div class="row"><label>turgor</label><input id="turgor" type="range" min="0" max="1.5" step="0.02" value="0.7"><span class="val" id="turgorv"></span></div>
  <div class="row"><label>stiffness</label><input id="stiff" type="range" min="0.1" max="1" step="0.02" value="0.55"><span class="val" id="stiffv"></span></div>
  <div class="row"><label>viscosity</label><input id="damp" type="range" min="0.5" max="0.98" step="0.01" value="0.9"><span class="val" id="dampv"></span></div>
  <div class="row"><label>grab hold</label><input id="hold" type="range" min="1" max="12" step="1" value="4"><span class="val" id="holdv"></span></div>
  <div class="toggles">
    <button id="bMem" class="on">membrane</button>
    <button id="bNuc" class="on">nucleus</button>
    <button id="bStrain" class="on">strain</button>
    <button id="bXlink">crosslinks</button>
    <button id="bNodes">nodes</button>
    <button id="bReset">reset</button>
  </div>
</div>
<div id="stat"></div>
<div id="leg"><div>fiber strain</div><div class="bar"></div><div class="lab"><span>−squish</span><span>rest</span><span>stretch+</span></div></div>
<div id="hint"><b class="k">grab</b> a point &amp; drag to stretch · <b>drag empty space</b> = rotate · <b>scroll</b> = zoom · release to spring back</div>

<script src="https://cdnjs.cloudflare.com/ajax/libs/three.js/r128/three.min.js"></script>
<script src="https://cdn.jsdelivr.net/npm/three@0.128.0/examples/js/controls/OrbitControls.js"></script>
<script>
const P = /*__PAYLOAD__*/;
document.title = P.title || "FF cell playground";
function decF32(b){const s=atob(b);const u=new Uint8Array(s.length);for(let i=0;i<s.length;i++)u[i]=s.charCodeAt(i);return new Float32Array(u.buffer);}
function decU32(b){const s=atob(b);const u=new Uint8Array(s.length);for(let i=0;i<s.length;i++)u[i]=s.charCodeAt(i);return new Uint32Array(u.buffer);}
function decU8(b){const s=atob(b);const u=new Uint8Array(s.length);for(let i=0;i<s.length;i++)u[i]=s.charCodeAt(i);return u;}
const R=P.R_um||7.5, K=P.skin_k;

// shared compartment topology (unit-sphere directions + faces), decoded once
const MEM_DIR=decF32(P.membrane.verts_b64), MEM_FACE=decU32(P.membrane.faces_b64), MEM_NV=P.membrane.nv, MEM_INFLATE=P.membrane.inflate;
const NUC_DIR=decF32(P.nucleus.verts_b64), NUC_FACE=decU32(P.nucleus.faces_b64), NUC_NV=P.nucleus.nv, NUC_SCALE=P.nucleus.scale;

// ---- three.js scene -------------------------------------------------------------------------
const isMobile=/Mobi|Android|iPhone|iPad|iPod/i.test(navigator.userAgent)||(navigator.maxTouchPoints>1&&Math.min(innerWidth,innerHeight)<820);
const renderer=new THREE.WebGLRenderer({antialias:!isMobile,powerPreference:"high-performance"});
renderer.setSize(innerWidth,innerHeight);renderer.setPixelRatio(Math.min(devicePixelRatio,isMobile?1.5:2));
document.body.appendChild(renderer.domElement);
const scene=new THREE.Scene();scene.background=new THREE.Color(0x0b0e14);
const cam=new THREE.PerspectiveCamera(45,innerWidth/innerHeight,0.05,5000);
cam.position.set(R*2.4,R*1.5,R*2.4);
const ctr=new THREE.OrbitControls(cam,renderer.domElement);ctr.target.set(0,0,0);
ctr.enableDamping=true;ctr.dampingFactor=0.12;ctr.update();
scene.add(new THREE.AmbientLight(0xffffff,0.75));
const dl=new THREE.DirectionalLight(0xffffff,0.55);dl.position.set(1,1.5,2);scene.add(dl);

// ---- per-tier state (rebuilt on density switch) ---------------------------------------------
let N,M,pos,prev,r0,V0,si,sj,srest,skind,winv,segIdx,xlIdx,nSeg,nXl;
let memIdx,memW,nucIdx,nucW;                     // compartment skin bindings (per tier)
let fibers,xlink,nodes,fgeo,fpos,fcol,xgeo,xpos,ngeo,npos,memMesh,memGeo,memPos,nucMesh,nucGeo,nucPos;
let cellCx=0,cellCy=0,cellCz=0;                  // cell centroid (for compartment reconstruction)
const gmat=new THREE.MeshBasicMaterial({color:0xffe14d});
const grabDot=new THREE.Mesh(new THREE.SphereGeometry(0.42,16,16),gmat);grabDot.visible=false;scene.add(grabDot);

function disposeTier(){
  for(const o of [fibers,xlink,nodes,memMesh,nucMesh]) if(o){scene.remove(o);o.geometry.dispose();o.material.dispose();}
}
function loadTier(idx){
  disposeTier();grabbed=-1;grabDot.visible=false;
  const T=P.tiers[idx];
  N=T.n_nodes;const pos0=decF32(T.pos_b64);pos=pos0.slice();prev=pos0.slice();
  si=decU32(T.si_b64);sj=decU32(T.sj_b64);srest=decF32(T.srest_b64);skind=decU8(T.skind_b64);M=srest.length;
  winv=new Float32Array(N).fill(1);
  memIdx=decU32(T.mem_idx_b64);memW=decF32(T.mem_w_b64);nucIdx=decU32(T.nuc_idx_b64);nucW=decF32(T.nuc_w_b64);
  r0=new Float32Array(N);V0=0;
  for(let i=0;i<N;i++){const x=pos0[3*i],y=pos0[3*i+1],z=pos0[3*i+2];const r=Math.hypot(x,y,z);r0[i]=r;V0+=r*r*r;}
  segIdx=[];xlIdx=[];for(let s=0;s<M;s++)(skind[s]===0?segIdx:xlIdx).push(s);nSeg=segIdx.length;nXl=xlIdx.length;

  fpos=new Float32Array(nSeg*6);fcol=new Float32Array(nSeg*6);
  fgeo=new THREE.BufferGeometry();
  fgeo.setAttribute('position',new THREE.BufferAttribute(fpos,3));fgeo.setAttribute('color',new THREE.BufferAttribute(fcol,3));
  fibers=new THREE.LineSegments(fgeo,new THREE.LineBasicMaterial({vertexColors:true,transparent:true,opacity:0.85}));
  fibers.frustumCulled=false;scene.add(fibers);

  xpos=new Float32Array(nXl*6);
  xgeo=new THREE.BufferGeometry();xgeo.setAttribute('position',new THREE.BufferAttribute(xpos,3));
  xlink=new THREE.LineSegments(xgeo,new THREE.LineBasicMaterial({color:0x39506e,transparent:true,opacity:0.3}));
  xlink.frustumCulled=false;xlink.visible=U.showXl;scene.add(xlink);

  npos=new Float32Array(N*3);
  ngeo=new THREE.BufferGeometry();ngeo.setAttribute('position',new THREE.BufferAttribute(npos,3));
  const psize=N>12000?0.1:(N>6000?0.13:0.16);
  nodes=new THREE.Points(ngeo,new THREE.PointsMaterial({color:0xffd98a,size:psize,sizeAttenuation:true,transparent:true,opacity:0.9}));
  nodes.frustumCulled=false;nodes.visible=U.showNodes;scene.add(nodes);

  // membrane: translucent skinned surface, just outside the fiber shell
  memPos=new Float32Array(MEM_NV*3);memGeo=new THREE.BufferGeometry();
  memGeo.setAttribute('position',new THREE.BufferAttribute(memPos,3));memGeo.setIndex(new THREE.BufferAttribute(MEM_FACE,1));
  memMesh=new THREE.Mesh(memGeo,new THREE.MeshStandardMaterial({color:0xbcd4ff,transparent:true,opacity:0.15,
    side:THREE.DoubleSide,depthWrite:false,roughness:0.35,metalness:0.0}));
  memMesh.frustumCulled=false;memMesh.renderOrder=3;memMesh.visible=U.showMem;scene.add(memMesh);
  // nucleus: blue skinned blob at ~0.55 R
  nucPos=new Float32Array(NUC_NV*3);nucGeo=new THREE.BufferGeometry();
  nucGeo.setAttribute('position',new THREE.BufferAttribute(nucPos,3));nucGeo.setIndex(new THREE.BufferAttribute(NUC_FACE,1));
  nucMesh=new THREE.Mesh(nucGeo,new THREE.MeshStandardMaterial({color:0x6f9bef,transparent:true,opacity:0.6,
    side:THREE.DoubleSide,depthWrite:false,roughness:0.5,metalness:0.0}));
  nucMesh.frustumCulled=false;nucMesh.renderOrder=2;nucMesh.visible=U.showNuc;scene.add(nucMesh);

  iters=5;frameN=0;updateBuffers();      // prime geometry so the first frame is correct
}

// ---- UI -------------------------------------------------------------------------------------
const U={turgor:0.7,stiff:0.55,damp:0.9,hold:4,strain:true,showXl:false,showNodes:false,showMem:true,showNuc:true};
function bind(id,key,fmt){const el=document.getElementById(id),v=document.getElementById(id+'v');
  const set=()=>{U[key]=parseFloat(el.value);if(v)v.textContent=fmt?fmt(U[key]):el.value;};el.oninput=set;set();}
bind('turgor','turgor');bind('stiff','stiff');bind('damp','damp');bind('hold','hold',x=>x.toFixed(0));
function tog(id,key,obj){const b=document.getElementById(id);b.classList.toggle('on',!!U[key]);
  b.onclick=()=>{U[key]=!U[key];b.classList.toggle('on',U[key]);if(obj)obj();};}
tog('bMem','showMem',()=>memMesh.visible=U.showMem);tog('bNuc','showNuc',()=>nucMesh.visible=U.showNuc);
tog('bStrain','strain');tog('bXlink','showXl',()=>xlink.visible=U.showXl);tog('bNodes','showNodes',()=>nodes.visible=U.showNodes);
document.getElementById('bReset').onclick=()=>{pos.set(decF32(P.tiers[curTier].pos_b64));prev.set(pos);grabbed=-1;grabDot.visible=false;winv.fill(1);};

const tsel=document.getElementById('tier');
P.tiers.forEach((t,i)=>{const o=document.createElement('option');o.value=i;o.textContent=`${t.name} · ${t.n_nodes.toLocaleString()} nodes`;tsel.appendChild(o);});
let curTier=P.default_tier;tsel.value=curTier;tsel.onchange=()=>{curTier=parseInt(tsel.value);loadTier(curTier);};

// ---- picking + grab -------------------------------------------------------------------------
const ray=new THREE.Raycaster();const ndc=new THREE.Vector2();
let grabbed=-1;const grabTarget=new THREE.Vector3();const grabPlane=new THREE.Plane();const hitP=new THREE.Vector3();
let neighbours=[];
function pointerNDC(e){ndc.x=(e.clientX/innerWidth)*2-1;ndc.y=-(e.clientY/innerHeight)*2+1;}
function pickNode(){
  ray.setFromCamera(ndc,cam);const o=ray.ray.origin,d=ray.ray.direction;let best=-1,bestT=1e9;const thr=0.9;
  for(let i=0;i<N;i++){const x=pos[3*i]-o.x,y=pos[3*i+1]-o.y,z=pos[3*i+2]-o.z;const t=x*d.x+y*d.y+z*d.z;if(t<0)continue;
    const px=x-t*d.x,py=y-t*d.y,pz=z-t*d.z;if(Math.hypot(px,py,pz)<thr && t<bestT){bestT=t;best=i;}}
  return best;
}
function ringOf(idx){const s=new Set();for(let e=0;e<M;e++){if(si[e]===idx)s.add(sj[e]);else if(sj[e]===idx)s.add(si[e]);}return [...s];}
addEventListener('pointerdown',e=>{                       // capture phase → pre-empts OrbitControls rotate
  pointerNDC(e);const hit=pickNode();
  if(hit>=0){grabbed=hit;neighbours=ringOf(hit);ctr.enabled=false;
    grabPlane.setFromNormalAndCoplanarPoint(cam.getWorldDirection(new THREE.Vector3()),new THREE.Vector3(pos[3*hit],pos[3*hit+1],pos[3*hit+2]));
    grabTarget.set(pos[3*hit],pos[3*hit+1],pos[3*hit+2]);grabDot.visible=true;}
},true);
addEventListener('pointermove',e=>{if(grabbed<0)return;pointerNDC(e);ray.setFromCamera(ndc,cam);if(ray.ray.intersectPlane(grabPlane,hitP))grabTarget.copy(hitP);});
function release(){if(grabbed<0)return;grabbed=-1;neighbours=[];ctr.enabled=true;grabDot.visible=false;winv.fill(1);}
addEventListener('pointerup',release);

// ---- PBD solver (turgor once/frame; iters adapts to frame cost) ------------------------------
let iters=5,frameN=0;const compEvery=isMobile?2:1;   // compartments (visual) update at half-rate on mobile
function step(){
  const damp=U.damp;
  for(let i=0;i<N;i++){if(i===grabbed)continue;for(let c=0;c<3;c++){const k=3*i+c;const cur=pos[k];const v=(cur-prev[k])*damp;prev[k]=cur;pos[k]=cur+v;}}
  if(grabbed>=0){winv[grabbed]=0;pos[3*grabbed]=grabTarget.x;pos[3*grabbed+1]=grabTarget.y;pos[3*grabbed+2]=grabTarget.z;
    prev[3*grabbed]=grabTarget.x;prev[3*grabbed+1]=grabTarget.y;prev[3*grabbed+2]=grabTarget.z;}
  const kFib=U.stiff,kXl=U.stiff*0.55,pull=Math.min(1,U.hold/12);
  for(let it=0;it<iters;it++){
    for(let s=0;s<M;s++){const a=si[s],b=sj[s];const wa=winv[a],wb=winv[b];const wsum=wa+wb;if(wsum===0)continue;
      const ax=pos[3*a],ay=pos[3*a+1],az=pos[3*a+2],bx=pos[3*b],by=pos[3*b+1],bz=pos[3*b+2];
      let dx=bx-ax,dy=by-ay,dz=bz-az;let d=Math.hypot(dx,dy,dz);if(d<1e-6)continue;
      const diff=(d-srest[s])/d*(skind[s]===0?kFib:kXl);const ca=diff*wa/wsum,cb=diff*wb/wsum;
      pos[3*a]+=dx*ca;pos[3*a+1]+=dy*ca;pos[3*a+2]+=dz*ca;pos[3*b]-=dx*cb;pos[3*b+1]-=dy*cb;pos[3*b+2]-=dz*cb;}
    if(grabbed>=0 && pull>0)for(const j of neighbours){if(winv[j]===0)continue;
      pos[3*j]+=(grabTarget.x-pos[3*j])*pull*0.15;pos[3*j+1]+=(grabTarget.y-pos[3*j+1])*pull*0.15;pos[3*j+2]+=(grabTarget.z-pos[3*j+2])*pull*0.15;}
  }
  if(U.turgor>0){
    let V=0;for(let i=0;i<N;i++){const x=pos[3*i],y=pos[3*i+1],z=pos[3*i+2];V+=Math.hypot(x,y,z)**3;}
    const g=(V0-V)/Math.max(V0,1e-6)*U.turgor*0.5,s2=1+g;
    for(let i=0;i<N;i++){if(winv[i]===0)continue;pos[3*i]*=s2;pos[3*i+1]*=s2;pos[3*i+2]*=s2;
      const r2=Math.hypot(pos[3*i],pos[3*i+1],pos[3*i+2])||1e-6;const t=0.03*(r0[i]-r2)/r2;pos[3*i]+=pos[3*i]*t;pos[3*i+1]+=pos[3*i+1]*t;pos[3*i+2]+=pos[3*i+2]*t;}
  }
  let cx=0,cy=0,cz=0;for(let i=0;i<N;i++){cx+=pos[3*i];cy+=pos[3*i+1];cz+=pos[3*i+2];}cx/=N;cy/=N;cz/=N;const rc=0.01;
  for(let i=0;i<N;i++){if(i===grabbed)continue;pos[3*i]-=cx*rc;pos[3*i+1]-=cy*rc;pos[3*i+2]-=cz*rc;}
  cellCx=cx*(1-rc);cellCy=cy*(1-rc);cellCz=cz*(1-rc);
}

// ---- render ---------------------------------------------------------------------------------
function strainColor(strain,out,off){const t=Math.max(-1,Math.min(1,strain/0.3));let r,g,b;
  if(!U.strain){r=1.0;g=0.83;b=0.42;}else if(t>=0){r=1.0;g=0.83-0.55*t;b=0.42-0.30*t;}else{const u=-t;r=1.0-0.68*u;g=0.83-0.31*u;b=0.42+0.58*u;}
  out[off]=r;out[off+1]=g;out[off+2]=b;}
function skinInto(dst,nv,idx,w,scale){    // dst[v] = centroid + (Σ w·node - centroid)*scale
  for(let v=0;v<nv;v++){let Sx=0,Sy=0,Sz=0;const base=v*K;
    for(let k=0;k<K;k++){const n=idx[base+k],ww=w[base+k];Sx+=ww*pos[3*n];Sy+=ww*pos[3*n+1];Sz+=ww*pos[3*n+2];}
    dst[3*v]=cellCx+(Sx-cellCx)*scale;dst[3*v+1]=cellCy+(Sy-cellCy)*scale;dst[3*v+2]=cellCz+(Sz-cellCz)*scale;}
}
function updateBuffers(){
  for(let n=0;n<nSeg;n++){const s=segIdx[n];const a=si[s],b=sj[s];
    const ax=pos[3*a],ay=pos[3*a+1],az=pos[3*a+2],bx=pos[3*b],by=pos[3*b+1],bz=pos[3*b+2];const o=6*n;
    fpos[o]=ax;fpos[o+1]=ay;fpos[o+2]=az;fpos[o+3]=bx;fpos[o+4]=by;fpos[o+5]=bz;
    const d=Math.hypot(bx-ax,by-ay,bz-az);const st=d/Math.max(srest[s],1e-6)-1;strainColor(st,fcol,o);strainColor(st,fcol,o+3);}
  fgeo.attributes.position.needsUpdate=true;fgeo.attributes.color.needsUpdate=true;
  if(U.showXl){for(let n=0;n<nXl;n++){const s=xlIdx[n];const a=si[s],b=sj[s];const o=6*n;
    xpos[o]=pos[3*a];xpos[o+1]=pos[3*a+1];xpos[o+2]=pos[3*a+2];xpos[o+3]=pos[3*b];xpos[o+4]=pos[3*b+1];xpos[o+5]=pos[3*b+2];}xgeo.attributes.position.needsUpdate=true;}
  if(U.showNodes){npos.set(pos);ngeo.attributes.position.needsUpdate=true;}
  if((frameN%compEvery)===0){                    // compartments are visual; half-rate on mobile
    if(U.showMem){skinInto(memPos,MEM_NV,memIdx,memW,MEM_INFLATE);memGeo.attributes.position.needsUpdate=true;memGeo.computeVertexNormals();}
    if(U.showNuc){skinInto(nucPos,NUC_NV,nucIdx,nucW,NUC_SCALE);nucGeo.attributes.position.needsUpdate=true;nucGeo.computeVertexNormals();}
  }
  if(grabbed>=0)grabDot.position.set(pos[3*grabbed],pos[3*grabbed+1],pos[3*grabbed+2]);
}
addEventListener('resize',()=>{cam.aspect=innerWidth/innerHeight;cam.updateProjectionMatrix();renderer.setSize(innerWidth,innerHeight);});
let workEwma=8,dtEwma=16,dtMax=16,dprLow=false,tPrev=performance.now();
const statEl=document.getElementById('stat');
function tick(now){requestAnimationFrame(tick);frameN++;
  const t0=performance.now();step();updateBuffers();const work=performance.now()-t0;workEwma=0.92*workEwma+0.08*work;
  ctr.update();renderer.render(scene,cam);
  // adapt on the ACTUAL frame interval — captures GPU/compositor cost (overdraw, DOM), not just CPU solve
  const dt=now-tPrev;tPrev=now;dtEwma=0.9*dtEwma+0.1*Math.min(dt,100);dtMax=Math.max(dt,dtMax*0.94);
  if(dtEwma>22 && iters>3) iters--; else if(dtEwma<15 && iters<7) iters++;
  if(!dprLow && dtEwma>30 && iters<=3){dprLow=true;renderer.setPixelRatio(Math.max(1,renderer.getPixelRatio()*0.7));}
  if((frameN%20)===0)                            // throttle DOM writes: per-frame innerHTML is a real mobile cost
    statEl.innerHTML=`<b>${N.toLocaleString()}</b> nodes · <b>${nSeg.toLocaleString()}</b> fibers<br><b>${(1000/dtEwma).toFixed(0)}</b> fps · <b>${(1000/dtMax).toFixed(0)}</b> min · ${workEwma.toFixed(1)} ms · ${iters} it${dprLow?" · dpr↓":""}`;
}
loadTier(curTier);requestAnimationFrame(tick);
</script></body></html>"""


def main() -> None:
    ap = argparse.ArgumentParser(description="Export the interactive FF cell 'grab & stretch' playground HTML.")
    ap.add_argument("--out", default="ffn_sim/outputs/ff/figs/ff_playground.html")
    ap.add_argument("--tiers", default="1400,2100,3200,4500", help="comma-separated filament counts")
    ap.add_argument("--default-tier", type=int, default=1)
    ap.add_argument("--cutoff", type=float, default=0.85)
    ap.add_argument("--k-per-node", type=int, default=3)
    ap.add_argument("--R-um", type=float, default=7.5)
    ap.add_argument("--beads", type=int, default=7)
    ap.add_argument("--seg-um", type=float, default=0.5)
    ap.add_argument("--length-dist", default="exponential", choices=["mono", "exponential"])
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--title", default="FF cell — grab & stretch")
    args = ap.parse_args()

    tiers = [int(x) for x in args.tiers.split(",") if x.strip()]
    payload = build_payload(tiers=tiers, default_tier=args.default_tier, cutoff=args.cutoff,
                            k_per_node=args.k_per_node, R_um=args.R_um, beads=args.beads,
                            seg_um=args.seg_um, length_dist=args.length_dist, seed=args.seed)
    build_html(payload, args.out, args.title)
    print(f"[ff_playground] wrote {args.out}  ({len(tiers)} tiers, default={payload['default_tier']}, "
          f"membrane nv={payload['membrane']['nv']}, nucleus nv={payload['nucleus']['nv']})")
    for i, t in enumerate(payload["tiers"]):
        star = " *" if i == payload["default_tier"] else "  "
        print(f" {star} {t['name']:>16}: nodes={t['n_nodes']:>6}  fibers={t['n_segments']:>6}  "
              f"crosslinks={t['n_crosslinks']:>6}  springs={t['n_springs']:>6}")


if __name__ == "__main__":
    main()
