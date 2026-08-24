#!/usr/bin/env python
r"""Render the C-1 field run to a self-contained HTML viewer — the setup check, not decoration.

WHY THIS EXISTS.  An experiment wired to the wrong array produces a perfectly clean exponent from the
wrong field, and no scalar in the run record says so.  ``pore_mass_rel_drift`` would still be 1e-9, the
``1/e`` crossing would still interpolate, and the log-log fit would still return 2.00.  The only cheap
check that catches a mis-wired setup is looking at it: is the blob centred, is it round, is it inside
the box, is the contact patch the size the record claims, and is the thing that decays the thing the
observers are reading.

The viewer draws the mid-plane ``z`` slice at each stored time, with the blob width and the observer's
contact patch overlaid at true scale, beside the three observer curves and their ``1/e`` crossings.
Colour is a fixed per-frame-normalised map so a spreading blob does not merely look like a dimming one.

Charter: every run writes its own figures beside its record, verified in a REAL BROWSER, not by
grepping the HTML.  Host-only; reads a record, writes a file, evaluates no physics.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np

#: Central crop of the mid-plane slice, in multiples of the blob width. The far field is empty by
#: construction and shipping it would make the page tens of megabytes without showing anything; the
#: crop factor is recorded in the page so a reader knows the box is larger than what is drawn.
CROP_IN_A: float = 6.0

#: Decimals kept per cell. The field is normalised to its own frame maximum before rounding, so this is
#: display precision only — the record keeps full precision and the scoring never reads this file.
ROUND_DP: int = 4


def _crop(slice_2d: np.ndarray, dx: float, a_um: float) -> tuple[np.ndarray, float]:
    """Return the central crop and its half-width [µm]."""
    n = slice_2d.shape[0]
    half = min(int(math.ceil(CROP_IN_A * a_um / dx)), (n - 1) // 2)
    c = n // 2
    return slice_2d[c - half:c + half + 1, c - half:c + half + 1], half * dx


def build(record_path: Path, out_path: Path) -> None:
    rec = json.loads(record_path.read_text())
    pts = rec["sweeps"]["a_um"]
    build_commit = (rec.get("build") or {}).get("commit", "unknown")

    frames = []
    for p in pts:
        a, dx = p["a_um"], p["dx_um"]
        fr = []
        for s in p["slices"]:
            arr = np.asarray(s["z_midplane"], dtype=float)
            crop, half_um = _crop(arr, dx, a)
            peak = float(np.max(np.abs(crop))) or 1.0
            fr.append({
                "t": s["t_s"],
                "peak": peak,
                "half_um": half_um,
                "g": np.round(crop / peak, ROUND_DP).tolist(),
            })
        frames.append({
            "a": a, "dx": dx, "L": p["L_um"],
            "patch_um": p["patch_radius_um"],
            "tau": {k: v["tau_measured_s"] for k, v in p["observers"].items()},
            "tau_oracle": {k: v["tau_oracle_s"] for k, v in p["observers"].items()},
            "box_frac": p["box_contamination"],
            "mass_drift": p["pore_mass_rel_drift"],
            "curve": p["curve"],
            "slices": fr,
        })

    payload = json.dumps({"build": build_commit, "crop_in_a": CROP_IN_A, "runs": frames},
                         separators=(",", ":"))
    out_path.write_text(_HTML.replace("__PAYLOAD__", payload))


_HTML = r"""<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<title>C-1 step 1 — Biot field setup verification</title>
<style>
  :root { --bg:#0f1115; --fg:#e6e8ee; --dim:#9aa3b2; --line:#262b36; --accent:#67b7ff; }
  body { margin:0; background:var(--bg); color:var(--fg);
         font:14px/1.5 ui-sans-serif,system-ui,-apple-system,"Segoe UI",sans-serif; }
  header { padding:18px 24px 6px; border-bottom:1px solid var(--line); }
  h1 { font-size:17px; margin:0 0 4px; font-weight:600; }
  .sub { color:var(--dim); font-size:12px; }
  .wrap { display:flex; gap:24px; padding:18px 24px; flex-wrap:wrap; }
  .panel { background:#151822; border:1px solid var(--line); border-radius:8px; padding:14px; }
  canvas { display:block; image-rendering:pixelated; border-radius:4px; }
  .ctl { display:flex; gap:16px; align-items:center; padding:0 24px 10px; flex-wrap:wrap; }
  label { color:var(--dim); font-size:12px; }
  input[type=range] { width:220px; accent-color:var(--accent); }
  table { border-collapse:collapse; font-size:12px; }
  td,th { padding:3px 10px 3px 0; text-align:left; }
  th { color:var(--dim); font-weight:500; }
  .mono { font-family:ui-monospace,SFMono-Regular,Menlo,monospace; }
  .note { color:var(--dim); font-size:11.5px; max-width:640px; padding:0 24px 20px; }
</style></head><body>
<header>
  <h1>C-1 step 1 — does the poroelastic clock carry a length scale?</h1>
  <div class="sub">Biot pore-pressure field alone: no mechanics, no contact, no cell.
    Mid-plane <span class="mono">z</span> slice at true scale. build <span class="mono" id="build"></span></div>
</header>
<div class="ctl">
  <label>blob width a <select id="run"></select></label>
  <label>time <input type="range" id="frame" min="0" value="0"> <span class="mono" id="tlab"></span></label>
</div>
<div class="wrap">
  <div class="panel">
    <canvas id="field" width="420" height="420"></canvas>
    <div class="sub" style="margin-top:8px">
      dashed = Gaussian width <span class="mono">a</span> &nbsp;·&nbsp; solid = contact patch of the
      <span class="mono">patch_mean</span> observer
    </div>
  </div>
  <div class="panel">
    <canvas id="curves" width="460" height="300"></canvas>
    <div class="sub" style="margin-top:8px">observer readings, normalised; horizontal line = 1/e</div>
  </div>
  <div class="panel">
    <table id="stats"></table>
  </div>
</div>
<p class="note">The far field is empty by construction, so only the central crop is drawn — the box is
larger than what is shown, and its size is in the table. Each frame is normalised to its own maximum,
so a blob that spreads looks like it spreads rather than merely dims.</p>
<script>
const DATA = __PAYLOAD__;
document.getElementById('build').textContent = DATA.build.slice(0,8);
const runSel = document.getElementById('run'), fr = document.getElementById('frame');
DATA.runs.forEach((r,i)=>{ const o=document.createElement('option'); o.value=i;
  o.textContent = `a = ${r.a} µm`; runSel.appendChild(o); });

function colour(v){                       // diverging, zero at mid grey; positive warm
  const t = Math.max(-1, Math.min(1, v));
  if (t >= 0){ const u=t; return [15+240*u, 17+150*u, 21+40*u]; }
  const u=-t; return [15+40*u, 17+110*u, 21+230*u];
}
function drawField(){
  const r = DATA.runs[+runSel.value], f = r.slices[+fr.value];
  const cv = document.getElementById('field'), cx = cv.getContext('2d');
  const g = f.g, n = g.length, px = cv.width/n;
  const img = cx.createImageData(n,n);
  for(let j=0;j<n;j++) for(let i=0;i<n;i++){
    const [R,G,B]=colour(g[i][j]); const k=4*(j*n+i);
    img.data[k]=R; img.data[k+1]=G; img.data[k+2]=B; img.data[k+3]=255;
  }
  const tmp=document.createElement('canvas'); tmp.width=n; tmp.height=n;
  tmp.getContext('2d').putImageData(img,0,0);
  cx.imageSmoothingEnabled=false; cx.clearRect(0,0,cv.width,cv.height);
  cx.drawImage(tmp,0,0,cv.width,cv.height);
  const perUm = cv.width/(2*f.half_um), c = cv.width/2;
  cx.strokeStyle='#ffffff'; cx.setLineDash([5,4]); cx.lineWidth=1.2;
  cx.beginPath(); cx.arc(c,c,r.a*perUm,0,7); cx.stroke();
  cx.setLineDash([]); cx.strokeStyle='#67b7ff'; cx.lineWidth=1.6;
  cx.beginPath(); cx.arc(c,c,r.patch_um*perUm,0,7); cx.stroke();
  document.getElementById('tlab').textContent = `${(f.t*1e3).toFixed(2)} ms`;
}
function drawCurves(){
  const r = DATA.runs[+runSel.value], cv=document.getElementById('curves'), cx=cv.getContext('2d');
  const W=cv.width, H=cv.height, m={l:48,r:12,t:12,b:34};
  cx.clearRect(0,0,W,H);
  const t = r.curve.t_s, tmax=t[t.length-1];
  const keys=['l2_norm','peak','patch_mean'], cols=['#67b7ff','#ffb86b','#9be27a'];
  const X=v=>m.l+(W-m.l-m.r)*v/tmax, Y=v=>m.t+(H-m.t-m.b)*(1-(Math.log10(Math.max(v,1e-3))+3)/3);
  cx.strokeStyle='#262b36'; cx.beginPath(); cx.moveTo(m.l,m.t); cx.lineTo(m.l,H-m.b);
  cx.lineTo(W-m.r,H-m.b); cx.stroke();
  cx.strokeStyle='#4a5568'; cx.setLineDash([4,4]); cx.beginPath();
  cx.moveTo(m.l,Y(1/Math.E)); cx.lineTo(W-m.r,Y(1/Math.E)); cx.stroke(); cx.setLineDash([]);
  keys.forEach((k,i)=>{
    const y=r.curve[k]; if(!y) return; const y0=y[0];
    cx.strokeStyle=cols[i]; cx.lineWidth=1.6; cx.beginPath();
    y.forEach((v,j)=>{ const p=[X(t[j]),Y(v/y0)]; j?cx.lineTo(...p):cx.moveTo(...p); }); cx.stroke();
    const tau=r.tau[k];
    if(tau!=null){ cx.strokeStyle=cols[i]; cx.globalAlpha=.45; cx.beginPath();
      cx.moveTo(X(tau),m.t); cx.lineTo(X(tau),H-m.b); cx.stroke(); cx.globalAlpha=1; }
  });
  cx.fillStyle='#9aa3b2'; cx.font='11px ui-monospace,monospace';
  cx.fillText('t [ms]', W/2-18, H-10);
  cx.fillText(`0`, m.l-6, H-m.b+14); cx.fillText(`${(tmax*1e3).toFixed(1)}`, W-m.r-22, H-m.b+14);
  keys.forEach((k,i)=>{ cx.fillStyle=cols[i]; cx.fillText(k, W-m.r-96, m.t+14+14*i); });
}
function drawStats(){
  const r = DATA.runs[+runSel.value];
  const rows = [
    ['a (Gaussian width)', `${r.a} µm`],
    ['box L', `${r.L} µm`],
    ['dx', `${r.dx} µm`],
    ['patch radius', `${r.patch_um.toFixed(3)} µm`],
    ['tau / tau_box', r.box_frac.toExponential(3)],
    ['pore-mass drift', r.mass_drift.toExponential(3)],
    ['&nbsp;', '&nbsp;'],
    ['tau l2_norm', `${r.tau.l2_norm.toExponential(5)} s`],
    ['&nbsp;&nbsp;oracle', `${r.tau_oracle.l2_norm.toExponential(5)} s`],
    ['tau peak', `${r.tau.peak.toExponential(5)} s`],
    ['&nbsp;&nbsp;oracle', `${r.tau_oracle.peak.toExponential(5)} s`],
    ['tau patch_mean', `${r.tau.patch_mean.toExponential(5)} s`],
    ['&nbsp;&nbsp;oracle', 'none (finite contact)'],
    ['&nbsp;', '&nbsp;'],
    ['<b>l2 / peak</b>', `<b>${(r.tau.l2_norm/r.tau.peak).toFixed(4)}</b>`],
    ['closed form', '2.9477'],
  ];
  document.getElementById('stats').innerHTML =
    '<tr><th>quantity</th><th>value</th></tr>' +
    rows.map(([a,b])=>`<tr><td>${a}</td><td class="mono">${b}</td></tr>`).join('');
}
function sync(){ const r=DATA.runs[+runSel.value];
  fr.max=r.slices.length-1; if(+fr.value>+fr.max) fr.value=fr.max;
  drawField(); drawCurves(); drawStats(); }
runSel.onchange=sync; fr.oninput=drawField;
sync();
</script></body></html>
"""


def main() -> None:
    ap = argparse.ArgumentParser(description="Render a C-1 field record to a self-contained HTML page.")
    ap.add_argument("--record", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args()
    build(a.record, a.out)
    print(f"[c1-render] wrote {a.out} ({a.out.stat().st_size / 1e6:.2f} MB)")


if __name__ == "__main__":
    main()
