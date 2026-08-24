"""FF ECM — interactive HTML viewer for the contact-guidance result (ROADMAP NEAR #6).

Renders the ALIGNED collagen-I microstructure across the nematic-order ladder (S = 0 → 0.83) as ONE
multi-scene, full-resolution interactive HTML (per the no-downsample rule), each scene labelled with its
measured S, the native ECM virial stress anisotropy R_σ=σ∥/σ⊥, and the library E∥/E⊥ — so the reader can SEE
the fibers rotate into the director as S rises and read off the emergent directional stiffness. A bright green
line marks the fixed in-plane director (the fiber-alignment axis). Complements the comparison figure
`figs/contact_guidance_anisotropy_cg_native.png` (the numbers) with the morphology (the mechanism).

Output (standard ecm_lib location): aleph/outputs/ff/ecm_lib/contact_guidance_viewer.html

Run:  python -m aleph.scripts.ff_contact_guidance_viewer            # default S ladder, ~40 µm REV
      python -m aleph.scripts.ff_contact_guidance_viewer --box 50   # bigger REV
"""

from __future__ import annotations

import argparse
import json
import os

import numpy as np

from aleph.laws import ecm_library as L
from aleph.scripts.ff_ecm_viewer import ecm_layers
from aleph.scripts.ff_viewer_html import build_viewer

OUT = "aleph/outputs/ff/ecm_lib"

# library-validated ladder (README): S → E∥/E⊥ = 1 / 3.1 / 10.3 / 63
E_ANISO_LIB = {0.0: 1.0, 0.30: 3.1, 0.59: 10.3, 0.83: 63.0}


def _native_rsigma() -> dict:
    """R_σ(S) mean from the committed native #6 JSON, keyed by rounded S (for scene annotation)."""
    path = f"{OUT}/contact_guidance_cg_native.json"
    out = {}
    if not os.path.exists(path):
        return out
    d = json.load(open(path))
    rows = d.get("rows", [])
    by = {}
    for r in rows:
        s = round(float(r["S_target"]), 2)
        v = r.get("R_sigma")
        if v is not None and np.isfinite(v):
            by.setdefault(s, []).append(float(v))
    for s, vs in by.items():
        out[s] = float(np.mean(vs))
    return out


def director_line(lo, hi, director=(1.0, 0.0, 0.0)) -> dict:
    """A bright line through the box centre along the (in-plane) director — the fiber-alignment axis."""
    c = 0.5 * (np.asarray(lo, float) + np.asarray(hi, float))
    d = np.asarray(director, float)
    d = d / (np.linalg.norm(d) or 1.0)
    half = 0.48 * float(np.min(np.asarray(hi, float) - np.asarray(lo, float)))
    seg = np.array([[c - half * d, c + half * d]], float)          # (1,2,3)
    return {"name": "director (fiber axis)", "kind": "lines", "verts": seg, "color": "#33ff88",
            "opacity": 1.0, "size": 5.0, "on_top": True}


def build(box: float = 40.0, conc: float = 1.5, S_list=(0.0, 0.30, 0.59, 0.83),
          director=(1.0, 0.0, 0.0), out: str | None = None) -> str:
    """Build the multi-scene contact-guidance viewer (aligned collagen across S)."""
    os.makedirs(OUT, exist_ok=True)
    out = out or f"{OUT}/contact_guidance_viewer.html"
    spec = L.get_spec("collagen_I")
    lo, hi = [0.0, 0.0, 0.0], [box, box, box]
    rsig = _native_rsigma()
    scenes = {}
    for S in S_list:
        rng = np.random.default_rng(7)
        ecm = L.build_fibrillar_ecm(spec, lo, hi, concentration=conc, dim=3, alignment_S=float(S),
                                    director=director, pin_faces=(), target_z=3.2, rng=rng)
        Sm = float(getattr(ecm, "S_measured", S))
        Elib = E_ANISO_LIB.get(round(float(S), 2))
        Rn = rsig.get(round(float(S), 2))
        tag = (f"collagen-I  S_target={S:g}  (measured {Sm:.2f})   "
               f"E∥/E⊥={'%.0f×' % Elib if Elib else '—'}   "
               f"native R_σ={('%.1f×' % Rn) if Rn else '—'}")
        layers = ecm_layers(ecm, name=f"collagen S={S:g}")
        layers.append(director_line(lo, hi, director))
        scenes[tag] = layers
    title = ("FF ECM NEAR #6 — contact guidance: aligned collagen-I microstructure vs nematic order S  "
             "(green = fiber-alignment director; R_σ = emergent ECM virial stress anisotropy, native)")
    build_viewer(scenes, out, title=title)
    print(f"wrote {out}  ({len(scenes)} scenes: S={list(S_list)})")
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--box", type=float, default=40.0, help="REV box side [µm] for the microstructure view")
    ap.add_argument("--conc", type=float, default=1.5, help="collagen concentration [mg/mL]")
    ap.add_argument("--S", default="0.0,0.30,0.59,0.83", help="comma-separated nematic order values")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    S_list = tuple(float(x) for x in a.S.split(",") if x.strip())
    build(box=a.box, conc=a.conc, S_list=S_list, out=a.out)


if __name__ == "__main__":
    main()
