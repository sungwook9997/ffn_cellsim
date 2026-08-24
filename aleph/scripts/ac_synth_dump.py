"""Synthetic full-graph v2 dump — proves the staged-build viewer's component/connector fan-out END-TO-END.

The real ``ac/engine`` graph is 14 components / 36 connectors, but only the incumbent slice (cortex + NMII +
nucleus + membrane + LINC) has a real gbook dump today. This generator writes a v2 dump for the ENTIRE graph so
the viewer's per-component isolation, per-connector load-path, and cumulative-composition scenes can be
rendered + browser-checked NOW — before gbook binds each slice. It is honest about being a **schema/render
proof, NOT physics**: the TOPOLOGY (component names, roles, connector families, endpoint components) is the REAL
``reference_cell_architecture()`` graph (CPU-safe, no CUDA); only the node GEOMETRY and the |F|/load fields are
synthetic placeholders. Node counts are deliberately small (fast local render), never a native claim.

Everything goes through :func:`dump_state.write_v2_dump` — the SAME writer the real gbook producers use — so a
synthetic dump is byte-schema-identical to a real one and the viewer never special-cases it.

    python -m aleph.scripts.ac_synth_dump --out aleph/outputs/ac/cell_assembled/synth_graph_v2.npz
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

from aleph.engine.contracts import EvidenceLabel, EvidenceRung, QuantitativeClaim
from aleph.scripts.dump_state import _c, _k, write_v2_dump
from aleph.scripts.ff_cell_morphology import uv_sphere

R = 7.5                                          # µm cell radius
RNG = np.random.default_rng(7)


def _radial_filaments(n: int, r0: float, r1: float, jitter: float = 0.15, seed: int = 0):
    """n short radial filaments (2 nodes each) between radii r0→r1 on random directions; + a radial |F| ramp."""
    rng = np.random.default_rng(seed)
    d = rng.normal(size=(n, 3)); d /= np.linalg.norm(d, axis=1, keepdims=True)
    a = d * (r0 + rng.uniform(-jitter, jitter, n)[:, None])
    b = d * (r1 + rng.uniform(-jitter, jitter, n)[:, None])
    pos = np.stack([a, b], axis=1).reshape(-1, 3)                       # (2n,3)
    seg = np.arange(2 * n).reshape(n, 2)
    rr = np.linalg.norm(pos, axis=1)
    fmag = 10.0 ** (np.interp(rr, [r0 - 1, r1 + 1], [-2.0, 2.5]))       # log ramp 0.01→~300 pN toward the shell
    aid = np.repeat(np.arange(n), 2)
    return pos.astype(np.float32), seg, fmag, aid


def _aster(n: int, r1: float, center=(0.0, 0.0, 0.0), seed: int = 0):
    """n rods from a central point outward (MT aster / filopodial bundle)."""
    return _radial_filaments(n, 0.3, r1, jitter=0.05, seed=seed)


def build_components():
    comps = []
    # membrane — outer Helfrich sheet (mesh)
    mv, mf = uv_sphere(np.zeros(3), R, nu=40, nv=22)
    comps.append(_c("membrane", "plasma membrane (Helfrich+γ_mem)", "mesh", mv, faces=mf, order=0,
                    meta={"n_actors": 1}))
    # cortex — 2000 synthetic radial F-actin near the shell
    pos, seg, fmag, aid = _radial_filaments(2000, R - 0.5, R - 0.05, seed=1)
    comps.append(_c("cortex", "cortex F-actin (synthetic)", "filaments", pos, seg=seg, fmag=fmag,
                    actor_id=aid, order=1, meta={"n_actors": 2000}))
    # nucleus — inner core (mesh)
    nv, nf = uv_sphere(np.zeros(3), 5.0, nu=28, nv=16)
    comps.append(_c("nucleus", "deformable nucleus", "mesh", nv, faces=nf, order=2, meta={"n_actors": 1}))
    # sf_arc — long ventral stress fibers crossing the cell
    rng = np.random.default_rng(3)
    a = rng.normal(size=(60, 3)); a /= np.linalg.norm(a, axis=1, keepdims=True); a *= (R - 1.0)
    b = -a + rng.normal(scale=0.5, size=(60, 3))
    sf_pos = np.stack([a, b], axis=1).reshape(-1, 3); sf_seg = np.arange(120).reshape(60, 2)
    comps.append(_c("sf_arc", "stress fibers / arcs (synthetic)", "filaments", sf_pos, seg=sf_seg,
                    fmag=np.full(120, 50.0), actor_id=np.repeat(np.arange(60), 2), order=3,
                    meta={"n_actors": 60}))
    # microtubule — aster from the centrosome
    mtp, mts, mtf, mta = _aster(120, R - 0.3, seed=4)
    comps.append(_c("microtubule", "microtubule aster (synthetic)", "filaments", mtp, seg=mts, fmag=mtf,
                    actor_id=mta, order=4, meta={"n_actors": 120}))
    # intermediate_filament — perinuclear cage
    ifp, ifs, iff, ifa = _radial_filaments(150, 5.0, 6.5, jitter=0.3, seed=5)
    comps.append(_c("intermediate_filament", "IF perinuclear cage (synthetic)", "filaments", ifp, seg=ifs,
                    fmag=iff, actor_id=ifa, order=5, meta={"n_actors": 150}))
    # lamellipodium — flat protrusion sheet near +x pole
    lp, ls, lf, la = _radial_filaments(200, R - 0.2, R + 1.2, jitter=0.4, seed=6)
    keep = lp.reshape(-1, 2, 3)[:, 1, 0] > R * 0.5                       # keep only the +x-leaning ones
    lp2 = lp.reshape(-1, 2, 3)[keep].reshape(-1, 3)
    comps.append(_c("lamellipodium", "lamellipodium (synthetic)", "filaments", lp2,
                    seg=np.arange(lp2.shape[0]).reshape(-1, 2), fmag=np.full(lp2.shape[0], 20.0),
                    actor_id=np.repeat(np.arange(lp2.shape[0] // 2), 2), order=6,
                    meta={"n_actors": lp2.shape[0] // 2}))
    # filopodium — a few finger bundles beyond the shell at +y
    fp = []
    for k in range(8):
        dirv = np.array([RNG.uniform(-0.3, 0.3), 1.0, RNG.uniform(-0.3, 0.3)]); dirv /= np.linalg.norm(dirv)
        fp += [dirv * (R - 0.3), dirv * (R + 2.0)]
    fp = np.asarray(fp, np.float32)
    comps.append(_c("filopodium", "filopodia (synthetic)", "filaments", fp,
                    seg=np.arange(fp.shape[0]).reshape(-1, 2), fmag=np.full(fp.shape[0], 8.0),
                    actor_id=np.repeat(np.arange(8), 2), order=7, meta={"n_actors": 8}))
    # ecm — collagen fibers OUTSIDE the cell
    ep, es, ef, ea = _radial_filaments(300, R + 0.5, R + 4.0, jitter=1.0, seed=8)
    comps.append(_c("ecm", "collagen ECM (synthetic)", "filaments", ep, seg=es,
                    fmag=np.full(ep.shape[0], 5.0), actor_id=ea, order=8, meta={"n_actors": 300}))
    return comps


def _nearest_pairs(pa, pb, n, seed=0):
    """n synthetic joints: for n random nodes of A, connect to the nearest node of B → world endpoints."""
    rng = np.random.default_rng(seed)
    ia = rng.choice(pa.shape[0], size=min(n, pa.shape[0]), replace=False)
    out = []
    for i in ia:
        j = int(np.argmin(np.linalg.norm(pb - pa[i], axis=1)))
        out += [pa[i], pb[j]]
    return np.asarray(out, np.float32)


def build_connectors(comps):
    by = {c["key"]: c["pos"] for c in comps}
    # (key, family, a, b, n_joints) — families/endpoints mirror reference_cell_architecture()
    specs = [
        ("membrane_erm_cortex", "ERM", "membrane", "cortex", 120),
        ("integrin_collagen_clutch", "FA-clutch", "cortex", "ecm", 60),
        ("actin_cap_linc", "LINC", "sf_arc", "nucleus", 40),
        ("mt_nucleus_linc", "LINC", "microtubule", "nucleus", 30),
        ("if_nucleus_linc", "LINC", "intermediate_filament", "nucleus", 30),
        ("nmii_sf_motor", "MOTOR", "sf_arc", "cortex", 40),
        ("mt_cortex_capture", "MOTOR", "microtubule", "cortex", 30),
        ("if_sf_plectin", "plectin", "intermediate_filament", "sf_arc", 30),
        ("mt_sf_spectraplakin", "spectraplakin", "microtubule", "sf_arc", 25),
        ("sf_cortex_transient", "transient-actin", "sf_arc", "cortex", 30),
        ("lamellipodium_cortex_seam", "transient-actin", "lamellipodium", "cortex", 25),
        ("filopodium_cortex_root", "transient-actin", "filopodium", "cortex", 12),
        ("sf_cytosol_transfer", "immersed-transfer", "sf_arc", "nucleus", 20),
        ("membrane_ecm_contact", "CONTACT", "membrane", "ecm", 40),
    ]
    cons = []
    for key, fam, a, b, n in specs:
        if a not in by or b not in by:
            continue
        ep = _nearest_pairs(by[a], by[b], n, seed=hash(key) % 2**31)
        load = RNG.uniform(1.0, 40.0, ep.shape[0] // 2)
        cons.append(_k(key, fam, key.replace("_", " "), a, b, ep, load=load))
    return cons


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=str,
                    default="aleph/outputs/ac/cell_assembled/synth_graph_v2.npz")
    a = ap.parse_args()
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    comps = build_components()
    cons = build_connectors(comps)
    # Counts come from what was actually built, never from a literal — the same drift that left
    # "13 comp / 32 conn" in a render title for three days after connectors 32 -> 35 landed.
    evidence = EvidenceLabel(
        rung=EvidenceRung.SCHEMA_PROOF,
        quantitative=QuantitativeClaim.BLOCKED,
        basis=f"topology only: {len(comps)} components / {len(cons)} connectors, no physics evaluated",
    )
    meta = {"stage": {"name": f"SYNTHETIC full-graph ({len(comps)} comp / real topology)",
                      **evidence.as_artifact_fields(),
                      "note": "topology real (reference_cell_architecture); geometry+forces synthetic — NOT physics"},
            "report": {}, "ledger": {}}
    out = write_v2_dump(a.out, comps, cons, meta)
    print(f"wrote {out}  ({len(comps)} components, {len(cons)} connectors)")


if __name__ == "__main__":
    main()
