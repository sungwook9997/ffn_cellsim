"""Viz-owned state dumper → the generic v2 component/connector schema the staged-build viewer reads.

This is the WRITER counterpart of :mod:`ac_viz_common` (the reader). It is the SINGLE source of the v2 dump
schema: :func:`write_v2_dump` serialises a list of components + connectors (plain numpy) into ONE ``.npz``
with a ``manifest_json`` that names every flat array. Both the real gbook producers below and the local
synthetic generator (``ac_synth_dump.py``) go through this one writer, so a synthetic dump and a real dump are
byte-schema-identical and the viewer never special-cases them.

Producers:

* :func:`legacy_npz_to_v2` — re-express the incumbent ``ac/cell/dump_state.py`` output
  (``assembled_state.npz``: cortex F-actin + NMII + deformable nucleus + membrane + LINC) as the generic
  component/connector model with **stable per-component keys + GLOBAL-UNIQUE actor IDs** (each physical
  filament belongs to exactly ONE component → no double-draw). Pure numpy — runs anywhere, incl. the dev Mac.
* :func:`dump_incumbent_cell` — gbook A5000: build the incumbent ``build_cell`` on CUDA, run one transactional
  ``--from-resting`` outer step via the frozen driver, dump the legacy ``.npz`` (reusing the tested
  ``ac.cell.dump_state``), then convert to v2. READ-ONLY over physics.
* :func:`dump_engine_actor` — gbook A5000, FORWARD path: iterate the ``ac/engine``
  ``reference_cell_architecture()`` graph and read each bound component runtime's ``position_d`` / ``force_d`` /
  segments / ``persistent_filament_id_d`` and each connector's resolved endpoints. Best-effort over whichever
  components/connectors are already bound (the engine is mid-binding). READ-ONLY.

I0-A: the CUDA producers run on the gbook A5000; the dev Mac cannot launch the kernels. All three emit a plain
``.npz`` the Mac renders with numpy only. This module edits no ``ac/`` / ``ff/`` / solver / scheduler file.

    # gbook (CUDA): incumbent build_cell → v2 dump
    ~/miniconda3/envs/ffn_sim/bin/python -m aleph.scripts.dump_state --incumbent \
        --out aleph/outputs/ac/cell_assembled/assembled_state_v2.npz
    # dev Mac: convert an existing legacy dump → v2 (no CUDA)
    python -m aleph.scripts.dump_state --from-legacy \
        aleph/outputs/ac/cell_assembled/assembled_state.npz --out .../assembled_state_v2.npz
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

__all__ = ["write_v2_dump", "legacy_npz_to_v2", "dump_incumbent_cell", "dump_engine_actor",
           "dump_ecm_network", "ACTOR_ID_STRIDE"]

# Global-unique actor-ID namespace: component of binding-order k owns IDs in [k·STRIDE, (k+1)·STRIDE). Guarantees
# each physical filament belongs to exactly ONE component (cortex/SF/… disjoint) so nothing is double-drawn.
ACTOR_ID_STRIDE = 1_000_000_000


def _c(key, label, kind, pos, *, seg=None, faces=None, fmag=None, actor_id=None, order=0,
       color=None, array_meta=None, meta=None) -> dict:
    """Build one component spec (numpy arrays kept alongside the manifest entry until write)."""
    return {"key": key, "label": label, "kind": kind, "order": int(order), "color": color,
            "pos": np.asarray(pos, np.float32).reshape(-1, 3),
            "seg": None if seg is None else np.asarray(seg, np.int64).reshape(-1, 2),
            "faces": None if faces is None else np.asarray(faces, np.int64).reshape(-1, 3),
            "fmag": None if fmag is None else np.asarray(fmag, np.float64).reshape(-1),
            "actor_id": None if actor_id is None else np.asarray(actor_id, np.int64).reshape(-1),
            "array_meta": array_meta or {}, "meta": meta or {}}


def _k(key, family, label, a, b, endpoints, *, load=None, color=None, meta=None) -> dict:
    """Build one connector spec (world-space endpoint pairs [a0,b0,a1,b1,…])."""
    return {"key": key, "family": family, "label": label, "a_component": a, "b_component": b, "color": color,
            "endpoints": np.asarray(endpoints, np.float32).reshape(-1, 3),
            "load": None if load is None else np.asarray(load, np.float64).reshape(-1), "meta": meta or {}}


def write_v2_dump(out_npz, components: list[dict], connectors: list[dict], meta: dict) -> str:
    """Serialise components + connectors into ONE v2 ``.npz`` (manifest_json + flat prefixed arrays).

    ``components`` items come from :func:`_c`, ``connectors`` from :func:`_k`. Global-unique actor IDs are
    assigned per binding order if a component omits them. Returns the written path.
    """
    arrays: dict[str, np.ndarray] = {}
    man_components: list[dict] = []
    for i, c in enumerate(sorted(components, key=lambda c: c["order"])):
        key = c["key"]
        pre = f"cmp_{key}"
        arrays[f"{pre}_pos"] = c["pos"]
        entry = {"key": key, "label": c["label"], "kind": c["kind"], "order": c["order"],
                 "color": c["color"], "pos_key": f"{pre}_pos", "meta": c["meta"]}
        if c["seg"] is not None:
            arrays[f"{pre}_seg"] = c["seg"]; entry["seg_key"] = f"{pre}_seg"
        if c["faces"] is not None:
            arrays[f"{pre}_faces"] = c["faces"]; entry["faces_key"] = f"{pre}_faces"
        if c["fmag"] is not None:
            arrays[f"{pre}_fmag"] = c["fmag"]; entry["fmag_key"] = f"{pre}_fmag"
        aid = c["actor_id"]
        if aid is None:                                            # no per-node actor map: namespace per node
            aid = i * ACTOR_ID_STRIDE + np.arange(c["pos"].shape[0], dtype=np.int64)
            # actor count is a declared aggregate (minifilaments, 1 surface, …), not the node count
            n_actors = int(c["meta"].get("n_actors", c["pos"].shape[0]))
        else:
            aid = i * ACTOR_ID_STRIDE + aid.astype(np.int64)       # namespace the local ids → global-unique
            n_actors = int(np.unique(aid).size)
        arrays[f"{pre}_actorid"] = aid; entry["actor_id_key"] = f"{pre}_actorid"
        entry["n_actors"] = n_actors
        am: dict[str, str] = {}
        for field_name, arr in c["array_meta"].items():
            ak = f"{pre}_{field_name}"
            arrays[ak] = np.asarray(arr); am[field_name] = ak
        if am:
            entry["array_meta"] = am
        entry["meta"] = {**c["meta"], "n_actors": n_actors}
        man_components.append(entry)

    man_connectors: list[dict] = []
    for c in connectors:
        pre = f"con_{c['key']}"
        arrays[f"{pre}_xyz"] = c["endpoints"]
        entry = {"key": c["key"], "family": c["family"], "label": c["label"], "color": c["color"],
                 "a_component": c["a_component"], "b_component": c["b_component"],
                 "endpoints_key": f"{pre}_xyz", "meta": c["meta"]}
        if c["load"] is not None:
            arrays[f"{pre}_load"] = c["load"]; entry["load_key"] = f"{pre}_load"
        man_connectors.append(entry)

    # disjointness assertion: no actor-id shared across components (the no-double-count invariant)
    seen: dict[int, str] = {}
    for e in man_components:
        for a in np.unique(arrays[e["actor_id_key"]]).tolist():
            if a in seen and seen[a] != e["key"]:
                raise ValueError(f"actor id {a} shared by {seen[a]} and {e['key']} — components must be disjoint")
            seen[a] = e["key"]

    manifest = {"schema_version": 2, "components": man_components, "connectors": man_connectors,
                "meta": meta, "report": meta.get("report", {}), "ledger": meta.get("ledger", {}),
                "stage": meta.get("stage", {})}
    out_npz = str(out_npz)
    np.savez_compressed(out_npz, manifest_json=np.array(json.dumps(manifest), dtype=object), **arrays)
    with open(out_npz.rsplit(".npz", 1)[0] + ".json", "w") as fh:
        json.dump({"stage": manifest["stage"], "report": manifest["report"],
                   "components": [{k: e[k] for k in ("key", "kind", "order", "n_actors")} for e in man_components],
                   "connectors": [{k: e[k] for k in ("key", "family", "a_component", "b_component")}
                                  for e in man_connectors]}, fh, indent=2, default=float)
    return out_npz


# ──────────────────────────────── producer 1: legacy npz → v2 (pure numpy) ────────────────────────────────
def legacy_npz_to_v2(legacy_npz, out_npz, *, stage: dict | None = None) -> str:
    """Convert the incumbent ``assembled_state.npz`` → the generic v2 schema (runs anywhere, no CUDA)."""
    from aleph.scripts.ac_viz_common import normalize_legacy
    d = np.load(str(legacy_npz), allow_pickle=True)
    dump = normalize_legacy(d)
    comps: list[dict] = []
    for c in dump.components:
        am = {}
        if "f_steric" in c.meta:
            am["f_steric"] = c.meta["f_steric"]
        if "head_verts" in c.meta:
            am["head_verts"] = c.meta["head_verts"]
        comps.append(_c(c.key, c.label, c.kind, c.pos, seg=c.seg, faces=c.faces, fmag=c.fmag,
                        actor_id=c.actor_id, order=c.order, color=c.color, array_meta=am,
                        meta={k: v for k, v in c.meta.items() if not isinstance(v, np.ndarray)}))
    cons = [_k(cn.key, cn.family, cn.label, cn.a_component, cn.b_component, cn.endpoints,
               load=cn.load, color=cn.color) for cn in dump.connectors]
    meta = dict(dump.meta)
    meta["stage"] = stage or dump.meta.get("stage", {"name": "incumbent build_cell", "schema": "legacy→v2"})
    return write_v2_dump(out_npz, comps, cons, meta)


# ──────────────────────────── producer 2: incumbent build_cell → v2 (gbook CUDA) ──────────────────────────
def dump_incumbent_cell(out_npz, *, n_filaments: int = 70686, n_inner: int = 60, dt_phys: float = 0.05,
                        with_myosin: bool = True, with_steric: bool = True, with_pressure: bool = True,
                        device: str | None = None, stage: dict | None = None) -> str:
    """gbook: build the incumbent composed cell on CUDA, run one resting step, dump v2 (READ-ONLY physics)."""
    import tempfile

    from aleph.components.incumbent.assemble import CellConfig
    from aleph.components.incumbent.dump_state import dump_state as _legacy_dump   # tested gbook step; writes legacy npz

    cfg = CellConfig(n_filaments=n_filaments, with_myosin=with_myosin, with_steric=with_steric,
                     with_pressure=with_pressure, device=device)
    tmp = Path(tempfile.mkdtemp(prefix="acdump_")) / "legacy.npz"
    report = _legacy_dump(cfg, str(tmp), n_inner=n_inner, dt_phys=dt_phys)
    st = stage or {"name": f"incumbent build_cell (N={n_filaments:,})",
                   "evidence": "CONNECTED", "accepted": report.get("outer_accepted")}
    return legacy_npz_to_v2(tmp, out_npz, stage=st)


# ─────────────────────────── producer 3: ac/engine actor graph → v2 (gbook CUDA) ──────────────────────────
def dump_engine_actor(actor, arch, out_npz, *, stage: dict | None = None) -> str:
    """gbook FORWARD path: read whatever ``ac/engine`` components/connectors are bound on ``actor``.

    Topology (component/connector names, families, endpoint components) comes from ``arch`` (CPU-safe
    ``reference_cell_architecture()``); per-node arrays come from each bound runtime's ``position_d`` /
    ``force_d`` / segments / ``persistent_filament_id_d`` via ``.numpy()``. Components/connectors that are not
    yet bound are skipped (best-effort — the engine is mid-binding), so this grows automatically as slices land.
    """
    def _np(x):
        return None if x is None else np.asarray(x.numpy())

    def _runtime(name):
        for attr in ("component_runtime", "component", "get_component"):
            fn = getattr(actor, attr, None)
            if callable(fn):
                try:
                    return fn(name)
                except Exception:
                    pass
        return None

    def _first(rt, names):
        for n in names:
            v = getattr(rt, n, None)
            if v is not None:
                return v
        return None

    comps: list[dict] = []
    node_world: dict[str, np.ndarray] = {}                          # for resolving connector endpoints
    for order, cc in enumerate(arch.components):
        rt = _runtime(cc.name)
        pos_d = _first(rt, ("position_d", "pos_d")) if rt is not None else None
        if pos_d is None:
            continue                                               # geometry-less (focal_adhesion/world) or unbound
        pos = _np(pos_d).reshape(-1, 3)
        node_world[cc.name] = pos
        force_d = _first(rt, ("force_d",))
        fmag = np.linalg.norm(_np(force_d), axis=1) if force_d is not None else None
        seg_d = _first(rt, ("segments_d",))
        seg = _np(seg_d).reshape(-1, 2) if seg_d is not None else None
        off_d = _first(rt, ("fiber_offset_d", "fiber_offsets_d"))
        if seg is None and off_d is not None:
            from aleph.scripts.ac_viz_common import seg_pairs_from_offsets
            seg = seg_pairs_from_offsets(_np(off_d), pos.shape[0])
        faces_d = _first(rt, ("faces_d",))
        faces = None
        if faces_d is not None:
            faces = _np(faces_d).reshape(-1, 3)
            faces = faces - int(faces.min()) if faces.size else faces
        aid_d = _first(rt, ("persistent_filament_id_d",))
        actor_id = _np(aid_d) if aid_d is not None else None
        kind = "mesh" if faces is not None else ("filaments" if seg is not None else "points")
        comps.append(_c(cc.name, cc.name, kind, pos, seg=seg, faces=faces, fmag=fmag, actor_id=actor_id,
                        order=order, meta={"role": getattr(cc.role, "name", str(cc.role))}))

    cons: list[dict] = []
    for cc in arch.connectors:
        a, b = cc.component_a, cc.component_b
        if a not in node_world or b not in node_world:
            continue
        joint = None
        for attr in ("connector_runtime", "connector", "get_connector"):
            fn = getattr(actor, attr, None)
            if callable(fn):
                try:
                    joint = fn(cc.name)
                    break
                except Exception:
                    pass
        ep = None
        load = None
        if joint is not None:
            ai = _first(joint, ("endpoint_a_node_d", "filament_actor_id_d", "a_node_d"))
            bi = _first(joint, ("endpoint_b_node_d", "nucleus_actor_id_d", "b_node_d"))
            if ai is not None and bi is not None:
                ia, ib = _np(ai).reshape(-1).astype(np.int64), _np(bi).reshape(-1).astype(np.int64)
                m = min(ia.size, ib.size, node_world[a].shape[0], node_world[b].shape[0])
                if m:
                    pa = node_world[a][np.clip(ia[:m], 0, node_world[a].shape[0] - 1)]
                    pb = node_world[b][np.clip(ib[:m], 0, node_world[b].shape[0] - 1)]
                    ep = np.stack([pa, pb], axis=1).reshape(-1, 3)
            ld = _first(joint, ("load_d",))
            load = _np(ld) if ld is not None else None
        if ep is None:
            continue
        fam = getattr(cc.family, "name", str(cc.family)).replace("_", "-").replace("FA-CLUTCH", "FA-clutch")
        cons.append(_k(cc.name, fam, cc.name, a, b, ep, load=load))

    meta = {"stage": stage or {"name": "ac/engine actor graph", "evidence": "CONNECTED"},
            "report": {}, "ledger": {}}
    return write_v2_dump(out_npz, comps, cons, meta)


# ─────────────────────── producer 4: native ECM Mikado network → v2 (gbook CUDA) ──────────────────────────
def dump_ecm_network(out_npz, *, n_fibers: int = 220, fiber_length_um: float = 12.0,
                     target_segment_um: float = 1.5, box_half_um: float = 9.0, material: str = "collagen_I",
                     shear: float = 0.05, device: str | None = None, stage: dict | None = None) -> str:
    """gbook A5000: build the native ECM Mikado collagen network + its constitutive force, dump v2 (READ-ONLY).

    The ECM topology is CUDA-only (``ac.ecm.device_schema.require_cuda_device`` rejects every non-CUDA runtime),
    so this producer runs ONLY on the gbook; the dev Mac renders the emitted ``.npz`` with numpy. It builds the
    real :class:`~aleph.components.ecm.MikadoTopologyBuilder` SoA, binds the SOURCED collagen constitutive force
    (:func:`~aleph.engine.ecm_mechanics.build_collagen_constitutive_force`, EA/κ from the fibrillar
    collagen card), and reads two force states over the SAME network so the viewer shows the constitutive law:

      * RELAXED (reference config): segment lengths == rest, straight rods → restoring |F| ≈ 0.
      * SHEARED (a small ``shear`` simple-shear ``x += γ·z`` diagnostic perturbation): the stretched/bent
        segments develop a restoring |F| — the collagen network resisting deformation.

    The emitted ``ecm`` component carries the REFERENCE geometry (active segments) coloured by the SHEARED
    restoring |F| (log-turbo pN); the relaxed field (≈0) is stored alongside for provenance. No modulus band was
    tuned — EA/κ come straight from the card, and the shear is a documented diagnostic, not a fit.

    This performs host readbacks ONLY after the two completed diagnostic accumulations (for plotting); it is not
    a physical-time step and decides no device event. Runs on gbook; pending a gbook dump on the dev Mac.
    """
    import warp as wp

    from aleph.components.ecm import MikadoInitConfig, MikadoTopologyBuilder
    from aleph.engine.ecm_mechanics import build_collagen_constitutive_force, collagen_force_provenance
    from aleph.laws.ecm_library import get_spec

    dev = str(wp.get_device(device))
    cfg = MikadoInitConfig(
        box_lo_um=(-box_half_um, -box_half_um, -box_half_um),
        box_hi_um=(box_half_um, box_half_um, box_half_um),
        n_fibers=int(n_fibers), fiber_length_um=float(fiber_length_um),
        target_segment_um=float(target_segment_um), crosslink_capture_um=0.4,
        pin_faces=("x_lo", "x_hi", "y_lo", "y_hi", "z_lo", "z_hi"), pin_margin_um=0.25,
        rng_seed=1234, max_refinement_level=1, persistent_id_base=80_000)
    topology = MikadoTopologyBuilder(cfg, device=dev).initialize()
    card = get_spec(material)
    force_pass = build_collagen_constitutive_force(topology, card, device=dev)

    pos0 = np.asarray(topology.position_d.numpy(), np.float64).reshape(-1, 3)
    n_nodes = pos0.shape[0]
    seg = np.asarray(topology.segments_d.numpy(), np.int64).reshape(-1, 2)
    seg_active = np.asarray(topology.segment_active_d.numpy()).reshape(-1).astype(bool)
    active_seg = seg[seg_active]
    node_fiber = np.asarray(topology.node_fiber_d.numpy(), np.int64).reshape(-1)
    n_active_fibers = int(np.unique(node_fiber[np.unique(active_seg)]).size) if active_seg.size else 0

    def _restoring_mag(pos_np: np.ndarray) -> np.ndarray:
        pos_d = wp.array(np.ascontiguousarray(pos_np), dtype=wp.vec3d, device=dev)
        force_d = wp.zeros(n_nodes, dtype=wp.vec3d, device=dev)
        force_pass.accumulate(pos_d, force_d)
        wp.synchronize_device(dev)
        return np.linalg.norm(np.asarray(force_d.numpy(), np.float64).reshape(-1, 3), axis=1)

    f_relaxed = _restoring_mag(pos0)
    pos_sheared = pos0.copy()
    pos_sheared[:, 0] += float(shear) * pos0[:, 2]                 # simple shear x += γ·z (diagnostic)
    f_sheared = _restoring_mag(pos_sheared)

    st = stage or {
        "name": f"native ECM Mikado {card.name} (n_fibers={n_fibers})", "evidence": "KERNEL_BOUND",
        "ecm_shear_gamma": float(shear),
        "ecm_relaxed_maxF_pN": float(f_relaxed.max()) if f_relaxed.size else 0.0,
        "ecm_sheared_maxF_pN": float(f_sheared.max()) if f_sheared.size else 0.0,
        "provenance": collagen_force_provenance(card)}
    comp = _c("ecm", f"ECM collagen network ({card.name}) · restoring |F| @ γ={shear:g}", "filaments",
              pos0, seg=active_seg, fmag=f_sheared, actor_id=node_fiber, order=0, color="#8d6e63",
              array_meta={"f_relaxed": f_relaxed}, meta={"n_actors": n_active_fibers, "material": card.key})
    meta = {"stage": st, "report": {"n_active_fibers": n_active_fibers,
                                     "n_active_segments": int(active_seg.shape[0])}, "ledger": {}}
    return write_v2_dump(out_npz, [comp], [], meta)


def _build_argparser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Dump Active-Cell state → generic v2 component/connector schema.")
    src = p.add_mutually_exclusive_group(required=True)
    src.add_argument("--incumbent", action="store_true", help="gbook CUDA: build build_cell + one resting step")
    src.add_argument("--from-legacy", type=str, metavar="NPZ", help="convert an existing legacy npz (no CUDA)")
    src.add_argument("--engine", action="store_true", help="gbook CUDA: read the ac/engine bound actor graph")
    src.add_argument("--ecm", action="store_true", help="gbook CUDA: build the native ECM Mikado network + force")
    p.add_argument("--out", type=str, required=True, help="output v2 .npz path")
    p.add_argument("--n-filaments", type=int, default=70686, help="cortex F-actin count (full native=70686)")
    p.add_argument("--n-inner", type=int, default=60)
    p.add_argument("--dt-phys", type=float, default=0.05)
    p.add_argument("--no-myosin", action="store_true")
    p.add_argument("--no-steric", action="store_true")
    p.add_argument("--no-pressure", action="store_true")
    p.add_argument("--ecm-n-fibers", type=int, default=220, help="ECM collagen fiber count (--ecm)")
    p.add_argument("--ecm-shear", type=float, default=0.05, help="ECM diagnostic simple-shear γ (--ecm)")
    p.add_argument("--ecm-material", type=str, default="collagen_I", help="ECM fibrillar card key (--ecm)")
    p.add_argument("--device", type=str, default=None)
    return p


def main() -> None:
    a = _build_argparser().parse_args()
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    if a.from_legacy:
        out = legacy_npz_to_v2(a.from_legacy, a.out)
    elif a.incumbent:
        out = dump_incumbent_cell(a.out, n_filaments=a.n_filaments, n_inner=a.n_inner, dt_phys=a.dt_phys,
                                  with_myosin=not a.no_myosin, with_steric=not a.no_steric,
                                  with_pressure=not a.no_pressure, device=a.device)
    elif a.ecm:
        out = dump_ecm_network(a.out, n_fibers=a.ecm_n_fibers, shear=a.ecm_shear, material=a.ecm_material,
                               device=a.device)
    else:  # --engine
        from aleph.engine import reference_cell_architecture
        try:
            from aleph.engine import CellActor  # type: ignore
        except Exception:
            CellActor = None
        arch = reference_cell_architecture()
        actor = CellActor() if CellActor is not None else None
        if actor is None:
            raise SystemExit("--engine needs a bound CellActor; construct+bind it on gbook then call "
                             "dump_engine_actor(actor, arch, out) directly.")
        out = dump_engine_actor(actor, arch, a.out)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
