"""Shared, data-driven scene model for the staged Active-Cell (ac/engine) visualisers.

The staged ``ac/engine`` build is a component/connector graph: components OWN a disjoint node population;
connectors are the ONLY mechanical joints (co-location in an array is NEVER a connection). This module turns a
state **dump** into interactive-viewer scenes *without hardcoding which components exist*, so the moment a new
component (SF, MT, IF, ECM, lamellipodium, filopodium, …) or a new connector family (ERM, FA-clutch, LINC,
MOTOR, plectin, spectraplakin, transient-actin, immersed-transfer) appears in the dump, it gets an isolation
scene, a coloured load-path scene, and a slot in the cumulative-composition view — with no edit here.

Two dump schemas are normalised to the same in-memory model:

* **v2 (generic, written by ``scripts/dump_state.py``):** a ``manifest_json`` lists components + connectors, each
  naming the flat ``.npz`` arrays that hold its nodes / faces / per-node |F| / global-unique actor IDs /
  connector world-space endpoint pairs. Fully data-driven — this module iterates the manifest.
* **legacy (written by ``ac/cell/dump_state.py``):** the incumbent ``assembled_state.npz`` (cortex F-actin + NMII
  + deformable nucleus + membrane + LINC). :func:`normalize_legacy` re-expresses it as the same component/
  connector model so the generic renderer works on the real incumbent dump that already exists on disk.

Rendering is numpy/matplotlib only (dev Mac, NO Warp/CUDA); it reuses :func:`ff_viewer_html.build_viewer`.
Colour integrity (project rule): perceptually-ordered **turbo**; |F| on a **log₁₀** axis (>4 decades: resting
floor ~0.02 pN vs excluded-volume hotspots 10²–10³ pN), smooth loads linear; every ramp is annotated with units
and no axis is truncated. Full native resolution — NO downsampling.
"""
from __future__ import annotations

import json
from dataclasses import dataclass, field

import matplotlib
import numpy as np

#: Evidence rungs a figure may NOT print as an achieved rung, because `STATE.md` (c) blocks them.
#: `CONNECTED` is (c) 4: the wiring is complete but nothing SOLVES, so no magnitude from a composed
#: candidate is physical — and a dump that self-stamps it (as `gate_b_dynamic_v2.npz` does, which
#: `STATE.md` (b) records as "wrong") hands that word to every frame rendered from it.
BLOCKED_RUNGS: frozenset[str] = frozenset({"CONNECTED"})

_TURBO = matplotlib.colormaps["turbo"]

# |F| display range [pN] (log₁₀): resting floor → excluded-volume hotspot ceiling. A display range, annotated —
# NOT tuned to an outcome (the real per-component min/max are always printed in the scene name + caption).
F_VMIN, F_VMAX = 1.0e-2, 2.0e3
HOTSPOT_THR = 1.0                      # |F| [pN] above which a node/segment is drawn as a force hotspot

# Categorical colours held CLEAR of the turbo force ramp (blue→red), so a structural actor or a connector is
# never misread as a force hotspot. Components default to a per-key colour; connectors colour by FAMILY.
COMPONENT_COLORS = {
    "cortex": "#7f8fb0", "sf": "#ff8c42", "stress_fiber": "#ff8c42", "arc": "#ffb454",
    "lamellipodium": "#37d67a", "filopodium": "#00e0c7", "microtubule": "#ffd166", "mt": "#ffd166",
    "if": "#c98cff", "intermediate_filament": "#c98cff", "keratin": "#c98cff", "vimentin": "#a06bff",
    "nucleus": "#ffb454", "membrane": "#7fd4ff", "ecm": "#8d6e63", "collagen": "#8d6e63",
    "nmii": "#eaeaea", "myosin": "#eaeaea",
}
CONNECTOR_COLORS = {           # by connector FAMILY (force PATH between two components)
    "ERM": "#ff5db1", "FA-clutch": "#ffd60a", "FA": "#ffd60a", "LINC": "#c98cff", "MOTOR": "#ff35d6",
    "plectin": "#4dd0e1", "spectraplakin": "#7cffb2", "transient-actin": "#ff9e6d",
    "immersed-transfer": "#9db4ff",
}
_FALLBACK_COMPONENT = "#9aa7bf"
_FALLBACK_CONNECTOR = "#e0e0e0"


def component_color(key: str) -> str:
    return COMPONENT_COLORS.get(key.lower(), _FALLBACK_COMPONENT)


def connector_color(family: str) -> str:
    return CONNECTOR_COLORS.get(family, CONNECTOR_COLORS.get(family.upper(), _FALLBACK_CONNECTOR))


# ───────────────────────────────────────── colour ramps ──────────────────────────────────────────
def log_rgb(fmag: np.ndarray, vmin: float = F_VMIN, vmax: float = F_VMAX) -> np.ndarray:
    """|F| [pN] → (N,3) uint8 turbo on a clamped log₁₀ axis."""
    f = np.clip(np.asarray(fmag, np.float64), vmin, vmax)
    t = (np.log10(f) - np.log10(vmin)) / (np.log10(vmax) - np.log10(vmin))
    return (_TURBO(t)[:, :3] * 255.0 + 0.5).astype(np.uint8)


def lin_rgb(val: np.ndarray, vmin: float, vmax: float) -> np.ndarray:
    """Scalar load → (N,3) uint8 turbo on a clamped LINEAR axis."""
    t = (np.clip(np.asarray(val, np.float64), vmin, vmax) - vmin) / (vmax - vmin + 1e-30)
    return (_TURBO(t)[:, :3] * 255.0 + 0.5).astype(np.uint8)


def cbar(label: str, vmin: float, vmax: float, *, unit: str = "pN", log: bool = True, n: int = 7) -> dict:
    """Turbo colour-bar spec for build_viewer (gradient bottom→top; ticks span [vmin,vmax])."""
    grad = [matplotlib.colors.to_hex(_TURBO(i / (n - 1))) for i in range(n)]
    lo, hi = (np.log10(vmin), np.log10(vmax)) if log else (vmin, vmax)
    return {"label": label, "unit": (("log₁₀ " + unit) if log else unit), "grad": grad,
            "lo": float(lo), "hi": float(hi)}


def seg_pairs_from_offsets(offsets: np.ndarray, n_nodes: int) -> np.ndarray:
    """Internal per-fiber consecutive-node segment pairs → (Nseg,2) (last node of each fiber skipped)."""
    is_last = np.zeros(n_nodes, bool)
    is_last[np.asarray(offsets, np.int64)[1:] - 1] = True
    starts = np.arange(n_nodes)[~is_last]
    return np.stack([starts, starts + 1], axis=1)


# ───────────────────────────────────── normalized data model ─────────────────────────────────────
@dataclass
class Component:
    """One state-owning component. Owns a DISJOINT node population (global-unique ``actor_id``)."""
    key: str
    label: str
    kind: str                                  # "filaments" | "mesh" | "points"
    pos: np.ndarray                            # (N,3) float32 world-space nodes
    color: str
    seg: np.ndarray | None = None              # (M,2) int segment node-index pairs (filaments)
    faces: np.ndarray | None = None            # (F,3) int LOCAL triangle indices (mesh)
    fmag: np.ndarray | None = None             # (N,) per-node |F| [pN]
    actor_id: np.ndarray | None = None         # (N,) global-unique filament/actor id (prevents double-draw)
    order: int = 0                             # binding order for the cumulative view
    meta: dict = field(default_factory=dict)

    @property
    def n_nodes(self) -> int:
        return int(self.pos.shape[0])

    @property
    def n_actors(self) -> int:
        # a component may DECLARE its aggregate actor count (minifilaments, 1 surface) even when it namespaces
        # a per-node actor id for the disjointness invariant; the declared count wins for honest reporting.
        if "n_actors" in self.meta:
            return int(self.meta["n_actors"])
        if self.actor_id is not None and self.actor_id.size:
            return int(np.unique(self.actor_id).size)
        return 0


@dataclass
class Connector:
    """One connector family = the explicit force PATH between two components (co-location ≠ connection)."""
    key: str
    family: str
    label: str
    a_component: str
    b_component: str
    endpoints: np.ndarray                       # (2M,3) world-space [a0,b0,a1,b1,…] — the drawn joints
    color: str
    load: np.ndarray | None = None              # (M,) per-connector load/|F| [pN]
    meta: dict = field(default_factory=dict)

    @property
    def n(self) -> int:
        return int(self.endpoints.shape[0] // 2)


@dataclass
class CellDump:
    components: list[Component]
    connectors: list[Connector]
    meta: dict                                  # report / ledger / stage / R_cell / centroid / title bits

    def component(self, key: str) -> Component | None:
        return next((c for c in self.components if c.key == key), None)


# ──────────────────────────────────────── schema loaders ─────────────────────────────────────────
def load_dump(npz_path) -> CellDump:
    """Load a state dump (v2 generic or legacy) → the normalized :class:`CellDump`."""
    d = np.load(str(npz_path), allow_pickle=True)
    if "manifest_json" in d:
        return _normalize_v2(d)
    return normalize_legacy(d)


def _get(d, key, default=None):
    return d[key] if key in d.files else default


def _normalize_v2(d) -> CellDump:
    """Generic v2 dump: iterate ``manifest_json``; each component/connector names its own flat arrays."""
    man = json.loads(str(d["manifest_json"]))
    comps: list[Component] = []
    for i, c in enumerate(man.get("components", [])):
        pos = np.asarray(d[c["pos_key"]], np.float32).reshape(-1, 3)
        seg = np.asarray(d[c["seg_key"]], np.int64).reshape(-1, 2) if c.get("seg_key") else None
        if seg is None and c.get("offsets_key"):
            seg = seg_pairs_from_offsets(d[c["offsets_key"]], pos.shape[0])
        faces = np.asarray(d[c["faces_key"]], np.int64).reshape(-1, 3) if c.get("faces_key") else None
        fmag = np.asarray(d[c["fmag_key"]], np.float64).reshape(-1) if c.get("fmag_key") else None
        aid = np.asarray(d[c["actor_id_key"]], np.int64).reshape(-1) if c.get("actor_id_key") else None
        meta = dict(c.get("meta", {}))
        for field_name, arr_key in (c.get("array_meta") or {}).items():   # extra per-component arrays
            if arr_key in d.files:                                        # e.g. NMII head_verts, cortex f_steric
                meta[field_name] = np.asarray(d[arr_key])
        comps.append(Component(
            key=c["key"], label=c.get("label", c["key"]), kind=c.get("kind", "filaments"),
            pos=pos, color=c.get("color") or component_color(c["key"]), seg=seg, faces=faces,
            fmag=fmag, actor_id=aid, order=int(c.get("order", i)), meta=meta))
    cons: list[Connector] = []
    for c in man.get("connectors", []):
        ep = np.asarray(d[c["endpoints_key"]], np.float32).reshape(-1, 3) if c.get("endpoints_key") else \
            np.zeros((0, 3), np.float32)
        load = np.asarray(d[c["load_key"]], np.float64).reshape(-1) if c.get("load_key") else None
        cons.append(Connector(
            key=c["key"], family=c.get("family", c["key"]), label=c.get("label", c["key"]),
            a_component=c["a_component"], b_component=c["b_component"], endpoints=ep,
            color=c.get("color") or connector_color(c.get("family", c["key"])), load=load,
            meta=c.get("meta", {})))
    meta = man.get("meta", {})
    meta.setdefault("report", man.get("report", {}))
    meta.setdefault("ledger", man.get("ledger", {}))
    meta.setdefault("stage", man.get("stage", {}))
    _finalize_meta(meta, comps)
    return CellDump(comps, cons, meta)


def normalize_legacy(d) -> CellDump:
    """Re-express the incumbent ``assembled_state.npz`` (ac/cell/dump_state) as the generic model."""
    n_actin = int(d["n_actin"])
    pos = np.asarray(d["pos_post"], np.float32).reshape(-1, 3)          # combined [actin | myosin | meshes]
    pos_actin = pos[:n_actin]
    offsets = np.asarray(d["fiber_offsets"], np.int64)
    report = json.loads(str(d["report_json"]))
    ledger = json.loads(str(d["ledger_json"]))
    comps: list[Component] = []
    cons: list[Connector] = []

    # cortex F-actin (owns actin ids 0..n_actin-1 via fiber membership)
    fib_of_node = np.zeros(n_actin, np.int64)
    fib_of_node[offsets[:-1]] = 1
    fib_of_node = np.cumsum(fib_of_node) - 1
    comps.append(Component(
        key="cortex", label="cortex F-actin", kind="filaments", pos=pos_actin,
        color=component_color("cortex"), seg=seg_pairs_from_offsets(offsets, n_actin),
        fmag=np.asarray(d["f_total_post"], np.float64), actor_id=fib_of_node, order=0,
        meta={"n_actors": int(report["n_fibers"]),
              "f_steric": np.asarray(d["f_steric_post"], np.float64)}))

    # NMII minifilaments (backbones + heads) — a structural component drawn as line pairs
    mf_bb = np.asarray(d["mf_backbone_bonds"], np.int64).reshape(-1, 2)
    mf_hb = np.asarray(d["mf_head_bonds"], np.int64).reshape(-1, 2)
    if mf_bb.size:
        bb_v = pos[mf_bb].reshape(-1, 3)
        comps.append(Component(key="nmii", label="NMII minifilaments (Stam-Hocky)", kind="lines_raw",
                               pos=bb_v, color=component_color("nmii"), order=1,
                               meta={"n_actors": int(ledger.get("myosin_n_minifilaments", 0)),
                                     "head_verts": pos[mf_hb].reshape(-1, 3) if mf_hb.size else
                                     np.zeros((0, 3), np.float32),
                                     "n_heads": int(ledger.get("myosin_n_heads", 0))}))

    # deformable nucleus + plasma membrane (surface meshes)
    for key, label, order in (("nucleus", "deformable nucleus (Helfrich+lamina+ν→½)", 2),
                              ("membrane", "plasma membrane (Helfrich+γ_mem)", 3)):
        off = int(d[f"{key}_off"]) if f"{key}_off" in d.files else -1
        if off >= 0 and int(d[f"{key}_nverts"]) > 0:
            nv = int(d[f"{key}_nverts"])
            comps.append(Component(key=key, label=label, kind="mesh", pos=pos[off:off + nv],
                                   color=component_color(key), faces=np.asarray(d[f"{key}_faces"], np.int64),
                                   order=order, meta={"n_actors": 1}))

    # LINC connector: nucleus ↔ cortex (global node index pairs → world endpoints)
    linc = np.asarray(d["linc"], np.int64).reshape(-1, 2) if "linc" in d.files and d["linc"].size else \
        np.zeros((0, 2), np.int64)
    if linc.size:
        cons.append(Connector(key="linc", family="LINC", label="LINC nucleus↔cortex", a_component="nucleus",
                              b_component="cortex", endpoints=pos[linc].reshape(-1, 3),
                              color=connector_color("LINC")))

    meta = {"report": report, "ledger": ledger, "stage": {"name": "incumbent build_cell", "schema": "legacy"}}
    _finalize_meta(meta, comps)
    return CellDump(comps, cons, meta)


def _finalize_meta(meta: dict, comps: list[Component]) -> None:
    all_pos = np.concatenate([c.pos for c in comps if c.pos.size], 0) if comps else np.zeros((1, 3))
    meta["centroid"] = all_pos.mean(0).astype(float).tolist()
    meta.setdefault("R_cell", 7.5)


# ─────────────────────────────────────── scene construction ──────────────────────────────────────
def _component_layers(c: Component, *, color_by_force: bool, opacity: float, size: float,
                      on_top: bool = False, clip: bool = True) -> list[dict]:
    """Render one component as viewer layers (kind-aware).

    Every layer carries ``group="compartment"`` and a ``swatch`` = the component's categorical colour so the
    viewer's per-layer checkbox panel (PI 2026-07-23) groups it correctly and shows a meaningful colour chip
    even when the drawn colour is white (force-coloured layers paint per-node |F|, not the flat colour).
    """
    if c.kind == "mesh" and c.faces is not None:
        return [{"name": c.label, "kind": "mesh", "verts": c.pos, "faces": c.faces,
                 "color": c.color, "swatch": c.color, "group": "compartment", "opacity": opacity, "clip": clip}]
    if c.kind == "lines_raw":                                   # pre-paired line verts (e.g. NMII bonds)
        out = [{"name": c.label, "kind": "lines", "verts": c.pos, "color": c.color, "swatch": c.color,
                "group": "compartment", "size": size, "opacity": opacity, "on_top": on_top}]
        hv = c.meta.get("head_verts")
        if hv is not None and len(hv):
            out.append({"name": c.label + " · motor heads", "kind": "lines", "verts": hv,
                        "color": "#ff35d6", "swatch": "#ff35d6", "group": "compartment",
                        "size": max(size - 1.0, 1.0), "opacity": opacity, "on_top": on_top})
        return out
    if c.seg is None:                                           # bare point cloud
        return [{"name": c.label, "kind": "points", "verts": c.pos, "color": c.color, "swatch": c.color,
                 "group": "compartment", "size": size, "opacity": opacity, "clip": clip}]
    verts = c.pos[c.seg].reshape(-1, 3)
    layer = {"name": c.label, "kind": "lines", "verts": verts, "color": c.color, "swatch": c.color,
             "group": "compartment", "size": size, "opacity": opacity, "clip": clip, "on_top": on_top}
    if color_by_force and c.fmag is not None:
        layer["color"] = "#ffffff"
        layer["color_frames"] = [log_rgb(c.fmag)[c.seg].reshape(-1, 3)]
    return [layer]


def _connector_layer(cn: Connector, *, color_by_load: bool, size: float = 2.0, opacity: float = 0.95) -> dict:
    layer = {"name": f"{cn.label} · {cn.family} ({cn.n:,} joints)", "kind": "lines", "verts": cn.endpoints,
             "color": cn.color, "swatch": cn.color, "group": "connector", "size": size, "opacity": opacity,
             "on_top": True}
    if color_by_load and cn.load is not None and cn.load.size:
        lo, hi = float(np.nanmin(cn.load)), float(max(np.nanmax(cn.load), 1e-6))
        layer["color"] = "#ffffff"
        layer["color_frames"] = [np.repeat(lin_rgb(cn.load, lo, hi), 2, axis=0)]
    return layer


def ref_sphere_layer(dump: CellDump, opacity: float = 0.04):
    """Faint R=R_cell scale-guide sphere (imported lazily to avoid a hard morphology dep at import)."""
    from aleph.scripts.ff_cell_morphology import uv_sphere
    c = np.asarray(dump.meta["centroid"], float)
    v, f = uv_sphere(c, float(dump.meta.get("R_cell", 7.5)), nu=48, nv=28)
    return {"name": f"R={dump.meta.get('R_cell', 7.5):g}µm reference sphere", "kind": "mesh",
            "verts": v, "faces": f, "color": "#5a6b8c", "swatch": "#5a6b8c", "group": "reference",
            "opacity": opacity, "clip": True}


def _gate_b_dynamic_tension_scene(dump: CellDump, faint_ref: dict) -> tuple[str, list[dict]] | None:
    """GATE-B milestone scene: the tensed cortex coloured by |F| with the BOUND NMII motors overlaid.

    Data-driven + guarded: built ONLY when the dump carries the dynamic cortex-motor signals — a ``cortex``
    component with a per-node |F| field and an ``nmii`` component carrying per-head ``head_bound`` / ``head_pos``
    (written by ``ac_gate_b_cortex_motor_native.py --dump-viz``).  For the incumbent assembled dump (no bound
    head map) it returns ``None`` and the scene simply does not appear.

    The scene makes "the cortical tension EMERGED from exactly these bound motors" legible: the cortex is
    coloured by the emergent per-node |F| (same log-turbo ramp/units as every other force scene), the BOUND
    heads are drawn as distinct bright markers, and the head→actin crossbridge force PATH is overlaid — so the
    reader can see the elevated-|F| cortex regions sitting under the bound crossbridges.
    """
    cortex = dump.component("cortex")
    nmii = dump.component("nmii")
    if cortex is None or cortex.fmag is None or cortex.seg is None:
        return None
    if nmii is None or "head_bound" not in nmii.meta or "head_pos" not in nmii.meta:
        return None

    stage = dump.meta.get("stage", {})
    bpct, gtot = stage.get("bound_pct"), stage.get("gamma_total_pn_per_um")
    # bound-% / γ go in the LAYER names (and the always-on global title), NOT the scene KEY — the key stays a
    # stable, scriptable identifier so `browser_check --scene "◆ GATE-B …"` selects it regardless of run values.
    gtag = (f" · γ={float(gtot):.3g} pN/µm" if gtot is not None else "")

    # cortex coloured by the emergent per-node |F| (the standard log-turbo force ramp; force_cbar registered).
    layers: list[dict] = _component_layers(cortex, color_by_force=True, opacity=0.7, size=1.3, clip=True)
    # faint minifilament backbones for structural context (NOT the head arms — the BOUND heads carry the story).
    if nmii.pos.size:
        layers.append({"name": f"NMII minifilament backbones (faint · {nmii.n_actors:,})", "kind": "lines",
                       "verts": nmii.pos, "color": nmii.color, "swatch": nmii.color, "group": "compartment",
                       "size": 1.0, "opacity": 0.22})
    # the BOUND heads as distinct bright markers — the motors that produced the tension.
    hp = np.asarray(nmii.meta["head_pos"], np.float32).reshape(-1, 3)
    hb = np.asarray(nmii.meta["head_bound"]).reshape(-1).astype(bool)
    n_bound = int(hb.sum())
    if hp.shape[0] == hb.shape[0] and n_bound:
        bpct_s = f"{float(bpct):.0f}%" if bpct is not None else f"{100.0 * n_bound / max(hb.size, 1):.0f}%"
        # SMALL depth-tested dots via the per-layer pt_size override: at the native bound count (thousands of
        # heads, ~10 µm⁻² on the shell) a large marker overdraws into a solid magenta COVER that HIDES the
        # emergent-|F| cortex — the whole point of the scene.  pt_size≈0.3 µm (~1/8 the global PT_SIZE=2.5) keeps
        # the γ field the STAR and makes the bound heads a legible magenta STIPPLE over it (verified on a
        # native-scale ~5k-head render: at 0.83 the field is covered, at ~0.3 it reads through).  Tunable — one
        # number: raise toward 0.5 if a sparser dump wants punchier markers, lower if a denser one still covers.
        # Depth-tested so the front shell occludes the far hemisphere (points ignore on_top).
        layers.append({"name": f"● BOUND NMII heads ({n_bound:,} · {bpct_s}) — the motors that made γ{gtag}",
                       "kind": "points", "verts": hp[hb], "color": "#ff35d6", "swatch": "#ff35d6",
                       "group": "connector", "pt_size": 0.3, "opacity": 1.0})
    # the explicit head→actin crossbridge force PATH (MOTOR connector), in the MOTOR colour.
    for cn in dump.connectors:
        if cn.family.upper() == "MOTOR" and {cn.a_component, cn.b_component} == {"nmii", "cortex"}:
            layers.append(_connector_layer(cn, color_by_load=False, size=2.6, opacity=1.0))
    layers.append(faint_ref)
    return "◆ GATE-B dynamic tension · cortex |F| + bound NMII heads", layers


def build_scenes(dump: CellDump) -> tuple[dict, dict]:
    """Assemble the full data-driven scene set + colour bars from a normalized dump.

    Produces, with NO per-component hardcoding:
      * one ISOLATION scene per component (each alone + a faint scale ref; force-coloured if |F| present),
      * one CONNECTOR load-path scene per connector family (the explicit joints between two components),
      * a CUMULATIVE composition family that grows one component at a time (+ the connectors that become
        drawable as both endpoints land), coloured by per-component |F|,
      * a COMPOSED whole-cell scene and a steric-hotspot scene (if a steric field is present).
    """
    scenes: dict[str, list] = {}
    cbars: dict[str, dict] = {}
    ref = ref_sphere_layer(dump)
    faint_ref = dict(ref, opacity=0.03)
    force_cbar = cbar("per-node |F|", F_VMIN, F_VMAX)

    comps = sorted(dump.components, key=lambda c: c.order)
    ctx = None
    cortex = dump.component("cortex")
    if cortex is not None and cortex.seg is not None:
        ctx = {"name": f"cortex context (faint · {cortex.n_actors:,} fil)", "kind": "lines",
               "verts": cortex.pos[cortex.seg].reshape(-1, 3), "color": "#39435a", "swatch": "#39435a",
               "group": "reference", "size": 0.7, "opacity": 0.06, "clip": True}

    # ---- COMPOSED whole cell (all components, cortex force-coloured) --------------------------------------
    whole: list[dict] = []
    for c in comps:
        whole += _component_layers(c, color_by_force=(c.key == "cortex"),
                                   opacity=0.45 if c.kind == "filaments" else
                                   (0.55 if c.key == "nucleus" else 0.08 if c.key == "membrane" else 0.9),
                                   size=1.15)
    for cn in dump.connectors:
        whole.append(_connector_layer(cn, color_by_load=False, opacity=0.5, size=1.4))
    whole.append(ref)
    k_whole = "composed cell · cortex |F| (log pN) · CUT slider"
    scenes[k_whole] = whole
    cbars[k_whole] = force_cbar

    # ---- GATE-B dynamic-tension scene (guarded: only for a dynamic cortex-motor dump) --------------------
    gb = _gate_b_dynamic_tension_scene(dump, faint_ref)
    if gb is not None:
        scenes[gb[0]] = gb[1]
        cbars[gb[0]] = force_cbar

    # ---- per-component ISOLATION scenes (data-driven: one per component in the dump) ---------------------
    for c in comps:
        has_f = c.fmag is not None
        layers = _component_layers(c, color_by_force=has_f, opacity=0.62 if c.kind != "mesh" else 0.72,
                                   size=1.25, on_top=(c.kind == "lines_raw"))
        if c.key != "cortex" and ctx is not None:
            layers = [ctx] + layers
        layers.append(faint_ref)
        key = f"◦ {c.key} only · {c.label}"
        scenes[key] = layers
        if has_f:
            cbars[key] = force_cbar

    # ---- per-connector LOAD-PATH scenes (the force PATH; co-location ≠ connection) -----------------------
    for cn in dump.connectors:
        a, b = dump.component(cn.a_component), dump.component(cn.b_component)
        layers = []
        for comp, op in ((a, 0.20), (b, 0.20)):
            if comp is not None:
                layers += _component_layers(comp, color_by_force=False, opacity=op, size=0.9, clip=True)
        color_load = cn.load is not None and cn.load.size > 0
        layers.append(_connector_layer(cn, color_by_load=color_load, size=2.4, opacity=1.0))
        layers.append(faint_ref)
        key = f"⇄ {cn.family} · {cn.a_component}↔{cn.b_component}"
        scenes[key] = layers
        if color_load:
            lo, hi = float(np.nanmin(cn.load)), float(max(np.nanmax(cn.load), 1e-6))
            cbars[key] = cbar(f"{cn.family} load", max(lo, 1e-6), hi, log=False)

    # ---- CUMULATIVE composition (grows one component at a time; connectors join when both ends land) -----
    present: set[str] = set()
    for i, c in enumerate(comps, 1):
        present.add(c.key)
        layers: list[dict] = []
        for c2 in comps:
            if c2.key not in present:
                continue
            newest = c2.key == c.key
            layers += _component_layers(
                c2, color_by_force=(c2.fmag is not None),
                opacity=(0.7 if newest else 0.28) if c2.kind != "mesh" else (0.6 if newest else 0.18),
                size=1.4 if newest else 1.0, on_top=(c2.kind == "lines_raw" and newest))
        for cn in dump.connectors:
            if cn.a_component in present and cn.b_component in present:
                layers.append(_connector_layer(cn, color_by_load=False, size=1.8, opacity=0.85))
        layers.append(faint_ref)
        key = f"▸ compose {i}/{len(comps)} · +{c.key}"
        scenes[key] = layers
        cbars[key] = force_cbar

    # ---- STERIC interpenetration hotspots (legacy cortex carries a steric-only field) -------------------
    if cortex is not None and "f_steric" in cortex.meta:
        fs = np.asarray(cortex.meta["f_steric"], np.float64)
        seg = cortex.seg
        hot = (fs[seg[:, 0]] > HOTSPOT_THR) | (fs[seg[:, 1]] > HOTSPOT_THR)
        n_interp = int(np.count_nonzero(fs > 0.0))
        s_rgb = log_rgb(fs, 1.0, F_VMAX)
        key = "steric interpenetration hotspots (missing 7 nm EV)"
        scenes[key] = [
            dict(ctx, opacity=0.10) if ctx else faint_ref,
            {"name": f"steric hotspots ({n_interp:,} nodes w/o 7nm EV · {int(hot.sum()):,} seg >{HOTSPOT_THR:g}pN)",
             "kind": "lines", "verts": cortex.pos[seg[hot]].reshape(-1, 3), "color": "#ffffff",
             "swatch": "#ff5db1", "group": "compartment", "size": 2.2,
             "opacity": 1.0, "on_top": True, "color_frames": [s_rgb[seg[hot]].reshape(-1, 3)]}]
        cbars[key] = cbar("steric-only |F|", 1.0, F_VMAX)

    return scenes, cbars


def compose_title(dump: CellDump) -> str:
    """One-line provenance/stat title for the viewer (units, per-component counts, stage, |F| range)."""
    rep = dump.meta.get("report", {})
    stage = dump.meta.get("stage", {})
    parts = [f"{c.key}:{c.n_actors:,}a/{c.n_nodes:,}n" for c in sorted(dump.components, key=lambda c: c.order)]
    conparts = [f"{cn.family}:{cn.n:,}" for cn in dump.connectors]
    cortex = dump.component("cortex")
    frange = ""
    if cortex is not None and cortex.fmag is not None:
        fm = cortex.fmag[np.isfinite(cortex.fmag)]
        if fm.size:
            frange = (f" · cortex |F|∈[{fm.min():.1e},{fm.max():.0f}] pN med "
                      f"{np.percentile(fm, 50):.3f} (log turbo)")
    txn = ""
    if "outer_accepted" in rep:
        txn = (" · ACCEPTED" if rep.get("outer_accepted") else
               " · REJECTED+ROLLED BACK" if rep.get("outer_rolled_back") else " · REJECTED")
    # ── the title may not print what the record itself says is not quotable ──────────────────────────
    # A figure is the strongest re-quotation vector this project has: a screenshot outlives the record
    # it came from, and STATE.md (c) 10 already lists two retired numbers "still live on the public
    # gh-pages gallery". The dump CARRIES `quantitative_claim_status`; until 2026-07-29 this function
    # ignored it and stamped both a (c) 4-blocked rung and a (c) 17-blocked γ into every frame.
    rung = str(stage.get("evidence") or "")
    if rung in BLOCKED_RUNGS:
        stg = f"{stage.get('name', 'ac/engine')} [rung self-stamped {rung} — NOT QUOTABLE, STATE.md (c) 4]"
    else:
        stg = f"{stage.get('name', 'ac/engine')}" + (f" [{rung}]" if rung else "")
    gb = ""
    blocked = str(stage.get("quantitative_claim_status") or "").upper() == "BLOCKED"
    if stage.get("bound_pct") is not None:
        # the bound FRACTION is a count, not a magnitude, and stays visible either way
        gb = f" · GATE-B bound {float(stage['bound_pct']):.0f}%"
        if stage.get("gamma_total_pn_per_um") is not None:
            gb += (" · γ MAGNITUDE WITHHELD (QuantitativeClaim: BLOCKED)" if blocked else
                   f" · γ_total {float(stage['gamma_total_pn_per_um']):.3g} pN/µm (source "
                   f"{float(stage.get('gamma_source_pn_per_um', 0.0)):.3g} / network "
                   f"{float(stage.get('gamma_network_pn_per_um', 0.0)):.3g})")
    return (f"AC engine staged build — {stg}{txn}{gb} · components {{{' '.join(parts)}}} · "
            f"connectors {{{' '.join(conparts) or '—'}}}{frange} · units µm · FULL-RES (no downsample)")
