"""Render the REAL composed ``ac/engine`` world (14 components / 36 connectors) in the per-compartment viewer.

This is the bridge the staged-build viewer was missing: it drives the ACTUAL composed-world path —
:func:`~aleph.engine.build_composed_cell_world` + :func:`~aleph.engine.dump_composed_world_state` —
to obtain the authoritative per-compartment CENSUS (node / face / endpoint counts, the disjoint global
actor / filament-ID namespace, and the composite α2β1–collagen FA-joint de-duplication), and then feeds that
census into the generic v2 dump schema (:func:`aleph.scripts.dump_state.write_v2_dump`) that
:mod:`aleph.scripts.ac_viz_common` renders.  The moment the GPU lane produces a real *device* census (with
real geometry) via ``dump_engine_actor``, the same bridge structure renders it with no viewer edit.

What is REAL here vs SYNTHETIC (honest, per the project rules):

* REAL — the composition itself: the 14-component / 36-connector topology, roles, connector families and
  endpoint components come from :func:`~aleph.engine.reference_cell_architecture` driven through the real
  ``build_composed_cell_world``; the per-compartment counts, the **global-unique disjoint filament-ID blocks**,
  and the composite-FA-joint-is-ONE-actor de-dup come from the real ``dump_composed_world_state`` census (its
  ``actor_ids_are_unique`` / ``filament_id_blocks_are_disjoint`` invariants are asserted and recorded).  The
  cortex census is the KB-sourced native count (70,686 filaments / 494,802 nodes; CLAUDE.md first baseline).
* SYNTHETIC — the node GEOMETRY and the |F| / load *fields* are placeholders (the real device dump comes from
  the A5000 lane later).  Non-cortex census counts are small synthetic render placeholders, NOT physiological
  claims (MT / IF / filopodium counts are PI-GAPs and are deliberately NOT asserted as native).  Cortex is
  drawn as a labelled subsample of its 70,686-filament block.

Every component becomes its own isolation layer and every connector its own explicit joint line-set
(co-location ≠ connection: the joints are drawn between the two endpoint clouds, not as mere proximity), all
independently toggleable in the viewer's grouped compartment / connector / reference checkbox panel.

    python -m aleph.scripts.ac_composed_world_dump \
        --out aleph/outputs/ac/cell_assembled/composed_world_v2.npz --html
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from aleph.engine import (
    ComposedCellRuntimes,
    FacadeDispatch,
    build_composed_cell_world,
    dump_composed_world_state,
    reference_cell_architecture,
)
from aleph.engine.contracts import (
    CellArchitecture,
    EvidenceLabel,
    EvidenceRung,
    QuantitativeClaim,
)
from aleph.engine.dispatch import canonical_facade_claims
from aleph.engine.ecm_world import ECMWorld
from aleph.engine.fluid_core import FluidCore
from aleph.engine.intermediate_filament_rig import IntermediateFilamentRig
from aleph.engine.microtubule_rig import MicrotubuleRig
from aleph.engine.nmii_actuator import NMIIActuator
from aleph.engine.protrusion import ProtrusionActors
from aleph.engine.stress_fiber import StressFiberActor
from aleph.engine.surface_body import SurfaceBody
from aleph.scripts.ac_synth_dump import _radial_filaments
from aleph.scripts.dump_state import _c, _k, write_v2_dump
from aleph.scripts.ff_cell_morphology import uv_sphere

R = 7.5          # µm cell radius (synthetic render placeholder geometry, consistent with ac_synth_dump)
R_NUC = 5.0      # µm nucleus radius (synthetic render placeholder; real value is a PI-GAP, not asserted)
RNG = np.random.default_rng(23)

# ── composed-world census (the counts dump_composed_world_state reads from each bound runtime) ─────────────
# cortex is the KB-sourced native first baseline (CLAUDE.md: 70,686 cortical F-actin filaments / 494,802 actin
# nodes). Every OTHER count is a small SYNTHETIC render placeholder — NOT a physiological claim (MT / IF /
# filopodium population sizes are PI-GAPs and are deliberately not asserted here).
_CENSUS: dict[str, dict[str, int]] = {
    "cortex": {"n_filaments": 70686, "n_nodes": 494802},          # KB-sourced native baseline
    "sf_arc": {"n_filaments": 60, "n_nodes": 120},                # synthetic placeholders below
    "microtubule": {"n_filaments": 120, "n_nodes": 240},
    "intermediate_filament": {"n_filaments": 150, "n_nodes": 300},
    "lamellipodium": {"n_filaments": 200, "n_nodes": 400},
    "filopodium": {"n_filaments": 8, "n_nodes": 16},
    "nmii": {"n_filaments": 100, "n_nodes": 200, "n_endpoints": 400},
    "ecm": {"n_filaments": 300, "n_nodes": 600},
    "cytosol": {"n_nodes": 343},
    "focal_adhesion": {"n_endpoints": 80},
}
_CORTEX_DRAW = 2000                                               # labelled subsample of the 70,686 block


# ─────────────────────────────── composed-world runtime doubles (CPU-safe) ────────────────────────────────
@dataclass(slots=True)
class _CensusRuntime:
    """A component/connector state-owner double carrying only its census + the full transaction/ledger API.

    The dev Mac has no CUDA, so the real Warp-resident owners cannot be built here; these doubles expose the
    host-side count attributes ``dump_composed_world_state`` reads and the transaction/mechanics/ledger hooks
    the composed-world builder validates.  They own no geometry — geometry is synthesised from the census.
    """

    name: str | None = None
    component_a: str | None = None
    component_b: str | None = None
    mechanical_group: str | None = None
    n_filaments: int = 0
    n_nodes: int = 0
    n_faces: int = 0
    n_endpoints: int = 0

    def accumulate(self, *args) -> None: ...
    def snapshot_candidate(self) -> None: ...
    def rollback(self, accepted) -> None: ...
    def commit_irreversible(self, accepted, dt_phys, rng_seed) -> None: ...
    def accumulate_ledger(self, ledger) -> None: ...

    # Required of every connector declaring `kinetics=True` (and every component declaring
    # `has_events=True`) by `CellActor.assert_fully_bound`. Its absence is why this script raised
    # `TypeError: connector 'membrane_erm_cortex' has an incomplete runtime API` and could not
    # regenerate the cumulative composition render — which made the per-stage visualization gate
    # unsatisfiable for EVERY track, not just for one component. The double proposes nothing; it
    # exists so the census world can be BUILT, and the engine's requirement is untouched.
    def propose_events(self, *args, **kwargs) -> None: ...


@dataclass(slots=True)
class _FacadeSpy:
    """Dispatch-facade double exposing every canonical orchestration method the manifest may claim."""

    calls: list = field(default_factory=list)

    def accumulate_mechanics(self, *a, **k) -> None: self.calls.append(("accumulate_mechanics",))
    def accumulate_candidate(self, *a, **k) -> None: self.calls.append(("accumulate_candidate",))
    def candidate_iteration(self, *a, **k) -> None: self.calls.append(("candidate_iteration",))


_FACADE_METHOD = {
    SurfaceBody: "accumulate_mechanics",
    FluidCore: "candidate_iteration",
    StressFiberActor: "accumulate_mechanics",
    MicrotubuleRig: "accumulate_mechanics",
    IntermediateFilamentRig: "accumulate_mechanics",
    ProtrusionActors: "accumulate_mechanics",
    NMIIActuator: "accumulate_candidate",
    ECMWorld: "accumulate_mechanics",
}
_FA_EDGES = ("fa_actin_anchor", "integrin_collagen_clutch")
_FA_GROUP = "alpha2beta1_collagen_series"


def build_census_world(architecture: CellArchitecture):
    """Drive the REAL composed-world builder over census-carrying doubles → a bound ``ComposedCellWorld``."""
    component_owners: dict[str, object] = {}
    for component in architecture.components:
        if not component.dynamically_evolving and not component.owns_geometry:
            component_owners[component.name] = object()          # world_boundary: inert reference frame
            continue
        component_owners[component.name] = _CensusRuntime(**_CENSUS.get(component.name, {}))

    fa_joint = _CensusRuntime(mechanical_group=_FA_GROUP, n_endpoints=_CENSUS["focal_adhesion"]["n_endpoints"])
    connector_runtimes: dict[str, object] = {}
    for connector in architecture.connectors:
        if connector.name in _FA_EDGES:
            connector_runtimes[connector.name] = fa_joint        # the composite joint is ONE object
            continue
        connector_runtimes[connector.name] = _CensusRuntime(
            name=connector.name, component_a=connector.component_a, component_b=connector.component_b,
        )

    facade_dispatch = {}
    for claim in canonical_facade_claims():
        method = _FACADE_METHOD[claim.owner_type]
        args = (object(), 0.05) if method == "candidate_iteration" else (object(),)
        facade_dispatch[claim.owner_type] = FacadeDispatch(_FacadeSpy(), args=args)

    runtimes = ComposedCellRuntimes(component_owners, connector_runtimes, facade_dispatch)
    return build_composed_cell_world(runtimes, architecture)


# ─────────────────────────────────── geometry archetypes (synthetic) ─────────────────────────────────────
def _fluid_grid(n_target: int) -> np.ndarray:
    """Interior fluid-field cell centres on a cubic lattice clipped to the cell interior (faint point cloud)."""
    m = max(int(round(n_target ** (1 / 3))), 2)
    lin = np.linspace(-R * 0.8, R * 0.8, m)
    gx, gy, gz = np.meshgrid(lin, lin, lin, indexing="ij")
    pts = np.stack([gx.ravel(), gy.ravel(), gz.ravel()], axis=1)
    return pts[np.linalg.norm(pts, axis=1) < R * 0.85].astype(np.float32)


def _fa_sites(n: int) -> np.ndarray:
    """n integrin-clutch sites clustered on the basal (−z) shell (focal adhesions are point-like clutches)."""
    d = RNG.normal(size=(n, 3))
    d[:, 2] = -np.abs(d[:, 2]) - 0.3                              # push onto the ventral/basal hemisphere
    d /= np.linalg.norm(d, axis=1, keepdims=True)
    return (d * (R - 0.1)).astype(np.float32)


def _component_specs(dump, architecture: CellArchitecture) -> tuple[list[dict], dict[str, np.ndarray]]:
    """Build one v2 component spec per census entry; return (specs, node_clouds) for connector resolution."""
    specs: list[dict] = []
    clouds: dict[str, np.ndarray] = {}
    order = 0

    def add(key, label, kind, pos, **kw):
        nonlocal order
        clouds[key] = np.asarray(pos, np.float32).reshape(-1, 3)
        specs.append(_c(key, label, kind, pos, order=order, **kw))
        order += 1

    def census_meta(name: str) -> dict:
        """Stamp the REAL composed-world census (counts + global-unique filament-ID block) into the layer."""
        try:
            e = dump.by_name(name)
        except KeyError:
            return {"census": "unbound"}
        return {"census_actor_id": e.actor_id, "census_filament_id_base": e.filament_id_base,
                "census_n_filaments": e.n_filaments, "census_n_nodes": e.n_nodes,
                "census_n_faces": e.n_faces, "census_n_endpoints": e.n_endpoints,
                "n_actors": int(e.n_filaments) or 1}

    # membrane — outer Helfrich sheet (mesh)
    mv, mf = uv_sphere(np.zeros(3), R, nu=40, nv=22)
    add("membrane", "plasma membrane (Helfrich) · SYNTH geom", "mesh", mv, faces=mf, meta=census_meta("membrane"))
    # cortex — KB 70,686-filament block, drawn as a labelled subsample near the shell
    pos, seg, fmag, aid = _radial_filaments(_CORTEX_DRAW, R - 0.5, R - 0.05, seed=1)
    cm = census_meta("cortex")
    cm.update({"n_actors": 70686, "drawn_subsample": _CORTEX_DRAW})
    add("cortex", f"cortex F-actin · 70,686 block ({_CORTEX_DRAW:,} drawn) · SYNTH geom", "filaments", pos,
        seg=seg, fmag=fmag, actor_id=aid, meta=cm)
    # cytosol — Biot/Darcy fluid field, drawn as faint interior grid-cell centres
    add("cytosol", "cytosol (Biot/Darcy field) · SYNTH grid", "points", _fluid_grid(_CENSUS["cytosol"]["n_nodes"]),
        meta=census_meta("cytosol"))
    # nucleus — inner deformable core (mesh)
    nv, nf = uv_sphere(np.zeros(3), R_NUC, nu=28, nv=16)
    add("nucleus", "deformable nucleus · SYNTH geom", "mesh", nv, faces=nf, meta=census_meta("nucleus"))
    # sf_arc — long ventral stress fibers crossing the cell
    a = RNG.normal(size=(_CENSUS["sf_arc"]["n_filaments"], 3)); a /= np.linalg.norm(a, axis=1, keepdims=True)
    a *= (R - 1.0); b = -a + RNG.normal(scale=0.5, size=a.shape)
    sf_pos = np.stack([a, b], axis=1).reshape(-1, 3); n_sf = a.shape[0]
    add("sf_arc", "stress fibers / arcs · SYNTH geom", "filaments", sf_pos,
        seg=np.arange(2 * n_sf).reshape(n_sf, 2), fmag=np.full(2 * n_sf, 50.0),
        actor_id=np.repeat(np.arange(n_sf), 2), meta=census_meta("sf_arc"))
    # focal_adhesion — point-like integrin clutches on the basal shell (owns no bulk geometry)
    add("focal_adhesion", "focal adhesions (α2β1 clutches) · SYNTH sites", "points",
        _fa_sites(_CENSUS["focal_adhesion"]["n_endpoints"]), meta=census_meta("fa_actin_anchor"))
    # microtubule — aster from the centrosome
    mtp, mts, mtf, mta = _radial_filaments(_CENSUS["microtubule"]["n_filaments"], 0.3, R - 0.3, jitter=0.05, seed=4)
    add("microtubule", "microtubule aster · SYNTH geom (count=PI-GAP)", "filaments", mtp, seg=mts, fmag=mtf,
        actor_id=mta, meta=census_meta("microtubule"))
    # intermediate_filament — perinuclear cage
    ifp, ifs, iff, ifa = _radial_filaments(_CENSUS["intermediate_filament"]["n_filaments"], R_NUC, 6.5, jitter=0.3, seed=5)
    add("intermediate_filament", "IF perinuclear cage · SYNTH geom (count=PI-GAP)", "filaments", ifp, seg=ifs,
        fmag=iff, actor_id=ifa, meta=census_meta("intermediate_filament"))
    # lamellipodium — flat protrusion sheet near the +x pole
    lp, _, _, _ = _radial_filaments(_CENSUS["lamellipodium"]["n_filaments"], R - 0.2, R + 1.2, jitter=0.4, seed=6)
    keep = lp.reshape(-1, 2, 3)[:, 1, 0] > R * 0.5
    lp2 = lp.reshape(-1, 2, 3)[keep].reshape(-1, 3); n_lp = lp2.shape[0] // 2
    add("lamellipodium", "lamellipodium · SYNTH geom", "filaments", lp2,
        seg=np.arange(2 * n_lp).reshape(n_lp, 2), fmag=np.full(2 * n_lp, 20.0),
        actor_id=np.repeat(np.arange(n_lp), 2), meta=census_meta("lamellipodium"))
    # filopodium — finger bundles beyond the shell at +y
    fp = []
    for _ in range(_CENSUS["filopodium"]["n_filaments"]):
        dv = np.array([RNG.uniform(-0.3, 0.3), 1.0, RNG.uniform(-0.3, 0.3)]); dv /= np.linalg.norm(dv)
        fp += [dv * (R - 0.3), dv * (R + 2.0)]
    fp = np.asarray(fp, np.float32); n_fp = fp.shape[0] // 2
    add("filopodium", "filopodia · SYNTH geom (count=PI-GAP)", "filaments", fp,
        seg=np.arange(2 * n_fp).reshape(n_fp, 2), fmag=np.full(2 * n_fp, 8.0),
        actor_id=np.repeat(np.arange(n_fp), 2), meta=census_meta("filopodium"))
    # nmii — bipolar minifilament rods (Stam-Hocky backbone), drawn as short crossing rods
    mp, ms, _, ma = _radial_filaments(_CENSUS["nmii"]["n_filaments"], R - 2.5, R - 0.5, jitter=0.5, seed=9)
    add("nmii", "NMII minifilaments (Stam-Hocky) · SYNTH geom", "filaments", mp, seg=ms,
        fmag=np.full(mp.shape[0], 30.0), actor_id=ma, meta=census_meta("nmii"))
    # ecm — collagen fibers OUTSIDE the cell
    ep, es, _, ea = _radial_filaments(_CENSUS["ecm"]["n_filaments"], R + 0.5, R + 4.0, jitter=1.0, seed=8)
    add("ecm", "collagen ECM · SYNTH geom", "filaments", ep, seg=es, fmag=np.full(ep.shape[0], 5.0),
        actor_id=ea, meta=census_meta("ecm"))
    # world_boundary — faint physiological far-field reference shell (owns no dynamic geometry)
    wv, _ = uv_sphere(np.zeros(3), R + 4.0, nu=24, nv=14)
    add("world_boundary", "world far-field reference frame · SYNTH shell", "points", wv,
        meta={"census": "reference-frame", "n_actors": 1})
    # extracellular_medium — the FREE face of the asymmetric world boundary (PI D5-A). Drawn as a shell
    # just outside the membrane so the viewer shows the medium exists and touches only the membrane.
    # DECLARED ONLY: there is no exterior solve behind it until T10, so the label says so rather than
    # letting a rendered layer read as a working compartment.
    mv, _ = uv_sphere(np.zeros(3), R + 0.25, nu=20, nv=12)
    add("extracellular_medium", "extracellular medium (free face) · DECLARED ONLY, no solve until T10",
        "points", mv, meta={"census": "declared-only", "n_actors": 1})
    return specs, clouds


def _nearest_pairs(pa: np.ndarray, pb: np.ndarray, n: int, seed: int) -> np.ndarray:
    """n explicit joints: n random A-nodes each connected to their nearest B-node → world-space endpoint pairs.

    This draws the ACTUAL joint (a→b line), not co-location: two clouds that merely overlap in space share no
    line unless a joint is explicitly emitted here.
    """
    if not pa.size or not pb.size:
        return np.zeros((0, 3), np.float32)
    rng = np.random.default_rng(seed)
    ia = rng.choice(pa.shape[0], size=min(n, pa.shape[0]), replace=False)
    out = []
    for i in ia:
        j = int(np.argmin(np.linalg.norm(pb - pa[i], axis=1)))
        out += [pa[i], pb[j]]
    return np.asarray(out, np.float32)


def _connector_specs(dump, architecture: CellArchitecture, clouds: dict[str, np.ndarray]) -> list[dict]:
    """One explicit joint line-set per REGISTERED connector edge, sized by the census endpoint count."""
    registered = set(dump_registered_connectors(dump, architecture))
    cons: list[dict] = []
    for cc in architecture.connectors:
        if cc.name not in registered:
            continue
        a, b = cc.component_a, cc.component_b
        pa, pb = clouds.get(a), clouds.get(b)
        if pa is None or pb is None:
            continue
        # census endpoint count for this edge (composite FA edges resolve to the one joint's endpoints).
        try:
            n_joints = max(int(dump.by_name(cc.name).n_endpoints), 8)
        except KeyError:
            n_joints = 24
        n_joints = min(n_joints, 60)                             # render budget (synthetic)
        if a == b:                                               # internal crosslink: join distinct nodes in one cloud
            ep = _nearest_pairs(pa[0::2], pa[1::2] if pa.shape[0] > 1 else pa, n_joints, seed=hash(cc.name) % 2**31)
        else:
            ep = _nearest_pairs(pa, pb, n_joints, seed=hash(cc.name) % 2**31)
        if not ep.size:
            continue
        fam = getattr(cc.family, "name", str(cc.family)).replace("_", "-")
        fam = {"FA-CLUTCH": "FA-clutch"}.get(fam, fam)
        load = RNG.uniform(1.0, 40.0, ep.shape[0] // 2)
        cons.append(_k(cc.name, fam, cc.name.replace("_", " "), a, b, ep, load=load))
    return cons


def dump_registered_connectors(dump, architecture: CellArchitecture) -> list[str]:
    """Connector names that appear (as a primary name or an alias) in the composed-world census."""
    names: set[str] = set()
    for e in dump.entries:
        if e.kind == "connector":
            names.add(e.name)
            names.update(e.aliases)
    return [c.name for c in architecture.connectors if c.name in names]


def build_dump(out_npz: str | Path) -> tuple[str, dict]:
    """Build the composed-world census, synthesise its geometry, and write the v2 dump. Returns (path, stats)."""
    architecture = reference_cell_architecture()
    world = build_census_world(architecture)
    dump = dump_composed_world_state(world)

    # The census invariants that make a per-compartment renderer safe (no double-draw) — asserted, not assumed.
    assert dump.actor_ids_are_unique(), "composed-world actor ids must be unique"
    assert dump.filament_id_blocks_are_disjoint(), "composed-world filament-ID blocks must be disjoint"

    comps, clouds = _component_specs(dump, architecture)
    cons = _connector_specs(dump, architecture, clouds)
    # Counts come FROM the architecture, never from a literal: the title read "32 conn" for three days
    # after connectors 32 -> 35 landed, while the panel beside it correctly drew 35.
    evidence = EvidenceLabel(
        rung=EvidenceRung.CENSUS_WIRED,
        quantitative=QuantitativeClaim.BLOCKED,
        basis=(
            f"census + topology for {len(architecture.components)} components / "
            f"{len(architecture.connectors)} connectors; geometry and forces are SYNTHETIC"
        ),
    )
    meta = {
        "stage": {
            "name": (
                f"COMPOSED ac/engine world ({len(architecture.components)} comp / "
                f"{len(architecture.connectors)} conn · real census)"
            ),
            **evidence.as_artifact_fields(),
            "note": "topology + counts + disjoint global filament-IDs from dump_composed_world_state; geometry "
                    "+ forces SYNTHETIC (real device dump = A5000 lane). cortex=KB native; MT/IF/filo count=PI-GAP.",
        },
        "report": {
            "n_actors": dump.n_actors, "n_filaments_total": dump.n_filaments_total,
            "n_nodes_total_census": dump.n_nodes_total,
            "actor_ids_unique": dump.actor_ids_are_unique(),
            "filament_blocks_disjoint": dump.filament_id_blocks_are_disjoint(),
            "n_components_drawn": len(comps), "n_connectors_drawn": len(cons),
        },
        "ledger": {},
    }
    path = write_v2_dump(out_npz, comps, cons, meta)
    stats = {"n_components": len(comps), "n_connectors": len(cons), "census_actors": dump.n_actors,
             "census_filaments_total": dump.n_filaments_total}
    return path, stats


def build_html(npz: str | Path, out_html: str | Path) -> tuple[str, dict]:
    """Render the composed-world v2 dump through the existing per-component/connector viewer."""
    from aleph.scripts.ac_cell_assembled_viz import build as _viz_build
    return _viz_build(Path(npz), Path(out_html))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    _root = Path(__file__).resolve().parents[1]
    ap.add_argument("--out", type=Path,
                    default=_root / "outputs" / "ac" / "cell_assembled" / "composed_world_v2.npz")
    ap.add_argument("--html", action="store_true", help="also render the interactive HTML viewer")
    ap.add_argument("--html-out", type=Path, default=None, help="HTML output path (default: alongside --out)")
    a = ap.parse_args()
    a.out.parent.mkdir(parents=True, exist_ok=True)
    path, stats = build_dump(a.out)
    print(f"wrote {path}")
    for k, v in stats.items():
        print(f"  {k:26s} {v:,}" if isinstance(v, int) else f"  {k:26s} {v}")
    if a.html:
        html_out = a.html_out or a.out.with_name("composed_world.html")
        hpath, counts = build_html(path, html_out)
        print(f"rendered {hpath}")
        for k, v in counts.items():
            print(f"  {k:34s} {v:,}" if isinstance(v, int) else f"  {k:34s} {v}")


if __name__ == "__main__":
    main()
