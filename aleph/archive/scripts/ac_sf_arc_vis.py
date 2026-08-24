r"""Milestone viz for the ``sf_arc`` stress-fiber component (DISJOINT population + KERNEL_BOUND mechanics).

This is the per-stage visualization checkpoint (PI 2026-07-22) for the two ``sf_arc`` landings this session:

* ``sf_population.build_sf_arc_population`` — the DISJOINT F-actin population (ventral / dorsal / transverse-arc /
  perinuclear-cap, each fiber owned by ``sf_arc`` and no other component; unique global IDs; cortex-disjoint).
* ``sf_mechanics`` (8ddfda15 / fc762243) — ``sf_arc`` SEAMED→KERNEL_BOUND: the component LAUNCHES its OWN
  rod-cable Warp kernels (``link_spring`` axial backbone + ``cytosim_bending`` NF2007), with the emergent
  contractile prestress coming from the discrete NMII cumsum (``stress_fiber_load_path``), NEVER a lumped k_SF.

Both builders are host NumPy (CPU-importable), so this viz builds the REAL native ``sf_arc`` geometry + its
emergent tension + the EXACT kernel-law restoring force ON THE DEV MAC — no gbook dump needed.  It emits ONE
self-contained interactive WebGL HTML (three.js, full native resolution, NO downsampling) via the shared
``ff_viewer_html.build_viewer``, reusing the ``ac_viz_common`` colour ramps / cbars / categorical colours so it
reads as one system with the assembled-cell viewer, and is verified in a REAL browser with
``scripts/browser_check.py`` (a grep of the HTML is NOT verification).

Scenes (each an independent per-layer-checkbox set):

* ``sf_arc population · disjoint (by class)`` — the four emergent SF classes, each its own colour + checkbox, so
  the reader SEES the disjoint population geometry (ventral basal, dorsal rising, arcs elevated, cap over the
  nucleus).  The class tag is a build-time GEOMETRY descriptor only — the firewall keeps it out of the force.
* ``sf_arc emergent contractile tension`` — every fiber coloured by the per-node emergent axial tension SHAPE
  from the NMII cumsum load path (``stress_fiber_load_path``), on a LINEAR ramp: − (cool) = compressive lobe at a
  dorsal free end, + (warm) = contractile tension.  The **magnitude is a PI-GAP** (per-head f_stall I0-B3 +
  engaged-head density I0-B6 unresolved), so the ramp is per-f_head (dimensionless) and labelled as such — this
  is the emergent SHAPE + sign, never a tuned number.  NMII bipolar seeds are overlaid as the motor stations.
* ``sf_arc connectors — the ONLY couplings`` — the four canonical ``sf_arc`` connector families' SF-side
  endpoints (FA anchor / cortex-transient / cap-LINC) + the INTERNAL dorsal↔arc crosslink drawn as real lines
  (both ends SF-owned).  Co-location in an array is NEVER a connection; the inter-component partner sides
  (focal_adhesion / cortex / nucleus) are owned by those components and appear on the composed gbook dump.
* ``sf_arc KERNEL_BOUND passive |F| (link+bend law)`` — every fiber coloured by the per-node restoring |F| from
  the EXACT NumPy reference of the launched kernels (``link_spring_force_reference`` + ``bending_force_reference``)
  under an illustrative small contraction, so the KERNEL_BOUND mechanics that landed is visible localizing at the
  anchors / motor stations.  Absolute magnitude is illustrative (axial k is a modelling GAP + the contraction is
  a diagnostic perturbation); the LOCALIZATION is the physics, log-turbo in pN.

    PYTHONPATH=/Users/sw1/ffn_cellsim python aleph/scripts/ac_sf_arc_vis.py
        -> aleph/outputs/ac/sf_arc/sf_arc.html  (+ browser-verified scene screenshots)
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

from aleph.engine.sf_mechanics import (
    bending_force_reference,
    build_sf_mechanics_topology,
    link_spring_force_reference,
)
from aleph.engine.sf_population import SFArcPopulation, SFClass, build_sf_arc_population
from aleph.scripts.ac_viz_common import cbar, lin_rgb, log_rgb
from aleph.scripts.ff_cell_morphology import uv_sphere
from aleph.scripts.ff_viewer_html import build_viewer

_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUT = _ROOT / "outputs" / "ac" / "sf_arc" / "sf_arc.html"

# Per-class categorical colours (held CLEAR of the turbo force ramp so a class is never misread as a hotspot).
_CLASS_COLOR = {
    SFClass.VENTRAL: "#ff8c42",          # basal, both-ends-FA (closed contractile dipole)
    SFClass.DORSAL: "#ffb454",           # basal FA → dorsal free end (opens into an arc)
    SFClass.TRANSVERSE_ARC: "#ffd166",   # elevated, no FA (couples dorsal free ends + cortex)
    SFClass.PERINUCLEAR_CAP: "#c98cff",  # arches over the nucleus, LINC-anchored
}
_CLASS_LABEL = {
    SFClass.VENTRAL: "ventral SF (FA↔FA, closed dipole)",
    SFClass.DORSAL: "dorsal SF (FA→free end)",
    SFClass.TRANSVERSE_ARC: "transverse arc (elevated, no FA)",
    SFClass.PERINUCLEAR_CAP: "perinuclear cap (LINC↔LINC over nucleus)",
}
# connector-family SF-side colours (mirror ac_viz_common.CONNECTOR_COLORS families where they exist).
_FA_COLOR, _LINC_COLOR, _CORTEX_COLOR, _ARC_COLOR = "#ffd60a", "#c98cff", "#ff5db1", "#4dd0e1"
_MOTOR_COLOR = "#ff35d6"

# emergent-tension display range [pN per pN-head] — the load path spans a small compressive lobe (~−0.75 at a
# dorsal free end) to full contractile tension (+1.0).  A DISPLAYED range, annotated; NOT tuned to an outcome.
_TENS_LO, _TENS_HI = -0.75, 1.0
# passive |F| display range [pN] for the KERNEL_BOUND restoring-force scene (log turbo). Interior colinear nodes
# cancel to ~0 (floor); the load is carried at ends / motor stations (warm) — the localization is the point.
_PF_LO, _PF_HI = 1.0e-2, 1.0e1
# Isotropic magnification for the single-sarcomere scene ONLY (a microscope objective, stated in the scene name):
# the viewer frames at the cell scale, so a 0.301 µm minifilament needs ~1e1 to be legible.  Uniform on all three
# axes ⇒ no axis truncation and no aspect distortion; every quoted length in the label is the TRUE µm value.
_ZOOM_MAGNIFICATION = 15.0


def _global_offsets(pop: SFArcPopulation) -> np.ndarray:
    """Flat per-fiber node offsets across the whole population (for consecutive-segment pairing)."""
    offs = [0]
    for b in pop.bundles:
        loc = np.asarray(b.fiber_offsets, np.int64)
        for f in range(loc.size - 1):
            offs.append(offs[-1] + int(loc[f + 1] - loc[f]))
    return np.asarray(offs, np.int64)


def _class_segments(pop: SFArcPopulation) -> dict[str, np.ndarray]:
    """Per-class (Nseg,2) consecutive-node global segment pairs (last node of each fiber skipped)."""
    by_class: dict[str, list[tuple[int, int]]] = {c: [] for c in SFClass.ALL}
    for b in pop.bundles:
        loc = np.asarray(b.fiber_offsets, np.int64)
        for f in range(loc.size - 1):
            lo = b.node_base + int(loc[f])
            hi = b.node_base + int(loc[f + 1])
            by_class[b.sf_class].extend((n, n + 1) for n in range(lo, hi - 1))
    return {c: np.asarray(v, np.int64).reshape(-1, 2) for c, v in by_class.items() if v}


def _node_tension_shape(pop: SFArcPopulation) -> np.ndarray:
    """Per-node emergent axial tension SHAPE (per f_head) from the NMII cumsum load path (magnitude = GAP).

    For each bundle the tension profile ``tension_per_fhead(s)`` is a function of the axial coordinate
    ``s = pos·axis``; each node is coloured by interpolating that profile at its own ``s`` — the same axial
    projection the load path uses.  Contractile (+), compressive lobe at a dorsal free end (−).
    """
    tens = np.zeros(pop.n_nodes, np.float64)
    for b, lp in zip(pop.bundles, pop.load_paths(f_head_pN=None), strict=True):
        lo, hi = b.node_base, b.node_base + b.pos.shape[0]
        s = pop.pos[lo:hi] @ np.asarray(lp.axis, np.float64)
        tens[lo:hi] = np.interp(s, lp.s_bins, lp.tension_per_fhead, left=0.0, right=0.0)
    return tens


def _passive_force_mag(pop: SFArcPopulation, *, k_axial_pn_per_um: float, contraction: float) -> np.ndarray:
    """Per-node restoring |F| [pN] from the EXACT launched-kernel law under an illustrative contraction.

    Builds the same disjoint rod-cable topology the KERNEL_BOUND ``SFFilamentMechanics`` uploads, contracts each
    bundle a small fraction toward its centroid (a labelled diagnostic — a stretched/contracted SF develops
    restoring force, the module's own sanity gate), and evaluates the NumPy references of the two launched
    kernels (byte-identical force law).  Absolute magnitude is illustrative (axial k is a modelling GAP); the
    localization at anchors / motor stations is the physics.
    """
    topo = build_sf_mechanics_topology(pop, k_axial_pn_per_um=k_axial_pn_per_um)
    pos = pop.pos.copy()
    for b in pop.bundles:
        lo, hi = b.node_base, b.node_base + b.pos.shape[0]
        centroid = pos[lo:hi].mean(0)
        pos[lo:hi] = centroid + (pos[lo:hi] - centroid) * (1.0 - contraction)
    f = link_spring_force_reference(pos, topo.links, topo.link_k, topo.link_r0)
    f = f + bending_force_reference(pos, topo.bend_triples, topo.bend_alpha)
    return np.linalg.norm(f, axis=1)


def _ref_layers(cell_r: float, nucleus_r: float) -> list[dict]:
    """Faint nucleus + R_cell scale-guide spheres (reference group)."""
    nv, nf = uv_sphere((0.0, 0.0, 0.0), nucleus_r, nu=36, nv=22)
    cv, cf = uv_sphere((0.0, 0.0, 0.0), cell_r, nu=48, nv=28)
    return [
        {"name": f"nucleus (R={nucleus_r:g}µm)", "kind": "mesh", "verts": nv, "faces": nf,
         "color": "#ffb454", "swatch": "#ffb454", "group": "reference", "opacity": 0.10, "clip": True},
        {"name": f"R_cell={cell_r:g}µm reference sphere", "kind": "mesh", "verts": cv, "faces": cf,
         "color": "#5a6b8c", "swatch": "#5a6b8c", "group": "reference", "opacity": 0.03, "clip": True},
    ]


def _lines_layer(name, verts, color, *, group="compartment", size=2.2, opacity=0.9, on_top=False,
                 color_frames=None, swatch=None) -> dict:
    layer = {"name": name, "kind": "lines", "verts": np.asarray(verts, np.float32).reshape(-1, 3),
             "color": color, "swatch": swatch or color, "group": group, "size": size,
             "opacity": opacity, "on_top": on_top, "clip": True}
    if color_frames is not None:
        layer["color"] = "#ffffff"
        layer["color_frames"] = [color_frames]
    return layer


def _points_layer(name, verts, color, *, group="connector", pt_size=0.35, opacity=1.0, swatch=None) -> dict:
    return {"name": name, "kind": "points", "verts": np.asarray(verts, np.float32).reshape(-1, 3),
            "color": color, "swatch": swatch or color, "group": group, "pt_size": pt_size, "opacity": opacity}


def build_scenes(pop: SFArcPopulation, *, cell_r: float, nucleus_r: float, k_axial: float,
                 contraction: float) -> tuple[dict, dict]:
    """Assemble the four ``sf_arc`` milestone scenes + their colour bars from the built population."""
    pos = np.asarray(pop.pos, np.float32)
    class_seg = _class_segments(pop)
    all_seg = np.concatenate([s for s in class_seg.values()], 0) if class_seg else np.zeros((0, 2), np.int64)
    scenes: dict[str, list] = {}
    cbars: dict[str, dict] = {}

    # NMII bipolar seed stations (the discrete motors whose cumsum IS the emergent prestress).
    myo_nodes = np.array(
        sorted({int(b.node_base + i) for b in pop.bundles for i in (*b.myo_i, *b.myo_j)}), np.int64)

    # ── scene 1: population by class (disjoint geometry) ─────────────────────────────────────────────
    s1: list[dict] = []
    for sf_class in SFClass.ALL:
        seg = class_seg.get(sf_class)
        if seg is None or not seg.size:
            continue
        n_fib = sum(b.global_fiber_ids.size for b in pop.bundles if b.sf_class == sf_class)
        s1.append(_lines_layer(f"{_CLASS_LABEL[sf_class]} · {n_fib} fil", pos[seg].reshape(-1, 3),
                               _CLASS_COLOR[sf_class], size=2.6, opacity=0.95))
    s1 += _ref_layers(cell_r, nucleus_r)
    scenes["sf_arc population · disjoint (by class)"] = s1

    # ── scene 2: emergent contractile tension (NMII cumsum SHAPE; magnitude PI-GAP) ──────────────────
    tens = _node_tension_shape(pop)
    tens_rgb = lin_rgb(tens, _TENS_LO, _TENS_HI)
    s2 = [_lines_layer("sf_arc fibers · emergent axial tension (− compress … + contractile)",
                       pos[all_seg].reshape(-1, 3), "#ffffff", size=3.0, opacity=1.0,
                       color_frames=tens_rgb[all_seg].reshape(-1, 3), swatch="#d62728")]
    s2.append(_points_layer(f"NMII bipolar seeds ({myo_nodes.size} motor stations)", pos[myo_nodes],
                            _MOTOR_COLOR, pt_size=0.5, group="compartment"))
    s2 += _ref_layers(cell_r, nucleus_r)
    k2 = "sf_arc emergent contractile tension (NMII cumsum · magnitude PI-GAP)"
    scenes[k2] = s2
    cbars[k2] = cbar("emergent axial tension / f_head", _TENS_LO, _TENS_HI, unit="pN·pN⁻¹", log=False)

    # ── scene 3: connectors — the ONLY couplings (co-location ≠ connection) ───────────────────────────
    s3 = [_lines_layer("sf_arc fibers (faint context)", pos[all_seg].reshape(-1, 3), "#39435a",
                       group="reference", size=0.8, opacity=0.14)]
    if pop.fa_sites.size:
        s3.append(_points_layer(f"fa_actin_anchor · SF↔focal_adhesion ({pop.fa_sites.size})",
                                pos[pop.fa_sites], _FA_COLOR, pt_size=0.6))
    if pop.linc_sites.size:
        s3.append(_points_layer(f"actin_cap_linc · SF↔nucleus ({pop.linc_sites.size})",
                                pos[pop.linc_sites], _LINC_COLOR, pt_size=0.6))
    if pop.cortex_sites.size:
        s3.append(_points_layer(f"sf_cortex_transient · SF↔cortex ({pop.cortex_sites.size})",
                                pos[pop.cortex_sites], _CORTEX_COLOR, pt_size=0.3))
    joints = np.asarray(pop.dorsal_arc_joints, np.int64).reshape(-1, 2)
    if joints.size:
        s3.append(_lines_layer(f"dorsal_arc_crosslink · SF↔SF internal ({joints.shape[0]})",
                               pos[joints].reshape(-1, 3), _ARC_COLOR, group="connector", size=3.4,
                               opacity=1.0, on_top=True))
    s3 += _ref_layers(cell_r, nucleus_r)
    scenes["sf_arc connectors — the ONLY couplings (co-location≠connection)"] = s3

    # ── scene 4: KERNEL_BOUND passive |F| (exact launched-kernel law; illustrative contraction) ───────
    pf = _passive_force_mag(pop, k_axial_pn_per_um=k_axial, contraction=contraction)
    pf_rgb = log_rgb(pf, _PF_LO, _PF_HI)
    s4 = [_lines_layer("sf_arc fibers · restoring |F| (link_spring + cytosim_bending)",
                       pos[all_seg].reshape(-1, 3), "#ffffff", size=3.0, opacity=1.0,
                       color_frames=pf_rgb[all_seg].reshape(-1, 3), swatch="#4dd0e1")]
    s4 += _ref_layers(cell_r, nucleus_r)
    k4 = "sf_arc KERNEL_BOUND passive |F| (link+bend law · illustrative contraction)"
    scenes[k4] = s4
    cbars[k4] = cbar("passive restoring |F|", _PF_LO, _PF_HI)

    # ── scenes 5–6: nmii_sf_motor — the bipolar straddle placement (only when the population carries it) ──
    motor = _motor_straddle_scene(pop, cell_r=cell_r, nucleus_r=nucleus_r)
    if motor is not None:
        scenes[motor[0]] = motor[1]
        zoom = _motor_zoom_scene(pop)
        if zoom is not None:
            scenes[zoom[0]] = zoom[1]
    return scenes, cbars


def _motor_zoom_scene(pop: SFArcPopulation) -> tuple[str, list[dict]] | None:
    """ONE sarcomere at its own scale — the only view at which the straddle geometry is actually legible.

    At whole-cell scale (R≈7.5 µm) a 0.3 µm minifilament is a dot, so the composed motor scene can show WHERE
    the motors sit but not WHETHER they straddle.  This scene contains a single sarcomere and its motor, so the
    viewer frames it at ~1 µm and the three load-bearing geometric facts become visible directly:

      * the two anti-parallel filaments run SIDE BY SIDE, separated by ``2·head_offset``, overlapping around the
        station (not meeting at a point, which is what the legacy geometry did);
      * the backbone sits BETWEEN them, and the ``+`` heads reach one filament while the ``−`` heads reach the
        other — so the motor pulls the pair together instead of pulling one filament two ways;
      * the barbed ends point OUTWARD (arrows at the fiber ends), which is what makes that pull contractile.
    """
    from aleph.engine.sf_nmii_population import (
        barbed_ward_tangents,
        build_sf_straddle_nmii_population,
    )
    from aleph.components.motor.minifilament_topology import MinifilamentTopology

    bundle = next((b for b in pop.bundles if b.motor_station_ready), None)
    if bundle is None:
        return None
    topology = MinifilamentTopology(
        n_bb=14, n_heads_per_side=10,
        backbone_length_um=float(bundle.sarcomere_overlap_um),
        head_offset_um=0.5 * float(bundle.sarcomere_lateral_um))
    nmii = build_sf_straddle_nmii_population(pop, topology, id_base=3_000_000)

    # the first motor + the bundle it straddles, re-centred on the station so the viewer frames ~1 µm.
    station = nmii.station_nodes[0]
    centre = np.asarray(pop.pos[station].mean(axis=0), np.float64)
    n_per_fiber = int(bundle.fiber_offsets[1])
    # Uniform ISOTROPIC magnification: the viewer frames at the cell scale (R≈7.5 µm), so a 0.3 µm motor is a
    # dot.  Scaling every coordinate of this one scene by the SAME factor is a microscope objective, not an axis
    # distortion (no axis is truncated, no aspect ratio changes); the factor is stated in the scene name.
    mag = float(_ZOOM_MAGNIFICATION)
    local = (np.asarray(bundle.pos, np.float64) - centre) * mag
    tangent = barbed_ward_tangents(bundle.pos, bundle.fiber_offsets, bundle.polarity)

    # Clip each filament to a WINDOW around the station, or the 13 µm sarcomere would frame the whole cell again
    # and the 0.4 µm separation would collapse to one line.  The window is 2x the backbone contour (so the whole
    # motor plus a margin is inside it) and the straight filaments are exact over it, so two endpoints suffice.
    window = 2.0 * float(bundle.sarcomere_overlap_um) * mag
    layers: list[dict] = []
    for f, (lo, hi, colour, label) in enumerate((
            (0, n_per_fiber, "#ff8c42", "filament A (barbed end OUTWARD ↙)"),
            (n_per_fiber, 2 * n_per_fiber, "#ffd166", "filament B (anti-parallel, barbed OUTWARD ↗)"))):
        axis = local[hi - 1] - local[lo]
        axis = axis / max(float(np.linalg.norm(axis)), 1e-12)
        # the point on THIS filament closest to the station (the station is the local origin), then ±window.
        base = local[lo] - float(local[lo] @ axis) * axis
        stub = np.vstack([base - window * axis, base + window * axis])
        layers.append(_lines_layer(label, stub.reshape(-1, 3), colour, size=4.0, opacity=1.0))
        # barbed-direction arrow stub, drawn AT the window edge on the barbed side (polarity decides which).
        t_barbed = tangent[lo if int(bundle.polarity[f]) < 0 else hi - 1]
        anchor = base + window * float(np.sign(t_barbed @ axis) or 1.0) * axis
        layers.append(_lines_layer(
            f"  ↳ barbed direction ({'A' if f == 0 else 'B'}) — OUTWARD",
            np.vstack([anchor, anchor + 0.35 * window * t_barbed]).reshape(-1, 3),
            "#ffffff", group="reference", size=2.4, opacity=0.9, on_top=True))

    mf_lo = int(nmii.minifilament_offset[0])
    mf_hi = int(nmii.minifilament_offset[1])
    mf_local = (np.asarray(nmii.pos[mf_lo:mf_hi], np.float64) - centre) * mag
    n_bb = int(topology.n_bb)
    bb_seg = np.column_stack([np.arange(n_bb - 1), np.arange(1, n_bb)])
    layers.append(_lines_layer("NMII backbone (Stam-Hocky, BETWEEN the two filaments)",
                               mf_local[bb_seg].reshape(-1, 3), _MOTOR_COLOR, group="compartment",
                               size=5.0, opacity=1.0, on_top=True))
    head_local = mf_local[np.asarray(nmii.head_node[:topology.n_heads], np.int64) - mf_lo]
    half = int(topology.n_heads_per_side)
    layers.append(_points_layer(f"'+' heads ON filament A ({half})", head_local[:half], _FA_COLOR,
                                pt_size=0.035, group="connector"))
    layers.append(_points_layer(f"'−' heads ON filament B ({half})", head_local[half:], _LINC_COLOR,
                                pt_size=0.035, group="connector"))
    # the head→backbone arms, so the reader sees each head is a real lever, not a free-floating dot.
    arms = np.concatenate([
        np.vstack([head_local[h], mf_local[int(nmii.head_bonds[h, 1]) - mf_lo]])
        for h in range(topology.n_heads)], 0)
    layers.append(_lines_layer("head↔backbone arms (r0 = head_offset, force-free at rest)",
                               arms.reshape(-1, 3), "#8899aa", group="connector", size=1.4, opacity=0.7))

    name = (f"nmii_sf_motor · ONE sarcomere ×{_ZOOM_MAGNIFICATION:g} isotropic (lateral "
            f"{bundle.sarcomere_lateral_um:g} µm = 2·head_offset · overlap "
            f"{bundle.sarcomere_overlap_um:g} µm = backbone contour)")
    return name, layers


def _motor_straddle_scene(
    pop: SFArcPopulation, *, cell_r: float, nucleus_r: float,
) -> tuple[str, list[dict]] | None:
    """The ``nmii_sf_motor`` scene: explicit minifilaments STRADDLING each straight sarcomere's pair.

    Guarded — it appears only for a population built with the straddle sarcomere geometry (the legacy
    end-to-end population has its two filaments co-located at ``mid`` and cannot host a motor).  What the
    reader must be able to SEE, because these are exactly the claims the native gate rests on:

      * each minifilament's backbone lies BETWEEN two anti-parallel filaments, its ``+`` heads touching one and
        its ``−`` heads the other (two separate colours) — the bipolar straddle, so the pair is pulled INWARD;
      * the heads sit ON the actin lines, not floating near them (they are drawn at full resolution, and the
        measured perpendicular residual is in the scene name);
      * the CURVED perinuclear cap carries NO motor and is drawn desaturated with its exclusion in the label —
        an arch's two filaments are not anti-parallel over the shared span, so a rigid bipolar minifilament
        cannot straddle it (measured bipolar dot −0.600 vs −1 for every straight class).
    """
    if not pop.straddle_ready:
        return None
    from aleph.engine.sf_nmii_population import build_sf_straddle_nmii_population
    from aleph.engine.sf_motor_slice import sf_sarcomere_geometry_for
    from aleph.components.motor.minifilament_topology import MinifilamentTopology

    bundle = pop.bundles[0]
    head_offset = 0.5 * float(bundle.sarcomere_lateral_um)
    topology = MinifilamentTopology(
        n_bb=14, n_heads_per_side=10,
        backbone_length_um=float(bundle.sarcomere_overlap_um), head_offset_um=head_offset)
    if abs(sf_sarcomere_geometry_for(topology)["sarcomere_lateral_um"]
           - float(bundle.sarcomere_lateral_um)) > 1e-9:
        return None
    nmii = build_sf_straddle_nmii_population(pop, topology, id_base=3_000_000)

    pos = np.asarray(pop.pos, np.float32)
    class_seg = _class_segments(pop)
    ready_classes = set(pop.motor_station_census()["ready_by_class"])
    layers: list[dict] = []
    for sf_class in SFClass.ALL:
        seg = class_seg.get(sf_class)
        if seg is None or not seg.size:
            continue
        hosts = sf_class in ready_classes
        label = (f"{_CLASS_LABEL[sf_class]} · hosts a motor" if hosts else
                 f"{_CLASS_LABEL[sf_class]} · NO motor (curved sarcomere — cannot be straddled)")
        layers.append(_lines_layer(label, pos[seg].reshape(-1, 3), _CLASS_COLOR[sf_class],
                                   size=2.6 if hosts else 1.4, opacity=0.95 if hosts else 0.22))

    mf_pos = np.asarray(nmii.pos, np.float32)
    n_bb = int(topology.n_bb)
    backbone_verts = np.concatenate([
        mf_pos[int(nmii.minifilament_offset[m]) + np.column_stack(
            [np.arange(n_bb - 1), np.arange(1, n_bb)]).reshape(-1)]
        for m in range(nmii.n_minifilaments)], 0) if nmii.n_minifilaments else np.zeros((0, 3), np.float32)
    layers.append(_lines_layer(f"NMII backbones ({nmii.n_minifilaments} minifilaments, Stam-Hocky bipolar)",
                               backbone_verts, _MOTOR_COLOR, group="compartment", size=3.6, opacity=1.0,
                               on_top=True))
    head_pos = mf_pos[np.asarray(nmii.head_node, np.int64)]
    side = np.asarray(nmii.head_side, np.int64)
    layers.append(_points_layer(f"'+' heads → filament A ({int((side == 0).sum())})",
                                head_pos[side == 0], _FA_COLOR, pt_size=0.42, group="connector"))
    layers.append(_points_layer(f"'−' heads → filament B ({int((side == 1).sum())})",
                                head_pos[side == 1], _LINC_COLOR, pt_size=0.42, group="connector"))
    layers.append(_points_layer(f"motor stations ({nmii.station_nodes.shape[0]} anti-parallel node pairs)",
                                pos[np.unique(nmii.station_nodes.reshape(-1))], "#ffffff",
                                pt_size=0.30, group="connector", opacity=0.85))
    layers += _ref_layers(cell_r, nucleus_r)

    excluded = nmii.station_census["excluded_by_class"] or {"none": 0}
    name = (f"nmii_sf_motor · bipolar straddle ({nmii.n_heads} heads on "
            f"{nmii.n_minifilaments} minifilaments · excluded: "
            f"{', '.join(f'{k}×{v}' for k, v in excluded.items())})")
    return name, layers


def compose_title(pop: SFArcPopulation, *, k_axial: float, contraction: float) -> str:
    c = pop.census()
    by = c["N_by_class"]
    return (
        f"sf_arc milestone — DISJOINT population + KERNEL_BOUND mechanics · "
        f"{c['N_unique_sf_fibers']} fibers ({c['N_unique_sf_fibers']}=unique, no cortex double-count) · "
        f"classes {{V:{by.get('ventral', 0)} D:{by.get('dorsal', 0)} "
        f"arc:{by.get('transverse_arc', 0)} cap:{by.get('perinuclear_cap', 0)}}} · "
        f"tension SHAPE per f_head (magnitude=PI-GAP I0-B3/B6) · passive |F| illustrative "
        f"(k_axial={k_axial:g} pN/µm GAP, {100 * contraction:g}% contraction) · units µm · FULL-RES (no downsample)"
    )


def build(out_html: Path, *, n_ventral: int, n_dorsal: int, n_arc: int, n_cap: int, n_per_fiber: int,
          cell_r: float, nucleus_r: float, k_axial: float, contraction: float,
          with_motor: bool = True) -> tuple[str, dict]:
    """Build the ``sf_arc`` population on the dev Mac and emit the interactive milestone viewer HTML."""
    straddle: dict[str, float] = {}
    if with_motor:
        from aleph.engine.sf_motor_slice import sf_sarcomere_geometry_for
        from aleph.components.motor.minifilament_topology import MinifilamentTopology

        # the sarcomere geometry is DERIVED from the minifilament that must sit in it, never chosen here.
        straddle = sf_sarcomere_geometry_for(MinifilamentTopology(
            n_bb=14, n_heads_per_side=10, backbone_length_um=0.301, head_offset_um=0.200))
    pop = build_sf_arc_population(
        n_ventral=n_ventral, n_dorsal=n_dorsal, n_arc=n_arc, n_cap=n_cap, n_per_fiber=n_per_fiber,
        cell_radius_um=cell_r, nucleus_radius_um=nucleus_r, **straddle)
    pop.assert_partitioned()                       # no fiber double-counted (each owned by exactly one component)
    scenes, cbars = build_scenes(
        pop, cell_r=cell_r, nucleus_r=nucleus_r, k_axial=k_axial, contraction=contraction)
    build_viewer(scenes, out=str(out_html), title=compose_title(pop, k_axial=k_axial, contraction=contraction),
                 cbars=cbars)
    census = pop.census()
    counts = {"n_fibers": census["N_unique_sf_fibers"], "n_nodes": pop.n_nodes,
              "n_scenes": len(scenes), "by_class": census["N_by_class"],
              "magnitude_status": census["magnitude_status"]}
    if with_motor:
        counts["motor_stations"] = pop.motor_station_census()
    return str(out_html), counts


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    ap.add_argument("--n-ventral", type=int, default=6)
    ap.add_argument("--n-dorsal", type=int, default=4)
    ap.add_argument("--n-arc", type=int, default=3)
    ap.add_argument("--n-cap", type=int, default=3)
    ap.add_argument("--n-per-fiber", type=int, default=7)
    ap.add_argument("--cell-radius", type=float, default=7.5)
    ap.add_argument("--nucleus-radius", type=float, default=3.0)
    # illustrative-only knobs for the KERNEL_BOUND passive-|F| scene (magnitude is a documented GAP).
    ap.add_argument("--k-axial", type=float, default=100.0, help="illustrative axial k [pN/µm] (modelling GAP)")
    ap.add_argument("--contraction", type=float, default=0.04, help="illustrative bundle contraction fraction")
    ap.add_argument("--no-motor", action="store_true",
                    help="omit the nmii_sf_motor straddle scene (legacy end-to-end sarcomere geometry)")
    a = ap.parse_args()
    a.out.parent.mkdir(parents=True, exist_ok=True)
    path, counts = build(
        a.out, n_ventral=a.n_ventral, n_dorsal=a.n_dorsal, n_arc=a.n_arc, n_cap=a.n_cap,
        n_per_fiber=a.n_per_fiber, cell_r=a.cell_radius, nucleus_r=a.nucleus_radius,
        k_axial=a.k_axial, contraction=a.contraction, with_motor=not a.no_motor)
    print(f"wrote {path}")
    for k, v in counts.items():
        print(f"  {k:20s} {v}")


if __name__ == "__main__":
    main()
