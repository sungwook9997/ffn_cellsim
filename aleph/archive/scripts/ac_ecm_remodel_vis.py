"""Visualize the real CUDA B2 ECM topology transaction.

This diagnostic runs the committed ``ECMRemodelRuntime`` on an NVIDIA Warp CUDA
device and renders a common-axis step-through of baseline, proposal, rejection,
accepted refinement, coarsening proposal, and accepted coarsening.  Host reads occur
only after each completed diagnostic transaction for plotting; they are not simulation
state and are never used to decide a device event.

The figure is deliberately a topology/remap visualization.  It contains no collagen
constitutive law, traction response, damage threshold, or FA chemistry.

Usage (gbook A5000):
    python -m aleph.scripts.ac_ecm_remodel_vis --device cuda:0
"""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import matplotlib
import numpy as np
import warp as wp

matplotlib.use("Agg")
import matplotlib.animation as animation  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

from aleph.components.ecm import (  # noqa: E402
    ECMRemodelRuntime,
    MikadoInitConfig,
    MikadoTopologyBuilder,
    RemodelEvent,
)


@wp.kernel
def _set_i32_kernel(array: wp.array(dtype=wp.int32), index: wp.int32, value: wp.int32) -> None:
    if wp.tid() == 0:
        array[index] = value


@wp.kernel
def _set_f64_kernel(array: wp.array(dtype=wp.float64), index: wp.int32, value: wp.float64) -> None:
    if wp.tid() == 0:
        array[index] = value


def _set_i32(array: wp.array, index: int, value: int, *, device: str) -> None:
    wp.launch(
        _set_i32_kernel,
        dim=1,
        inputs=[array, wp.int32(index), wp.int32(value)],
        device=device,
    )


def _set_f64(array: wp.array, index: int, value: float, *, device: str) -> None:
    wp.launch(
        _set_f64_kernel,
        dim=1,
        inputs=[array, wp.int32(index), wp.float64(value)],
        device=device,
    )


@dataclass(frozen=True, slots=True)
class Projection:
    """Target-fiber-aligned orthographic projection."""

    origin: np.ndarray
    axis_x: np.ndarray
    axis_y: np.ndarray

    def apply(self, position: np.ndarray) -> np.ndarray:
        centered = np.asarray(position, dtype=np.float64) - self.origin
        return np.column_stack((centered @ self.axis_x, centered @ self.axis_y))


def _projection_from_target(topology: Any, fiber: int) -> Projection:
    positions = topology.position_d.numpy()
    segments = topology.segments_d.numpy()
    active = topology.segment_active_d.numpy().astype(bool)
    owners = topology.segment_fiber_d.numpy()
    previous = topology.segment_prev_d.numpy()
    following = topology.segment_next_d.numpy()
    head = int(np.flatnonzero(active & (owners == fiber) & (previous < 0))[0])
    tail = head
    while following[tail] >= 0:
        tail = int(following[tail])
    node_a = int(segments[head, 0])
    node_b = int(segments[tail, 1])
    axis_x = positions[node_b] - positions[node_a]
    axis_x /= np.linalg.norm(axis_x)
    helper = np.array([0.0, 0.0, 1.0])
    if abs(float(axis_x @ helper)) > 0.9:
        helper = np.array([0.0, 1.0, 0.0])
    axis_y = np.cross(helper, axis_x)
    axis_y /= np.linalg.norm(axis_y)
    origin = 0.5 * (positions[node_a] + positions[node_b])
    return Projection(origin=origin, axis_x=axis_x, axis_y=axis_y)


def _chain_slots(topology: Any, fiber: int) -> list[int]:
    active = topology.segment_active_d.numpy().astype(bool)
    owners = topology.segment_fiber_d.numpy()
    previous = topology.segment_prev_d.numpy()
    following = topology.segment_next_d.numpy()
    heads = np.flatnonzero(active & (owners == fiber) & (previous < 0))
    if heads.size != 1:
        raise RuntimeError(f"fiber {fiber} must have exactly one active chain head")
    chain: list[int] = []
    segment = int(heads[0])
    while segment >= 0:
        if segment in chain:
            raise RuntimeError("cycle detected in ECM chain")
        chain.append(segment)
        segment = int(following[segment])
    return chain


def _longest_fiber(topology: Any) -> int:
    """Choose the least-clipped initialized fiber using post-build geometry only."""
    positions = topology.position_d.numpy()
    segments = topology.segments_d.numpy()
    lengths: list[float] = []
    for fiber in range(topology.n_fiber_capacity):
        chain = _chain_slots(topology, fiber)
        node_a = int(segments[chain[0], 0])
        node_b = int(segments[chain[-1], 1])
        lengths.append(float(np.linalg.norm(positions[node_b] - positions[node_a])))
    return int(np.argmax(lengths))


def _material_point(
    positions: np.ndarray,
    segments: np.ndarray,
    segment: int,
    u: float,
) -> np.ndarray:
    node_a, node_b = segments[segment]
    return (1.0 - u) * positions[node_a] + u * positions[node_b]


def _capture_state(
    *,
    label: str,
    subtitle: str,
    topology: Any,
    projection: Projection,
    target_fiber: int,
    endpoint_specs: list[tuple[str, int, float]],
    highlight_slots: list[int],
    plan: dict[str, Any] | None = None,
) -> dict[str, Any]:
    wp.synchronize_device(topology.device)
    positions_3d = topology.position_d.numpy()
    positions_2d = projection.apply(positions_3d)
    segments = topology.segments_d.numpy()
    active = topology.segment_active_d.numpy().astype(bool)
    owners = topology.segment_fiber_d.numpy()
    levels = topology.segment_refinement_level_d.numpy()
    damage = topology.segment_damage_d.numpy()
    persistent = topology.persistent_segment_id_d.numpy()
    generations = topology.segment_generation_d.numpy()
    segment_rows: list[dict[str, Any]] = []
    for slot in np.flatnonzero(active):
        node_a, node_b = segments[slot]
        segment_rows.append(
            {
                "slot": int(slot),
                "fiber": int(owners[slot]),
                "a": positions_2d[node_a].round(9).tolist(),
                "b": positions_2d[node_b].round(9).tolist(),
                "level": int(levels[slot]),
                "damage": float(damage[slot]),
                "persistent_id": int(persistent[slot]),
                "generation": int(generations[slot]),
            }
        )
    endpoints: list[dict[str, Any]] = []
    for endpoint_label, segment, u in endpoint_specs:
        point = _material_point(positions_3d, segments, segment, u)
        point_2d = projection.apply(point[None, :])[0]
        endpoints.append(
            {
                "label": endpoint_label,
                "segment": int(segment),
                "u": float(u),
                "xy": point_2d.round(9).tolist(),
            }
        )
    return {
        "label": label,
        "subtitle": subtitle,
        "segments": segment_rows,
        "target_chain": _chain_slots(topology, target_fiber),
        "highlight_slots": [int(slot) for slot in highlight_slots],
        "endpoints": endpoints,
        "active_population": topology.active_population_d.numpy().astype(int).tolist(),
        "free_population": topology.free_population_d.numpy().astype(int).tolist(),
        "next_persistent_id": int(topology.next_persistent_id_d.numpy()[0]),
        "topology_epoch": int(topology.topology_epoch_d.numpy()[0]),
        "topology_dirty": int(topology.topology_dirty_d.numpy()[0]),
        "plan": plan,
    }


def _proposal_plan(runtime: ECMRemodelRuntime, source: int) -> dict[str, Any]:
    remap = runtime.endpoint_remap
    return {
        "source_slot": source,
        "source_persistent_id": int(remap.source_persistent_id_d.numpy()[source]),
        "new_node_slot": int(runtime.plan.new_node_slot_d.numpy()[source]),
        "target_a_slot": int(remap.target_a_segment_d.numpy()[source]),
        "target_b_slot": int(remap.target_b_segment_d.numpy()[source]),
        "target_a_persistent_id": int(remap.target_a_persistent_id_d.numpy()[source]),
        "target_b_persistent_id": int(remap.target_b_persistent_id_d.numpy()[source]),
        "candidate_valid": int(runtime.candidate_valid_d.numpy()[0]),
    }


def _ack_split(runtime: ECMRemodelRuntime, source: int, *, device: str) -> None:
    _set_i32(runtime.endpoint_remap.ack_count_d, source, 3, device=device)
    _set_i32(runtime.endpoint_remap.target_a_refcount_d, source, 2, device=device)
    _set_i32(runtime.endpoint_remap.target_b_refcount_d, source, 1, device=device)


def _ack_merge(runtime: ECMRemodelRuntime, left: int, right: int, *, device: str) -> None:
    _set_i32(runtime.endpoint_remap.ack_count_d, left, 2, device=device)
    _set_i32(runtime.endpoint_remap.target_a_refcount_d, left, 2, device=device)
    _set_i32(runtime.endpoint_remap.ack_count_d, right, 1, device=device)
    _set_i32(runtime.endpoint_remap.target_a_refcount_d, right, 1, device=device)


def _submit_refine(
    runtime: ECMRemodelRuntime,
    source: int,
    *,
    device: str,
) -> None:
    runtime.snapshot_candidate()
    _set_i32(runtime.proposal.segment_event_d, source, RemodelEvent.REFINE, device=device)
    _set_f64(runtime.proposal.damage_increment_d, source, 0.25, device=device)
    runtime.prepare_candidate()


def _submit_coarsen(
    runtime: ECMRemodelRuntime,
    left: int,
    right: int,
    *,
    device: str,
) -> None:
    runtime.snapshot_candidate()
    _set_i32(runtime.proposal.segment_event_d, left, RemodelEvent.COARSEN, device=device)
    _set_f64(runtime.proposal.damage_increment_d, left, 0.10, device=device)
    _set_f64(runtime.proposal.damage_increment_d, right, 0.30, device=device)
    runtime.prepare_candidate()


def _plot_state(
    ax: plt.Axes,
    state: dict[str, Any],
    *,
    limits: tuple[float, float, float, float],
    show_axis_labels: bool,
    show_endpoint_labels: bool,
) -> None:
    target_chain = set(state["target_chain"])
    highlights = set(state["highlight_slots"])
    for segment in state["segments"]:
        start = segment["a"]
        end = segment["b"]
        slot = segment["slot"]
        if slot in highlights:
            color, width, alpha = "#d62728", 4.2, 1.0
        elif slot in target_chain:
            color, width, alpha = "#007f8b", 2.8, 1.0
        else:
            color, width, alpha = "#8d99ae", 0.7, 0.28
        ax.plot(
            (start[0], end[0]),
            (start[1], end[1]),
            color=color,
            linewidth=width,
            alpha=alpha,
            solid_capstyle="round",
            zorder=2 if slot in target_chain else 1,
        )
    label_offsets = ((-18, 14), (0, 30), (18, 14))
    for endpoint_index, endpoint in enumerate(state["endpoints"]):
        ax.scatter(
            endpoint["xy"][0],
            endpoint["xy"][1],
            s=58,
            facecolor="#ffbf00",
            edgecolor="#3b2f00",
            linewidth=0.9,
            zorder=6,
        )
        if show_endpoint_labels:
            ax.annotate(
                f"{endpoint['label']}\nu={endpoint['u']:.2f}",
                endpoint["xy"],
                xytext=label_offsets[endpoint_index % len(label_offsets)],
                textcoords="offset points",
                ha="center",
                va="bottom",
                fontsize=7,
                arrowprops={"arrowstyle": "-", "color": "#6b5b00", "linewidth": 0.6},
                zorder=7,
            )
    if state["plan"] is not None and state["plan"].get("prospective_xy") is not None:
        planned = state["plan"]["prospective_xy"]
        ax.plot(
            (planned["a"][0], planned["mid"][0], planned["b"][0]),
            (planned["a"][1], planned["mid"][1], planned["b"][1]),
            color="#9467bd",
            linewidth=5.5,
            linestyle=(0, (3, 2)),
            zorder=4,
        )
        ax.scatter(
            planned["mid"][0],
            planned["mid"][1],
            marker="D",
            s=52,
            color="#9467bd",
            edgecolor="#3d2457",
            zorder=5,
        )
    if state.get("accepted_midpoint_xy") is not None:
        midpoint = state["accepted_midpoint_xy"]
        ax.scatter(
            midpoint[0],
            midpoint[1],
            marker="D",
            s=56,
            color="#9467bd",
            edgecolor="#3d2457",
            zorder=5,
        )
    ax.set_xlim(limits[0], limits[1])
    ax.set_ylim(limits[2], limits[3])
    ax.set_aspect("equal", adjustable="box")
    ax.grid(alpha=0.16, linewidth=0.5)
    ax.set_title(f"{state['label']}\n{state['subtitle']}", fontsize=10, fontweight="bold")
    if show_axis_labels:
        ax.set_xlabel("target-fiber axis [µm]")
        ax.set_ylabel("orthogonal projection [µm]")
    populations = state["active_population"]
    ax.text(
        0.015,
        0.018,
        (
            f"active N/S/B={populations[1]}/{populations[2]}/{populations[3]}  "
            f"epoch={state['topology_epoch']}  dirty={state['topology_dirty']}"
        ),
        transform=ax.transAxes,
        fontsize=7.5,
        ha="left",
        va="bottom",
        bbox={"boxstyle": "round,pad=0.25", "facecolor": "white", "alpha": 0.82},
    )


def _common_limits(states: list[dict[str, Any]]) -> tuple[float, float, float, float]:
    points = np.array(
        [
            point
            for state in states
            for segment in state["segments"]
            for point in (segment["a"], segment["b"])
        ],
        dtype=np.float64,
    )
    lo = points.min(axis=0)
    hi = points.max(axis=0)
    span = np.maximum(hi - lo, 1.0)
    pad = 0.04 * span
    return (
        float(lo[0] - pad[0]),
        float(hi[0] + pad[0]),
        float(lo[1] - pad[1]),
        float(hi[1] + pad[1]),
    )


def _target_limits(states: list[dict[str, Any]]) -> tuple[float, float, float, float]:
    points: list[list[float]] = []
    for state in states:
        target_chain = set(state["target_chain"])
        for segment in state["segments"]:
            if segment["slot"] in target_chain:
                points.extend((segment["a"], segment["b"]))
    array = np.asarray(points, dtype=np.float64)
    lo = array.min(axis=0)
    hi = array.max(axis=0)
    x_pad = max(0.35, 0.08 * (hi[0] - lo[0]))
    y_center = 0.5 * (lo[1] + hi[1])
    y_half = max(0.85, 0.15 * (hi[0] - lo[0]))
    return (
        float(lo[0] - x_pad),
        float(hi[0] + x_pad),
        float(y_center - y_half),
        float(y_center + y_half),
    )


def _render_outputs(states: list[dict[str, Any]], output_dir: Path) -> None:
    figure_dir = output_dir / "figs"
    figure_dir.mkdir(parents=True, exist_ok=True)
    global_limits = _common_limits(states)
    target_limits = _target_limits(states)
    fig, axes = plt.subplots(2, 3, figsize=(16, 10), constrained_layout=True)
    for index, (ax, state) in enumerate(zip(axes.flat, states, strict=True)):
        _plot_state(
            ax,
            state,
            limits=global_limits if index == 0 else target_limits,
            show_axis_labels=index % 3 == 0,
            show_endpoint_labels=index != 0,
        )
    legend = (
        Line2D([0], [0], color="#8d99ae", lw=1.0, alpha=0.5, label="other ECM fibers"),
        Line2D([0], [0], color="#007f8b", lw=3.0, label="tracked fiber"),
        Line2D([0], [0], color="#d62728", lw=4.0, label="edited segment"),
        Line2D([0], [0], color="#9467bd", lw=5.0, ls="--", label="planned split"),
        Line2D(
            [0],
            [0],
            color="#ffbf00",
            marker="o",
            markeredgecolor="#3b2f00",
            lw=0,
            label="graph endpoint",
        ),
    )
    fig.legend(handles=legend, loc="lower center", ncol=5, frameon=False)
    fig.suptitle(
        "ECM B2 — real Warp-CUDA topology transaction (A5000)\n"
        "topology/remap only; no constitutive traction or damage-threshold claim",
        fontsize=15,
        fontweight="bold",
    )
    fig.savefig(figure_dir / "ecm_b2_transaction.png", dpi=180)
    plt.close(fig)

    animation_fig, animation_ax = plt.subplots(figsize=(11, 7), constrained_layout=True)

    def update(frame: int) -> tuple[plt.Axes]:
        animation_ax.clear()
        _plot_state(
            animation_ax,
            states[frame],
            limits=target_limits,
            show_axis_labels=True,
            show_endpoint_labels=True,
        )
        animation_fig.suptitle(
            "ECM B2 transaction — actual CUDA states\n"
            "yellow endpoints retain their physical material points across remap",
            fontsize=13,
            fontweight="bold",
        )
        return (animation_ax,)

    movie = animation.FuncAnimation(
        animation_fig,
        update,
        frames=len(states),
        interval=1350,
        repeat=True,
        blit=False,
    )
    movie.save(figure_dir / "ecm_b2_transaction.gif", writer=animation.PillowWriter(fps=0.75))
    plt.close(animation_fig)


def _write_report(
    *,
    output_dir: Path,
    states: list[dict[str, Any]],
    reject_max_delta_um: float,
    source_rest_before_um: float,
    source_rest_after_um: float,
) -> None:
    baseline = states[0]
    refined = states[3]
    coarsened = states[-1]
    report = f"""# ECM B2 transaction visualization

## Scope

This is an RTX A5000 execution of the Warp-CUDA ECM topology transaction. It visualizes
proposal planning, explicit rejection, accepted refinement with graph endpoint remap,
and accepted coarsening. It is not a collagen constitutive, traction, FA-chemistry, or
production-material validation.

## Sanity measurements

- Rejected maximum position change: `{reject_max_delta_um:.17g} µm`.
- Source rest length before refine: `{source_rest_before_um:.17g} µm`.
- Rest length after refine→coarsen: `{source_rest_after_um:.17g} µm`.
- Active node/segment/bend populations: baseline `{baseline["active_population"][1:]}`;
  refined `{refined["active_population"][1:]}`; coarsened `{coarsened["active_population"][1:]}`.
- Topology epochs: `{[state["topology_epoch"] for state in states]}`.
- Final endpoint population is conserved at 3.

## Figures

- `figs/ecm_b2_transaction.png` — common-axis six-state transaction audit.
- `figs/ecm_b2_transaction.gif` — animated step-through of the same actual CUDA states.
- `ecm_b2_transaction.json` — projected segment topology, endpoint coordinates, IDs,
  generations, populations, epochs, and proposal plan used by the figures.

## Held

The collagen modulus conflict (30–100 Pa versus 5–100 Pa) remains HELD. No force law or
damage threshold was introduced to make this visualization.
"""
    (output_dir / "REPORT.md").write_text(report, encoding="utf-8")


def run(*, device: str, output_dir: Path) -> None:
    cfg = MikadoInitConfig(
        box_lo_um=(-7.0, -7.0, -7.0),
        box_hi_um=(7.0, 7.0, 7.0),
        n_fibers=18,
        fiber_length_um=6.0,
        target_segment_um=1.5,
        crosslink_capture_um=0.4,
        pin_faces=("x_lo", "x_hi", "y_lo", "y_hi", "z_lo", "z_hi"),
        pin_margin_um=0.15,
        rng_seed=741,
        max_refinement_level=2,
        persistent_id_base=90_000,
    )
    topology = MikadoTopologyBuilder(cfg, device=device).initialize()
    runtime = ECMRemodelRuntime(topology, max_refinement_level=cfg.max_refinement_level)
    target_fiber = _longest_fiber(topology)
    source = target_fiber * cfg.segment_slots_per_fiber + (cfg.nodes_per_fiber - 2) // 2
    projection = _projection_from_target(topology, target_fiber)
    _set_i32(topology.endpoint_refcount_d, source, 3, device=device)
    baseline_position = topology.position_d.numpy().copy()
    source_rest_before = float(topology.segment_rest_length_d.numpy()[source])
    baseline_endpoints = [("E1", source, 0.20), ("E2", source, 0.50), ("E3", source, 0.85)]
    states: list[dict[str, Any]] = [
        _capture_state(
            label="0 · Baseline",
            subtitle="3 graph endpoints on one material segment",
            topology=topology,
            projection=projection,
            target_fiber=target_fiber,
            endpoint_specs=baseline_endpoints,
            highlight_slots=[source],
        )
    ]

    _submit_refine(runtime, source, device=device)
    plan = _proposal_plan(runtime, source)
    positions = topology.position_d.numpy()
    segments = topology.segments_d.numpy()
    node_a, node_b = segments[source]
    prospective = projection.apply(
        np.vstack(
            (positions[node_a], 0.5 * (positions[node_a] + positions[node_b]), positions[node_b])
        )
    )
    plan["prospective_xy"] = {
        "a": prospective[0].round(9).tolist(),
        "mid": prospective[1].round(9).tolist(),
        "b": prospective[2].round(9).tolist(),
    }
    states.append(
        _capture_state(
            label="1 · CUDA proposal",
            subtitle="fresh slots/IDs planned; topology still unchanged",
            topology=topology,
            projection=projection,
            target_fiber=target_fiber,
            endpoint_specs=baseline_endpoints,
            highlight_slots=[source],
            plan=plan,
        )
    )
    _ack_split(runtime, source, device=device)
    runtime.validate_remap_acknowledgements()
    rejected = wp.array([0], dtype=wp.int32, device=device)
    runtime.rollback(rejected)
    runtime.commit_irreversible(rejected, dt_phys=0.01, rng_seed=1)
    wp.synchronize_device(device)
    reject_max_delta = float(np.max(np.abs(topology.position_d.numpy() - baseline_position)))
    states.append(
        _capture_state(
            label="2 · Rejected",
            subtitle=f"bit-exact rollback: max |Δx|={reject_max_delta:.1e} µm",
            topology=topology,
            projection=projection,
            target_fiber=target_fiber,
            endpoint_specs=baseline_endpoints,
            highlight_slots=[source],
        )
    )

    _submit_refine(runtime, source, device=device)
    _ack_split(runtime, source, device=device)
    runtime.validate_remap_acknowledgements()
    right = int(runtime.endpoint_remap.target_b_segment_d.numpy()[source])
    accepted = wp.array([1], dtype=wp.int32, device=device)
    runtime.rollback(accepted)
    runtime.commit_irreversible(accepted, dt_phys=0.01, rng_seed=2)
    refined_positions = topology.position_d.numpy()
    refined_segments = topology.segments_d.numpy()
    refined_midpoint = int(refined_segments[source, 1])
    refined_midpoint_xy = projection.apply(refined_positions[refined_midpoint][None, :])[0]
    refined_endpoints = [("E1", source, 0.40), ("E2", source, 1.00), ("E3", right, 0.70)]
    refined_state = _capture_state(
        label="3 · Accepted refine",
        subtitle="left/right generations and IDs are fresh; endpoints conserved 2+1",
        topology=topology,
        projection=projection,
        target_fiber=target_fiber,
        endpoint_specs=refined_endpoints,
        highlight_slots=[source, right],
    )
    refined_state["accepted_midpoint_xy"] = refined_midpoint_xy.round(9).tolist()
    states.append(refined_state)

    _submit_coarsen(runtime, source, right, device=device)
    merge_plan = _proposal_plan(runtime, source)
    merge_positions = topology.position_d.numpy()
    merge_segments = topology.segments_d.numpy()
    merge_nodes = np.asarray((merge_segments[source, 0], merge_segments[right, 1]))
    merge_xy = projection.apply(merge_positions[merge_nodes])
    merge_plan["prospective_xy"] = {
        "a": merge_xy[0].round(9).tolist(),
        "mid": refined_midpoint_xy.round(9).tolist(),
        "b": merge_xy[1].round(9).tolist(),
    }
    states.append(
        _capture_state(
            label="4 · Coarsen proposal",
            subtitle="sibling root/path verified; graph remaps both sources to left slot",
            topology=topology,
            projection=projection,
            target_fiber=target_fiber,
            endpoint_specs=refined_endpoints,
            highlight_slots=[source, right],
            plan=merge_plan,
        )
    )
    _ack_merge(runtime, source, right, device=device)
    runtime.validate_remap_acknowledgements()
    runtime.rollback(accepted)
    runtime.commit_irreversible(accepted, dt_phys=0.01, rng_seed=3)
    final_endpoints = [("E1", source, 0.20), ("E2", source, 0.50), ("E3", source, 0.85)]
    states.append(
        _capture_state(
            label="5 · Accepted coarsen",
            subtitle="material interval/rest length restored; new persistent identity",
            topology=topology,
            projection=projection,
            target_fiber=target_fiber,
            endpoint_specs=final_endpoints,
            highlight_slots=[source],
        )
    )
    source_rest_after = float(topology.segment_rest_length_d.numpy()[source])
    output_dir.mkdir(parents=True, exist_ok=True)
    payload = {
        "device": str(wp.get_device(device)),
        "config": {
            "n_fibers": cfg.n_fibers,
            "nodes_per_fiber": cfg.nodes_per_fiber,
            "max_refinement_level": cfg.max_refinement_level,
            "rng_seed": cfg.rng_seed,
        },
        "target_fiber": target_fiber,
        "source_segment": source,
        "reject_max_delta_um": reject_max_delta,
        "source_rest_before_um": source_rest_before,
        "source_rest_after_um": source_rest_after,
        "states": states,
    }
    (output_dir / "ecm_b2_transaction.json").write_text(
        json.dumps(payload, indent=2), encoding="utf-8"
    )
    _render_outputs(states, output_dir)
    _write_report(
        output_dir=output_dir,
        states=states,
        reject_max_delta_um=reject_max_delta,
        source_rest_before_um=source_rest_before,
        source_rest_after_um=source_rest_after,
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--device", default=None, help="Warp CUDA device; no CPU path exists")
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("aleph/outputs/ac_ecm"),
    )
    args = parser.parse_args()
    run(device=str(wp.get_device(args.device)), output_dir=args.output_dir)


if __name__ == "__main__":
    main()
