"""Per-step wall-time PROFILER for the GPU-friendly DCM stack (gbook A5000).

Answers the one honest question the PI asked: of the ~30-35 ms/step in the live
GPU DCM sim, how many ms are spent in PYTHON callback / gpu_local context
overhead vs ACTUAL GPU kernel compute vs HOOMD-internal machinery — so we know
whether FUSING the custom forces (Python-bound) helps, or whether a compiled
C++/CUDA HOOMD plugin is the only lever (HOOMD-internal-bound).

How it measures (GPU async-aware):
  * Each ``md.force.Custom`` subclass's ``set_forces`` and the BAOAB ``act`` are
    wrapped at the CLASS level BEFORE the sim is built (HOOMD's
    ``CustomForceCompute`` captures the bound ``self.set_forces`` at attach, so an
    instance-level patch after build would be missed — the class patch is picked
    up because the bound method resolves through the wrapped class method).
  * Each wrapper brackets the call with ``cupy.cuda.runtime.deviceSynchronize()``
    on a GPU device, so the measured time includes kernel COMPLETION not just the
    async launch.
  * TOTAL per-step wall is measured by ``sim.run(1)`` in a loop with a device-sync
    each step; "HOOMD internals" is the residual (total − Σ custom − baoab).
  * An "empty gpu_local context" probe times the per-force context enter/exit +
    asarray + zero-write floor (everything a fused force pays ONCE/step instead of
    once-per-force) — this separates irreducible Python/context overhead from
    kernel compute.

Run on the gbook A5000:
    PYTHONPATH=. ~/miniconda3/envs/ffn_sim/bin/python scripts/dcm_gpu_profile.py \
        --n-cells 200 --steps 600 --warmup 100 --active
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np

import hoomd

import ffn_sim.cell.dcm_gpu_forces as gforces
import ffn_sim.integrator.baoab as baoab_mod
from ffn_sim.cell.dcm_gpu_build import (
    ResolvedGpuDCM,
    build_gpu_dcm_simulation,
)

_OUT = Path("ffn_sim/outputs/h_dcm_gpu_lod/profile")

_ON_GPU = False  # set in main()


def _sync():
    if _ON_GPU:
        import cupy as cp
        cp.cuda.runtime.deviceSynchronize()


class Timer:
    def __init__(self) -> None:
        self.sum: dict[str, float] = {}
        self.cnt: dict[str, int] = {}

    def add(self, label: str, dt: float) -> None:
        self.sum[label] = self.sum.get(label, 0.0) + dt
        self.cnt[label] = self.cnt.get(label, 0) + 1


_TIMER = Timer()
_TIMING_ON = False  # gate so warmup steps aren't counted


def _make_timed_set_forces(orig, label: str):
    def timed(self, timestep):
        if not _TIMING_ON:
            return orig(self, timestep)
        _sync()
        t0 = time.perf_counter()
        try:
            return orig(self, timestep)
        finally:
            _sync()
            _TIMER.add(label, time.perf_counter() - t0)
    return timed


def _make_timed_act(orig, label: str):
    def timed(self, timestep):
        if not _TIMING_ON:
            return orig(self, timestep)
        _sync()
        t0 = time.perf_counter()
        try:
            return orig(self, timestep)
        finally:
            _sync()
            _TIMER.add(label, time.perf_counter() - t0)
    return timed


def _patch_classes() -> None:
    """Wrap each force class's set_forces + BAOAB act at the CLASS level."""
    gforces.DcmTurgorForceGPU.set_forces = _make_timed_set_forces(
        gforces.DcmTurgorForceGPU.set_forces, "set_forces_turgor")
    gforces.DcmTentContactGPU.set_forces = _make_timed_set_forces(
        gforces.DcmTentContactGPU.set_forces, "set_forces_tent")
    gforces.DcmSubstrateForceGPU.set_forces = _make_timed_set_forces(
        gforces.DcmSubstrateForceGPU.set_forces, "set_forces_substrate")
    gforces.DcmActiveRimTractionGPU.set_forces = _make_timed_set_forces(
        gforces.DcmActiveRimTractionGPU.set_forces, "set_forces_active")
    baoab_mod.LeimkuhlerMatthewsBAOAB.act = _make_timed_act(
        baoab_mod.LeimkuhlerMatthewsBAOAB.act, "baoab_act")


def _empty_context_cost(h, n_probe: int) -> None:
    """Per-force gpu_local context enter/exit + asarray + zero-write floor."""
    turgor = h["turgor"]
    d = turgor._dispatch()
    xp = d.xp
    for _ in range(n_probe):
        _sync()
        t0 = time.perf_counter()
        with d.snapshot() as snap:
            tag = xp.asarray(snap.particles.tag)
            pos = xp.asarray(snap.particles.position, dtype=xp.float64)
            perm = xp.argsort(tag)
            _ = pos[perm]
            F = xp.zeros_like(pos)
        with d.force_arrays() as arr:
            arr.force[:] = F
        _sync()
        _TIMER.add("empty_ctx_per_force", time.perf_counter() - t0)


def main() -> None:
    global _ON_GPU, _TIMING_ON
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--n-cells", type=int, default=200)
    ap.add_argument("--steps", type=int, default=600)
    ap.add_argument("--warmup", type=int, default=100)
    ap.add_argument("--active", action="store_true")
    ap.add_argument("--cpu", action="store_true")
    args = ap.parse_args()

    _OUT.mkdir(parents=True, exist_ok=True)
    _patch_classes()  # BEFORE building so the bound methods resolve to wrappers

    p = ResolvedGpuDCM()
    device = hoomd.device.CPU(notice_level=0) if args.cpu else None
    h = build_gpu_dcm_simulation(p, args.n_cells, device=device,
                                 active=args.active)
    sim = h["sim"]
    _ON_GPU = isinstance(sim.device, hoomd.device.GPU)
    print(f"[profile] device={type(sim.device).__name__} N={args.n_cells} "
          f"nodes={sim.state.N_particles} active={args.active}")

    # warm up (timing OFF): triggers lazy DeviceDispatch + JIT + buffer alloc.
    sim.run(args.warmup)
    _sync()

    _TIMING_ON = True
    t_total0 = time.perf_counter()
    for _ in range(args.steps):
        sim.run(1)
        _sync()
    total_wall = time.perf_counter() - t_total0
    _TIMING_ON = False

    _empty_context_cost(h, n_probe=200)

    steps = args.steps
    ms = lambda s: 1e3 * s
    per_step = {k: _TIMER.sum[k] / max(1, _TIMER.cnt[k]) for k in _TIMER.sum}

    custom_labels = [k for k in per_step if k.startswith("set_forces_")]
    sum_custom = sum(per_step[k] for k in custom_labels)
    baoab = per_step.get("baoab_act", 0.0)
    total_per_step = total_wall / steps
    hoomd_internal = total_per_step - sum_custom - baoab
    empty_floor = per_step.get("empty_ctx_per_force", float("nan"))

    rows = []
    for k in sorted(custom_labels):
        rows.append((k, per_step[k]))
    rows.append(("baoab_act", baoab))
    rows.append(("--- SUM custom set_forces", sum_custom))
    rows.append(("--- HOOMD internals (residual)", hoomd_internal))
    rows.append(("=== TOTAL per step", total_per_step))
    rows.append(("(probe) empty gpu_local ctx/force", empty_floor))

    print(f"\n=== PER-STEP BREAKDOWN  ({steps} steps, sync-bracketed, "
          f"{'GPU' if _ON_GPU else 'CPU'}) ===")
    print(f"{'component':45s} {'ms/step':>10s} {'% total':>9s}")
    for label, val in rows:
        pct = 100.0 * val / total_per_step if total_per_step > 0 else 0.0
        print(f"{label:45s} {ms(val):10.4f} {pct:8.1f}%")

    n_forces = len(custom_labels)
    py_overhead_separate = empty_floor * n_forces
    print(f"\n[estimate] empty gpu_local ctx floor   = {ms(empty_floor):.4f} ms/force")
    print(f"[estimate] x {n_forces} custom forces           = "
          f"{ms(py_overhead_separate):.4f} ms/step  "
          f"({100.0 * py_overhead_separate / total_per_step:.1f}% of total)")
    print(f"[estimate] kernel-compute in custom    ~ "
          f"{ms(sum_custom - py_overhead_separate):.4f} ms/step "
          f"(SUM custom - ctx floor x {n_forces})")

    out = dict(
        device=type(sim.device).__name__, on_gpu=_ON_GPU,
        n_cells=args.n_cells, nodes=int(sim.state.N_particles),
        active=args.active, steps=steps, warmup=args.warmup,
        total_wall_s=total_wall, total_ms_per_step=ms(total_per_step),
        per_step_ms={k: ms(v) for k, v in per_step.items()},
        per_step_cnt={k: _TIMER.cnt[k] for k in _TIMER.sum},
        sum_custom_ms=ms(sum_custom), baoab_ms=ms(baoab),
        hoomd_internal_ms=ms(hoomd_internal),
        empty_ctx_floor_ms=ms(empty_floor),
        n_custom_forces=n_forces,
        py_overhead_separate_ms=ms(py_overhead_separate),
        kernel_compute_ms=ms(sum_custom - py_overhead_separate),
    )
    act_tag = "active" if args.active else "passive"
    tag = f"n{args.n_cells}_{act_tag}_{type(sim.device).__name__}"
    fp = _OUT / f"profile_{tag}.json"
    fp.write_text(json.dumps(out, indent=2))
    print(f"\n[json] {fp}")


if __name__ == "__main__":
    main()
