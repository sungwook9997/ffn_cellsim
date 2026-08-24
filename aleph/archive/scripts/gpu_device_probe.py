"""Probe: does Warp come up on this workstation's GPU, and can it run what ffn_cellsim needs?

Answers five questions with evidence, in order of how badly a NO would hurt:
  1. Does `wp.init()` resolve a CUDA device at all under WSL2?  (CLAUDE.md forbids a CPU fallback,
     so a non-CUDA resolution must be reported as a HARD FAIL, never worked around.)
  2. Does a kernel COMPILE for this card's compute capability?  Ada (8.9) has never been targeted
     by this project; every measurement to date is on 8.6.
  3. Does it RUN and produce the right answer?  A kernel that compiles and returns garbage is the
     failure mode a smoke test exists to catch.
  4. Do float64 atomics work, and what do they cost relative to float32?  The inner solve and the
     force assembly are float64, and consumer boards rate-limit FP64 hard — this is the number that
     decides whether these cards are usable for the real workload at all.
  5. Where is the kernel cache?  Mixing 8.9 and 8.6 over one cache directory is a known
     false-positive generator in this repo.

This is a SMOKE PROBE, not a physics measurement: the timings below are microbenchmarks of a
synthetic kernel and may not be quoted as engine performance.
"""

from __future__ import annotations

import argparse
import json
import os
import platform
import time

import numpy as np
import warp as wp

N = 1 << 22          # 4.19 M elements — big enough to be bandwidth/ALU bound, small enough to be quick
REPEATS = 5


@wp.kernel
def spring_accum_f32(pos: wp.array(dtype=wp.vec3), idx: wp.array(dtype=wp.int32),
                     k: wp.float32, out: wp.array(dtype=wp.vec3)):
    """Hookean pair force with an ATOMIC scatter — the shape of our force assembly, in float32."""
    t = wp.tid()
    i = idx[2 * t]
    j = idx[2 * t + 1]
    d = pos[j] - pos[i]
    f = k * d
    wp.atomic_add(out, i, f)
    wp.atomic_sub(out, j, f)


@wp.kernel
def spring_accum_f64(pos: wp.array(dtype=wp.vec3d), idx: wp.array(dtype=wp.int32),
                     k: wp.float64, out: wp.array(dtype=wp.vec3d)):
    """The same force in float64 — what the engine actually runs."""
    t = wp.tid()
    i = idx[2 * t]
    j = idx[2 * t + 1]
    d = pos[j] - pos[i]
    f = k * d
    wp.atomic_add(out, i, f)
    wp.atomic_sub(out, j, f)


def _time_then_measure(kernel, dim, inputs, out, device) -> float:
    """Return best-of-REPEATS seconds, then leave ``out`` holding EXACTLY ONE launch's result.

    The two phases must not share an accumulator.  The first version of this probe timed six
    launches into the same array and then compared it against a one-launch reference, which reports
    a relative error of exactly 5.0 — the kernel being right six times over, read as a failure.  The
    final zero-and-relaunch is therefore load-bearing, not tidiness.
    """
    wp.launch(kernel, dim=dim, inputs=inputs, device=device)   # warmup: compile + first-launch cost
    wp.synchronize_device(device)
    best = float("inf")
    for _ in range(REPEATS):
        t0 = time.perf_counter()
        wp.launch(kernel, dim=dim, inputs=inputs, device=device)
        wp.synchronize_device(device)
        best = min(best, time.perf_counter() - t0)
    out.zero_()
    wp.launch(kernel, dim=dim, inputs=inputs, device=device)
    wp.synchronize_device(device)
    return best


def main() -> int:
    print("=" * 70)
    print("ffn_cellsim GPU probe")
    print("=" * 70)
    print(f"CUDA_VISIBLE_DEVICES = {os.environ.get('CUDA_VISIBLE_DEVICES', '(unset)')!r}")
    print(f"SLURM_JOB_ID         = {os.environ.get('SLURM_JOB_ID', '(none)')!r}")

    wp.init()
    print(f"warp {wp.config.version}")
    print(f"kernel cache         = {wp.config.kernel_cache_dir}")

    devices = wp.get_devices()
    print(f"devices              = {[str(d) for d in devices]}")
    cuda = [d for d in devices if d.is_cuda]
    if not cuda:
        print("\nHARD FAIL: no CUDA device resolved. CLAUDE.md forbids a CPU path; not proceeding.")
        return 2

    dev = cuda[0]
    print(f"\n[1] CUDA device       : {dev}")
    print(f"    name              : {dev.name}")
    print(f"    arch (CC*10)      : {dev.arch}")
    print(f"    total memory      : {dev.total_memory / 1e9:.2f} GB")
    print(f"    free memory       : {dev.free_memory / 1e9:.2f} GB")
    print(f"    is_uva / is_mempool: {dev.is_uva} / {dev.is_mempool_supported}")

    rng = np.random.default_rng(0)
    pos_np = rng.standard_normal((N, 3))
    idx_np = rng.integers(0, N, size=2 * N, dtype=np.int32)

    # ---- float32 -------------------------------------------------------------------------------
    print("\n[2/3] compiling + running float32 atomic scatter ...")
    pos32 = wp.array(pos_np.astype(np.float32), dtype=wp.vec3, device=dev)
    idx = wp.array(idx_np, dtype=wp.int32, device=dev)
    out32 = wp.zeros(N, dtype=wp.vec3, device=dev)
    t32 = _time_then_measure(spring_accum_f32, N, [pos32, idx, wp.float32(2.5), out32], out32, dev)
    got32 = out32.numpy()

    # ---- float64 -------------------------------------------------------------------------------
    print("[4]   compiling + running float64 atomic scatter ...")
    pos64 = wp.array(pos_np, dtype=wp.vec3d, device=dev)
    out64 = wp.zeros(N, dtype=wp.vec3d, device=dev)
    t64 = _time_then_measure(spring_accum_f64, N, [pos64, idx, wp.float64(2.5), out64], out64, dev)
    got64 = out64.numpy()

    # ---- correctness against a host reference ---------------------------------------------------
    ref = np.zeros((N, 3))
    d = pos_np[idx_np[1::2]] - pos_np[idx_np[0::2]]
    np.add.at(ref, idx_np[0::2], 2.5 * d)
    np.subtract.at(ref, idx_np[1::2], 2.5 * d)

    err64 = np.max(np.abs(got64 - ref))
    err32 = np.max(np.abs(got32 - ref))
    scale = np.max(np.abs(ref))
    print("\n[3] CORRECTNESS (vs numpy reference, atomics are order-dependent so exactness is not expected)")
    print(f"    max|err| float64  : {err64:.3e}   (relative {err64 / scale:.3e})")
    print(f"    max|err| float32  : {err32:.3e}   (relative {err32 / scale:.3e})")
    ok64 = err64 / scale < 1e-12
    ok32 = err32 / scale < 1e-4
    print(f"    float64 verdict   : {'PASS' if ok64 else 'FAIL'}")
    print(f"    float32 verdict   : {'PASS' if ok32 else 'FAIL'}")

    print("\n[4] FP64 COST (microbenchmark of THIS kernel only — not an engine number)")
    print(f"    float32 best      : {t32 * 1e3:.3f} ms  ({N / t32 / 1e9:.2f} G-pairs/s)")
    print(f"    float64 best      : {t64 * 1e3:.3f} ms  ({N / t64 / 1e9:.2f} G-pairs/s)")
    print(f"    f64/f32 ratio     : {t64 / t32:.2f}x")

    print(f"\n[5] kernel cache dir  : {wp.config.kernel_cache_dir}")
    import pathlib as _pl
    cache = _pl.Path(wp.config.kernel_cache_dir)
    entries = sorted(p.name for p in cache.iterdir()) if cache.is_dir() else []
    print(f"    entries            : {entries[:12]}")
    for sub in entries:
        d = cache / sub
        if d.is_dir():
            print(f"      {sub}/ -> {sorted(q.name for q in d.iterdir())[:8]}")
    tagged = any(str(dev.arch) in e for e in entries)
    print(f"    arch {dev.arch} appears in a cache path name: {tagged}")
    print("    (if False, a cache shared between sm_89 and sm_86 must be proven safe before mixing cards)")

    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", default="", help="write a JSON record here (the ffn_gpu.py artifact)")
    out_path = ap.parse_args().out

    print("\n" + "=" * 70)
    verdict = "PASS" if (ok64 and ok32) else "FAIL"
    if out_path:
        record = {
            "schema": "gpu-device-probe@1",
            "verdict": verdict,
            "host": platform.node(),
            "cuda_visible_devices": os.environ.get("CUDA_VISIBLE_DEVICES"),
            "slurm_job_id": os.environ.get("SLURM_JOB_ID"),
            "warp_version": wp.config.version,
            "device": {"name": dev.name, "arch": dev.arch,
                       "total_memory_gb": dev.total_memory / 1e9,
                       "mempool": bool(dev.is_mempool_supported)},
            "kernel": {"name": "spring_accum atomic scatter", "n_pairs": N, "repeats": REPEATS},
            "timing_s": {"f32_best": t32, "f64_best": t64, "f64_over_f32": t64 / t32},
            "correctness": {"f64_rel_err": float(err64 / scale), "f32_rel_err": float(err32 / scale),
                            "f64_pass": bool(ok64), "f32_pass": bool(ok32)},
            "kernel_cache_dir": str(wp.config.kernel_cache_dir),
            "NOT_QUOTABLE_AS": "engine performance — synthetic kernel, no population, no gate",
        }
        with open(out_path, "w", encoding="utf-8") as fh:
            json.dump(record, fh, indent=2)
        print(f"wrote record -> {out_path}")
    print(f"PROBE VERDICT: {verdict} — CUDA resolved, kernels compiled for arch {dev.arch}, float64 {'verified' if ok64 else 'WRONG'}")
    print("=" * 70)
    return 0 if verdict == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
