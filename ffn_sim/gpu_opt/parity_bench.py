"""Parity + benchmark harness for the DCM GPU kernel port.

Compares the reference numpy kernels (dcm_kernels_cpu.py) against the GPU cupy kernels
(dcm_kernels_gpu.py — written by the optimizer) on synthetic many-cell DCM arrays, and
times both. The port is ACCEPTED when every kernel matches to rtol≈1e-9 AND the GPU
version is faster on a large system.

Usage:
    python gpu_opt/parity_bench.py            # parity (+ bench if cupy/GPU available)
    python gpu_opt/parity_bench.py --n-cells 200 --n-per-cell 162   # large bench

The GPU module must expose the SAME function names/signatures as dcm_kernels_cpu.py but
accept/return cupy arrays (or numpy-on-GPU). The harness moves inputs to device, calls the
GPU kernel, moves the result back, and compares to the CPU reference.
"""

from __future__ import annotations

import argparse
import time

import numpy as np

from ffn_sim.gpu_opt import dcm_kernels_cpu as K


def _try_gpu():
    try:
        from ffn_sim.gpu_opt import dcm_kernels_gpu as G  # written by the optimizer
        import cupy as cp  # noqa: F401
        return G
    except Exception as e:  # noqa: BLE001
        print(f"[parity] GPU kernels not available yet ({e}); CPU-only run.")
        return None


KERNELS = [
    ("turgor", lambda M, mod: mod.turgor_forces(
        M["pos"], M["faces"], M["face_cell"], M["n_cells"], M["V0"], M["turgor_dP0"], M["K_vol"])),
    ("substrate", lambda M, mod: mod.substrate_forces(
        M["pos"], M["z0"], M["k_well"], M["rng_subst"])),
    ("fa_clutch", lambda M, mod: mod.fa_clutch_forces(
        M["pos"], M["pairs"], M["r0"], M["k_fa"])),
    ("cell_cell_adhesion", lambda M, mod: mod.cell_cell_adhesion_forces(
        M["pos"], M["cell_of_node"], M["sigma"], M["r_cut"], M["k_core"], M["F_adh0"])),
]


def _to_np(x):
    try:
        import cupy as cp
        if isinstance(x, cp.ndarray):
            return cp.asnumpy(x)
    except Exception:  # noqa: BLE001
        pass
    return np.asarray(x)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-cells", type=int, default=40)
    ap.add_argument("--n-per-cell", type=int, default=42)
    ap.add_argument("--reps", type=int, default=20)
    ap.add_argument("--rtol", type=float, default=1e-9)
    args = ap.parse_args()

    M = K.make_synthetic(n_cells=args.n_cells, n_per_cell=args.n_per_cell)
    N = M["pos"].shape[0]
    print(f"[parity] synthetic system: {args.n_cells} cells × {args.n_per_cell} = {N} nodes")
    G = _try_gpu()

    all_ok = True
    for name, call in KERNELS:
        Fcpu = call(M, K)
        t0 = time.perf_counter()
        for _ in range(args.reps):
            call(M, K)
        t_cpu = (time.perf_counter() - t0) / args.reps
        line = f"[{name:20s}] CPU {t_cpu*1e3:8.2f} ms/call"
        if G is not None:
            try:
                Fgpu = _to_np(call(M, G))
                ok = np.allclose(Fcpu, Fgpu, rtol=args.rtol, atol=1e-18)
                # warm + time GPU
                import cupy as cp
                call(M, G); cp.cuda.Stream.null.synchronize()
                t0 = time.perf_counter()
                for _ in range(args.reps):
                    call(M, G)
                cp.cuda.Stream.null.synchronize()
                t_gpu = (time.perf_counter() - t0) / args.reps
                speed = t_cpu / t_gpu if t_gpu > 0 else float("nan")
                maxerr = float(np.max(np.abs(Fcpu - Fgpu)))
                line += f" | GPU {t_gpu*1e3:8.2f} ms ({speed:5.1f}×) | parity {'OK' if ok else 'FAIL'} (maxabs={maxerr:.2e})"
                all_ok = all_ok and ok
            except Exception as e:  # noqa: BLE001
                line += f" | GPU kernel ERROR: {e}"
                all_ok = False
        print(line)

    print(f"\n[parity] {'ALL KERNELS PASS' if (G is not None and all_ok) else ('CPU-only (no GPU kernels)' if G is None else 'PARITY FAILURES — see above')}")
    return 0 if (G is None or all_ok) else 1


if __name__ == "__main__":
    raise SystemExit(main())
