"""Parity + benchmark for the cupy port of kernels_cpu.py.

Compares the NumPy reference kernels (kernels_cpu.py) against the cupy kernels
(kernels_gpu.py — written by the optimizer) on synthetic grouped point-mesh arrays, and
times both. Accept when every kernel matches to rtol≈1e-9 AND the cupy version is faster
on a large system.

    python -m ffn_sim.gpu_opt.bench                              # parity (+bench if cupy)
    python -m ffn_sim.gpu_opt.bench --n-groups 200 --n-per 162  # large-N speedup
"""

from __future__ import annotations

import argparse
import time

import numpy as np

from ffn_sim.gpu_opt import kernels_cpu as K


def _try_gpu():
    try:
        from ffn_sim.gpu_opt import kernels_gpu as G   # written by the optimizer
        import cupy as cp  # noqa: F401
        return G
    except Exception as e:  # noqa: BLE001
        print(f"[bench] cupy kernels not available yet ({e}); CPU-only run.")
        return None


KERNELS = [
    ("mesh_pressure", lambda M, m: m.mesh_pressure_forces(
        M["pos"], M["tris"], M["tri_group"], M["n_groups"], M["V0"], M["p0"], M["K"])),
    ("plane_well", lambda M, m: m.plane_well_forces(M["pos"], M["z0"], M["k_well"], M["w"])),
    ("bond_spring", lambda M, m: m.bond_spring_forces(M["pos"], M["bonds"], M["r0"], M["k_bond"])),
    ("group_pair", lambda M, m: m.group_pair_forces(
        M["pos"], M["group_id"], M["sigma"], M["r_cut"], M["k_core"], M["f_well0"])),
    ("tent_contact", lambda M, m: m.tent_contact_forces(
        M["pos"], M["group_id"], M["r_contact"], M["c_adh"], M["rep_strength"],
        M["adh_strength"], M["patch_area"], M["force_cap"])),
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
    ap.add_argument("--n-groups", type=int, default=40)
    ap.add_argument("--n-per", type=int, default=42)
    ap.add_argument("--reps", type=int, default=20)
    ap.add_argument("--rtol", type=float, default=1e-9)
    args = ap.parse_args()

    M = K.make_test_arrays(n_groups=args.n_groups, n_per_group=args.n_per)
    N = M["pos"].shape[0]
    print(f"[bench] synthetic system: {args.n_groups} groups × {args.n_per} = {N} particles")
    G = _try_gpu()

    all_ok = True
    for name, call in KERNELS:
        Fc = call(M, K)
        t0 = time.perf_counter()
        for _ in range(args.reps):
            call(M, K)
        tc = (time.perf_counter() - t0) / args.reps
        line = f"[{name:14s}] CPU {tc*1e3:8.2f} ms/call"
        if G is not None:
            try:
                Fg = _to_np(call(M, G))
                ok = np.allclose(Fc, Fg, rtol=args.rtol, atol=1e-18)
                import cupy as cp
                call(M, G); cp.cuda.Stream.null.synchronize()
                t0 = time.perf_counter()
                for _ in range(args.reps):
                    call(M, G)
                cp.cuda.Stream.null.synchronize()
                tg = (time.perf_counter() - t0) / args.reps
                line += (f" | GPU {tg*1e3:8.2f} ms ({tc/tg:5.1f}×) | "
                         f"parity {'OK' if ok else 'FAIL'} (maxabs={np.max(np.abs(Fc-Fg)):.2e})")
                all_ok = all_ok and ok
            except Exception as e:  # noqa: BLE001
                line += f" | GPU ERROR: {e}"; all_ok = False
        print(line)

    print(f"\n[bench] {'ALL KERNELS PASS' if (G is not None and all_ok) else ('CPU-only (no cupy kernels yet)' if G is None else 'PARITY FAILURES')}")
    return 0 if (G is None or all_ok) else 1


if __name__ == "__main__":
    raise SystemExit(main())
