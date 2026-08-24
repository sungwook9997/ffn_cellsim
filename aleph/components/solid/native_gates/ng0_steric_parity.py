"""NG-0 (I2b steric) — Warp-CUDA <-> NumPy device-parity gate (run on the gbook A5000 by the lead).

Certifies that the device hash-grid WCA excluded-volume kernel IS its CPU-verified oracle: it drives
``steric_warp.StericForce`` (the device ``wca_steric_kernel`` over a ``wp.HashGrid``) and the pure-NumPy
brute-force O(N^2) reference ``steric_reference.steric_force_bruteforce`` with the IDENTICAL fixed-seed node
cloud (positions, fiber_ids, sigma, k_ev) and compares the per-node steric force to round-off:

    NG-0.steric  wca_steric_kernel (hash-grid)  vs  steric_force_bruteforce (numpy)

The hash grid is a pure ACCELERATOR — it prunes only pairs beyond ``r_c`` that contribute exactly 0 — so
the accelerated per-node force must equal the brute-force sum to summation round-off. Both paths get the
SAME f64 positions; the grid uses an f32 position copy only to bin/prune. PASS iff the relative residual
max|warp - ref| / max|ref| < 1e-10. To also exercise the Magic-Number-Block grid-invariance statement, the
gate reruns with a second hash-table dimension (the number of hash cells must NOT change the physics).

This is a ROUND-OFF parity check, not a magnitude verdict: sigma / k_ev are test values (any positive value
gives the same verdict). Do NOT tune, loosen the tolerance, or edit the kernel to force a pass — a genuine
NG-0 failure is a real on-device bug the CPU oracle could not catch, and must be reported honestly.

I0-A: constructing ``StericForce`` (device HashGrid + Warp arrays) REQUIRES a CUDA device, so this runner
only executes on the gbook A5000 (never the dev Mac).

Run (from ~/ffn_ac_native on gbook):
    PYTHONPATH=. ~/miniconda3/envs/ffn_sim/bin/python aleph/components/solid/native_gates/ng0_steric_parity.py
"""

from __future__ import annotations

import sys

import numpy as np
import warp as wp

from aleph.components.solid import steric_reference
from aleph.components.solid.steric_warp import StericForce

# ----------------------------------------------------------------------------------------------------
# Fixed test case (deterministic; STRUCTURAL — magnitudes irrelevant to a round-off parity verdict).
# ----------------------------------------------------------------------------------------------------
DEVICE: str | None = None  # resolve the current CUDA device; the hardware contract forbids fixed ordinals
SEED = 20260717
RTOL = 1.0e-10          # round-off pass threshold (relative)

N_NODES = 160           # small enough for a trivial O(N^2) brute force, crowded enough for many overlaps
N_FIBERS = 8            # per-node filament id (same-filament pairs excluded in BOTH paths)
SIGMA = 0.5             # steric diameter [um] (test value; r_c = 2^(1/6) sigma ~= 0.561 um)
K_EV = 100.0            # contact stiffness [pN/um] (Magic-Number-Block scale; test value here)
BOX = 2.5              # cube side [um] -> mean spacing ~0.47 um < r_c => plenty of within-cutoff pairs


def _residual(warp_arr: np.ndarray, ref_arr: np.ndarray) -> tuple[float, float]:
    """Return (max_abs, rel) = (max|warp-ref|, max|warp-ref| / max|ref|) with a zero-ref guard."""
    warp_arr = np.asarray(warp_arr, dtype=np.float64)
    ref_arr = np.asarray(ref_arr, dtype=np.float64)
    max_abs = float(np.max(np.abs(warp_arr - ref_arr)))
    scale = float(np.max(np.abs(ref_arr)))
    rel = max_abs / scale if scale > 0.0 else max_abs
    return max_abs, rel


def _report(name: str, max_abs: float, rel: float) -> bool:
    passed = rel < RTOL
    verdict = "PASS" if passed else "FAIL"
    print(f"NG-0 {name}: max_abs={max_abs:.3e} rel={rel:.3e} -> {verdict}")
    return passed


def _device_force(
    pos_np: np.ndarray, fiber_id: np.ndarray, grid_dim: tuple[int, int, int], device: str,
) -> tuple[np.ndarray, float]:
    """Per-node WCA steric force from the device hash-grid kernel (StericForce._accumulate_pos)."""
    steric = StericForce(fiber_id, SIGMA, K_EV, device, grid_dim=grid_dim)
    pos = wp.array(np.ascontiguousarray(pos_np), dtype=wp.vec3d, device=device)
    out = wp.zeros(pos_np.shape[0], dtype=wp.vec3d, device=device)
    steric._accumulate_pos(pos, out)                                       # grid build + kernel launch
    wp.synchronize()
    return out.numpy(), steric.epsilon


def gate_steric(device: str) -> bool:
    """One WCA steric evaluation: device hash-grid kernel vs the numpy brute-force sum, same input."""
    rng = np.random.default_rng(SEED)
    pos_np = (rng.random((N_NODES, 3)) * BOX).astype(np.float64)           # crowded node cloud
    fiber_id = rng.integers(0, N_FIBERS, size=N_NODES).astype(np.int32)    # per-node filament id

    # --- Warp-CUDA path (default hash-table dim) ---
    f_warp, epsilon = _device_force(pos_np, fiber_id, grid_dim=(64, 64, 64), device=device)

    # --- NumPy reference: brute-force O(N^2) WCA, SAME positions / fiber_id / sigma / epsilon ---
    f_ref = steric_reference.steric_force_bruteforce(
        pos_np, SIGMA, epsilon, fiber_id=fiber_id, exclude_same_fiber=True
    )

    max_abs, rel = _residual(f_warp, f_ref)
    ok = _report("wca_steric_kernel", max_abs, rel)

    # --- grid-invariance: a DIFFERENT hash-table dimension must give the identical physics (accelerator) ---
    f_warp2, _ = _device_force(pos_np, fiber_id, grid_dim=(128, 128, 128), device=device)
    inv_abs, inv_rel = _residual(f_warp2, f_ref)
    ok = _report("wca_steric_kernel[grid_dim=128^3]", inv_abs, inv_rel) and ok

    # Diagnostics: pair count within r_c (excl. same-fiber), momentum closure, |max force|.
    r_c = steric_reference.wca_cutoff(SIGMA)
    d = np.linalg.norm(pos_np[:, None, :] - pos_np[None, :, :], axis=2)
    same = fiber_id[:, None] == fiber_id[None, :]
    npair = int(np.sum((d > 0.0) & (d < r_c) & (~same))) // 2
    net = float(np.linalg.norm(f_ref.sum(axis=0)))
    print(f"    diag: within-r_c cross-fiber pairs={npair}  r_c={r_c:.4f} um  "
          f"epsilon={epsilon:.4e} pN.um  |F|_max={np.max(np.linalg.norm(f_ref, axis=1)):.4e} pN  "
          f"|net F|={net:.3e} pN (Newton-3rd)")
    return ok


def main() -> int:
    wp.init()
    dev = wp.get_device(DEVICE)
    if not dev.is_cuda:
        raise RuntimeError(f"NG-0 steric parity requires CUDA; CPU execution is forbidden (resolved {dev})")
    device = str(dev)
    print(f"NG-0 I2b steric parity gate | device={dev} | warp={wp.config.version} | "
          f"N={N_NODES} fibers={N_FIBERS} sigma={SIGMA} k_ev={K_EV} seed={SEED} rtol={RTOL:.0e}")
    ok = gate_steric(device)
    print(f"NG-0 STERIC OVERALL: {'PASS' if ok else 'FAIL'}")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
