"""Parity test for kernels_gpu.py on a CPU-only box, via a numpy-backed cupy shim.

Installs fake 'cupy'/'cupyx' modules (backed by numpy) into sys.modules BEFORE
importing kernels_gpu, then checks all 4 kernels against kernels_cpu on small and
large synthetic systems plus edge cases. Exits nonzero on any failure.
"""

from __future__ import annotations

import sys
import time
import types

import numpy as np

# --------------------------------------------------------------------------
# numpy-backed cupy/cupyx shim (must be installed before importing kernels_gpu)
# --------------------------------------------------------------------------
fake_cp = types.ModuleType("cupy")
for _k in dir(np):
    if not _k.startswith("_"):
        setattr(fake_cp, _k, getattr(np, _k))
fake_cp.asnumpy = np.asarray


def _scatter_add(a, idx, val):
    np.add.at(a, idx, val)


fake_cupyx = types.ModuleType("cupyx")
fake_cupyx.scatter_add = _scatter_add
fake_cp.scatter_add = _scatter_add  # in case of `from cupy import ...` style use

sys.modules["cupy"] = fake_cp
sys.modules["cupyx"] = fake_cupyx

import kernels_cpu as KC  # noqa: E402
import kernels_gpu as KG  # noqa: E402

RTOL, ATOL = 1e-9, 1e-18
results = []  # (name, ok, maxabs, note)


def check(name, Fc, Fg, note=""):
    Fg = np.asarray(Fg)
    maxabs = float(np.max(np.abs(Fc - Fg))) if Fc.size else 0.0
    ok = np.allclose(Fc, Fg, rtol=RTOL, atol=ATOL) and Fc.shape == Fg.shape
    results.append((name, ok, maxabs, note))


def run_case(tag, M, do_k4=True, k4_note=""):
    Fc = KC.mesh_pressure_forces(M["pos"], M["tris"], M["tri_group"], M["n_groups"],
                                 M["V0"], M["p0"], M["K"])
    Fg = KG.mesh_pressure_forces(M["pos"], M["tris"], M["tri_group"], M["n_groups"],
                                 M["V0"], M["p0"], M["K"])
    check(f"K1 mesh_pressure {tag}", Fc, Fg)

    Fc = KC.plane_well_forces(M["pos"], M["z0"], M["k_well"], M["w"])
    Fg = KG.plane_well_forces(M["pos"], M["z0"], M["k_well"], M["w"])
    check(f"K2 plane_well    {tag}", Fc, Fg)

    Fc = KC.bond_spring_forces(M["pos"], M["bonds"], M["r0"], M["k_bond"])
    Fg = KG.bond_spring_forces(M["pos"], M["bonds"], M["r0"], M["k_bond"])
    check(f"K3 bond_spring   {tag}", Fc, Fg)

    if do_k4:
        t0 = time.perf_counter()
        Fc = KC.group_pair_forces(M["pos"], M["group_id"], M["sigma"], M["r_cut"],
                                  M["k_core"], M["f_well0"])
        t_ref = time.perf_counter() - t0
        Fg = KG.group_pair_forces(M["pos"], M["group_id"], M["sigma"], M["r_cut"],
                                  M["k_core"], M["f_well0"])
        check(f"K4 group_pair   {tag}", Fc, Fg, f"cpu ref {t_ref:.1f}s{k4_note}")
    else:
        results.append((f"K4 group_pair   {tag}", True, float("nan"),
                        "SKIPPED (CPU O(N^2) ref too slow; covered by (a)+(c))"))


# (a) small system ----------------------------------------------------------
Ma = KC.make_test_arrays(n_groups=40, n_per_group=42)
run_case("(a) 40x42  ", Ma)

# edge cases on (a): masked plane well, empty bonds --------------------------
mask = (np.arange(Ma["pos"].shape[0]) % 3) == 0
Fc = KC.plane_well_forces(Ma["pos"], Ma["z0"], Ma["k_well"], Ma["w"], mask=mask)
Fg = KG.plane_well_forces(Ma["pos"], Ma["z0"], Ma["k_well"], Ma["w"], mask=mask)
check("K2 plane_well   mask      ", Fc, Fg)

empty = np.zeros((0, 2), dtype=np.int64)
Fc = KC.bond_spring_forces(Ma["pos"], empty, np.zeros(0), Ma["k_bond"])
Fg = KG.bond_spring_forces(Ma["pos"], empty, np.zeros(0), Ma["k_bond"])
check("K3 bond_spring  empty     ", Fc, Fg)

# (c) K4 randomized dense stress case ---------------------------------------
rng = np.random.default_rng(42)
Nc = 2000
sigma, r_cut = Ma["sigma"], Ma["r_cut"]
pos_c = rng.uniform(0.0, 4.0 * r_cut, size=(Nc, 3))
gid_c = rng.integers(-1, 8, size=Nc).astype(np.int64)  # includes some -1
Fc = KC.group_pair_forces(pos_c, gid_c, sigma, r_cut, Ma["k_core"], Ma["f_well0"])
Fg = KG.group_pair_forces(pos_c, gid_c, sigma, r_cut, Ma["k_core"], Ma["f_well0"])
check("K4 group_pair   (c) stress", Fc, Fg, "2000 pts dense box, gid incl -1")

# (b) large system -----------------------------------------------------------
Mb = KC.make_test_arrays(n_groups=200, n_per_group=162)
# Estimate K4 CPU-reference cost from case (a)'s O(N^2) scaling before committing.
Na, Nb = Ma["pos"].shape[0], Mb["pos"].shape[0]
t0 = time.perf_counter()
KC.group_pair_forces(Ma["pos"], Ma["group_id"], Ma["sigma"], Ma["r_cut"],
                     Ma["k_core"], Ma["f_well0"])
ta = time.perf_counter() - t0
est_b = ta * (Nb / Na) ** 2
run_case("(b) 200x162", Mb, do_k4=(est_b <= 120.0),
         k4_note=f" (est {est_b:.0f}s)")

# ---------------------------------------------------------------------------
print(f"\n{'kernel / case':<28s} {'status':<6s} {'max abs err':<12s} note")
print("-" * 78)
all_ok = True
for name, ok, maxabs, note in results:
    all_ok &= ok
    ma = "-" if maxabs != maxabs else f"{maxabs:.2e}"
    print(f"{name:<28s} {'PASS' if ok else 'FAIL':<6s} {ma:<12s} {note}")
print("-" * 78)
print(f"OVERALL: {'PASS' if all_ok else 'FAIL'}")
sys.exit(0 if all_ok else 1)
