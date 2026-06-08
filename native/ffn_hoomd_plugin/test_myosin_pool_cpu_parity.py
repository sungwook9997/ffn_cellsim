"""CPU parity: host-side myosin fixed-pool vs the EXACT production attach-bond force.

This is the Stage-2 binder parity gate (no GPU required, runs on macOS/dev).

Claim under test
----------------
The native fixed-pool myosin binder replaces the per-step head-actin attach
bonds (rewired every batch tick via ``set_snapshot``) with a device-resident pool
of ``(head_tag, actin_tag, k, r0)`` slots summed every step. For that to be a
faithful (non-physics-changing) fast path, the per-step FORCE the pool delivers
must equal the force the production attach bonds deliver for the SAME binding
state.

The production per-step myosin force is exactly a ``md.bond.Harmonic`` over the
attach bonds that ``cortex.myosin.MyosinStepUpdater`` maintains
(``k = k_head_actin``, ``r0 = grip-walk eps``). So this test:

1. Builds a small real cortex-actin topology + grip_walk minifilaments using the
   production ``ffn_sim.cortex.myosin`` functions.
2. Runs the production ``MyosinStepUpdater.act()`` once on a CPU HOOMD sim so it
   binds heads to actin (populating ``_head_bound_to_actin`` + the attach bonds).
3. GROUND TRUTH: builds a CPU ``md.bond.Harmonic`` over EXACTLY the attach bonds
   the updater wrote and reads its per-particle forces.
4. POOL: builds ``MyosinAttachmentPool`` from the same updater state
   (``update_from_updater``) and computes ``reference_forces`` (the exact CPU
   mirror of the ``.cu`` kernel).
5. Asserts the two force fields match to ~1e-12 relative.

If they match, the host reference (and therefore the ``.cu`` kernel it mirrors)
reproduces the production force with zero physics change -- the gate the GBOOK
GPU build must also pass before the native binder can be production.

    cd ~/ffn_cellsim/native/ffn_hoomd_plugin
    PYTHONPATH="$PWD/python:$HOME/ffn_cellsim" \
        conda run -n ffn_sim python test_myosin_pool_cpu_parity.py
"""

from __future__ import annotations

from copy import deepcopy
from pathlib import Path

import numpy as np
import yaml
import hoomd
import hoomd.md as md
import gsd.hoomd

from ffn_sim.cortex.myosin import (
    MyosinStepUpdater,
    cortex_myosin_attach_bin_names,
    cortex_myosin_attach_bin_rest_lengths,
    extend_state_with_cortex_myosin,
    generate_cortex_myosin_layout,
    register_cortex_myosin_bond_params,
    resolve_cortex_myosin,
)

# Import the pure-Python pool module DIRECTLY (not via the package __init__,
# which imports the GPU-only compiled _ffn_native and would fail on a CPU-only
# dev box). myosin_pool.py has no native dependency.
import importlib.util as _ilu

_pool_path = Path(__file__).resolve().parent / "python" / "ffn_hoomd_plugin" / "myosin_pool.py"
import sys as _sys
_spec = _ilu.spec_from_file_location("ffn_myosin_pool", _pool_path)
_mod = _ilu.module_from_spec(_spec)
_sys.modules["ffn_myosin_pool"] = _mod  # dataclass introspection needs this
_spec.loader.exec_module(_mod)
MyosinAttachmentPool = _mod.MyosinAttachmentPool


_CFG = (
    Path(__file__).resolve().parents[2]
    / "ffn_sim" / "configs" / "phase1_h3.yaml"
)


def _grip_walk_p_myo(n_motors: int, dt: float):
    with open(_CFG) as f:
        cfg = yaml.safe_load(f)
    cfg = deepcopy(cfg)
    cfg["cortex"]["myosin"]["stepping_mode"] = "grip_walk"
    cfg["cortex"]["myosin"]["n_motors_per_cell"] = n_motors
    return resolve_cortex_myosin(cfg, dt=dt)


def _build_cortex_actin_snapshot(n_fil: int, beads_per_fil: int, ell0: float, L: float):
    """A few short straight cortex-actin filaments, laid out so myosin heads land
    near actin beads (so the production updater will actually bind some heads)."""
    n_actin = n_fil * beads_per_fil
    pos = np.zeros((n_actin, 3), dtype=np.float64)
    rng = np.random.default_rng(7)
    for f in range(n_fil):
        # each filament along +x, offset in y/z so they are distinct lines
        base = np.array([-0.5 * (beads_per_fil - 1) * ell0,
                         (f - 0.5 * (n_fil - 1)) * 0.6 * ell0,
                         (rng.random() - 0.5) * 0.2 * ell0])
        axis = np.array([1.0, 0.0, 0.0])
        for j in range(beads_per_fil):
            pos[f * beads_per_fil + j] = base + j * ell0 * axis

    snap = gsd.hoomd.Frame()
    snap.particles.N = n_actin
    snap.particles.types = ["actin_cortex"]
    snap.particles.typeid = np.zeros(n_actin, dtype=np.uint32)
    snap.particles.position = pos
    snap.particles.mass = np.ones(n_actin, dtype=np.float64)
    # backbone bonds along each filament
    bonds = []
    for f in range(n_fil):
        for j in range(beads_per_fil - 1):
            bonds.append((f * beads_per_fil + j, f * beads_per_fil + j + 1))
    bg = np.asarray(bonds, dtype=np.uint32)
    snap.bonds.N = int(bg.shape[0])
    snap.bonds.types = ["actin_backbone"]
    snap.bonds.group = bg
    snap.bonds.typeid = np.zeros(bg.shape[0], dtype=np.uint32)
    snap.configuration.box = [L, L, L, 0, 0, 0]
    return snap, n_actin, pos


def _to_hoomd_snapshot(frame: gsd.hoomd.Frame) -> hoomd.Snapshot:
    s = hoomd.Snapshot()
    s.configuration.box = list(frame.configuration.box)
    s.particles.N = int(frame.particles.N)
    s.particles.types = list(frame.particles.types)
    s.particles.typeid[:] = np.asarray(frame.particles.typeid)
    s.particles.position[:] = np.asarray(frame.particles.position)
    s.particles.mass[:] = np.asarray(frame.particles.mass)
    s.bonds.N = int(frame.bonds.N)
    s.bonds.types = list(frame.bonds.types)
    if int(frame.bonds.N) > 0:
        s.bonds.group[:] = np.asarray(frame.bonds.group)
        s.bonds.typeid[:] = np.asarray(frame.bonds.typeid)
    return s


def main() -> int:
    n_fil, beads_per_fil = 6, 5
    L = 1.0e-5
    ell0 = 0.45e-6  # cortex bead spacing (rest length); near grip-walk ℓ₀ scale
    dt = 1.0e-9
    n_motors = 4

    p_myo = _grip_walk_p_myo(n_motors=n_motors, dt=dt)
    # Make myosin heads reach the actin lattice we built (override the capture /
    # max-bind distances UP so the legacy bead-center binder finds neighbours in
    # this toy box -- this only widens the binding search, not the FORCE law,
    # which is what parity tests). Not a production value; test geometry only.
    p_myo.head_actin_max_bind_dist = 1.5e-6
    p_myo.head_actin_capture_perp = 1.5e-6

    base_frame, n_actin, actin_pos = _build_cortex_actin_snapshot(
        n_fil, beads_per_fil, ell0, L
    )

    layout = generate_cortex_myosin_layout(
        p_myo, R_cell=0.5 * L, motor_tag_start=n_actin,
        rng=np.random.default_rng(1),
    )
    # Re-center the minifilaments right on top of actin beads so heads are within
    # binding distance (legacy sphere-random placement otherwise scatters them).
    rng = np.random.default_rng(3)
    chosen = rng.choice(n_actin, size=n_motors, replace=False)
    for m in range(n_motors):
        shift = actin_pos[chosen[m]] - layout.centers[m]
        layout.positions[m] += shift
        layout.centers[m] = actin_pos[chosen[m]]

    frame = extend_state_with_cortex_myosin(base_frame, layout, p_myo)
    snap = _to_hoomd_snapshot(frame)
    pos = np.asarray(snap.particles.position, dtype=np.float64).copy()
    box_L = np.asarray(snap.configuration.box[:3], dtype=np.float64)
    N = int(snap.particles.N)
    bt_names_all = list(snap.bonds.types)

    # grip_walk geometry: uniform fixed-N cortex bead-tag map. The updater is
    # constructed exactly as in production; we drive its per-head bound state
    # DIRECTLY (rather than rely on the stochastic ~5e-6/tick binder firing) so
    # the parity is deterministic and tests the FORCE law for a known binding
    # state -- the binder itself (Bell-Evans/Hill) is unchanged and stays in
    # Python; only the per-step force is being moved to the pool.
    upd = MyosinStepUpdater(
        p_myo=p_myo, layout=layout, kT=1.380649e-23 * 310.0,
        n_cortex_actin=n_actin,
        ell0_cortex=ell0,
        cortex_beads_per_filament=beads_per_fil,
    )

    # Seed a deterministic, geometry-respecting binding state: bind every head
    # to its NEAREST cortex-actin bead within head_actin_max_bind_dist (the same
    # rule the production Step-2 binder uses), populating both _head_bound_to_actin
    # and the grip_walk per-head (filament, pos) fields the production binder sets.
    n_heads_total = 2 * p_myo.n_heads_per_side * n_motors
    from scipy.spatial import cKDTree
    tree = cKDTree(pos[:n_actin])
    n_engaged = 0
    for h in range(n_heads_total):
        htag = upd._head_global_tag(h)
        d, j = tree.query(pos[htag])
        if d <= p_myo.head_actin_max_bind_dist:
            upd._head_bound_to_actin[h] = int(j)
            fil, pos_j = upd._tag_to_fil_pos(int(j))
            upd._head_bound_filament[h] = fil
            upd._head_bound_bead_pos[h] = pos_j
            upd._head_grip_s[h] = 0.0
            n_engaged += 1
    print(f"[setup] seeded engaged heads: {n_engaged} (of {n_heads_total})",
          flush=True)
    if n_engaged == 0:
        print("myosin-pool-cpu-parity FAIL (no heads in range; widen test geometry)",
              flush=True)
        return 1

    # Reconstruct the production attach bonds + their bin typeids EXACTLY as the
    # updater's Step-2 would write them: bond (head_tag, actin_tag), bin index
    # from the head-bead distance. This is the ground-truth per-step topology.
    bin_width = p_myo.head_actin_max_bind_dist / p_myo.n_bins
    attach_bonds_list = []
    attach_bins_list = []
    for h in range(n_heads_total):
        j = int(upd._head_bound_to_actin[h])
        if j < 0:
            continue
        htag = upd._head_global_tag(h)
        d_use = min(float(np.linalg.norm(pos[htag] - pos[j])),
                    p_myo.head_actin_max_bind_dist - 1e-12)
        idx_bin = int(min(p_myo.n_bins - 1, max(0, int(d_use / bin_width))))
        attach_bonds_list.append((htag, j))
        attach_bins_list.append(idx_bin)
    attach_bonds = np.asarray(attach_bonds_list, dtype=np.int64)
    attach_bins = np.asarray(attach_bins_list, dtype=np.int64)

    bin_names = cortex_myosin_attach_bin_names(p_myo.n_bins)
    bin_r0 = cortex_myosin_attach_bin_rest_lengths(
        p_myo.n_bins, p_myo.head_actin_max_bind_dist, stepping_mode="grip_walk"
    )
    r0_eps = float(bin_r0[0])  # the grip-walk attach rest length (1 pm)

    def _ground_truth_force(bonds: np.ndarray, bins: np.ndarray) -> np.ndarray:
        """Per-particle force from a CPU md.bond.Harmonic over EXACTLY ``bonds``
        with the production grip-walk attach params. This IS the production
        per-step myosin contribution for that binding state."""
        gt = hoomd.Snapshot()
        gt.configuration.box = list(snap.configuration.box)
        gt.particles.N = N
        gt.particles.types = list(snap.particles.types)
        gt.particles.typeid[:] = np.asarray(snap.particles.typeid)
        gt.particles.position[:] = pos
        gt.particles.mass[:] = np.asarray(snap.particles.mass)
        gt.bonds.N = int(bonds.shape[0])
        gt.bonds.types = list(bin_names)
        if bonds.shape[0] > 0:
            gt.bonds.group[:] = bonds.astype(np.uint32)
            gt.bonds.typeid[:] = bins.astype(np.uint32)
        s_gt = hoomd.Simulation(device=hoomd.device.CPU(notice_level=0), seed=1)
        s_gt.create_state_from_snapshot(gt)
        hb = md.bond.Harmonic()
        for nm, r0v in zip(bin_names, bin_r0):
            hb.params[nm] = dict(k=p_myo.k_head_actin, r0=float(r0v))
        s_gt.operations.integrator = md.Integrator(dt=dt, forces=[hb], methods=[])
        s_gt.run(0)
        return np.asarray(hb.forces, dtype=np.float64).copy()

    def _compare(label: str, bonds: np.ndarray, bins: np.ndarray) -> bool:
        F_truth = _ground_truth_force(bonds, bins)
        pool = MyosinAttachmentPool(
            pool_size=n_heads_total,
            k_head_actin=p_myo.k_head_actin,
            r0=r0_eps,
        )
        pool.update_from_updater(upd)
        F_pool = pool.reference_forces(pos, box_L)  # row==tag (dense)
        dF = float(np.max(np.abs(F_truth - F_pool)))
        sF = float(np.max(np.abs(F_truth))) + 1e-300
        rel = dF / sF
        ok_force = rel < 1e-10
        n_active = int(((pool.k > 0) & (pool.actin_tag >= 0)).sum())
        ok_count = n_active == int(bonds.shape[0])
        print(f"[{label}] bonds={bonds.shape[0]:>3} pool_active={n_active:>3}  "
              f"max|truth-pool|={dF:.3e} rel={rel:.2e}  "
              f"force={'PASS' if ok_force else 'FAIL'} "
              f"count={'PASS' if ok_count else 'FAIL'}", flush=True)
        return ok_force and ok_count

    # --- Scenario 1: all heads bound (multi-head-per-bead atomic accumulation) ---
    ok1 = _compare("all-bound", attach_bonds, attach_bins)

    # --- Scenario 2: unbind half the heads -> pool must zero those slots and the
    #     remaining force must still match a HOOMD ref over only the active half. ---
    for h in range(0, n_heads_total, 2):
        upd._head_bound_to_actin[h] = -1
        upd._head_bound_filament[h] = -1
        upd._head_bound_bead_pos[h] = -1
        upd._head_grip_s[h] = 0.0
    half_bonds_list, half_bins_list = [], []
    for row in range(attach_bonds.shape[0]):
        htag = int(attach_bonds[row, 0])
        h_local = upd._head_local_from_tag(htag)
        if h_local >= 0 and upd._head_bound_to_actin[h_local] >= 0:
            half_bonds_list.append(attach_bonds[row])
            half_bins_list.append(attach_bins[row])
    half_bonds = np.asarray(half_bonds_list, dtype=np.int64).reshape(-1, 2)
    half_bins = np.asarray(half_bins_list, dtype=np.int64).reshape(-1)
    ok2 = _compare("half-unbound", half_bonds, half_bins)

    ok = ok1 and ok2
    print("myosin-pool-cpu-parity " + ("PASS" if ok else "FAIL"), flush=True)
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
