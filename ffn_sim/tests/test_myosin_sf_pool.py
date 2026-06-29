"""Actin-pool generalization tests for MyosinStepUpdater (KU-3.5 SF NMII).

The cortical updater assumed "actin occupies global tags [0, n_cortex_actin)".
For ventral-stress-fiber NMII the bindable actin is the ``sf_actin`` chain, whose
global tags are NON-contiguous (they sit after the cortex/nucleus/... blocks). The
updater now takes an explicit ``actin_pool_tags`` (default ``arange`` →
byte-identical). These tests assert the indirection is TRANSPARENT:

* **Parity** — a chain at tags ``[0, nb)`` with the default pool, vs the SAME
  geometry shifted to tags ``[pad, pad+nb)`` with ``actin_pool_tags`` passed,
  produce IDENTICAL head→actin binding (the bound global tags differ only by the
  ``pad`` offset; the same heads engage the same chain positions).
* **Confinement** — every head binds a tag INSIDE the pool, never a pad bead or a
  myosin bead.

The binding is driven for real (segment-projection KDTree path), so the generalized
pool / global-tag mapping is exercised — not hand-placed.
"""

from __future__ import annotations

from copy import deepcopy
from pathlib import Path

import numpy as np
import pytest
import yaml

import hoomd

from ffn_sim.archive.hoomd_legacy.cortex.cortex import resolve_h3_derived
from ffn_sim.archive.hoomd_legacy.cortex.myosin import (
    CortexMyosinLayout,
    MyosinStepUpdater,
    cortex_myosin_attach_bin_names,
    resolve_cortex_myosin,
)

CONFIG_PATH = Path(__file__).resolve().parents[1] / "configs" / "phase1_h3.yaml"


def _load_cfg() -> dict:
    with open(CONFIG_PATH) as f:
        return yaml.safe_load(f)


def _demo(cfg: dict) -> dict:
    cfg = deepcopy(cfg)
    cfg["cortex"]["n_filaments"] = 50
    cfg["cortex"]["demo_mode"] = True
    return cfg


def _build_pool_fixture(pad: int, use_pool: bool):
    """Minimal chain + one minifilament; chain at tags [pad, pad+nb).

    Returns (sim, updater, pool_tags, motor_tag_start, nb).
    """
    cfg_full = _load_cfg()
    pc = resolve_h3_derived(_demo(cfg_full))
    ell0 = pc.rest_length
    nb = 7

    cfg = _demo(cfg_full)
    cfg["cortex"]["myosin"]["n_motors_per_cell"] = 1
    cfg["cortex"]["myosin"]["n_heads_per_side"] = 2
    cfg["cortex"]["myosin"]["stepping_mode"] = "binned_r0"   # isolate the pool path
    # Drive de-novo binding deterministically: k_on huge (p_bind≈1), k_off0 tiny.
    cfg["cortex"]["myosin"]["head_actin_k_on"] = 1.0e12
    cfg["cortex"]["myosin"]["head_actin_k_off0"] = 1.0e-12
    p = resolve_cortex_myosin(cfg, dt=pc.dt_cfl)

    N = p.n_backbone
    H = p.n_heads_per_side
    per_motor = p.n_particles_per_motor

    # actin chain along +x at y=0, tags [pad, pad+nb)
    bead_x = np.arange(nb) * ell0
    actin_pos = np.zeros((nb, 3), dtype=np.float64)
    actin_pos[:, 0] = bead_x

    # minifilament backbone along +x centred over the chain, heads ±perp (close)
    perp = 50.0e-9
    bb_center_x = bead_x[nb // 2]
    seg = p.backbone_segment_length
    bb_pos = np.zeros((N, 3), dtype=np.float64)
    bb_pos[:, 0] = bb_center_x + (np.arange(N) - (N - 1) / 2.0) * seg
    # heads distributed along the backbone span at y=±perp.
    head_axis = np.linspace(bead_x[0], bead_x[nb - 1], H)
    plus_heads = np.stack([head_axis, np.full(H, perp), np.zeros(H)], axis=1)
    minus_heads = np.stack([head_axis, np.full(H, -perp), np.zeros(H)], axis=1)
    motor_pos = np.vstack([bb_pos, plus_heads, minus_heads])

    # pad dummy particles (off to the side, far from everything) BEFORE the chain.
    pad_pos = np.zeros((pad, 3), dtype=np.float64)
    if pad > 0:
        pad_pos[:, 0] = -10.0e-6 - np.arange(pad) * 1.0e-6

    all_pos = np.vstack([pad_pos, actin_pos, motor_pos])
    n_total = all_pos.shape[0]
    chain_tag0 = pad
    motor_tag_start = pad + nb

    p_types = ["pad", "cortex_actin", "cortex_myosin_backbone", "cortex_myosin_head"]
    typeid = np.zeros(n_total, dtype=np.uint32)
    typeid[:pad] = 0
    typeid[pad:pad + nb] = 1
    typeid[motor_tag_start:motor_tag_start + N] = 2
    typeid[motor_tag_start + N:] = 3

    attach_names = cortex_myosin_attach_bin_names(p.n_bins)
    bond_types = ["cortex_actin_backbone", "cortex_myosin_backbone",
                  "cortex_myosin_head_backbone"] + attach_names
    bonds: list[tuple[int, int]] = []
    btids: list[int] = []
    seg_groups: list[tuple[int, int]] = []
    for i in range(nb - 1):
        a, b = chain_tag0 + i, chain_tag0 + i + 1
        bonds.append((a, b)); btids.append(0)
        seg_groups.append((a, b))
    for i in range(N - 1):
        bonds.append((motor_tag_start + i, motor_tag_start + i + 1)); btids.append(1)
    for i in range(2 * H):
        bonds.append((motor_tag_start + N + i, motor_tag_start + N // 2)); btids.append(2)

    box_L = 80.0e-6
    snap = hoomd.Snapshot()
    snap.particles.N = n_total
    snap.particles.types = p_types
    snap.particles.typeid[:] = typeid
    snap.particles.position[:] = all_pos
    snap.particles.velocity[:] = np.zeros((n_total, 3))
    snap.particles.mass[:] = np.ones(n_total)
    snap.configuration.box = [box_L, box_L, box_L, 0, 0, 0]
    snap.bonds.N = len(bonds)
    snap.bonds.types = bond_types
    snap.bonds.group[:] = np.asarray(bonds, dtype=np.int64)
    snap.bonds.typeid[:] = np.asarray(btids, dtype=np.uint32)

    dev = hoomd.device.CPU(notice_level=0)
    sim = hoomd.Simulation(device=dev, seed=1)
    sim.create_state_from_snapshot(snap)

    layout = CortexMyosinLayout(
        centers=np.array([[bb_center_x, 0.0, 0.0]]),
        axes=np.array([[1.0, 0.0, 0.0]]),
        positions=motor_pos[None, :, :].copy(),
        motor_tag_start=motor_tag_start,
    )
    pool_tags = np.arange(chain_tag0, chain_tag0 + nb, dtype=np.int64)
    upd = MyosinStepUpdater(
        p_myo=p, layout=layout, kT=pc.kT, n_cortex_actin=nb,
        cortex_bond_groups=np.asarray(seg_groups, dtype=np.int64),
        actin_pool_tags=(pool_tags if use_pool else None),
    )
    upd._sim_ref = sim
    return sim, upd, pool_tags, motor_tag_start, nb


def _bound_tags(upd) -> np.ndarray:
    return upd._head_bound_to_actin.copy()


def test_shifted_pool_parity_with_default():
    """Default pool [0,nb) vs the SAME chain shifted to [pad,pad+nb) bind identically."""
    pad = 13
    sim0, upd0, pool0, mts0, nb = _build_pool_fixture(pad=0, use_pool=False)
    simP, updP, poolP, mtsP, _ = _build_pool_fixture(pad=pad, use_pool=True)

    upd0.act(0)
    updP.act(0)

    b0 = _bound_tags(upd0)
    bP = _bound_tags(updP)
    engaged0 = b0 >= 0
    engagedP = bP >= 0
    # Same heads engaged.
    assert np.array_equal(engaged0, engagedP)
    assert int(engaged0.sum()) > 0, "no heads bound — fixture geometry broke"
    # Bound chain tags differ ONLY by the pad offset (same chain position).
    assert np.array_equal(bP[engagedP], b0[engaged0] + pad)


def test_bound_tags_confined_to_pool():
    """Every bound tag is INSIDE the actin pool (never a pad / myosin bead)."""
    pad = 13
    sim, upd, pool, mts, nb = _build_pool_fixture(pad=pad, use_pool=True)
    upd.act(0)
    bound = upd._head_bound_to_actin
    for t in bound[bound >= 0]:
        assert int(t) in set(pool.tolist()), f"bound tag {t} outside the pool {pool}"
        assert int(t) < mts, "bound tag is a myosin bead, not actin"


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q"]))
