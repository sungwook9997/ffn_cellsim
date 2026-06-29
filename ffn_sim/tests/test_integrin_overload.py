"""H.7 integrin-overload break-before-overflow hardening (KU-2.5 updater).

Diagnosis: ``H7_FA_INTEGRIN_OVERLOAD_FIX_2026-06-07.md``. In the rigid
(constrained) FA build an over-stretched ``integrin_ligand`` bond can be dragged
past the Pereverzev slip range (|F|/F_c ≈ 940), where ``pereverzev_k_off`` raises
``FloatingPointError`` (its EXP_ARG_GUARD) — a crash. The secondary hardening
converts such an over-range bond into a FORCED UNBIND (the physically-correct
outcome: a bond past the slip range is already broken), instead of raising.

This test exercises that path directly: it engages one integrin at a huge
separation (force far above the overflow ceiling) and asserts ``act()`` does NOT
raise, the bond breaks, and the forced-break diagnostic increments — while
leaving the normal O(pN) regime untouched.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
import yaml

import hoomd

from ffn_sim.archive.hoomd_legacy.bridge.fa import build_h4_simulation, resolve_h4
from ffn_sim.validation.pereverzev import EXP_ARG_GUARD

CONFIG_PATH = Path(__file__).resolve().parents[1] / "configs" / "phase1_h4.yaml"


@pytest.fixture(scope="module")
def resolved():
    with open(CONFIG_PATH) as f:
        cfg = yaml.safe_load(f)
    return resolve_h4(cfg)


def _inject_overstretched_bond(sim, action, layouts, *, separation_m: float):
    """Move integrin 0 to ``separation_m`` from its ligand + add an engaged bond.

    Returns the integrin-ligand bond type id. Mutates the sim state in place.
    """
    L = layouts[0]
    int_tag = int(L.integrin_tag_start)
    lig_tag = int(L.ligand_tag)

    read = sim.state.get_snapshot()
    pos = np.asarray(read.particles.position, dtype=np.float64).copy()
    # Snapshot is tag-ordered → row == tag. Place the integrin separation_m away
    # from its ligand along +x (huge stretch → F = k_int·separation ≫ ceiling).
    pos[int_tag] = pos[lig_tag] + np.array([separation_m, 0.0, 0.0])

    bond_types = list(read.bonds.types)
    int_bt = bond_types.index("integrin_ligand")
    old_bg = np.asarray(read.bonds.group, dtype=np.int64)
    old_bt = np.asarray(read.bonds.typeid, dtype=np.uint32)
    new_bg = np.concatenate(
        [old_bg.reshape(-1, 2), np.array([[int_tag, lig_tag]], dtype=np.int64)],
        axis=0,
    )
    new_bt = np.concatenate([old_bt, np.array([int_bt], dtype=np.uint32)])

    write = hoomd.Snapshot()
    write.particles.N = int(read.particles.N)
    write.particles.types = list(read.particles.types)
    write.particles.typeid[:] = np.asarray(read.particles.typeid)
    write.particles.position[:] = pos
    write.particles.velocity[:] = np.asarray(read.particles.velocity)
    write.particles.mass[:] = np.asarray(read.particles.mass)
    write.particles.image[:] = np.asarray(read.particles.image)
    write.configuration.box = list(read.configuration.box)
    write.bonds.N = int(new_bg.shape[0])
    write.bonds.types = bond_types
    write.bonds.group[:] = new_bg.astype(np.uint32)
    write.bonds.typeid[:] = new_bt
    for grp in ("angles", "dihedrals", "impropers"):
        src = getattr(read, grp)
        if int(src.N) > 0:
            dst = getattr(write, grp)
            dst.N = int(src.N)
            dst.types = list(src.types)
            dst.group[:] = np.asarray(src.group)
            dst.typeid[:] = np.asarray(src.typeid)
    sim.state.set_snapshot(write)

    # Mark the integrin engaged in the updater's local state.
    local = action._integrin_tag_to_local[int_tag]
    action._engaged[local] = True
    return int_bt


def _build():
    with open(CONFIG_PATH) as f:
        cfg = yaml.safe_load(f)
    p = resolve_h4(cfg)
    sim, layouts, meta = build_h4_simulation(
        p, with_motors=False, with_baoab=True, with_integrin_updater=True,
    )
    sim.run(0)
    return p, sim, layouts, meta


def test_overstretched_bond_force_breaks_without_raising():
    p, sim, layouts, meta = _build()
    action = meta["integrin_action"]

    # Separation that puts the bond force WELL above the overflow ceiling
    # F_overflow = EXP_ARG_GUARD·min(F_s, F_c). r = 2·F_overflow/k_int.
    F_overflow = EXP_ARG_GUARD * min(p.pereverzev.F_s, p.pereverzev.F_c)
    separation = 2.0 * F_overflow / p.k_int_bare
    _inject_overstretched_bond(sim, action, layouts, separation_m=separation)

    n_forced0 = action.n_forced_break
    n_break0 = action.n_break_total
    # Directly invoke the action — must NOT raise FloatingPointError.
    action.act(0)

    # The over-range bond was force-broken (not crashed).
    assert action.n_forced_break == n_forced0 + 1
    assert action.n_break_total == n_break0 + 1
    # The bond is gone from the topology.
    snap = sim.state.get_snapshot()
    bt = np.asarray(snap.bonds.typeid)
    int_bt = list(snap.bonds.types).index("integrin_ligand")
    assert int(np.count_nonzero(bt == int_bt)) == 0
    # And the integrin is marked unbound again.
    L = layouts[0]
    local = action._integrin_tag_to_local[int(L.integrin_tag_start)]
    assert not action._engaged[local]


def test_overflow_ceiling_matches_pereverzev_guard():
    """The ceiling is exactly EXP_ARG_GUARD·min(F_s,F_c) (no magic number)."""
    p, _, _, _ = _build()
    F_overflow = EXP_ARG_GUARD * min(p.pereverzev.F_s, p.pereverzev.F_c)
    # At the KU-2.18 anchor (F_c = 7 pN) this is 700·7 pN = 4.9 nN — far above
    # the O(pN) catch-bond operating point, so normal dynamics are untouched.
    assert F_overflow == pytest.approx(EXP_ARG_GUARD * p.pereverzev.F_c)
    assert F_overflow > 1.0e-9  # nN-scale, ≫ the pN operating range
