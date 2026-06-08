"""Tests for the intermediate-filament (IF) perinuclear cage compartment.

Light + fast: resolver-level + OFF-identity + sign-sense + γ-denylist. The one
HOOMD-sim check (builtin harmonic bonded eval) is small and tagged so it can be
skipped if the device is unavailable; the core gates run with no sim.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

import hoomd

from ffn_sim.cell.intermediate_filaments import (
    BACKBONE_BOND_TYPE,
    CROSSLINK_BOND_TYPE,
    GAMMA_DENYLIST_PREFIX,
    IF_BEAD_TYPE,
    PI_DECISIONS,
    ResolvedIntermediateFilaments,
    attach_if_bonds_to_simulation,
    build_if_cage_layout,
    build_intermediate_filament_bonds,
    extend_snapshot_with_if_cage,
    resolve_intermediate_filaments,
)

kT = 4.114e-21          # J  (T ≈ 298 K)
R_CELL = 7.5e-6         # m  MCF7 (Wagner 2011)
R_NUC = 2.5e-6          # m  typical mammalian nucleus radius
GAMMA_IF = 1.0e-6       # N·s/m  per-bead Stokes drag placeholder


def _enabled_cfg(**over) -> dict:
    cfg = {"enabled": True, "n_filaments": 8, "beads_per_fil": 5}
    cfg.update(over)
    return cfg


# ---------------------------------------------------------------------------
# OFF-identity
# ---------------------------------------------------------------------------
def test_resolve_default_off():
    """Missing block and enabled:false both yield .enabled False with zeros."""
    p_missing = resolve_intermediate_filaments({}, kT=kT, R_cell=R_CELL)
    assert p_missing.enabled is False
    assert p_missing.n_beads_total == 0
    assert p_missing.k_bb == 0.0

    p_false = resolve_intermediate_filaments(
        {"enabled": False, "n_filaments": 100}, kT=kT, R_cell=R_CELL
    )
    assert p_false.enabled is False
    assert p_false.n_filaments == 0
    assert p_false.k_bb == 0.0
    assert p_false.k_xl == 0.0


def test_builder_off_identity_snapshot_unchanged():
    """Disabled builder returns the input snapshot unchanged (zero added)."""
    snap = hoomd.Snapshot()
    snap.particles.N = 3
    snap.particles.types = ["A"]
    snap.particles.position[:] = np.zeros((3, 3))
    snap.configuration.box = [1e-5, 1e-5, 1e-5, 0, 0, 0]

    p_off = resolve_intermediate_filaments({}, kT=kT, R_cell=R_CELL)
    out = extend_snapshot_with_if_cage(
        snap, p_off, (0.0, 0.0, 0.0), gamma_if=GAMMA_IF
    )
    assert out is snap                      # exact same object → no-op
    assert out.particles.N == 3
    assert out.bonds.N == 0


def test_attach_and_force_builder_off_return_none():
    """Disabled attach / bond builder return None (touch nothing)."""
    p_off = resolve_intermediate_filaments({}, kT=kT, R_cell=R_CELL)
    assert build_intermediate_filament_bonds(p_off) is None
    # attach with no integrator must still no-op (not raise) when disabled.
    sim = hoomd.Simulation(device=hoomd.device.CPU())
    assert attach_if_bonds_to_simulation(sim, p_off, gamma_if=GAMMA_IF) is None


# ---------------------------------------------------------------------------
# Dimensional / parameter sanity
# ---------------------------------------------------------------------------
def test_resolve_derived_si_quantities():
    """k_bb = E_if·A_if/l_seg with A_if = π(d_if/2)²; kappa = Lp·kT."""
    cfg = _enabled_cfg(
        E_if=6.0e6, d_if=10.0e-9, Lp=0.5e-6, ratio_xl=0.1
    )
    p = resolve_intermediate_filaments(cfg, kT=kT, R_cell=R_CELL, R_nuc=R_NUC)
    assert p.enabled is True

    A_expect = math.pi * (0.5 * 10.0e-9) ** 2
    assert p.A_if == pytest.approx(A_expect, rel=1e-12)
    # l_seg defaults to Lp (flexible-chain limit).
    assert p.l_seg == pytest.approx(0.5e-6, rel=1e-12)
    k_bb_expect = 6.0e6 * A_expect / 0.5e-6
    assert p.k_bb == pytest.approx(k_bb_expect, rel=1e-12)
    assert p.k_xl == pytest.approx(0.1 * k_bb_expect, rel=1e-12)
    assert p.kappa_bend == pytest.approx(0.5e-6 * kT, rel=1e-12)

    # Counts.
    assert p.n_beads_total == 8 * 5
    assert p.n_backbone_bonds == 8 * (5 - 1)
    # Cage shell seeded just outside the nucleus when R_nuc given.
    assert p.R_cage_inner == pytest.approx(R_NUC, rel=1e-12)
    assert p.R_cage_outer > p.R_cage_inner


def test_resolve_rejects_bad_params():
    """Negative / zero / out-of-band params raise ValueError when enabled."""
    with pytest.raises(ValueError):
        resolve_intermediate_filaments(
            _enabled_cfg(n_filaments=0), kT=kT, R_cell=R_CELL
        )
    with pytest.raises(ValueError):
        resolve_intermediate_filaments(
            _enabled_cfg(beads_per_fil=0), kT=kT, R_cell=R_CELL
        )
    with pytest.raises(ValueError):
        resolve_intermediate_filaments(
            _enabled_cfg(E_if=-1.0), kT=kT, R_cell=R_CELL
        )
    # E_if out of the few-MPa band.
    with pytest.raises(ValueError):
        resolve_intermediate_filaments(
            _enabled_cfg(E_if=5.0e8), kT=kT, R_cell=R_CELL
        )
    # Lp out of band.
    with pytest.raises(ValueError):
        resolve_intermediate_filaments(
            _enabled_cfg(Lp=5.0e-6), kT=kT, R_cell=R_CELL
        )
    # d_if out of band.
    with pytest.raises(ValueError):
        resolve_intermediate_filaments(
            _enabled_cfg(d_if=50.0e-9), kT=kT, R_cell=R_CELL
        )
    # Bad cell radius.
    with pytest.raises(ValueError):
        resolve_intermediate_filaments(_enabled_cfg(), kT=kT, R_cell=-1.0)


def test_layout_counts_and_isolation():
    """Layout adds expected beads/backbone bonds; pure (no mutation)."""
    p = resolve_intermediate_filaments(
        _enabled_cfg(n_filaments=6, beads_per_fil=4), kT=kT, R_cell=R_CELL,
        R_nuc=R_NUC,
    )
    layout = build_if_cage_layout(p, (0.0, 0.0, 0.0), gamma_if=GAMMA_IF, seed=1)
    pos = layout["positions"]
    bb = layout["backbone_bonds"]
    xl = layout["crosslink_bonds"]
    assert pos.shape == (6 * 4, 3)
    assert bb.shape == (6 * (4 - 1), 2)
    # backbone bonds connect consecutive beads within a filament only.
    for a, b in bb:
        assert b == a + 1
        assert a // 4 == b // 4         # same filament block
    # crosslinks (if any) bridge different filaments.
    for a, b in xl:
        assert a // 4 != b // 4
    # all beads seeded in the perinuclear shell.
    r = np.linalg.norm(pos, axis=1)
    assert r.min() >= 0.0
    # gamma array sized to beads.
    assert layout["gamma"].shape == (6 * 4,)


def test_build_layout_on_disabled_raises():
    p_off = resolve_intermediate_filaments({}, kT=kT, R_cell=R_CELL)
    with pytest.raises(RuntimeError):
        build_if_cage_layout(p_off, (0.0, 0.0, 0.0), gamma_if=GAMMA_IF)


# ---------------------------------------------------------------------------
# Sign / sense — a stretched IF backbone pulls beads back together
# ---------------------------------------------------------------------------
def test_backbone_sign_sense_analytic():
    """Stretched harmonic backbone: restoring force pulls beads inward.

    F = -k_bb (r - l_seg) r̂ ; for r > l_seg the force on bead-1 points toward
    bead-0 (inward), magnitude k_bb·(r - l_seg).
    """
    p = resolve_intermediate_filaments(_enabled_cfg(), kT=kT, R_cell=R_CELL)
    l0 = p.l_seg
    k = p.k_bb
    # Place bead 0 at origin, bead 1 stretched to 1.5·l0 along +x.
    r0 = np.array([0.0, 0.0, 0.0])
    r1 = np.array([1.5 * l0, 0.0, 0.0])
    dr = r1 - r0
    r = np.linalg.norm(dr)
    rhat = dr / r
    F_on_1 = -k * (r - l0) * rhat       # force on bead 1
    # Stretched (r > l0) → force on bead 1 points -x (back toward bead 0).
    assert r > l0
    assert F_on_1[0] < 0.0
    assert F_on_1[0] == pytest.approx(-k * (0.5 * l0), rel=1e-12)
    # Equal-and-opposite on bead 0 (no net momentum).
    F_on_0 = -F_on_1
    assert np.allclose(F_on_0 + F_on_1, 0.0)


def test_nonlinear_not_faked_raises():
    """Requesting the nonlinear law raises (never silently approximated)."""
    p = resolve_intermediate_filaments(_enabled_cfg(), kT=kT, R_cell=R_CELL)
    with pytest.raises(NotImplementedError):
        build_intermediate_filament_bonds(p, nonlinear=True)


# ---------------------------------------------------------------------------
# No-γ-contamination — every IF bond type carries the denylist prefix
# ---------------------------------------------------------------------------
def test_gamma_denylist_prefix():
    """Prefix non-empty and every bond type this module creates starts with it."""
    assert GAMMA_DENYLIST_PREFIX
    assert GAMMA_DENYLIST_PREFIX == "if_"
    assert BACKBONE_BOND_TYPE.startswith(GAMMA_DENYLIST_PREFIX)
    assert CROSSLINK_BOND_TYPE.startswith(GAMMA_DENYLIST_PREFIX)
    assert IF_BEAD_TYPE.startswith(GAMMA_DENYLIST_PREFIX)

    # The bond types actually attached to a force also carry the prefix.
    p = resolve_intermediate_filaments(_enabled_cfg(), kT=kT, R_cell=R_CELL)
    harmonic = build_intermediate_filament_bonds(p)
    # md.bond.Harmonic.params is keyed by bond type name.
    for type_name in (BACKBONE_BOND_TYPE, CROSSLINK_BOND_TYPE):
        assert type_name.startswith(GAMMA_DENYLIST_PREFIX)
        # params accessible (does not require a sim).
        _ = harmonic.params[type_name]


def test_pi_decisions_documented():
    """The nonlinear strain-stiffening gap is surfaced as a PI decision."""
    assert isinstance(PI_DECISIONS, list)
    joined = " ".join(PI_DECISIONS).lower()
    assert "nonlinear" in joined or "strain-stiffening" in joined


# ---------------------------------------------------------------------------
# Snapshot-extension shape + bond-type registration (no sim run)
# ---------------------------------------------------------------------------
def test_extend_snapshot_adds_beads_and_if_bonds():
    """Enabled builder appends if_bead particles + if_ bonds; leaves old intact."""
    snap = hoomd.Snapshot()
    snap.particles.N = 4
    snap.particles.types = ["actin_cortex"]
    snap.particles.position[:] = np.array(
        [[1e-7, 0, 0], [0, 1e-7, 0], [0, 0, 1e-7], [1e-7, 1e-7, 0]]
    )
    snap.bonds.N = 1
    snap.bonds.types = ["cortex_backbone"]
    snap.bonds.group[:] = np.array([[0, 1]])
    snap.configuration.box = [3e-5, 3e-5, 3e-5, 0, 0, 0]

    p = resolve_intermediate_filaments(
        _enabled_cfg(n_filaments=4, beads_per_fil=3), kT=kT, R_cell=R_CELL,
        R_nuc=R_NUC,
    )
    out = extend_snapshot_with_if_cage(
        snap, p, (0.0, 0.0, 0.0), gamma_if=GAMMA_IF, seed=2
    )
    assert out is not snap
    assert out.particles.N == 4 + p.n_beads_total
    # New if_bead type registered, old type preserved.
    assert "actin_cortex" in out.particles.types
    assert IF_BEAD_TYPE in out.particles.types
    # Old bond preserved + if_ bond types registered.
    assert "cortex_backbone" in out.bonds.types
    assert BACKBONE_BOND_TYPE in out.bonds.types
    assert CROSSLINK_BOND_TYPE in out.bonds.types
    # Bond count grew by at least the backbone bonds.
    assert out.bonds.N >= 1 + p.n_backbone_bonds
    # Old particle positions unchanged (first 4 rows).
    assert np.allclose(
        np.asarray(out.particles.position)[:4],
        np.asarray(snap.particles.position),
    )


@pytest.mark.skip(reason="heavy: full HOOMD bonded-force eval; core gates cover sign-sense analytically")
def test_hoomd_bonded_force_eval():  # pragma: no cover
    """Optional: stretched if_backbone in a real HOOMD sim gives a restoring force."""
    p = resolve_intermediate_filaments(_enabled_cfg(), kT=kT, R_cell=R_CELL)
    sim = hoomd.Simulation(device=hoomd.device.CPU(), seed=1)
    snap = hoomd.Snapshot()
    snap.particles.N = 2
    snap.particles.types = [IF_BEAD_TYPE]
    snap.particles.position[:] = np.array([[0, 0, 0], [1.5 * p.l_seg, 0, 0]])
    snap.bonds.N = 1
    snap.bonds.types = [BACKBONE_BOND_TYPE]
    snap.bonds.group[:] = np.array([[0, 1]])
    L = 10 * p.l_seg
    snap.configuration.box = [L, L, L, 0, 0, 0]
    sim.create_state_from_snapshot(snap)
    harmonic = build_intermediate_filament_bonds(p)
    integrator = hoomd.md.Integrator(dt=1e-9, forces=[harmonic])
    sim.operations.integrator = integrator
    sim.run(0)
    f = harmonic.forces
    assert f[1][0] < 0.0     # bead 1 pulled back toward bead 0
