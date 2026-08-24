r"""Resting bound-myosin SETPOINT gate — the fine-grained-faithful resting cortical-tension source.

Host (dev-Mac) tests exercise the pure-NumPy planner + the source-gated contract; the CUDA-only tests apply the
seed to a real device segment-hand state and build a small composed cell with the setpoint ON.  The seed is a
build-time construction (like ERM/LINC pairing), so the planner runs anywhere and only the device write / full
build require a GPU.
"""

from __future__ import annotations

import numpy as np
import pytest

from aleph.components.motor.resting_setpoint import (
    RestingBoundMyosinSetpoint,
    RestingBoundPlan,
    plan_resting_bound_heads,
    resting_tangential_load_reference,
)


def _straight_actin_patch():
    """6 short actin segments along x + 2 heads above each midpoint (within 0.05 µm capture)."""
    seg_mid = np.array([[i * 0.5, 0.0, 0.0] for i in range(6)], np.float64)
    acts, sa, sb, pol = [], [], [], []
    for i in range(6):
        acts.append(seg_mid[i] + [-0.15, 0.0, 0.0])
        acts.append(seg_mid[i] + [0.15, 0.0, 0.0])
        sa.append(2 * i)
        sb.append(2 * i + 1)
        pol.append(1)
    acts = np.array(acts, np.float64)
    heads = np.array([seg_mid[i] + off for i in range(6) for off in ([0, 0.03, 0], [0, 0.04, 0])], np.float64)
    pos = np.vstack([acts, heads])
    head_node = np.arange(acts.shape[0], acts.shape[0] + heads.shape[0], dtype=np.int64)
    return pos, head_node, np.array(sa, np.int64), np.array(sb, np.int64), np.array(pol, np.int64)


# ── source-gated contract ────────────────────────────────────────────────────────────────────────
def test_setpoint_contract_rejects_partial_or_defaulted() -> None:
    with pytest.raises(ValueError):
        RestingBoundMyosinSetpoint(fraction=1.5, per_head_force_pn=0.5, source="s")     # fraction > 1
    with pytest.raises(ValueError):
        RestingBoundMyosinSetpoint(fraction=0.0, per_head_force_pn=0.5, source="s")     # fraction 0
    with pytest.raises(ValueError):
        RestingBoundMyosinSetpoint(fraction=0.1, per_head_force_pn=-1.0, source="s")    # force <= 0
    with pytest.raises(ValueError):
        RestingBoundMyosinSetpoint(fraction=0.1, per_head_force_pn=0.5, source="   ")   # no source
    ok = RestingBoundMyosinSetpoint(fraction=0.1, per_head_force_pn=0.5, source="TEST; PI-GAP")
    assert ok.fraction == 0.1 and ok.per_head_force_pn == 0.5


# ── planner (pure NumPy, dev-Mac) ──────────────────────────────────────────────────────────────────
def test_planner_binds_requested_fraction_and_realizes_f_head() -> None:
    pos, head_node, sa, sb, pol = _straight_actin_patch()
    sp = RestingBoundMyosinSetpoint(fraction=0.5, per_head_force_pn=0.5, source="TEST; PI-GAP",
                                    capture_radius_um=0.05)
    plan = plan_resting_bound_heads(pos, head_node, sa, sb, pol, sp, k_xb=1000.0, r0_xb=0.0,
                                    capture_radius_um=0.05)
    assert plan.n_bound == 6                              # round(0.5 * 12)
    assert plan.diagnostics["n_target"] == 6
    # every seeded head carries EXACTLY the requested per-head tangential tension (nearest-point construction)
    loads = resting_tangential_load_reference(pos, head_node, plan, k_xb=1000.0, r0_xb=0.0)
    assert np.allclose(loads, 0.5, atol=1e-9)
    # the abscissa is a PHYSICAL, non-negative power-stroke displacement ≈ f_head/k_xb (no bond pre-strain)
    assert np.all(plan.abscissa >= 0.0)
    assert np.allclose(plan.abscissa, 0.5 / 1000.0, atol=1e-6)


def test_planner_is_deterministic_and_seed_free() -> None:
    pos, head_node, sa, sb, pol = _straight_actin_patch()
    sp = RestingBoundMyosinSetpoint(fraction=0.5, per_head_force_pn=0.4, source="TEST; PI-GAP",
                                    capture_radius_um=0.05)
    a = plan_resting_bound_heads(pos, head_node, sa, sb, pol, sp, k_xb=1000.0, r0_xb=0.0, capture_radius_um=0.05)
    b = plan_resting_bound_heads(pos, head_node, sa, sb, pol, sp, k_xb=1000.0, r0_xb=0.0, capture_radius_um=0.05)
    assert np.array_equal(a.head_index, b.head_index)
    assert np.array_equal(a.abscissa, b.abscissa)


def test_planner_skips_heads_outside_capture() -> None:
    """Heads beyond the capture reach are ineligible; the shortfall is recorded, never faked up."""
    pos, head_node, sa, sb, pol = _straight_actin_patch()
    # push the heads far from the actin so NONE are within a tiny capture radius
    pos[head_node] += np.array([0.0, 2.0, 0.0])
    sp = RestingBoundMyosinSetpoint(fraction=1.0, per_head_force_pn=0.5, source="TEST; PI-GAP",
                                    capture_radius_um=0.05)
    plan = plan_resting_bound_heads(pos, head_node, sa, sb, pol, sp, k_xb=1000.0, r0_xb=0.0,
                                    capture_radius_um=0.05)
    assert plan.n_bound == 0
    assert plan.diagnostics["n_eligible"] == 0
    assert plan.diagnostics["eligible_shortfall"] is True


def test_planner_fraction_one_binds_all_eligible() -> None:
    pos, head_node, sa, sb, pol = _straight_actin_patch()
    sp = RestingBoundMyosinSetpoint(fraction=1.0, per_head_force_pn=0.5, source="TEST; PI-GAP",
                                    capture_radius_um=0.05)
    plan = plan_resting_bound_heads(pos, head_node, sa, sb, pol, sp, k_xb=1000.0, r0_xb=0.0,
                                    capture_radius_um=0.05)
    assert plan.n_bound == 12
    loads = resting_tangential_load_reference(pos, head_node, plan, k_xb=1000.0, r0_xb=0.0)
    assert np.allclose(loads, 0.5, atol=1e-9)


def test_build_cell_setpoint_is_all_or_none_before_any_cuda() -> None:
    """A partial setpoint is rejected at config validation (before wp.init) — runs on the dev Mac."""
    from aleph.components.incumbent.assemble import CellConfig, build_cell

    with pytest.raises(ValueError, match="all-or-none"):
        build_cell(CellConfig(resting_bound_myosin_fraction=0.1))          # missing per-head force
    with pytest.raises(ValueError, match="all-or-none"):
        build_cell(CellConfig(resting_bound_myosin_force_pn=0.5))          # missing fraction


# ── device apply + full build (CUDA only) ─────────────────────────────────────────────────────────
def _cuda_device():
    wp = pytest.importorskip("warp")
    wp.init()
    dev = next((d for d in wp.get_devices() if d.is_cuda), None)
    if dev is None:
        pytest.skip("I0-A: resting bound-myosin device gate requires a CUDA GPU")
    return wp, str(dev)


def test_apply_writes_device_segment_hand_state() -> None:
    wp, device = _cuda_device()
    from aleph.components.motor.resting_setpoint import apply_resting_bound_heads
    from aleph.components.motor.segment_motor import allocate_segment_hand_state

    n_heads = 8
    state = allocate_segment_hand_state(n_heads, device=device)
    plan = RestingBoundPlan(
        head_index=np.array([1, 4, 6], np.int64),
        seg_id=np.array([0, 2, 3], np.int64),
        seg_a=np.array([10, 12, 14], np.int64),
        seg_b=np.array([11, 13, 15], np.int64),
        bary_t=np.array([0.5, 0.25, 0.75], np.float64),
        walk_dir=np.tile(np.array([1.0, 0.0, 0.0]), (3, 1)),
        abscissa=np.array([5e-4, 6e-4, 7e-4], np.float64),
        diagnostics={},
    )
    n = apply_resting_bound_heads(state, plan, device=device)
    assert n == 3
    bound = state["bound"].numpy()
    assert bound.sum() == 3
    assert bound[1] == 1 and bound[4] == 1 and bound[6] == 1
    assert np.array_equal(np.nonzero(bound)[0], np.array([1, 4, 6]))
    seg_a = state["seg_a"].numpy()
    assert seg_a[1] == 10 and seg_a[4] == 12 and seg_a[6] == 14
    assert np.allclose(state["abscissa"].numpy()[[1, 4, 6]], [5e-4, 6e-4, 7e-4])


def test_build_cell_resting_setpoint_on_vs_off() -> None:
    wp, device = _cuda_device()
    from aleph.components.incumbent.assemble import CellConfig, build_cell

    base = dict(n_filaments=400, with_nucleus=False, with_membrane=True, with_pressure=False,
                membrane_subdivisions=2, device=device)
    off = build_cell(CellConfig(**base))
    assert off.ledger["resting_bound_myosin_status"] == "OFF_UNBOUND_AT_REST"
    assert off.ledger["resting_bound_myosin_n_bound"] == 0
    assert "PI-GAP" in off.ledger["resting_bound_myosin_gap"]
    # OFF parity: the production segment motor is unbound at rest.
    assert off.myosin.segment_runtime.state["bound"].numpy().sum() == 0

    on = build_cell(CellConfig(
        resting_bound_myosin_fraction=0.2, resting_bound_myosin_force_pn=0.5,
        resting_bound_myosin_source="TEST fraction — NOT physiological; PI-GAP",
        resting_bound_myosin_capture_um=0.6, **base))   # generous TEST reach (heads sit ~0.2 µm off the shell)
    assert on.ledger["resting_bound_myosin_status"] == "ON_SEEDED_ISOMETRIC"
    n_bound = int(on.ledger["resting_bound_myosin_n_bound"])
    assert n_bound > 0
    assert on.myosin.segment_runtime.state["bound"].numpy().sum() == n_bound
    # every seeded head carries the requested resting per-head tension (device load kernel readback)
    assert on.ledger["resting_bound_myosin_realized_load_err_pn"] < 1e-6
    on.myosin.compute_loads(on.pos_d)
    wp.synchronize_device(device)
    loads = on.myosin.loads.numpy()
    bound = on.myosin.segment_runtime.state["bound"].numpy()
    assert np.allclose(loads[bound == 1], 0.5, atol=1e-4)
    # ON engages radial ERM pairing so the tension can reach the membrane
    assert on.ledger["resting_bound_myosin_erm_pairing"] == "RADIAL"
    assert on.membrane.erm_pairing_mode.startswith("RADIAL")
