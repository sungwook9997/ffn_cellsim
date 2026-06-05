"""Resume-gating tests for the mcf7_fullcell_stage1 checkpoint/resume.

The driver now skips the warm-up + all completed samples on a re-run of the
SAME config (no redo-every-time). Correctness hinges on the fingerprint: it
MUST be deterministic for identical params and MUST change when any
physics/schedule param changes, so a resume can never splice an incompatible
run. The end-to-end skip-warm-up behaviour is exercised by running the driver
twice (see the session log); these are the fast pure-function guards.
"""
from __future__ import annotations

import pytest

from ffn_sim.scripts.mcf7_fullcell_stage1 import (
    DEFAULT_TURGOR_PA,
    _require_physiological_pressure,
    _require_production_device,
    _run_fingerprint,
)

_BASE = dict(
    stepping_mode="grip_walk", n_fil=1000, n_motors=200, n_xl=2000,
    force_scaling=True, v0_accel=1.0, n_warmup=8000, n_sample_eff=6,
    interval_eff=20000, backbone_nm=700, kon_scale=1.0, bind_scale=1.0,
    connected_mesh=True, turgor_pa=DEFAULT_TURGOR_PA, motors_off=False,
    compartments_on=True, only=None, aggregation=True, seed=1,
)


def test_fingerprint_deterministic():
    assert _run_fingerprint(**_BASE) == _run_fingerprint(**_BASE)


def test_fingerprint_sensitive_to_every_param():
    base = _run_fingerprint(**_BASE)
    for key, alt in [
        ("n_fil", 400), ("n_motors", 500), ("n_xl", 800),
        ("v0_accel", 2.0), ("n_warmup", 4000), ("n_sample_eff", 10),
        ("interval_eff", 10000), ("backbone_nm", 300.0), ("kon_scale", 2.0),
        ("bind_scale", 6.0), ("connected_mesh", False), ("turgor_pa", None),
        ("motors_off", True), ("compartments_on", False),
        ("only", "p_cytoplasm"), ("aggregation", False), ("seed", 2),
        ("stepping_mode", "binned_r0"), ("force_scaling", False),
    ]:
        changed = dict(_BASE, **{key: alt})
        assert _run_fingerprint(**changed) != base, f"{key} did not change the fingerprint"


def test_turgor_none_vs_zero_vs_physiological_distinct():
    # None/0.0 Pa are dev/attribution states, while DEFAULT_TURGOR_PA is the
    # physiological production baseline. They must never resume into each other.
    fp_none = _run_fingerprint(**dict(_BASE, turgor_pa=None))
    fp_zero = _run_fingerprint(**dict(_BASE, turgor_pa=0.0))
    fp_phys = _run_fingerprint(**dict(_BASE, turgor_pa=DEFAULT_TURGOR_PA))
    assert fp_none != fp_zero
    assert fp_none != fp_phys
    assert fp_zero != fp_phys


def test_cpu_requires_explicit_dev_escape_hatch():
    with pytest.raises(ValueError, match="dev-only"):
        _require_production_device("cpu", allow_cpu_dev=False)
    _require_production_device("cpu", allow_cpu_dev=True)


def test_full_cell_requires_positive_turgor_unless_dev_escape_hatch():
    with pytest.raises(ValueError, match="baseline osmotic turgor"):
        _require_physiological_pressure(
            compartments_on=True,
            only=None,
            turgor_pa=0.0,
            allow_unpressurized_dev=False,
        )
    _require_physiological_pressure(
        compartments_on=True,
        only=None,
        turgor_pa=DEFAULT_TURGOR_PA,
        allow_unpressurized_dev=False,
    )
    _require_physiological_pressure(
        compartments_on=True,
        only=None,
        turgor_pa=0.0,
        allow_unpressurized_dev=True,
    )
