"""Resume-gating tests for the mcf7_fullcell_stage1 checkpoint/resume.

The driver now skips the warm-up + all completed samples on a re-run of the
SAME config (no redo-every-time). Correctness hinges on the fingerprint: it
MUST be deterministic for identical params and MUST change when any
physics/schedule param changes, so a resume can never splice an incompatible
run. The end-to-end skip-warm-up behaviour is exercised by running the driver
twice (see the session log); these are the fast pure-function guards.
"""
from __future__ import annotations

from ffn_sim.scripts.mcf7_fullcell_stage1 import _run_fingerprint

_BASE = dict(
    stepping_mode="grip_walk", n_fil=1000, n_motors=200, n_xl=2000,
    force_scaling=True, v0_accel=1.0, n_warmup=8000, n_sample_eff=6,
    interval_eff=20000, backbone_nm=700, kon_scale=1.0, bind_scale=1.0,
    connected_mesh=True, turgor_pa=None, motors_off=False,
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
        ("bind_scale", 6.0), ("connected_mesh", False), ("turgor_pa", 40.0),
        ("motors_off", True), ("compartments_on", False),
        ("only", "p_cytoplasm"), ("aggregation", False), ("seed", 2),
        ("stepping_mode", "binned_r0"), ("force_scaling", False),
    ]:
        changed = dict(_BASE, **{key: alt})
        assert _run_fingerprint(**changed) != base, f"{key} did not change the fingerprint"


def test_turgor_none_vs_zero_distinct():
    # None (unpressurised default) and 0.0 Pa must be distinguishable so a
    # baseline run never resumes from a turgor run and vice-versa.
    fp_none = _run_fingerprint(**dict(_BASE, turgor_pa=None))
    fp_zero = _run_fingerprint(**dict(_BASE, turgor_pa=0.0))
    assert fp_none != fp_zero
