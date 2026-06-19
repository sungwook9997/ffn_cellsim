"""NVIDIA Warp port of the ffn_cellsim runtime (Phase B — parity spike).

Each module here mirrors a committed HOOMD reference piece (``integrator/``,
``native/ffn_hoomd_plugin/``, ``cell/`` forces) as a ``warp.kernel`` and is
gated by a committed bit-parity test in ``ffn_sim/tests/warp_port/`` against a
committed HOOMD reference fixture. Nothing here is "done" without a passing
committed parity test (see
``docs/v2_audit/WARP_SPIKE_THEN_DECOHESION_PLAN_2026-06-19.md`` §Phase B).

Parity policy (the anti-drift contract):
  * Warp has a CPU backend → parity runs on CPU (Warp-CPU vs the numpy HOOMD
    reference). The gbook A5000 is only for an optional speed benchmark.
  * kT=0 → bit-for-bit (< 1e-9). kT>0 → identical injected noise, drift < 1e-7
    (float-op-ordering only — the RNG implementation is a separate concern).
"""
