"""NVIDIA Warp port of the ffn_cellsim runtime — the Phase-C GPU-resident, differentiable DCM engine.

**Adopted as the going-forward DCM engine (PI 2026-06-21).** Each module here mirrors a
committed HOOMD reference piece (``integrator/``, ``native/ffn_hoomd_plugin/``, ``cell/``
forces) as a ``warp.kernel``, gated by a committed bit-parity test in
``ffn_sim/tests/warp_port/`` against a committed HOOMD reference fixture. Nothing is
"done" without a passing committed parity test.

See ``ENGINE.md`` (this directory) for the engine architecture + entry points, and
``docs/v2_audit/PHASE_C_WARP_DCM_ENGINE_PLAN_2026-06-21.md`` for the adoption plan
(``WARP_SPIKE_THEN_DECOHESION_PLAN_2026-06-19.md`` is the original spike plan).

Parity status (the anti-drift contract):
  * CPU backend = the gate (Warp-CPU vs the numpy/HOOMD reference): force-laws machine-eps,
    iterative (M-SHAKE) / atomic-reduce kernels at a documented tol.
  * GPU backend = **VERIFIED on the gbook A5000** (G1, 2026-06-21): the SAME kernels reproduce
    parity on Warp's CUDA backend — force-laws machine-eps, atomic-reduce < 1e-8, reverse-mode
    autodiff clean. See ``docs/v2_audit/PHASE_C_GPU_PARITY_2026-06-21.md``.
  * kT=0 → bit-for-bit on CPU (< 1e-9); on GPU machine-eps (FMA reordering), still < 1e-9.
"""
