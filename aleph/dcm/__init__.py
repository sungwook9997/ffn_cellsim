"""NVIDIA Warp port of the ffn_cellsim runtime — the Phase-C GPU-resident, differentiable DCM engine.

**Adopted as the going-forward DCM engine (PI 2026-06-21).** Each module here mirrors a
committed HOOMD reference piece (all DELETED 2026-07-29; git history only) (``integrator/``, ``cell/``
forces) as a ``warp.kernel``, once gated by bit-parity tests in ``aleph/tests/dcm/``.
**Those tests and their fixtures were DELETED 2026-07-29 with the HOOMD tree**, so the
parity contract below is HISTORY, not a live gate: per ``STATE.md`` (c) 13 the 25
``warp-*`` claims are permanently uncheckable, because regenerating them would require
the runtime this repo forbids. This tree is PARKED and outside the native import closure.

See ``ENGINE.md`` (this directory) for the engine architecture + entry points, and
``docs/v2_audit/_historical/PHASE_C_WARP_DCM_ENGINE_PLAN_2026-06-21.md`` for the adoption plan
(``_historical/WARP_SPIKE_THEN_DECOHESION_PLAN_2026-06-19.md`` is the original spike plan).

Parity status (the anti-drift contract):
  * CPU backend = the gate (Warp-CPU vs the numpy/HOOMD reference): force-laws machine-eps,
    iterative (M-SHAKE) / atomic-reduce kernels at a documented tol.
  * GPU backend = **VERIFIED on the gbook A5000** (G1, 2026-06-21): the SAME kernels reproduce
    parity on Warp's CUDA backend — force-laws machine-eps, atomic-reduce < 1e-8, reverse-mode
    autodiff clean. See ``docs/v2_audit/PHASE_C_GPU_PARITY_2026-06-21.md``.
  * kT=0 → bit-for-bit on CPU (< 1e-9); on GPU machine-eps (FMA reordering), still < 1e-9.
"""
