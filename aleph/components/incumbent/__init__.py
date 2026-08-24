"""ac/cell — the lead-owned cell assembly + one-physical-clock driver (P5).

Net-new lead infrastructure that COMPOSES the per-track force primitives of the new Active Cell engine
into ONE runnable native cell under a single outer-physical / inner-mechanical scheduler, RETIRING the
lumped ``ff/`` mechanisms by NON-USE (the new driver never launches them). See ``assemble.py`` (builds the
physiological t0 cell) and ``driver.py`` (the outer/inner loop + the ``--from-resting`` stability run).

Runtime: NVIDIA Warp on CUDA GPU only (I0-A). Authored on the dev Mac; every device kernel runs on the
gbook A5000. HOOMD is never imported.
"""
