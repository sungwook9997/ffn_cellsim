# acs_kb — v1 frozen as validation reference

As of 2026-05-19, `acs_kb/` is **frozen** as the v1 validation oracle for the
HOOMD-blue full-fidelity rewrite.

- v1 (numpy + single-chain cortex) is preserved here in entirety: rate laws,
  KU sanity gates, parameter dictionaries, regression tests, and notebooks
  all continue to run and remain the ground truth for v2.
- v2 work begins in a new top-level package, `acs_hoomd/`. Do not add new
  features to `acs_kb/`; bug fixes that affect a KU contract should land
  there and be mirrored into `acs_hoomd/` once it exists.
- The audit that classified every module as **reuse / port / archive** lives
  at `acs_kb/outputs/v2_audit/CODEBASE_AUDIT.md`. Start there before
  touching any file in this package.

**Reference**: Notion *🚀 Simulation Build Plan v2 — HOOMD-blue Full-Fidelity*
— <https://www.notion.so/365120daec5d81799efefcf078f2039e>

**Modules explicitly archived** (replaced by HOOMD-native equivalents):

- `acs_kb/cell/cortex.py` — single-chain circular cortex; v2 uses ~1000
  effective filaments per cell as HOOMD bonded particles.
- `acs_kb/cell/force_balance.py` — overdamped Euler over cortex beads;
  replaced by `hoomd.md.methods.Brownian`.
- `acs_kb/ecm/integrator.py` — Euler-Maruyama numpy step; replaced by
  HOOMD's built-in overdamped integrators.
- `acs_kb/ecm/shear_lees_edwards.py` — hand-rolled Lees-Edwards MI;
  replaced by HOOMD-native sheared box deformation.
