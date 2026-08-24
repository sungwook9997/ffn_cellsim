"""Closed-form acceptance oracles for passive single-cell/sphere mechanics.

These modules support the hierarchical single-cell → spheroid mechanics
validation flow (``docs/v2_audit/HIERARCHICAL_MECHANICS_VALIDATION_PLAN_2026-07-14.md``,
advisor flow, Stages S1–S2). They wrap the classical closed forms a passive
homogeneous elastic/viscoelastic body must satisfy so that the FF emergent
force–displacement / stress / relaxation response can be gated against an
analytic ground truth BEFORE a run is accepted — never fit to a run afterwards.

VALIDATION oracles (runtime-import-forbidden), same discipline as the rest of
``aleph.validation.oracles``: pure functional forms + geometry, no physics
magic constants, SI units throughout.

Modules
-------
- ``hertz``                : forward Hertz contact (sphere–plane, sphere-between-plates),
                             the S1 distributed-load force–displacement oracle.
- (planned) ``thin_shell`` : pressurized-shell inflation / indentation (turgor regime).
- (planned) ``viscoelastic_relaxation`` : Maxwell / SLS / Kelvin-Voigt G(t), J(t).
- (planned) ``boussinesq`` : point-load half-space stress decay (S2 spatial-decay gate).
"""
