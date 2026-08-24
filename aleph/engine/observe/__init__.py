r"""Engine observables — built ONCE here so seven tracks cannot each grow their own instrument.

The 2026-07-25 compartment-track design left ~118 observables unbuilt, and seven of the tracks
independently specify per-connector two-sided resultants plus adjoint-work closure while four
independently specify per-component transmission matrices.  Built per track, that is seven
implementations of one measurement, free to disagree — and a disagreement between two instruments
measuring the same thing is indistinguishable, in an artifact, from a physics difference.  Hence one
package, imported by drivers, with the physics left where it lives.

Nothing in here computes a force, owns state, or decides acceptance.  Every module takes a callable
(a matrix-free tangent, a force evaluation, an observable) plus arrays, and returns a measurement — so
each is exercisable against a pure-NumPy reference with no device, while the driver hands it the real
Warp-CUDA launches.

Modules:
    * :mod:`~aleph.engine.observe.operator_probe` — explicit assembly of a matrix-free tangent,
      and ``‖K − Kᵀ‖`` as the magnitude-free double-count guard.
    * :mod:`~aleph.engine.observe.spectrum` — the eigen-spectrum of ``K = ∇²U`` as a field, the
      true ``λ_max`` behind the Gershgorin CFL bound, the null-space count, and mode localisation
      (the stress decay length); plus the matrix-free Lanczos path, which is the only way to reach
      the operator at native population, since column probing refuses above 8,192 DOF by design.
    * :mod:`~aleph.engine.observe.sensitivity` — the log-log Jacobian, the Fisher information
      matrix, and the sloppiness spectrum that bounds how many parameters are inferable at all.
    * :mod:`~aleph.engine.observe.energy` — closed-loop work (an energy test needing no energy
      function), the accepted-step energy balance built out of the step's own force launches, and the
      mobility scaling that says whether a non-closing balance is the integrator or a leak.
    * :mod:`~aleph.engine.observe.artifact` — the self-stamping, two-axis, gate-classified run
      record every driver writes.
"""

from __future__ import annotations

from aleph.engine.observe.artifact import (
    ARTIFACT_SCHEMA,
    build_stamp,
    config_hash,
    observation_artifact,
    timing_block,
    write_artifact,
)

__all__ = [
    "ARTIFACT_SCHEMA",
    "CENTRAL_DIFFERENCE_RELATIVE_STEP",
    "BalanceConvergence",
    "ClosedLoopWork",
    "EnergyBalance",
    "FisherReport",
    "LanczosSpectrum",
    "LogSensitivity",
    "LoopConvergence",
    "ModeLocalization",
    "NodeMatvec",
    "StepEnergyLedger",
    "StiffnessSpectrum",
    "SymmetryReport",
    "accumulate_squared_displacement_kernel",
    "assemble_dense_operator",
    "balance_convergence",
    "build_stamp",
    "closed_loop_work",
    "config_hash",
    "dissipated_work",
    "fisher_report",
    "lanczos_extremal_spectrum",
    "lanczos_ritz",
    "log_log_jacobian",
    "loop_convergence",
    "mode_localization",
    "observation_artifact",
    "paired_path_work",
    "path_work",
    "step_energy_balance",
    "stiffness_spectrum",
    "symmetry_report",
    "timing_block",
    "write_artifact",
]


#: Analysis modules resolved ON DEMAND (PEP 562).  `artifact` stays eager because every driver writes a
#: run record; these four are instruments a particular measurement reaches for, and importing them here
#: put ~2,000 lines of eigensolvers, Fisher information and energy ledgers into EVERY native run — the
#: same eager-package defect fixed one level up in `ac/engine/__init__` on 2026-07-28, missed because
#: fixing a package is not the same as looking at its subpackages (2026-07-29).
_LAZY_MODULES: dict[str, tuple[str, ...]] = {
    "aleph.engine.observe.energy": ('BalanceConvergence', 'ClosedLoopWork', 'EnergyBalance', 'LoopConvergence', 'StepEnergyLedger', 'accumulate_squared_displacement_kernel', 'balance_convergence', 'closed_loop_work', 'dissipated_work', 'loop_convergence', 'paired_path_work', 'path_work', 'step_energy_balance'),
    "aleph.engine.observe.operator_probe": ('NodeMatvec', 'SymmetryReport', 'assemble_dense_operator', 'symmetry_report'),
    "aleph.engine.observe.sensitivity": ('CENTRAL_DIFFERENCE_RELATIVE_STEP', 'FisherReport', 'LogSensitivity', 'fisher_report', 'log_log_jacobian'),
    "aleph.engine.observe.spectrum": ('LanczosSpectrum', 'ModeLocalization', 'StiffnessSpectrum', 'lanczos_extremal_spectrum', 'lanczos_ritz', 'mode_localization', 'stiffness_spectrum'),
}
_LAZY_NAMES: dict[str, str] = {
    name: module for module, names in _LAZY_MODULES.items() for name in names
}


def __getattr__(name: str):
    """Resolve a deferred analysis symbol, so `from ac.engine.observe import X` keeps working."""
    module = _LAZY_NAMES.get(name)
    if module is None:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    import importlib

    return getattr(importlib.import_module(module), name)


def __dir__() -> list[str]:
    """Include the deferred names so tab-completion and `dir()` still show the full surface."""
    return sorted(set(globals()) | set(_LAZY_NAMES))
