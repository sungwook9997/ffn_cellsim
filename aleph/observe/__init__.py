"""Aleph's observation layer: what was measured, how, and whether it may be believed.

Eight modules, one discipline — a number may not leave this package without the protocol that
produced it and an honest statement of its own precision.

``manifest``
    :class:`~aleph.observe.manifest.ProtocolManifest`: the thirty-four things the manuscript's §8.1
    says must be declared before a measurement means anything, with ``NOT_APPLICABLE`` and
    ``NOT_RECORDED`` kept as distinct typed values because the difference between them decides
    whether two datasets may be pooled.

``operator``
    The observation seam. ``raw_observable()`` and ``apparent_quantity()`` stay separate methods,
    because a simulator that computes the true quantity directly and compares it against a published
    one has skipped every step where the disagreement lives. Refusal is a return value here, never
    an exception and never a NaN.

``stationarity``
    Domain-independent numerical-experiment hygiene. ``sem = sigma*sqrt(2*tau_int/T)``, never
    ``sigma/sqrt(N)`` — on a correlated series the naive form is six times too small and a nominal
    95% interval covers the truth 22% of the time.

``fluctuation``
    The first concrete operator: the passive membrane thermal fluctuation spectrum, reported only
    over the wavevector band the mesh can actually support.

``spectral_observer``
    Eigen-spectrum observables of a symmetric operator, with a zero tolerance derived from backward
    stability rather than chosen by hand.

``tether``
    The AFM membrane-tether forward operator. Composes ``sigma_bilayer`` and the membrane-cortex
    attachment energy ``W`` itself, rather than accepting one lumped tension, so that it *can* be
    asked whether the coupling reproduces a measured force. ``f = 2*pi*sqrt(2*kappa*sigma_app)``.

``traction``
    Traction force per cell from TFM — and the first operator here that **refuses on a runtime
    precondition**. Two defects are open on the load path it reads (the discarded
    ``ProtrusionActinTarget`` wrapper, and a delivery probe that scores a row from the ledger rather
    than from motion), so it returns ``UPSTREAM_DEFECT_OPEN`` and names the lanes that own the fix.
    It begins reporting by itself when the audit arrays say the path is sound — no edit required.

``kinematics``
    Spreading rate, volume flux and PIV field speed, and the unit discipline they are exported
    under. Where an external record states no unit, a *magnitude* comparison is a refusal and a
    *direction* comparison is not, because a sign survives an unknown positive scale and a magnitude
    does not. ``DirectionAgreement`` carries ``is_acceptance_criterion = False`` as a field.

This package must never import ``validation/**``. The analytic oracles judge the observation code;
observation code that can reach them can validate itself. ``tests/firewall/`` enforces it.

Nothing here imports :mod:`aleph.state` or :mod:`aleph.units` either. That is a deliberate posture
toward concurrently-developing lanes: the state is reached through the
:class:`~aleph.observe.operator.AcceptedState` Protocol and the mode basis through
:class:`~aleph.observe.fluctuation.ModeBasis`, so neither lane blocks the other and neither can
break the other by settling its API.
"""

from __future__ import annotations

__all__ = [
    "fluctuation",
    "kinematics",
    "manifest",
    "operator",
    "spectral_observer",
    "stationarity",
    "tether",
    "traction",
]
