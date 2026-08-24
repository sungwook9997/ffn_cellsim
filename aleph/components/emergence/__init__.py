"""Label-blind emergence detector (increment I5 detector-core) — the Active Cell acceptance layer.

The unified actomyosin network is supposed to *self-organize* into cortex / ventral+dorsal SF /
transverse arcs / perinuclear cap / filopodium (P3, STRESS_FIBER_TARGET). Whether it actually does is
the *falsifiable* emergence gate: bundles must be shown to CONDENSE from an isotropic seed, and the
verdict must be read by a detector that NEVER sees the construction labels (region / type). A detector
that read the labels would confirm the architecture it was handed — not measure emergence.

This package is that detector-core. It consumes ONLY the FF fiber-array geometry contract
(``pos`` / ``fiber_offsets`` / ``n_fibers`` — which PREDATES the I4 unified weave), reuses the
validated nematic order parameter ``S`` (the Q-tensor largest-eigenvalue form from
``ff/architecture_metrics.parallel_order_parameter``), and turns the VACUOUS ``bundle_count = n_fibers``
into a real, localizing condensation measure with a pre-registered effect-size rule against the
finite-N isotropic null.

Modules:
  * ``nematic``      — the validated Q-tensor S + fiber axes/centroids (global order).
  * ``null_model``   — the finite-N isotropic null band (analytic ``E[sum lambda^2] = 3/(2N)`` +
                       Monte-Carlo band) and the pre-registered effect-size z-rule.
  * ``condensation`` — local nematic + density fields and label-blind bundle LOCALIZATION (which
                       spatial patch condenses while the background stays isotropic).
  * ``detector``     — the top-level ``EmergenceDetector.detect()`` structured report (the drop-in
                       replacement for the vacuous ``bundle_count``).
  * ``synthetic``    — pure-NumPy synthetic config generators (the oracle fixtures).

Scope of THIS increment (I5 detector-core) = the detector + synthetic-config oracles ONLY. The
4-toggle ablation-control harness (myosin-off / unanchored / KMC-off / angle-off) and the native
emergence PROOF need I1c + I3 + I4 hooks that do not exist yet — they are specified, not built, in
``INTEGRATION.md``.

Anti-coupling firewall (HARD): no function in this package accepts a region / type / label argument.
The detector is a pure function of geometry. A flat result is a FINDING to PI + a SEEDED-labeled
scaffold fallback — NEVER re-tuned to force condensation (falsifiability contract, §5 / hard-truth #1).

Runtime note: this is the host-side (dev-Mac, pure NumPy/SciPy) acceptance layer. It imports no Warp
and no HOOMD; the native emergence run is gated AGAINST these oracles, never the other way round.
"""
