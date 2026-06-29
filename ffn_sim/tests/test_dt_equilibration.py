"""STATIC + analytical sanity-gate tests for the PI-ratified B1 + B2 items.

Covers:

- B1 — global CFL dt reconciliation (``ffn_sim/cell/dt_reconcile.py``):
    * §1 Dimensional  (TestB1Dimensional — dt>0 in s; cortex-only == dt_cfl;
                       per-element breakdown keys/sums sane)
    * §2 Monotonicity (TestB1Monotone — a stiffer subsystem only LOWERS dt_min;
                       dt_min ≤ safety·τ_e for every active element)
    * §4 Numerical    (TestB1Degenerate — zero/negative/∞ k raises)
    * Builder off-path no-op (TestB1Builder.test_off_path_dt_unchanged)
    * Builder on-path  (TestB1Builder.test_on_path_dt_drops)
- B2 — cell equilibration prelude (``ffn_sim/cell/equilibration.py``):
    * §1 Boundary      (TestB2.test_negative_budget_raises / missing-sim)
    * §3 Non-increase  (TestB2.test_soft_phase_does_not_increase_max_force)
    * Numerical        (TestB2.test_no_nan_after_prelude)
    * No-op fast path   (TestB2.test_noop_both_phases_off)

Mirrors the style of ``ffn_sim/tests/test_enclosed_volume.py`` +
``ffn_sim/tests/test_turnover.py``: small demo cortex (CI-fast), the
closed-form CFL relation ``dt = safety·γ/k`` as the ACCEPTANCE ORACLE, the
runtime being :func:`compute_global_cfl_dt` / :func:`equilibrate_cell`.

The B1 ``p_erm`` channel is exercised with a minimal duck-typed stub whose
only read attribute is ``k_ERM`` — that is the EXACT contract
``compute_global_cfl_dt`` and ``build_cortex_full_simulation`` use for ERM
(the builder reads ``p_erm`` ONLY to feed dt reconciliation; it does not
attach an ERM force in this path). Using a stub keeps the stiff-subsystem
driver parameter-free and grid-invariant (no magic numbers — the stiffness
is derived from ``γ_b`` and the cortex ``τ_min`` so the expected dt is a
closed form ``safety · τ_min / factor``).

The B1 builder ON-path test (``test_on_path_dt_drops``) is currently
marked ``xfail`` because of a pre-existing bug in the B1 *wiring* inside
``ffn_sim/cell/cell.py`` — see that test's docstring + the report. The B1
*module* itself is correct and fully exercised here directly via
:func:`compute_global_cfl_dt`.

Multi-minute checks would be gated behind ``B1B2_PRODUCTION=1`` per the
H3_*_PRODUCTION pattern; the default suite is CI-fast.
"""

from __future__ import annotations

import math
import os
from copy import deepcopy
from pathlib import Path

import numpy as np
import pytest
import yaml

from ffn_sim.archive.hoomd_legacy.cell.cell import build_cortex_full_simulation
from ffn_sim.archive.hoomd_legacy.cell.dt_reconcile import (
    CFLElement,
    GlobalCFLResult,
    compute_global_cfl_dt,
)
from ffn_sim.archive.hoomd_legacy.cell.equilibration import (
    DEFAULT_N_BAOAB,
    DEFAULT_N_SOFTSTART,
    equilibrate_cell,
)
from ffn_sim.archive.hoomd_legacy.cortex.cortex import resolve_h3_derived


CONFIG_PATH = (
    Path(__file__).resolve().parents[1] / "configs" / "phase1_h3.yaml"
)


def _load_cfg() -> dict:
    with open(CONFIG_PATH) as f:
        return yaml.safe_load(f)


def _demo_cortex_cfg(n_filaments: int = 20) -> dict:
    """Small cortex for fast CI runs of the B1/B2 sanity gates."""
    cfg = deepcopy(_load_cfg())
    cfg["cortex"]["n_filaments"] = n_filaments
    cfg["cortex"]["demo_mode"] = True
    return cfg


class _StubStiffSubsystem:
    """Minimal duck-typed ``p_erm`` carrying ONLY ``k_ERM``.

    ``compute_global_cfl_dt`` reads the ERM channel solely through
    ``p_erm.k_ERM``, so this is the exact contract — no need to construct
    the full :class:`ffn_sim.archive.hoomd_legacy.cortex.erm.ResolvedERM` dataclass.
    """

    def __init__(self, k_ERM: float) -> None:
        self.k_ERM = k_ERM


@pytest.fixture(scope="module")
def resolved_cortex():
    return resolve_h3_derived(_demo_cortex_cfg())


def _stiff_erm_for(p_cortex, *, factor: float = 10.0) -> _StubStiffSubsystem:
    """A stiff ERM stub whose τ = γ_b / k_ERM = τ_min / ``factor``.

    Grid-invariant + magic-number-free: the stiffness is derived from the
    cortex drag and τ_min, so the reconciled dt is the closed form
    ``safety · τ_min / factor``.
    """
    k = p_cortex.gamma_b / (p_cortex.tau_min / factor)
    return _StubStiffSubsystem(k_ERM=k)


# ---------------------------------------------------------------------------
# B1 — §1 Dimensional analysis
# ---------------------------------------------------------------------------
class TestB1Dimensional:
    def test_returns_positive_dt_in_seconds(self, resolved_cortex):
        """dt_min is finite, positive, and of seconds magnitude (cortex
        dt_cfl is O(1e-8 s) for the demo cortex)."""
        r = compute_global_cfl_dt(resolved_cortex)
        assert isinstance(r, GlobalCFLResult)
        assert math.isfinite(r.dt_min)
        assert r.dt_min > 0.0
        # Sane seconds scale for an overdamped cortex bead (ps–ms).
        assert 1e-12 < r.dt_min < 1e-3

    def test_cortex_only_equals_dt_cfl(self, resolved_cortex):
        """Boundary: with ONLY the cortex active, dt_min == the cortex
        ``cfl_safety_factor · tau_min`` (== ResolvedH3.dt_cfl). This is the
        acceptance oracle for the cortex-only case."""
        r = compute_global_cfl_dt(resolved_cortex)
        expected = resolved_cortex.cfl_safety_factor * resolved_cortex.tau_min
        assert math.isclose(r.dt_min, expected, rel_tol=1e-12)
        assert math.isclose(r.dt_min, resolved_cortex.dt_cfl, rel_tol=1e-12)
        assert r.binding.subsystem == "cortex"
        assert r.safety == resolved_cortex.cfl_safety_factor

    def test_dt_bound_is_safety_times_tau(self, resolved_cortex):
        """Each element's dt_bound == safety · τ (dimensional consistency:
        [s] = [dimensionless]·[s])."""
        r = compute_global_cfl_dt(resolved_cortex)
        for e in r.elements:
            if math.isfinite(e.tau):
                assert math.isclose(
                    e.dt_bound, r.safety * e.tau, rel_tol=1e-12
                )

    def test_breakdown_keys_and_sums_sane(self, resolved_cortex):
        """Per-element breakdown is non-empty, the binding element is the
        smallest-τ active one, dt_min == safety·τ_binding, and
        format_breakdown mentions dt_min + the binding subsystem."""
        erm = _stiff_erm_for(resolved_cortex)
        r = compute_global_cfl_dt(resolved_cortex, p_erm=erm)
        assert len(r.elements) >= 2  # cortex + erm
        for e in r.elements:
            assert isinstance(e, CFLElement)
            assert e.subsystem and e.label
        active = [
            e for e in r.elements
            if not e.is_constraint_removed and not e.is_informational
        ]
        tau_min_active = min(e.tau for e in active)
        assert math.isclose(r.binding.tau, tau_min_active, rel_tol=1e-12)
        assert math.isclose(r.dt_min, r.safety * tau_min_active, rel_tol=1e-12)
        text = r.format_breakdown()
        assert "dt_min" in text
        assert r.binding.subsystem in text


# ---------------------------------------------------------------------------
# B1 — §2/§3 Monotonicity
# ---------------------------------------------------------------------------
class TestB1Monotone:
    def test_stiffer_subsystem_lowers_dt(self, resolved_cortex):
        """Adding a stiffer subsystem can only LOWER dt_min (monotone min);
        the reconciled value is the closed form safety·γ_b/k_ERM."""
        base = compute_global_cfl_dt(resolved_cortex).dt_min
        erm = _stiff_erm_for(resolved_cortex, factor=10.0)
        r = compute_global_cfl_dt(resolved_cortex, p_erm=erm)
        assert r.dt_min < base
        assert r.binding.subsystem == "erm"
        expected = (
            resolved_cortex.cfl_safety_factor
            * (resolved_cortex.gamma_b / erm.k_ERM)
        )
        assert math.isclose(r.dt_min, expected, rel_tol=1e-12)
        assert math.isclose(
            r.dt_min,
            resolved_cortex.cfl_safety_factor * resolved_cortex.tau_min / 10.0,
            rel_tol=1e-12,
        )

    def test_softer_subsystem_does_not_change_dt(self, resolved_cortex):
        """A subsystem SOFTER than the cortex leaves dt_min at the cortex
        value (the cortex stays binding)."""
        base = compute_global_cfl_dt(resolved_cortex).dt_min
        # τ = 10·τ_min → softer than cortex; k = γ_b/(10·τ_min).
        soft = _StubStiffSubsystem(
            k_ERM=resolved_cortex.gamma_b / (10.0 * resolved_cortex.tau_min)
        )
        r = compute_global_cfl_dt(resolved_cortex, p_erm=soft)
        assert math.isclose(r.dt_min, base, rel_tol=1e-12)
        assert r.binding.subsystem == "cortex"

    def test_dt_min_le_safety_tau_for_every_active_element(self, resolved_cortex):
        """Monotonicity invariant: dt_min ≤ safety·τ_e for every active
        (non-removed, non-informational) element."""
        erm = _stiff_erm_for(resolved_cortex)
        r = compute_global_cfl_dt(resolved_cortex, p_erm=erm)
        for e in r.elements:
            if e.is_constraint_removed or e.is_informational:
                continue
            assert r.dt_min <= e.dt_bound + 1e-30, (
                f"dt_min {r.dt_min:.3e} > dt_bound of "
                f"{e.subsystem}.{e.label} ({e.dt_bound:.3e})"
            )


# ---------------------------------------------------------------------------
# B1 — §4 Numerical sanity (degenerate inputs raise)
# ---------------------------------------------------------------------------
class TestB1Degenerate:
    def test_zero_stiffness_raises(self, resolved_cortex):
        with pytest.raises(ValueError, match="k must be finite-positive"):
            compute_global_cfl_dt(
                resolved_cortex, p_erm=_StubStiffSubsystem(k_ERM=0.0)
            )

    def test_negative_stiffness_raises(self, resolved_cortex):
        with pytest.raises(ValueError, match="k must be finite-positive"):
            compute_global_cfl_dt(
                resolved_cortex, p_erm=_StubStiffSubsystem(k_ERM=-1.0)
            )

    def test_nonfinite_stiffness_raises(self, resolved_cortex):
        with pytest.raises(ValueError, match="k must be finite-positive"):
            compute_global_cfl_dt(
                resolved_cortex, p_erm=_StubStiffSubsystem(k_ERM=math.inf)
            )

    def test_bad_safety_override_raises(self, resolved_cortex):
        with pytest.raises(ValueError, match="cfl_safety_factor"):
            compute_global_cfl_dt(resolved_cortex, safety=-0.1)


# ---------------------------------------------------------------------------
# B1 — builder off-path no-op (default) + on-path drop
# ---------------------------------------------------------------------------
class TestB1Builder:
    def test_off_path_dt_unchanged(self, resolved_cortex):
        """reconcile_dt=False (default) → integrator dt is bit-for-bit the
        cortex dt_cfl. A default build and an explicit reconcile_dt=False
        build are identical (in-process dt-value compare, no git)."""
        h_default = build_cortex_full_simulation(
            resolved_cortex, with_baoab=True,
        )
        h_explicit_off = build_cortex_full_simulation(
            resolved_cortex, with_baoab=True, reconcile_dt=False,
        )
        dt_default = float(h_default["sim"].operations.integrator.dt)
        dt_off = float(h_explicit_off["sim"].operations.integrator.dt)
        assert dt_default == dt_off
        assert dt_default == resolved_cortex.dt_cfl

    def test_on_path_dt_drops(self, resolved_cortex):
        """reconcile_dt=True with a deliberately stiffer subsystem → the
        integrator dt drops to safety·γ_b/k_ERM (B1 wiring now fixed: cell.py
        defines the module logger _LOG, so the reconcile branch runs cleanly).

        Also asserts the BAOAB Action sees the lowered dt and that the B1
        result is surfaced in the handles (cfl_result + dt_used)."""
        erm = _stiff_erm_for(resolved_cortex, factor=10.0)
        h = build_cortex_full_simulation(
            resolved_cortex, with_baoab=True, reconcile_dt=True, p_erm=erm,
        )
        dt = float(h["sim"].operations.integrator.dt)
        expected = (
            resolved_cortex.cfl_safety_factor
            * (resolved_cortex.gamma_b / erm.k_ERM)
        )
        assert dt < resolved_cortex.dt_cfl
        # HOOMD stores the integrator dt in float32, so compare to the
        # float64 oracle at float32 precision.
        assert math.isclose(dt, expected, rel_tol=1e-6)
        # dt_used is the raw float the builder passed through (no float32
        # round-trip) → it equals the float64 cfl_result.dt_min exactly.
        assert math.isclose(h["dt_used"], expected, rel_tol=1e-12)
        # The BAOAB Action runs at the same (lowered) dt as the integrator.
        assert math.isclose(
            float(h["baoab_action"].dt), dt, rel_tol=1e-6
        )
        # B1 results surfaced in handles.
        assert h["cfl_result"] is not None
        assert math.isclose(h["cfl_result"].dt_min, h["dt_used"], rel_tol=1e-12)

    def test_on_path_within_bound_keeps_dt(self, resolved_cortex):
        """reconcile_dt=True with NO stiffer subsystem → the cortex stays
        binding, so dt is unchanged (== cortex dt_cfl); the builder takes the
        'within bound' branch without error and surfaces cfl_result."""
        h = build_cortex_full_simulation(
            resolved_cortex, with_baoab=True, reconcile_dt=True,
        )
        dt = float(h["sim"].operations.integrator.dt)
        # integrator dt is float32-stored; dt_used is the raw float passed in.
        assert math.isclose(dt, resolved_cortex.dt_cfl, rel_tol=1e-6)
        assert math.isclose(h["dt_used"], resolved_cortex.dt_cfl, rel_tol=1e-12)
        assert h["cfl_result"] is not None
        assert h["cfl_result"].binding.subsystem == "cortex"

    def test_off_path_surfaces_none_cfl_result(self, resolved_cortex):
        """reconcile_dt=False (default) → cfl_result is None and dt_used is
        the cortex dt_cfl (the B1 surface is present but inert off-path)."""
        h = build_cortex_full_simulation(resolved_cortex, with_baoab=True)
        assert h["cfl_result"] is None
        assert math.isclose(h["dt_used"], resolved_cortex.dt_cfl, rel_tol=1e-12)


# ---------------------------------------------------------------------------
# B2 — equilibration helper sanity gates (via the dedicated public function)
# ---------------------------------------------------------------------------
class TestB2:
    def test_negative_budget_raises(self, resolved_cortex):
        """§1 Boundary: a negative step budget raises ValueError."""
        h = build_cortex_full_simulation(resolved_cortex, with_baoab=True)
        with pytest.raises(ValueError, match="n_softstart"):
            equilibrate_cell(
                {"sim": h["sim"], "baoab_updater": h["baoab_updater"]},
                n_softstart=-1, n_baoab=0,
                rest_length=resolved_cortex.rest_length,
                gamma_b=resolved_cortex.gamma_b,
            )

    def test_missing_sim_raises(self, resolved_cortex):
        with pytest.raises(KeyError, match="sim"):
            equilibrate_cell(
                {},
                n_softstart=1, n_baoab=0,
                rest_length=resolved_cortex.rest_length,
                gamma_b=resolved_cortex.gamma_b,
            )

    def test_bad_gamma_raises(self, resolved_cortex):
        h = build_cortex_full_simulation(resolved_cortex, with_baoab=True)
        with pytest.raises(ValueError, match="gamma_b"):
            equilibrate_cell(
                {"sim": h["sim"], "baoab_updater": h["baoab_updater"]},
                n_softstart=1, n_baoab=0,
                rest_length=resolved_cortex.rest_length,
                gamma_b=-1.0,
            )

    def test_noop_both_phases_off(self, resolved_cortex):
        """Both phases 0 → no-op fast path: before == after == after_baoab,
        and the reported step counts are 0 (build identical / unchanged)."""
        h = build_cortex_full_simulation(resolved_cortex, with_baoab=True)
        diag = equilibrate_cell(
            {"sim": h["sim"], "baoab_updater": h["baoab_updater"]},
            n_softstart=0, n_baoab=0,
            rest_length=resolved_cortex.rest_length,
            gamma_b=resolved_cortex.gamma_b,
        )
        assert diag["n_softstart"] == 0.0
        assert diag["n_baoab"] == 0.0
        assert diag["max_force_before"] == diag["max_force_after_soft"]
        assert diag["max_force_after_soft"] == diag["max_force_after_baoab"]
        assert math.isfinite(diag["max_force_before"])

    def test_soft_phase_does_not_increase_max_force(self, resolved_cortex):
        """§3 Non-increase: the clipped-Brownian soft phase drains (steepest
        descent); it does not INJECT mechanical energy. Post-soft max|F| must
        not meaningfully exceed the construction max|F|.

        Tolerance note: a freshly-built cortex shell is already at its rest
        configuration, so max|F| sits at the thermal/round-off floor
        (~1e-24 N). At that floor the clipped-Brownian update's per-particle
        ±round-off can nudge the reported max by O(1e-26 N); a small relative
        tolerance keeps the gate honest (it still catches a real force
        INJECTION, which would be many orders of magnitude larger)."""
        h = build_cortex_full_simulation(resolved_cortex, with_baoab=True)
        diag = equilibrate_cell(
            {"sim": h["sim"], "baoab_updater": h["baoab_updater"]},
            n_softstart=50, n_baoab=50,
            rest_length=resolved_cortex.rest_length,
            gamma_b=resolved_cortex.gamma_b,
        )
        assert math.isfinite(diag["max_force_before"])
        assert math.isfinite(diag["max_force_after_soft"])
        assert math.isfinite(diag["max_force_after_baoab"])
        before = diag["max_force_before"]
        after = diag["max_force_after_soft"]
        # No force INJECTION: post-soft ≤ construction (+1% floor tolerance).
        assert after <= before * 1.01 + 1e-30, (
            f"soft phase increased max|F| from {before:.3e} to {after:.3e} N "
            "(more than the round-off floor — possible energy injection)"
        )

    def test_no_nan_after_prelude(self, resolved_cortex):
        """Numerical: positions stay finite through a short soft+BAOAB
        prelude (the int32 image guard would fire first otherwise)."""
        h = build_cortex_full_simulation(resolved_cortex, with_baoab=True)
        equilibrate_cell(
            {"sim": h["sim"], "baoab_updater": h["baoab_updater"]},
            n_softstart=30, n_baoab=30,
            rest_length=resolved_cortex.rest_length,
            gamma_b=resolved_cortex.gamma_b,
        )
        with h["sim"].state.cpu_local_snapshot as s:
            pos = np.asarray(s.particles.position).copy()
        assert np.isfinite(pos).all()


# ---------------------------------------------------------------------------
# B2 — builder-path wiring (equilibrate kwargs consumed by the builder)
# ---------------------------------------------------------------------------
class TestB2BuilderPath:
    def test_builder_equilibrate_populates_diagnostics(self, resolved_cortex):
        """equilibrate_steps>0 → the builder runs the prelude in-line and
        stashes the diagnostics dict under handles['equilibrate_diagnostics'];
        positions stay finite (no NaN), and the reported step counts match."""
        h = build_cortex_full_simulation(
            resolved_cortex, with_baoab=True,
            equilibrate=True, equilibrate_steps=20,
            equilibrate_softstart_steps=10,
        )
        diag = h["equilibrate_diagnostics"]
        assert diag is not None
        assert diag["n_softstart"] == 10.0
        assert diag["n_baoab"] == 20.0
        for key in (
            "max_force_before", "max_force_after_soft", "max_force_after_baoab",
        ):
            assert math.isfinite(diag[key])
        with h["sim"].state.cpu_local_snapshot as s:
            pos = np.asarray(s.particles.position).copy()
        assert np.isfinite(pos).all()

    def test_builder_equilibrate_steps_zero_default_no_key(self, resolved_cortex):
        """Default (equilibrate=False, equilibrate_steps=0) → the prelude is
        NEVER run: equilibrate_diagnostics is None and the integrator dt is
        the plain cortex dt_cfl (bit-for-bit pre-B2 behaviour)."""
        h = build_cortex_full_simulation(resolved_cortex, with_baoab=True)
        assert h["equilibrate_diagnostics"] is None
        # integrator dt is float32-stored; compare at float32 precision.
        assert math.isclose(
            float(h["sim"].operations.integrator.dt),
            resolved_cortex.dt_cfl, rel_tol=1e-6,
        )

    def test_builder_equilibrate_flag_only(self, resolved_cortex):
        """equilibrate=True with equilibrate_steps=0 still triggers the
        prelude (the helper's no-op fast path runs and returns a diagnostics
        dict with zero step counts)."""
        h = build_cortex_full_simulation(
            resolved_cortex, with_baoab=True,
            equilibrate=True, equilibrate_steps=0,
            equilibrate_softstart_steps=0,
        )
        diag = h["equilibrate_diagnostics"]
        assert diag is not None
        assert diag["n_softstart"] == 0.0
        assert diag["n_baoab"] == 0.0


# ---------------------------------------------------------------------------
# Opt-in production-scale prelude (multi-minute; gated like H3_*_PRODUCTION)
# ---------------------------------------------------------------------------
@pytest.mark.skipif(
    not bool(int(os.environ.get("B1B2_PRODUCTION", "0"))),
    reason=(
        "Full-budget (DEFAULT_N_SOFTSTART + DEFAULT_N_BAOAB) prelude on a "
        "production-scale cell; opt-in via B1B2_PRODUCTION=1. Multi-minute "
        "wall (same rationale as the H3_*_PRODUCTION opt-ins)."
    ),
)
class TestB1B2Production:
    def test_full_budget_prelude_completes(self, resolved_cortex):
        h = build_cortex_full_simulation(resolved_cortex, with_baoab=True)
        diag = equilibrate_cell(
            {"sim": h["sim"], "baoab_updater": h["baoab_updater"]},
            n_softstart=DEFAULT_N_SOFTSTART, n_baoab=DEFAULT_N_BAOAB,
            rest_length=resolved_cortex.rest_length,
            gamma_b=resolved_cortex.gamma_b,
        )
        assert math.isfinite(diag["max_force_after_baoab"])
        with h["sim"].state.cpu_local_snapshot as s:
            pos = np.asarray(s.particles.position).copy()
        assert np.isfinite(pos).all()
