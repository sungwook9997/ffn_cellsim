"""Self-test of the top-level detector + the anti-coupling firewall (pure NumPy, no Warp/CUDA).

Proves the detector is a pure function of geometry (no region/type/label input anywhere in its API — the
I5 firewall), that its report is honest at the synthetic corners, and that ``n_emerged_bundles`` is the
real replacement for the vacuous ``bundle_count = n_fibers``.
"""

from __future__ import annotations

import inspect
import types

import numpy as np
import pytest

from aleph.components.emergence import condensation, detector, nematic, null_model, synthetic
from aleph.components.emergence import synthetic as syn
from aleph.components.emergence.detector import EmergenceDetector, detect

# Words that would betray a construction label leaking into the geometry-only detector API.
_FORBIDDEN_PARAM_TOKENS = ("label", "region", "type", "cortex", "bundle_id", "kind", "class")
_PUBLIC_MODULES = (nematic, condensation, detector, null_model, synthetic)


def test_no_public_function_accepts_a_label_argument() -> None:
    """ANTI-COUPLING FIREWALL: no public detector function takes a region/type/label parameter.

    The detector must read geometry only; if it could read the construction labels it would confirm the
    architecture it was handed instead of measuring emergence. Enforced structurally over every public
    signature in the package.
    """
    offenders: list[str] = []
    for mod in _PUBLIC_MODULES:
        for name in getattr(mod, "__all__", []):
            obj = getattr(mod, name)
            if not callable(obj):
                continue
            try:
                params = inspect.signature(obj).parameters
            except (ValueError, TypeError):
                continue
            for pname in params:
                low = pname.lower()
                if any(tok in low for tok in _FORBIDDEN_PARAM_TOKENS):
                    offenders.append(f"{mod.__name__}.{name}({pname})")
    assert not offenders, "label-like parameter leaked into the detector API: " + ", ".join(offenders)


def test_detect_signature_is_geometry_only() -> None:
    """``EmergenceDetector.detect`` takes exactly (pos, fiber_offsets) — the FF fiber-array contract."""
    params = list(inspect.signature(EmergenceDetector.detect).parameters)
    assert params == ["self", "pos", "fiber_offsets"]


def test_report_aligned_is_ordered_one_bundle() -> None:
    """A globally aligned network: ordered, one emerged bundle covering ~all fibers."""
    pos, off = syn.aligned_network(300, seed=1, jitter_deg=6.0)
    rep = EmergenceDetector(null_mc=1500).detect(pos, off)
    assert rep.is_ordered_global
    assert rep.z_global > null_model.EFFECT_SIZE_Z_CRIT
    assert rep.n_emerged_bundles == 1
    assert rep.condensed_fraction > 0.8


def test_report_isotropic_is_flat() -> None:
    """A homogeneous isotropic network: NOT ordered, no emerged bundle, ~zero condensed fraction (the
    falsifiable NEGATIVE — a flat result the detector must be willing to return)."""
    pos, off = syn.isotropic_network(440, seed=2)
    rep = EmergenceDetector(null_mc=1500).detect(pos, off)
    assert not rep.is_ordered_global
    assert rep.n_emerged_bundles == 0
    assert rep.condensed_fraction == 0.0


def test_localized_bundle_flies_under_the_global_radar() -> None:
    """A small localized bundle need NOT lift global S past the bar, yet is localized: this is exactly
    why a global-only order parameter is insufficient and the local detector is required."""
    pos, off, mask, center = syn.planted_bundle(400, 40, seed=3)
    rep = EmergenceDetector(null_mc=1500).detect(pos, off)
    assert rep.n_emerged_bundles >= 1
    # the honest bundle count is O(1), NOT the vacuous n_fibers
    assert rep.n_emerged_bundles < rep.n_fibers
    assert 0.0 < rep.condensed_fraction < 0.3


def test_honest_bundle_count_replaces_vacuous_one() -> None:
    """``n_emerged_bundles`` is a measured O(1) count, unlike ``bundle_count = n_fibers`` (which would
    return the whole fiber count regardless of whether anything condensed)."""
    from aleph.laws.architecture_metrics import bundle_count

    pos, off = syn.isotropic_network(440, seed=4)
    net = types.SimpleNamespace(pos=pos, fiber_offsets=off, n_fibers=off.shape[0] - 1)
    rep = EmergenceDetector(null_mc=1500).detect(pos, off)
    assert bundle_count(net) == 440          # vacuous: "440 filaments => a 440-fiber bundle"
    assert rep.n_emerged_bundles == 0        # honest: nothing condensed out of the isotropic network


def test_detect_network_duck_types_fiber_array() -> None:
    """``detect_network`` runs on anything exposing pos/fiber_offsets (e.g. a FiberNetwork)."""
    pos, off = syn.aligned_network(200, seed=5, jitter_deg=5.0)
    net = types.SimpleNamespace(pos=pos, fiber_offsets=off)
    rep = EmergenceDetector(null_mc=1200).detect_network(net)
    assert rep.is_ordered_global and rep.n_fibers == 200


def test_module_level_detect_matches_class() -> None:
    """The convenience ``detect(...)`` equals ``EmergenceDetector(...).detect(...)`` (same seed/config)."""
    pos, off, mask, center = syn.planted_bundle(300, 30, seed=6)
    a = detect(pos, off)
    b = EmergenceDetector().detect(pos, off)
    assert a.s_global == pytest.approx(b.s_global, abs=1e-12)
    assert a.n_emerged_bundles == b.n_emerged_bundles


def test_report_fields_are_finite_and_typed() -> None:
    """The report carries finite scalars, a unit director, and a populated null band (no NaNs)."""
    pos, off, mask, center = syn.planted_bundle(300, 40, seed=7)
    rep = EmergenceDetector(null_mc=1200).detect(pos, off)
    assert np.isfinite(rep.s_global) and np.isfinite(rep.z_global)
    assert rep.director.shape == (3,) and abs(np.linalg.norm(rep.director) - 1.0) < 1e-9
    assert rep.null_band.n_fibers == rep.n_fibers and rep.null_band.std > 0.0
    assert rep.radius > 0.0
