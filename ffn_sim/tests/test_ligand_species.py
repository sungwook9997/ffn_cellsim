"""Static tests for bridge/ligand_species.py (ECM ligand identity → clutch kinetics).

Verifies the literature-anchored ligand registry and the Bell-Evans↔Pereverzev
mapping WITHOUT building a simulation (pure-Python kinetics math — fast, no HOOMD).
This module is import-isolated (no runtime file imports ligand_species), so these
tests cannot perturb the live build; they regression-guard the prep module against
the values in docs/LIGAND_IDENTITY_FA_SPEC.md.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from ffn_sim.archive.hoomd_legacy.bridge.ligand_species import (
    CATCH_DISABLED_FRACTION,
    DEFAULT_LIGAND_FOR_CONDITION,
    K_BT,
    LIGAND_REGISTRY,
    bell_evans_k_off,
    catch_peak_force_for,
    pereverzev_params_for,
)
from ffn_sim.validation.pereverzev import PereverzevParams, pereverzev_k_off

PN = 1.0e-12  # one piconewton, in Newtons

SLIP_SPECIES = ["collagen_I", "collagen_I_locked_open", "laminin_111"]
ALL_SPECIES = SLIP_SPECIES + ["fibronectin"]


def test_registry_constructs_and_is_complete():
    """Every registry entry is a valid LigandSpecies; conditions map to real keys."""
    assert set(ALL_SPECIES) <= set(LIGAND_REGISTRY)
    for key, lig in LIGAND_REGISTRY.items():
        assert lig.name == key
        assert lig.bond_class in ("slip", "catch-slip")
        assert lig.integrin.startswith("alpha")
    for cond, key in DEFAULT_LIGAND_FOR_CONDITION.items():
        assert key in LIGAND_REGISTRY, f"condition {cond!r} → unknown ligand {key!r}"


@pytest.mark.parametrize(
    "species, expected_pN",
    [
        ("collagen_I", 18.6),          # x_β = 0.23 nm
        ("collagen_I_locked_open", 6.1),  # x_β = 0.70 nm
        ("laminin_111", 15.3),         # x_β = 0.28 nm
    ],
)
def test_F_s_equals_kT_over_xbeta(species, expected_pN):
    """F_s = k_BT / x_β exactly, and lands at the documented pN value."""
    lig = LIGAND_REGISTRY[species]
    assert lig.x_beta is not None
    assert lig.F_s == pytest.approx(K_BT / lig.x_beta, rel=1e-12)
    assert lig.F_s / PN == pytest.approx(expected_pN, abs=0.5)


def test_laminin_F_s_matches_reported_f_b():
    """Internal cross-check: the proxy F_s ≈ reported f_b 14.1±1.3 pN (α7β1-invasin)."""
    lig = LIGAND_REGISTRY["laminin_111"]
    assert lig.F_s / PN == pytest.approx(14.1, abs=2.0)


@pytest.mark.parametrize("species", SLIP_SPECIES)
def test_slip_maps_to_pereverzev_with_catch_disabled(species):
    """Slip ligands → PereverzevParams with k_s=k_off0, F_s=kT/x_β, catch disabled."""
    lig = LIGAND_REGISTRY[species]
    p = pereverzev_params_for(species)
    assert isinstance(p, PereverzevParams)
    assert p.k_s == pytest.approx(lig.k_off0)
    assert p.F_s == pytest.approx(lig.F_s)
    # catch pathway negligible (k_c/k_s == the disabled fraction)
    assert p.k_c / p.k_s == pytest.approx(CATCH_DISABLED_FRACTION, rel=1e-9)


@pytest.mark.parametrize("species", SLIP_SPECIES)
def test_slip_has_no_catch_peak(species):
    """A slip bond has no catch peak → F* is NaN (KU-2.5 gate must skip slip)."""
    assert math.isnan(catch_peak_force_for(species))


def test_fibronectin_is_full_catch_slip_with_peak():
    """FN keeps the KU-2.18 catch-slip default and a finite catch peak ~7 pN."""
    lig = LIGAND_REGISTRY["fibronectin"]
    assert lig.bond_class == "catch-slip"
    p = pereverzev_params_for("fibronectin")
    assert p is lig.catch_params
    F_star = catch_peak_force_for("fibronectin")
    assert math.isfinite(F_star)
    assert 1.0 < F_star / PN < 30.0  # catch-slip transition in the tensile regime


@pytest.mark.parametrize("species", SLIP_SPECIES)
def test_bell_evans_matches_analytic_slip(species):
    """bell_evans_k_off reproduces k_off0·exp(F·x_β/kT); equals k_off0 at F=0."""
    lig = LIGAND_REGISTRY[species]
    assert bell_evans_k_off(0.0, species) == pytest.approx(lig.k_off0)
    for F_pN in (5.0, 20.0, 50.0):
        F = F_pN * PN
        expected = lig.k_off0 * math.exp(F * lig.x_beta / K_BT)
        assert bell_evans_k_off(F, species) == pytest.approx(expected, rel=1e-9)


@pytest.mark.parametrize("species", SLIP_SPECIES)
def test_bell_evans_consistent_with_pereverzev_mapping(species):
    """The mapped Pereverzev k_off equals the Bell-Evans slip (catch term negligible)."""
    p = pereverzev_params_for(species)
    for F_pN in (0.0, 10.0, 30.0):
        F = F_pN * PN
        assert pereverzev_k_off(F, p) == pytest.approx(
            bell_evans_k_off(F, species), rel=1e-6
        )


def test_bell_evans_rejects_catch_slip():
    """bell_evans_k_off is slip-only; FN (catch-slip) must raise."""
    with pytest.raises(ValueError):
        bell_evans_k_off(10.0 * PN, "fibronectin")


@pytest.mark.parametrize("F_pN", [5.0, 10.0, 20.0])
def test_laminin_weaker_than_fibronectin(F_pN):
    """Rel-traction anchor: LN-111 (α6β1) clutch unbinds FASTER than FN (α5β1)."""
    F = F_pN * PN
    k_lam = pereverzev_k_off(F, pereverzev_params_for("laminin_111"))
    k_fn = pereverzev_k_off(F, pereverzev_params_for("fibronectin"))
    assert float(k_lam) > float(k_fn), (
        f"laminin should unbind faster (weaker clutch) at {F_pN} pN: "
        f"k_lam={float(k_lam):.3g} vs k_fn={float(k_fn):.3g}"
    )


def test_proxy_flags():
    """Only laminin (direct kinetics absent in lit) is proxy-flagged."""
    assert LIGAND_REGISTRY["laminin_111"].proxy is True
    assert LIGAND_REGISTRY["collagen_I"].proxy is False
    assert LIGAND_REGISTRY["fibronectin"].proxy is False
    # confidence strings present
    for lig in LIGAND_REGISTRY.values():
        assert lig.confidence and lig.source


def test_kT_value():
    """k_BT at 310.15 K ≈ 4.2816e-21 J (anchors the F_s conversion)."""
    assert K_BT == pytest.approx(4.2816e-21, rel=1e-3)
