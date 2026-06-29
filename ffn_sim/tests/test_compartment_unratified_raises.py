"""Platform SAFETY-CONTRACT canary: un-ratified paths MUST refuse to run.

The companion of ``test_compartment_enabled_canary.py``. The positive canary
proves each compartment's enabled path runs *with* SMOKE-ONLY candidate constants;
this negative canary proves the PRODUCTION path *refuses* (raises) when a
None-gated / PI-pending constant is absent — the no-magic-number "honesty over
completeness" guarantee the whole platform rests on. A refactor that silently
makes an un-ratified path runnable (re-introducing an invented placeholder) is
caught here. One place, platform-wide.
"""

from __future__ import annotations

import numpy as np
import pytest

import hoomd.md as md

_KT = 4.28e-21


def test_stress_fibers_backbone_halts_without_mu_SF_or_k_actin():
    from ffn_sim.archive.hoomd_legacy.cell.stress_fibers import (
        resolve_backbone_k_bond,
        resolve_stress_fibers,
    )

    p = resolve_stress_fibers(
        {"stress_fibers": {"enabled": True, "n_SF": 1, "n_beads_per_SF": 12}}
    )
    assert p.mu_SF is None and p.k_actin is None
    with pytest.raises(NotImplementedError):
        resolve_backbone_k_bond(p, 1.0e-7)


def test_linc_bond_potential_halts_without_k_linc():
    from ffn_sim.archive.hoomd_legacy.cell.linc import configure_linc_bond_potential, resolve_linc

    p = resolve_linc({"linc": {"enabled": True}}, R_cell=7.5e-6, R_nuc=3.0e-6)
    assert p.k_linc is None
    with pytest.raises(NotImplementedError):
        configure_linc_bond_potential(md.bond.Harmonic(), p)


def test_intermediate_filaments_nonlinear_law_not_faked():
    from ffn_sim.archive.hoomd_legacy.cell.intermediate_filaments import (
        register_if_bond_params,
        resolve_intermediate_filaments,
    )

    p = resolve_intermediate_filaments(
        {"intermediate_filaments": {"enabled": True, "n_filaments": 4,
                                    "beads_per_fil": 4}},
        kT=_KT, R_cell=7.5e-6,
    )
    with pytest.raises(NotImplementedError):
        register_if_bond_params(md.bond.Harmonic(), p, nonlinear=True)


def test_microtubules_di_updater_unimplemented():
    from ffn_sim.archive.hoomd_legacy.cell.microtubules import MTDynamicInstability, resolve_microtubules

    p = resolve_microtubules(
        {"microtubules": {"enabled": True, "n_mt": 4, "beads_per_mt": 5,
                          "L_mt": 2.0e-6, "dynamic_instability": True}},
        kT=_KT, gamma_b=1.0e-5,
    )
    with pytest.raises(NotImplementedError):
        MTDynamicInstability(p, dt=1.0e-9)


def test_microtubules_require_L_mt_host_geometry():
    from ffn_sim.archive.hoomd_legacy.cell.microtubules import resolve_microtubules

    with pytest.raises(ValueError):
        resolve_microtubules(
            {"microtubules": {"enabled": True, "n_mt": 4, "beads_per_mt": 5}},
            kT=_KT, gamma_b=1.0e-5,
        )


def test_membrane_reservoir_bleb_updater_unconditionally_blocked():
    from ffn_sim.archive.hoomd_legacy.cell.membrane_reservoir import (
        MembraneTetherLayout,
        MembraneTetherUpdater,
        resolve_membrane_reservoir,
    )

    p = resolve_membrane_reservoir(
        {"membrane_reservoir": {"enabled": True}}, R_cell=7.5e-6
    )
    assert p.sigma_crit_bleb is None and p.f_excess is None
    layout = MembraneTetherLayout(
        tether_pairs=np.empty((0, 2), dtype=np.int64),
        tether_r0=np.empty((0,), dtype=np.float64),
        n_tether=0, rupture_force=0.0,
    )
    with pytest.raises(NotImplementedError):
        MembraneTetherUpdater(p, layout, kT=_KT)


def test_cadherin_junction_no_invented_cadherin_count():
    from ffn_sim.archive.hoomd_legacy.junction.cadherin import resolve_cadherin_junction

    with pytest.raises(ValueError):
        resolve_cadherin_junction(
            {"junction": {"cadherin": {"enabled": True}}},
            kT=_KT, dt=6.98e-7, contact_zone_width=1.0e-8,
        )


def test_junctional_actin_stub_blocks_build_and_laws():
    from ffn_sim.archive.hoomd_legacy.junction.junctional_actin import (
        catch_off_rate,
        coupling_force_magnitude,
        extend_snapshot_with_junctional_actin,
        resolve_junctional_actin,
    )

    p = resolve_junctional_actin(
        {"junction": {"junctional_actin": {"enabled": True}}},
        kT=_KT, dt=6.98e-7, n_cadherin=10,
    )
    assert not getattr(p, "is_anchored", True)
    with pytest.raises(NotImplementedError):
        coupling_force_magnitude(p, 1.0e-9)
    with pytest.raises(NotImplementedError):
        catch_off_rate(p, 5.0e-12)
    with pytest.raises(NotImplementedError):
        extend_snapshot_with_junctional_actin(object(), p)
