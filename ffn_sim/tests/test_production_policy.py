"""Unit tests for executable production-run guardrails."""

from __future__ import annotations

from types import SimpleNamespace

import pytest

from ffn_sim.common.production_policy import (
    DEFAULT_MCF7_TURGOR_PA,
    require_full_cell_physiological_baseline,
    require_physiological_cytoplasm,
    require_physiological_turgor,
    require_production_device,
)


def _hoomd_stub(*, gpu_enabled: bool):
    return SimpleNamespace(version=SimpleNamespace(gpu_enabled=gpu_enabled))


def _cyto(*, eta_eff=65.9, eta_water=6.913e-4, cell_type="MCF7", is_noop=False):
    return SimpleNamespace(
        eta_eff=eta_eff,
        eta_water=eta_water,
        cell_type=cell_type,
        is_noop=is_noop,
    )


def _ev(*, turgor=DEFAULT_MCF7_TURGOR_PA):
    return SimpleNamespace(turgor_dP0=turgor)


def test_device_policy_requires_gpu_or_explicit_cpu_dev():
    with pytest.raises(ValueError, match="dev-only"):
        require_production_device(
            "cpu", allow_cpu_dev=False, hoomd_module=_hoomd_stub(gpu_enabled=False)
        )
    require_production_device(
        "cpu", allow_cpu_dev=True, hoomd_module=_hoomd_stub(gpu_enabled=False)
    )
    with pytest.raises(ValueError, match="gpu_enabled=False"):
        require_production_device(
            "gpu", allow_cpu_dev=False, hoomd_module=_hoomd_stub(gpu_enabled=False)
        )
    require_production_device(
        "gpu", allow_cpu_dev=False, hoomd_module=_hoomd_stub(gpu_enabled=True)
    )


def test_cytoplasm_policy_rejects_missing_or_water_drag():
    with pytest.raises(ValueError, match="p_cytoplasm=None"):
        require_physiological_cytoplasm(None)
    with pytest.raises(ValueError, match="water/no-op"):
        require_physiological_cytoplasm(_cyto(eta_eff=6.913e-4, cell_type="water", is_noop=True))
    require_physiological_cytoplasm(_cyto())


def test_turgor_policy_rejects_missing_zero_or_none_pressure():
    with pytest.raises(ValueError, match="positive"):
        require_physiological_turgor(None)
    with pytest.raises(ValueError, match="baseline osmotic turgor"):
        require_physiological_turgor(_ev(turgor=0.0))
    with pytest.raises(ValueError, match="baseline osmotic turgor"):
        require_physiological_turgor(_ev(turgor=None))
    require_physiological_turgor(_ev())


def test_full_cell_policy_requires_both_cytoplasm_and_turgor():
    require_full_cell_physiological_baseline(
        {"p_cytoplasm": _cyto(), "p_enclosed_volume": _ev()}
    )
    with pytest.raises(ValueError, match="cytoplasm"):
        require_full_cell_physiological_baseline(
            {"p_cytoplasm": None, "p_enclosed_volume": _ev()}
        )
    with pytest.raises(ValueError, match="turgor"):
        require_full_cell_physiological_baseline(
            {"p_cytoplasm": _cyto(), "p_enclosed_volume": _ev(turgor=0.0)}
        )
