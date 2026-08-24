"""Acceptance and CUDA gates for unilateral ERM membrane--cortex tethers."""

from __future__ import annotations

import numpy as np
import pytest

from aleph.components.incumbent.erm_tether import ERMBellKinetics
from aleph.components.incumbent.erm_tether_analytic import erm_tether_energy, erm_tether_tension


def test_unilateral_tether_is_force_free_in_compression_and_at_rest() -> None:
    length = np.array([0.0, 0.9, 1.0, 1.2], dtype=np.float64)
    tension = erm_tether_tension(length, k_erm=10.0, rest_length=1.0)
    energy = erm_tether_energy(length, k_erm=10.0, rest_length=1.0)
    np.testing.assert_array_equal(tension[:3], 0.0)
    np.testing.assert_array_equal(energy[:3], 0.0)
    assert tension[3] == pytest.approx(2.0)
    assert energy[3] == pytest.approx(0.2)


def test_tension_is_positive_extension_energy_gradient() -> None:
    length = 1.2
    step = 1.0e-6
    derivative = (
        erm_tether_energy(length + step, k_erm=10.0, rest_length=1.0)
        - erm_tether_energy(length - step, k_erm=10.0, rest_length=1.0)
    ) / (2.0 * step)
    expected = erm_tether_tension(length, k_erm=10.0, rest_length=1.0)
    assert float(derivative) == pytest.approx(float(expected), rel=1.0e-10)


def test_unilateral_oracle_rejects_nonphysical_inputs() -> None:
    with pytest.raises(ValueError, match="k_erm"):
        erm_tether_tension(1.0, k_erm=0.0, rest_length=1.0)
    with pytest.raises(ValueError, match="rest_length"):
        erm_tether_tension(1.0, k_erm=1.0, rest_length=-1.0)
    with pytest.raises(ValueError, match="length"):
        erm_tether_tension([-0.1], k_erm=1.0, rest_length=0.0)


def test_bell_contract_has_no_unsourced_or_partial_defaults() -> None:
    kinetics = ERMBellKinetics(
        k_on_s=2.0,
        k_off0_s=0.5,
        bell_force_pn=3.0,
        capture_radius_um=0.02,
        source="TEST_ONLY",
    )
    assert kinetics.rebind_rest_policy == "formation_length"
    with pytest.raises(ValueError, match="positive"):
        ERMBellKinetics(0.0, 0.5, 3.0, 0.02, "TEST_ONLY")
    with pytest.raises(ValueError, match="provenance"):
        ERMBellKinetics(2.0, 0.5, 3.0, 0.02, "")
    with pytest.raises(ValueError, match="formation_length"):
        ERMBellKinetics(2.0, 0.5, 3.0, 0.02, "TEST_ONLY", "fixed_rest")


def test_mcf7_production_requires_complete_bell_contract_before_cuda_build() -> None:
    from aleph.components.incumbent.assemble import CellConfig, build_cell

    with pytest.raises(ValueError, match="all-or-none"):
        build_cell(CellConfig(erm_k_on_s=1.0))
    with pytest.raises(ValueError, match="complete source-gated Bell"):
        build_cell(CellConfig(
            erm_density_per_um2=10.0,
            erm_density_source="TEST_ONLY",
            erm_density_mcf7_production=True,
        ))


def test_cuda_force_and_commit_contract() -> None:
    """Compression/rest are force-free; tension is Newton-3rd; only accepted overload ruptures."""
    wp = pytest.importorskip("warp")
    from aleph.components.incumbent.erm_tether import erm_tether_force_kernel, erm_tether_rupture_kernel

    wp.init()
    cuda_device = next((device for device in wp.get_devices() if device.is_cuda), None)
    if cuda_device is None:
        pytest.skip("I0-A: ERM kernel gate requires a CUDA GPU")
    device = str(cuda_device)

    # Four disjoint pairs: compressed, resting, tensile below rupture, tensile above rupture.
    pos = np.array([
        [0.9, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [1.2, 0.0, 0.0],
        [1.5, 0.0, 0.0],
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
    ], dtype=np.float64)
    membrane_idx = np.arange(4, dtype=np.int32)
    cortex_idx = np.arange(4, 8, dtype=np.int32)
    rest = np.ones(4, dtype=np.float64)
    bound_d = wp.ones(4, dtype=wp.int32, device=device)
    pos_d = wp.array(pos, dtype=wp.vec3d, device=device)
    membrane_d = wp.array(membrane_idx, dtype=wp.int32, device=device)
    cortex_d = wp.array(cortex_idx, dtype=wp.int32, device=device)
    rest_d = wp.array(rest, dtype=wp.float64, device=device)
    force_d = wp.zeros(8, dtype=wp.vec3d, device=device)

    wp.launch(
        erm_tether_force_kernel,
        dim=4,
        inputs=[pos_d, membrane_d, cortex_d, bound_d, wp.float64(10.0), rest_d, wp.float64(3.0)],
        outputs=[force_d],
        device=device,
    )
    measured = force_d.numpy()
    np.testing.assert_array_equal(measured[[0, 1, 3, 4, 5, 7]], 0.0)
    np.testing.assert_allclose(measured[2], [-2.0, 0.0, 0.0], rtol=0.0, atol=2.0e-15)
    np.testing.assert_allclose(measured[6], [2.0, 0.0, 0.0], rtol=0.0, atol=2.0e-15)
    np.testing.assert_allclose(measured.sum(axis=0), 0.0, rtol=0.0, atol=2.0e-15)

    accepted_d = wp.zeros(1, dtype=wp.int32, device=device)
    wp.launch(
        erm_tether_rupture_kernel,
        dim=4,
        inputs=[
            pos_d, membrane_d, cortex_d, bound_d, wp.float64(10.0), rest_d, wp.float64(3.0), accepted_d,
        ],
        device=device,
    )
    np.testing.assert_array_equal(bound_d.numpy(), 1)
    accepted_d.fill_(1)
    wp.launch(
        erm_tether_rupture_kernel,
        dim=4,
        inputs=[
            pos_d, membrane_d, cortex_d, bound_d, wp.float64(10.0), rest_d, wp.float64(3.0), accepted_d,
        ],
        device=device,
    )
    np.testing.assert_array_equal(bound_d.numpy(), [1, 1, 1, 0])


def test_cuda_bell_force_and_transactional_on_off_kmc() -> None:
    """Bell mode carries full load, rejection is exact, and accepted detach/rebind uses formation length."""
    wp = pytest.importorskip("warp")
    from aleph.components.incumbent.erm_tether import (
        erm_bell_force_kernel,
        erm_bell_kmc_kernel,
        increment_erm_epoch_if_accepted_kernel,
    )

    wp.init()
    cuda_device = next((device for device in wp.get_devices() if device.is_cuda), None)
    if cuda_device is None:
        pytest.skip("I0-A: ERM Bell kernel gate requires a CUDA GPU")
    device = str(cuda_device)

    # Pair 0 is bound and strongly tensile. Pair 1 is unbound inside capture. Pair 2 is unbound outside.
    pos = np.array([
        [2.0, 0.0, 0.0], [0.2, 0.0, 0.0], [1.0, 0.0, 0.0],
        [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0],
    ], dtype=np.float64)
    membrane_d = wp.array(np.arange(3, dtype=np.int32), dtype=wp.int32, device=device)
    cortex_d = wp.array(np.arange(3, 6, dtype=np.int32), dtype=wp.int32, device=device)
    pos_d = wp.array(pos, dtype=wp.vec3d, device=device)
    bound_d = wp.array(np.array([1, 0, 0], dtype=np.int32), dtype=wp.int32, device=device)
    rest_d = wp.array(np.ones(3, dtype=np.float64), dtype=wp.float64, device=device)
    force_d = wp.zeros(6, dtype=wp.vec3d, device=device)
    detach_d = wp.zeros(1, dtype=wp.int32, device=device)
    attach_d = wp.zeros(1, dtype=wp.int32, device=device)
    rng_epoch_d = wp.zeros(1, dtype=wp.int32, device=device)

    # Full 10 pN load remains present; Bell detachment is a later physical-time event, not a force cutoff.
    wp.launch(
        erm_bell_force_kernel,
        dim=3,
        inputs=[pos_d, membrane_d, cortex_d, bound_d, wp.float64(10.0), rest_d],
        outputs=[force_d],
        device=device,
    )
    force = force_d.numpy()
    np.testing.assert_allclose(force[0], [-10.0, 0.0, 0.0], rtol=0.0, atol=1.0e-14)
    np.testing.assert_allclose(force[3], [10.0, 0.0, 0.0], rtol=0.0, atol=1.0e-14)
    np.testing.assert_allclose(force.sum(axis=0), 0.0, rtol=0.0, atol=1.0e-14)

    accepted_d = wp.zeros(1, dtype=wp.int32, device=device)
    kernel_inputs = [
        pos_d, membrane_d, cortex_d, bound_d, wp.float64(10.0), rest_d,
        wp.float64(1.0e6), wp.float64(1.0e6), wp.float64(1.0), wp.float64(0.5),
        wp.float64(1.0), wp.int32(1234), rng_epoch_d, accepted_d, detach_d, attach_d,
    ]
    wp.launch(erm_bell_kmc_kernel, dim=3, inputs=kernel_inputs, device=device)
    wp.launch(
        increment_erm_epoch_if_accepted_kernel,
        dim=1,
        inputs=[accepted_d, rng_epoch_d],
        device=device,
    )
    np.testing.assert_array_equal(bound_d.numpy(), [1, 0, 0])
    np.testing.assert_array_equal(rest_d.numpy(), 1.0)
    assert int(detach_d.numpy()[0]) == 0
    assert int(attach_d.numpy()[0]) == 0
    assert int(rng_epoch_d.numpy()[0]) == 0

    accepted_d.fill_(1)
    wp.launch(erm_bell_kmc_kernel, dim=3, inputs=kernel_inputs, device=device)
    wp.launch(
        increment_erm_epoch_if_accepted_kernel,
        dim=1,
        inputs=[accepted_d, rng_epoch_d],
        device=device,
    )
    np.testing.assert_array_equal(bound_d.numpy(), [0, 1, 0])
    np.testing.assert_allclose(rest_d.numpy(), [1.0, 0.2, 1.0], rtol=0.0, atol=1.0e-15)
    assert int(detach_d.numpy()[0]) == 1
    assert int(attach_d.numpy()[0]) == 1
    assert int(rng_epoch_d.numpy()[0]) == 1
