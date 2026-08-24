"""Deterministic observation-operator tests for synthetic microscopy."""

from __future__ import annotations

import numpy as np
import pytest

from aleph.virtual_cell.synthetic_microscopy import (
    CameraConfig,
    MicroscopyConfig,
    SyntheticObservationProvenance,
    apply_camera_model,
    render_multichannel_points,
    render_point_channel,
)


def _config() -> MicroscopyConfig:
    return MicroscopyConfig(
        height_px=32,
        width_px=32,
        pixel_size_um=0.25,
        psf_sigma_um=0.4,
        background=0.0,
    )


def _centroid_x(image: np.ndarray) -> float:
    xs = np.arange(image.shape[1], dtype=np.float64)
    return float(np.sum(image * xs[None, :]) / np.sum(image))


def test_point_renderer_conserves_emitted_intensity_away_from_frame_edges() -> None:
    image = render_point_channel([[4.0, 4.0]], [2.0], _config())
    assert float(np.sum(image)) == pytest.approx(2.0, rel=1e-12, abs=1e-12)
    assert image.flags.writeable is False


def test_frame_edge_loses_out_of_view_photons_instead_of_creating_ghost() -> None:
    config = _config()
    edge = render_point_channel([[0.05, 4.0]], [1.0], config)
    outside = render_point_channel([[-1.75, 4.0]], [1.0], config)
    assert 0.0 < float(np.sum(edge)) < 1.0
    assert 0.0 < float(np.sum(outside)) < 1e-4


def test_translation_in_physical_space_translates_image_centroid() -> None:
    config = _config()
    left = render_point_channel([[3.0, 4.0]], [1.0], config)
    right = render_point_channel([[4.0, 4.0]], [1.0], config)
    expected_shift_px = 1.0 / config.pixel_size_um
    assert _centroid_x(right) - _centroid_x(left) == pytest.approx(expected_shift_px, abs=0.03)


def test_multichannel_renderer_keeps_channels_independent() -> None:
    rendered = render_multichannel_points(
        {
            "nucleus": (np.array([[2.0, 2.0]]), np.array([5.0])),
            "actin": (np.array([[5.0, 5.0]]), np.array([7.0])),
        },
        _config(),
    )
    assert set(rendered) == {"nucleus", "actin"}
    assert float(np.sum(rendered["nucleus"])) == pytest.approx(5.0)
    assert float(np.sum(rendered["actin"])) == pytest.approx(7.0)
    assert not np.array_equal(rendered["nucleus"], rendered["actin"])


def test_world_y_up_maps_larger_y_to_smaller_image_row() -> None:
    config = _config()
    lower = render_point_channel([[4.0, 3.0]], [1.0], config)
    upper = render_point_channel([[4.0, 4.0]], [1.0], config)
    rows = np.arange(config.height_px, dtype=np.float64)
    lower_row = float(np.sum(lower * rows[:, None]) / np.sum(lower))
    upper_row = float(np.sum(upper * rows[:, None]) / np.sum(upper))
    assert lower_row - upper_row == pytest.approx(
        1.0 / config.pixel_size_um,
        abs=0.03,
    )


def test_out_of_frame_point_contributes_nothing() -> None:
    image = render_point_channel([[1_000.0, 1_000.0]], [2.0], _config())
    assert float(np.sum(image)) == 0.0


@pytest.mark.parametrize(
    ("points", "intensities", "message"),
    [
        ([[1.0, 2.0, 3.0]], [1.0], "shape"),
        ([[1.0, 2.0]], [1.0, 2.0], "intensities"),
        ([[1.0, 2.0]], [-1.0], "nonnegative"),
    ],
)
def test_invalid_point_channels_are_rejected(points, intensities, message) -> None:
    with pytest.raises(ValueError, match=message):
        render_point_channel(points, intensities, _config())


def test_camera_model_is_seed_replayable() -> None:
    expected = np.full((8, 9), 0.4)
    config = CameraConfig(
        photons_per_intensity=100.0,
        read_noise_std_photons=2.0,
        saturation_photons=100.0,
        bit_depth=12,
    )
    first = apply_camera_model(expected, config, seed=17)
    second = apply_camera_model(expected, config, seed=17)
    np.testing.assert_array_equal(first, second)
    assert first.dtype == np.uint16
    assert not first.flags.writeable


def test_zero_signal_and_zero_read_noise_stays_zero() -> None:
    expected = np.zeros((4, 5))
    observed = apply_camera_model(
        expected,
        CameraConfig(
            photons_per_intensity=10.0,
            read_noise_std_photons=0.0,
            saturation_photons=20.0,
            bit_depth=8,
        ),
        seed=2,
    )
    assert np.count_nonzero(observed) == 0


def test_camera_saturation_is_bounded_by_bit_depth() -> None:
    expected = np.full((100, 100), 1e6)
    observed = apply_camera_model(
        expected,
        CameraConfig(
            photons_per_intensity=1.0,
            read_noise_std_photons=0.0,
            saturation_photons=10.0,
            bit_depth=8,
        ),
        seed=4,
    )
    assert np.all(observed == 255)


def test_camera_rejects_negative_expected_intensity() -> None:
    with pytest.raises(ValueError, match="nonnegative"):
        apply_camera_model(
            [[0.0, -1.0]],
            CameraConfig(
                photons_per_intensity=1.0,
                read_noise_std_photons=0.0,
                saturation_photons=10.0,
            ),
            seed=0,
        )


def test_camera_rejects_photon_expectation_overflow() -> None:
    with pytest.raises(ValueError, match="overflowed"):
        apply_camera_model(
            np.full((2, 2), 1e308),
            CameraConfig(
                photons_per_intensity=1e308,
                read_noise_std_photons=0.0,
                saturation_photons=10.0,
            ),
            seed=0,
        )


def test_observation_provenance_binds_optics_camera_seed_and_images() -> None:
    microscopy = _config()
    camera = CameraConfig(
        photons_per_intensity=10.0,
        read_noise_std_photons=1.0,
        saturation_photons=20.0,
    )
    provenance = SyntheticObservationProvenance(
        source_state_artifact_sha256="a" * 64,
        microscopy_config_sha256=microscopy.config_hash,
        camera_config_sha256=camera.config_hash,
        expected_image_sha256="b" * 64,
        observed_sidecar_sha256="c" * 64,
        channel_emitter_mapping_sha256="d" * 64,
        world_to_camera_sha256="e" * 64,
        renderer_id="gaussian-point-psf@1",
        rng_algorithm="numpy-pcg64",
        camera_seed=19,
    )
    assert len(provenance.record_hash) == 64


def test_multichannel_names_cannot_collide_after_normalization() -> None:
    with pytest.raises(ValueError, match="collide"):
        render_multichannel_points(
            {
                "actin": (np.array([[2.0, 2.0]]), np.array([1.0])),
                " actin ": (np.array([[3.0, 3.0]]), np.array([1.0])),
            },
            _config(),
        )
