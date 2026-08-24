"""Deterministic state-to-image observation oracle for the image-AI lane.

The first image milestone is not a generative neural network.  It is a checkable observation operator
that converts point-like simulated labels to microscopy-like channels through a Gaussian PSF while
preserving declared intensity.  Neural encoders, domain randomisation, camera noise, and differentiable
renderers can replace this backend later while retaining the same state/observation contract.
"""

from __future__ import annotations

import hashlib
import json
import math
import re
from collections.abc import Mapping
from dataclasses import dataclass
from types import MappingProxyType
from typing import Any

import numpy as np

__all__ = [
    "CameraConfig",
    "MicroscopyConfig",
    "SyntheticObservationProvenance",
    "apply_camera_model",
    "render_multichannel_points",
    "render_point_channel",
]


@dataclass(frozen=True, slots=True)
class MicroscopyConfig:
    """Geometry and PSF of a deterministic 2-D synthetic microscope."""

    height_px: int
    width_px: int
    pixel_size_um: float
    psf_sigma_um: float
    background: float = 0.0
    psf_radius_sigma: float = 4.0
    world_y_direction: str = "up"

    def __post_init__(self) -> None:
        for name in ("height_px", "width_px"):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, int) or value < 1:
                raise ValueError(f"{name} must be a positive integer")
        for name in ("pixel_size_um", "psf_sigma_um", "psf_radius_sigma"):
            value = getattr(self, name)
            if (
                isinstance(value, bool)
                or not isinstance(value, (int, float))
                or not math.isfinite(value)
                or value <= 0.0
            ):
                raise ValueError(f"{name} must be positive-finite")
        if (
            isinstance(self.background, bool)
            or not isinstance(self.background, (int, float))
            or not math.isfinite(self.background)
            or self.background < 0.0
        ):
            raise ValueError("background must be finite and nonnegative")
        if self.world_y_direction not in {"up", "down"}:
            raise ValueError("world_y_direction must be 'up' or 'down'")

    @property
    def config_hash(self) -> str:
        payload = json.dumps(
            {
                "height_px": self.height_px,
                "width_px": self.width_px,
                "pixel_size_um": self.pixel_size_um,
                "psf_sigma_um": self.psf_sigma_um,
                "background": self.background,
                "psf_radius_sigma": self.psf_radius_sigma,
                "world_y_direction": self.world_y_direction,
                "origin_semantics": (
                    "lower-left-pixel-edge"
                    if self.world_y_direction == "up"
                    else "upper-left-pixel-edge"
                ),
            },
            sort_keys=True,
            separators=(",", ":"),
        ).encode()
        return hashlib.sha256(payload).hexdigest()


@dataclass(frozen=True, slots=True)
class CameraConfig:
    """Photon/read-noise, saturation, and quantization parameters for a CPU camera oracle."""

    photons_per_intensity: float
    read_noise_std_photons: float
    saturation_photons: float
    bit_depth: int = 16

    def __post_init__(self) -> None:
        for name in ("photons_per_intensity", "saturation_photons"):
            value = getattr(self, name)
            if (
                isinstance(value, bool)
                or not isinstance(value, (int, float))
                or not math.isfinite(value)
                or value <= 0.0
            ):
                raise ValueError(f"{name} must be positive-finite")
        if (
            isinstance(self.read_noise_std_photons, bool)
            or not isinstance(self.read_noise_std_photons, (int, float))
            or not math.isfinite(self.read_noise_std_photons)
            or self.read_noise_std_photons < 0.0
        ):
            raise ValueError("read_noise_std_photons must be finite and nonnegative")
        if (
            isinstance(self.bit_depth, bool)
            or not isinstance(self.bit_depth, int)
            or not 1 <= self.bit_depth <= 32
        ):
            raise ValueError("bit_depth must be an integer in [1, 32]")

    @property
    def config_hash(self) -> str:
        payload = json.dumps(
            {
                "photons_per_intensity": self.photons_per_intensity,
                "read_noise_std_photons": self.read_noise_std_photons,
                "saturation_photons": self.saturation_photons,
                "bit_depth": self.bit_depth,
            },
            sort_keys=True,
            separators=(",", ":"),
        ).encode()
        return hashlib.sha256(payload).hexdigest()


_SHA256 = re.compile(r"^[0-9a-f]{64}$")


@dataclass(frozen=True, slots=True)
class SyntheticObservationProvenance:
    """Hash bindings needed to replay one synthetic observation."""

    source_state_artifact_sha256: str
    microscopy_config_sha256: str
    camera_config_sha256: str
    expected_image_sha256: str
    observed_sidecar_sha256: str
    channel_emitter_mapping_sha256: str
    world_to_camera_sha256: str
    renderer_id: str
    rng_algorithm: str
    camera_seed: int

    def __post_init__(self) -> None:
        for field_name in (
            "source_state_artifact_sha256",
            "microscopy_config_sha256",
            "camera_config_sha256",
            "expected_image_sha256",
            "observed_sidecar_sha256",
            "channel_emitter_mapping_sha256",
            "world_to_camera_sha256",
        ):
            if not _SHA256.fullmatch(getattr(self, field_name)):
                raise ValueError(f"{field_name} must be a lowercase SHA-256 digest")
        if not isinstance(self.renderer_id, str) or not self.renderer_id.strip():
            raise ValueError("renderer_id must be a non-empty string")
        if not isinstance(self.rng_algorithm, str) or not self.rng_algorithm.strip():
            raise ValueError("rng_algorithm must be a non-empty string")
        if (
            isinstance(self.camera_seed, bool)
            or not isinstance(self.camera_seed, int)
            or self.camera_seed < 0
        ):
            raise ValueError("camera_seed must be a nonnegative integer")

    @property
    def record_hash(self) -> str:
        payload = json.dumps(
            {
                "source_state_artifact_sha256": self.source_state_artifact_sha256,
                "microscopy_config_sha256": self.microscopy_config_sha256,
                "camera_config_sha256": self.camera_config_sha256,
                "expected_image_sha256": self.expected_image_sha256,
                "observed_sidecar_sha256": self.observed_sidecar_sha256,
                "channel_emitter_mapping_sha256": self.channel_emitter_mapping_sha256,
                "world_to_camera_sha256": self.world_to_camera_sha256,
                "renderer_id": self.renderer_id.strip(),
                "rng_algorithm": self.rng_algorithm.strip(),
                "camera_seed": self.camera_seed,
            },
            sort_keys=True,
            separators=(",", ":"),
        ).encode()
        return hashlib.sha256(payload).hexdigest()


def apply_camera_model(
    expected_intensity: Any,
    config: CameraConfig,
    *,
    seed: int,
) -> np.ndarray:
    """Sample photons/read noise, then saturate and quantize with an explicit replay seed."""
    if not isinstance(config, CameraConfig):
        raise TypeError("config must be a CameraConfig")
    expected = np.asarray(expected_intensity, dtype=np.float64)
    if expected.ndim < 2 or expected.size == 0:
        raise ValueError("expected_intensity must be a non-empty image or image stack")
    if not np.all(np.isfinite(expected)) or np.any(expected < 0.0):
        raise ValueError("expected_intensity must be finite and nonnegative")
    if isinstance(seed, bool) or not isinstance(seed, int) or seed < 0:
        raise ValueError("seed must be a nonnegative integer")

    generator = np.random.default_rng(seed)
    maximum_expected = float(np.max(expected))
    if maximum_expected > np.finfo(np.float64).max / float(config.photons_per_intensity):
        raise ValueError("camera photon expectation overflowed; rescale the declared observation")
    expected_photons = expected * float(config.photons_per_intensity)
    if not np.all(np.isfinite(expected_photons)):
        raise ValueError("camera photon expectation overflowed; rescale the declared observation")
    observed_photons = generator.poisson(expected_photons).astype(np.float64)
    if config.read_noise_std_photons > 0.0:
        observed_photons += generator.normal(
            loc=0.0,
            scale=float(config.read_noise_std_photons),
            size=expected.shape,
        )
    observed_photons = np.clip(observed_photons, 0.0, float(config.saturation_photons))
    digital_max = (1 << config.bit_depth) - 1
    digital = np.rint(observed_photons * (digital_max / float(config.saturation_photons)))
    dtype = np.uint16 if config.bit_depth <= 16 else np.uint32
    result = digital.astype(dtype)
    result.setflags(write=False)
    return result


def _validated_points(points_xy_um: Any, intensities: Any) -> tuple[np.ndarray, np.ndarray]:
    points = np.asarray(points_xy_um, dtype=np.float64)
    weights = np.asarray(intensities, dtype=np.float64)
    if points.ndim != 2 or points.shape[1] != 2:
        raise ValueError(f"points_xy_um must have shape (n, 2), got {points.shape}")
    if weights.ndim != 1 or weights.shape[0] != points.shape[0]:
        raise ValueError(f"intensities must have shape ({points.shape[0]},), got {weights.shape}")
    if not np.all(np.isfinite(points)) or not np.all(np.isfinite(weights)):
        raise ValueError("points and intensities must be finite")
    if np.any(weights < 0.0):
        raise ValueError("intensities must be nonnegative")
    return points, weights


def render_point_channel(
    points_xy_um: Any,
    intensities: Any,
    config: MicroscopyConfig,
    *,
    origin_xy_um: tuple[float, float] = (0.0, 0.0),
) -> np.ndarray:
    """Render point labels through a truncated Gaussian PSF.

    ``origin_xy_um`` is the lower-left image edge when world y points up, and the upper-left edge when
    it points down.  The truncated PSF is normalized before frame clipping, so photons outside the
    field of view are lost rather than renormalized into a bright boundary ghost.
    """
    if not isinstance(config, MicroscopyConfig):
        raise TypeError("config must be a MicroscopyConfig")
    points, weights = _validated_points(points_xy_um, intensities)
    origin = np.asarray(origin_xy_um, dtype=np.float64)
    if origin.shape != (2,) or not np.all(np.isfinite(origin)):
        raise ValueError("origin_xy_um must be a finite (x, y) pair")

    image = np.full(
        (config.height_px, config.width_px),
        float(config.background),
        dtype=np.float64,
    )
    sigma_px = config.psf_sigma_um / config.pixel_size_um
    radius_px = max(1, int(math.ceil(config.psf_radius_sigma * sigma_px)))

    for point, intensity in zip(points, weights, strict=True):
        if intensity == 0.0:
            continue
        physical = (point - origin) / config.pixel_size_um
        cx = float(physical[0] - 0.5)
        cy = (
            float(config.height_px - physical[1] - 0.5)
            if config.world_y_direction == "up"
            else float(physical[1] - 0.5)
        )
        full_x0 = int(math.floor(cx)) - radius_px
        full_x1 = int(math.floor(cx)) + radius_px + 2
        full_y0 = int(math.floor(cy)) - radius_px
        full_y1 = int(math.floor(cy)) + radius_px + 2
        x0 = max(0, full_x0)
        x1 = min(config.width_px, full_x1)
        y0 = max(0, full_y0)
        y1 = min(config.height_px, full_y1)
        if x0 >= x1 or y0 >= y1:
            continue
        xs = np.arange(full_x0, full_x1, dtype=np.float64)
        ys = np.arange(full_y0, full_y1, dtype=np.float64)
        dx2 = np.square(xs - cx)[None, :]
        dy2 = np.square(ys - cy)[:, None]
        full_patch = np.exp(-0.5 * (dx2 + dy2) / (sigma_px * sigma_px))
        patch_sum = float(np.sum(full_patch))
        if patch_sum == 0.0:
            continue
        full_patch *= float(intensity) / patch_sum
        patch = full_patch[
            y0 - full_y0 : y1 - full_y0,
            x0 - full_x0 : x1 - full_x0,
        ]
        image[y0:y1, x0:x1] += patch
    image.setflags(write=False)
    return image


def render_multichannel_points(
    channels: Mapping[str, tuple[Any, Any]],
    config: MicroscopyConfig,
    *,
    origin_xy_um: tuple[float, float] = (0.0, 0.0),
) -> Mapping[str, np.ndarray]:
    """Render independent named point channels with a shared microscope configuration."""
    if not isinstance(channels, Mapping):
        raise TypeError("channels must map channel name to (points_xy_um, intensities)")
    rendered: dict[str, np.ndarray] = {}
    for name, payload in channels.items():
        if not isinstance(name, str) or not name.strip():
            raise ValueError("channel names must be non-empty strings")
        normalized_name = name.strip()
        if normalized_name in rendered:
            raise ValueError(
                f"channel names must not collide after whitespace normalization: {name!r}"
            )
        if not isinstance(payload, tuple) or len(payload) != 2:
            raise ValueError(f"channel {name!r} must contain (points_xy_um, intensities)")
        rendered[normalized_name] = render_point_channel(
            payload[0], payload[1], config, origin_xy_um=origin_xy_um
        )
    return MappingProxyType(rendered)
