"""One CPU-only cross-lane rehearsal without claiming native or biological evidence."""

from __future__ import annotations

import hashlib
import io

import numpy as np
import pytest

from aleph.virtual_cell.contracts import (
    ArtifactEnvelope,
    CellStateManifest,
    EvidenceSource,
    RepresentationDescriptor,
    RepresentationKind,
    UncertaintySource,
)
from aleph.virtual_cell.observation import SyntheticObservationArtifact
from aleph.virtual_cell.reduction import select_array_representation
from aleph.virtual_cell.sidecar import ArrayAxis, ArraySidecarDescriptor
from aleph.virtual_cell.synthetic_microscopy import (
    CameraConfig,
    MicroscopyConfig,
    SyntheticObservationProvenance,
    apply_camera_model,
    render_point_channel,
)


def test_manifest_to_reduced_state_to_synthetic_observation_keeps_provenance() -> None:
    manifest = CellStateManifest(
        manifest_id="generic-contractile-oracle",
        species="Homo sapiens",
        lineage="generic-contractile",
        identity={"source": "analytic fixture"},
        biological_state={"regime": "pulsatile"},
        environment={"geometry": "2D sandbox"},
        components=("membrane", "cortex"),
        connectors=("membrane_cortex",),
        observations=("synthetic_actin_channel",),
        missing_components=("nucleus",),
        missing_laws=("cell-type-specific-regulation",),
        provenance=("analytic-oracle",),
    )

    state_tensor = np.einsum(
        "t,n,d->tnd",
        np.linspace(1.0, 2.0, 5),
        np.linspace(0.5, 1.5, 12),
        np.array([1.0, 2.0, 3.0]),
    )
    reduction = select_array_representation(
        state_tensor,
        relative_tolerance=1e-12,
    )
    assert reduction.kind == "tensor-train"

    microscopy_config = MicroscopyConfig(
        height_px=24,
        width_px=24,
        pixel_size_um=0.25,
        psf_sigma_um=0.4,
    )
    expected_image = render_point_channel(
        [[2.0, 2.0], [4.0, 3.0]],
        [1.0, 2.0],
        microscopy_config,
    )
    camera_config = CameraConfig(
        photons_per_intensity=100.0,
        read_noise_std_photons=1.0,
        saturation_photons=100.0,
        bit_depth=12,
    )
    observed_image = apply_camera_model(
        expected_image,
        camera_config,
        seed=23,
    )
    buffer = io.BytesIO()
    np.save(buffer, observed_image, allow_pickle=False)
    payload = buffer.getvalue()
    sidecar = ArraySidecarDescriptor(
        uri="objects/synthetic-actin.npy",
        sha256=hashlib.sha256(payload).hexdigest(),
        media_type="application/x-npy",
        schema_id="ffn.virtual-cell/synthetic-image@1",
        dtype=observed_image.dtype.str,
        value_semantics="camera digital number",
        value_unit="ADU",
        axes=(
            ArrayAxis(name="y", length=24, semantics="image row", unit="pixel"),
            ArrayAxis(name="x", length=24, semantics="image column", unit="pixel"),
        ),
        byte_size=len(payload),
    )
    np.testing.assert_array_equal(sidecar.verify_npy_array(payload), observed_image)
    observation_provenance = SyntheticObservationProvenance(
        source_state_artifact_sha256="e" * 64,
        microscopy_config_sha256=microscopy_config.config_hash,
        camera_config_sha256=camera_config.config_hash,
        expected_image_sha256=hashlib.sha256(expected_image.tobytes(order="C")).hexdigest(),
        observed_sidecar_sha256=sidecar.sha256,
        channel_emitter_mapping_sha256="f" * 64,
        world_to_camera_sha256="a" * 64,
        renderer_id="gaussian-point-psf@1",
        rng_algorithm="numpy-pcg64",
        camera_seed=23,
    )
    assert len(observation_provenance.record_hash) == 64

    envelope = ArtifactEnvelope(
        artifact_id="synthetic-observation-oracle",
        cell_state_manifest_hash=manifest.manifest_hash,
        evidence_source=EvidenceSource.SYNTHETIC,
        representation=RepresentationDescriptor(
            kind=RepresentationKind.SYNTHETIC_IMAGE,
            uncertainty_sources=(
                UncertaintySource.MEASUREMENT,
                UncertaintySource.REPRESENTATION,
            ),
            axes=("y", "x"),
        ),
        source_artifact_ids=("analytic-state-fixture",),
        observation_provenance_sha256=observation_provenance.record_hash,
    )
    assert envelope.evidence_source is EvidenceSource.SYNTHETIC
    composite = SyntheticObservationArtifact(
        envelope=envelope,
        sidecar=sidecar,
        provenance=observation_provenance,
        manifest=manifest,
    )
    assert len(composite.artifact_hash) == 64


def test_cross_lane_composite_rejects_sidecar_swap() -> None:
    manifest = CellStateManifest(
        manifest_id="swap-test",
        species="Homo sapiens",
        lineage="generic",
        identity={"fixture": True},
        biological_state={"state": "oracle"},
        environment={"geometry": "2D"},
        components=("membrane",),
        connectors=(),
    )
    microscopy = MicroscopyConfig(8, 8, 0.5, 0.5)
    camera = CameraConfig(10.0, 0.0, 20.0)
    payload = b"payload-a"
    sidecar = ArraySidecarDescriptor(
        uri="objects/a.npy",
        sha256=hashlib.sha256(payload).hexdigest(),
        media_type="application/x-npy",
        schema_id="ffn.virtual-cell/synthetic-image@1",
        dtype="<u2",
        value_semantics="camera digital number",
        value_unit="ADU",
        axes=(
            ArrayAxis("y", 8, "row", "pixel"),
            ArrayAxis("x", 8, "column", "pixel"),
        ),
        byte_size=len(payload),
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
        camera_seed=0,
    )
    envelope = ArtifactEnvelope(
        artifact_id="swap-test",
        cell_state_manifest_hash=manifest.manifest_hash,
        evidence_source=EvidenceSource.SYNTHETIC,
        representation=RepresentationDescriptor(
            kind=RepresentationKind.SYNTHETIC_IMAGE,
            uncertainty_sources=(UncertaintySource.REPRESENTATION,),
            axes=("y", "x"),
        ),
        source_artifact_ids=("state-a",),
        observation_provenance_sha256=provenance.record_hash,
    )
    with pytest.raises(ValueError, match="sidecar"):
        SyntheticObservationArtifact(envelope, sidecar, provenance, manifest)
