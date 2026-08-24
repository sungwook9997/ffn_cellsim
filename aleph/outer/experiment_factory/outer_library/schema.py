"""Canonical, deliberately non-authoritative observation schema."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Any


@dataclass(frozen=True, slots=True)
class Observation:
    observation_id: str
    dataset_id: str
    lab_group: str
    provider: str
    accession: str
    license_id: str
    source_sha256: str
    source_locator: str
    sample_id: str
    biological_replicate: str
    technical_replicate: str
    condition: str
    cell_type: str
    cell_state: str
    modality: str
    observable: str
    value: float
    unit: str
    unit_source_sha256: str = ""
    unit_source_locator: str = ""
    coordinate_name: str = ""
    coordinate_value: float | None = None
    coordinate_unit: str = ""
    label_authority: str = "author_reported"
    aleph_authority: str = "none"

    def as_dict(self) -> dict[str, Any]:
        return asdict(self)

    def validate(self) -> None:
        if self.label_authority != "author_reported":
            raise ValueError("only author-reported labels enter the training manifest")
        if self.aleph_authority != "none":
            raise ValueError("outer-library observations cannot carry Aleph authority")
        required = (
            self.observation_id,
            self.dataset_id,
            self.lab_group,
            self.source_sha256,
            self.source_locator,
            self.sample_id,
            self.condition,
            self.modality,
            self.observable,
            self.unit,
        )
        if any(not str(value).strip() for value in required):
            raise ValueError("traceability fields must be non-empty")
        if bool(self.unit_source_sha256) != bool(self.unit_source_locator):
            raise ValueError("unit source hash and locator must be recorded together")
        if self.coordinate_name:
            if self.coordinate_value is None or not self.coordinate_unit:
                raise ValueError("coordinates require a value and unit")
        elif self.coordinate_value is not None or self.coordinate_unit:
            raise ValueError("coordinate value/unit require a coordinate name")
