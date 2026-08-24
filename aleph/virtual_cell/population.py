"""Heterotypic multicell population and interaction-graph manifests."""

from __future__ import annotations

import hashlib
import json
import math
import re
from collections.abc import Mapping
from dataclasses import dataclass, field
from types import MappingProxyType

from aleph.virtual_cell.contracts import CellStateManifest, JSONScalar

__all__ = [
    "CellInstanceManifest",
    "CellInteraction",
    "CellPopulationManifest",
]

_SHA256 = re.compile(r"^[0-9a-f]{64}$")


def _nonempty(value: str, *, what: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{what} must be a non-empty string")
    return value.strip()


def _scalar_mapping(
    values: Mapping[str, JSONScalar],
    *,
    what: str,
) -> Mapping[str, JSONScalar]:
    if not isinstance(values, Mapping):
        raise TypeError(f"{what} must be a mapping")
    result: dict[str, JSONScalar] = {}
    for raw_key, value in values.items():
        key = _nonempty(raw_key, what=f"{what} key")
        if key in result:
            raise ValueError(f"{what} keys must not collide after normalization")
        if not isinstance(value, (str, int, float, bool, type(None))):
            raise TypeError(f"{what}[{key!r}] must be a JSON scalar")
        if isinstance(value, float) and not math.isfinite(value):
            raise ValueError(f"{what}[{key!r}] must be finite")
        result[key] = value
    return MappingProxyType(result)


@dataclass(frozen=True, slots=True)
class CellInstanceManifest:
    """One stable cell identity and transform within a population."""

    cell_id: str
    archetype_manifest_hash: str
    translation_um: tuple[float, float, float]
    state_overrides: Mapping[str, JSONScalar] = field(default_factory=dict)

    def __post_init__(self) -> None:
        object.__setattr__(self, "cell_id", _nonempty(self.cell_id, what="cell_id"))
        if not _SHA256.fullmatch(self.archetype_manifest_hash):
            raise ValueError("archetype_manifest_hash must be a lowercase SHA-256 digest")
        translation = tuple(self.translation_um)
        if len(translation) != 3 or any(
            isinstance(value, bool)
            or not isinstance(value, (int, float))
            or not math.isfinite(value)
            for value in translation
        ):
            raise ValueError("translation_um must be a finite (x, y, z) tuple")
        object.__setattr__(
            self,
            "translation_um",
            tuple(float(value) for value in translation),
        )
        object.__setattr__(
            self,
            "state_overrides",
            _scalar_mapping(self.state_overrides, what="state_overrides"),
        )


@dataclass(frozen=True, slots=True)
class CellInteraction:
    """One namespaced cell-cell interaction whose endpoint cells remain separate state owners."""

    interaction_id: str
    interaction_kind: str
    endpoint_cell_ids: tuple[str, ...]
    law_id: str

    def __post_init__(self) -> None:
        for field_name in ("interaction_id", "interaction_kind", "law_id"):
            object.__setattr__(
                self,
                field_name,
                _nonempty(getattr(self, field_name), what=field_name),
            )
        endpoints = tuple(
            _nonempty(endpoint, what="interaction endpoint") for endpoint in self.endpoint_cell_ids
        )
        if len(endpoints) < 2 or len(set(endpoints)) != len(endpoints):
            raise ValueError("cell-cell interaction needs at least two distinct endpoint cells")
        object.__setattr__(self, "endpoint_cell_ids", endpoints)


@dataclass(frozen=True, slots=True)
class CellPopulationManifest:
    """Heterotypic cell instances, their interaction graph, and shared environment."""

    population_id: str
    archetypes: tuple[CellStateManifest, ...]
    cells: tuple[CellInstanceManifest, ...]
    interactions: tuple[CellInteraction, ...]
    shared_environment: Mapping[str, JSONScalar]
    shared_components: tuple[str, ...] = ()
    missing_collective_laws: tuple[str, ...] = ()
    provenance: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "population_id",
            _nonempty(self.population_id, what="population_id"),
        )
        archetypes = tuple(self.archetypes)
        if not archetypes or any(
            not isinstance(archetype, CellStateManifest) for archetype in archetypes
        ):
            raise ValueError("archetypes must contain at least one CellStateManifest")
        archetype_hashes = tuple(archetype.manifest_hash for archetype in archetypes)
        if len(set(archetype_hashes)) != len(archetype_hashes):
            raise ValueError("archetype manifests must be unique")
        object.__setattr__(self, "archetypes", archetypes)

        cells = tuple(self.cells)
        if not cells or any(not isinstance(cell, CellInstanceManifest) for cell in cells):
            raise ValueError("cells must contain at least one CellInstanceManifest")
        cell_ids = tuple(cell.cell_id for cell in cells)
        if len(set(cell_ids)) != len(cell_ids):
            raise ValueError("cell IDs must be unique")
        unknown_archetypes = {
            cell.archetype_manifest_hash
            for cell in cells
            if cell.archetype_manifest_hash not in archetype_hashes
        }
        if unknown_archetypes:
            raise ValueError("cell instance references an unknown archetype manifest")
        object.__setattr__(self, "cells", cells)

        interactions = tuple(self.interactions)
        interaction_ids = tuple(interaction.interaction_id for interaction in interactions)
        if any(not isinstance(interaction, CellInteraction) for interaction in interactions):
            raise TypeError("interactions must contain CellInteraction values")
        if len(set(interaction_ids)) != len(interaction_ids):
            raise ValueError("interaction IDs must be unique")
        known_cells = set(cell_ids)
        for interaction in interactions:
            unknown_cells = set(interaction.endpoint_cell_ids) - known_cells
            if unknown_cells:
                raise ValueError(
                    "interaction references unknown cell IDs: " + ", ".join(sorted(unknown_cells))
                )
        object.__setattr__(self, "interactions", interactions)
        object.__setattr__(
            self,
            "shared_environment",
            _scalar_mapping(self.shared_environment, what="shared_environment"),
        )
        for field_name in (
            "shared_components",
            "missing_collective_laws",
            "provenance",
        ):
            values = tuple(
                _nonempty(value, what=f"{field_name} entry") for value in getattr(self, field_name)
            )
            if len(set(values)) != len(values):
                raise ValueError(f"{field_name} must not contain duplicates")
            object.__setattr__(self, field_name, values)

    def to_dict(self) -> dict[str, object]:
        return {
            "population_id": self.population_id,
            "archetypes": [
                archetype.to_dict()
                for archetype in sorted(
                    self.archetypes,
                    key=lambda item: item.manifest_hash,
                )
            ],
            "cells": [
                {
                    "cell_id": cell.cell_id,
                    "archetype_manifest_hash": cell.archetype_manifest_hash,
                    "translation_um": list(cell.translation_um),
                    "state_overrides": dict(sorted(cell.state_overrides.items())),
                }
                for cell in sorted(self.cells, key=lambda item: item.cell_id)
            ],
            "interactions": [
                {
                    "interaction_id": interaction.interaction_id,
                    "interaction_kind": interaction.interaction_kind,
                    "endpoint_cell_ids": list(interaction.endpoint_cell_ids),
                    "law_id": interaction.law_id,
                }
                for interaction in sorted(
                    self.interactions,
                    key=lambda item: item.interaction_id,
                )
            ],
            "shared_environment": dict(sorted(self.shared_environment.items())),
            "shared_components": sorted(self.shared_components),
            "missing_collective_laws": sorted(self.missing_collective_laws),
            "provenance": sorted(self.provenance),
        }

    @property
    def population_hash(self) -> str:
        payload = json.dumps(
            self.to_dict(),
            sort_keys=True,
            separators=(",", ":"),
            ensure_ascii=False,
        ).encode()
        return hashlib.sha256(payload).hexdigest()
