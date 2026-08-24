from __future__ import annotations

import pytest

from aleph.virtual_cell.contracts import CellStateManifest
from aleph.virtual_cell.population import (
    CellInstanceManifest,
    CellInteraction,
    CellPopulationManifest,
)


def _archetype(
    manifest_id: str,
    lineage: str,
    components: tuple[str, ...],
    *,
    missing_laws: tuple[str, ...] = (),
) -> CellStateManifest:
    return CellStateManifest(
        manifest_id=manifest_id,
        species="Homo sapiens",
        lineage=lineage,
        identity={"source": "analytic fixture"},
        biological_state={"state": "generic"},
        environment={"context": "co-culture"},
        components=components,
        connectors=(),
        missing_laws=missing_laws,
    )


def test_heterotypic_population_namespaces_same_component_names() -> None:
    endothelial = _archetype(
        "endothelial-flow",
        "endothelial",
        ("membrane", "cortex", "nucleus", "junction"),
    )
    immune = _archetype(
        "amoeboid-immune",
        "immune",
        ("membrane", "cortex", "nucleus"),
        missing_laws=("fast-amoeboid-remodelling",),
    )
    population = CellPopulationManifest(
        population_id="endothelial-immune-contact",
        archetypes=(endothelial, immune),
        cells=(
            CellInstanceManifest(
                "endo-001",
                endothelial.manifest_hash,
                (0.0, 0.0, 0.0),
            ),
            CellInstanceManifest(
                "immune-001",
                immune.manifest_hash,
                (12.0, 0.0, 0.0),
                state_overrides={"activation": "rolling"},
            ),
        ),
        interactions=(
            CellInteraction(
                "contact-001",
                "heterotypic-adhesion",
                ("endo-001", "immune-001"),
                "selectin-integrin-law-missing",
            ),
        ),
        shared_environment={"flow": "present", "matrix": "basement membrane"},
        shared_components=("extracellular-medium", "ecm"),
        missing_collective_laws=("selectin-integrin-law-missing",),
    )
    assert len(population.population_hash) == 64
    assert population.cells[0].cell_id != population.cells[1].cell_id
    assert "membrane" in endothelial.components
    assert "membrane" in immune.components


def test_anucleate_archetype_does_not_force_nucleus_component() -> None:
    erythrocyte = _archetype(
        "erythrocyte",
        "erythroid",
        ("membrane", "cortex"),
        missing_laws=("spectrin-bilayer-law",),
    )
    assert "nucleus" not in erythrocyte.components
    assert "spectrin-bilayer-law" in erythrocyte.missing_laws


def test_population_rejects_interaction_with_unknown_cell() -> None:
    archetype = _archetype("generic", "generic", ("membrane",))
    with pytest.raises(ValueError, match="unknown cell IDs"):
        CellPopulationManifest(
            population_id="bad-contact",
            archetypes=(archetype,),
            cells=(
                CellInstanceManifest(
                    "cell-001",
                    archetype.manifest_hash,
                    (0.0, 0.0, 0.0),
                ),
            ),
            interactions=(
                CellInteraction(
                    "contact",
                    "junction",
                    ("cell-001", "cell-missing"),
                    "junction-law",
                ),
            ),
            shared_environment={},
        )
