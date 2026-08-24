#!/usr/bin/env python3
"""Build experiment-hierarchical features from the official MULTIMOT R object.

The source contains 472,173 cell-time rows, but those rows are not independent
training samples.  This builder reproduces the paper's declared 18-variable
complete-case panel and emits one robust median per technical replicate.  Lab,
person, independent experiment, condition, and technical replicate remain
explicit, so downstream splitting cannot mistake cells or frames for biological
replication.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from schema import Observation


ACCESSION = "10.17044/scilifelab.21407402"
DATASET_ID = "scilifelab-21407402-multisite-live-cell-2d"
SOURCE_SHA256 = "2c8c6a7d263dcff2e2749afbed93f6f7b916d5cbb87429081a83136f4ddd76d6"
PAPER_SELECTED_OBSERVABLES = (
    "Cells_Area",
    "Cells_Compactness",
    "Cells_Eccentricity",
    "Cells_MajorAxisLength",
    "Cells_MaxFeretDiameter",
    "Cells_MaximumRadius",
    "Cells_MeanRadius",
    "Cells_MinFeretDiameter",
    "Cells_MinorAxisLength",
    "Cells_Perimeter",
    "Cells_Solidity",
    "Cells_Protrusions",
    "Cells_Retractions",
    "Cells_ShortLivedRegions",
    "Cells_ICS",
    "Nuc_Area",
    "Nuc_INS",
    "Nuc_Perimeter",
)
GEOMETRY_AREA = {"Cells_Area", "Nuc_Area"}
GEOMETRY_LENGTH = {
    "Cells_MajorAxisLength", "Cells_MaxFeretDiameter", "Cells_MaximumRadius",
    "Cells_MeanRadius", "Cells_MinFeretDiameter", "Cells_MinorAxisLength",
    "Cells_Perimeter", "Nuc_Perimeter",
}
DIMENSIONLESS = {"Cells_Compactness", "Cells_Eccentricity", "Cells_Solidity"}
AUTHOR_SPEED_SCALE = 0.8260495552435883 / 6.0


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _unit(observable: str) -> str:
    if observable in GEOMETRY_AREA:
        return "square_micrometre"
    if observable in GEOMETRY_LENGTH:
        return "micrometre"
    if observable in DIMENSIONLESS:
        return "dimensionless"
    if observable in {"Cells_ICS", "Nuc_INS"}:
        # The released analysis applies this scale but does not attach a unit.
        # Retaining that uncertainty prevents an unsafe physical-axis mapping.
        return "author_scaled_source_unit"
    return "author_processed_source_unit"


def _condition(value: str) -> tuple[str, str]:
    mapping = {
        "C": ("control", "untreated_migrating_HT1080"),
        "T": ("ROCK_inhibited", "ROCK_inhibited_migrating_HT1080"),
    }
    if value not in mapping:
        raise ValueError(f"unexpected author condition: {value}")
    return mapping[value]


def build(source: Path) -> list[Observation]:
    digest = _sha256(source)
    if digest != SOURCE_SHA256:
        raise ValueError(f"official MULTIMOT dat_v2.R hash mismatch: {digest}")
    try:
        import pyreadr
    except ImportError as exc:  # pragma: no cover - dependency error is explicit
        raise RuntimeError("build_multisite_migration_manifest requires pyreadr") from exc
    objects = pyreadr.read_r(str(source))
    if set(objects) != {"dat"}:
        raise ValueError(f"expected only R object 'dat', found {list(objects)}")
    frame = objects["dat"]
    required = {
        *PAPER_SELECTED_OBSERVABLES, "Lab", "Person", "Experiment",
        "Technical_replicate", "Condition",
    }
    missing = sorted(required - set(frame.columns))
    if missing:
        raise ValueError(f"MULTIMOT source is missing columns: {missing}")
    complete = frame.loc[frame[list(PAPER_SELECTED_OBSERVABLES)].notna().all(axis=1)].copy()
    complete.loc[:, ["Cells_ICS", "Nuc_INS"]] *= AUTHOR_SPEED_SCALE
    group_columns = [
        "Lab", "Person", "Experiment", "Technical_replicate", "Condition"
    ]
    rows: list[Observation] = []
    for keys, group in complete.groupby(group_columns, observed=True, sort=True):
        lab, person, experiment, technical, author_condition = map(str, keys)
        condition, state = _condition(author_condition)
        biological = f"lab_{lab}:person_{person}:experiment_{experiment}"
        technical_number = technical[-1]
        if technical[0] != author_condition or technical_number not in {"1", "2", "3"}:
            raise ValueError(f"condition/technical mismatch: {keys}")
        sample = f"{DATASET_ID}:{biological}:{condition}:technical_{technical_number}"
        medians = group[list(PAPER_SELECTED_OBSERVABLES)].median(axis=0)
        for observable in PAPER_SELECTED_OBSERVABLES:
            transform = (
                f"*{AUTHOR_SPEED_SCALE:.17g}" if observable in {"Cells_ICS", "Nuc_INS"}
                else "identity"
            )
            locator = (
                "R:dat;complete_case=paper_18_variable_panel;"
                f"Lab={lab};Person={person};Experiment={experiment};"
                f"Technical_replicate={technical};Condition={author_condition};"
                f"column={observable};statistic=median;n={len(group)};transform={transform}"
            )
            identifier = hashlib.sha256(f"{digest}:{locator}".encode()).hexdigest()[:20]
            row = Observation(
                observation_id=f"multimot-2d:{identifier}",
                dataset_id=DATASET_ID,
                lab_group=f"lab_{lab}",
                provider="SciLifeLab Data Repository (Figshare)",
                accession=ACCESSION,
                license_id="cc-by-4.0",
                source_sha256=digest,
                source_locator=locator,
                sample_id=sample,
                biological_replicate=biological,
                technical_replicate=f"technical_{technical_number}",
                condition=condition,
                cell_type="HT1080 H2B-EGFP Lifeact-mCherry",
                cell_state=state,
                modality="high_content_live_cell_imaging_derived",
                observable=observable,
                value=float(medians[observable]),
                unit=_unit(observable),
            )
            row.validate()
            rows.append(row)
    if len(rows) != 162 * len(PAPER_SELECTED_OBSERVABLES):
        raise ValueError(f"expected 2,916 observations, found {len(rows)}")
    if len({row.sample_id for row in rows}) != 162:
        raise ValueError("expected 162 technical-replicate samples")
    return rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    rows = build(args.source)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8") as handle:
        for row in rows:
            handle.write(json.dumps(row.as_dict(), sort_keys=True) + "\n")
    print(json.dumps({
        "observations": len(rows),
        "technical_replicates": len({row.sample_id for row in rows}),
        "independent_experiments": len({row.biological_replicate for row in rows}),
        "labs": sorted({row.lab_group for row in rows}),
        "observables": len(PAPER_SELECTED_OBSERVABLES),
    }, sort_keys=True))


if __name__ == "__main__":
    main()
