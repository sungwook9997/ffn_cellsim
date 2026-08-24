#!/usr/bin/env python3
"""Build a leakage-safe HPA IF embedding tensor from the official v25.1 export.

The split unit is a connected component of the gene--antibody bipartite graph.
Consequently neither a gene nor an antibody can cross train, validation,
calibration, and test.  A fixed, predeclared subset of the provider's 1,024
image features is retained; no labels or holdout rows select dimensions.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import zipfile

import numpy as np


SCHEMA = "aleph.outer_library.hpa_if_embedding_tensor.v1"
DATASET_ID = "hpa-v25.1-subcellular-image-embeddings"
LAB_GROUP = "human-protein-atlas-subcellular-resource"
OFFICIAL_ARCHIVE_SIZE = 632_515_574
OFFICIAL_ARCHIVE_SHA256 = (
    "c23170d2fc1107d66c87c5a181a10b97f19cb4c2b2c5e2f2bc369ba7d547f92d"
)
MEMBER = "subcell_image_umap_features.tsv"
EXPECTED_ROWS = 81_007
FEATURE_INDEXES = np.arange(0, 1024, 8, dtype=np.int64)  # fixed 128-D view


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


class _DisjointSet:
    def __init__(self) -> None:
        self.parent: dict[str, str] = {}
        self.size: dict[str, int] = {}

    def find(self, value: str) -> str:
        if value not in self.parent:
            self.parent[value] = value
            self.size[value] = 1
        root = value
        while self.parent[root] != root:
            root = self.parent[root]
        while self.parent[value] != value:
            nxt = self.parent[value]
            self.parent[value] = root
            value = nxt
        return root

    def union(self, left: str, right: str) -> None:
        left, right = self.find(left), self.find(right)
        if left == right:
            return
        if self.size[left] < self.size[right]:
            left, right = right, left
        self.parent[right] = left
        self.size[left] += self.size[right]


def _tokens(value: str) -> list[str]:
    return [item.strip() for item in value.replace(";", ",").split(",") if item.strip()]


def _split(component: str) -> str:
    bucket = int(hashlib.sha256(component.encode("utf-8")).hexdigest()[:8], 16) % 100
    if bucket < 65:
        return "train"
    if bucket < 80:
        return "validation"
    if bucket < 90:
        return "calibration"
    return "test"


def _open_member(archive: Path):
    container = zipfile.ZipFile(archive)
    if container.namelist() != [MEMBER]:
        container.close()
        raise ValueError("HPA archive member contract changed")
    return container, container.open(MEMBER)


def build(archive: Path, output: Path) -> dict[str, object]:
    if archive.stat().st_size != OFFICIAL_ARCHIVE_SIZE:
        raise ValueError("HPA archive size differs from the official v25.1 receipt")
    source_sha256 = _sha256(archive)
    if source_sha256 != OFFICIAL_ARCHIVE_SHA256:
        raise ValueError("HPA archive SHA-256 differs from the official v25.1 receipt")

    dsu = _DisjointSet()
    antibodies: list[str] = []
    genes_text: list[str] = []
    container, handle = _open_member(archive)
    try:
        header = handle.readline().decode("utf-8").rstrip("\n").split("\t")
        if len(header) != 1034 or header[10] != "image feature 1" or header[-1] != "image feature 1024":
            raise ValueError("HPA embedding header contract changed")
        for raw in handle:
            fields = raw.decode("utf-8").split("\t", 5)
            if len(fields) != 6:
                raise ValueError("malformed HPA metadata row")
            antibody = fields[2]
            genes = _tokens(fields[3])
            if not antibody or not genes:
                raise ValueError("HPA row lacks antibody or gene provenance")
            antibody_node = f"antibody:{antibody}"
            for gene in genes:
                dsu.union(antibody_node, f"gene:{gene}")
            antibodies.append(antibody)
            genes_text.append(",".join(genes))
    finally:
        handle.close()
        container.close()
    if len(antibodies) != EXPECTED_ROWS:
        raise ValueError(f"expected {EXPECTED_ROWS} HPA rows, found {len(antibodies)}")

    component_ids = np.asarray(
        [dsu.find(f"antibody:{antibody}") for antibody in antibodies], dtype="U64"
    )
    split = np.asarray([_split(component) for component in component_ids], dtype="U11")
    X = np.empty((EXPECTED_ROWS, len(FEATURE_INDEXES)), dtype=np.float16)
    file_prefix: list[str] = []
    cell_line: list[str] = []
    locations: list[str] = []
    umap3d = np.empty((EXPECTED_ROWS, 3), dtype=np.float32)

    container, handle = _open_member(archive)
    try:
        handle.readline()
        for row_index, raw in enumerate(handle):
            fields = raw.decode("utf-8").rstrip("\n").split("\t")
            if len(fields) != 1034:
                raise ValueError(f"HPA row {row_index} has {len(fields)} fields")
            file_prefix.append(fields[0])
            cell_line.append(fields[1])
            locations.append(",".join(_tokens(fields[4])))
            umap3d[row_index] = np.asarray(fields[7:10], dtype=np.float32)
            features = np.fromstring("\t".join(fields[10:]), sep="\t", dtype=np.float32)
            if features.shape != (1024,) or not np.isfinite(features).all():
                raise ValueError(f"invalid HPA feature vector at row {row_index}")
            X[row_index] = features[FEATURE_INDEXES]
    finally:
        handle.close()
        container.close()

    # Explicit invariants make accidental sample-level leakage a hard failure.
    entity_splits: dict[str, set[str]] = {}
    for antibody, genes, row_split in zip(antibodies, genes_text, split):
        for entity in [f"antibody:{antibody}", *[f"gene:{g}" for g in _tokens(genes)]]:
            entity_splits.setdefault(entity, set()).add(str(row_split))
    overlap = [entity for entity, splits in entity_splits.items() if len(splits) != 1]
    if overlap:
        raise ValueError(f"gene/antibody leakage across splits: {overlap[:3]}")

    output.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        output,
        X=X,
        feature_indexes=FEATURE_INDEXES,
        file_prefix=np.asarray(file_prefix, dtype="U96"),
        cell_line=np.asarray(cell_line, dtype="U32"),
        antibody=np.asarray(antibodies, dtype="U32"),
        genes=np.asarray(genes_text, dtype="U128"),
        locations=np.asarray(locations, dtype="U256"),
        component_id=component_ids,
        split=split,
        umap3d=umap3d,
        source_sha256=np.asarray([source_sha256], dtype="U64"),
        dataset_id=np.asarray([DATASET_ID], dtype="U64"),
        lab_group=np.asarray([LAB_GROUP], dtype="U64"),
    )
    split_counts = {name: int((split == name).sum()) for name in np.unique(split)}
    component_counts = {
        name: len(set(component_ids[split == name])) for name in np.unique(split)
    }
    return {
        "schema": SCHEMA,
        "dataset_id": DATASET_ID,
        "lab_group": LAB_GROUP,
        "provider_version": "HPA v25.1",
        "license": "CC BY 4.0",
        "official_source": "https://www.proteinatlas.org/about/download",
        "archive_sha256": source_sha256,
        "archive_bytes": archive.stat().st_size,
        "archive_member": MEMBER,
        "image_count": EXPECTED_ROWS,
        "gene_count": len({g for value in genes_text for g in _tokens(value)}),
        "antibody_count": len(set(antibodies)),
        "cell_line_count": len(set(cell_line)),
        "author_location_label_count": len({x for value in locations for x in _tokens(value)}),
        "retained_embedding_dimensions": int(len(FEATURE_INDEXES)),
        "feature_selection": "fixed every eighth provider feature, starting at feature 1",
        "split_unit": "connected_component_of_gene_antibody_bipartite_graph",
        "split_counts": split_counts,
        "component_counts": component_counts,
        "gene_or_antibody_cross_split_overlap_count": 0,
        "independent_lab_holdout": False,
        "aleph_authority": "none",
        "may_select_aleph_parameter": False,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--archive", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--receipt", type=Path, required=True)
    args = parser.parse_args()
    receipt = build(args.archive, args.output)
    args.receipt.parent.mkdir(parents=True, exist_ok=True)
    args.receipt.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(receipt, sort_keys=True))


if __name__ == "__main__":
    main()
