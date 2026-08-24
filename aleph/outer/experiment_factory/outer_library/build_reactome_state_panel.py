#!/usr/bin/env python3
"""Build a provenance-locked state-gene panel from Reactome participants."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
from pathlib import Path
import re
import urllib.request


PATHWAYS = {
    "cell_cycle": "R-HSA-1640170",
    "apoptosis": "R-HSA-109581",
    "cellular_senescence": "R-HSA-2559583",
    "TGF_beta_EMT": "R-HSA-2173791",
    "hypoxia_response": "R-HSA-1234174",
    "DNA_repair": "R-HSA-73894",
    "unfolded_protein_response": "R-HSA-381119",
    "innate_immune_system": "R-HSA-168249",
    "interferon_signaling": "R-HSA-913531",
}
BASE = "https://reactome.org/ContentService/data/participants"


def _gene_annotations(path: Path) -> tuple[dict[str, list[int]], int]:
    by_symbol: dict[str, list[int]] = {}
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        header = next(handle).split()
        if header != ["id", "gene_short_name"]:
            raise ValueError("sci-Plex gene annotation schema changed")
        count = 0
        for index, line in enumerate(handle, start=1):
            _, symbol = line.rstrip("\n").split(" ", 1)
            by_symbol.setdefault(symbol, []).append(index)
            count += 1
    return by_symbol, count


def build(gene_annotation_path: Path) -> dict[str, object]:
    gene_indices, annotation_count = _gene_annotations(gene_annotation_path)
    pathways = []
    union: set[str] = set()
    for axis, stable_id in PATHWAYS.items():
        url = f"{BASE}/{stable_id}"
        request = urllib.request.Request(url, headers={"User-Agent": "Project-Aleph-research/1.0"})
        with urllib.request.urlopen(request, timeout=60) as response:
            payload = response.read()
        participants = json.loads(payload)
        symbols: set[str] = set()
        for participant in participants:
            for entity in participant.get("refEntities", []):
                if entity.get("schemaClass") not in {"ReferenceGeneProduct", "ReferenceDNASequence"}:
                    continue
                display = entity.get("displayName", "")
                match = re.search(r"(?:UniProt|ENSEMBL):\S+\s+([^\s\[]+)", display)
                if match:
                    symbols.add(match.group(1))
        matched = sorted(symbol for symbol in symbols if symbol in gene_indices)
        union.update(matched)
        pathways.append({
            "axis": axis,
            "reactome_stable_id": stable_id,
            "official_url": f"https://reactome.org/content/detail/{stable_id}",
            "participants_api_url": url,
            "participants_response_sha256": hashlib.sha256(payload).hexdigest(),
            "participant_record_count": len(participants),
            "reference_gene_symbol_count": len(symbols),
            "sciplex_matched_gene_count": len(matched),
            "sciplex_matched_genes": matched,
        })
    panel = sorted(union)
    matrix_rows = {
        symbol: gene_indices[symbol]
        for symbol in panel
    }
    return {
        "schema": "aleph.outer_library.reactome_state_panel.v1",
        "retrieved_at": "2026-08-06",
        "source": "Reactome_ContentService_curated_human_pathway_participants",
        "sciplex_gene_annotation_count": annotation_count,
        "pathway_count": len(pathways),
        "panel_gene_symbol_count": len(panel),
        "panel_gene_symbols": panel,
        "sciplex_matrix_rows_1_based": matrix_rows,
        "pathways": pathways,
        "selection_rule": "all Reactome reference gene products intersecting exact sci-Plex gene symbols",
        "manual_marker_additions": [],
        "training_target_authority": "pathway_feature_panel_not_cell_state_ground_truth",
        "aleph_authority": "none",
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--gene-annotations", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = build(args.gene_annotations)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({
        "pathway_count": result["pathway_count"],
        "panel_gene_symbol_count": result["panel_gene_symbol_count"],
    }, sort_keys=True))


if __name__ == "__main__":
    main()
