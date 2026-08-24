#!/usr/bin/env python3
"""Download a pinned, paper-linked M-TRACK A549 EMT trajectory pilot.

The Science Advances paper links the public M-TRACK repository, whose README
links the GUI input folder on Google Drive.  The released GUI example is one
TGF-beta-treated ``xy01`` trajectory.  It is deliberately sampled sparsely and
must never be counted as ten independent biological samples.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import hashlib
import json
from pathlib import Path
import time
import urllib.error
import urllib.request


PAPER_DOI = "10.1126/sciadv.aba9319"
PAPER_URL = "https://pmc.ncbi.nlm.nih.gov/articles/PMC7473671/"
REPOSITORY_URL = "https://github.com/opnumten/M-TRACK"
PUBLIC_FOLDER_ID = "1KwVwR9moblQlA0U2aDPWHNeKDNONrk9J"

PINNED_PAIRS = (
    (3, "12HIwjd7H4D69m-J0sE3mKVlosvnscJEL", 10062,
     "1cTuJPB-rR2MKAuo6eg1LCyIkudDin3y0", 7645628),
    (13, "14a-kfz6O3sfPgS2c6AAcRlmvB3rPdGAH", 10057,
     "1opUR6kiFM3j7cidvH-jWFGfmP-0VuXjT", 7645628),
    (25, "1CFjeD0JAKq3JbL6bNsqBiuT2F_56VC4I", 10459,
     "1jJo-qXpl3YI8TpicGiLGrYnemOD3jXQH", 7645628),
    (37, "140eSF3imy7AUbaX9s-r88QKzAs8Om7a9", 11386,
     "1ARKU63bAEsVVnY4Vmy5NkSG1ybNRpjSY", 7645628),
    (61, "1a8VZTvNJFptJ3ZBPGpqU4jqL5AkY9iRN", 11997,
     "1j68vSfZ0T3--Z5UUAqd3ZO45qBLyV10K", 7645628),
    (77, "1yQ_D8C_1mUaNE_l9RHqdATs2IhO8ZW_M", 12817,
     "1klXq6MVAXkOIrvFO_Zct1VYPxKUye11z", 7645628),
    (93, "1sCysARaFDCEdwoyn3GxKa_A3ziFjwdrP", 12328,
     "1POa5UYz1RiT6QCHBZNkzDyVb6lYsgqAP", 7645628),
    (107, "1BXk5SWxraOnE2Zwypza_PINwVCtvY3FW", 12548,
     "1MzAoTHnEuH0cvt8TZglVwdrlDuu-oTNR", 7645628),
    (123, "1JC5p4C7OQddBxZPd9_4oG8Ateclm2nb8", 12953,
     "1llwaI0tcZMzhoFTdh0iLYs6yXrU7AkTi", 7645628),
    (145, "1tfLrkAs9K9su0e4baj0-HY0dzXl91qz_", 14483,
     "1Mrr0ZzxPbw9v22EidXZqygtSwSzJLDJY", 7645628),
)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _download(file_id: str, destination: Path, expected_size: int) -> dict[str, object]:
    url = (
        "https://drive.usercontent.google.com/download"
        f"?id={file_id}&export=download&confirm=t"
    )
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.is_file() and destination.stat().st_size == expected_size:
        return {
            "google_drive_file_id": file_id,
            "local_name": destination.name,
            "size": expected_size,
            "sha256": _sha256(destination),
        }
    temporary = destination.with_suffix(destination.suffix + ".part")
    headers = {"User-Agent": "Project-Aleph-research/1.0"}
    for attempt in range(6):
        try:
            request = urllib.request.Request(url, headers=headers)
            with urllib.request.urlopen(request, timeout=180) as response, temporary.open("wb") as handle:
                while True:
                    block = response.read(1024 * 1024)
                    if not block:
                        break
                    handle.write(block)
            if temporary.stat().st_size != expected_size:
                raise ValueError(
                    f"Google Drive size changed for {file_id}: "
                    f"{temporary.stat().st_size} != {expected_size}"
                )
            temporary.replace(destination)
            return {
                "google_drive_file_id": file_id,
                "local_name": destination.name,
                "size": expected_size,
                "sha256": _sha256(destination),
            }
        except (urllib.error.URLError, TimeoutError):
            if attempt == 5:
                raise
            time.sleep(min(30, 2 ** attempt))
    raise RuntimeError("unreachable")


def acquire(output_dir: Path, workers: int = 6) -> dict[str, object]:
    jobs = []
    for time_index, mask_id, mask_size, fluor_id, fluor_size in PINNED_PAIRS:
        jobs.extend((
            (time_index, "segmentation_mask", mask_id, mask_size,
             output_dir / f"seg_a549_tgf4ngxy01t{time_index:03d}c2.png"),
            (time_index, "VIM_RFP", fluor_id, fluor_size,
             output_dir / f"a549_tgf4ngxy01t{time_index:03d}c1.tif"),
        ))
    records: list[dict[str, object]] = []
    with ThreadPoolExecutor(max_workers=max(1, workers)) as executor:
        futures = {
            executor.submit(_download, file_id, path, size): (time_index, kind)
            for time_index, kind, file_id, size, path in jobs
        }
        for future in as_completed(futures):
            time_index, kind = futures[future]
            record = future.result()
            record.update({"time_index": time_index, "kind": kind})
            records.append(record)
    records.sort(key=lambda row: (int(row["time_index"]), str(row["kind"])))
    return {
        "schema": "aleph.outer_library.mtrack_emt_acquisition.v1",
        "source": {
            "paper_doi": PAPER_DOI,
            "paper_url": PAPER_URL,
            "repository_url": REPOSITORY_URL,
            "public_google_drive_folder_id": PUBLIC_FOLDER_ID,
            "release_role": "paper_linked_GUI_input_example",
            "license_id": "not_stated_for_image_data",
            "connector_visibility_status": (
                "access_not_verified_but_paper_repository_links_public_folder"
            ),
        },
        "selection": {
            "trajectory": "A549_TGFbeta_4ng_per_mL_xy01",
            "time_indices": [row[0] for row in PINNED_PAIRS],
            "selection_frozen_before_feature_evaluation": True,
            "all_frames_form_one_split_group": True,
        },
        "files": records,
        "counts": {"timepoints": len(PINNED_PAIRS), "files": len(records)},
        "limitations": [
            "the GUI release exposes one treated xy01 trajectory, not a biological replicate panel",
            "no untreated control trajectory is present in the linked GUI example",
            "the image-data license is not stated in the repository or public folder metadata",
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--receipt", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=6)
    args = parser.parse_args()
    receipt = acquire(args.output_dir, args.workers)
    args.receipt.parent.mkdir(parents=True, exist_ok=True)
    args.receipt.write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(receipt["counts"], sort_keys=True))


if __name__ == "__main__":
    main()
