"""Hash the untracked raw acquisition data, so it can be verified — or re-acquired after deletion.

About 12.7 GB of downloaded microscopy and matrix files were never in git, in this tree or in
`Project_Aleph`'s. `aleph/outputs/outer/PROVENANCE.md` records that retiring that tree must not
delete them "until the raw data has somewhere else to live". Whichever way the PI decides that, both
answers need the same thing first: a per-file record of what is on disk, and a check that it matches
what the acquisition records claim.

What this writes is that record — `raw_acquisition_manifest.jsonl`, one line per file, small enough
to commit. It is not a copy and it does not make one. It makes deletion *checkable* rather than
safe: a manifest says what would have to come back, not that it can.

Each file is checked against the digest its own acquisition record already pinned, where one exists.
Two things came out of building this and both are recorded rather than smoothed over:

- **The AllenCell filename prefix is not a digest of the file.** It looks exactly like one — 64 hex
  characters — and a first version of this module checked against it and reported all 64 files
  corrupt. It is AllenCell's S3 logical key. A check that fails on 100% of its inputs is evidence
  about the check, and the rule this tree already writes for a threshold applies to a verifier too:
  it may not be re-aimed after seeing the data without saying so, which is what this paragraph is.
- **The 25 files of the `downloads` block have no recorded digest anywhere.** They are named by GEO
  and SCP accession, so they can be fetched again, but nothing pins the bytes that were used. They
  are recorded here as `verified_against: null` — unverifiable is a finding, not a blank.

Sanity Gate
-----------
- Dimensions: bytes in, hex digest out; sizes are `st_size`, no unit conversion anywhere.
- Boundary: an absent root is reported and skipped, never silently counted as zero files.
- Invariant: `--check` re-reads the manifest and recomputes; equal digests or a non-zero exit.
- Sign sense: `verified` True/False/None are three states, not two — matched, contradicted, and
  never pinned. Collapsing the third into either of the others is how 1.9 GB with no digest would
  come to read as checked.
- Measurement protocol: digests over whole files, 1 MiB blocks, no sampling.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

DEFAULT_ROOT = Path.home() / (
    "Project_Aleph/.claude/worktrees/external-training-corpus"
    "/corpus/external_training/experiment_factory/outer_library"
)
MANIFEST = Path(__file__).with_name("raw_acquisition_manifest.jsonl")

#: Each block is (name, path under the root, how it is re-acquired). The third field is a pointer to
#: a record that already exists — this module does not restate provenance, it locates it.
BLOCKS: list[tuple[str, str, str]] = [
    (
        "allencell_label_free_confirmation",
        "data/allencell_label_free_confirmation",
        "results/allencell_label_free_acquisition.json (per-file url + sha256 + etag)",
    ),
    (
        "downloads",
        "data/corpus/external_training/experiment_factory/outer_library/downloads",
        "download_*.py per GEO/SCP accession — no digest was recorded for these",
    ),
]

#: Acquisition records that pin a digest, as (record path, how to read it into {name: sha256}).
#: `allencell_label_free_acquisition.json` keys by path within the block; the receipts key by the
#: downloaded filename. Nothing else in the tree pins a raw payload.
PINNED = "results/allencell_label_free_acquisition.json"
RECEIPTS = "acquisition_receipts.jsonl"


def sha256(path: Path) -> str:
    """Digest a file in 1 MiB blocks, so a 9 GB TIFF does not become a 9 GB allocation."""
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def pinned_digests(root: Path) -> dict[str, str]:
    """Collect every digest the acquisition records already pin, keyed by basename.

    Keyed by basename rather than by path because the two records disagree about what a path is —
    one stores `crop_raw/<file>` relative to its block, the other stores a bare filename — and the
    basenames are unique across all 89 files, which is asserted rather than assumed.
    """
    pinned: dict[str, str] = {}
    record = root / PINNED
    if record.is_file():
        for entry in json.loads(record.read_text()).get("files", []):
            pinned[Path(entry["path"]).name] = entry["sha256"]
    receipts = root / RECEIPTS
    if receipts.is_file():
        for line in receipts.read_text().splitlines():
            entry = json.loads(line)
            if entry.get("sha256"):
                pinned[Path(entry["downloaded_file"]).name] = entry["sha256"]
    return pinned


def scan(root: Path) -> list[dict]:
    """Walk every block under ``root`` and return one record per file, sorted by path."""
    pinned = pinned_digests(root)
    records: list[dict] = []
    for name, rel, reacquire in BLOCKS:
        base = root / rel
        if not base.is_dir():
            print(f"  ! block {name!r} absent at {base}", file=sys.stderr)
            continue
        files = sorted(p for p in base.rglob("*") if p.is_file())
        print(f"  {name}: {len(files)} files", file=sys.stderr)
        for path in files:
            digest = sha256(path)
            expected = pinned.get(path.name)
            records.append(
                {
                    "block": name,
                    "path": str(path.relative_to(root)),
                    "size_bytes": path.stat().st_size,
                    "sha256": digest,
                    # None = no record pinned this file. False = a record pinned it and disagrees.
                    "verified": None if expected is None else expected == digest,
                    "verified_against": None if expected is None else PINNED,
                    "reacquire_via": reacquire,
                }
            )
    names = [Path(r["path"]).name for r in records]
    if len(set(names)) != len(names):
        raise ValueError("basenames are not unique — keying digests by basename would mis-verify")
    return records


def report(records: list[dict]) -> int:
    """Print the totals and return the number of files that contradict a pinned digest."""
    total = sum(r["size_bytes"] for r in records)
    bad = [r for r in records if r["verified"] is False]
    unpinned = [r for r in records if r["verified"] is None]
    print(f"{len(records)} files, {total / 1e9:.2f} GB", file=sys.stderr)
    print(
        f"  verified {sum(1 for r in records if r['verified'])}"
        f" · contradicted {len(bad)}"
        f" · never pinned {len(unpinned)}"
        f" ({sum(r['size_bytes'] for r in unpinned) / 1e9:.2f} GB)",
        file=sys.stderr,
    )
    for r in bad:
        print(f"  MISMATCH {r['path']}", file=sys.stderr)
    return len(bad)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--root", type=Path, default=DEFAULT_ROOT)
    ap.add_argument(
        "--check",
        action="store_true",
        help="recompute against the committed manifest instead of writing it; "
        "run this before anyone deletes the source tree",
    )
    args = ap.parse_args()

    records = scan(args.root)
    if not records:
        print("no files found — refusing to write an empty manifest", file=sys.stderr)
        return 2
    mismatched = report(records)

    if args.check:
        if not MANIFEST.exists():
            print(f"no manifest at {MANIFEST}", file=sys.stderr)
            return 2
        stored = {json.loads(line)["path"]: json.loads(line) for line in MANIFEST.read_text().splitlines()}
        drift = [r["path"] for r in records if stored.get(r["path"], {}).get("sha256") != r["sha256"]]
        missing = sorted(set(stored) - {r["path"] for r in records})
        for p in drift:
            print(f"  DRIFT   {p}", file=sys.stderr)
        for p in missing:
            print(f"  MISSING {p}", file=sys.stderr)
        print(f"check: {len(drift)} drifted, {len(missing)} missing", file=sys.stderr)
        return 1 if (drift or missing) else 0

    MANIFEST.write_text("".join(json.dumps(r, sort_keys=True) + "\n" for r in records))
    print(f"wrote {MANIFEST}", file=sys.stderr)
    # A mismatch is recorded, not hidden — but it is not a reason to withhold the manifest.
    return 1 if mismatched else 0


if __name__ == "__main__":
    raise SystemExit(main())
