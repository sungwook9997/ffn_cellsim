"""Command-line entry point for the external RAG--TAG--CAG adapter."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from aleph.harness import DataRoot

from aleph.outer.ragtagcag.ingest import verify_idempotent_ingest


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-root", required=True, help="isolated Aleph data root outside Git")
    parser.add_argument("--manifest", action="append", required=True)
    parser.add_argument("--receipt", action="append", default=[])
    parser.add_argument("--summary", help="write a compact JSON validation summary")
    args = parser.parse_args()

    root = DataRoot.resolve(args.data_root)
    result = verify_idempotent_ingest(args.manifest, root, receipt_paths=args.receipt)
    rendered = json.dumps(result, indent=2, sort_keys=True, ensure_ascii=False) + "\n"
    if args.summary:
        Path(args.summary).write_text(rendered, encoding="utf-8")
    print(rendered, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
