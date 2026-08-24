"""CLI for the combined first/second-pass projection."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from aleph.harness import DataRoot
from aleph.outer.ragtagcag.second_pass.combined import run_combined_validation


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-root", required=True)
    parser.add_argument("--initial-manifest", action="append", required=True)
    parser.add_argument("--second-pass-manifest", required=True)
    parser.add_argument("--receipt", action="append", required=True)
    parser.add_argument("--summary")
    args = parser.parse_args()
    result = run_combined_validation(
        initial_manifests=args.initial_manifest,
        second_pass_manifest=args.second_pass_manifest,
        receipt_paths=args.receipt,
        data_root=DataRoot.resolve(args.data_root),
    )
    rendered = json.dumps(result, indent=2, sort_keys=True, ensure_ascii=False) + "\n"
    if args.summary:
        Path(args.summary).write_text(rendered, encoding="utf-8")
    print(rendered, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
