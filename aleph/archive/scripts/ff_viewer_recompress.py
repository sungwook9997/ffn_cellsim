#!/usr/bin/env python3
"""Re-encode an existing "FF viewer" HTML through the current packer — same physics, smaller file.

The viewers written before buffer pooling (commit 6a5f67a4) and int16 position transport
(45bdfdd6) store every layer's geometry as its own uncompressed float32 base64 string, so a
buffer shared by six scenes is embedded six times. Regenerating them normally would mean
re-running the original simulation, and for most of these the exact script+arguments are no
longer recoverable from the filename.

That is unnecessary: the HTML already contains the full geometry. This decodes the embedded
payload back into arrays and feeds it through :func:`ff_viewer_html.build_viewer`, which pools
duplicate buffers and quantizes positions. Nothing is resampled — every filament, node and
segment survives; only the transport encoding changes (positions gain a bounded ~1e-4 um
quantization error, ~15x below one filament radius).

Safety: the rebuilt file is written to a temporary path and only replaces the original after it
is verified non-empty, so a failed conversion never destroys an artifact.

Usage:
    python -m aleph.scripts.ff_viewer_recompress outputs/ff/figs/native_cortex_morph.html
    python -m aleph.scripts.ff_viewer_recompress --dry-run outputs/ff/figs/*.html
"""
from __future__ import annotations

import argparse
import base64
import json
import re
import shutil
import tempfile
from pathlib import Path

import numpy as np

from aleph.scripts.ff_viewer_html import build_viewer

_ANCHOR = re.compile(r"const P\s*=\s*")


def extract_payload(html: str) -> dict:
    """Extract the embedded ``const P = {...}`` viewer payload from a viewer HTML.

    Args:
        html: Full HTML source.

    Returns:
        The decoded payload dict.

    Raises:
        ValueError: If no payload object can be located or parsed.
    """
    m = _ANCHOR.search(html)
    if not m:
        raise ValueError("no `const P =` payload anchor")
    start = html.index("{", m.end())
    decoder = json.JSONDecoder()
    payload, _ = decoder.raw_decode(html[start:])
    return payload


def _blob(payload: dict, ref, dtype) -> np.ndarray:
    """Decode one geometry buffer, accepting both payload generations.

    Older viewers store the base64 string inline on the layer; newer ones store an integer
    index into a shared ``blobs`` pool (and may carry int16 + ``quant`` dequantization meta).
    """
    if isinstance(ref, str):
        raw = base64.b64decode(ref)
        return np.frombuffer(raw, dtype=dtype)
    blobs = payload.get("blobs") or []
    raw = base64.b64decode(blobs[int(ref)])
    quant = (payload.get("quant") or {}).get(str(int(ref)))
    if quant is not None and dtype == np.float32:
        q = np.frombuffer(raw, dtype=np.int16).reshape(-1, 3).astype(np.float32)
        origin = np.asarray(quant[:3], dtype=np.float32)
        span = np.asarray(quant[3:], dtype=np.float32)
        return (origin + (q + 32767.0) / 65534.0 * span).reshape(-1)
    return np.frombuffer(raw, dtype=dtype)


def payload_to_scenes(payload: dict) -> dict:
    """Rebuild the ``scenes`` dict that :func:`build_viewer` consumes from a decoded payload."""
    scenes: dict[str, list[dict]] = {}
    for sname, layers in (payload.get("scenes") or {}).items():
        out: list[dict] = []
        for layer in layers:
            kind = layer.get("kind", "lines")
            if kind == "plates":
                out.append({"name": layer.get("name", "plates"), "kind": "plates",
                            "verts": layer.get("z", []), "half_xy": layer.get("half_xy", 12.0),
                            "color": layer.get("color", "#888888")})
                continue
            entry = {k: layer[k] for k in
                     ("name", "kind", "color", "size", "opacity", "on_top", "clip",
                      "group", "swatch", "pt_size", "cbar") if k in layer}
            entry["kind"] = kind
            entry["verts"] = _blob(payload, layer["b64"], np.float32).reshape(-1, 3)
            if "faces_b64" in layer:
                entry["faces"] = _blob(payload, layer["faces_b64"], np.uint32).reshape(-1, 3)
            if layer.get("frames"):
                entry["frames"] = [_blob(payload, f, np.float32).reshape(-1, 3)
                                   for f in layer["frames"]]
            if layer.get("color_frames"):
                entry["color_frames"] = [_blob(payload, c, np.uint8).reshape(-1, 3)
                                         for c in layer["color_frames"]]
            out.append(entry)
        scenes[sname] = out
    return scenes


def recompress(path: Path, *, dry_run: bool = False, force: bool = False) -> tuple[int, int]:
    """Re-encode one viewer in place (atomically). Returns ``(bytes_before, bytes_after)``.

    ``bytes_after == bytes_before`` means the file was left untouched.

    Args:
        path: Viewer HTML to rewrite.
        dry_run: Report the result without writing.
        force: Replace even when the rebuild is not smaller. Needed to roll an already-packed
            viewer onto an improved template (e.g. fat lines, scale bar), where the payload is
            already minimal and the win is in the rendering code rather than the byte count.
    """
    before = path.stat().st_size
    html = path.read_text(errors="replace")
    payload = extract_payload(html)
    scenes = payload_to_scenes(payload)
    if not scenes:
        return before, before
    with tempfile.NamedTemporaryFile(suffix=".html", delete=False) as tmp:
        tmp_path = Path(tmp.name)
    build_viewer(scenes, out=str(tmp_path), title=payload.get("title", "FF viewer"),
                 cbars=payload.get("cbars") or {})
    after = tmp_path.stat().st_size
    if after <= 0:
        tmp_path.unlink(missing_ok=True)
        raise ValueError("rebuilt viewer is empty")
    if dry_run or (after >= before and not force):
        tmp_path.unlink(missing_ok=True)
        return before, before
    shutil.move(str(tmp_path), str(path))
    return before, after


def main() -> None:
    """CLI entry point."""
    ap = argparse.ArgumentParser(description=__doc__.split("\n", 1)[0])
    ap.add_argument("paths", nargs="+", type=Path)
    ap.add_argument("--dry-run", action="store_true", help="report the gain, write nothing")
    ap.add_argument("--force", action="store_true",
                    help="rewrite even if not smaller — use to roll viewers onto an improved template")
    args = ap.parse_args()

    total_before = total_after = 0
    for path in args.paths:
        if not path.is_file():
            continue
        try:
            before, after = recompress(path, dry_run=args.dry_run, force=args.force)
        except Exception as exc:                     # keep going; the original is untouched
            print(f"  SKIP {path.name:44s} {type(exc).__name__}: {exc}")
            continue
        total_before += before
        total_after += after
        tag = "same" if after == before else f"{before / max(after, 1):.1f}x"
        print(f"  {'DRY ' if args.dry_run else 'OK  '} {path.name:44s} "
              f"{before / 1e6:8.1f} MB -> {after / 1e6:7.1f} MB  {tag}")
    if total_before:
        print(f"\n  total {total_before / 1e6:.0f} MB -> {total_after / 1e6:.0f} MB "
              f"({(1 - total_after / total_before) * 100:.0f}% smaller)")


if __name__ == "__main__":
    main()
