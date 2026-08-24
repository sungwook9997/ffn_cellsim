"""Export the standing cell to a compact binary the native viewer opens — FULL population, no thinning.

**Why this exists: the HTML viewer cannot hold the cell.** A `render_cell_view` page of eight
populations is 78–124 MB of markup and the browser has to parse every line as a DOM-adjacent object
before it draws anything. The cell is 4.6 M nodes and ~4.3 M segments; that is not a large amount of
DATA — it is 86 MB as raw float32 and int32 — it is only large as text. So the fix is not to thin the
cell, which this project forbids, and not to compress the page: it is to stop shipping geometry as a
document.

**Never thinned. That is the whole point.** `feedback-viewer-no-downsample` and
`feedback-viewer-fat-lines-no-dots` are standing instructions, and every viewer this project has
shipped has honoured them; this one keeps them at a population where a page cannot.

**The format is deliberately boring** — a JSON header and raw little-endian arrays, so it can be read
by anything, including a future viewer nobody has written yet:

    magic      b"ALEPHCEL"                      8 bytes
    header_len uint32                            4 bytes, little-endian
    header     UTF-8 JSON                        header_len bytes
    blocks     raw arrays, in header order, contiguous

`positions` is float32 despite the arena being float64. **That is a stated approximation, not an
oversight**: the cell is ~15 µm across and float32 holds ~7 decimal digits, so the worst placement
error is ~1e-6 µm — a picometre, four orders below the finest structure in the cell (the 0.05 µm
cortex mesh). The header records it so a reader is never guessing what it holds. Nothing quantitative
is computed from this file; it exists to be looked at.

Sanity Gate:
    * dimensions — positions µm, radii µm; no unit conversion happens here.
    * boundary cases — a population with zero nodes is skipped and SAID so rather than written empty;
      an empty export is refused, because a viewer opening a file with nothing in it should learn
      that from the exporter, not from a black window.
    * conservation — every index written is checked to land inside the block it addresses, so a
      viewer cannot read outside a population and draw a line to an unrelated structure.
    * precision — float32 positions, error bound stated above and recorded in the header.
    * sign sense — not applicable: this writes geometry and computes no force.
    * measurement protocol — the header records the node and segment count per population, so what
      was drawn is checkable against the run record without opening the binary.

Runtime: reads a CUDA arena. Refuses a non-CUDA device, like every other driver here.
"""

from __future__ import annotations

import argparse
import json
import struct
import time
from pathlib import Path

MAGIC = b"ALEPHCEL"
#: Bumped when the layout changes in a way an existing reader would get wrong.
FORMAT_VERSION = 1


def _write(out: Path, header: dict, blocks: list) -> int:
    """Write the container. Returns bytes written."""
    payload = b"".join(b.tobytes() for b in blocks)
    head = json.dumps(header, separators=(",", ":")).encode("utf-8")
    with out.open("wb") as fh:
        fh.write(MAGIC)
        fh.write(struct.pack("<I", len(head)))
        fh.write(head)
        fh.write(payload)
    return len(MAGIC) + 4 + len(head) + len(payload)


def write_sequence(out: Path, pops: dict, frames: list, header_extra: dict | None = None) -> int:
    """Write a POSITION SEQUENCE: topology once, positions per frame.

    The saving is the whole design. A frame of 4.59 M positions is 55 MB; the topology is 39 MB and
    does not move, so writing it once per file rather than once per frame is the difference between
    a sequence being affordable and not. It also keeps the viewer honest: a per-frame topology would
    imply the connectivity is time-dependent, and it is not — the arena's claims are fixed for the
    life of a run and that invariant is what `assert_partitioned` protects.

    Args:
        out: the ``.alephcell`` path.
        pops: population name -> builder result, for the topology.
        frames: ``[(step_index, {population: (N, 3) float32}), ...]`` in order.
        header_extra: merged into the header.

    Returns:
        Bytes written.

    Raises:
        ValueError: on an empty frame list, or a frame whose populations differ from the first —
            a sequence whose membership changes between frames is not one cell over time.
    """
    import numpy as np

    if not frames:
        raise ValueError("refused: a sequence with no frames. An empty file teaches a viewer nothing.")
    first = set(frames[0][1])
    for i, (_step, f) in enumerate(frames):
        if set(f) != first:
            raise ValueError(
                f"frame {i} holds populations {sorted(set(f))} against frame 0's {sorted(first)}. A "
                "sequence whose membership changes between frames is not one cell over time.")

    entries: list[dict] = []
    blocks: list = []
    offset = 0

    def _add(name: str, arr, kind: str, frame: int | None = None) -> None:
        nonlocal offset
        a = np.ascontiguousarray(arr)
        e = {"population": name, "kind": kind, "dtype": str(a.dtype), "shape": list(a.shape),
             "offset": offset, "nbytes": int(a.nbytes)}
        if frame is not None:
            e["frame"] = int(frame)
        entries.append(e)
        blocks.append(a)
        offset += int(a.nbytes)

    # topology first, once
    for name, obj in pops.items():
        lo, hi = (obj["claims"]["node"] if isinstance(obj, dict) else (obj.nodes.lo, obj.nodes.hi))
        n = int(hi - lo)
        if n == 0:
            continue
        if not isinstance(obj, dict) and hasattr(obj, "face_idx"):
            fi = obj.face_idx
            f = np.asarray(fi.numpy() if hasattr(fi, "numpy") else fi, np.int64).reshape(-1, 3) - lo
            _add(name, f.astype(np.int32), "faces")
            continue
        if isinstance(obj, dict) and "nodes_per_strand" in obj:
            n_str, n_per = int(obj["n_strands"]), int(obj["nodes_per_strand"])
            k = np.arange(n_str, dtype=np.int64).repeat(n_per - 1)
            i = np.tile(np.arange(n_per - 1, dtype=np.int64), n_str)
            pairs = np.stack([k * n_per + i, k * n_per + i + 1], axis=1)
        elif not isinstance(obj, dict) and getattr(obj, "seg_node", None) is not None:
            sn = obj.seg_node
            pairs = np.asarray(sn.numpy() if hasattr(sn, "numpy") else sn, np.int64) - lo
        else:
            continue                       # a population with no drawable topology is skipped, not faked
        _add(name, pairs.astype(np.int32), "segments")

    for idx, (_step, f) in enumerate(frames):
        for name in sorted(f):
            _add(name, f[name], "positions_um", frame=idx)

    header = {
        "format": "aleph-cell-binary", "version": FORMAT_VERSION,
        "n_frames": len(frames), "frame_steps": [int(s) for s, _ in frames],
        "thinning": "NONE. Every node of every population is present in every frame.",
        "topology_written": "ONCE — the arena's claims are fixed for the life of a run, so the "
                            "connectivity is not time-dependent and is not re-written per frame",
        "positions_dtype": "float32",
        **(header_extra or {}),
        "blocks": entries,
    }
    out.parent.mkdir(parents=True, exist_ok=True)
    return _write(out, header, blocks)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--device", default="cuda:0")
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--nmii-heads-per-side", type=int, default=None,
                    help="draw NMII at this head count. NO DEFAULT — decision-queue item 9's band is "
                         "unratified, so the file states which test point it is a picture OF.")
    ap.add_argument("--seg-um", type=float, default=0.05)
    args = ap.parse_args()

    import numpy as np
    import warp as wp

    from aleph.world.arena import Kind, WorldArena
    from aleph.world.bond import BondCount, SourceClass
    from aleph.world.build import build_all, cortex_shell
    from aleph.world.populations import build_remaining_populations
    from aleph.world.build.filopodium import build_filopodia
    from aleph.world.build.intermediate_filament import build_intermediate_filaments
    from aleph.world.build.lamellipodium import build_lamellipodium
    from aleph.world.build.microtubule import build_microtubules
    from aleph.world.build.nmii import _anchor_beads, build_nmii
    from aleph.world.build.stress_fiber import build_stress_fibers
    from aleph.world.geometry import assert_inside_membrane, footprint

    wp.init()
    dev = wp.get_device(args.device)
    if not dev.is_cuda:
        raise RuntimeError(f"native-only; {args.device!r} resolved to {dev}, not CUDA")

    t0 = time.perf_counter()

    def gap(n: float, what: str, basis: str = "explicit") -> BondCount:
        return BondCount(value=float(n), basis=basis, source_class=SourceClass.PI_GAP,
                         scope=f"NO VALUE EXISTS — per-cell structure count absent for {what}",
                         provenance="PI-GAP; a size placeholder, not a claim")

    arena = WorldArena(capacity={Kind.NODE: 12_000_000, Kind.SEGMENT: 12_000_000,
                                 Kind.ANGLE3: 12_000_000, Kind.ANGLE4: 2_000_000,
                                 Kind.FACE: 1_500_000, Kind.STRAND: 200_000,
                                 Kind.BOND: 1_000_000, Kind.GRID_CELL: 4_000_000},
                       device=args.device)
    seg = args.seg_um
    cell = build_all(args.device, arena=arena)
    fp = footprint(7.5)
    # ⚠ ONE copy of the population build, shared with the PHASE 1 and PHASE 4 drivers.
    # This file used to carry its own, and the two copies drifted from the cell they were
    # placing — 95.5% of the stress fibres ended up outside the membrane while
    # assert_partitioned reported success, because partition is an ID-range property.
    pops: dict[str, object] = {
        "cortex": cell.cortex, "membrane": cell.membrane,
        "nuclear_envelope": cell.envelope,
        **build_remaining_populations(arena, fp, seg_um=seg),
    }

    if args.nmii_heads_per_side is not None:
        pops["nmii"] = build_nmii(
            arena, count=BondCount(value=0.625, basis="areal", source_class=SourceClass.PI_GAP,
                                   scope="MCF7 cortical shell, resting",
                                   provenance="Nie 2015 areal density 0.625 um^-2"),
            support=706.86, n_bb=14, n_heads_per_side=int(args.nmii_heads_per_side),
            backbone_length_um=0.301, head_offset_um=0.200,
            radius_um=cortex_shell()[0], thickness_um=cortex_shell()[1])

    arena.assert_partitioned()
    placement = assert_inside_membrane(arena.node_arrays["position"].numpy(), pops, fp.r_cell_um)
    build_s = time.perf_counter() - t0

    pos_host = arena.node_arrays["position"].numpy()
    blocks: list = []
    entries: list[dict] = []
    offset = 0
    skipped: list[str] = []
    unrepresentable: list[dict] = []

    def _add(name: str, arr, kind: str) -> None:
        nonlocal offset
        a = np.ascontiguousarray(arr)
        entries.append({"population": name, "kind": kind, "dtype": str(a.dtype),
                        "shape": list(a.shape), "offset": offset, "nbytes": int(a.nbytes)})
        blocks.append(a)
        offset += int(a.nbytes)

    for name, obj in pops.items():
        lo, hi = (obj["claims"]["node"] if isinstance(obj, dict) else (obj.nodes.lo, obj.nodes.hi))
        n = int(hi - lo)
        if n == 0:
            skipped.append(name)          # said, not written empty
            continue
        _add(name, pos_host[lo:hi].astype(np.float32), "positions_um")

        # Segment pairs, LOCAL to this population's block so a viewer cannot address across
        # populations — the conservation check in the Sanity Gate.
        if name == "nmii":
            n_mf, n_bb = int(obj["n_minifilaments"]), int(obj["n_bb"])
            h, n_per = int(obj["n_heads_per_side"]), int(obj["nodes_per_minifilament"])
            plus, minus = _anchor_beads(n_bb, h)
            local = ([(i, i + 1) for i in range(n_bb - 1)]
                     + [(plus[i], n_bb + i) for i in range(h)]
                     + [(minus[i], n_bb + h + i) for i in range(h)])
            base = np.arange(n_mf, dtype=np.int64) * n_per
            pairs = (base[:, None, None] + np.asarray(local, np.int64)[None]).reshape(-1, 2)
        elif not isinstance(obj, dict) and hasattr(obj, "face_idx"):
            fi = obj.face_idx
            f = np.asarray(fi.numpy() if hasattr(fi, "numpy") else fi, np.int64).reshape(-1, 3) - lo
            _add(name, f.astype(np.int32), "faces")      # a surface is a MESH, not a wire cage
            continue
        elif isinstance(obj, dict) and "n_strands" not in obj:
            # ⚠ A dict is NOT necessarily a strand population, and this branch assumed it was.
            # `lamina` and `chromatin` carry no `n_strands`/`nodes_per_strand` — they describe a
            # shell of filaments (`n_filaments_implied`) and a chain of subunits (`n_subunits`,
            # `n_backbone_segments`) — so the export died with `KeyError: 'n_strands'` on the
            # ELEVENTH population. Found 2026-08-21 by exporting the full cell for the first time:
            # every export before it had ten populations, and the branch had never met the other two.
            #
            # ⚠ REFUSED, not guessed. Inventing a topology for them would put lines in the viewer
            # that the build does not have, and a viewer that invents connectivity can be quoted
            # about structure. And NOT merged with `populations_skipped_empty`: "there was nothing
            # to draw" and "this tool cannot represent what is there" are different facts with
            # different owners, and collapsing them sends the fix to the wrong place.
            unrepresentable.append({
                "population": name, "n_nodes": int(n),
                "why": ("POSITIONS ARE WRITTEN; TOPOLOGY IS NOT. The builder's record carries no "
                        "n_strands/nodes_per_strand, and this exporter can only derive segment "
                        "pairs from those. The nodes are in the file and nothing connects them."),
                "keys_it_does_have": sorted(k for k in obj if k.startswith("n_")),
            })
            continue
        elif isinstance(obj, dict):
            n_str, n_per = int(obj["n_strands"]), int(obj["nodes_per_strand"])
            k = np.arange(n_str, dtype=np.int64).repeat(n_per - 1)
            i = np.tile(np.arange(n_per - 1, dtype=np.int64), n_str)
            pairs = np.stack([k * n_per + i, k * n_per + i + 1], axis=1)
        else:
            sn = obj.seg_node
            pairs = np.asarray(sn.numpy() if hasattr(sn, "numpy") else sn, np.int64) - lo

        if pairs.size and (pairs.min() < 0 or pairs.max() >= n):
            raise RuntimeError(
                f"{name}: a segment index leaves the population block "
                f"([{pairs.min()}, {pairs.max()}] against {n} nodes). A viewer reading this would draw "
                "a line to an unrelated structure, which is worse than not drawing it."
            )
        _add(name, pairs.astype(np.int32), "segments")

    if not entries:
        raise RuntimeError("refused: nothing to export. An empty file teaches a viewer nothing.")
    if unrepresentable:
        names = ", ".join(u["population"] for u in unrepresentable)
        print(f"[export] ⚠ {len(unrepresentable)} population(s) written WITHOUT TOPOLOGY — {names}. "
              "Their positions are in the file; nothing connects them, so a viewer can only draw "
              "them as points and this project does not draw points. The header names them and the "
              "keys they do carry. A picture from this file is not the whole cell.", flush=True)

    header = {
        "format": "aleph-cell-binary",
        "version": FORMAT_VERSION,
        "device": str(dev),
        # ⚠ This was `max(...)` — the LARGEST population, not the total, under the name
        # `n_nodes_total`. It reported 4,197,654 (the cortex alone) while the summary line
        # printed by the same command said 4,588,667, and nothing reconciled the two. Found
        # 2026-08-21 by reading the viewer's own output. The max is kept, under a name that
        # says what it is, because a reader sizing a buffer does want it.
        "n_nodes_total": int(sum(e["shape"][0] for e in entries if e["kind"] == "positions_um")),
        "n_nodes_largest_population": int(
            max(e["shape"][0] for e in entries if e["kind"] == "positions_um")),
        "positions_dtype": "float32",
        "positions_precision_note": (
            "the arena is float64; positions are stored float32. The cell is ~15 um across and "
            "float32 holds ~7 decimal digits, so the worst placement error is ~1e-6 um — four orders "
            "below the 0.05 um cortex mesh. Nothing quantitative is computed from this file."),
        "thinning": "NONE. Every node and every segment of every population is present.",
        "populations_skipped_empty": skipped,
        # ⚠ A SEPARATE list from `skipped`, on purpose. Empty means "nothing was built there";
        # unrepresentable means "something is built there and this tool cannot write it". A reader
        # who merges them concludes the cell has fewer parts than it has.
        "populations_unrepresentable": unrepresentable,
        "completeness": ("EVERY population has positions AND topology" if not unrepresentable else
                         f"⚠ INCOMPLETE — {len(unrepresentable)} population(s) have POSITIONS BUT NO "
                         "TOPOLOGY; see populations_unrepresentable"),
        "nmii_heads_per_side": args.nmii_heads_per_side,
        "nmii_head_count_status": (
            None if args.nmii_heads_per_side is None else
            "UNRATIFIED test point — decision-queue item 9. Neither candidate traces to a primary "
            "measurement passing the citation audit."),
        "placement_rule": "every node lies inside the built membrane (PI 2026-08-21)",
        "placement_check": placement,
        "placement_envelope": fp.as_record(),
        "build_s": round(build_s, 3),
        "blocks": entries,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    n_bytes = _write(args.out, header, blocks)
    n_nodes = sum(e["shape"][0] for e in entries if e["kind"] == "positions_um")
    n_lines = sum(e["shape"][0] for e in entries if e["kind"] == "segments")
    n_faces = sum(e["shape"][0] for e in entries if e["kind"] == "faces")
    print(f"[export] {n_nodes:,} nodes · {n_lines:,} segments · {n_faces:,} faces · "
          f"{n_bytes / 1e6:.1f} MB · none thinned · wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
