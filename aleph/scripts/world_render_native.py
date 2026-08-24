"""Render the standing native cell to a browser page — every node, nothing sampled away.

WHY A DRIVER RATHER THAN A CHANGE TO THE VIEWER.  ``world/viewer.py`` takes ``Strand`` objects, one per
filament, which is the model the arena had before PHASE 1 replaced per-strand claims with per-population
ones. Its reader needs only ``population`` / ``nodes`` / ``seg_node``, and a ``StrandPopulation`` carries
all three — so this adapts rather than editing a module two sessions are writing into.

⚠ THE PAGE IS THE LIMIT, NEVER A STRIDE.  ``viewer.py``'s own rule: *"When a scene grows past what one
page can hold, the honest limit is the PAGE — a byte budget, stated — and never a stride through the
population."* A 4.4 M-node cortex at 16 B/node is ~300 MB of base64, which no browser tab opens. So this
driver draws WHOLE POPULATIONS or none: ``--only`` selects which populations are in the page, the page
states which were left out by name, and no population is ever thinned. A picture of three complete
populations is a true picture of three things; a picture of every population at one node in ten is a
picture of nothing.

Sanity Gate (recorded per CLAUDE.md, before first execution):
  * dimensional — positions [µm] throughout; the viewer derives its scale bar from the payload extent.
  * boundary — a request naming a population that was not built is refused by name rather than skipped;
    an empty selection is refused rather than emitting an empty page.
  * conservation/invariant — the node count written into the page is read back from the arena census,
    not from what this driver believes it passed.
  * CFL/precision — nothing is integrated. float64 positions narrowed to float32 for transport, ~1e-7
    relative, four orders below the µm features drawn, and the viewer states it on the page.
  * sign sense — not applicable.
  * measurement protocol — ONE host readback of the built arrays, between builds, never inside a loop.
    There is no loop here: the cell is standing, not stepping.

engine units: length µm.  Runtime: builds on CUDA, renders on the host.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--device", default="cuda:0")
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--only", default="membrane,nuclear_envelope,microtubule,intermediate_filament",
                    help="comma-separated populations to draw IN FULL. Everything not named is left "
                         "out entirely and said so on the page — never thinned.")
    ap.add_argument("--nmii-heads-per-side", type=int, default=None,
                    help="draw the NMII minifilaments at this head count. NO DEFAULT: it is "
                         "decision-queue item 9 and unratified, so the picture states which "
                         "test point it is a picture OF. Omitted, NMII is not drawn.")
    ap.add_argument("--seg-um", type=float, default=0.05)
    args = ap.parse_args()

    import numpy as np
    import warp as wp

    from aleph.world.arena import Kind, WorldArena
    from aleph.world.bond import BondCount, SourceClass
    from aleph.world.geometry import assert_inside_membrane, footprint
    from aleph.world.build import build_all, cortex_shell
    from aleph.world.populations import build_remaining_populations
    from aleph.world.build.filopodium import build_filopodia
    from aleph.world.build.intermediate_filament import build_intermediate_filaments
    from aleph.world.build.lamellipodium import build_lamellipodium
    from aleph.world.build.microtubule import build_microtubules
    from aleph.world.build.nmii import _anchor_beads, build_nmii
    from aleph.world.build.stress_fiber import build_stress_fibers
    from aleph.world.viewer import render_cell_view

    wp.init()
    if not wp.get_device(args.device).is_cuda:
        raise RuntimeError(f"{args.device!r} is not CUDA; this renders the NATIVE cell only")

    def gap(n: float, what: str, basis: str = "explicit") -> BondCount:
        return BondCount(value=float(n), basis=basis, source_class=SourceClass.PI_GAP,
                         scope=f"NO VALUE EXISTS — per-cell structure count absent for {what}",
                         provenance="PI-GAP 2026-08-20; a size placeholder, not a claim")

    def derived(n: int, what: str) -> BondCount:
        """A count that FOLLOWS from geometry, so it is not a gap and must not be labelled one."""
        return BondCount(value=float(n), basis="explicit", source_class=SourceClass.DERIVED,
                         scope="MCF7 basal footprint AS BUILT — R_cell 7.5 um, substrate z -7.0 um",
                         provenance=f"DERIVED from the contact disc: {what} = "
                                    "floor(2 r_footprint / pitch) + 1, the pitch held at its "
                                    "physiological value so the sourced quantity stays sourced")

    arena = WorldArena(capacity={Kind.NODE: 12_000_000, Kind.SEGMENT: 12_000_000,
                                 Kind.ANGLE3: 12_000_000, Kind.ANGLE4: 2_000_000,
                                 Kind.FACE: 1_500_000, Kind.STRAND: 200_000,
                                 Kind.BOND: 1_000_000}, device=args.device)
    seg = args.seg_um
    cell = build_all(args.device, arena=arena)
    fp = footprint(7.5)      # placements are DERIVED — see aleph/world/geometry.py
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
            arena,
            # ⚠ the minifilament COUNT is not the gap here — Nie 2015's areal density is the datum and
            # 442 falls out of it. What is unratified is n_heads_per_side, which the flag carries.
            count=BondCount(value=0.625, basis="areal", source_class=SourceClass.PI_GAP,
                            scope="MCF7 cortical shell, resting",
                            provenance="Nie 2015 areal density 0.625 um^-2 — the density the frozen "
                                       "incumbent already carried"),
            support=706.86,
            n_bb=14, n_heads_per_side=int(args.nmii_heads_per_side),
            backbone_length_um=0.301, head_offset_um=0.200,
            radius_um=cortex_shell()[0], thickness_um=cortex_shell()[1])
    arena.assert_partitioned()

    want = [w.strip() for w in args.only.split(",") if w.strip()]
    unknown = [w for w in want if w not in pops]
    if unknown:
        raise SystemExit(f"refused: {unknown} were never built. Built: {sorted(pops)}")
    if not want:
        raise SystemExit("refused: an empty selection would emit an empty page")
    left_out = sorted(set(pops) - set(want))

    class _Range:
        __slots__ = ("lo", "hi")

        def __init__(self, lo, hi):
            self.lo, self.hi = int(lo), int(hi)

    class _StrandView:
        """The three fields viewer.py reads off a strand, served from a population.

        Two shapes arrive here: A1's builders return an object with ``.nodes`` and a device
        ``.seg_node``; A2's return a dict and build NO topology array at all, because their strands are
        uniform chains and segment ``s`` of strand ``k`` is ``lo + k*n_per + i`` by arithmetic. That is
        a real saving at 5 M nodes — an int64 pair per segment is not stored — so the adapter DERIVES
        the pairs for drawing rather than asking the builder to have materialised them.
        """

        def __init__(self, name, obj, pos_host):
            self.population = name
            if isinstance(obj, dict):
                lo, hi = obj["claims"]["node"]
                self.nodes = _Range(lo, hi)
                n_str, n_per = int(obj["n_strands"]), int(obj["nodes_per_strand"])
                k = np.arange(n_str, dtype=np.int64).repeat(n_per - 1)
                i = np.tile(np.arange(n_per - 1, dtype=np.int64), n_str)
                a = lo + k * n_per + i
                self.seg_node = np.stack([a, a + 1], axis=1)
            else:
                self.nodes = obj.nodes
                sn = obj.seg_node
                self.seg_node = np.asarray(sn.numpy() if hasattr(sn, "numpy") else sn)
            self.position = pos_host[self.nodes.lo:self.nodes.hi]
            self.n_nodes = self.nodes.hi - self.nodes.lo

    class _NmiiView:
        """NMII is the one population whose strand is NOT a chain, so it gets its own adapter.

        ⚠ A bipolar minifilament is ``[bb_0 … bb_{n_bb-1}, head+_0 … head+_{H-1}, head-_0 … head-_{H-1}]``.
        Nodes past ``n_bb`` are head ARMS hanging off backbone beads, not a continuation of the chain.
        Feeding this to ``_StrandView``'s ``(a, a+1)`` derivation — correct for every other population
        here — would draw a spurious segment from the last backbone bead into the first head and then
        thread a line through all twenty heads. **The picture would be of a structure that does not
        exist**, and a picture is what gets looked at.

        The anchor arithmetic is NOT restated here: :func:`aleph.world.build.nmii._anchor_beads` is the
        module's own host mirror of the kernel, checked against it there. Two copies of index
        arithmetic is how a renderer drifts from the thing it renders.
        """

        def __init__(self, obj, pos_host):
            self.population = "nmii"
            lo, hi = obj["claims"]["node"]
            self.nodes = _Range(lo, hi)
            n_mf, n_bb = int(obj["n_minifilaments"]), int(obj["n_bb"])
            h, n_per = int(obj["n_heads_per_side"]), int(obj["nodes_per_minifilament"])
            plus, minus = _anchor_beads(n_bb, h)
            local = [(i, i + 1) for i in range(n_bb - 1)]                       # the backbone chain
            local += [(plus[i], n_bb + i) for i in range(h)]                    # + side head arms
            local += [(minus[i], n_bb + h + i) for i in range(h)]               # - side head arms
            base = lo + np.arange(n_mf, dtype=np.int64) * n_per
            pair = np.asarray(local, dtype=np.int64)
            self.seg_node = (base[:, None, None] + pair[None, :, :]).reshape(-1, 2)
            self.position = pos_host[lo:hi]
            self.n_nodes = hi - lo

    class _SurfaceView:
        """What viewer.py reads off a surface, served from a ClosedSurface.

        The population keeps ``face_idx`` on the DEVICE and holds no host ``position`` — the positions
        live in the arena's shared array, which is the whole point of the arena. So this reads both
        back once, here, between builds and never inside a loop.
        """

        def __init__(self, obj, pos_host):
            self.population = obj.population
            self.nodes = obj.nodes
            fi = obj.face_idx
            self.face_idx = np.asarray(fi.numpy() if hasattr(fi, "numpy") else fi)
            self.position = pos_host[self.nodes.lo:self.nodes.hi]
            self.n_nodes = self.nodes.hi - self.nodes.lo

    pos_host = arena.node_arrays["position"].numpy()
    strands, surfaces = [], []
    for name in want:
        obj = pops[name]
        if name == "nmii":
            strands.append(_NmiiView(obj, pos_host))
        elif not isinstance(obj, dict) and hasattr(obj, "face_idx"):
            surfaces.append(_SurfaceView(obj, pos_host))
        else:
            strands.append(_StrandView(name, obj, pos_host))

    note = (f"PHASE 1, {args.device}. Drawn IN FULL: {', '.join(want)}. "
            + (f"LEFT OUT ENTIRELY (not thinned): {', '.join(left_out)}. " if left_out else "")
            + "Five of the eight populations carry PI_GAP counts — no sourced per-cell structure count "
              "exists for MT, IF, filopodium, lamellipodium or stress fibre. Geometry only: no force, "
              "no law, nothing has moved.")
    html = render_cell_view(arena, strands=strands, surfaces=surfaces,
                            title="Aleph arena — native cell, PHASE 1", note=note)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(html)
    census = arena.census()
    print(json.dumps({"page_bytes": len(html), "drawn": want, "left_out": left_out,
                      "census": census}, indent=1, default=str), flush=True)
    print(f"\nwrote {args.out}  ({len(html)/2**20:.1f} MB)", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
