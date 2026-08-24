r"""Stress-fiber prestress + traction-loop wiring — the emergent contractile bundle load path — I4 / P3.

THE SF WIRING GAP (FILAMENT_SUBSYSTEM_PLAN §3 G1/G2, [[project-filament-subsystem-plan]]). ``ac.weave.regions``
already SEEDS ventral / dorsal SF + transverse arc + perinuclear cap as isotropic-conditional bundles (fibers +
alpha-actinin + NMII motor seeds + FA/LINC anchor sites) and ``weave_cell`` concatenates them into the one
network. What was still missing — the "biggest structural gap" of the walkthrough — is that a stress fiber is
not merely a bundle of fibers: it is a **pre-stressed contractile load path**. Its defining property (the reason
it is an SF and not a random crosslinked patch) is the AXIAL PRESTRESS that the NMII bipolar motors generate and
that the FA anchors react as a traction dipole into the substrate. This module is that representation + the
hand-off that puts the SF into the whole-cell loop.

EMERGENT, NOT LUMPED (CLAUDE.md hard rule; [[project-unified-actin-architecture]]). The axial tension is NOT a
lumped ``k_SF`` bundle stiffness. It is computed from the DISCRETE NMII head reactions (each head walks toward
its actin barbed end -> Newton-3rd reaction on the actin; an antiparallel bipolar minifilament -> net inward
contraction) accumulated along the bundle axis as the load-path integral ``T(s) = -sum_{s'>s} f_axial(s')`` (the
same "actin axial tension = cumulative sum of external forces" pattern the gamma harness uses,
[[project-ac-magnitude-integration-blockers]]). The bundle CONDENSATION that sharpens this prestress is the I5
native proof; here we deliver the buildable load-path representation + the magnitude-independent SHAPE gate.

EMERGENCE FIREWALL (build-plan §5-confirm-3). This module NEVER reads ``region_id`` (the SF/cap/arc diagnostic
label). A stress fiber is discovered label-blind: a **motorized, anchored** connected bundle. The emergent
taxonomy falls out of geometry + mechanics (INV-4/INV-5), it is not a hand-set class:
  * contractile (has NMII motors) + FA-anchored  -> ventral / dorsal SF
  * contractile + LINC-anchored                  -> perinuclear cap SF (I7 nucleus load path)
  * anchored, NO motor                           -> filopodium (fascin bundle, non-contractile)
An un-anchored motorized region (the cortex) is NOT a stress fiber and is excluded (no FA/LINC anchor).

MAGNITUDE IS A GAP — report-not-tune (CLAUDE.md; [[project-sf-nmii-forcescale-result]]). The per-head stall
force and the number of engaged heads per cross-section (the bundle density) are I0-B3 / I0-B6 GAPs; the single-
SF ~5-6 nN scale is density-floored and NOT closable from the coarse seed. So ``f_head_pN`` defaults to ``None``
and the load path returns its magnitude-independent SHAPE (axis, contractility sign, cross-section ordering, peak
location) with ``tension = None`` and ``magnitude_status = "GAP"``. A caller may pass a PROVISIONAL per-head
value to get a density-floored FINDING profile — explicitly labeled coarse-seed, never tuned to a traction band.

engine units: length um, force pN, stiffness pN/um. Host NumPy only (no Warp import); the native force assembly
(reusing the I3 ``ac.motor`` Hill-NMII kernels) is the lead's on the gbook A5000 — see ``INTEGRATION.md`` §3g.

Sanity Gate (self-tested in tests/ac/weave/test_stress_fiber_oracle.py):
  * dimensional/convention: axis is a unit vector; ``tension_per_fhead`` is dimensionless (pN per pN-head);
    ``tension`` (when ``f_head`` given) is in pN; the cumulative load path is exact (endpoint tension = 0 for a
    balanced isolated SF, machine precision).
  * sign-sense: an antiparallel FA<->FA bundle with barbed ends OUTWARD is CONTRACTILE (peak internal tension
    > 0, dipole pulls the two anchors together); a parallel (non-sarcomeric) bundle is not net-contractile.
  * measurement-protocol: the ventral (both-ends) SF prestress peaks in the SHAFT (interior); a dorsal
    (one-end) SF is asymmetric (peak shifted to the anchored end) — the tip/shaft/basal load-path partition.
  * firewall: the load path is a pure function of geometry + motors + anchors and never reads ``region_id``
    (grep-checked in the test); ``extract_sf_bundles`` groups label-blind by anchor+motor connectivity.
  * ledger: unique SF-fiber IDs, no cortex/arc/SF double count; bundles partition the anchored-motorized set.
  * GAP guard: ``f_head_pN=None`` -> ``tension is None`` and ``magnitude_status == "GAP"`` (never a tuned number).
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import numpy.typing as npt

from aleph.components.emergence.nematic import fiber_axes, nematic_order_and_director
from aleph.components.weave.walk_dir import barbed_end_node

__all__ = [
    "GAP_STATUS",
    "OK_STATUS",
    "StressFiberLoadPath",
    "WovenCellSFWiring",
    "stress_fiber_load_path",
    "extract_sf_bundles",
    "wire_stress_fibers",
]

_EPS = 1.0e-12

#: magnitude-status tags (report-not-tune): the shape gate is magnitude-independent; the magnitude is GAP until
#: the I0-B3 per-head stall + I0-B6 bundle density close (surfaced to PI, never tuned to a traction band).
GAP_STATUS = "GAP — NMII per-head f_stall (I0-B3) + engaged-head density (I0-B6) unresolved; surface to PI"
OK_STATUS = "FINDING (density-floored coarse seed; NOT the native SF magnitude) — report-not-tune"


@dataclass(slots=True)
class StressFiberLoadPath:
    """The contractile load path of ONE label-blind bundle (SF / cap SF), computed from its discrete motors.

    All arrays index the bundle's own fibers/nodes (local). ``tension`` is ``None`` (GAP) unless a provisional
    per-head force was supplied; ``tension_per_fhead`` is the magnitude-INDEPENDENT shape (always populated).
    """

    axis: npt.NDArray[np.float64]                       # (3,) unit contractile axis (nematic director, sign-fixed)
    nematic_order: float                                # bundle alignment S in [0,1] (0=isotropic seed, ->1 condensed)
    n_fibers: int
    n_motors: int
    contractile: bool                                   # has NMII motors (vs a non-contractile fascin bundle)
    anchor_kind: str                                    # "fa" | "linc" | "fa+linc" | "none"
    s_bins: npt.NDArray[np.float64]                     # (B,) cross-section axial coordinate [um]
    tension_per_fhead: npt.NDArray[np.float64]          # (B,) internal axial tension / f_head [pN per pN] (shape)
    peak_frac: float                                    # argmax(tension) position along the axis in [0,1]
    dipole_per_fhead: float                             # traction-dipole magnitude / f_head [pN per pN]
    dipole_axis: npt.NDArray[np.float64]                # (3,) unit dipole direction (anchor_a -> anchor_b)
    balance_residual_per_fhead: float                   # |sum of end reactions| / f_head (0 => closed both-end SF)
    magnitude_status: str                               # GAP_STATUS or OK_STATUS
    f_head_pN: float | None = None                      # provisional per-head stall passed IN (never chosen here)
    tension: npt.NDArray[np.float64] | None = None      # (B,) internal axial tension [pN] — None if GAP
    dipole_pN: float | None = None                      # traction-dipole magnitude [pN] — None if GAP

    def load_path_thirds(self) -> dict[str, float]:
        """Mean prestress (per f_head) in the three axial thirds of the bundle — the load-path partition.

        A closed both-ends (ventral) SF carries prestress across its SHAFT (interior third); a one-end (dorsal)
        SF is asymmetric (prestress funnels toward its single anchored end). This is the geometric prestress
        distribution — the tip/shaft/basal ADHESION-maturation classification it feeds is I6's (``fa_maturation``),
        not set here.
        """
        b = self.tension_per_fhead
        n = b.shape[0]
        third = max(n // 3, 1)
        return {
            "near_a": float(b[:third].mean()),
            "shaft": float(b[third:n - third].mean()) if n - 2 * third > 0 else float(b.mean()),
            "near_b": float(b[-third:].mean()),
        }

    @property
    def closed_dipole(self) -> bool:
        """True if the two anchored ends balance (a self-equilibrated contractile SF, e.g. ventral FA<->FA)."""
        return self.dipole_per_fhead > _EPS and self.balance_residual_per_fhead <= 0.5 * self.dipole_per_fhead


def _resolve_axis(
    pos: npt.NDArray[np.float64], fiber_offsets: npt.NDArray[np.int64],
    anchor_nodes: npt.NDArray[np.int64],
) -> tuple[npt.NDArray[np.float64], float]:
    """The bundle contractile axis = nematic director, sign-fixed toward the anchor spread (label-blind).

    The nematic director is n -> -n symmetric; a stress fiber has a well-defined contractile ORIENTATION but
    the tip/shaft/basal ordering needs a sign, so we orient the axis along the anchor cloud's principal spread
    (FA<->FA for a ventral SF). With < 2 anchors the sign is arbitrary (kept as the director's).
    """
    axes = fiber_axes(pos, fiber_offsets)
    s, director = nematic_order_and_director(axes)
    if anchor_nodes.size >= 2:
        a = pos[anchor_nodes]
        spread = a - a.mean(axis=0)
        # principal spread direction of the anchor cloud
        _, evec = np.linalg.eigh(spread.T @ spread)
        anchor_dir = evec[:, -1]
        if float(director @ anchor_dir) < 0.0:
            director = -director
    return director, s


def _motor_node_reactions(
    pos: npt.NDArray[np.float64], myo_i: npt.NDArray[np.int64], myo_j: npt.NDArray[np.int64],
    node_fiber: npt.NDArray[np.int64], barbed_node: npt.NDArray[np.int64],
) -> npt.NDArray[np.float64]:
    """Per-node UNIT-f_head reaction force of every NMII motor (Newton-3rd of the barbed-directed power stroke).

    A head bound to actin node ``a`` on a filament with barbed end ``b`` walks toward ``b``; the reaction on the
    actin is ``-f_head * unit(pos[b]-pos[a])`` (the actin is pulled the opposite way — toward the minifilament
    centre). A bipolar minifilament straddling two anti-parallel filaments (barbed ends outward) therefore pulls
    both pointed ends inward -> net contraction. Returns forces in units of ``f_head`` (magnitude factored out).
    """
    n = pos.shape[0]
    f = np.zeros((n, 3), np.float64)
    for a, bnode_side in ((myo_i, myo_j), (myo_j, myo_i)):
        if a.size == 0:
            continue
        b = barbed_node[node_fiber[a]]
        v = pos[b] - pos[a]
        nv = np.linalg.norm(v, axis=1)
        good = nv > _EPS
        wdir = np.zeros_like(v)
        wdir[good] = v[good] / nv[good, None]
        np.add.at(f, a, -wdir)          # reaction: actin pulled opposite the head's walk direction
    return f


def stress_fiber_load_path(
    pos: npt.NDArray[np.float64],
    fiber_offsets: npt.NDArray[np.int64],
    polarity: npt.NDArray[np.int64],
    myo_i: npt.NDArray[np.int64],
    myo_j: npt.NDArray[np.int64],
    fa_sites: npt.NDArray[np.int64],
    linc_sites: npt.NDArray[np.int64],
    *,
    f_head_pN: float | None = None,
    n_bins: int = 24,
) -> StressFiberLoadPath:
    """Compute a bundle's emergent contractile load path (label-blind; magnitude GAP-guarded).

    Args:
        pos: (N,3) bundle node positions [um].
        fiber_offsets: (F+1,) contiguous offsets (fiber f = nodes [off[f], off[f+1])).
        polarity: (F,) barbed-end flag per fiber (+1 last node, -1 first node).
        myo_i, myo_j: (M,) NMII motor-seed endpoint node pairs (empty => non-contractile bundle).
        fa_sites, linc_sites: substrate / nucleus anchor node sets (the traction / LINC endpoints).
        f_head_pN: PROVISIONAL per-head stall [pN] to scale the density-floored FINDING; ``None`` => GAP (the
            magnitude-independent SHAPE only). NEVER chosen here to hit a band (report-not-tune).
        n_bins: number of axial cross-sections for the tension profile.

    Returns:
        The :class:`StressFiberLoadPath` (shape always populated; magnitude None unless ``f_head_pN`` given).
    """
    pos = np.asarray(pos, np.float64)
    fiber_offsets = np.asarray(fiber_offsets, np.int64)
    n_fib = fiber_offsets.shape[0] - 1
    anchors = np.concatenate([np.asarray(fa_sites, np.int64), np.asarray(linc_sites, np.int64)]) \
        if (fa_sites.size or linc_sites.size) else np.zeros(0, np.int64)

    axis, s_order = _resolve_axis(pos, fiber_offsets, anchors)
    barbed = barbed_end_node(fiber_offsets, np.asarray(polarity, np.int64))
    node_fiber = np.repeat(np.arange(n_fib), np.diff(fiber_offsets)).astype(np.int64)

    contractile = myo_i.size > 0
    # per-node motor reaction (unit f_head), then its axial projection = the external axial load on the actin
    fvec = _motor_node_reactions(pos, np.asarray(myo_i, np.int64), np.asarray(myo_j, np.int64), node_fiber, barbed)
    f_axial = fvec @ axis                                   # (N,) per-node axial force / f_head
    s_node = (pos - pos.mean(axis=0)) @ axis

    # ---- PER-FIBER load path: tension is carried by the actin BACKBONE between a motor attachment and the
    #      fiber's anchor, NOT across a free-body cut of loose nodes (a single bipolar motor's two heads sit at
    #      one axial station and would spuriously cancel). For each fiber the internal tension at a node =
    #      sum of the axial motor forces on the OUTBOARD side (away from the anchor); the anchor supplies the
    #      balancing reaction (= the traction that fiber delivers). Bundle prestress = sum over fibers.
    s_lo, s_hi = float(s_node.min()), float(s_node.max())
    span = max(s_hi - s_lo, _EPS)
    s_bins = np.linspace(s_lo, s_hi, n_bins)
    tension_per_fhead = np.zeros(n_bins, np.float64)
    fiber_reaction = np.zeros(n_fib, np.float64)            # axial reaction each fiber hands its anchor (/f_head)
    anchor_set = set(int(a) for a in anchors)
    for f in range(n_fib):
        lo, hi = int(fiber_offsets[f]), int(fiber_offsets[f + 1])
        sf = s_node[lo:hi]
        ff = f_axial[lo:hi]
        order = np.argsort(sf)
        sfs, ffs = sf[order], ff[order]
        # anchor reference = the fiber's anchored end (low-|Δ| to an anchor node); else the low-s end
        anchor_local = [k for k in range(lo, hi) if k in anchor_set]
        anchor_at_high = bool(anchor_local) and (s_node[anchor_local[0]] >= 0.5 * (sfs[0] + sfs[-1]))
        # tension = outboard axial force projected on the OUTWARD-from-anchor direction (positive = stretching).
        # anchor at the low end: outward = +axis, T = sum of forces on the high side (cumsum from high).
        # anchor at the high end: outward = -axis, T = -(sum of forces on the low side) (negated cumsum from low).
        t_node = -np.cumsum(ffs) if anchor_at_high else np.cumsum(ffs[::-1])[::-1]
        # accumulate into the shared bundle cross-sections this fiber spans
        tension_per_fhead += np.interp(s_bins, sfs, t_node, left=0.0, right=0.0)
        fiber_reaction[f] = float(ffs.sum())               # net axial load -> balanced by the anchor reaction

    peak_frac = float((s_bins[int(np.argmax(np.abs(tension_per_fhead)))] - s_lo) / span)

    # ---- traction dipole: the net axial reaction delivered to each of the two anchor ends -----------------
    dipole_axis, dipole_per_fhead, balance_residual = _anchor_dipole(
        pos, anchors, fiber_reaction, axis, node_fiber)

    if fa_sites.size and linc_sites.size:
        anchor_kind = "fa+linc"
    elif fa_sites.size:
        anchor_kind = "fa"
    elif linc_sites.size:
        anchor_kind = "linc"
    else:
        anchor_kind = "none"

    status = OK_STATUS if f_head_pN is not None else GAP_STATUS
    tension = tension_per_fhead * f_head_pN if f_head_pN is not None else None
    dipole_pN = dipole_per_fhead * f_head_pN if f_head_pN is not None else None

    return StressFiberLoadPath(
        axis=axis, nematic_order=s_order, n_fibers=n_fib, n_motors=int(myo_i.size), contractile=contractile,
        anchor_kind=anchor_kind, s_bins=s_bins, tension_per_fhead=tension_per_fhead, peak_frac=peak_frac,
        dipole_per_fhead=dipole_per_fhead, dipole_axis=dipole_axis,
        balance_residual_per_fhead=balance_residual, magnitude_status=status,
        f_head_pN=f_head_pN, tension=tension, dipole_pN=dipole_pN)


def _anchor_dipole(
    pos: npt.NDArray[np.float64], anchors: npt.NDArray[np.int64], fiber_reaction: npt.NDArray[np.float64],
    axis: npt.NDArray[np.float64], node_fiber: npt.NDArray[np.int64],
) -> tuple[npt.NDArray[np.float64], float, float]:
    """Traction-dipole unit axis + per-f_head magnitude + balance residual = the SF's entry to the whole-cell loop.

    Split the anchored fibers into the two AXIAL ends (by their anchor coordinate about the median) and sum the
    per-fiber reaction each end carries. ``dipole`` = the larger end's contractile pull (the traction the SF
    delivers to the substrate/nucleus). ``residual`` = ``|R_lo + R_hi|`` = the net force handed to the rest of the
    network: ~0 for a self-equilibrated both-ends (ventral) SF, ~= the dipole for a one-end (dorsal) SF whose free
    end is carried by the network.
    """
    if anchors.size < 1:
        return axis.copy(), 0.0, 0.0
    centroid = pos.mean(axis=0)
    sa = (pos[anchors] - centroid) @ axis
    med = float(np.median(sa))
    hi_nodes = anchors[sa >= med]
    lo_nodes = anchors[sa < med]
    fibers_hi = np.unique(node_fiber[hi_nodes]) if hi_nodes.size else np.zeros(0, np.int64)
    fibers_lo = np.unique(node_fiber[lo_nodes]) if lo_nodes.size else np.zeros(0, np.int64)
    r_hi = float(fiber_reaction[fibers_hi].sum()) if fibers_hi.size else 0.0
    r_lo = float(fiber_reaction[fibers_lo].sum()) if fibers_lo.size else 0.0
    dipole = max(abs(r_hi), abs(r_lo))
    residual = abs(r_hi + r_lo)
    if hi_nodes.size and lo_nodes.size:
        d = pos[hi_nodes].mean(axis=0) - pos[lo_nodes].mean(axis=0)
        nd = float(np.linalg.norm(d))
        daxis = d / nd if nd > _EPS else axis.copy()
    else:
        daxis = axis.copy()
    return daxis, dipole, residual


# ── whole-cell wiring: label-blind bundle extraction from a WovenCell ─────────────────────────────
def extract_sf_bundles(cell) -> list[dict]:
    """Discover the stress-fiber bundles in a :class:`~ac.weave.woven_cell.WovenCell` — LABEL-BLIND (no region_id).

    A stress fiber is a **motorized + anchored** connected bundle. We build the intra-bundle connectivity graph
    from crosslinks + motors (NOT the named cross-region ``xbond`` bridges), keep only fibers that carry an FA or
    LINC anchor, and return each connected component. The emergent taxonomy (ventral/dorsal SF vs cap SF vs
    filopodium) falls out of anchor-kind + contractility — it is never read from a label.

    Returns:
        One dict per bundle with the bundle's global fiber indices + the local arrays needed by
        :func:`stress_fiber_load_path` (pos/offsets/polarity/motors/anchors remapped to local indices).
    """
    n_fib = cell.n_fibers
    node_fiber = cell.node_fiber

    anchor_nodes = np.concatenate([cell.fa_sites, cell.linc_sites]) if (cell.fa_sites.size or cell.linc_sites.size) \
        else np.zeros(0, np.int64)
    anchored_fiber = np.zeros(n_fib, bool)
    if anchor_nodes.size:
        anchored_fiber[np.unique(node_fiber[anchor_nodes])] = True
    fa_fiber = np.zeros(n_fib, bool)
    linc_fiber = np.zeros(n_fib, bool)
    if cell.fa_sites.size:
        fa_fiber[np.unique(node_fiber[cell.fa_sites])] = True
    if cell.linc_sites.size:
        linc_fiber[np.unique(node_fiber[cell.linc_sites])] = True

    # union-find over intra-bundle edges (crosslinks + motors) restricted to anchored fibers
    parent = np.arange(n_fib)

    def find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(x: int, y: int) -> None:
        rx, ry = find(x), find(y)
        if rx != ry:
            parent[rx] = ry

    for ei, ej in ((cell.xl_i, cell.xl_j), (cell.myo_i, cell.myo_j)):
        if ei.size == 0:
            continue
        fi, fj = node_fiber[ei], node_fiber[ej]
        # vectorized pre-filter: only edges whose BOTH fibers are anchored can be intra-SF (keeps the python
        # union loop tiny even with the full 70,686-fiber cortex, whose fibers carry no FA/LINC anchor).
        keep = anchored_fiber[fi] & anchored_fiber[fj]
        for u, v in zip(fi[keep], fj[keep]):
            union(int(u), int(v))

    # collect components among anchored fibers
    comps: dict[int, list[int]] = {}
    for f in np.where(anchored_fiber)[0]:
        comps.setdefault(find(int(f)), []).append(int(f))

    bundles: list[dict] = []
    for fibs in comps.values():
        fibs = sorted(fibs)
        bundles.append(_slice_bundle(cell, np.asarray(fibs, np.int64), fa_fiber, linc_fiber))
    return bundles


def _slice_bundle(cell, fibs: npt.NDArray[np.int64], fa_fiber, linc_fiber) -> dict:
    """Remap the global fibers ``fibs`` of ``cell`` to a local, contiguous bundle array set."""
    node_lists, foff = [], [0]
    gnode_to_local: dict[int, int] = {}
    for f in fibs:
        lo, hi = int(cell.fiber_offsets[f]), int(cell.fiber_offsets[f + 1])
        for g in range(lo, hi):
            gnode_to_local[g] = len(gnode_to_local)
        node_lists.append(cell.pos[lo:hi])
        foff.append(foff[-1] + (hi - lo))
    pos = np.concatenate(node_lists, axis=0)
    fiber_offsets = np.asarray(foff, np.int64)
    polarity = cell.polarity[fibs]

    def _local_pairs(gi, gj):
        keep = np.array([(int(a) in gnode_to_local and int(b) in gnode_to_local) for a, b in zip(gi, gj)], bool) \
            if gi.size else np.zeros(0, bool)
        if not keep.any():
            return np.zeros(0, np.int64), np.zeros(0, np.int64)
        li = np.array([gnode_to_local[int(a)] for a in gi[keep]], np.int64)
        lj = np.array([gnode_to_local[int(b)] for b in gj[keep]], np.int64)
        return li, lj

    myo_li, myo_lj = _local_pairs(cell.myo_i, cell.myo_j)

    def _local_nodes(gsites):
        if gsites.size == 0:
            return np.zeros(0, np.int64)
        loc = [gnode_to_local[int(g)] for g in gsites if int(g) in gnode_to_local]
        return np.asarray(loc, np.int64)

    fa_loc = _local_nodes(cell.fa_sites)
    linc_loc = _local_nodes(cell.linc_sites)
    return {
        "global_fibers": fibs, "pos": pos, "fiber_offsets": fiber_offsets, "polarity": polarity,
        "myo_i": myo_li, "myo_j": myo_lj, "fa_sites": fa_loc, "linc_sites": linc_loc,
    }


@dataclass(slots=True)
class WovenCellSFWiring:
    """The whole-cell SF wiring product: per-bundle load paths + a firewall-safe unique-ID ledger."""

    load_paths: list[StressFiberLoadPath]
    global_fibers: list[npt.NDArray[np.int64]]      # per-bundle global fiber IDs (partition; no double count)
    f_head_pN: float | None
    magnitude_status: str

    def ledger(self) -> dict[str, object]:
        """Unique-active SF ledger — no cortex/arc/SF double count (each fiber in at most one bundle)."""
        all_fibs = np.concatenate(self.global_fibers) if self.global_fibers else np.zeros(0, np.int64)
        contractile = [lp for lp in self.load_paths if lp.contractile]
        return {
            "N_sf_bundles": len(self.load_paths),
            "N_sf_bundles_contractile": len(contractile),
            "N_unique_sf_fibers": int(np.unique(all_fibs).size),
            "N_sf_fibers_summed": int(all_fibs.size),        # == unique iff bundles truly partition (asserted)
            "anchor_kinds": [lp.anchor_kind for lp in self.load_paths],
            "closed_dipole": [lp.closed_dipole for lp in self.load_paths],
            "magnitude_status": self.magnitude_status,
            "dipole_pN": [lp.dipole_pN for lp in self.load_paths],
        }


def wire_stress_fibers(cell, *, f_head_pN: float | None = None, n_bins: int = 24) -> WovenCellSFWiring:
    """Wire every emergent stress fiber of a woven cell into a contractile load path (label-blind, GAP-guarded).

    This is the SF's entry into the whole-cell loop (FILAMENT_SUBSYSTEM_PLAN G1/G2): each contractile FA-anchored
    bundle becomes a traction dipole the (lead-owned) ff clutch/ECM loop reacts; each LINC-anchored cap SF is the
    I7 nucleus load path. Discovery is label-blind (:func:`extract_sf_bundles`); magnitude is GAP unless a
    provisional per-head force is supplied (report-not-tune).

    Args:
        cell: the :class:`~ac.weave.woven_cell.WovenCell` (consumed read-only; ``region_id`` is NOT read).
        f_head_pN: provisional per-head NMII stall [pN] or ``None`` (GAP — the default).
        n_bins: cross-sections per bundle tension profile.

    Returns:
        The :class:`WovenCellSFWiring`.
    """
    bundles = extract_sf_bundles(cell)
    load_paths, gfibs = [], []
    for b in bundles:
        lp = stress_fiber_load_path(
            b["pos"], b["fiber_offsets"], b["polarity"], b["myo_i"], b["myo_j"],
            b["fa_sites"], b["linc_sites"], f_head_pN=f_head_pN, n_bins=n_bins)
        load_paths.append(lp)
        gfibs.append(b["global_fibers"])
    status = OK_STATUS if f_head_pN is not None else GAP_STATUS
    return WovenCellSFWiring(load_paths=load_paths, global_fibers=gfibs, f_head_pN=f_head_pN, magnitude_status=status)
