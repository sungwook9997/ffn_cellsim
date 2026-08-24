"""Population builders for the arena, and the one driver that stands them up together and measures.

WHAT LIVES HERE.  :mod:`~aleph.world.build.cortex`, :mod:`~aleph.world.build.membrane` and
:mod:`~aleph.world.build.envelope` each build ONE population into a :class:`~aleph.world.arena.WorldArena`
using Warp CUDA kernels.  This module holds the driver, because the number PHASE 1 exists to take —
peak GPU bytes — is a property of the three standing in one process at once and is not any single
builder's to report.

THE DRIVER IS WHERE THE AXES ARE DECLARED, AND THAT IS DELIBERATE.  Every builder refuses to default a
physiological value; something has to supply one, and that something should be the thing that also
writes the record.  :data:`NATIVE_AXES` is that declaration: each entry is a value, a unit, and a
provenance label copied from where the value came from.  A reader of the emitted record can see the
CONVENIENCE and PI-GAP labels without having to find this file.

⚠ **SCOPE.**  ``build_all`` stands up the three populations THIS package owns — cortex, plasma
membrane, nuclear envelope.  MT, IF, filopodium, lamellipodium and stress fibres are being built in
parallel and are not imported here; the whole-cell total is assembled above this layer, not inside it.
Any figure this module prints is therefore a figure for three populations and says so.

WHAT MAY NOT BE CLAIMED FROM WHAT THIS PRODUCES.  No physics of any kind: no force is computed, no law
is bound, nothing is relaxed and nothing steps.  Not a step throughput either — there is no step.  The
claim available is that the populations STAND, that the arena's claims partition its live prefix, and
what the standing cell costs in bytes and in wall time.

⚠ **AND THE PEAK-BYTES FIGURE IS TWO NUMBERS, NOT ONE.**  ``warp_mempool_high_bytes`` is Warp's
high-water mark over its own allocator, which is every byte this build allocates.
``device_used_delta_bytes`` is the change in device-wide used memory across the build, which also
carries the CUDA context and any other process on the card.  Neither is
``exact_peak_gpu_bytes`` in the sense the ledger uses that name: that one is NVML's per-process
``maxMemoryUsage``.  ⚠ **2026-08-22 — this module now READS that third number.**  The implementation
was ``components/incumbent/gpu_memory_accounting.py``, which ``world/`` may not import
(``test_layer_directions``), so this record carried ``exact_peak_gpu_bytes: None`` and a status
string calling the port an open item.  The port happened: it is :mod:`aleph.world.gpu_memory`.

⚠ **Being read is not being exact.**  NVML accounting mode must be enabled by the driver BEFORE the
process starts, which needs administrator authority this repository never takes, so on a host where
it is off the probe returns an explicit non-exact status and the field stays ``None``.  All three
numbers are reported under names that say what they are, and ``exact_peak_gpu_bytes_status`` now
carries what the probe actually said rather than a description of a missing import.

Sanity Gate (recorded per CLAUDE.md, before first execution):
  * dimensional — lengths µm, areal density µm⁻², memory bytes, time seconds.
  * boundary — a non-CUDA device raises through the arena rather than falling back; the arena capacity
    is computed from the same arithmetic the builders use, so an over-capacity claim is a bug in that
    arithmetic and is allowed to raise rather than being padded around.
  * conservation/invariant — ``assert_partitioned`` on every kind and ``assert_disjoint_spans`` on
    every node array, after the builds and before any number is reported.
  * CFL/precision — no integration; float64 throughout.
  * sign sense — the only signed quantity is each closed surface's enclosed volume, asserted positive
    by its own builder.
  * measurement protocol — wall time is taken around a ``wp.synchronize()`` so an asynchronous launch
    cannot be reported as a fast build; memory is sampled after that synchronise.

engine units: length µm.  Runtime: NVIDIA Warp on CUDA.
"""

from __future__ import annotations

import json
import math
import time
from dataclasses import dataclass

import warp as wp

from aleph.laws.cell_geometry import resolve_cell_geometry
from aleph.laws.cortex_assembly import CortexParams
from aleph.world.arena import Kind, WorldArena
from aleph.world.gpu_memory import query_process_peak
from aleph.world.build.cortex import (
    StrandPopulation,
    build_cortex,
    build_strand_population,
    filaments_for_density,
)
from aleph.world.build.envelope import ENVELOPE_MESH_PI_GAP, build_envelope, nuclear_radius_um
from aleph.world.build.membrane import (
    ClosedSurface,
    build_closed_surface,
    build_membrane,
    icosphere_counts,
    mean_edge_um,
    subdivisions_for_mesh,
)

__all__ = [
    "ClosedSurface",
    "ENVELOPE_MESH_PI_GAP",
    "NATIVE_AXES",
    "StandingCell",
    "StrandPopulation",
    "build_all",
    "build_closed_surface",
    "build_cortex",
    "build_envelope",
    "build_membrane",
    "build_strand_population",
    "filaments_for_density",
    "icosphere_counts",
    "mean_edge_um",
    "measure_native",
    "native_counts",
    "verify_standing_cell",
    "nuclear_radius_um",
    "subdivisions_for_mesh",
]


# ── the declared axes ───────────────────────────────────────────────────────────────────────────
#: Cortical / membrane-skeleton mesh band [µm], SOURCED: Morone 2006, Bovellan 2014,
#: Chugh & Paluch 2018.  Both ends are used, for different quantities, and which end is used where is
#: itself a decision — see the two entries below.
CORTEX_MESH_BAND_UM = (0.050, 0.100)


def _axes() -> dict[str, dict[str, object]]:
    """The declared physiological axes, each with its provenance label.

    Read from ``aleph.laws`` wherever a law already carries the value, so a swept or corrected law
    moves this driver with it instead of leaving a second copy behind.
    """
    geom = resolve_cell_geometry("mcf7")
    cortex_law = CortexParams()
    return {
        "R_cell_um": {
            "value": float(geom.R_cell_um), "unit": "um", "provenance": "SOURCED",
            "source": "laws/cell_geometry.MCF7_GEOMETRY — MCF7 cell volume 1760 um^3 (BNID 115154, "
                      "Gamcsik 1995 via Wagner 2011) -> sphere-equivalent R 7.49 um; confirmed by "
                      "photoacoustic sizing n=37.",
        },
        "h_cortex_um": {
            "value": 0.200, "unit": "um", "provenance": "SOURCED",
            "source": "h_cortex ~ 200 nm (KB-3.1 / KB-3.5; Salbreux 2012 TCB, Charras 2008 BJ, "
                      "Clark 2013 10.1016/j.bpj.2013.05.057), quoted in laws/cortex_assembly.py:61.",
        },
        "cortex_areal_density_um2": {
            "value": float(cortex_law.areal_density_um2), "unit": "um^-2", "provenance": "SOURCED",
            "source": "actin filament areal density ~100 um^-2 (KB-3.18); read back from "
                      "laws/cortex_assembly.CortexParams rather than retyped.",
        },
        "cortex_contour_um": {
            "value": float(cortex_law.L_filament_um), "unit": "um", "provenance": "CONVENIENCE",
            "source": "laws/cortex_assembly.CortexParams.L_filament_um. ⚠ PARAM_PROVENANCE_AUDIT "
                      "2026-07-24 rates 3.0 um CONVENIENCE (coarse): sourced formin ~1 um / Arp2/3 "
                      "~0.1 um means 3-30x TOO LONG, and the Arp2/3 short-filament population (~1/3 "
                      "of cortical actin by mass) is a MISSING ARCHITECTURE rather than a wrong "
                      "number. Carried forward unchanged so this phase changes one axis, not two.",
        },
        "cortex_seg_um": {
            "value": CORTEX_MESH_BAND_UM[0], "unit": "um", "provenance": "SOURCED (band floor)",
            "source": f"cortical mesh / crosslink spacing {CORTEX_MESH_BAND_UM[0] * 1e3:.0f}-"
                      f"{CORTEX_MESH_BAND_UM[1] * 1e3:.0f} nm (Morone 2006, Bovellan 2014, Chugh & "
                      "Paluch 2018); the FLOOR of that band, directed by the PI 2026-08-20. It "
                      "replaces the incumbent's 0.5 um, which the audit rates CONVENIENCE with 0 KB "
                      "rows and 5-10x TOO COARSE. ⚠ The band CEILING (0.100 um) is equally sourced "
                      "and halves the cortex node count; which end is used is a PI decision, not a "
                      "derivation, and it is recorded here as one.",
        },
        "membrane_mesh_um": {
            "value": CORTEX_MESH_BAND_UM[1], "unit": "um", "provenance": "SOURCED (band ceiling)",
            "source": "same band, read as an upper bound: the membrane mesh must be no COARSER than "
                      "100 nm, so the subdivision is the coarsest level lying inside the band. That "
                      "rule returns 7 at R_cell, which is the level WORLD_PORT_PLAN and PI 2026-07-22 "
                      "both name — reached from the spacing rather than typed.",
        },
        "envelope_mesh_um": {
            "value": CORTEX_MESH_BAND_UM[1], "unit": "um", "provenance": "PI-GAP",
            "source": ENVELOPE_MESH_PI_GAP,
        },
        "nc_ratio": {
            "value": float(geom.nc_ratio), "unit": "1", "provenance": "SOURCED",
            "source": "nucleus:cell RADIUS ratio 0.68 +/- 0.08 — imaging flow cytometry n=2164 "
                      "(PMC7000884) and photoacoustic n=37 (10.1007/s10765-016-2129-y), both 0.68. "
                      "⚠ CONFLICT: components/incumbent/compartments.py carries r_eq_um = 5.0 "
                      "labelled 'provisional geometry; PI GAP'. The sourced ratio gives 5.10 um. The "
                      "laws/ value is used; the conflict is reported, not closed.",
        },
        "seed": {
            "value": 0, "unit": "1", "provenance": "NUMERICAL",
            "source": "RNG seed for the cortical filament frames. Not physiological: it selects one "
                      "realisation of the quenched disorder, and a conclusion that depends on it is a "
                      "conclusion about one draw.",
        },
    }


#: The declared axes, evaluated once at import.  Each entry carries value, unit, provenance and source.
NATIVE_AXES: dict[str, dict[str, object]] = _axes()


def _v(name: str) -> float:
    """The value of a declared axis."""
    return float(NATIVE_AXES[name]["value"])  # type: ignore[arg-type]


# ── the arithmetic, before any device exists ────────────────────────────────────────────────────
def native_counts() -> dict[str, object]:
    """Every count the native build will claim, from :data:`NATIVE_AXES`, as pure arithmetic.

    Available without a card, which is what lets the arena be sized exactly rather than padded, and
    what lets the derivation be checked in a test that does not need a GPU.
    """
    r_cell = _v("R_cell_um")
    r_shell = r_cell - 0.5 * _v("h_cortex_um")
    n_fil = filaments_for_density(_v("cortex_areal_density_um2"), r_shell)
    n_seg_per = int(round(_v("cortex_contour_um") / _v("cortex_seg_um")))
    n_per = n_seg_per + 1

    mem_level = subdivisions_for_mesh(r_cell, _v("membrane_mesh_um"))
    mem_v, mem_e, mem_f = icosphere_counts(mem_level)
    r_nuc = nuclear_radius_um("mcf7")
    env_level = subdivisions_for_mesh(r_nuc, _v("envelope_mesh_um"))
    env_v, env_e, env_f = icosphere_counts(env_level)

    return {
        "cortex": {
            "shell_radius_um": r_shell, "n_filaments": n_fil, "nodes_per_filament": n_per,
            "node": n_fil * n_per, "segment": n_fil * n_seg_per, "angle3": n_fil * (n_per - 2),
        },
        "membrane": {
            "radius_um": r_cell, "subdivisions": mem_level,
            "mesh_um_realised": mean_edge_um(r_cell, mem_level),
            "node": mem_v, "face": mem_f, "angle4": mem_e,
        },
        "nuclear_envelope": {
            "radius_um": r_nuc, "subdivisions": env_level,
            "mesh_um_realised": mean_edge_um(r_nuc, env_level),
            "node": env_v, "face": env_f, "angle4": env_e,
        },
        "total": {
            "node": n_fil * n_per + mem_v + env_v,
            "segment": n_fil * n_seg_per,
            "angle3": n_fil * (n_per - 2),
            "face": mem_f + env_f,
            "angle4": mem_e + env_e,
        },
    }


@dataclass(frozen=True, slots=True)
class StandingCell:
    """The three populations this package owns, standing in one arena.

    Attributes:
        arena: the world they were claimed from.
        cortex / membrane / envelope: the built populations.
    """

    arena: WorldArena
    cortex: StrandPopulation
    membrane: ClosedSurface
    envelope: ClosedSurface

    def record(self, *, with_stats: bool = True) -> dict[str, object]:
        """Host-side census of the three populations plus the arena's own."""
        return {
            "populations": [
                self.cortex.record(with_stats=with_stats),
                self.membrane.record(),
                self.envelope.record(),
            ],
            "arena": self.arena.census(),
            "topology_bytes_total": (
                self.cortex.topology_bytes + self.membrane.topology_bytes + self.envelope.topology_bytes
            ),
        }


def cortex_shell() -> tuple[float, float]:
    """The cortical shell as ``(radius_um, thickness_um)`` — the ONE place it is derived.

    ⚠ **Added 2026-08-22 to remove a second derivation, not to add a helper.** `geometry.py` used to
    carry its own `nmii_radius_um = r_cell - tip_clearance - 0.30`, which was a different formula for
    the same shell and drifted from this one by 0.45 µm — enough that the NMII population stood clear
    of the cortex it exists to pull and **0 of 442 minifilaments could station** (`STATE.md` (c) 20).
    PI queue 14 named the cause exactly: *nothing owned the relationship between the two shells.*
    Now something does.

    ⚠ **And the old formula could not have been repaired by changing its constant.** `r_cell -
    tip_clearance` alone already reaches 7.25 with the cortex shell at 7.40, so no subtraction gets
    there: the SHAPE was wrong. The NMII shell is not the membrane pulled in by a clearance, it is the
    cortical shell — which is what `build/nmii.py`'s own self-check had been using all along.

    Returns:
        ``(radius_um, thickness_um)`` in µm, both derived — neither is typed here.
    """
    counts = native_counts()
    return (float(counts["cortex"]["shell_radius_um"]),  # type: ignore[index]
            float(_v("h_cortex_um")))


#: The ONE builder whose module name is not the population name it stands. ⚠ This is a REPORTING
#: accommodation, not a ruling: `test_families_census.py` carries an OPEN PI disagreement over whether
#: the endpoint is `nucleus` or `nuclear_envelope`, marked xfail(strict) so that ruling it fails loudly.
#: Mapping it here decides nothing — it stops :func:`builders_not_standing` from reporting a population
#: absent because two names for it disagree, which would bury the real absences under a naming argument.
_BUILDER_ALIAS = {"envelope": "nuclear_envelope"}


def builders_not_standing(populations: object) -> tuple[str, ...]:
    """Which builder modules under ``build/`` produced NO population in the cell that was stood.

    ⚠ **WHY THIS EXISTS.** The PHASE 4 cell the 2026-08-21 tau runs were measured on stands eleven
    populations, and `build/nmii.py` is not among the builders that made it — **there is no motor in
    that cell**, and `build/cytosol.py` never ran either, which is why `grid_cell` live is 0. Neither
    fact appears in the run record. A reader had to know the population list by heart to notice, and
    the record's own `placement_envelope` carried an `nmii_radius_um` while containing no NMII.
    `test_observe_gamma.py` had already written the finding down, in the docstring of a test defending
    a *different* thing — `gamma_source` must be `None` rather than `0.0`, because "the motors
    contributed nothing" and "there are no motors" are different findings. The distinction stood; what
    it was distinguishing never reached a record.

    ⚠ **DERIVED from the filesystem, never typed.** A builder added tomorrow is covered with no edit
    here. A typed list is the same defect one level up: it would report the absences someone remembered
    to list, which is exactly the failure this closes.

    ⚠ **It reports BUILDER modules, not populations.** A builder may stand several populations or, like
    `cytosol`, none that carry a name in the census. So a name here means *"this module contributed
    nothing to this cell"*, which is the question worth asking, and never *"this population is empty"*.

    Args:
        populations: any iterable of the population names the cell actually stood.

    Returns:
        Builder module stems, sorted, that match no standing population. Empty when the cell is whole.
    """
    from pathlib import Path as _Path

    stood = {str(p) for p in populations}  # type: ignore[union-attr]
    here = _Path(__file__).resolve().parent
    absent = []
    for f in sorted(here.glob("*.py")):
        if f.stem.startswith("_"):
            continue
        if _BUILDER_ALIAS.get(f.stem, f.stem) not in stood:
            absent.append(f.stem)
    return tuple(absent)


def build_all(device: str = "cuda:0", *, arena: WorldArena | None = None,
              seed: int | None = None) -> StandingCell:
    """Stand cortex, plasma membrane and nuclear envelope in one arena, on the device.

    Args:
        device: the CUDA device string.  Ignored when ``arena`` is given.
        arena: an existing arena to claim from.  When ``None`` one is allocated with exactly the
            capacity :func:`native_counts` says is needed — exactly, so that an over-capacity claim
            surfaces an error in the arithmetic instead of being absorbed by padding.

    Returns:
        The :class:`StandingCell`.

    Raises:
        RuntimeError: if ``device`` is not CUDA — raised by the arena, which has no CPU path.
    """
    counts = native_counts()
    if arena is None:
        total = counts["total"]  # type: ignore[index]
        arena = WorldArena(
            capacity={
                Kind.NODE: total["node"], Kind.SEGMENT: total["segment"],
                Kind.ANGLE3: total["angle3"], Kind.FACE: total["face"], Kind.ANGLE4: total["angle4"],
            },
            device=device,
        )

    cx = build_cortex(
        arena,
        radius_um=counts["cortex"]["shell_radius_um"],  # type: ignore[index]
        thickness_um=_v("h_cortex_um"),
        areal_density_um2=_v("cortex_areal_density_um2"),
        contour_um=_v("cortex_contour_um"),
        seg_um=_v("cortex_seg_um"),
        density_provenance=str(NATIVE_AXES["cortex_areal_density_um2"]["source"]),
        # ⚠ `seed` overrides the declared axis, and it exists so INDEPENDENT REPLICATES are
        # possible. Added 2026-08-21: `STATE.md` (f) carries a standing rule that gamma is judged
        # against SEED SCATTER rather than within-run `sem` — three replicates give a spread 67x the
        # sem — and (e) 1's D-2 calls a criterion built on within-run variance "the single most likely
        # way to write the amendment and still be measuring nothing". Without this argument, running
        # the same driver three times at three `--seed` values produced the SAME CELL, and offering
        # those as a scatter would have been that mistake with extra steps.
        seed=int(_v("seed") if seed is None else seed),
    )
    mem = build_membrane(
        arena,
        radius_um=_v("R_cell_um"),
        mesh_um=_v("membrane_mesh_um"),
        mesh_provenance=str(NATIVE_AXES["membrane_mesh_um"]["source"]),
    )
    env = build_envelope(
        arena,
        radius_um=nuclear_radius_um("mcf7"),
        mesh_um=_v("envelope_mesh_um"),
        mesh_provenance=ENVELOPE_MESH_PI_GAP,
    )
    return StandingCell(arena=arena, cortex=cx, membrane=mem, envelope=env)


# ── the measurement ─────────────────────────────────────────────────────────────────────────────
def _device_used_bytes(device: wp.Device) -> int | None:
    """Device-wide used memory [bytes], or ``None`` if the driver does not report it.

    Device-WIDE: it includes the CUDA context and anything else on the card, which is exactly why the
    figure derived from it is reported under a name that says ``device_used`` rather than ``peak``.
    """
    total = getattr(device, "total_memory", None)
    free = getattr(device, "free_memory", None)
    if total is None or free is None:
        return None
    return int(total) - int(free)


def measure_native(device: str = "cuda:0", out: str | None = None) -> dict[str, object]:
    """Stand the three populations natively and measure what standing them costs.

    Measures build wall time per population (around a synchronise, so an asynchronous launch cannot
    be reported as a fast build), Warp's allocator high-water mark, the device-wide used-memory delta,
    and bytes per node.  Asserts the arena's partition and byte-span invariants before reporting.

    **Claims nothing about physics and nothing about a step.**  Nothing steps.

    Args:
        device: the CUDA device string.
        out: path to write the record to as JSON.  ``None`` writes nothing and only returns it.

    Returns:
        The record.

    Raises:
        RuntimeError: if ``device`` is not CUDA.
        AssertionError: on any arena invariant.
    """
    wp.init()
    dev = wp.get_device(device)
    if not dev.is_cuda:
        raise RuntimeError(
            f"measure_native resolved {device!r} to {dev!r}, which is not CUDA. Warp CUDA is the only "
            "simulation runtime and a number measured on CPU is not a weak result; it is not a result."
        )

    counts = native_counts()
    used_before = _device_used_bytes(dev)
    mempool = bool(wp.is_mempool_enabled(dev)) if wp.is_mempool_supported(dev) else False

    t0 = time.perf_counter()
    total = counts["total"]  # type: ignore[index]
    arena = WorldArena(
        capacity={
            Kind.NODE: total["node"], Kind.SEGMENT: total["segment"], Kind.ANGLE3: total["angle3"],
            Kind.FACE: total["face"], Kind.ANGLE4: total["angle4"],
        },
        device=device,
    )
    wp.synchronize_device(dev)
    t_arena = time.perf_counter() - t0

    timings: dict[str, float] = {"arena_alloc_s": t_arena}
    t = time.perf_counter()
    cx = build_cortex(
        arena, radius_um=counts["cortex"]["shell_radius_um"],  # type: ignore[index]
        thickness_um=_v("h_cortex_um"), areal_density_um2=_v("cortex_areal_density_um2"),
        contour_um=_v("cortex_contour_um"), seg_um=_v("cortex_seg_um"),
        density_provenance=str(NATIVE_AXES["cortex_areal_density_um2"]["source"]),
        seed=int(_v("seed")),
    )
    wp.synchronize_device(dev)
    timings["cortex_s"] = time.perf_counter() - t

    t = time.perf_counter()
    mem = build_membrane(arena, radius_um=_v("R_cell_um"), mesh_um=_v("membrane_mesh_um"),
                         mesh_provenance=str(NATIVE_AXES["membrane_mesh_um"]["source"]))
    wp.synchronize_device(dev)
    timings["membrane_s"] = time.perf_counter() - t

    t = time.perf_counter()
    env = build_envelope(arena, radius_um=nuclear_radius_um("mcf7"), mesh_um=_v("envelope_mesh_um"),
                         mesh_provenance=ENVELOPE_MESH_PI_GAP)
    wp.synchronize_device(dev)
    timings["envelope_s"] = time.perf_counter() - t
    timings["total_s"] = time.perf_counter() - t0

    # The invariants, before any number is reported.
    arena.assert_partitioned()
    for name in sorted(arena.node_arrays):
        arena.assert_disjoint_spans(name)

    used_after = _device_used_bytes(dev)
    high = int(wp.get_mempool_used_mem_high(dev)) if mempool else None
    standing = StandingCell(arena=arena, cortex=cx, membrane=mem, envelope=env)

    n_nodes = arena.n_live(Kind.NODE)
    census = arena.census()
    node_bytes = int(census["node_bytes_arithmetic"])  # type: ignore[arg-type]
    topo_bytes = int(standing.record(with_stats=False)["topology_bytes_total"])  # type: ignore[arg-type]

    verification = verify_standing_cell(standing)

    # Read AFTER the build, per the probe's own Sanity Gate: it queries a driver statistic and
    # introduces no simulation-state D2H. A probe that cannot answer must say so in the record
    # rather than take the run down with it — the same reporting rule the drivers use.
    try:
        _nvml = query_process_peak(device)
        nvml_peak = _nvml.exact_peak_bytes if _nvml.exact else None
        nvml_status = _nvml.probe_status
    except Exception as exc:  # noqa: BLE001 — an unavailable probe is reported, never swallowed
        nvml_peak = None
        nvml_status = f"PROBE_RAISED: {type(exc).__name__}: {exc}"

    record: dict[str, object] = {
        "kind": "diagnostic",
        "phase": "WORLD_PORT PHASE 1a — geometry only",
        "scope": "cortex + plasma membrane + nuclear envelope ONLY. MT, IF, filopodium, "
                 "lamellipodium and stress fibres are built elsewhere and are NOT in these totals.",
        "may_not_claim": [
            "any physics — no force is computed and no law is bound",
            "any step throughput — nothing steps",
            "that the geometry is correct because it built",
        ],
        "device": {
            "requested": device, "resolved": str(dev), "name": getattr(dev, "name", None),
            "arch": getattr(dev, "arch", None), "total_memory_bytes": getattr(dev, "total_memory", None),
            "mempool_enabled": mempool,
        },
        "warp_version": wp.config.version,
        "axes": NATIVE_AXES,
        "counts": counts,
        "populations": standing.record(),
        "timings_s": timings,
        "memory": {
            "warp_mempool_high_bytes": high,
            "warp_mempool_high_note": "Warp allocator high-water mark. Every byte this build "
                                      "allocates goes through it; the CUDA context does not.",
            "device_used_before_bytes": used_before,
            "device_used_after_bytes": used_after,
            "device_used_delta_bytes": (None if used_before is None or used_after is None
                                        else used_after - used_before),
            "device_used_note": "device-WIDE, so it carries the CUDA context and any other process "
                                "on the card. Not a per-process figure.",
            "arena_node_bytes_arithmetic": node_bytes,
            "population_topology_bytes": topo_bytes,
            "bytes_per_node_arena": node_bytes / n_nodes if n_nodes else None,
            "bytes_per_node_with_topology": (node_bytes + topo_bytes) / n_nodes if n_nodes else None,
            "bytes_per_node_measured": (None if high is None or not n_nodes else high / n_nodes),
            "exact_peak_gpu_bytes": nvml_peak,
            "exact_peak_gpu_bytes_status": nvml_status,
            "exact_peak_gpu_bytes_note": (
                "NVML per-process maxMemoryUsage via aleph.world.gpu_memory, read AFTER the build. "
                "It is null unless the driver had accounting mode enabled before this process "
                "started, which needs administrator authority this repository never takes — so a "
                "null here is a HOST configuration fact, not a measurement that failed. The two "
                "figures above are what this module can measure without it."
            ),
        },
        "invariants": {
            "assert_partitioned": "PASS",
            "assert_disjoint_spans": sorted(arena.node_arrays),
            "verification": verification,
            "verification_note": "read back from what the kernels wrote; the topology entry is an "
                                 "EXACT comparison against strand.build_strand, which is what makes "
                                 "the claim-model change bookkeeping rather than geometry.",
        },
    }

    if out is not None:
        with open(out, "w", encoding="utf-8") as fh:
            json.dump(record, fh, indent=2, ensure_ascii=False)

    return record


# ── verification ────────────────────────────────────────────────────────────────────────────────
def verify_standing_cell(cell: StandingCell) -> dict[str, object]:
    """Check on the DEVICE'S OWN OUTPUT that the kernels built what they claim to have built.

    Every check above this point runs without a card and therefore checks arithmetic, not kernels.
    This one reads back what the kernels wrote and tests it.  Four things can go wrong in a build like
    this and none of them raises on its own: an off-by-one in a stride writes a filament's nodes over
    its neighbour's; a wrong arc formula leaves nodes off the shell; a wrong chord leaves the
    population pre-strained at t0; and a mis-globalised index makes a face address another population.
    All four are silent, and all four are visible here.

    ⚠ **This reads the full node array back to the host** — at native that is ~100 MB.  It is a
    build-time transfer, taken once, between accepted steps by construction because nothing has
    stepped; it is not free and it is not in any loop.  Call it once after a build, never per step.

    The topology check is EXACT rather than approximate, and that is the point of it: because claims
    are contiguous and per-kind, laying ``n`` filaments through ``strand.build_strand`` produces the
    SAME global node, segment and angle ids as one population claim of the same total — so the two
    claim models must agree bit for bit on every index.  If they do, the change this package makes is
    bookkeeping and not geometry, which is exactly the claim being made.

    Args:
        cell: a :class:`StandingCell` that has been built on a device.

    Returns:
        A dict of the measured deviations, for the record.  Every entry is a number that SHOULD be
        zero or bounded, so a reader can see the margin rather than a bare PASS.

    Raises:
        AssertionError: on any violated invariant, naming the deviation and the bound.
    """
    import numpy as np

    from aleph.world.strand import build_strand

    out: dict[str, object] = {}
    cx = cell.cortex
    n_fil, n_per = cx.n_strands, cx.nodes_per_strand
    n_seg_per = cx.segments_per_strand
    centre = np.asarray(cx.centre_um, np.float64)

    pos = cell.arena.download_nodes(cx.nodes, "position").reshape(n_fil, n_per, 3)
    assert np.isfinite(pos).all(), "the cortex kernel wrote a non-finite position"

    # 1. Every node of a filament sits on THAT filament's own radius. A wrong arc formula, a wrong
    #    stride, or a frame that is not orthonormal all break this and nothing else would notice.
    r = np.linalg.norm(pos - centre, axis=2)
    r_fil = r[:, 0:1]
    radius_dev = float(np.max(np.abs(r - r_fil) / r_fil))
    assert radius_dev < 1e-12, f"nodes drift off their own shell radius by {radius_dev:.3e} (relative)"
    out["max_relative_radius_drift_within_filament"] = radius_dev

    # 2. Those radii fill the shell thickness and do not leave it. A dispersion bug shows up as a
    #    span of zero (every filament on one sphere) or as a span past the layer.
    lo, hi = cx.radius_um - 0.5 * cx.thickness_um, cx.radius_um + 0.5 * cx.thickness_um
    r_min, r_max = float(r.min()), float(r.max())
    assert lo - 1e-12 <= r_min <= r_max <= hi + 1e-12, \
        f"filament radii span [{r_min}, {r_max}], outside the shell [{lo}, {hi}]"
    out["shell_radius_span_um"] = [r_min, r_max]
    out["shell_fill_fraction"] = (r_max - r_min) / cx.thickness_um if cx.thickness_um else None

    # 3. The rest length is the built chord, and the built chord is what the arc geometry predicts.
    #    This is the check that the population starts UNSTRAINED: seg_rest is measured from positions,
    #    so agreeing with the independent analytic chord means both are right rather than jointly wrong.
    rest = cx.seg_rest_um.numpy().reshape(n_fil, n_seg_per)
    predicted = 2.0 * r_fil * np.sin(cx.seg_um_realised / (2.0 * r_fil))
    chord_dev = float(np.max(np.abs(rest - predicted) / predicted))
    assert chord_dev < 1e-9, f"built chord differs from the arc prediction by {chord_dev:.3e}"
    out["max_relative_chord_deviation"] = chord_dev
    out["chord_vs_requested_step"] = float(np.mean(rest) / cx.seg_um_requested)

    # 4. The topology, against the per-strand builder it replaces — exact, for the reason in the
    #    docstring. A handful of filaments is enough: the arithmetic is uniform, so an error in it is
    #    an error at every stride, and doing this at native would be a slower way to learn the same fact.
    probe_fil = min(8, n_fil)
    probe = WorldArena(capacity={
        Kind.NODE: probe_fil * n_per, Kind.SEGMENT: probe_fil * n_seg_per,
        Kind.ANGLE3: probe_fil * (n_per - 2),
    })
    ref_seg, ref_ang = [], []
    for _ in range(probe_fil):
        st = build_strand(probe, "cortex", start=(0.0, 0.0, 0.0), direction=(1.0, 0.0, 0.0),
                          contour_um=cx.contour_um_requested, seg_um=cx.seg_um_requested)
        ref_seg.append(st.seg_node)
        ref_ang.append(st.angle_idx)
    got_seg = cx.seg_node.numpy()[: probe_fil * n_seg_per] - cx.nodes.lo
    got_ang = cx.angle_idx.numpy()[: probe_fil * (n_per - 2)] - cx.nodes.lo
    assert np.array_equal(got_seg, np.concatenate(ref_seg)), \
        "segment topology differs from strand.build_strand — the claim model changed the GEOMETRY"
    assert np.array_equal(got_ang, np.concatenate(ref_ang)), \
        "angle topology differs from strand.build_strand — the claim model changed the GEOMETRY"
    assert len(probe.claims(population="cortex")) == 3 * probe_fil
    assert len(cell.arena.claims(population=cx.population)) == 3
    out["topology_identical_to_per_strand_builder"] = True
    out["claims_per_strand_builder"] = 3 * probe_fil
    out["claims_this_builder"] = 3

    # 5. Nothing addresses outside its own claim. A shared node is a weld, not a connection, and it is
    #    the one arrangement the architecture forbids outright.
    for name, arr in (("seg_node", cx.seg_node), ("angle_idx", cx.angle_idx)):
        idx = arr.numpy()
        assert int(idx.min()) >= cx.nodes.lo and int(idx.max()) < cx.nodes.hi, \
            f"cortex {name} addresses outside [{cx.nodes.lo}, {cx.nodes.hi})"
    sid = cell.arena.download_nodes(cx.nodes, "strand_id").reshape(n_fil, n_per)
    assert np.array_equal(sid, np.repeat(np.arange(n_fil, dtype=sid.dtype)[:, None], n_per, axis=1)), \
        "strand_id does not label each node with its own filament"

    # 6. The two closed surfaces: on their radius, wound outward, closed, and inside their claims.
    for surf in (cell.membrane, cell.envelope):
        v = cell.arena.download_nodes(surf.nodes, "position")
        assert np.isfinite(v).all(), f"{surf.population} kernel wrote a non-finite vertex"
        rad = np.linalg.norm(v - np.asarray(surf.centre_um), axis=1)
        dev = float(np.max(np.abs(rad - surf.radius_um)) / surf.radius_um)
        assert dev < 1e-12, f"{surf.population} vertices are off the sphere by {dev:.3e}"
        faces = surf.face_idx.numpy()
        hinges = surf.hinge_idx.numpy()
        for name, idx in (("face_idx", faces), ("hinge_idx", hinges)):
            assert int(idx.min()) >= surf.nodes.lo and int(idx.max()) < surf.nodes.hi, \
                f"{surf.population} {name} addresses outside its node claim"
        assert hinges.shape[0] == 3 * faces.shape[0] // 2, "a closed mesh has 3F/2 hinges, not F"
        assert (hinges[:, 0] < hinges[:, 3]).all(), \
            "hinge opposites are not ordered — the build depends on kernel scheduling"
        # Every vertex is used, and used by the six-ish faces an icosphere gives it. An unreferenced
        # vertex means the subdivision dropped one; the twelve original poles keep five.
        used = np.bincount(faces.reshape(-1) - surf.nodes.lo, minlength=surf.nodes.count)
        top_valence = 6 if surf.subdivisions >= 1 else 5  # the bare icosahedron is 5 everywhere
        assert used.min() == 5 and used.max() == top_valence, (
            f"{surf.population} vertex valence spans [{used.min()}, {used.max()}], not the "
            f"5/{top_valence} an icosphere has — the subdivision lost or duplicated a vertex"
        )
        assert int((used == 5).sum()) == 12, "an icosphere has exactly twelve five-valent vertices"
        # The discrete volume is BELOW the analytic sphere and approaches it with refinement. Being
        # above it would mean the winding or the volume sum is wrong, not that the mesh is coarse.
        analytic = 4.0 / 3.0 * math.pi * surf.radius_um ** 3
        gap = (analytic - surf.volume0_um3) / analytic
        assert 0.0 < gap < 1e-3, f"{surf.population} discrete/analytic volume gap {gap:.3%}"
        out[f"{surf.population}_max_relative_radius_deviation"] = dev
        out[f"{surf.population}_discrete_analytic_volume_gap"] = gap
        out[f"{surf.population}_area0_total_um2"] = surf.area0_total_um2

    cell.arena.assert_partitioned()
    for name in sorted(cell.arena.node_arrays):
        cell.arena.assert_disjoint_spans(name)
    out["assert_partitioned"] = "PASS"
    return out


def _demo() -> None:
    """Self-check: the axes carry provenance and the counts derive, both without a card."""
    for name, axis in NATIVE_AXES.items():
        assert axis["provenance"], f"{name} has no provenance label"
        assert axis["source"], f"{name} has no source"
        assert isinstance(axis["value"], (int, float)), f"{name} has a non-numeric value"
    assert NATIVE_AXES["envelope_mesh_um"]["provenance"] == "PI-GAP", "the gap must stay labelled"
    assert NATIVE_AXES["cortex_contour_um"]["provenance"] == "CONVENIENCE"

    c = native_counts()
    # Every count is derived; none is typed. The relations are what a reader can check.
    assert c["cortex"]["shell_radius_um"] == _v("R_cell_um") - 0.5 * _v("h_cortex_um")
    assert c["cortex"]["n_filaments"] == filaments_for_density(
        _v("cortex_areal_density_um2"), c["cortex"]["shell_radius_um"])
    assert c["cortex"]["nodes_per_filament"] == int(round(
        _v("cortex_contour_um") / _v("cortex_seg_um"))) + 1
    assert c["cortex"]["node"] == c["cortex"]["n_filaments"] * c["cortex"]["nodes_per_filament"]
    assert c["cortex"]["segment"] == c["cortex"]["node"] - c["cortex"]["n_filaments"]
    assert c["cortex"]["angle3"] == c["cortex"]["node"] - 2 * c["cortex"]["n_filaments"]
    for surf in ("membrane", "nuclear_envelope"):
        v, e, f = icosphere_counts(c[surf]["subdivisions"])
        assert (c[surf]["node"], c[surf]["angle4"], c[surf]["face"]) == (v, e, f)
        assert c[surf]["angle4"] == 3 * c[surf]["face"] // 2, "a closed mesh has 3F/2 hinges"
        assert c[surf]["mesh_um_realised"] <= _v(f"{'membrane' if surf == 'membrane' else 'envelope'}_mesh_um")
    assert c["total"]["node"] == c["cortex"]["node"] + c["membrane"]["node"] + c["nuclear_envelope"]["node"]
    assert math.isclose(c["nuclear_envelope"]["radius_um"], _v("nc_ratio") * _v("R_cell_um"))

    # The dominant population is the cortex, and by how much is the fact the phase turns on.
    assert c["cortex"]["node"] / c["total"]["node"] > 0.9

    print(
        "build self-check OK — "
        f"cortex {c['cortex']['node']:,} + membrane {c['membrane']['node']:,} (subdiv "
        f"{c['membrane']['subdivisions']}) + envelope {c['nuclear_envelope']['node']:,} (subdiv "
        f"{c['nuclear_envelope']['subdivisions']}, PI-GAP) = {c['total']['node']:,} nodes, "
        "3 populations, 9 claims"
    )


if __name__ == "__main__":
    _demo()
