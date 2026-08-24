"""The cytosol: a FIELD, not a node population — and the arena has no primitive for it yet.

WHAT THE CYTOSOL IS, AND WHY IT IS NOT NODES.  Three independent pieces of the frozen port source
say the same thing, and none of them is a preference:

1. ``components/fluid/field_grid.py:9-11`` — *"the Eulerian Biot/RAD PDE fields live on this
   structured, conservative finite-volume grid"* — and the next line separates it from the Warp
   ``HashGrid`` used for filament neighbour search: *"a neighbour structure, NOT the PDE field grid.
   They may share origin/cell size but are distinct allocations."*  The distinction a node population
   would collapse is one the port source drew on purpose.
2. ``engine/fluid_core.py:215-243`` — ``FluidVolumeStateOwner`` validates ``pressure_d`` and
   ``mask_d`` at ``ndim=3`` and refuses anything else.  The Biot update walks a 6-point masked
   Laplacian over ``p[i, j, k]`` (``biot_substrate.py:80``), so the stencil IS the data structure.
3. ``engine/immersed_transfer.py:141-152`` — the six ``*_cytosol_transfer`` edges scatter
   ``-alpha * V_node * grad(p)`` onto the SOLID's own force array through a Peskin-4 interpolation.
   The cytosol end of that edge is a STENCIL AROUND A POSITION, not an ID.  There is no cytosol node
   to bond to, and minting one would replace an interpolation with a lumped stand-in.

``Kind.GRID_CELL`` — ADDED 2026-08-20, PI-approved.  ``arena.Kind`` had seven members and called the
set *"closed by argument rather than by convenience"*, then rejected a species pool with an argument
that named an eighth by name: *"a species pool is a named channel on a ``GRID_CELL`` field plus one
conserved scalar per consuming range."*  The argument that closed the set already presupposed a
member that had never been added, and the cytosol is what forced it.  :class:`ArenaKindGapError`
stays, and so does :func:`grid_cell_kind`: a rebase onto a tree without the member must still produce
a self-describing refusal rather than an ``AttributeError``.

The arena allocates device arrays **only** for ``Kind.NODE``
(:meth:`~aleph.world.arena.WorldArena.__post_init__`).  ``SEGMENT``, ``ANGLE3``, ``FACE`` and now
``GRID_CELL`` are pure bookkeeping — capacity, claim, partition invariant — and the kind-specific
arrays belong to the builder: ``build/cortex.py`` owns ``seg_node``/``seg_rest_um``,
``build/membrane.py`` owns ``face_idx``/``hinge_idx``, and this module owns the fields.  The 1-D claim
over a 3-D array is not new either: a claim counts ``S`` segments while the array is ``(S, 2)``; here
it counts ``nx*ny*nz`` cells while the arrays are ``(nx, ny, nz)``.

⚠ **AND ``dx`` IS NOT A DENSITY — IT ENTERS COST AS dx^-5.**  A FLUID_VOLUME has no per-cell
structure count to source.  It has a RESOLUTION, and that is a discretisation choice bounded on both
sides by sourced physics (:data:`CYTOSOL_AXES`).  Cells scale as dx^-3 and the explicit consolidation
CFL ``safety*S*dx^2/(2*d*mobility)`` (``field_grid.py:125-127``) makes subcycles per outer step scale
as dx^-2.  :func:`cytosol_counts` computes both, without a card, so the decision is arithmetic rather
than taste.  At the cortex's own sourced 0.05 µm the grid is ~28 M cells against ~4.9 M solid nodes.

⚠ **PI 2026-08-20 set** :data:`PI_DX_UM` **= 0.25 µm, and it is a TEST POINT, not a derivation.**  It
is inside the band, one step finer than the incumbent's 0.5 (which ``assemble.py:333`` labels
*"numerical sizing; not literature"*), and it costs 28x that — while deliberately leaving 0.1
(2,480x) as a question the built cell can ANSWER by measurement once ``k_xl`` lands.  Choosing 0.05
now would have spent that measurement.  **The value is therefore still passed as an argument and
never defaulted**: a build must state which resolution it ran at, the band is still enforced at both
ends, and the provenance travels into the record.  A ``k_xl`` decision that lands in the field band
reopens this one.

⚠ **TWO DIFFERENT QUANTITIES ARE BOTH CALLED ``mobility`` AND THEY MUST NOT BE CROSSED.**  Here it
is the DARCY mobility ``k / mu_pore`` [µm²/(Pa·s)] — a permeability over the INTERSTITIAL PORE
viscosity, which is 2-8x water (KB-3.31, whose own text says *"Wire pore viscosity into the fluid
channel … poroelastic, not one lumped viscosity"*).  The arena's per-node ``mobility``
[µm/(pN·s)] is something else entirely: the parameter of a DRAG law, which ``build/membrane.py:49-50``
leaves at zero on purpose because *"a builder that filled it in would be inventing"*.  A bulk
cytoplasm viscosity such as ``laws/units.ETA_CYTOPLASM`` belongs to neither of them, and session F
found it standing in for a nucleoplasm value on the drag side (``laws/motility_warp.py:429``).
Nothing here reads it, and nothing here may.

⚠ **KB-3.21's ``V_cyto ~ 4200 um^3`` IS NOT THIS CELL'S AND IS NOT QUOTED HERE.**  The row is
``verified``/``High`` and it is not wrong — 4,200 µm³ is a sphere of radius 10.0 µm, and the arena
builds MCF7 at ``R_cell`` 7.50 µm, whose ``V_cell - V_nuc`` is 1,211 µm³, a factor 3.47 away.  A
value that is perfectly sourced and scoped to a different cell is ``bond.SourceClass.UNRATIFIED_PROXY``,
not a source.  So the enclosed-volume check here is self-contained: ``sum_FLUID dV`` measured from the
mask the build produces, against the enclosed volumes the membrane and envelope builders already
measured from what THEY built (``ClosedSurface.volume0_um3``), never against an external row.

Sanity Gate (recorded per CLAUDE.md, before first execution):
  * dimensional — ``dx_um``, ``margin_um``, radii and the box are lengths [µm]; ``D_p`` is [µm²/s];
    ``dt_phys_s`` is [s]; ``mask`` is a class code and carries no unit; every field is allocated and
    left at zero, because PHASE 1 computes no physics.
  * boundary — ``dx_um`` has NO default and the build refuses without it; a ``dx`` outside the sourced
    band is refused with both ends named; a grid that does not strictly contain the membrane is
    refused, because ``live_mesh_domain`` requires the surface to stay inside the box for the winding
    query bound to hold; a non-CUDA device raises through the arena.
  * conservation/invariant — one claim of ``nx*ny*nz`` cells, so ``assert_partitioned`` holds by
    construction; ``sum_FLUID dV`` is reported against ``volume0_um3(membrane) -
    volume0_um3(envelope)`` and the discretisation error is reported, never absorbed.
    ``assert_disjoint_spans`` does not apply: like the cortex's ``seg_*`` arrays these are separate
    allocations, not views into an arena array.
  * CFL/precision — nothing integrates here.  Fields are float64 to match the port source.  The
    winding-number query is float32 (``wp.Mesh`` is a float32 structure), which is why the
    surface-tie epsilon is carried explicitly rather than assumed away.
  * sign sense — a watertight outward-wound surface gives ``query.sign < 0`` INSIDE; the envelope is
    applied after the membrane and overwrites FLUID with NUCLEUS.
  * measurement protocol — every count in :func:`cytosol_counts` is host arithmetic over declared
    axes and needs no device; the mask reduction is a device reduction read once, between builds.

engine units: length µm, time s.  Runtime: NVIDIA Warp on CUDA only.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from inspect import Parameter, signature

import warp as wp

from aleph.world.arena import Claim, Kind, WorldArena

__all__ = [
    "CYTOSOL_AXES",
    "PI_DX_UM",
    "GRID_CELL_DECISION",
    "ArenaKindGapError",
    "CytosolField",
    "build_cytosol",
    "cytosol_counts",
    "grid_cell_kind",
]

# Mask codes, identical to ``components/fluid/fv_reference`` so the arena's classification and the
# port source's discretisation cannot silently disagree about what a cell class means.
OUTSIDE = 0
FLUID = 1
NUCLEUS = 2

_OUTSIDE = wp.constant(0)
_FLUID = wp.constant(1)
_NUCLEUS = wp.constant(2)

#: float64 fields this population owns, per grid cell.  ``mask`` is int32; every other is 8 B.  Named
#: rather than counted so a reader can see WHICH arrays the byte figure is for.
_FIELD_NAMES = ("p", "p_new", "p_bar", "s_water", "s_membrane", "s_total", "div_vs")
_BYTES_PER_CELL = 8 * len(_FIELD_NAMES) + 4

GRID_CELL_DECISION = (
    "arena.Kind has no GRID_CELL member. The cytosol is a FIELD — the port source is Eulerian "
    "(field_grid.py:9-11), FluidVolumeStateOwner validates ndim=3 (fluid_core.py:215-243), and the "
    "six *_cytosol_transfer edges couple through a Peskin-4 stencil, not a node id "
    "(immersed_transfer.py:141-152). arena.py:77-78 already NAMES `GRID_CELL` inside the argument "
    "that closes the seven-member set, so an eighth was assumed and never added. Adding it costs one "
    "enum member: the arena allocates arrays only for Kind.NODE (arena.py:170-181) and every other "
    "kind is bookkeeping whose arrays the builder owns. Lead owns arena.py and the decision; this "
    "module refuses until it lands. See docs/v2_audit/CYTOSOL_ARENA_REPRESENTATION_2026-08-20.md."
)


class ArenaKindGapError(NotImplementedError):
    """A population cannot be claimed because the arena has no primitive for what it IS.

    The sibling of :class:`aleph.world.families.ConnectorGapError`, and it exists for the same reason:
    a refusal that carries its question is a product, while a partially-built object that looks built
    is a defect.  This one is about the arena's vocabulary rather than a missing count — no value
    would close it, only a decision.

    Args:
        population: the population that cannot be claimed.
        kind_name: the primitive that does not exist.
        why: the argument for adding it, and where the decision is recorded.
    """

    def __init__(self, population: str, kind_name: str, why: str) -> None:
        self.population = population
        self.kind_name = kind_name
        self.why = why
        super().__init__(
            f"{population}: arena.Kind has no {kind_name!r}, so this population cannot be claimed. "
            f"This is a PENDING DECISION, not an incomplete file. {why}"
        )


def grid_cell_kind() -> Kind:
    """The arena's grid primitive, or a refusal that names the decision it waits on.

    Returns:
        ``Kind.GRID_CELL`` once it exists.

    Raises:
        ArenaKindGapError: while it does not.  Deliberately not ``AttributeError``: a caller and a
            reader must both be able to tell a pending decision from a broken import.
    """
    kind = getattr(Kind, "GRID_CELL", None)
    if kind is None:
        raise ArenaKindGapError("cytosol", "GRID_CELL", GRID_CELL_DECISION)
    return kind


# ── the declared axes ───────────────────────────────────────────────────────────────────────────
#: Continuum FLOOR for a pore-pressure field [µm]: below the cytoskeletal pore/mesh size a Darcy
#: continuum pressure is not defined, and the intrinsic permeability k ~ xi^2 is set by it.
PORE_MESH_BAND_UM = (0.014, 0.040)

#: The resolution PI set 2026-08-20 [µm].  **A TEST POINT INSIDE THE BAND, NOT A DERIVATION**, and it
#: is deliberately NOT a default: :func:`build_cytosol` still requires ``dx_um`` as an argument so a
#: record always says which resolution it ran at.  Chosen one step finer than the incumbent's 0.5 µm
#: (28x its fluid work) while leaving 0.1 µm — 2,480x, and the next question — answerable BY
#: MEASUREMENT on the built cell once ``k_xl`` lands.  Picking 0.05 now would have spent that.
PI_DX_UM = 0.25

CYTOSOL_AXES: dict[str, dict[str, object]] = {
    "dx_um": {
        "value": PI_DX_UM, "unit": "um", "provenance": "PI-DECISION 2026-08-20 (test point)",
        "source": "PI set 0.25 um in session, transcript d380041d-06a6-49ea-9daf-1e4e57e76bf0. NOT "
                  "DERIVED and not a threshold: a FLUID_VOLUME has no per-cell structure count — it "
                  "has a RESOLUTION, and dx is a discretisation choice, not a physiological density. "
                  "It is still passed as an argument and never defaulted, so a record states the "
                  "resolution it ran at rather than inheriting one. Reopened if k_xl lands in the "
                  "field band. The "
                  f"band is sourced at both ends: FLOOR = pore/mesh size {PORE_MESH_BAND_UM[0] * 1e3:.0f}"
                  f"-{PORE_MESH_BAND_UM[1] * 1e3:.0f} nm (KB-DRAFT-3.B-27, KB-DRAFT-3.B-18, both "
                  "DRAFT), below which a continuum pore pressure is undefined; CEILING = the "
                  "poroelastic diffusion length sqrt(D_p * dt_phys), which must be resolved rather "
                  "than stepped over. The incumbent's 0.5 um sits inside the band and was not chosen "
                  "from it — assemble.py:333 labels it 'numerical sizing; not literature', the same "
                  "class as the cortex seg_um=0.5 this phase already replaced. ⚠ dx enters cost as "
                  "dx^-5 (cells dx^-3 x CFL subcycles dx^-2), so the choice is worth more than the "
                  "k_xl lever it runs against.",
    },
    "poroelastic_diffusion_um2_s": {
        "value": None, "unit": "um^2/s", "provenance": "SOURCED (band, cross-cell-type)",
        "source": "KB-3.B3.2 (verified): D_p ~ 40-60 um^2/s — HeLa 41+/-11, HT1080 40+/-10, MDCK "
                  "61+/-10. ⚠ The row states in its own text that this is a cross-cell-type anchor "
                  "and NOT breast-specific, so for MCF7 it is a proxy whose scope is recorded. It "
                  "sets the CFL and therefore the dx ceiling; it is required, never defaulted.",
    },
    "dt_phys_s": {
        "value": None, "unit": "s", "provenance": "NUMERICAL",
        "source": "The outer physical-time step the fluid must be subcycled across. Not "
                  "physiological: it is the driver's clock. Required so that the subcycle count "
                  "reported beside a dx is the one that dx would actually pay.",
    },
}


# ── the arithmetic, before any device or any Kind exists ────────────────────────────────────────
def cytosol_counts(
    *,
    dx_um: float,
    radius_um: float,
    poroelastic_diffusion_um2_s: float,
    dt_phys_s: float,
    centre_um: tuple[float, float, float] = (0.0, 0.0, 0.0),
    cfl_safety: float = 0.9,
) -> dict[str, object]:
    """Every count the cytosol build would claim, plus what that resolution costs to step.

    Available with no card and no ``Kind.GRID_CELL``, which is the point: the resolution decision is
    arithmetic over declared axes and can be taken before anything is built.

    The box is the smallest cube that strictly contains the membrane with one OUTSIDE cell layer on
    every side.  One layer is DERIVED, not chosen: ``live_mesh_domain`` requires the live surface to
    stay inside the field box for its winding-query distance bound to hold, and a mask boundary
    condition cannot be expressed on a face that is also the surface.  ⚠ One layer is the PHASE 1
    minimum because nothing moves in PHASE 1; a membrane that moves needs more, and that is PHASE 4's
    to size against a measured excursion rather than this function's to guess.

    Args:
        dx_um: uniform cell size [µm].
        radius_um: the membrane radius the box must contain [µm].
        poroelastic_diffusion_um2_s: ``D_p = mobility / S`` [µm²/s], which sets the CFL.
        dt_phys_s: the outer step the fluid is subcycled across [s].
        centre_um: the membrane centre [µm].
        cfl_safety: the safety factor ``field_grid.cfl_dt`` defaults to.

    Returns:
        Shape, cell count, byte footprint split into fields and mask, the CFL step, the subcycles one
        outer step costs, and the diffusion length the resolution is being judged against.

    Raises:
        ValueError: on a non-finite or non-positive argument, or a ``dx_um`` below the sourced
            continuum floor — below the pore size the field being discretised does not exist.
    """
    for label, value in (("dx_um", dx_um), ("radius_um", radius_um),
                         ("poroelastic_diffusion_um2_s", poroelastic_diffusion_um2_s),
                         ("dt_phys_s", dt_phys_s), ("cfl_safety", cfl_safety)):
        if not math.isfinite(value) or value <= 0.0:
            raise ValueError(f"{label} must be finite and positive; got {value!r}")
    if dx_um < PORE_MESH_BAND_UM[0]:
        raise ValueError(
            f"dx_um={dx_um} um is below the sourced continuum floor {PORE_MESH_BAND_UM[0]} um "
            f"(pore/mesh size {PORE_MESH_BAND_UM[0] * 1e3:.0f}-{PORE_MESH_BAND_UM[1] * 1e3:.0f} nm). "
            "A Darcy pore pressure is not defined below the pore size, so this is not a finer "
            "discretisation of the same field — it is a different field."
        )

    diffusion_length_um = math.sqrt(poroelastic_diffusion_um2_s * dt_phys_s)
    if dx_um > diffusion_length_um:
        raise ValueError(
            f"dx_um={dx_um} um exceeds the poroelastic diffusion length {diffusion_length_um:.4f} um "
            f"(sqrt(D_p={poroelastic_diffusion_um2_s} * dt_phys={dt_phys_s})). One cell would then be "
            "coarser than the transport one outer step carries, so the grid cannot represent the "
            "physics it is being stepped through. This is the band's CEILING, not a tuned threshold: "
            "no number was chosen to make anything pass, and how many cells per diffusion length are "
            "ENOUGH is a separate question this refusal deliberately does not answer — "
            "`cells_per_diffusion_length` is reported so it can be judged."
        )

    half = radius_um + dx_um                      # one OUTSIDE layer, derived (see the docstring)
    n = 2 * int(math.ceil(half / dx_um)) + 1      # odd, so a cell centre sits on the membrane centre
    cells = n ** 3
    lo = tuple(c - 0.5 * (n - 1) * dx_um for c in centre_um)

    cfl_dt_s = cfl_safety * dx_um * dx_um / (2.0 * 3.0 * poroelastic_diffusion_um2_s)
    return {
        "shape": (n, n, n),
        "n_cells": cells,
        "dx_um": dx_um,
        "origin_um": lo,
        "box_half_um": half,
        "field_bytes": 8 * len(_FIELD_NAMES) * cells,
        "mask_bytes": 4 * cells,
        "total_bytes": _BYTES_PER_CELL * cells,
        "cfl_dt_s": cfl_dt_s,
        "subcycles_per_outer_step": math.ceil(dt_phys_s / cfl_dt_s),
        "diffusion_length_um": diffusion_length_um,
        "cells_per_diffusion_length": diffusion_length_um / dx_um,
    }


# ── the build (device; waits on Kind.GRID_CELL) ─────────────────────────────────────────────────
@dataclass(slots=True)
class CytosolField:
    """The cytosol population: ONE grid claim, and the device fields this module owns.

    Attributes:
        population: the name the range is attributed to.
        cells: the arena claim for the ``GRID_CELL`` range.
        shape / dx_um / origin_um: the uniform conservative grid.
        mask: ``(nx, ny, nz)`` int32 OUTSIDE / FLUID / NUCLEUS.
        fields: the float64 3-D fields, by name — allocated and ZERO.  PHASE 1 computes no physics.
        dx_provenance: where ``dx_um`` came from.  Free text, carried into the record so a PI-GAP
            travels inside the artifact rather than beside it.
        fluid_volume_um3: ``sum_FLUID dV`` measured from the mask this build produced.
        reference_volume_um3: ``volume0_um3(membrane) - volume0_um3(envelope)``, measured by those
            builders from what THEY built.  The check is these two against each other; no external
            row is quoted (see the module docstring on KB-3.21).
        field_bytes: device bytes this population's own arrays hold.
    """

    population: str
    cells: Claim
    shape: tuple[int, int, int]
    dx_um: float
    origin_um: tuple[float, float, float]
    mask: wp.array = field(repr=False)
    fields: dict[str, wp.array] = field(repr=False)
    dx_provenance: str
    fluid_volume_um3: float
    reference_volume_um3: float
    field_bytes: int

    @property
    def volume_error(self) -> float:
        """Relative discretisation error of the masked volume against the built surfaces."""
        return abs(self.fluid_volume_um3 - self.reference_volume_um3) / self.reference_volume_um3


def build_cytosol(
    arena: WorldArena,
    *,
    dx_um: float,
    poroelastic_diffusion_um2_s: float,
    dt_phys_s: float,
    membrane: object,
    envelope: object,
    dx_provenance: str,
    population: str = "cytosol",
) -> CytosolField:
    """Claim and stand the cytosol field, classified against the surfaces that were BUILT.

    Not against ideal spheres: ``build/cortex.py``'s rule is that a downstream quantity reads what was
    built rather than what was asked for, and the membrane here is an icosphere whose vertices are
    already in the arena.  The classifier is therefore a winding-number query against a ``wp.Mesh``
    over the live vertices — the algorithm ``components/incumbent/live_mesh_domain.py`` settled and
    the only one in this repository that classifies a grid against a real surface.

    Args:
        arena: the arena to claim from.  Must already hold the membrane and envelope node ranges.
        dx_um: uniform cell size [µm].  No default: see :data:`CYTOSOL_AXES`.
        poroelastic_diffusion_um2_s: ``D_p`` [µm²/s], for the reported CFL.
        dt_phys_s: the outer step [s], for the reported subcycle count.
        membrane / envelope: the built :class:`~aleph.world.build.membrane.ClosedSurface` populations.
        dx_provenance: where ``dx_um`` came from, carried into the artifact.
        population: the claim's name.

    Returns:
        The standing :class:`CytosolField`.

    Raises:
        ArenaKindGapError: while ``arena.Kind`` has no ``GRID_CELL``.  **This is the current state.**
        ValueError: on an axis outside its sourced band, or a grid that does not contain the membrane.
        RuntimeError: if the arena holds no device.
    """
    kind = grid_cell_kind()          # refuses here, before anything is claimed or allocated

    counts = cytosol_counts(
        dx_um=dx_um,
        radius_um=float(membrane.radius_um),
        poroelastic_diffusion_um2_s=poroelastic_diffusion_um2_s,
        dt_phys_s=dt_phys_s,
        centre_um=tuple(float(c) for c in membrane.centre_um),
    )
    if arena.device is None:
        raise RuntimeError(
            "the cytosol field is a device population; the arena holds no device. Warp CUDA is the "
            "only simulation runtime and there is no host path."
        )

    shape = counts["shape"]                                          # type: ignore[assignment]
    cells = arena.claim(population, kind, int(counts["n_cells"]))    # type: ignore[arg-type]
    arena.assert_partitioned()

    with wp.ScopedDevice(arena.device):
        mask = wp.zeros(shape, dtype=wp.int32)
        fields = {name: wp.zeros(shape, dtype=wp.float64) for name in _FIELD_NAMES}
        _classify(arena, mask, shape, counts, membrane, envelope)
        fluid_um3 = _fluid_volume(mask, shape, dx_um)

    reference = float(membrane.volume0_um3) - float(envelope.volume0_um3)
    return CytosolField(
        population=population,
        cells=cells,
        shape=shape,
        dx_um=dx_um,
        origin_um=counts["origin_um"],                               # type: ignore[arg-type]
        mask=mask,
        fields=fields,
        dx_provenance=dx_provenance,
        fluid_volume_um3=fluid_um3,
        reference_volume_um3=reference,
        field_bytes=int(counts["total_bytes"]),                      # type: ignore[arg-type]
    )


@wp.kernel
def _mesh_points_kernel(
    position: wp.array(dtype=wp.vec3d), node_lo: wp.int32, points: wp.array(dtype=wp.vec3)
) -> None:
    """Narrow the arena's float64 vertices to the float32 points ``wp.Mesh`` requires."""
    v = wp.tid()
    p = position[node_lo + v]
    points[v] = wp.vec3(wp.float32(p[0]), wp.float32(p[1]), wp.float32(p[2]))


@wp.kernel
def _local_indices_kernel(
    face_idx: wp.array2d(dtype=wp.int32), node_lo: wp.int32, out: wp.array(dtype=wp.int32)
) -> None:
    """Rebase GLOBAL arena node ids to mesh-local vertex ids, flattened as ``wp.Mesh`` wants them."""
    f = wp.tid()
    for c in range(3):
        out[3 * f + c] = face_idx[f, c] - node_lo


@wp.kernel
def _classify_kernel(
    mesh_id: wp.uint64,
    origin: wp.vec3,
    dx: wp.float32,
    max_dist: wp.float32,
    surface_epsilon: wp.float32,
    inside_code: wp.int32,
    require_fluid: wp.int32,
    mask: wp.array3d(dtype=wp.int32),
) -> None:
    """Write ``inside_code`` where the cell centre is enclosed by the mesh.

    One kernel serves both surfaces.  For the membrane ``require_fluid`` is 0 and OUTSIDE cells are
    written, establishing the classification.  For the envelope it is 1, so only cells already FLUID
    become NUCLEUS and a nucleus lying (impossibly) outside the membrane cannot silently create fluid.

    A cell whose centre lands ON the surface within ``surface_epsilon`` counts as inside: the float32
    winding sign is not trustworthy there, and biasing consistently is what keeps the classification
    from depending on rounding.
    """
    i, j, k = wp.tid()
    if require_fluid == 1 and mask[i, j, k] != _FLUID:
        return
    point = origin + dx * wp.vec3(wp.float32(i), wp.float32(j), wp.float32(k))
    query = wp.mesh_query_point_sign_winding_number(
        mesh_id, point, max_dist, wp.float32(2.0), wp.float32(0.5))
    if not query.result:
        if require_fluid == 0:
            mask[i, j, k] = _OUTSIDE
        return
    closest = wp.mesh_eval_position(mesh_id, query.face, query.u, query.v)
    inside = query.sign < wp.float32(0.0) or wp.length(closest - point) <= surface_epsilon
    if inside:
        mask[i, j, k] = inside_code
    elif require_fluid == 0:
        mask[i, j, k] = _OUTSIDE


@wp.kernel
def _fluid_cells_kernel(mask: wp.array3d(dtype=wp.int32), out: wp.array(dtype=wp.int32)) -> None:
    """Count FLUID cells into one device scalar."""
    i, j, k = wp.tid()
    if mask[i, j, k] == _FLUID:
        wp.atomic_add(out, 0, 1)


def _surface_mesh(arena: WorldArena, surface: object) -> wp.Mesh:
    """A winding-number-capable ``wp.Mesh`` over one built surface's live arena vertices."""
    node_lo = int(surface.nodes.lo)
    n_verts = int(surface.nodes.count)
    n_faces = int(surface.faces.count)
    points = wp.zeros(n_verts, dtype=wp.vec3)
    indices = wp.zeros(3 * n_faces, dtype=wp.int32)
    wp.launch(_mesh_points_kernel, dim=n_verts,
              inputs=[arena.node_arrays["position"], wp.int32(node_lo)], outputs=[points])
    wp.launch(_local_indices_kernel, dim=n_faces,
              inputs=[surface.face_idx, wp.int32(node_lo)], outputs=[indices])
    return wp.Mesh(points=points, indices=indices, support_winding_number=True)


def _classify(
    arena: WorldArena, mask: wp.array, shape: tuple[int, int, int], counts: dict, membrane, envelope
) -> None:
    """Membrane first (FLUID vs OUTSIDE), then the envelope overwrites FLUID with NUCLEUS."""
    lo = counts["origin_um"]
    dx = float(counts["dx_um"])
    origin = wp.vec3(float(lo[0]), float(lo[1]), float(lo[2]))
    # The box diagonal bounds any point-to-surface distance, because the box was sized to contain the
    # surface with a margin. float32 is wp.Mesh's precision, so the tie band is stated, not assumed.
    max_dist = float(math.sqrt(3.0) * (shape[0] - 1) * dx)
    eps = 8.0 * 1.1920929e-07 * max_dist          # 8 float32 ulps over the query span
    for surface, code, require_fluid in ((membrane, _FLUID, 0), (envelope, _NUCLEUS, 1)):
        mesh = _surface_mesh(arena, surface)
        wp.launch(
            _classify_kernel, dim=shape,
            inputs=[mesh.id, origin, wp.float32(dx), wp.float32(max_dist), wp.float32(eps),
                    wp.int32(code), wp.int32(require_fluid)],
            outputs=[mask],
        )


def _fluid_volume(mask: wp.array, shape: tuple[int, int, int], dx_um: float) -> float:
    """``sum_FLUID dV`` [µm³] — one device reduction, read once, between builds and never in a loop."""
    counter = wp.zeros(1, dtype=wp.int32)
    wp.launch(_fluid_cells_kernel, dim=shape, inputs=[mask], outputs=[counter])
    return float(counter.numpy()[0]) * dx_um ** 3


# ── self-check ──────────────────────────────────────────────────────────────────────────────────
def _demo() -> None:
    """Assert the band, the axis discipline and the dx^-5 arithmetic. No device, nothing launched."""
    # 1. The refusal is typed, names the primitive, and is not an AttributeError.
    if not hasattr(Kind, "GRID_CELL"):
        try:
            grid_cell_kind()
        except ArenaKindGapError as exc:
            assert exc.kind_name == "GRID_CELL"
            assert "PENDING DECISION" in str(exc)
        else:                                                   # pragma: no cover - Kind landed
            raise AssertionError("Kind.GRID_CELL is absent but grid_cell_kind() did not refuse")

    # 2. dx has a PI value now, but it is recorded as a decision and is NOT a default: the builder
    #    still requires it, so a record cannot inherit a resolution it never stated.
    assert CYTOSOL_AXES["dx_um"]["value"] == PI_DX_UM
    assert "test point" in str(CYTOSOL_AXES["dx_um"]["provenance"])
    assert signature(build_cytosol).parameters["dx_um"].default is Parameter.empty
    #    The two axes with no PI point still carry none.
    assert CYTOSOL_AXES["poroelastic_diffusion_um2_s"]["value"] is None
    assert CYTOSOL_AXES["dt_phys_s"]["value"] is None

    # 3. Below the continuum floor is refused — a finer grid there is a different field.
    common = dict(radius_um=7.5, poroelastic_diffusion_um2_s=50.0, dt_phys_s=0.05)
    try:
        cytosol_counts(dx_um=0.005, **common)
    except ValueError as exc:
        assert "continuum floor" in str(exc)
    else:
        raise AssertionError("a sub-pore dx was accepted")

    # 3b. Above the diffusion length is refused at the other end of the same band.
    try:
        cytosol_counts(dx_um=3.0, **common)                     # sqrt(50*0.05) = 1.58 um
    except ValueError as exc:
        assert "diffusion length" in str(exc)
    else:
        raise AssertionError("a dx coarser than one step's transport was accepted")

    # 3c. The PI point sits inside the band and resolves the diffusion length several times over.
    pi = cytosol_counts(dx_um=PI_DX_UM, **common)
    assert PORE_MESH_BAND_UM[0] <= PI_DX_UM <= pi["diffusion_length_um"]
    assert pi["cells_per_diffusion_length"] > 1.0

    # 4. dx^-5: halving dx must multiply cells by ~8 and subcycles by ~4, i.e. work by ~32.
    coarse = cytosol_counts(dx_um=0.2, **common)
    fine = cytosol_counts(dx_um=0.1, **common)
    assert 7.0 < fine["n_cells"] / coarse["n_cells"] < 9.0, fine["n_cells"] / coarse["n_cells"]
    assert math.isclose(coarse["cfl_dt_s"] / fine["cfl_dt_s"], 4.0, rel_tol=1e-9)
    work = (fine["n_cells"] * fine["subcycles_per_outer_step"]) / (
        coarse["n_cells"] * coarse["subcycles_per_outer_step"])
    assert 28.0 < work < 36.0, work

    # 5. The box strictly contains the membrane with at least one cell to spare on every side.
    assert coarse["box_half_um"] >= 7.5 + coarse["dx_um"]
    assert coarse["shape"][0] % 2 == 1                          # a centre cell exists

    # 6. The byte split adds up and names its arrays.
    assert coarse["field_bytes"] + coarse["mask_bytes"] == coarse["total_bytes"]

    if hasattr(Kind, "GRID_CELL"):
        print(f"cytosol: dx = {PI_DX_UM} um (PI 2026-08-20, TEST POINT, not derived) -> "
              f"{pi['shape'][0]}^3 = {pi['n_cells']:,} cells, {pi['total_bytes'] / 1e9:.3f} GB, "
              f"{pi['subcycles_per_outer_step']:,} subcycles/step")
    else:                                                       # pragma: no cover - Kind absent
        print(f"cytosol: REFUSING — {GRID_CELL_DECISION.splitlines()[0]}")
    for dx in (0.5, 0.25, 0.1, 0.05):
        c = cytosol_counts(dx_um=dx, **common)
        print(f"  dx={dx:<5} shape={c['shape'][0]:>4}^3  cells={c['n_cells']:>12,}  "
              f"{c['total_bytes'] / 1e9:6.3f} GB  cfl_dt={c['cfl_dt_s']:.2e} s  "
              f"subcycles={c['subcycles_per_outer_step']:>6,}")
    print("self-check OK")


if __name__ == "__main__":
    _demo()
