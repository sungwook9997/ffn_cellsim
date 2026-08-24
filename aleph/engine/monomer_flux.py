"""G-actin conserved-pool CHEMICAL_FLUX connector — whole-cell slice-3 (owns the Cytosol MonomerField seam).

Spec §4 (``WHOLE_CELL_COMMON_CONTRACTS_SPEC_2026-07-23.md``) + roadmap §3 item 7.  The cytosolic **G-actin**
monomer field is owned by ONE component (Cytosol, backed by :class:`~aleph.components.fluid.transport.MonomerField`).
Cortex / SF / lamellipodium / filopodium each **draw and return** monomer through an EXPLICIT per-consumer
chemical-flux connector — never a shared array — and the whole-cell conservation

    A_total  =  integral phi c dV  +  sum_over_consumers  bound_c

is asserted every accepted step on the ledger **mass channel** (:meth:`GlobalCellLedger.add_mass`).

Three lanes, mirroring ``ac/engine/population.py`` (host accounting) + ``ac/fluid/transport_reference.py``
(NumPy physics) + a device kernel (CUDA), so the whole invariant is gated on the dev Mac:

* :class:`MonomerFluxAccount` — a PURE HOST accounting structure (CPU-importable, no device): the field's
  available-monomer scalar plus one bound pool per consumer, with the conservative ``draw``/``release``
  discipline and the disjoint-consumer guarantee.  This is what the device kernel mirrors and what the
  structural gate exercises directly (no CUDA needed).
* :class:`MonomerPoolOwner` — the Cytosol-owned conserved-pool transaction participant + ledger contributor
  (CUDA lane): wraps the ``MonomerField`` device state, exposes it to consumers, and pushes ``A_total`` to the
  mass channel.  It owns the field; consumers never touch ``u`` directly.
* :class:`MonomerFluxConnector` — one CHEMICAL_FLUX connector runtime per consumer (CUDA lane): ``propose_events``
  raises a candidate flux RATE only (no force, no host state); ``commit_irreversible`` materialises the
  conservative scalar exchange under the accepted predicate (``-drawn`` off the field pool, ``+drawn`` onto the
  consumer's bound pool — ``A_total`` invariant); ``rollback`` restores the snapshot on a rejected step.

Device work is CUDA-only (I0-A); the Warp kernels here are import-safe SOURCE (lazy JIT) and construct only
against CUDA arrays.  The flux RATE is a sourced-GAP: it defaults to zero (no flux) until ``CellState`` supplies
a KB-sourced polymerisation/depolymerisation rate — no magic number is invented here.

Sanity Gate:
    * ownership: the field pool and each consumer bound pool are DISTINCT scalars; a draw names exactly one
      consumer and moves monomer only between the field and that consumer's pool (no shared array).
    * conservation: every ``draw``/``release`` is mass-neutral to ``A_total``; ``commit_irreversible`` applies
      ``-drawn``/``+drawn`` under one device predicate so an accepted step conserves ``A_total`` to round-off.
    * boundary/sign: a rejected step leaves the field and every bound pool at the snapshot (device-gated
      rollback); the flux clock/epoch never advances on a rejected candidate.
"""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field

import warp as wp

from aleph.engine.contracts import ConnectorFamily

__all__ = [
    "MonomerFluxAccount",
    "MonomerPoolOwner",
    "MonomerFluxConnector",
    "monomer_field_content_reduce_kernel",
    "monomer_flux_snapshot_kernel",
    "monomer_flux_rollback_kernel",
    "monomer_flux_commit_kernel",
]

_FLUID = wp.constant(1)
CYTOSOL_COMPONENT = "cytosol"


def _device_is_cuda(array: object) -> bool:
    """Return whether an array-like object declares a CUDA device."""
    device = getattr(array, "device", None)
    return bool(getattr(device, "is_cuda", False))


def _storage_key(array: object) -> tuple[object, ...]:
    """Return a conservative storage identity without copying device data."""
    ptr = getattr(array, "ptr", None)
    if ptr is None:
        return ("object", id(array))
    return ("device-ptr", str(getattr(array, "device", None)), int(ptr))


# ======================================================================================================
# Host accounting lane (CPU-importable; the device kernel mirrors this)
# ======================================================================================================
@dataclass(slots=True)
class MonomerFluxAccount:
    """Pure-host conserved G-actin pool: one field scalar + one bound pool per consumer.

    The field pool holds the free-monomer mass available in the Cytosol; each consumer owns a disjoint bound
    pool (polymerised monomer).  ``draw`` moves monomer field->consumer (polymerisation); ``release`` moves it
    back (depolymerisation).  Both are exact and mass-neutral to :attr:`a_total`, and a consumer can only ever
    touch its OWN pool — the "no shared array" property, enforced by name lookup.

    Args:
        field_monomer: initial free-monomer mass in the Cytosol pool (a GAP for the native magnitude — the
            sourced ``c0`` x volume — but exact for this accounting gate); must be finite and nonnegative.
        consumers: the disjoint consumer names owning bound pools.
    """

    field_monomer: float
    consumers: Sequence[str]
    _bound: dict[str, float] = field(default_factory=dict, init=False, repr=False)
    _a_total0: float = field(default=0.0, init=False, repr=False)

    def __post_init__(self) -> None:
        if not math.isfinite(self.field_monomer) or self.field_monomer < 0.0:
            raise ValueError("field_monomer must be finite and nonnegative")
        names = tuple(self.consumers)
        if len(set(names)) != len(names):
            raise ValueError("consumer names must be unique (each owns a disjoint bound pool)")
        if not names:
            raise ValueError("need at least one monomer consumer")
        self.consumers = names
        self._bound = {name: 0.0 for name in names}
        self._a_total0 = float(self.field_monomer)

    @property
    def bound(self) -> Mapping[str, float]:
        """Read-only view of each consumer's bound-monomer count."""
        return dict(self._bound)

    @property
    def a_total(self) -> float:
        """Conserved total actin ``field_monomer + sum_c bound_c``."""
        return float(self.field_monomer + sum(self._bound.values()))

    def draw(self, consumer: str, quanta: float) -> None:
        """Polymerise: move ``quanta`` (>= 0) from the field pool to ``consumer``'s bound pool.

        Raises ``KeyError`` for an unknown consumer (no implicit shared pool) and ``ValueError`` if the field
        pool cannot cover the draw (monomer is finite — a consumer cannot conjure polymer from an empty pool).
        """
        if consumer not in self._bound:
            raise KeyError(f"unknown monomer consumer {consumer!r}; known: {tuple(self._bound)}")
        if not math.isfinite(quanta) or quanta < 0.0:
            raise ValueError("draw quanta must be finite and nonnegative")
        if quanta > self.field_monomer + 1e-12:
            raise ValueError(
                f"{consumer!r} draw {quanta!r} exceeds the free-monomer pool {self.field_monomer!r} "
                "(finite G-actin — the pool is conserved, never resampled)"
            )
        self.field_monomer -= quanta
        self._bound[consumer] += quanta

    def release(self, consumer: str, quanta: float) -> None:
        """Depolymerise: move ``quanta`` (>= 0) from ``consumer``'s bound pool back to the field pool."""
        if consumer not in self._bound:
            raise KeyError(f"unknown monomer consumer {consumer!r}; known: {tuple(self._bound)}")
        if not math.isfinite(quanta) or quanta < 0.0:
            raise ValueError("release quanta must be finite and nonnegative")
        if quanta > self._bound[consumer] + 1e-12:
            raise ValueError(f"{consumer!r} release {quanta!r} exceeds its bound pool {self._bound[consumer]!r}")
        self._bound[consumer] -= quanta
        self.field_monomer += quanta

    def assert_conserved(self, atol: float = 1e-9) -> None:
        """Assert ``A_total`` has not drifted from its initial value (mass-channel conservation gate)."""
        if abs(self.a_total - self._a_total0) > atol:
            raise AssertionError(
                f"A_total drifted {self.a_total!r} != {self._a_total0!r}: a flux was not mass-neutral "
                "(the field<->bound exchange leaked)"
            )


# ======================================================================================================
# Device kernels (Warp-CUDA SOURCE; import-safe, construct only on CUDA)
# ======================================================================================================
@wp.kernel
def monomer_field_content_reduce_kernel(
    u: wp.array3d(dtype=wp.float64),
    mask: wp.array3d(dtype=wp.int32),
    cell_vol: wp.float64,
    out_content: wp.array(dtype=wp.float64),
) -> None:
    """Atomically reduce ``sum_FLUID u * cell_vol`` (the Cytosol free-monomer field content) into ``out[0]``."""
    i, j, k = wp.tid()
    if mask[i, j, k] == _FLUID:
        wp.atomic_add(out_content, 0, u[i, j, k] * cell_vol)


@wp.kernel
def monomer_flux_snapshot_kernel(
    field_pool: wp.array(dtype=wp.float64),
    consumer_bound: wp.array(dtype=wp.float64),
    net_drawn: wp.array(dtype=wp.float64),
    field_pool_snap: wp.array(dtype=wp.float64),
    consumer_bound_snap: wp.array(dtype=wp.float64),
    net_drawn_snap: wp.array(dtype=wp.float64),
) -> None:
    """Snapshot the reversible flux scalars before candidate mutation (one thread)."""
    field_pool_snap[0] = field_pool[0]
    consumer_bound_snap[0] = consumer_bound[0]
    net_drawn_snap[0] = net_drawn[0]


@wp.kernel
def monomer_flux_rollback_kernel(
    accepted: wp.array(dtype=wp.int32),
    field_pool_snap: wp.array(dtype=wp.float64),
    consumer_bound_snap: wp.array(dtype=wp.float64),
    net_drawn_snap: wp.array(dtype=wp.float64),
    field_pool: wp.array(dtype=wp.float64),
    consumer_bound: wp.array(dtype=wp.float64),
    net_drawn: wp.array(dtype=wp.float64),
) -> None:
    """Restore the snapshot iff ``accepted[0] == 0`` — a rejected candidate re-exposes the prior pools."""
    if accepted[0] == wp.int32(0):
        field_pool[0] = field_pool_snap[0]
        consumer_bound[0] = consumer_bound_snap[0]
        net_drawn[0] = net_drawn_snap[0]


@wp.kernel
def monomer_flux_commit_kernel(
    accepted: wp.array(dtype=wp.int32),
    drawn_candidate: wp.array(dtype=wp.float64),
    field_pool: wp.array(dtype=wp.float64),
    consumer_bound: wp.array(dtype=wp.float64),
    net_drawn: wp.array(dtype=wp.float64),
) -> None:
    """Move the candidate draw field->bound iff ``accepted[0] != 0`` (one thread; A_total invariant).

    ``drawn_candidate[0]`` is the accepted-step monomer flux (polymerisation positive, depolymerisation
    negative), proposed as a RATE * dt by :meth:`MonomerFluxConnector.propose_events`.  The exchange subtracts
    from the field pool and adds to the consumer's bound pool by the SAME amount, so the whole-cell
    ``A_total = field + sum_c bound`` is conserved on the device with no host read of the predicate.
    """
    if accepted[0] != wp.int32(0):
        d = drawn_candidate[0]
        field_pool[0] = field_pool[0] - d
        consumer_bound[0] = consumer_bound[0] + d
        net_drawn[0] = net_drawn[0] + d


# ======================================================================================================
# Cytosol conserved-pool owner (CUDA lane)
# ======================================================================================================
@dataclass(frozen=True, slots=True)
class MonomerPoolOwner:
    """Cytosol-owned G-actin conserved pool: the ``MonomerField`` device state + the mass-channel gate.

    Args:
        monomer_field: the :class:`~aleph.components.fluid.transport.MonomerField` (owns ``u``/``reaction`` on the
            grid); this owner exposes it and reports ``A_total`` — consumers never write ``u`` directly.
        field_pool_d: ``(1,)`` float64 CUDA scalar — the free-monomer mass available for polymerisation flux
            (the conserved counterpart of the consumers' bound pools; ``A_total = field_pool + sum bound``).
        a_total_d: ``(1,)`` float64 CUDA scalar — scratch the mass reduce writes before the ledger push.
    """

    monomer_field: object
    field_pool_d: wp.array
    a_total_d: wp.array

    def __post_init__(self) -> None:
        grid = getattr(self.monomer_field, "grid", None)
        if grid is None or not hasattr(self.monomer_field, "u"):
            raise TypeError("MonomerPoolOwner requires an ac/fluid MonomerField")
        for label, arr in (("field_pool_d", self.field_pool_d), ("a_total_d", self.a_total_d)):
            if getattr(arr, "dtype", None) != wp.float64:
                raise TypeError(f"MonomerPoolOwner.{label} must be wp.float64")
            if getattr(arr, "shape", None) != (1,):
                raise ValueError(f"MonomerPoolOwner.{label} must be a single-entry scalar")
            if not _device_is_cuda(arr):
                raise ValueError(f"MonomerPoolOwner.{label} must be a Warp CUDA device array")
        if _storage_key(self.field_pool_d) == _storage_key(self.a_total_d):
            raise ValueError("field_pool_d and a_total_d must own distinct storage")

    @property
    def name(self) -> str:
        """The owning component (always the Cytosol)."""
        return CYTOSOL_COMPONENT

    @property
    def device(self) -> str:
        """CUDA device the pool lives on."""
        return str(self.field_pool_d.device)

    def accumulate_ledger(self, ledger: object) -> None:
        """Reduce the field content + free-monomer pool into the ledger's mass channel (A_total invariant).

        The consumers' bound pools are reduced by each :class:`MonomerFluxConnector` (below); together the
        mass channel accumulates ``A_total`` for the global conservation gate.  No host read of device state.
        """
        grid = self.monomer_field.grid
        cell_vol = float(grid.dx) ** int(grid.dim)
        self.a_total_d.zero_()
        wp.launch(
            monomer_field_content_reduce_kernel,
            dim=grid.shape,
            inputs=[self.monomer_field.u, grid.mask, wp.float64(cell_vol)],
            outputs=[self.a_total_d],
            device=self.device,
        )
        add_mass = getattr(ledger, "add_mass", None)
        if callable(add_mass):
            add_mass(self.a_total_d)   # field content ...
            add_mass(self.field_pool_d)  # ... + free-monomer pool (bound pools added by the connectors)


# ======================================================================================================
# CHEMICAL_FLUX connector (CUDA lane)
# ======================================================================================================
@dataclass(frozen=True, slots=True)
class MonomerFluxConnector:
    """One explicit G-actin chemical-flux connector: a consumer component <-> the Cytosol monomer pool.

    ``family`` is :attr:`ConnectorFamily.CHEMICAL_FLUX`; ``kinetics``/``commit_on_accept`` are True (the count
    exchange is irreversible and commits only on an accepted step).  The connector owns the reversible flux
    scalars and mirrors :class:`MonomerFluxAccount` on the device: ``propose_events`` writes the candidate
    ``drawn`` from a RATE (no force, no host state); ``commit_irreversible`` moves it field->bound under the
    accepted predicate; ``rollback`` restores the snapshot on rejection.

    Args:
        name: the connector name (e.g. ``"cortex_gactin_flux"``).
        component_a: the consumer component (``"cortex"``/``"sf_arc"``/``"lamellipodium"``/``"filopodium"``).
        pool: the Cytosol :class:`MonomerPoolOwner` whose ``field_pool_d`` this connector draws from.
        consumer_bound_d: ``(1,)`` float64 CUDA scalar — the consumer's own bound-monomer pool (distinct from
            the field pool and from every other consumer's pool; the "no shared array" invariant).
        drawn_candidate_d: ``(1,)`` float64 CUDA scalar — this step's proposed flux (rate * dt), set by
            ``propose_events``.
        net_drawn_d: ``(1,)`` float64 CUDA scalar — running committed flux (diagnostic / event log).
        rate_d: ``(1,)`` float64 CUDA scalar — the local polymerisation rate CellState sets (defaults zero:
            a sourced GAP, never an invented magic number).
        snapshot_d: ``(3,)`` float64 CUDA scratch — the (field, bound, net) snapshot for rollback.
    """

    name: str
    component_a: str
    pool: MonomerPoolOwner
    consumer_bound_d: wp.array
    drawn_candidate_d: wp.array
    net_drawn_d: wp.array
    rate_d: wp.array
    snapshot_d: wp.array
    component_b: str = CYTOSOL_COMPONENT
    family: ConnectorFamily = ConnectorFamily.CHEMICAL_FLUX

    def __post_init__(self) -> None:
        if not self.name.strip():
            raise ValueError("chemical-flux connector needs a name")
        if self.component_a == CYTOSOL_COMPONENT:
            raise ValueError("a monomer-flux consumer must be a distinct component from the cytosol")
        scalars = (
            ("consumer_bound_d", self.consumer_bound_d, (1,)),
            ("drawn_candidate_d", self.drawn_candidate_d, (1,)),
            ("net_drawn_d", self.net_drawn_d, (1,)),
            ("rate_d", self.rate_d, (1,)),
            ("snapshot_d", self.snapshot_d, (3,)),
        )
        for label, arr, shape in scalars:
            if getattr(arr, "dtype", None) != wp.float64:
                raise TypeError(f"MonomerFluxConnector.{label} must be wp.float64")
            if getattr(arr, "shape", None) != shape:
                raise ValueError(f"MonomerFluxConnector.{label} must have shape {shape}")
            if not _device_is_cuda(arr):
                raise ValueError(f"MonomerFluxConnector.{label} must be a Warp CUDA device array")
        # No shared array: the consumer bound pool must not alias the Cytosol field pool.
        if _storage_key(self.consumer_bound_d) == _storage_key(self.pool.field_pool_d):
            raise ValueError(
                "consumer_bound_d aliases the Cytosol field pool — a chemical-flux connector must move "
                "monomer between DISTINCT device arrays, never share one array"
            )
        if str(self.consumer_bound_d.device) != self.pool.device:
            raise ValueError("consumer bound pool and Cytosol field pool must share one CUDA device")

    @property
    def device(self) -> str:
        """CUDA device shared by the connector scalars and the Cytosol pool."""
        return self.pool.device

    def propose_events(self, rates: object, dt_phys: float, rng_seed: int, neighbors: object) -> None:
        """Write the candidate monomer flux ``drawn = rate * dt`` (RATE only; no force, no host state).

        The rate lives in :attr:`rate_d` (set by ``CellState`` in the biology phase; zero until sourced).  A
        deterministic Poisson/tau-leap draw folds ``rng_seed`` + the accepted-step epoch in the native lane;
        the mean-flux candidate here is ``rate * dt`` and materialises only inside ``commit_irreversible``.
        """
        if not math.isfinite(dt_phys) or dt_phys <= 0.0:
            raise ValueError("dt_phys must be finite and positive")
        wp.launch(
            _propose_monomer_flux_kernel,
            dim=1,
            inputs=[self.rate_d, wp.float64(dt_phys)],
            outputs=[self.drawn_candidate_d],
            device=self.device,
        )

    def snapshot_candidate(self) -> None:
        """Snapshot the reversible flux scalars before candidate mutation."""
        wp.launch(
            monomer_flux_snapshot_kernel,
            dim=1,
            inputs=[self.pool.field_pool_d, self.consumer_bound_d, self.net_drawn_d],
            outputs=[
                self.snapshot_d[0:1],
                self.snapshot_d[1:2],
                self.snapshot_d[2:3],
            ],
            device=self.device,
        )

    def rollback(self, accepted: wp.array) -> None:
        """Restore the snapshot under the scheduler-owned rejected predicate (device-gated)."""
        wp.launch(
            monomer_flux_rollback_kernel,
            dim=1,
            inputs=[
                accepted,
                self.snapshot_d[0:1],
                self.snapshot_d[1:2],
                self.snapshot_d[2:3],
            ],
            outputs=[self.pool.field_pool_d, self.consumer_bound_d, self.net_drawn_d],
            device=self.device,
        )

    def commit_irreversible(self, accepted: wp.array, dt_phys: float, rng_seed: int) -> None:
        """Move the candidate draw field->bound under the accepted predicate (A_total invariant)."""
        wp.launch(
            monomer_flux_commit_kernel,
            dim=1,
            inputs=[accepted, self.drawn_candidate_d],
            outputs=[self.pool.field_pool_d, self.consumer_bound_d, self.net_drawn_d],
            device=self.device,
        )

    def accumulate_ledger(self, ledger: object) -> None:
        """Reduce this consumer's bound pool into the mass channel (completes ``A_total`` with the field)."""
        add_mass = getattr(ledger, "add_mass", None)
        if callable(add_mass):
            add_mass(self.consumer_bound_d)

    def accumulate(self, pos: wp.array, force: wp.array) -> None:
        """Mechanical coupling hook — the monomer transfer's drag rides the existing IMMERSED_TRANSFER stencil.

        The chemical count exchange is mass, not force; the connector adds no mechanical force here (the drag
        of transported monomer on the pore fluid is the ``*_cytosol_transfer`` IMMERSED_TRANSFER connector's
        job).  The hook exists so the actor's mechanics-binding check is satisfied without a phantom force.
        """
        return None


@wp.kernel
def _propose_monomer_flux_kernel(
    rate: wp.array(dtype=wp.float64),
    dt_phys: wp.float64,
    drawn_candidate: wp.array(dtype=wp.float64),
) -> None:
    """Candidate flux ``drawn = rate * dt`` (one thread; rate is CellState-set, zero until sourced)."""
    drawn_candidate[0] = rate[0] * dt_phys
