"""Ensemble validation harness + per-run event-log schema for the stochastic whole-cell runtime.

Spec §7 (``WHOLE_CELL_COMMON_CONTRACTS_SPEC_2026-07-23.md``) — spine item #4 of the ratified commit
order, built **before** any stochastic subsystem.  It lands **zero new physics**: no constants, no
rates, no forces.  It is orchestration + statistics only.

Two things the spec requires:

1. **Per-run event-log schema** (:class:`EventRecord` / :class:`EventLog`).  Each event carries
   ``{type, position|region, start_time, lifetime, amplitude|rate_modifier, mechanistic_trigger}`` plus
   the run's :class:`~aleph.engine.cell_state.CellState` tag, written **out of the hot loop**.  The
   two either/or pairs (position XOR region, amplitude XOR rate_modifier) are validated at construction.

2. **Ensemble harness** (:class:`EnsembleRunner` → :class:`EnsembleSummary`).  A stochastic runtime is
   validated by **steady-state DISTRIBUTIONS over N realizations** (connectivity %, γ, filament count,
   length distribution) — never single-run asserts, never a fit.  Each realization (trajectory) carries
   exactly ONE CellState; :func:`mix_ensembles` forms the experimental distribution by MIXING
   state-conditioned ensembles at their observed proportions.

The harness never does physics.  Mirroring
:class:`~aleph.engine.transaction.CellTransaction` (caller-owned ``solve``), the caller supplies a
``simulate(rng, cell_state, realization_index) -> RealizationResult`` callable; the runner only builds the
independent per-realization RNG, invokes it, and aggregates.  This module imports no ``warp`` and runs no
CUDA — it is CPU-importable and drives the structural gates directly.

Sanity Gate:
    * reproducibility: :meth:`EnsembleRunner.seeds` is a pure function of ``base_seed`` +
      ``n_realizations`` (``np.random.SeedSequence(base_seed).spawn(...)``); same base_seed → identical
      seed list and identical summary numbers (when ``simulate`` uses the passed rng).
    * independence: spawned SeedSequence children give distinct, decorrelated per-realization streams.
    * event validation: :class:`EventRecord` rejects both/neither of position|region and
      amplitude|rate_modifier, and rejects negative / non-finite times.
    * band fraction: :meth:`ObservableDistribution.within_band` reports the fraction of samples in a
      literature band — the gate compares a DISTRIBUTION to a band, never a single value.
    * mixture: :func:`mix_ensembles` requires proportions > 0 summing to 1 (within 1e-9) and mixes only
      observables present in ALL summaries, by deterministic proportional allocation (no hidden RNG).
    * no physics: the runner asserts each returned result's tag + index match this ensemble and performs
      no mechanical/topological computation itself.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from types import MappingProxyType
from typing import Any, Callable

import numpy as np

from aleph.engine.cell_state import CellState

__all__ = [
    "EventRecord",
    "EventLog",
    "RealizationResult",
    "ObservableDistribution",
    "EnsembleSummary",
    "EnsembleRunner",
    "mix_ensembles",
]

# Number of axes in a CellState tag (see cell_state.CellState.axes).
_N_AXES = 6


def _is_real(value: object) -> bool:
    """Return whether ``value`` is a real (non-bool) number."""
    return not isinstance(value, bool) and isinstance(value, (int, float))


def _is_finite(value: float) -> bool:
    """Return whether ``value`` is a finite real number (rejects NaN and inf)."""
    return _is_real(value) and float(value) == float(value) and abs(float(value)) != float("inf")


def _coerce_tag(tag: object) -> tuple[str, ...]:
    """Return a validated 6-tuple CellState tag from a :class:`CellState` or a raw 6-tuple of strings.

    Args:
        tag: a :class:`CellState` (its ``axes`` are used) or a length-6 sequence of non-empty strings.

    Returns:
        The immutable 6-tuple of axis labels.

    Raises:
        TypeError: if ``tag`` is neither a CellState nor a sequence of strings.
        ValueError: if the tag does not have exactly 6 non-empty string entries.
    """
    if isinstance(tag, CellState):
        return tuple(tag.axes)
    if isinstance(tag, str) or not isinstance(tag, Sequence):
        raise TypeError("cell_state_tag must be a CellState or a 6-tuple of non-empty strings")
    values = tuple(tag)
    if len(values) != _N_AXES:
        raise ValueError(f"cell_state_tag must have {_N_AXES} axes, got {len(values)}")
    for entry in values:
        if not isinstance(entry, str) or not entry.strip():
            raise ValueError("every cell_state_tag axis must be a non-empty string")
    return values


def _require_nonneg_int(value: object, *, what: str) -> int:
    """Return ``value`` as a nonnegative int or raise."""
    if isinstance(value, bool) or not isinstance(value, int):
        raise TypeError(f"{what} must be an int")
    if value < 0:
        raise ValueError(f"{what} must be nonnegative")
    return value


def _as_1d_float_array(values: object, *, what: str) -> np.ndarray:
    """Return a defensively-copied 1-D float ndarray of finite values or raise."""
    array = np.asarray(values, dtype=float)
    if array.ndim != 1:
        raise ValueError(f"{what} must be 1-D, got shape {array.shape}")
    if array.size and not np.all(np.isfinite(array)):
        raise ValueError(f"{what} must be finite (no NaN/inf)")
    return np.array(array, dtype=float, copy=True)


@dataclass(frozen=True, slots=True)
class EventRecord:
    """One logged stochastic event, written out of the hot loop (spec §7).

    Exactly one of ``position`` / ``region`` must be set (where the event happened), and exactly one of
    ``amplitude`` / ``rate_modifier`` must be set (how strong it was) — the other of each pair stays
    ``None``.  These mutual-exclusions are enforced at construction.

    Args:
        event_type: the event kind (e.g. ``"nucleation"``, ``"myosin_pulse"``); non-empty.
        start_time: event onset [s]; finite and ``>= 0``.
        lifetime: event duration [s]; finite and ``>= 0``.
        mechanistic_trigger: the mechanistic cause (e.g. ``"g_actin_binding"``); non-empty.
        position: ``(x, y, z)`` site [µm], XOR ``region``.
        region: named region label (e.g. ``"cortex_apical"``), XOR ``position``.
        amplitude: event magnitude in its own units, XOR ``rate_modifier``.
        rate_modifier: dimensionless multiplier on a local rate, XOR ``amplitude``.
    """

    event_type: str
    start_time: float
    lifetime: float
    mechanistic_trigger: str
    position: tuple[float, float, float] | None = None
    region: str | None = None
    amplitude: float | None = None
    rate_modifier: float | None = None

    def __post_init__(self) -> None:
        if not isinstance(self.event_type, str) or not self.event_type.strip():
            raise ValueError("event_type must be a non-empty string")
        if not isinstance(self.mechanistic_trigger, str) or not self.mechanistic_trigger.strip():
            raise ValueError("mechanistic_trigger must be a non-empty string")
        for label, value in (("start_time", self.start_time), ("lifetime", self.lifetime)):
            if not _is_finite(value):
                raise ValueError(f"{label} must be finite (no NaN/inf)")
            if value < 0:
                raise ValueError(f"{label} must be nonnegative")

        # Location: exactly one of position / region.
        if (self.position is None) == (self.region is None):
            raise ValueError("EventRecord needs exactly one of position or region (not both, not neither)")
        if self.position is not None:
            coords = tuple(self.position)
            if len(coords) != 3 or not all(_is_finite(c) for c in coords):
                raise ValueError("position must be a finite (x, y, z) triple")
            object.__setattr__(self, "position", (float(coords[0]), float(coords[1]), float(coords[2])))
        if self.region is not None and (not isinstance(self.region, str) or not self.region.strip()):
            raise ValueError("region must be a non-empty string")

        # Magnitude: exactly one of amplitude / rate_modifier.
        if (self.amplitude is None) == (self.rate_modifier is None):
            raise ValueError(
                "EventRecord needs exactly one of amplitude or rate_modifier (not both, not neither)"
            )
        for label, value in (("amplitude", self.amplitude), ("rate_modifier", self.rate_modifier)):
            if value is not None:
                if not _is_finite(value):
                    raise ValueError(f"{label} must be finite (no NaN/inf)")
                object.__setattr__(self, label, float(value))

    def to_dict(self) -> dict[str, Any]:
        """Return a plain json-serializable dict for out-of-hot-loop persistence."""
        return {
            "event_type": self.event_type,
            "start_time": self.start_time,
            "lifetime": self.lifetime,
            "mechanistic_trigger": self.mechanistic_trigger,
            "position": list(self.position) if self.position is not None else None,
            "region": self.region,
            "amplitude": self.amplitude,
            "rate_modifier": self.rate_modifier,
        }


@dataclass(slots=True)
class EventLog:
    """Per-realization event history tagged with its CellState, seed, and realization index (spec §7).

    The log is appended to out of the hot loop; :meth:`to_dict` serializes the whole run for persistence.

    Args:
        cell_state_tag: the run's CellState tag (a :class:`CellState` or a 6-tuple of non-empty strings).
        base_seed: the ensemble base seed this realization descends from; nonnegative int.
        realization_index: this realization's index within the ensemble; nonnegative int.
    """

    cell_state_tag: Any
    base_seed: int
    realization_index: int
    _records: list[EventRecord] = field(default_factory=list, init=False, repr=False)

    def __post_init__(self) -> None:
        object.__setattr__(self, "cell_state_tag", _coerce_tag(self.cell_state_tag))
        self.base_seed = _require_nonneg_int(self.base_seed, what="base_seed")
        self.realization_index = _require_nonneg_int(self.realization_index, what="realization_index")

    def append(self, record: EventRecord) -> None:
        """Append one :class:`EventRecord` to the log."""
        if not isinstance(record, EventRecord):
            raise TypeError("EventLog.append expects an EventRecord")
        self._records.append(record)

    @property
    def records(self) -> tuple[EventRecord, ...]:
        """Read-only snapshot of the logged events in append order."""
        return tuple(self._records)

    def to_dict(self) -> dict[str, Any]:
        """Return a plain json-serializable dict of the tagged event history."""
        return {
            "cell_state_tag": list(self.cell_state_tag),
            "base_seed": self.base_seed,
            "realization_index": self.realization_index,
            "records": [record.to_dict() for record in self._records],
        }


@dataclass(frozen=True, slots=True)
class RealizationResult:
    """One realization's steady-state observables — the unit the caller's ``simulate`` returns.

    Args:
        cell_state_tag: the realization's CellState tag (a :class:`CellState` or a 6-tuple of strings).
        seed: the concrete per-realization seed this trajectory ran under; nonnegative int.
        realization_index: this realization's index within the ensemble; nonnegative int.
        scalars: scalar observables, e.g. ``{"connectivity_pct": .., "gamma": .., "filament_count": ..}``;
            each value a finite float.
        distributions: 1-D sample distributions, e.g. ``{"length_um": array}``; each a 1-D float ndarray.
        event_log: optional :class:`EventLog` for this realization.
    """

    cell_state_tag: Any
    seed: int
    realization_index: int
    scalars: Mapping[str, float] = field(default_factory=dict)
    distributions: Mapping[str, np.ndarray] = field(default_factory=dict)
    event_log: EventLog | None = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "cell_state_tag", _coerce_tag(self.cell_state_tag))
        object.__setattr__(self, "seed", _require_nonneg_int(self.seed, what="seed"))
        object.__setattr__(
            self, "realization_index", _require_nonneg_int(self.realization_index, what="realization_index")
        )
        scalars: dict[str, float] = {}
        for name, value in dict(self.scalars).items():
            if not isinstance(name, str) or not name.strip():
                raise ValueError("scalar observable names must be non-empty strings")
            if not _is_finite(value):
                raise ValueError(f"scalar {name!r} must be a finite float")
            scalars[name] = float(value)
        object.__setattr__(self, "scalars", MappingProxyType(scalars))
        dists: dict[str, np.ndarray] = {}
        for name, values in dict(self.distributions).items():
            if not isinstance(name, str) or not name.strip():
                raise ValueError("distribution observable names must be non-empty strings")
            dists[name] = _as_1d_float_array(values, what=f"distribution {name!r}")
        object.__setattr__(self, "distributions", MappingProxyType(dists))
        if self.event_log is not None and not isinstance(self.event_log, EventLog):
            raise TypeError("event_log must be an EventLog or None")


@dataclass(frozen=True, slots=True)
class ObservableDistribution:
    """A 1-D empirical distribution of one observable, compared to a literature BAND (never a fit).

    Args:
        name: the observable name.
        values: 1-D float samples (per-realization scalars, or pooled distribution samples).
    """

    name: str
    values: np.ndarray

    def __post_init__(self) -> None:
        if not isinstance(self.name, str) or not self.name.strip():
            raise ValueError("ObservableDistribution needs a non-empty name")
        object.__setattr__(self, "values", _as_1d_float_array(self.values, what=f"{self.name!r} values"))

    @property
    def count(self) -> int:
        """Number of samples."""
        return int(self.values.size)

    @property
    def mean(self) -> float:
        """Sample mean (``nan`` if empty)."""
        if self.values.size == 0:
            return float("nan")
        return float(np.mean(self.values))

    @property
    def std(self) -> float:
        """Sample standard deviation (ddof=1 when n>1, else 0.0)."""
        if self.values.size <= 1:
            return 0.0
        return float(np.std(self.values, ddof=1))

    def quantiles(self, levels: Sequence[float] = (0.05, 0.5, 0.95)) -> tuple[float, ...]:
        """Return the empirical quantiles at ``levels`` (each in ``[0, 1]``)."""
        levels = tuple(levels)
        for level in levels:
            if not _is_finite(level) or not 0.0 <= level <= 1.0:
                raise ValueError("quantile levels must lie in [0, 1]")
        if self.values.size == 0:
            return tuple(float("nan") for _ in levels)
        return tuple(float(q) for q in np.quantile(self.values, levels))

    def histogram(self, bins: int | Sequence[float] = 10) -> tuple[np.ndarray, np.ndarray]:
        """Return ``(counts, bin_edges)`` — a plain numpy histogram of the samples."""
        counts, edges = np.histogram(self.values, bins=bins)
        return counts, edges

    def within_band(self, lo: float, hi: float) -> float:
        """Return the fraction of samples in the closed literature band ``[lo, hi]``.

        This is how a gate scores a DISTRIBUTION against a literature range — the acceptance is a
        fraction-in-band, never an equality on a single value or a fitted parameter.
        """
        if not (_is_finite(lo) and _is_finite(hi)):
            raise ValueError("band bounds must be finite")
        if lo > hi:
            raise ValueError("band lower bound must not exceed the upper bound")
        if self.values.size == 0:
            return 0.0
        inside = np.count_nonzero((self.values >= lo) & (self.values <= hi))
        return float(inside) / float(self.values.size)


@dataclass(frozen=True, slots=True)
class EnsembleSummary:
    """Aggregated steady-state statistics over N realizations of ONE CellState.

    Args:
        cell_state_tag: the shared CellState tag of every realization.
        n_realizations: the number of realizations.
        seeds: the per-realization seeds, in realization order.
        scalar_distributions: one :class:`ObservableDistribution` per scalar observable, whose values are
            that scalar ACROSS the N realizations (per-realization spread).
        pooled_distributions: one :class:`ObservableDistribution` per distribution observable, whose
            values are all realizations' samples POOLED.
        results: the underlying realization results, in order.
    """

    cell_state_tag: tuple[str, ...]
    n_realizations: int
    seeds: tuple[int, ...]
    scalar_distributions: Mapping[str, ObservableDistribution]
    pooled_distributions: Mapping[str, ObservableDistribution]
    results: tuple[RealizationResult, ...]

    def __post_init__(self) -> None:
        object.__setattr__(self, "cell_state_tag", _coerce_tag(self.cell_state_tag))
        object.__setattr__(
            self, "scalar_distributions", MappingProxyType(dict(self.scalar_distributions))
        )
        object.__setattr__(
            self, "pooled_distributions", MappingProxyType(dict(self.pooled_distributions))
        )
        object.__setattr__(self, "seeds", tuple(self.seeds))
        object.__setattr__(self, "results", tuple(self.results))


@dataclass(slots=True)
class EnsembleRunner:
    """Runs N reproducible, independent realizations of one CellState and aggregates a summary.

    The runner owns ZERO physics.  It builds an independent per-realization
    :class:`numpy.random.Generator` and hands it to the caller-supplied ``simulate`` (mirroring
    :class:`~aleph.engine.transaction.CellTransaction`'s caller-owned ``solve``), then collects and
    aggregates the returned :class:`RealizationResult` objects.

    Args:
        cell_state: the :class:`CellState` every realization runs under.
        base_seed: the ensemble base seed; nonnegative int.  Same base_seed → identical seed list.
        n_realizations: the number of realizations; a positive int.
    """

    cell_state: CellState
    base_seed: int
    n_realizations: int

    def __post_init__(self) -> None:
        if not isinstance(self.cell_state, CellState):
            raise TypeError("EnsembleRunner.cell_state must be a CellState")
        self.base_seed = _require_nonneg_int(self.base_seed, what="base_seed")
        self.n_realizations = _require_nonneg_int(self.n_realizations, what="n_realizations")
        if self.n_realizations < 1:
            raise ValueError("n_realizations must be at least 1")

    def _seed_sequences(self) -> list[np.random.SeedSequence]:
        """Spawn ``n_realizations`` independent child SeedSequences from ``base_seed``.

        ``SeedSequence.spawn`` is the numpy-canonical way to derive many decorrelated, reproducible child
        streams from one parent seed: the children are independent (distinct spawn keys) yet fully
        reproducible (a pure function of ``base_seed`` + index), so a rerun with the same ``base_seed``
        replays every realization identically.
        """
        return list(np.random.SeedSequence(self.base_seed).spawn(self.n_realizations))

    def seeds(self) -> tuple[int, ...]:
        """Return the reproducible, independent per-realization integer seeds.

        Each seed is ``int(child.generate_state(1)[0])`` of a spawned child SeedSequence — reproducible in
        ``base_seed`` AND mutually independent.  Same ``base_seed`` yields an identical tuple; a different
        ``base_seed`` yields a (with overwhelming probability) different one.
        """
        return tuple(int(child.generate_state(1)[0]) for child in self._seed_sequences())

    def run(
        self,
        simulate: Callable[[np.random.Generator, CellState, int], RealizationResult],
    ) -> EnsembleSummary:
        """Run every realization through ``simulate`` and aggregate an :class:`EnsembleSummary`.

        Args:
            simulate: caller-owned stochastic runtime,
                ``simulate(rng, cell_state, realization_index) -> RealizationResult``.  It receives an
                independent :class:`numpy.random.Generator` per realization and returns that realization's
                steady-state observables.  The runner does NO physics itself.

        Returns:
            The aggregated :class:`EnsembleSummary`: per-realization scalar distributions and pooled
            sample distributions across the N realizations.

        Raises:
            TypeError: if ``simulate`` is not callable or returns a non-:class:`RealizationResult`.
            ValueError: if a returned result's CellState tag or realization index does not match.
        """
        if not callable(simulate):
            raise TypeError("simulate must be callable(rng, cell_state, realization_index)")
        tag = tuple(self.cell_state.axes)
        sequences = self._seed_sequences()
        seeds = tuple(int(child.generate_state(1)[0]) for child in sequences)

        results: list[RealizationResult] = []
        for index, child in enumerate(sequences):
            rng = np.random.default_rng(child)
            result = simulate(rng, self.cell_state, index)
            if not isinstance(result, RealizationResult):
                raise TypeError("simulate must return a RealizationResult")
            if tuple(result.cell_state_tag) != tag:
                raise ValueError(
                    f"realization {index} returned CellState tag {tuple(result.cell_state_tag)} "
                    f"!= ensemble tag {tag}"
                )
            if result.realization_index != index:
                raise ValueError(
                    f"realization {index} returned realization_index {result.realization_index}"
                )
            results.append(result)

        scalar_names = sorted({name for r in results for name in r.scalars})
        scalar_distributions = {
            name: ObservableDistribution(
                name=name,
                values=np.array([r.scalars[name] for r in results if name in r.scalars], dtype=float),
            )
            for name in scalar_names
        }
        dist_names = sorted({name for r in results for name in r.distributions})
        pooled_distributions = {}
        for name in dist_names:
            chunks = [r.distributions[name] for r in results if name in r.distributions]
            pooled = np.concatenate(chunks) if chunks else np.empty(0, dtype=float)
            pooled_distributions[name] = ObservableDistribution(name=name, values=pooled)

        return EnsembleSummary(
            cell_state_tag=tag,
            n_realizations=self.n_realizations,
            seeds=seeds,
            scalar_distributions=scalar_distributions,
            pooled_distributions=pooled_distributions,
            results=tuple(results),
        )


def mix_ensembles(
    weighted: Sequence[tuple[EnsembleSummary, float]],
) -> Mapping[str, ObservableDistribution]:
    """Form the state-conditioned experimental MIXTURE of several ensembles at observed proportions.

    The experimental distribution mixes state-conditioned ensembles at their observed proportions
    (spec §7).  Each ``(summary, proportion)`` pair contributes its observable samples weighted by its
    proportion.

    Mixing is **deterministic proportional allocation** (no hidden RNG): for each mixed observable, each
    ensemble contributes ``round(proportion * total_pool)`` of its own samples (evenly strided so the
    sub-sample spans the ensemble, not just its head), and the strided pieces are concatenated.  This
    reproduces the mixture composition exactly given fixed inputs, with no seed to thread.

    Only observables present in ALL summaries are mixed (an observable missing from any state has no
    defined mixture proportion), pooling each summary's scalar_distributions and pooled_distributions.

    Args:
        weighted: a sequence of ``(EnsembleSummary, proportion)`` pairs.  Every proportion must be
            ``> 0`` and they must sum to 1 within ``1e-9``.

    Returns:
        A mapping ``observable_name -> ObservableDistribution`` of the mixture, for every observable
        shared across all summaries.

    Raises:
        ValueError: if there are no summaries, a proportion is non-positive, or the proportions do not
            sum to 1 within tolerance.
    """
    pairs = list(weighted)
    if not pairs:
        raise ValueError("mix_ensembles needs at least one (summary, proportion) pair")
    proportions: list[float] = []
    for summary, proportion in pairs:
        if not isinstance(summary, EnsembleSummary):
            raise TypeError("mix_ensembles expects (EnsembleSummary, proportion) pairs")
        if not _is_finite(proportion) or proportion <= 0.0:
            raise ValueError("every mixture proportion must be a finite positive number")
        proportions.append(float(proportion))
    if abs(sum(proportions) - 1.0) > 1e-9:
        raise ValueError(f"mixture proportions must sum to 1 (got {sum(proportions)!r})")

    def _observables(summary: EnsembleSummary) -> dict[str, ObservableDistribution]:
        merged = dict(summary.scalar_distributions)
        merged.update(summary.pooled_distributions)
        return merged

    per_summary = [_observables(summary) for summary, _ in pairs]
    shared = set(per_summary[0])
    for observables in per_summary[1:]:
        shared &= set(observables)

    mixture: dict[str, ObservableDistribution] = {}
    for name in sorted(shared):
        available = [observables[name].values for observables in per_summary]
        total = sum(values.size for values in available)
        pieces: list[np.ndarray] = []
        for proportion, values in zip(proportions, available, strict=True):
            take = int(round(proportion * total))
            if take <= 0 or values.size == 0:
                continue
            take = min(take, values.size)
            # Evenly strided pick so the sub-sample spans the ensemble deterministically.
            idx = np.linspace(0, values.size - 1, take).round().astype(int)
            pieces.append(values[idx])
        pooled = np.concatenate(pieces) if pieces else np.empty(0, dtype=float)
        mixture[name] = ObservableDistribution(name=name, values=pooled)
    return MappingProxyType(mixture)
