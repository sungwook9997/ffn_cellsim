"""Declaring a sweep before running it.

A sweep that is specified in the same breath as it is launched tends to grow an extra axis every
time a result looks interesting, and by the end nobody can say which comparisons were planned and
which were found. So a campaign here is a value: it is written down, hashed, and the hash travels
with every result it produced.

Two things this module deliberately does not do. It does not decide how many points are enough —
that depends on identifiability, which is the inference layer's business. And it does not treat a
grid as the default shape; a grid over five axes is thirty-two thousand points and almost always
the wrong instrument. Latin hypercube and explicit point lists are first-class here for that reason.
"""

from __future__ import annotations

import hashlib
import itertools
import json
from collections.abc import Iterator, Mapping, Sequence
from dataclasses import dataclass, field
from enum import StrEnum
from typing import Any

import numpy as np
from numpy.typing import NDArray


class Sampling(StrEnum):
    """How points are drawn from the declared axes."""

    GRID = "grid"
    """Full Cartesian product. Honest for two axes, exponential past three."""

    LATIN_HYPERCUBE = "latin_hypercube"
    """Stratified: n points, each axis's range covered once. Scales with n, not with dimension."""

    EXPLICIT = "explicit"
    """A hand-written list of points. Use when the interesting points are already known."""


class Spacing(StrEnum):
    LINEAR = "linear"
    LOG = "log"


@dataclass(frozen=True, slots=True)
class Axis:
    """One swept parameter, with its units attached because a bare number is not a parameter."""

    name: str
    low: float
    high: float
    unit: str
    spacing: Spacing = Spacing.LINEAR
    count: int = 5

    def __post_init__(self) -> None:
        if self.high <= self.low:
            raise ValueError(f"axis {self.name!r}: high ({self.high}) must exceed low ({self.low})")
        if self.count < 1:
            raise ValueError(f"axis {self.name!r}: count must be at least 1, got {self.count}")
        if self.spacing is Spacing.LOG and self.low <= 0.0:
            raise ValueError(
                f"axis {self.name!r}: log spacing needs a positive lower bound, got {self.low}. "
                "A log axis through zero is not a rescaling, it is a different question."
            )

    def values(self) -> NDArray[np.float64]:
        if self.spacing is Spacing.LOG:
            return np.logspace(np.log10(self.low), np.log10(self.high), self.count)
        return np.linspace(self.low, self.high, self.count)


@dataclass(frozen=True, slots=True)
class SweepSpec:
    """A predeclared campaign. Hash it, run it, and quote the hash next to every number."""

    name: str
    axes: tuple[Axis, ...]
    sampling: Sampling = Sampling.GRID
    n_points: int = 0
    """Only used by LATIN_HYPERCUBE. Ignored for GRID (the axes decide) and EXPLICIT."""

    explicit_points: tuple[Mapping[str, float], ...] = ()
    seed: int = 0
    """Base seed. See `aleph.campaign.rng` — repeats share a stream, points do not."""

    repeats: int = 1
    """Independent stochastic realizations per point, for a variance estimate."""

    common_random_numbers: bool = True
    """Reuse the same noise across points so that differences between points are not swamped by
    sampling noise. This is a large variance reduction for comparisons and a subtle trap for
    absolute values, so it is recorded in the spec rather than chosen at the call site."""

    notes: str = ""
    metadata: Mapping[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if not self.name:
            raise ValueError("a campaign needs a name; unnamed results cannot be cited")
        names = [axis.name for axis in self.axes]
        if len(set(names)) != len(names):
            raise ValueError(f"duplicate axis names: {names}")
        if self.sampling is Sampling.LATIN_HYPERCUBE and self.n_points < 1:
            raise ValueError("latin hypercube sampling needs n_points >= 1")
        if self.sampling is Sampling.EXPLICIT and not self.explicit_points:
            raise ValueError("explicit sampling needs at least one point")
        if self.repeats < 1:
            raise ValueError("repeats must be at least 1")

    def points(self) -> list[dict[str, float]]:
        """Every parameter point this campaign will visit, in a deterministic order."""
        if self.sampling is Sampling.EXPLICIT:
            return [dict(point) for point in self.explicit_points]

        if self.sampling is Sampling.GRID:
            per_axis = [axis.values() for axis in self.axes]
            return [
                {axis.name: float(value) for axis, value in zip(self.axes, combination, strict=True)}
                for combination in itertools.product(*per_axis)
            ]

        # Latin hypercube: one stratum per point per axis, permuted independently.
        rng = np.random.default_rng(self.seed)
        n = self.n_points
        out: list[dict[str, float]] = [{} for _ in range(n)]
        for axis in self.axes:
            strata = (rng.permutation(n) + rng.random(n)) / n
            if axis.spacing is Spacing.LOG:
                lo, hi = np.log10(axis.low), np.log10(axis.high)
                drawn = 10.0 ** (lo + strata * (hi - lo))
            else:
                drawn = axis.low + strata * (axis.high - axis.low)
            for index in range(n):
                out[index][axis.name] = float(drawn[index])
        return out

    def total_runs(self) -> int:
        return len(self.points()) * self.repeats

    def canonical(self) -> dict[str, Any]:
        """A stable dict for hashing. Float formatting is pinned so the hash is reproducible."""
        return {
            "name": self.name,
            "sampling": str(self.sampling),
            "n_points": self.n_points,
            "repeats": self.repeats,
            "seed": self.seed,
            "common_random_numbers": self.common_random_numbers,
            "axes": [
                {
                    "name": axis.name,
                    "low": f"{axis.low:.17g}",
                    "high": f"{axis.high:.17g}",
                    "unit": axis.unit,
                    "spacing": str(axis.spacing),
                    "count": axis.count,
                }
                for axis in self.axes
            ],
            "explicit_points": [
                {key: f"{value:.17g}" for key, value in sorted(point.items())}
                for point in self.explicit_points
            ],
            "metadata": dict(sorted(self.metadata.items())),
        }

    def spec_hash(self) -> str:
        blob = json.dumps(self.canonical(), sort_keys=True, separators=(",", ":"))
        return hashlib.sha256(blob.encode("utf-8")).hexdigest()

    def iter_runs(self) -> Iterator[tuple[int, int, dict[str, float]]]:
        """Yield (point_index, repeat_index, parameters) for every run in the campaign."""
        for point_index, parameters in enumerate(self.points()):
            for repeat_index in range(self.repeats):
                yield point_index, repeat_index, parameters


def axes_from_pairs(pairs: Sequence[tuple[str, float, float, str]]) -> tuple[Axis, ...]:
    """Convenience for the common case of linear axes: (name, low, high, unit)."""
    return tuple(Axis(name=n, low=lo, high=hi, unit=u) for n, lo, hi, u in pairs)
