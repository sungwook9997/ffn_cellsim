"""Reproducible, independent random streams for a campaign.

The rule this module exists to enforce: **a run's randomness is a function of its coordinates, never
of the order in which runs happened to execute.** Draw seeds from a single generator in a loop and
the campaign becomes unresumable and unparallelizable — restart it halfway and the second half gets
different noise than it would have, and no artifact records that it did.

`numpy.random.SeedSequence.spawn_key` solves this directly: a stream is addressed by its coordinates,
so run (point 7, repeat 2) gets the same bits whether it runs first, last, or on another machine.

On common random numbers. When comparing points, using the *same* noise at every point removes the
sampling noise from the difference, which is usually a large variance reduction and occasionally a
trap — a single unlucky noise realization biases every point the same way, so an absolute value
computed under CRN carries a correlated error that the per-point spread does not reveal. Both modes
are therefore explicit and recorded in the spec, and `describe()` states which one produced a result.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True, slots=True)
class StreamAddress:
    """Where a random stream sits in the campaign, independent of execution order."""

    point_index: int
    repeat_index: int
    purpose: int = 0
    """Distinguishes independent uses within one run — thermal noise vs initial condition vs
    measurement noise should not share a stream, or a change to one silently perturbs the others."""


class CampaignRNG:
    """Addressable random streams for one campaign."""

    def __init__(self, base_seed: int, *, common_random_numbers: bool) -> None:
        self._base_seed = int(base_seed)
        self._crn = bool(common_random_numbers)

    @property
    def common_random_numbers(self) -> bool:
        return self._crn

    def seed_sequence(self, address: StreamAddress) -> np.random.SeedSequence:
        if self._crn:
            # Under CRN the point index is deliberately dropped from the address, so every point
            # in a given repeat sees identical noise and differences between points are signal.
            key = (address.repeat_index, address.purpose)
        else:
            key = (address.point_index, address.repeat_index, address.purpose)
        return np.random.SeedSequence(self._base_seed, spawn_key=key)

    def generator(self, address: StreamAddress) -> np.random.Generator:
        return np.random.default_rng(self.seed_sequence(address))

    def describe(self, address: StreamAddress) -> dict[str, object]:
        """What goes into the run record, so the noise is reconstructable from the artifact alone."""
        sequence = self.seed_sequence(address)
        return {
            "base_seed": self._base_seed,
            "common_random_numbers": self._crn,
            "spawn_key": list(sequence.spawn_key),
            "entropy": sequence.entropy,
            "point_index": address.point_index,
            "repeat_index": address.repeat_index,
            "purpose": address.purpose,
        }


def streams_are_independent(
    rng: CampaignRNG, a: StreamAddress, b: StreamAddress, *, draws: int = 4096
) -> tuple[bool, float]:
    """A cheap smoke test that two addresses really do give different bits.

    Returns (independent, |correlation|). This is a smoke test, not a randomness test — it catches
    the realistic failure (two addresses collapsing to the same stream, giving correlation exactly
    1.0), and says nothing about the quality of the underlying generator, which is not ours to
    judge.
    """
    left = rng.generator(a).standard_normal(draws)
    right = rng.generator(b).standard_normal(draws)
    if np.allclose(left, right):
        return False, 1.0
    correlation = float(abs(np.corrcoef(left, right)[0, 1]))
    return correlation < 5.0 / np.sqrt(draws), correlation
