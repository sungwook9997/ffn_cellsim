"""Top-level label-blind emergence detector — the honest replacement for ``bundle_count``.

``detect()`` reads ONLY the FF fiber-array geometry contract (``pos`` / ``fiber_offsets``) and returns a
structured verdict: the global nematic order S, its effect size against the pre-registered finite-N
isotropic null, and the list of spatially LOCALIZED condensed bundles. No argument names or reads a
region / type / construction label — that is the I5 anti-coupling firewall, enforced structurally by
``tests/ac/emergence/test_detector_contract.py`` (the signature carries geometry + statistical config
only).

The emergence-relevant scalars this exposes (all measured, none handed in):
  * ``n_emerged_bundles`` — how many ordered+dense patches CONDENSED (the honest ``bundle_count``);
  * ``condensed_fraction`` — fraction of fibers inside any recovered bundle;
  * ``is_ordered_global`` — whether the whole network beat the isotropic null at the pre-registered z.

Falsifiability (hard-truth #1 / §5): a flat result (``S`` inside the null band, no bundle recovered) is a
genuine NEGATIVE — surfaced to PI as a FINDING with a SEEDED-labeled scaffold fallback, NEVER re-tuned to
force a positive.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Protocol

import numpy as np
import numpy.typing as npt

from aleph.components.emergence.condensation import (
    NEIGHBORHOOD_SCALE,
    Bundle,
    compute_local_fields,
    localize_bundles,
)
from aleph.components.emergence.nematic import fiber_axes, nematic_order_and_director
from aleph.components.emergence.null_model import (
    EFFECT_SIZE_Z_CRIT,
    NullBand,
    effect_size,
    nematic_null_band,
)

__all__ = ["FiberArrayLike", "EmergenceReport", "EmergenceDetector", "detect"]


class FiberArrayLike(Protocol):
    """Structural type for anything carrying the FF fiber-array contract (e.g. ``FiberNetwork``)."""

    pos: npt.NDArray[np.float64]
    fiber_offsets: npt.NDArray[np.int64]


@dataclass(frozen=True, slots=True)
class EmergenceReport:
    """The detector's structured, label-free verdict.

    Attributes:
        n_fibers: fiber count.
        s_global: global nematic order S = lambda_max(Q) (the validated order parameter).
        director: (3,) global nematic director.
        null_band: the finite-N isotropic :class:`NullBand` S was scored against.
        z_global: effect size (S - mu_null) / sigma_null.
        is_ordered_global: whether ``z_global`` beat the pre-registered threshold.
        radius: neighborhood radius used for the local condensation fields.
        bundles: recovered localized bundles, largest-first.
    """

    n_fibers: int
    s_global: float
    director: npt.NDArray[np.float64]
    null_band: NullBand
    z_global: float
    is_ordered_global: bool
    radius: float
    bundles: list[Bundle]

    @property
    def n_emerged_bundles(self) -> int:
        """The honest bundle count: number of ordered+dense patches that CONDENSED (not ``n_fibers``)."""
        return len(self.bundles)

    @property
    def condensed_fraction(self) -> float:
        """Fraction of fibers inside any recovered bundle."""
        if self.n_fibers == 0:
            return 0.0
        n_in = int(sum(b.size for b in self.bundles))
        return n_in / self.n_fibers


@dataclass(frozen=True, slots=True)
class EmergenceDetector:
    """Label-blind nematic / bundle-condensation detector.

    Attributes:
        z_crit: pre-registered effect-size threshold (default :data:`EFFECT_SIZE_Z_CRIT`).
        scale: neighborhood multiplier for the local fields (default :data:`NEIGHBORHOOD_SCALE`).
        null_mc: Monte-Carlo draws for the global null band.
        null_seed: RNG seed for the global null band (reproducible).
    """

    z_crit: float = EFFECT_SIZE_Z_CRIT
    scale: float = NEIGHBORHOOD_SCALE
    null_mc: int = 4000
    null_seed: int = 20260716

    def detect(
        self, pos: npt.NDArray[np.float64], fiber_offsets: npt.NDArray[np.int64]
    ) -> EmergenceReport:
        """Run the full detector on a fiber-array (geometry only — NO labels).

        Args:
            pos: (N, 3) node positions.
            fiber_offsets: (F+1,) fiber offsets.

        Returns:
            The :class:`EmergenceReport` verdict.
        """
        axes = fiber_axes(pos, fiber_offsets)
        n_fibers = axes.shape[0]
        s_global, director = nematic_order_and_director(axes)
        band = nematic_null_band(n_fibers, n_mc=self.null_mc, seed=self.null_seed)
        z = effect_size(s_global, band)
        fields = compute_local_fields(pos, fiber_offsets, scale=self.scale)
        bundles = localize_bundles(fields, z_crit=self.z_crit)
        return EmergenceReport(
            n_fibers=n_fibers,
            s_global=s_global,
            director=director,
            null_band=band,
            z_global=z,
            is_ordered_global=bool(z > self.z_crit),
            radius=fields.radius,
            bundles=bundles,
        )

    def detect_network(self, net: FiberArrayLike) -> EmergenceReport:
        """Convenience: run :meth:`detect` on anything exposing ``pos`` / ``fiber_offsets``."""
        return self.detect(net.pos, net.fiber_offsets)


def detect(
    pos: npt.NDArray[np.float64],
    fiber_offsets: npt.NDArray[np.int64],
    z_crit: float = EFFECT_SIZE_Z_CRIT,
    scale: float = NEIGHBORHOOD_SCALE,
) -> EmergenceReport:
    """Module-level convenience wrapper around :meth:`EmergenceDetector.detect` (geometry only)."""
    return EmergenceDetector(z_crit=z_crit, scale=scale).detect(pos, fiber_offsets)
