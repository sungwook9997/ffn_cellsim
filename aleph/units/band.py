"""Acceptance bands that carry their source, and refuse rather than clamp. `ALEPH-PORT-4001`.

Two commitments, and they are the whole content of this module. The code under them is ordinary.

**A band cannot exist without a citation.** :class:`LiteratureBand` takes ``citation`` as a required
constructor argument and refuses an empty one. A bound with no source is somebody's memory of a
bound, and it cannot be argued with: when the guard fires, the reader has no way to decide whether
the band is wrong or the run is. The citation is the difference between a guard that ends an
argument and a guard that starts one.

**A guard raises and returns the value untouched. It never clamps.** A quantity that leaves the band
its sources support is a finding. Clamping it to the boundary destroys the finding and leaves a run
that looks healthy — which is how a number gets produced, quoted, and only later found to have been
outside what supported it. This repository's own record of that is `STATE.md` (c) and the withdrawn
headlines in `ff/ENGINE.md`.

What this module is not
-----------------------
**It is not a unit system.** ``unit`` here is a string carried into the error message so the reader
knows what they are looking at. It is *not* checked, and a value passed in nm against a band in µm
guards green. This is stated rather than hidden because a reader who assumes otherwise will be
wrong in the direction that matters. Dimensional enforcement needs a quantity type, which
`ALEPH-PORT-4001` §1 deliberately excludes; this repository's unit system is `ff/units.py` and a
second one would be a second answer to what the numbers mean.

**It does not belong in the inner loop.** Bands are checked when a value is declared, when a run's
inputs are assembled, and when a result is read back between accepted steps. A guard inside the
mechanical solve would be authoritative per-step host state, which the charter forbids. Nothing
enforces that yet — see `ALEPH-PORT-4001` §14.3.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

__all__ = [
    "LiteratureBand",
    "OutOfBandError",
    "band_guard",
]


class OutOfBandError(ValueError):
    """A quantity fell outside a cited acceptance band, or was not a finite number at all."""


def _bracket(inclusive: bool) -> tuple[str, str]:
    return ("[", "]") if inclusive else ("(", ")")


@dataclass(frozen=True, slots=True)
class LiteratureBand:
    """A named acceptance band with its unit and its source attached.

    Declare the band once and guard by name, so the bound and its citation are written in exactly
    one place::

        CORTICAL_TENSION = LiteratureBand(
            lo=1e-2, hi=1.0, unit="pN/um",
            citation="KU-3.1 -- cortical tension, adherent cells",
        )
        CORTICAL_TENSION.guard("sigma_cortex", measured)

    Attributes:
        lo: Lower bound. Must be finite and ``<= hi``.
        hi: Upper bound. Must be finite and ``>= lo``.
        unit: Unit symbol, carried into the error message. **Not enforced** — see the module
            docstring.
        citation: Where the bound comes from. Refused if empty or whitespace.
        inclusive: Whether the bounds themselves are inside the band.

    Raises:
        ValueError: If the band is malformed, or has no citation.
    """

    lo: float
    hi: float
    unit: str
    citation: str
    inclusive: bool = True

    def __post_init__(self) -> None:
        lo, hi = float(self.lo), float(self.hi)
        # A NaN bound compares false against everything, so a band carrying one looks like a check
        # and is not. One infinite bound is legitimate — a strictly-positive quantity is `(0, inf)`
        # and inventing a large finite ceiling for it would be a number with no source, which is
        # the thing this whole module exists to refuse. Two infinite bounds are vacuous.
        if math.isnan(lo) or math.isnan(hi):
            raise ValueError(f"a band bound cannot be NaN, got [{self.lo!r}, {self.hi!r}]")
        if math.isinf(lo) and math.isinf(hi):
            raise ValueError(
                "a band infinite on both sides admits every finite value, so it is not a check. "
                "Use a one-sided band, or none."
            )
        if lo > hi:
            raise ValueError(f"malformed band: lo={lo!r} is above hi={hi!r}")
        if not isinstance(self.citation, str) or not self.citation.strip():
            raise ValueError(
                f"the band [{lo!r}, {hi!r}] {self.unit!r} has no citation. A bound without a "
                "source is somebody's memory of a bound: when this guard fires, nobody can tell "
                "whether the band is wrong or the run is."
            )
        if not isinstance(self.unit, str) or not self.unit.strip():
            raise ValueError(f"the band [{lo!r}, {hi!r}] needs a unit symbol")

    def contains(self, value: float) -> bool:
        """Test membership without raising.

        A non-finite value — ``NaN`` or either infinity — is **outside every band**. ``NaN`` in
        particular is not a small error, and the two obvious spellings of the interval test
        disagree about it: ``lo <= x <= hi`` is ``False`` while ``not (x < lo or x > hi)`` is
        ``True``. This is the first, and it is checked explicitly rather than relied on.
        """
        numeric = float(value)
        if not math.isfinite(numeric):
            return False
        if self.inclusive:
            return self.lo <= numeric <= self.hi
        return self.lo < numeric < self.hi

    def guard(self, name: str, value: float) -> float:
        """Check one scalar. Returns it unchanged, or raises.

        Args:
            name: What is being guarded, for the message. Not a lookup key.
            value: The quantity.

        Returns:
            ``value``, unchanged and unclamped.

        Raises:
            OutOfBandError: If the value is outside the band or not finite.
        """
        if not self.contains(value):
            raise OutOfBandError(f"{name} = {float(value)!r} {self.unit} is outside {self}")
        return value

    def guard_all(self, name: str, values: np.ndarray) -> np.ndarray:
        """Check every element of an array. Returns the **same object**, unchanged.

        Returning the input object rather than a copy is deliberate: a guard that quietly hands back
        a different array is a guard that can be mistaken for a transform.

        Raises:
            OutOfBandError: If any element is outside the band or not finite.
            TypeError, ValueError: If the input cannot be read as a float array. It is never
                coerced to ``NaN``, because that would turn a wrong dtype into an out-of-band
                report and hide which of the two happened.
        """
        array = np.asarray(values, dtype=float)
        if self.inclusive:
            inside = (array >= self.lo) & (array <= self.hi)
        else:
            inside = (array > self.lo) & (array < self.hi)
        inside &= np.isfinite(array)
        if not bool(inside.all()):
            outside = np.flatnonzero(~inside)
            worst = int(outside[int(np.argmax(np.abs(array.ravel()[outside])))])
            shown = ", ".join(str(int(i)) for i in outside[:8])
            more = f", … ({outside.size} total)" if outside.size > 8 else ""
            raise OutOfBandError(
                f"{name}: {outside.size} of {array.size} elements are outside {self} "
                f"(worst: index {worst} = {array.ravel()[worst]!r}; indices {shown}{more})"
            )
        return values

    def __str__(self) -> str:
        left, right = _bracket(self.inclusive)
        return f"{left}{self.lo:g}, {self.hi:g}{right} {self.unit} ({self.citation})"


def band_guard(
    name: str,
    value: float,
    lo: float,
    hi: float,
    unit: str,
    citation: str,
    *,
    inclusive: bool = True,
) -> float:
    """Guard a scalar against a band declared inline.

    For a band checked from exactly one place. If the same bound is checked from two, declare a
    :class:`LiteratureBand` instead — two inline copies of a bound are two bounds as soon as one of
    them is edited.
    """
    return LiteratureBand(lo=lo, hi=hi, unit=unit, citation=citation, inclusive=inclusive).guard(
        name, value
    )
