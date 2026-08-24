"""Scalar channels and colour bars: what a rendered colour is allowed to be read as.

Manuscript §11: *"Color bars record units, actual extrema, and display clamps; visualization never
turns a clipped field into a quantitative claim."*

The hazard is specific and it is not hypothetical. A field is rendered with a display range chosen so
the interesting structure is visible. Values outside that range saturate to the end colours. Somebody
later reads a number off the picture — "the peak stress is about 40 pN" — and that number is the
*clamp*, not the peak. The picture cannot distinguish 40 from 4000, and neither can anyone reading it.

So :class:`ColourBar` carries three things that are usually conflated into one "range":

* ``units`` — without which the number is not a physical quantity at all;
* the **actual** extrema, measured over the active population of the buffer being rendered;
* the **display** clamps, which are a presentation choice.

And the reads are typed. A normalised coordinate in the saturated band refuses. An extremum whose end
is clamped refuses. The refusal is :class:`ClampedChannelError` and it names which end saturated and
by how much, because the useful next action is almost always "widen the display range and re-render",
which needs the actual extremum to do.

Note the asymmetry that keeps this usable: the *interior* of a clamped bar is still exact. Linear
interpolation between two unclamped display bounds recovers the value it was built from. Refusing the
whole bar because one end saturated would make the guard so blunt that someone would turn it off.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from enum import StrEnum, unique
from typing import Any

__all__ = [
    "ClampedChannelError",
    "ColourBar",
    "ScalarChannel",
    "ScaleKind",
    "colour_bar_for",
]


class ClampedChannelError(ValueError):
    """A quantitative read was attempted on a saturated part of a clamped channel."""


@unique
class ScaleKind(StrEnum):
    """How display coordinates map to values."""

    LINEAR = "linear"
    LOG10 = "log10"

    @property
    def requires_positive(self) -> bool:
        """Whether the mapping is undefined at or below zero."""
        return self is ScaleKind.LOG10


@dataclass(frozen=True, slots=True)
class ScalarChannel:
    """The one scalar field a frame is showing, named and dimensioned.

    Attributes:
        name: Channel identity, for example ``"force_magnitude"``.
        units: Physical units, for example ``"pN"``. Required.
        owner: Which authority the channel was derived from.
        derivation: What reduction produced it. Two channels with the same name and different
            reductions are different observables, and the frame has to say which one it is.
    """

    name: str
    units: str
    owner: str
    derivation: str

    def __post_init__(self) -> None:
        for field_name in ("name", "units", "owner", "derivation"):
            if not str(getattr(self, field_name)).strip():
                raise ValueError(
                    f"ScalarChannel.{field_name} is required; an unlabelled channel renders a "
                    "picture rather than a measurement"
                )

    def as_json_obj(self) -> dict[str, Any]:
        """Return the JSON-able channel block."""
        return {
            "name": self.name,
            "units": self.units,
            "owner": self.owner,
            "derivation": self.derivation,
        }


@dataclass(frozen=True, slots=True)
class ColourBar:
    """A display mapping, with the clamp it applied recorded alongside what it clamped.

    Attributes:
        channel: The channel this bar renders.
        actual_min: Smallest value actually present in the active population.
        actual_max: Largest value actually present.
        display_min: Lower end of the display range.
        display_max: Upper end of the display range.
        scale: How display coordinates map to values.
        palette: Name of the colour map, recorded so a re-render is reproducible.
    """

    channel: ScalarChannel
    actual_min: float
    actual_max: float
    display_min: float
    display_max: float
    scale: ScaleKind = ScaleKind.LINEAR
    palette: str = "viridis"

    def __post_init__(self) -> None:
        if not isinstance(self.channel, ScalarChannel):
            raise TypeError("ColourBar.channel must be a ScalarChannel")
        for field_name in ("actual_min", "actual_max", "display_min", "display_max"):
            value = getattr(self, field_name)
            if isinstance(value, bool) or not isinstance(value, (int, float)):
                raise TypeError(f"ColourBar.{field_name} must be a real number; got {value!r}")
            if not math.isfinite(float(value)):
                raise ValueError(
                    f"ColourBar.{field_name} is {value!r}. A non-finite bound means the field itself "
                    "carries a NaN or an infinity, which is a physics finding and not a display "
                    "setting; it must not be smoothed over by a renderer"
                )
            object.__setattr__(self, field_name, float(value))
        if self.actual_max < self.actual_min:
            raise ValueError(
                f"actual extrema are inverted: [{self.actual_min}, {self.actual_max}]"
            )
        if not (self.display_max > self.display_min):
            raise ValueError(
                f"display range must be non-degenerate; got [{self.display_min}, {self.display_max}]"
            )
        if self.scale.requires_positive and self.display_min <= 0.0:
            raise ValueError(
                f"a {self.scale.value} colour bar needs a positive display_min; got "
                f"{self.display_min}"
            )
        if not isinstance(self.palette, str) or not self.palette.strip():
            raise ValueError("ColourBar.palette is required so a render can be reproduced")

    # -- what the clamp did ------------------------------------------------------------------
    @property
    def clamped_low(self) -> bool:
        """Whether values below the display range exist and were saturated."""
        return self.actual_min < self.display_min

    @property
    def clamped_high(self) -> bool:
        """Whether values above the display range exist and were saturated."""
        return self.actual_max > self.display_max

    @property
    def is_clamped(self) -> bool:
        """Whether either end saturated."""
        return self.clamped_low or self.clamped_high

    @property
    def clamp_excess_low(self) -> float:
        """How far below ``display_min`` the field actually went. Zero when not clamped low."""
        return max(0.0, self.display_min - self.actual_min)

    @property
    def clamp_excess_high(self) -> float:
        """How far above ``display_max`` the field actually went. Zero when not clamped high."""
        return max(0.0, self.actual_max - self.display_max)

    def clamp_report(self) -> str:
        """One line naming which ends saturated and by how much."""
        if not self.is_clamped:
            return "not clamped: the display range covers the measured extrema"
        parts: list[str] = []
        if self.clamped_low:
            parts.append(
                f"low end saturated: actual_min {self.actual_min:.6g} is "
                f"{self.clamp_excess_low:.6g} {self.channel.units} below display_min "
                f"{self.display_min:.6g}"
            )
        if self.clamped_high:
            parts.append(
                f"high end saturated: actual_max {self.actual_max:.6g} is "
                f"{self.clamp_excess_high:.6g} {self.channel.units} above display_max "
                f"{self.display_max:.6g}"
            )
        return "; ".join(parts)

    # -- the reads ---------------------------------------------------------------------------
    def normalised_for(self, value: float) -> float:
        """Display coordinate in ``[0, 1]`` for ``value``, saturating at the clamps.

        This is the renderer's own direction and it is always allowed: clamping is what a renderer
        does. It is the *inverse* that has to refuse.
        """
        v = float(value)
        if self.scale is ScaleKind.LOG10:
            lo, hi = math.log10(self.display_min), math.log10(self.display_max)
            v = math.log10(v) if v > 0.0 else lo
        else:
            lo, hi = self.display_min, self.display_max
        return min(1.0, max(0.0, (v - lo) / (hi - lo)))

    def value_at_normalised(self, t: float, *, saturation_band: float = 0.0) -> float:
        """Value corresponding to display coordinate ``t``, refusing in the saturated band.

        Args:
            t: Display coordinate in ``[0, 1]``.
            saturation_band: Width of the band at each end treated as saturated. ``0.0`` means only
                exactly-at-the-end coordinates refuse. A real renderer with 8-bit colour should pass
                something like ``1/255``, because the last quantisation step is indistinguishable
                from the clamp.

        Returns:
            The value the display coordinate came from.

        Raises:
            ValueError: If ``t`` is outside ``[0, 1]``.
            ClampedChannelError: If ``t`` lies in a saturated band of a clamped end. The colour there
                is the same for every value beyond the clamp, so the inverse does not exist and
                returning ``display_min``/``display_max`` would be manufacturing a measurement.
        """
        t = float(t)
        if not (0.0 <= t <= 1.0):
            raise ValueError(f"normalised display coordinate must lie in [0, 1]; got {t}")
        band = float(saturation_band)
        if band < 0.0 or band >= 0.5:
            raise ValueError(f"saturation_band must lie in [0, 0.5); got {saturation_band}")
        if self.clamped_low and t <= band:
            raise ClampedChannelError(
                f"channel {self.channel.name!r} is clamped at the low end, so display coordinate "
                f"{t} has no unique value: every value at or below {self.display_min:.6g} "
                f"{self.channel.units} renders identically, and the field goes down to "
                f"{self.actual_min:.6g}. Widen display_min and re-render, or read the recorded "
                "actual extremum instead of the picture."
            )
        if self.clamped_high and t >= 1.0 - band:
            raise ClampedChannelError(
                f"channel {self.channel.name!r} is clamped at the high end, so display coordinate "
                f"{t} has no unique value: every value at or above {self.display_max:.6g} "
                f"{self.channel.units} renders identically, and the field goes up to "
                f"{self.actual_max:.6g}. Widen display_max and re-render, or read the recorded "
                "actual extremum instead of the picture."
            )
        if self.scale is ScaleKind.LOG10:
            lo, hi = math.log10(self.display_min), math.log10(self.display_max)
            return float(10.0 ** (lo + t * (hi - lo)))
        return float(self.display_min + t * (self.display_max - self.display_min))

    def quantitative_extremum(self, which: str) -> float:
        """The measured extremum, refusing when that end was clamped in the render.

        A clamped render is still allowed to *carry* the actual extremum — that is the whole point of
        recording it — but a caller asking this object for a number to quote is asking about what the
        render supports, and a saturated end supports nothing. Read
        :attr:`actual_min`/:attr:`actual_max` directly when the intent is "what did the field do",
        and call this when the intent is "what may I claim from this view".

        Args:
            which: ``"min"`` or ``"max"``.

        Returns:
            The measured extremum at that end.

        Raises:
            ClampedChannelError: If that end saturated.
            ValueError: If ``which`` is neither ``"min"`` nor ``"max"``.
        """
        if which not in ("min", "max"):
            raise ValueError(f"which must be 'min' or 'max'; got {which!r}")
        if which == "min" and self.clamped_low:
            raise ClampedChannelError(
                f"refusing to report a minimum from a low-clamped view of "
                f"{self.channel.name!r}: {self.clamp_report()}"
            )
        if which == "max" and self.clamped_high:
            raise ClampedChannelError(
                f"refusing to report a maximum from a high-clamped view of "
                f"{self.channel.name!r}: {self.clamp_report()}"
            )
        return self.actual_min if which == "min" else self.actual_max

    def require_quantitative(self) -> None:
        """Raise unless this bar may support a quantitative claim at all.

        Raises:
            ClampedChannelError: When either end is clamped.
        """
        if self.is_clamped:
            raise ClampedChannelError(
                f"channel {self.channel.name!r} may not support a quantitative claim: "
                f"{self.clamp_report()}"
            )

    def as_json_obj(self) -> dict[str, Any]:
        """Return the JSON-able colour-bar block, clamp state included.

        ``is_clamped`` is written even though a reader could derive it, because a reader who has to
        derive it is a reader who might not.
        """
        return {
            "channel": self.channel.as_json_obj(),
            "units": self.channel.units,
            "actual_min": self.actual_min,
            "actual_max": self.actual_max,
            "display_min": self.display_min,
            "display_max": self.display_max,
            "scale": self.scale.value,
            "palette": self.palette,
            "is_clamped": self.is_clamped,
            "clamped_low": self.clamped_low,
            "clamped_high": self.clamped_high,
            "clamp_excess_low": self.clamp_excess_low,
            "clamp_excess_high": self.clamp_excess_high,
            "clamp_report": self.clamp_report(),
            "supports_quantitative_claim": not self.is_clamped,
        }


def colour_bar_for(
    buffer: Any,
    *,
    display_min: float | None = None,
    display_max: float | None = None,
    scale: ScaleKind = ScaleKind.LINEAR,
    palette: str = "viridis",
) -> ColourBar:
    """Build a colour bar from a :class:`~aleph.viz.buffers.DerivedBuffer`.

    The actual extrema come from the buffer's *active* population, never from its full allocation.
    Reading extrema over an over-allocated buffer's unwritten tail is how a padded zero becomes the
    reported minimum of a field that never goes near zero.

    Args:
        buffer: A derived buffer with ``extrema``, ``units``, ``owner``, ``name``, ``derivation``.
        display_min: Lower display clamp; defaults to the actual minimum, i.e. no clamping.
        display_max: Upper display clamp; defaults to the actual maximum.
        scale: Display mapping.
        palette: Colour map name.

    Returns:
        The colour bar.
    """
    actual_min, actual_max = buffer.extrema
    lo = actual_min if display_min is None else float(display_min)
    hi = actual_max if display_max is None else float(display_max)
    if not (hi > lo):
        # A constant field has no range to display. Widen symmetrically by a unit rather than
        # refusing: a uniform field is a legitimate physical state and must still be renderable.
        span = max(abs(lo), 1.0)
        lo, hi = lo - 0.5 * span, hi + 0.5 * span
    return ColourBar(
        channel=ScalarChannel(
            name=buffer.name,
            units=buffer.units,
            owner=buffer.owner,
            derivation=buffer.derivation,
        ),
        actual_min=actual_min,
        actual_max=actual_max,
        display_min=lo,
        display_max=hi,
        scale=scale,
        palette=palette,
    )
