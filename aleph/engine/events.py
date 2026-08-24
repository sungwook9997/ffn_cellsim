"""Device-resident event clock + reproducible RNG epoch for the whole-cell transaction.

Spec §6 (``WHOLE_CELL_COMMON_CONTRACTS_SPEC_2026-07-23.md``): one device-resident event clock,
advanced **only** on an accepted physical step, plus a per-epoch reproducible RNG so a rejected
candidate re-draws identically.  Inner mechanical iterations are not biological time — they neither
advance the clock nor the epoch.

The advance is device-gated on the acceptance predicate: both ``event_time`` and the accepted-step
index are incremented by ``dt_phys`` / ``1`` inside a Warp kernel *iff* ``accepted[0] != 0``, so no
acceptance value is ever read on the host.  ``base_seed`` is the run's fixed host seed; the per-epoch
draw derives on-device from ``base_seed`` and the accepted-step index (a consumer's ``propose_events``
kernel folds :attr:`WarpEventClock.accepted_step_index_d` into ``wp.rand_init``), so an unaccepted
retry sees the same index and reproduces its candidate draw.

This module defines the Warp kernel (import-safe on a CUDA-free host) and :class:`WarpEventClock`,
which is only constructible against CUDA arrays and is therefore exercised on the native lane.  The
CPU structural gates drive :class:`~aleph.engine.transaction.CellTransaction` with a spy clock
implementing the same :class:`~aleph.engine.runtime.EventClock` surface.

Sanity Gate:
    * ownership: the clock owns exactly two ``(1,)`` device scalars on one CUDA device; distinct storage.
    * dimensional: ``event_time_d`` [s], ``accepted_step_index_d`` [count]; ``dt_phys`` [s] > 0.
    * boundary/sign: advance is gated by ``accepted[0] != 0`` — a rejected step leaves both scalars fixed.
    * numerical: no host read of the predicate; the gate is a device branch, never a Python ``bool()``.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import warp as wp

__all__ = ["advance_event_clock_kernel", "WarpEventClock", "make_event_clock"]


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


@wp.kernel
def advance_event_clock_kernel(
    accepted: wp.array(dtype=wp.int32),
    dt_phys: wp.float64,
    event_time_d: wp.array(dtype=wp.float64),
    accepted_step_index_d: wp.array(dtype=wp.int64),
) -> None:
    """Advance biological time and the accepted-step index iff ``accepted[0] != 0`` (device-only).

    One thread (``dim=1``).  The acceptance predicate is consumed on the device; a rejected step
    (``accepted[0] == 0``) leaves both scalars untouched so the RNG epoch re-draws identically.
    """
    if accepted[0] != wp.int32(0):
        event_time_d[0] = event_time_d[0] + dt_phys
        accepted_step_index_d[0] = accepted_step_index_d[0] + wp.int64(1)


@dataclass(frozen=True, slots=True)
class WarpEventClock:
    """Device-resident event clock; advanced only on an accepted step (native lane).

    Args:
        base_seed: fixed host RNG seed for the whole run (a nonnegative int); the per-epoch draw
            derives on-device from this and :attr:`accepted_step_index_d`.
        event_time_d: ``(1,)`` float64 CUDA scalar — accumulated biological time [s].
        accepted_step_index_d: ``(1,)`` int64 CUDA scalar — number of accepted physical steps (the epoch).

    Neither scalar is ever read on the host inside the physical-time loop; :meth:`advance` gates on the
    device acceptance predicate.
    """

    base_seed: int
    event_time_d: wp.array
    accepted_step_index_d: wp.array

    def __post_init__(self) -> None:
        if isinstance(self.base_seed, bool) or not isinstance(self.base_seed, int) or self.base_seed < 0:
            raise ValueError("base_seed must be a nonnegative integer")
        checks = (
            ("event_time_d", self.event_time_d, wp.float64),
            ("accepted_step_index_d", self.accepted_step_index_d, wp.int64),
        )
        for label, array, dtype in checks:
            if not _device_is_cuda(array):
                raise ValueError(f"event clock {label} must be a Warp CUDA device array")
            if getattr(array, "dtype", None) != dtype:
                raise TypeError(f"event clock {label} must have dtype {dtype}")
            if getattr(array, "shape", None) != (1,):
                raise ValueError(f"event clock {label} must be a single-entry scalar")
        if str(self.event_time_d.device) != str(self.accepted_step_index_d.device):
            raise ValueError("event clock scalars must share one CUDA device")
        if _storage_key(self.event_time_d) == _storage_key(self.accepted_step_index_d):
            raise ValueError("event clock scalars must own distinct storage")

    @property
    def device(self) -> str:
        """CUDA device shared by the clock scalars."""
        return str(self.event_time_d.device)

    def advance(self, accepted: wp.array, dt_phys: float) -> None:
        """Device-gated advance of biological time and the accepted-step epoch.

        ``accepted`` is the same ``(1,)`` int32 predicate the transaction forwards to every
        participant's ``rollback``/``commit_irreversible``; this method never reads it on the host.
        """
        if not math.isfinite(dt_phys) or dt_phys <= 0.0:
            raise ValueError("dt_phys must be finite and positive")
        wp.launch(
            advance_event_clock_kernel,
            dim=1,
            inputs=[accepted, wp.float64(dt_phys), self.event_time_d, self.accepted_step_index_d],
            device=self.device,
        )


def make_event_clock(*, base_seed: int, device: str) -> WarpEventClock:
    """Allocate a zeroed :class:`WarpEventClock` on ``device`` (CUDA lane convenience constructor).

    Not importable-safe to call on a CUDA-free host (it allocates device memory); the CPU structural
    gates drive the transaction with a spy clock instead.  ``device`` is required — never a hard-coded
    id in production, the caller passes the resolved device (e.g. ``str(wp.get_device())``).
    """
    if isinstance(base_seed, bool) or not isinstance(base_seed, int) or base_seed < 0:
        raise ValueError("base_seed must be a nonnegative integer")
    return WarpEventClock(
        base_seed=base_seed,
        event_time_d=wp.zeros(1, dtype=wp.float64, device=device),
        accepted_step_index_d=wp.zeros(1, dtype=wp.int64, device=device),
    )
