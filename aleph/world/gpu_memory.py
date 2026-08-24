"""Exact process-lifetime CUDA peak-memory accounting through the installed NVML driver library.

NVML's accounting statistic ``maxMemoryUsage`` is the maximum total device memory, in bytes, ever allocated
by one compute/graphics process.  Unlike Warp's mempool high-water mark or an endpoint free-memory sample, it
includes allocations outside Warp's pool.  The NVIDIA driver must have accounting mode enabled *before* the
simulation process starts; enabling it requires administrator authority and is intentionally never attempted
by this module.

The binding uses only Python's standard-library :mod:`ctypes` and ``libnvidia-ml.so.1`` already supplied by the
driver.  Device selection is by Warp's GPU UUID, never by a hard-coded or potentially remapped ordinal.

Reference:
    NVIDIA NVML API Reference, ``nvmlAccountingStats_t`` and Accounting Statistics:
    https://docs.nvidia.com/deploy/nvml-api/group__nvmlAccountingStats.html

Why it lives in ``world/`` and not where it was written.  It was
``components/incumbent/gpu_memory_accounting.py``, and the arena — CANONICAL since 2026-08-20 —
may not import ``components/`` (``test_layer_directions``).  So ``world/build/__init__.py`` reported
``exact_peak_gpu_bytes: None`` with a status string naming this module and calling the port "an open
item", while ``scripts/world_phase1_native.py`` reached back into the frozen tree for it anyway —
outside the guard's scope, which is the package.  This is that port.  The module binds only
``ctypes`` and ``warp``, neither of which is a layer this may not see, so nothing about it needed to
change; only where it sits did.

No re-export shim was left behind, and that is a departure from how ``PROVENANCE`` was moved on
2026-08-20.  The reason is the test: ``test_gpu_memory_accounting`` monkeypatches ``_load_nvml`` and
``_NvmlAccountingStats`` on the module object, and a re-exporting shim would take the patch on
itself while ``query_process_peak`` resolved those names in the real module — an inert patch
against live NVML, which on a machine with no driver is a failure and on a machine with one is
worse.  There are four import sites in the whole tree; rewriting four lines is smaller than an
aliasing trick and it leaves nothing pointing at the old home.  The ledger KEY names
(``gpu_memory_accounting_*``) are artifact fields read by ``assemble.py`` and ``ng2_ng3_ng6.py`` and
are deliberately NOT renamed: they name a record format, not a module.

⚠ Importing this module does NOT enable accounting mode and does not create a CUDA context.  The
ordering constraint the drivers observe — request accounting before the context exists — is about
``wp.init()``, not about this import, and it is unchanged by the move.

Sanity Gate:
    * Units: ``maxMemoryUsage`` is returned by NVML in bytes and stored without scaling or sampling.
    * Scope: the value is the current process's lifetime maximum on the selected UUID; the PID and UUID are
      recorded alongside it so concurrent-device or process ambiguity is visible.
    * Boundary cases: disabled/unsupported accounting, missing process stats, unavailable sentinel, a missing
      NVML library, and every NVML error remain explicit non-exact statuses.
    * Authority: this module never calls ``nvmlDeviceSetAccountingMode`` and never invokes ``sudo``.
    * Timing: callers query only after the physical-time loop; no simulation-state D2H transfer is introduced.
"""

from __future__ import annotations

import argparse
import ctypes
import ctypes.util
import json
import os
from dataclasses import asdict, dataclass

import warp as wp

__all__ = ["NvmlAccountingResult", "query_process_peak"]

_NVML_SUCCESS = 0
_NVML_FEATURE_ENABLED = 1
_NVML_VALUE_NOT_AVAILABLE_ULL = (1 << 64) - 1


class _NvmlAccountingStats(ctypes.Structure):
    """ABI layout of ``nvmlAccountingStats_t`` from NVIDIA's public NVML header."""

    _fields_ = [
        ("gpuUtilization", ctypes.c_uint),
        ("memoryUtilization", ctypes.c_uint),
        ("maxMemoryUsage", ctypes.c_ulonglong),
        ("time", ctypes.c_ulonglong),
        ("startTime", ctypes.c_ulonglong),
        ("isRunning", ctypes.c_uint),
        ("reserved", ctypes.c_uint * 5),
    ]


@dataclass(frozen=True)
class NvmlAccountingResult:
    """One post-loop NVML accounting result suitable for the native population ledger."""

    probe_status: str
    device_uuid: str
    pid: int
    accounting_enabled: bool | None
    exact_peak_bytes: int | None = None
    error: str | None = None

    @property
    def exact(self) -> bool:
        """Whether NVML returned a usable process-lifetime maximum in bytes."""
        return self.probe_status == "EXACT_PROCESS_LIFETIME_PEAK" and self.exact_peak_bytes is not None

    def ledger_fields(self) -> dict[str, object]:
        """Return namespaced fields without relabelling a failed probe as exact accounting."""
        fields: dict[str, object] = {
            "gpu_memory_accounting_backend": "NVML_PROCESS_LIFETIME_ACCOUNTING",
            "gpu_memory_accounting_probe_status": self.probe_status,
            "gpu_memory_accounting_device_uuid": self.device_uuid,
            "gpu_memory_accounting_pid": self.pid,
            "gpu_memory_accounting_mode_enabled": self.accounting_enabled,
        }
        if self.exact_peak_bytes is not None:
            fields["gpu_bytes_peak_exact"] = self.exact_peak_bytes
        if self.error:
            fields["gpu_memory_accounting_probe_error"] = self.error
        return fields


def _load_nvml() -> ctypes.CDLL:
    """Load the driver-owned NVML shared library without adding a Python package dependency."""
    path = ctypes.util.find_library("nvidia-ml") or "libnvidia-ml.so.1"
    return ctypes.CDLL(path)


def _bind(lib: ctypes.CDLL) -> None:
    """Declare only the stable NVML functions used by this probe."""
    lib.nvmlInit_v2.argtypes = []
    lib.nvmlInit_v2.restype = ctypes.c_int
    lib.nvmlShutdown.argtypes = []
    lib.nvmlShutdown.restype = ctypes.c_int
    lib.nvmlErrorString.argtypes = [ctypes.c_int]
    lib.nvmlErrorString.restype = ctypes.c_char_p
    lib.nvmlDeviceGetHandleByUUID.argtypes = [ctypes.c_char_p, ctypes.POINTER(ctypes.c_void_p)]
    lib.nvmlDeviceGetHandleByUUID.restype = ctypes.c_int
    lib.nvmlDeviceGetAccountingMode.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_uint)]
    lib.nvmlDeviceGetAccountingMode.restype = ctypes.c_int
    lib.nvmlDeviceGetAccountingStats.argtypes = [
        ctypes.c_void_p,
        ctypes.c_uint,
        ctypes.POINTER(_NvmlAccountingStats),
    ]
    lib.nvmlDeviceGetAccountingStats.restype = ctypes.c_int


def _error_string(lib: ctypes.CDLL, code: int) -> str:
    raw = lib.nvmlErrorString(code)
    return raw.decode("utf-8", errors="replace") if raw else f"NVML error {code}"


def query_process_peak(device: str | None = None, *, pid: int | None = None) -> NvmlAccountingResult:
    """Query the exact NVML process-lifetime peak, or return an explicit non-exact status.

    Args:
        device: Warp CUDA device alias. ``None`` resolves Warp's current CUDA device.
        pid: Process to query; defaults to this Python process. Exposed only for diagnostic tests.

    Returns:
        A result whose :attr:`NvmlAccountingResult.exact` property is true only for a successful, enabled,
        non-sentinel NVML accounting value.
    """
    process_id = os.getpid() if pid is None else int(pid)
    try:
        dev = wp.get_device(device)
    except Exception as exc:  # noqa: BLE001 - result must preserve the exact unavailable reason
        return NvmlAccountingResult(
            "DEVICE_UNAVAILABLE", "", process_id, None, error=f"{type(exc).__name__}: {exc}")
    device_uuid = str(dev.uuid)
    if not dev.is_cuda:
        return NvmlAccountingResult("NON_CUDA_DEVICE", device_uuid, process_id, None)

    try:
        lib = _load_nvml()
        _bind(lib)
    except Exception as exc:  # noqa: BLE001 - missing driver library/symbol is an expected probe outcome
        return NvmlAccountingResult(
            "NVML_UNAVAILABLE", device_uuid, process_id, None, error=f"{type(exc).__name__}: {exc}")

    initialized = False
    try:
        code = int(lib.nvmlInit_v2())
        if code != _NVML_SUCCESS:
            return NvmlAccountingResult(
                "NVML_INIT_FAILED", device_uuid, process_id, None, error=_error_string(lib, code))
        initialized = True

        handle = ctypes.c_void_p()
        code = int(lib.nvmlDeviceGetHandleByUUID(device_uuid.encode("ascii"), ctypes.byref(handle)))
        if code != _NVML_SUCCESS:
            return NvmlAccountingResult(
                "DEVICE_HANDLE_FAILED", device_uuid, process_id, None, error=_error_string(lib, code))

        mode = ctypes.c_uint(0)
        code = int(lib.nvmlDeviceGetAccountingMode(handle, ctypes.byref(mode)))
        if code != _NVML_SUCCESS:
            return NvmlAccountingResult(
                "ACCOUNTING_MODE_QUERY_FAILED", device_uuid, process_id, None, error=_error_string(lib, code))
        if mode.value != _NVML_FEATURE_ENABLED:
            return NvmlAccountingResult("ACCOUNTING_DISABLED", device_uuid, process_id, False)

        stats = _NvmlAccountingStats()
        code = int(lib.nvmlDeviceGetAccountingStats(handle, ctypes.c_uint(process_id), ctypes.byref(stats)))
        if code != _NVML_SUCCESS:
            return NvmlAccountingResult(
                "ACCOUNTING_STATS_QUERY_FAILED", device_uuid, process_id, True, error=_error_string(lib, code))
        if stats.maxMemoryUsage == _NVML_VALUE_NOT_AVAILABLE_ULL:
            return NvmlAccountingResult("ACCOUNTING_VALUE_NOT_AVAILABLE", device_uuid, process_id, True)
        return NvmlAccountingResult(
            "EXACT_PROCESS_LIFETIME_PEAK", device_uuid, process_id, True,
            exact_peak_bytes=int(stats.maxMemoryUsage),
        )
    except Exception as exc:  # noqa: BLE001 - ABI/symbol failures must stay visible in the ledger
        return NvmlAccountingResult(
            "NVML_QUERY_FAILED", device_uuid, process_id, None, error=f"{type(exc).__name__}: {exc}")
    finally:
        if initialized:
            lib.nvmlShutdown()


def main() -> None:
    """Print a standalone probe; useful immediately after the administrator enables accounting mode."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--device", default=None, help="Warp CUDA device alias (default: current CUDA device)")
    args = parser.parse_args()
    result = query_process_peak(args.device)
    print(json.dumps(asdict(result) | {"exact": result.exact}, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
