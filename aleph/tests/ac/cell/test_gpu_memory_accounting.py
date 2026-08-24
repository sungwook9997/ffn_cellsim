"""Pure-host contracts for exact NVML process-lifetime peak-memory accounting."""

from __future__ import annotations

import ctypes
from types import SimpleNamespace

from aleph.world import gpu_memory as accounting


class _FakeFunction:
    """ctypes-like callable whose ABI metadata can be assigned by ``_bind``."""

    def __init__(self, callback):
        self.callback = callback
        self.argtypes = None
        self.restype = None

    def __call__(self, *args):
        return self.callback(*args)


class _FakeNvml:
    def __init__(self, *, enabled: bool, peak_bytes: int = 0):
        self.stats_calls = 0
        self.nvmlInit_v2 = _FakeFunction(lambda: 0)
        self.nvmlShutdown = _FakeFunction(lambda: 0)
        self.nvmlErrorString = _FakeFunction(lambda _code: b"fake NVML error")
        self.nvmlDeviceGetHandleByUUID = _FakeFunction(self._handle)
        self.nvmlDeviceGetAccountingMode = _FakeFunction(
            lambda _handle, mode: self._write_mode(mode, enabled))
        self.nvmlDeviceGetAccountingStats = _FakeFunction(
            lambda _handle, _pid, stats: self._write_stats(stats, peak_bytes))

    @staticmethod
    def _handle(_uuid, handle) -> int:
        ctypes.cast(handle, ctypes.POINTER(ctypes.c_void_p))[0] = ctypes.c_void_p(123)
        return 0

    @staticmethod
    def _write_mode(mode, enabled: bool) -> int:
        ctypes.cast(mode, ctypes.POINTER(ctypes.c_uint))[0] = int(enabled)
        return 0

    def _write_stats(self, stats, peak_bytes: int) -> int:
        self.stats_calls += 1
        out = ctypes.cast(stats, ctypes.POINTER(accounting._NvmlAccountingStats))[0]
        out.maxMemoryUsage = peak_bytes
        out.isRunning = 1
        return 0


def test_nvml_accounting_stats_abi_layout() -> None:
    """The stable C ABI has no accidental Python padding or field reordering."""
    stats = accounting._NvmlAccountingStats
    assert ctypes.sizeof(stats) == 56
    assert stats.gpuUtilization.offset == 0
    assert stats.memoryUtilization.offset == 4
    assert stats.maxMemoryUsage.offset == 8
    assert stats.time.offset == 16
    assert stats.startTime.offset == 24
    assert stats.isRunning.offset == 32
    assert stats.reserved.offset == 36


def test_exact_peak_requires_enabled_successful_non_sentinel_nvml(monkeypatch) -> None:
    fake = _FakeNvml(enabled=True, peak_bytes=987_654_321)
    monkeypatch.setattr(accounting, "_load_nvml", lambda: fake)
    monkeypatch.setattr(
        accounting.wp, "get_device",
        lambda _device: SimpleNamespace(is_cuda=True, uuid="GPU-test-uuid"),
    )
    result = accounting.query_process_peak("cuda", pid=42)
    assert result.exact
    assert result.exact_peak_bytes == 987_654_321
    assert result.ledger_fields()["gpu_bytes_peak_exact"] == 987_654_321
    assert result.device_uuid == "GPU-test-uuid"
    assert result.pid == 42
    assert fake.stats_calls == 1


def test_disabled_accounting_stays_explicitly_non_exact(monkeypatch) -> None:
    fake = _FakeNvml(enabled=False)
    monkeypatch.setattr(accounting, "_load_nvml", lambda: fake)
    monkeypatch.setattr(
        accounting.wp, "get_device",
        lambda _device: SimpleNamespace(is_cuda=True, uuid="GPU-test-uuid"),
    )
    result = accounting.query_process_peak("cuda", pid=42)
    assert result.probe_status == "ACCOUNTING_DISABLED"
    assert result.accounting_enabled is False
    assert not result.exact
    assert "gpu_bytes_peak_exact" not in result.ledger_fields()
    assert fake.stats_calls == 0
