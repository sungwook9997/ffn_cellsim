r"""Arp2/3 70° branch-angle harmonic is wired into the ac/cell inner force-assembly loop.

GPU-only contract (CLAUDE.md I0-A, PI 2026-07-16): ``ac`` kernels are Warp-CUDA and are **never launched on
the CPU device**, not even in tests. So — following the sanctioned dev-box pattern (Mac = kernel SOURCE +
NumPy reference; native launch = the Lead on the A5000) — this module proves the wiring + physics WITHOUT ever
launching the kernel or allocating a device array:

  * launch-plan / gating (pure Python control flow, ``wp.launch`` MOCKED — no device work): ``driver
    ._accumulate_all`` launches ``branch_angle_kernel`` when ``AssembledCell.branch_triples_d`` is populated
    (mixed formin+Arp2/3 cortex) and launches NOTHING extra when it is ``None`` (the formin-only default) — so
    the default build stays bit-identical.
  * force MATH: a pure-NumPy reference of the harmonic-angle force, read faithfully off the ``@wp.kernel``
    body, is exercised on hand-built triples (at θ₀, over-opened, over-closed) — a junction sitting AT θ₀ is
    force-free, off θ₀ the force RESTORES toward θ₀ (correct sign both ways) with torque ``k_θ|θ-θ₀|``, and the
    three nodal forces sum to zero (internal potential). The reference also matches the shipped host oracle
    ``ac.weave.branch_angle.angle_forces`` bit-for-formula.
  * kernel-only properties (dormant-daughter gating + the atomic scatter into the GLOBAL force accumulator +
    the degenerate-arm guard) are asserted at the kernel SOURCE level — native-validated by the Lead.

The kernel/constants are reused from the lamellipodium path (``branch_angle_kernel`` + Faessler-2020 anchors),
not reinvented; this test proves the ac/cell wiring gates them and that the bound expression is the real
angle-harmonic restoring force.
"""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

# Importing the driver / kernel modules requires warp importable (they carry @wp.kernel bodies), but NOTHING
# here launches a kernel or allocates a device array (GPU-only contract: no CPU-Warp launches in ac tests).
pytest.importorskip("warp")

import aleph.components.incumbent.driver as driver  # noqa: E402
import aleph.components.weave.branch_angle_warp as _bw  # noqa: E402
from aleph.components.weave.branch_angle import (  # noqa: E402
    ARP23_K_THETA,
    ARP23_THETA0_RAD,
    angle_forces,
)
from aleph.components.weave.branch_angle_warp import branch_angle_kernel  # noqa: E402

_KERNEL_SRC = Path(_bw.__file__).read_text(encoding="utf-8")


# ----------------------------------------------------------------------------------------------------
# Wiring: launch-plan gating of driver._accumulate_all (pure control flow, wp.launch MOCKED — no device work).
# ----------------------------------------------------------------------------------------------------


def _bare_cell(**over: object) -> SimpleNamespace:
    """A minimal duck-typed AssembledCell with every optional primitive OFF except the branch fields.

    ``device`` is left ``None``: the sole consumer (``wp.launch(..., device=cell.device)``) is MOCKED in these
    wiring tests, so no device is ever selected and no kernel is launched — this exercises only the Python gate.
    """
    base = dict(
        device=None, n_total=3, n_tri=0, n_xl=0,
        myosin=None, nucleus=None, membrane=None, membrane_pressure=None, steric=None, pressure=None,
        branch_triples_d=None, branch_active_d=None, n_branch=0,
    )
    base.update(over)
    return SimpleNamespace(**base)


def _record_launches(monkeypatch) -> list[object]:
    """Capture the kernel object of every ``wp.launch`` the driver issues (mock — no real device work)."""
    seen: list[object] = []

    def _fake_launch(kernel, *_, **__):
        seen.append(kernel)

    monkeypatch.setattr(driver.wp, "launch", _fake_launch)
    return seen


def test_formin_only_default_launches_no_branch_kernel(monkeypatch):
    """branch_triples_d is None (formin-only cortex) ⇒ the branch kernel is NEVER launched (bit-identical)."""
    seen = _record_launches(monkeypatch)
    cell = _bare_cell(branch_triples_d=None, branch_active_d=None, n_branch=0)
    driver._accumulate_all(cell, pos=None, f=None)
    assert branch_angle_kernel not in seen
    # the force accumulator is still zeroed (the only launch on this stripped cell)
    assert seen == [driver._zero]


def test_arp23_cortex_launches_branch_kernel(monkeypatch):
    """branch_triples_d populated (mixed Arp2/3 cortex) ⇒ _accumulate_all launches the angle-harmonic kernel.

    The triples/active buffers are opaque sentinels: with ``wp.launch`` mocked the driver only reads the gate
    predicate (``branch_triples_d is not None and n_branch``), never the buffer contents — so no device array
    is allocated."""
    seen = _record_launches(monkeypatch)
    cell = _bare_cell(branch_triples_d=object(), branch_active_d=object(), n_branch=1)
    driver._accumulate_all(cell, pos=None, f=None)
    assert branch_angle_kernel in seen


# ----------------------------------------------------------------------------------------------------
# Force MATH: pure-NumPy reference of the harmonic-angle force, read faithfully off the @wp.kernel body.
# ----------------------------------------------------------------------------------------------------


def _branch_forces_np(theta: float, d1: float = 1.0, d2: float = 1.0):
    """NumPy mirror of ``branch_angle_kernel``: nodal forces of ``U = ½ k_θ (θ - θ₀)²`` on one triple.

    Faithful transcription of the kernel body (``pref = k_θ (θ - θ₀)/sin θ``;
    ``F_i = pref(r₂/(d₁d₂) - c r₁/d₁²)``, ``F_k = pref(r₁/(d₁d₂) - c r₂/d₂²)``, ``F_j = -(F_i+F_k)``). The
    triple is the mother arm along +x and the daughter arm at angle ``theta`` from it in the xy-plane, vertex
    ``j`` at the origin. Returns ``(F_mother, F_branch, F_daughter)`` [pN]. The device kernel then atomically
    SCATTERS these into the global accumulator — a kernel-only step asserted at the source level below."""
    j = np.zeros(3)
    i = np.array([d1, 0.0, 0.0])
    k = np.array([d2 * np.cos(theta), d2 * np.sin(theta), 0.0])
    r1 = i - j
    r2 = k - j
    dd1 = float(np.linalg.norm(r1))
    dd2 = float(np.linalg.norm(r2))
    c = float(np.clip(np.dot(r1, r2) / (dd1 * dd2), -1.0, 1.0))
    s = float(np.sqrt(max(1.0 - c * c, 1e-12)))
    th = float(np.arccos(c))
    pref = ARP23_K_THETA * (th - ARP23_THETA0_RAD) / s
    f_i = pref * (r2 / (dd1 * dd2) - c * r1 / (dd1 * dd1))
    f_k = pref * (r1 / (dd1 * dd2) - c * r2 / (dd2 * dd2))
    f_j = -(f_i + f_k)
    return f_i, f_j, f_k


def test_branch_force_is_zero_at_rest():
    """A junction sitting AT θ₀ is force-free (the harmonic potential is at its minimum, ∇U = 0)."""
    f_i, f_j, f_k = _branch_forces_np(ARP23_THETA0_RAD)
    assert np.linalg.norm(f_i) < 1e-12
    assert np.linalg.norm(f_j) < 1e-12
    assert np.linalg.norm(f_k) < 1e-12


def test_branch_force_restores_toward_theta0_both_ways():
    """Off θ₀ the force RESTORES: an over-opened junction (θ>θ₀) is pulled closed, an over-closed one is pushed
    open. The generalized force conjugate to θ on the daughter (``F_k·t̂``, t̂ the +θ tangent) has the sign of
    ``θ₀-θ`` and magnitude ``k_θ|θ-θ₀|/d₂`` — i.e. a restoring torque ``k_θ(θ-θ₀)``, exactly ``-dU/dθ``."""
    for dtheta, d2 in ((+0.30, 1.0), (-0.30, 0.7), (+0.15, 1.3), (-0.20, 0.5)):
        theta = ARP23_THETA0_RAD + dtheta
        _f_i, _f_j, f_k = _branch_forces_np(theta, d1=1.0, d2=d2)
        assert np.linalg.norm(f_k) > 1e-6  # a real, non-zero restoring force — the branch is mechanically live
        t_hat = np.array([-np.sin(theta), np.cos(theta), 0.0])  # +θ tangent at the daughter node
        generalized = float(np.dot(f_k, t_hat))
        # Sign: restoring toward θ₀ (opposite the displacement).
        assert np.sign(generalized) == np.sign(-dtheta)
        # Magnitude: |F_k·t̂|·d₂ is the restoring torque k_θ|θ-θ₀| (=|dU/dθ|).
        torque = abs(generalized) * d2
        assert np.isclose(torque, ARP23_K_THETA * abs(dtheta), rtol=1e-9, atol=1e-12)


def test_branch_force_sums_to_zero_and_matches_host_oracle():
    """Newton's 3rd law: an internal angle potential exerts no NET force on the junction (ΣF = 0). And the
    NumPy mirror reproduces the shipped host oracle ``angle_forces`` bit-for-formula (the same expression the
    device kernel is a port of), so the reference faithfully stands in for the kernel math."""
    for theta in (ARP23_THETA0_RAD + 0.3, ARP23_THETA0_RAD - 0.4, ARP23_THETA0_RAD + 0.9):
        f_i, f_j, f_k = _branch_forces_np(theta, d1=1.1, d2=0.8)
        assert np.linalg.norm(f_i + f_j + f_k) < 1e-12
        # Cross-check against the module's documented NumPy oracle (mother, branch, daughter positions).
        p_mother = np.array([1.1, 0.0, 0.0])
        p_branch = np.zeros(3)
        p_daughter = np.array([0.8 * np.cos(theta), 0.8 * np.sin(theta), 0.0])
        o_i, o_j, o_k = angle_forces(p_mother, p_branch, p_daughter, ARP23_K_THETA, ARP23_THETA0_RAD)
        assert np.allclose(f_i, o_i, atol=1e-12)
        assert np.allclose(f_j, o_j, atol=1e-12)
        assert np.allclose(f_k, o_k, atol=1e-12)


# ----------------------------------------------------------------------------------------------------
# Kernel-only properties (dormant gating + atomic scatter into the GLOBAL accumulator + degenerate guard):
# no NumPy analogue, asserted at the kernel SOURCE level — native-validated by the Lead on the A5000.
# ----------------------------------------------------------------------------------------------------


def test_dormant_gating_and_atomic_scatter_are_in_the_kernel_source():
    """The device kernel skips DORMANT (inactive) branch triples and atomically SCATTERS the three nodal forces
    into the GLOBAL node-force accumulator — the force-neutral-dormant-pool + accumulate behaviour that only the
    running kernel produces. Asserted at the source level; native-validated by the Lead (device launch)."""
    # Dormant daughters contribute no force (the N-fixed dormant pool stays force-neutral until nucleation).
    assert "if active[t] == 0:" in _KERNEL_SRC
    # The bound prefactor is the real harmonic-angle restoring coefficient k_θ(θ-θ₀)/sinθ.
    assert "pref = k_theta * (theta - theta0) / s" in _KERNEL_SRC
    # Atomic scatter of the three nodal forces into the shared global accumulator (Newton's 3rd on j).
    assert "wp.atomic_add(force, i, fi)" in _KERNEL_SRC
    assert "wp.atomic_add(force, k, fk)" in _KERNEL_SRC
    assert "wp.atomic_add(force, j, -(fi + fk))" in _KERNEL_SRC
    # Degenerate-arm guard (coincident nodes ⇒ no force) — kernel-only numerical safety.
    assert "if d1 < wp.float64(1e-12) or d2 < wp.float64(1e-12):" in _KERNEL_SRC
