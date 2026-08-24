"""Interface-residual measurement — the read that must precede any coupled-solve accelerator.

Two halves, tested where each can be tested. The STRUCTURAL half (which cuts exist, and what an
unmeasurable cut reports) is CPU-testable and is where the dangerous bug lives: a cut that is silently
omitted makes an incomplete world look balanced. The PHYSICAL half (does the gate actually fail on a
one-sided scatter) needs CUDA and carries the negative control.
"""

from __future__ import annotations

import numpy as np
import pytest
import warp as wp

from aleph.engine.contracts import reference_cell_architecture
from aleph.engine.coupled_solve import (
    CutResidual,
    InterfaceResidualReport,
    measure_interface_residual,
)

wp.init()
_CUDA_DEVICE = next((str(d) for d in wp.get_devices() if d.is_cuda), None)

#: The structural cases below never reach device code — every cut short-circuits on UNBOUND before a
#: kernel or an allocation. This is passed where the signature wants a device string, and it is NOT
#: "cpu": naming the CPU as a device in an ac test is forbidden outright, and rightly, because the
#: rule cannot tell a path that would launch there from one that never gets the chance.
_NO_DEVICE = "no-device-is-touched-on-this-path"
ARCH = reference_cell_architecture()


# ── structural: absence must never read as balance ────────────────────────────────────────────────

def test_an_unmeasurable_cut_is_reported_not_omitted() -> None:
    """A world binding nothing must report every cut UNBOUND — never an empty, balanced-looking list."""
    report = measure_interface_residual(
        world=None, architecture=ARCH, device=_NO_DEVICE, component_runtime=lambda _n: None,
    )
    assert len(report.cuts) == len(ARCH.balance_cuts()) == 27
    assert len(report.measured) == 0
    assert len(report.imbalanced) == 0
    assert report.worst is None
    assert {c.status for c in report.cuts} == {"UNBOUND"}
    assert all(c.residual_pn is None for c in report.cuts)


def test_a_component_without_a_force_array_is_unbound_not_zero() -> None:
    """A bound runtime that owns no force array cannot be scored; 0.0 pN would be a lie."""

    class _NoForce:
        pass

    report = measure_interface_residual(
        world=None, architecture=ARCH, device=_NO_DEVICE, component_runtime=lambda _n: _NoForce(),
    )
    assert {c.status for c in report.cuts} == {"UNBOUND"}


def test_the_report_keeps_the_pair_that_cannot_be_attributed() -> None:
    """cortex/membrane carries two power ports, so its reading names the PAIR, not one edge."""
    report = measure_interface_residual(
        world=None, architecture=ARCH, device=_NO_DEVICE, component_runtime=lambda _n: None,
    )
    shared = [c for c in report.cuts if not c.attributable]
    assert len(shared) == 1
    assert shared[0].components == ("cortex", "membrane")
    assert set(shared[0].connectors) == {"membrane_erm_cortex", "membrane_cortex_contact"}
    assert sum(1 for c in report.cuts if c.attributable) == 26


def test_as_dict_states_what_the_number_is_not() -> None:
    """The artifact must carry its own scope — a residual is not a convergence claim."""
    payload = measure_interface_residual(
        world=None, architecture=ARCH, device=_NO_DEVICE, component_runtime=lambda _n: None,
    ).as_dict()
    assert payload["n_cuts_declared"] == 27
    assert payload["n_unbound"] == 27
    assert "NOT a convergence claim" in payload["scope"]
    assert "no position was updated" in payload["scope"]


def test_worst_ignores_unbound_cuts() -> None:
    report = InterfaceResidualReport((
        CutResidual(("a", "b"), ("x",), "UNBOUND", None, None, None),
        CutResidual(("c", "d"), ("y",), "BALANCED", 1e-14, 1.0, 1e-16),
        CutResidual(("e", "f"), ("z",), "IMBALANCED", 3.0, 10.0, 1e-16),
    ))
    assert report.worst is not None and report.worst.residual_pn == 3.0
    assert len(report.measured) == 2
    assert len(report.imbalanced) == 1


# ── physical: the gate must FAIL on a one-sided scatter ───────────────────────────────────────────

def _pair_runtimes(force_a: np.ndarray, force_b: np.ndarray, device: str):
    """Bind two force arrays at the endpoints of one real cut and leave every other cut unbound."""
    cut = next(iter(ARCH.attributable_cuts()))

    class _Body:
        def __init__(self, arr):
            self.force_d = wp.array(np.ascontiguousarray(arr), dtype=wp.vec3d, device=device)

    bodies = {cut[0]: _Body(force_a), cut[1]: _Body(force_b)}
    return cut, (lambda name: bodies.get(name))


@pytest.mark.skipif(_CUDA_DEVICE is None, reason="I0-A: the balance gate requires a CUDA GPU")
def test_an_equal_and_opposite_pair_balances() -> None:
    """The positive reading: a genuine Newton pair cancels within the derived floor."""
    f = np.array([[1.0, -2.0, 0.5], [0.25, 0.0, -3.0]])
    cut, lookup = _pair_runtimes(f, -f, _CUDA_DEVICE)
    report = measure_interface_residual(
        world=None, architecture=ARCH, device=_CUDA_DEVICE, component_runtime=lookup,
    )
    scored = next(c for c in report.cuts if c.components == cut)
    assert scored.status == "BALANCED", scored
    assert report.imbalanced == ()


@pytest.mark.skipif(_CUDA_DEVICE is None, reason="I0-A: the balance gate requires a CUDA GPU")
def test_negative_control_a_one_sided_scatter_is_caught() -> None:
    """THE control. Without it, "everything balanced" and "the gate cannot fail" are one artifact.

    One endpoint keeps its half of the pair; the other is short by a finite force. That is exactly the
    one-sided scatter `ledger.py` says the gate exists to catch, and it must be caught.
    """
    f = np.array([[1.0, -2.0, 0.5], [0.25, 0.0, -3.0]])
    broken = -f.copy()
    broken[0, 0] += 0.5  # one node's reaction never scattered back
    cut, lookup = _pair_runtimes(f, broken, _CUDA_DEVICE)
    report = measure_interface_residual(
        world=None, architecture=ARCH, device=_CUDA_DEVICE, component_runtime=lookup,
    )
    scored = next(c for c in report.cuts if c.components == cut)
    assert scored.status == "IMBALANCED", scored
    assert scored.residual_pn is not None and scored.residual_pn == pytest.approx(0.5, rel=1e-9)
    assert report.worst is not None and report.worst.components == cut


@pytest.mark.skipif(_CUDA_DEVICE is None, reason="I0-A: the balance gate requires a CUDA GPU")
def test_a_sign_flip_is_caught_too() -> None:
    """A scatter with the wrong sign doubles instead of cancelling — the other failure mode."""
    f = np.array([[1.0, -2.0, 0.5], [0.25, 0.0, -3.0]])
    cut, lookup = _pair_runtimes(f, f, _CUDA_DEVICE)  # NOT negated
    report = measure_interface_residual(
        world=None, architecture=ARCH, device=_CUDA_DEVICE, component_runtime=lookup,
    )
    scored = next(c for c in report.cuts if c.components == cut)
    assert scored.status == "IMBALANCED", scored


# ── a BALANCED verdict that means nothing must be visible as such ─────────────────────────────────

def test_cancellation_is_the_reading_that_says_whether_balanced_means_anything() -> None:
    """A genuine pair cancels (ratio -> 0). Nothing-to-cancel does not (ratio -> 1)."""
    genuine = CutResidual(("a", "b"), ("x",), "BALANCED", 1e-16, 4.0, 1e-16, ("x",))
    assert genuine.cancellation is not None and genuine.cancellation < 1e-10

    nothing = CutResidual(("a", "b"), ("x",), "BALANCED", 4.99e-25, 4.99e-25, 1e-16, ("x",))
    assert nothing.cancellation == pytest.approx(1.0), "residual == scale means nothing cancelled"


def test_an_unwired_cut_is_flagged_however_balanced_it_looks() -> None:
    """Two arrays with no connector between them cannot fail to balance. That is not evidence."""
    unwired = CutResidual(("a", "b"), ("x",), "BALANCED", 0.0, 0.0, 1e-16, ())
    assert unwired.unwired
    wired = CutResidual(("a", "b"), ("x",), "BALANCED", 1e-16, 4.0, 1e-16, ("x",))
    assert not wired.unwired

    report = InterfaceResidualReport((unwired, wired))
    assert report.scoreable == (wired,), "an unwired cut must never be counted as scoreable"


def test_scoreable_excludes_the_forceless_cut() -> None:
    """No force means no cancellation ratio, so there is nothing to judge — excluded structurally."""
    forceless = CutResidual(("a", "b"), ("x",), "BALANCED", 0.0, 0.0, 1e-16, ("x",))
    assert forceless.cancellation is None
    assert InterfaceResidualReport((forceless,)).scoreable == ()


def test_the_scope_tells_the_reader_to_read_cancellation_first() -> None:
    payload = InterfaceResidualReport((
        CutResidual(("a", "b"), ("x",), "BALANCED", 1.0, 1.0, 1e-16, ("x",)),
    )).as_dict()
    assert payload["n_scoreable"] == 1
    assert payload["cuts"][0]["cancellation"] == pytest.approx(1.0)
    assert "cancellation is near 1 cancelled nothing" in payload["scope"]
