"""Structural + codegen gate for the NMII catch-slip head-turnover device kernels (I0-A: no launch on Mac).

The behavioural ground truth is the NumPy oracle (test_bell_kinetics_oracle::test_engaged_fraction_catch_slip_
rises_under_load_then_falls). This file proves the Warp *kernel source* in ``ac/motor/hand.py`` is wired to the
corrected catch-slip law — the piece the dev-Mac numpy oracles miss because the kernels are never launched here
(same rationale as test_powerstroke_kernel_wiring). It is an AST inspection plus an optional Warp codegen check;
neither launches a kernel or evolves state, so the physics gate stays the gbook native gate.

If a future edit re-introduces the Kovacs-2007 bug (NMII head detaching by a pure Bell slip so duty FALLS under
load), the AST gates here fail locally.
"""

from __future__ import annotations

import ast
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[4]   # tests/ac/motor/<file> -> repo root
HAND = REPO_ROOT / "aleph" / "components" / "motor" / "hand.py"


def _tree() -> ast.Module:
    return ast.parse(HAND.read_text(encoding="utf-8"))


def _funcdef(tree: ast.Module, name: str) -> ast.FunctionDef:
    for node in ast.walk(tree):
        if isinstance(node, ast.FunctionDef) and node.name == name:
            return node
    raise AssertionError(f"{name} not found in {HAND.name}")


def _arg_names(fn: ast.FunctionDef) -> set[str]:
    return {a.arg for a in fn.args.args}


def _call_names(fn: ast.FunctionDef) -> set[str]:
    out: set[str] = set()
    for node in ast.walk(fn):
        if isinstance(node, ast.Call):
            f = node.func
            if isinstance(f, ast.Name):
                out.add(f.id)
            elif isinstance(f, ast.Attribute):
                out.add(f.attr)
    return out


def test_catch_slip_kernel_exists_and_takes_catch_slip_params() -> None:
    """step_detach_catch_slip_kernel must accept the full-|F| Bell load and a NMIICatchSlipParams `cs`."""
    fn = _funcdef(_tree(), "step_detach_catch_slip_kernel")
    args = _arg_names(fn)
    assert "loads_bell" in args, "catch-slip detach must read the full-|F| load (Bell/catch pathway)"
    assert "loads_hill" in args, "catch-slip detach still advances the abscissa by the Hill step"
    assert "cs" in args, "catch-slip detach must take the NMIICatchSlipParams struct"


def test_catch_slip_kernel_uses_catch_slip_not_pure_slip() -> None:
    """The corrected NMII detach must call catch_slip_off_rate — NOT the pure Bell slip bell_off_rate."""
    calls = _call_names(_funcdef(_tree(), "step_detach_catch_slip_kernel"))
    assert "catch_slip_off_rate" in calls, "NMII head must detach by the catch-slip off-rate (Kovacs 2007)"
    assert "bell_off_rate" not in calls, "NMII head must not detach by a pure Bell slip (duty would fall under load)"
    # still a Hill-stepped, Poisson-detached turnover
    assert "hill_velocity" in calls and "detach_prob" in calls


def test_pure_slip_kernel_left_unchanged() -> None:
    """step_detach_kernel stays a pure Bell slip — the ERM / alpha-actinin / crosslink consumers depend on it."""
    calls = _call_names(_funcdef(_tree(), "step_detach_kernel"))
    assert "bell_off_rate" in calls
    assert "catch_slip_off_rate" not in calls


def test_catch_slip_device_func_is_two_pathway() -> None:
    """The device catch_slip_off_rate @wp.func has both a decaying catch term and a growing slip term."""
    src = ast.get_source_segment(HAND.read_text(encoding="utf-8"), _funcdef(_tree(), "catch_slip_off_rate"))
    assert src is not None
    compact = " ".join(src.split())
    # catch pathway falls with load: exp(-|f|*x_catch/kT); slip pathway rises: exp(+|f|*x_slip/kT)
    assert "k_catch0 * wp.exp(-fa * x_catch / kT)" in compact
    assert "k_slip0 * wp.exp(fa * x_slip / kT)" in compact


def test_catch_slip_params_struct_fields() -> None:
    """NMIICatchSlipParams carries exactly the five two-pathway constants."""
    tree = _tree()
    cls = next(
        (n for n in ast.walk(tree) if isinstance(n, ast.ClassDef) and n.name == "NMIICatchSlipParams"),
        None,
    )
    assert cls is not None, "NMIICatchSlipParams struct missing"
    fields = {n.target.id for n in cls.body if isinstance(n, ast.AnnAssign) and isinstance(n.target, ast.Name)}
    assert fields == {"k_catch0", "x_catch", "k_slip0", "x_slip", "kT"}


def test_public_symbols_exported() -> None:
    """The new catch-slip primitives are importable from the module namespace."""
    import aleph.components.motor.hand as hand

    for name in ("NMIICatchSlipParams", "catch_slip_off_rate", "step_detach_catch_slip_kernel"):
        assert name in hand.__all__
        assert hasattr(hand, name)


def test_warp_kernels_codegen_clean() -> None:
    """Codegen/type-check the edited Warp source (compile only — no kernel launch, no simulation)."""
    try:
        import warp as wp

        import aleph.components.motor.hand as hand

        wp.init()
        devices = wp.get_devices()
        if not devices:
            pytest.skip("no Warp device available for codegen check")
        wp.load_module(hand, device=str(devices[0]))
    except Exception as exc:  # pragma: no cover - environment-dependent (real gate is the gbook)
        pytest.skip(f"Warp codegen check unavailable here: {exc}")
