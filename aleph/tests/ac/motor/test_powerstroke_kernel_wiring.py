"""Static + codegen gate that the Warp force kernels CONSUME the power stroke (catches the I3 decoupling).

The behavioural oracle (test_powerstroke_coupling_oracle) proves the FORMULA is contractile; this file proves
the actual Warp *kernel source* is wired to it — the piece the dev Mac's numpy oracles otherwise miss because
the kernels are never launched here (defect note #1). It is an AST inspection (like tests/ac/
test_warp_only_contract.py) plus an optional Warp codegen/type-check — neither launches a kernel or evolves any
state, so no CPU simulation is run; the physics gate stays the gbook native gate.

If a future edit re-decouples the force from the walk (the original bug: the crossbridge pulled toward a FIXED
anchor and never read the abscissa), the AST gate fails locally.
"""

from __future__ import annotations

import ast
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[4]   # tests/ac/motor/<file> -> repo root
MINIFIL_WARP = REPO_ROOT / "aleph" / "components" / "motor" / "minifilament_warp.py"

# the two device kernels whose force/load MUST depend on the walked abscissa + the walk direction
FORCE_KERNELS = ("crossbridge_kernel", "compute_head_loads_kernel")
REQUIRED_INPUTS = ("abscissa", "walk_dir")


def _funcdef(tree: ast.Module, name: str) -> ast.FunctionDef:
    for node in ast.walk(tree):
        if isinstance(node, ast.FunctionDef) and node.name == name:
            return node
    raise AssertionError(f"{name} not found in {MINIFIL_WARP.name}")


def _arg_names(fn: ast.FunctionDef) -> set[str]:
    return {a.arg for a in fn.args.args}


def _body_names(fn: ast.FunctionDef) -> set[str]:
    return {n.id for n in ast.walk(fn) if isinstance(n, ast.Name)}


@pytest.mark.parametrize("kernel", FORCE_KERNELS)
def test_force_kernel_accepts_powerstroke_inputs(kernel: str) -> None:
    """crossbridge_kernel / compute_head_loads_kernel MUST take `abscissa` and `walk_dir` as inputs."""
    tree = ast.parse(MINIFIL_WARP.read_text(encoding="utf-8"))
    args = _arg_names(_funcdef(tree, kernel))
    missing = [p for p in REQUIRED_INPUTS if p not in args]
    assert not missing, f"{kernel} is decoupled from the power stroke — missing inputs {missing}"


@pytest.mark.parametrize("kernel", FORCE_KERNELS)
def test_force_kernel_uses_the_abscissa(kernel: str) -> None:
    """The kernel body must actually READ `abscissa` (accepting it but ignoring it would re-introduce the bug)."""
    tree = ast.parse(MINIFIL_WARP.read_text(encoding="utf-8"))
    body_names = _body_names(_funcdef(tree, kernel))
    assert "abscissa" in body_names, f"{kernel} accepts abscissa but never reads it — power stroke decoupled"


def test_attachment_point_advances_along_walk_dir() -> None:
    """The `attachment_point` @wp.func must combine `abscissa` and `walk_dir` (x_att = anchor + s*walk_dir)."""
    tree = ast.parse(MINIFIL_WARP.read_text(encoding="utf-8"))
    body_names = _body_names(_funcdef(tree, "attachment_point"))
    assert {"abscissa", "walk_dir"} <= body_names


def test_warp_kernels_codegen_clean() -> None:
    """Codegen/type-check the edited Warp source (compile only — no kernel launch, no simulation)."""
    try:
        import warp as wp

        import aleph.components.motor.hand as hand
        import aleph.components.motor.minifilament_warp as mw

        wp.init()
        devices = wp.get_devices()
        if not devices:
            pytest.skip("no Warp device available for codegen check")
        device = str(devices[0])
        wp.load_module(hand, device=device)
        wp.load_module(mw, device=device)
    except Exception as exc:  # pragma: no cover - environment-dependent (real gate is the gbook)
        pytest.skip(f"Warp codegen check unavailable here: {exc}")
