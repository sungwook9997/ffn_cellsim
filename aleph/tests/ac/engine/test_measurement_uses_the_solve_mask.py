"""Every force assembly in the driver must run under the SAME omit mask as the solve it describes.

WHAT THIS CAUGHT.  On 2026-07-29 the `--engine-cortex` A/B failed its declared criterion: with the engine
component supplying the cortex, γ_total differed from the incumbent arm by a CONSTANT 0.41526 pN/µm while
``max_f_cortex_pn`` and ``n_bound`` agreed to the last digit at every sample.  Identical dynamics with a
constant offset in one observable is not a physics disagreement — it is a measurement reading a different
force field from the one being integrated.

The cause: ``make_inner_solve(..., omit=omit_channels)`` masked the incumbent's cortex channels for the
SOLVE, but ``_measure_gamma`` called ``_accumulate_all`` with no mask.  The engine's cortex arrives through
the ``cell.myosin`` hook, which is not an omittable channel, so the unmasked measurement summed the
incumbent's bending/crosslink AND the engine's — double-counting in the estimator only.

WHY A STATIC TEST.  The failure is invisible to every runtime check available here: no exception, no NaN,
a plausible γ, and a *bit-identical* trajectory that positively argues the swap was clean.  It surfaced
only because an A/B was run and its criterion was declared in advance.  A grep-level invariant is the
cheapest thing that fails the moment someone adds a fourth assembly site and forgets the mask.
"""

from __future__ import annotations

import ast
from pathlib import Path

import pytest

#: `tests/ac/engine/<this>` → parents[3] is the `ffn_sim` package root.  The existence test below is not
#: ceremony: the first draft used parents[4] and every AST assertion passed VACUOUSLY on an empty file list.
DRIVER = Path(__file__).resolve().parents[3] / "scripts" / "ac_gate_b_cortex_motor_native.py"


def _tree() -> ast.Module:
    return ast.parse(DRIVER.read_text(encoding="utf-8"))


def _calls_named(tree: ast.Module, name: str) -> list[ast.Call]:
    return [n for n in ast.walk(tree)
            if isinstance(n, ast.Call) and isinstance(n.func, ast.Name) and n.func.id == name]


def test_the_driver_exists_where_this_test_looks_for_it() -> None:
    """A path that stopped resolving would make every assertion below vacuously pass."""
    assert DRIVER.is_file(), DRIVER


def test_every_accumulate_all_call_passes_an_omit_mask() -> None:
    """The mask is the whole point: an unmasked assembly beside a masked solve reads a different field."""
    unmasked = [c for c in _calls_named(_tree(), "_accumulate_all")
                if not any(k.arg == "omit" for k in c.keywords)]
    assert not unmasked, (
        "_accumulate_all called without omit= at line(s) "
        f"{[c.lineno for c in unmasked]} — if the solve is masked and this is not, the measurement "
        "double-counts whatever a component supplies through a non-omittable hook")


def test_every_measure_gamma_call_passes_an_omit_mask() -> None:
    """Same invariant one level up, so a caller cannot rely on the default and silently drop the mask."""
    unmasked = [c for c in _calls_named(_tree(), "_measure_gamma")
                if not any(k.arg == "omit" for k in c.keywords)]
    assert not unmasked, (
        f"_measure_gamma called without omit= at line(s) {[c.lineno for c in unmasked]}")


def test_the_dump_viz_helper_forwards_a_mask_rather_than_defaulting() -> None:
    """The figure prints γ in its title; a mismatched mask there mislabels a rendered scene."""
    unmasked = [c for c in _calls_named(_tree(), "_dump_viz")
                if not any(k.arg == "omit" for k in c.keywords)]
    assert not unmasked, f"_dump_viz called without omit= at line(s) {[c.lineno for c in unmasked]}"


@pytest.mark.parametrize("func", ["_measure_gamma", "_dump_viz"])
def test_the_helpers_accept_the_mask_at_all(func: str) -> None:
    """Guard against the guard being satisfied by a **kwargs sink that ignores what it is handed."""
    defs = [n for n in ast.walk(_tree())
            if isinstance(n, ast.FunctionDef) and n.name == func]
    assert defs, f"{func} not found — this test's target moved"
    args = defs[0].args
    named = [a.arg for a in list(args.args) + list(args.kwonlyargs)]
    assert "omit" in named, f"{func} does not declare an `omit` parameter; it has {named}"
