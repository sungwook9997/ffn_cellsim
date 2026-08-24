"""Static enforcement for the PI-ratified Warp-CUDA-only runtime contract."""

from __future__ import annotations

import ast
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
# The 2026-08-09 restructure split `ac/` into `engine/` + `components/`. Naming the old path
# here would have scanned nothing and passed — the HOOMD contract is exactly the kind of rule that
# is worth nothing the moment its scan surface goes empty, so `test_the_scan_surface_is_not_empty`
# below asserts these roots hold files.
ACTIVE_ROOTS = (
    REPO_ROOT / "aleph" / "engine",
    REPO_ROOT / "aleph" / "components",
    REPO_ROOT / "aleph" / "laws",
    REPO_ROOT / "aleph" / "tests" / "ac",
)


def test_the_scan_surface_is_not_empty() -> None:
    """A contract test whose roots do not exist passes by scanning nothing."""
    for root in ACTIVE_ROOTS:
        assert root.is_dir(), f"{root} does not exist; this contract would pass vacuously"
        assert any(root.rglob("*.py")), f"{root} holds no Python files"


def _is_hoomd_module(name: str | None) -> bool:
    """Return whether an import name resolves inside the forbidden HOOMD package."""

    return bool(name and (name == "hoomd" or name.startswith("hoomd.")))


def _forbidden_imports(path: Path) -> list[str]:
    """Find direct or literal dynamic HOOMD imports in one Python source file."""

    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    violations: list[str] = []

    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            violations.extend(
                f"{path.relative_to(REPO_ROOT)}:{node.lineno}: import {alias.name}"
                for alias in node.names
                if _is_hoomd_module(alias.name)
            )
        elif isinstance(node, ast.ImportFrom) and _is_hoomd_module(node.module):
            violations.append(
                f"{path.relative_to(REPO_ROOT)}:{node.lineno}: from {node.module} import ..."
            )
        elif isinstance(node, ast.Call) and node.args:
            callee = node.func
            is_dynamic_import = (
                isinstance(callee, ast.Name)
                and callee.id == "__import__"
                or isinstance(callee, ast.Attribute)
                and isinstance(callee.value, ast.Name)
                and callee.value.id == "importlib"
                and callee.attr == "import_module"
            )
            module_arg = node.args[0]
            if (
                is_dynamic_import
                and isinstance(module_arg, ast.Constant)
                and isinstance(module_arg.value, str)
                and _is_hoomd_module(module_arg.value)
            ):
                violations.append(
                    f"{path.relative_to(REPO_ROOT)}:{node.lineno}: dynamic import {module_arg.value}"
                )

    return violations


def test_active_engine_contains_no_hoomd_import() -> None:
    """The active runtime/test scope must never import or dynamically load HOOMD."""

    violations = [
        violation
        for root in ACTIVE_ROOTS
        for path in root.rglob("*.py")
        for violation in _forbidden_imports(path)
    ]
    assert not violations, "Forbidden HOOMD dependency in active scope:\n" + "\n".join(violations)
