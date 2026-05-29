#!/usr/bin/env python
"""Pre-run module-integrity guard for ffn_sim production drivers.

Syncthing (the gbook sync path) occasionally leaves a runtime module as a
stale, zero-byte copy, or drops a ``*.sync-conflict-*`` shadow file beside the
real one.  Either can silently break a long unattended run — a zero-byte module
imports as an empty namespace (``AttributeError`` deep in the physics), and a
shadow file is at best clutter and at worst imported by mistake.

This module performs a cheap pre-flight check that a driver runs *before* it
spends CPU on a build:

#. **Import** every ffn_sim runtime submodule (ecm / cortex / bridge / cell /
   integrator / common / junction) and report any failure.
#. **Scan** the ``ffn_sim/`` tree for

   * zero-byte ``.py`` files — *excluding* ``__init__.py`` (legitimately empty
     package markers are never flagged), and
   * ``*.sync-conflict-*`` Syncthing shadow files.

Severity policy
---------------
* Import failures and zero-byte source modules are **fatal**.
* ``*.sync-conflict-*`` files are **fatal only under source directories**
  (``ecm/ cortex/ bridge/ cell/ integrator/ common/ junction/ scripts/``); the
  same shadow under ``ffn_sim/outputs/`` is a **warning** (cosmetic data
  clutter, cannot affect import).

Callers typically invoke :func:`assert_clean` (raises :class:`IntegrityError`
on any fatal finding) best-effort, downgrading the exception to a warning so a
genuinely valid run is never aborted by a cosmetic issue.
"""
from __future__ import annotations

import importlib
import logging
from dataclasses import dataclass, field
from pathlib import Path

LOGGER = logging.getLogger("ffn_sim.common.integrity")

#: Package root (``ffn_sim/``) inferred from this file's location.
_PKG_ROOT = Path(__file__).resolve().parents[1]

#: Source directories whose contents are import-critical; a sync-conflict or
#: zero-byte module here is fatal.
SOURCE_DIRS: tuple[str, ...] = (
    "ecm",
    "cortex",
    "bridge",
    "cell",
    "integrator",
    "common",
    "junction",
    "scripts",
)

#: Subpackages whose modules are imported as the runtime importability probe.
#: ``scripts`` is intentionally excluded — drivers have heavy side-effecting
#: ``__main__`` blocks and argparse, so importing them all is neither cheap nor
#: meaningful; the file-scan still covers them.
_RUNTIME_SUBPACKAGES: tuple[str, ...] = (
    "ecm",
    "cortex",
    "bridge",
    "cell",
    "integrator",
    "common",
    "junction",
)


class IntegrityError(RuntimeError):
    """Raised by :func:`assert_clean` when a fatal integrity issue is found."""


@dataclass
class IntegrityReport:
    """Structured result of a :func:`check_module_integrity` scan.

    Attributes:
        import_failures: ``(module, error-message)`` pairs for modules that
            failed to import.
        zero_byte: Paths of zero-byte ``.py`` files (excluding ``__init__.py``).
        sync_conflicts_source: ``*.sync-conflict-*`` paths under source dirs
            (fatal).
        sync_conflicts_output: ``*.sync-conflict-*`` paths under ``outputs/``
            (warning only).
    """

    import_failures: list[tuple[str, str]] = field(default_factory=list)
    zero_byte: list[Path] = field(default_factory=list)
    sync_conflicts_source: list[Path] = field(default_factory=list)
    sync_conflicts_output: list[Path] = field(default_factory=list)

    @property
    def fatal(self) -> bool:
        """``True`` if any fatal-severity issue was recorded."""
        return bool(
            self.import_failures
            or self.zero_byte
            or self.sync_conflicts_source
        )

    def summary(self) -> str:
        """Return a human-readable multi-line description of all findings."""
        lines: list[str] = []
        for mod, err in self.import_failures:
            lines.append(f"  [import]  {mod}: {err}")
        for p in self.zero_byte:
            lines.append(f"  [0-byte]  {p}")
        for p in self.sync_conflicts_source:
            lines.append(f"  [conflict/source]  {p}")
        for p in self.sync_conflicts_output:
            lines.append(f"  [conflict/output:WARN]  {p}")
        if not lines:
            return "module integrity: clean"
        return "module integrity issues:\n" + "\n".join(lines)


def _is_under(path: Path, root: Path) -> bool:
    """Return ``True`` if ``path`` is ``root`` or nested beneath it.

    Args:
        path: Candidate path (assumed resolved).
        root: Directory to test containment against (assumed resolved).
    """
    try:
        path.relative_to(root)
        return True
    except ValueError:
        return False


def _scan_files(report: IntegrityReport) -> None:
    """Populate ``report`` with file-system findings under :data:`_PKG_ROOT`.

    Args:
        report: Report object to mutate in place.
    """
    outputs_root = (_PKG_ROOT / "outputs").resolve()
    source_roots = [(_PKG_ROOT / d).resolve() for d in SOURCE_DIRS]

    # Zero-byte .py files (excluding __init__.py package markers).
    for py in sorted(_PKG_ROOT.rglob("*.py")):
        if "__pycache__" in py.parts:
            continue
        if py.name == "__init__.py":
            continue
        try:
            if py.stat().st_size == 0:
                report.zero_byte.append(py)
        except OSError:  # pragma: no cover - racing fs
            continue

    # Syncthing shadow files: *.sync-conflict-*  (any extension).
    for shadow in sorted(_PKG_ROOT.rglob("*.sync-conflict-*")):
        if "__pycache__" in shadow.parts:
            continue
        rp = shadow.resolve()
        if _is_under(rp, outputs_root):
            report.sync_conflicts_output.append(shadow)
        elif any(_is_under(rp, sr) for sr in source_roots):
            report.sync_conflicts_source.append(shadow)
        else:
            # Under the package tree but outside known source / outputs roots
            # (e.g. docs/, configs/) — treat conservatively as a warning.
            report.sync_conflicts_output.append(shadow)


def _scan_imports(report: IntegrityReport) -> None:
    """Populate ``report`` with import-failure findings.

    Imports every ``*.py`` module (excluding dunder files) in the runtime
    subpackages; records the first error per module.

    Args:
        report: Report object to mutate in place.
    """
    for sub in _RUNTIME_SUBPACKAGES:
        subdir = _PKG_ROOT / sub
        if not subdir.is_dir():
            continue
        for py in sorted(subdir.glob("*.py")):
            if py.name.startswith("__"):
                continue
            mod_name = f"ffn_sim.{sub}.{py.stem}"
            try:
                importlib.import_module(mod_name)
            except Exception as exc:  # noqa: BLE001  (report, do not raise)
                report.import_failures.append((mod_name, repr(exc)))


def check_module_integrity(*, scan_imports: bool = True) -> IntegrityReport:
    """Run the full integrity scan and return a structured report.

    Args:
        scan_imports: If ``True`` (default), additionally attempt to import
            every runtime submodule.  Set ``False`` to limit the check to the
            cheap file-system scan (no import side effects).

    Returns:
        An :class:`IntegrityReport`.  Use :attr:`IntegrityReport.fatal` to test
        severity or :func:`assert_clean` to raise on fatal findings.
    """
    report = IntegrityReport()
    _scan_files(report)
    if scan_imports:
        _scan_imports(report)
    return report


def assert_clean(*, scan_imports: bool = True) -> IntegrityReport:
    """Run :func:`check_module_integrity` and raise on any fatal finding.

    Output-directory sync-conflict files are logged as a warning but do not
    raise.  All fatal findings (import failures, zero-byte source modules,
    source-dir sync-conflicts) are collected into one clear error message.

    Args:
        scan_imports: Forwarded to :func:`check_module_integrity`.

    Returns:
        The (clean-enough-to-proceed) :class:`IntegrityReport` when no fatal
        issue is present.

    Raises:
        IntegrityError: If any fatal-severity issue is detected.
    """
    report = check_module_integrity(scan_imports=scan_imports)

    for p in report.sync_conflicts_output:
        LOGGER.warning("cosmetic sync-conflict under outputs/: %s", p)

    if report.fatal:
        raise IntegrityError(report.summary())

    LOGGER.info("module integrity check passed")
    return report
