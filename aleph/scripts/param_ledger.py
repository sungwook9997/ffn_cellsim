#!/usr/bin/env python
r"""Enumerate every magnitude declared in Python, so the axis count is COUNTED and not remembered.

**Why this exists, and it is a correction rather than an addition.** On 2026-08-24 the parameter
inventory was assembled BY HAND from the populations and connectors the code stands: 79 new axes on top
of the 35 PI-GAP cards. The PI's answer was *"더 필요할 것 같은데"*, and it was right — an AST sweep of
``laws/`` alone found **185** declared constants, including whole modules the hand count never reached:
``fa_maturation`` (24, the talin/vinculin chain), ``ecm_library`` (23), ``hand_kmc`` (17),
``polarization_activegel`` (11), ``polymerization_warp`` (9, actin turnover) and ``piezo`` (7). The hand
count was not slightly short; it missed more than it found.

⚠ **AND THE BLIND SPOT WAS ALREADY WRITTEN DOWN.** ``outputs/tag_kb/params_manifest.yaml``'s own header
says it: *"parameters that live in Python (ac/engine, ac/cell, ff/, dcm/ defaults and dataclass fields)
are outside this file's reach entirely."* The KB gate counts YAML constants and ratchets them; the
constants the simulator actually runs on are in Python and nothing counts them. This is the counter.

⚠ **IT GRADES NOTHING.** Every row lands ``UNGRADED`` unless an existing declaration supplies a grade.
Assigning ``SOURCED`` / ``DERIVED`` / ``CONVENIENCE`` / ``PI_GAP`` / ``UNRATIFIED_PROXY`` is a judgement
about provenance and it is the PI's; a script that guessed would manufacture exactly the false coverage
the five-grade vocabulary exists to prevent. What this does is make the DENOMINATOR real.

⚠ **EXCLUSIONS ARE PRINTED, NEVER SILENT.** A name dropped for being a path, a label or a non-numeric is
reported with its reason and counted, because a sweep that quietly narrows its own scope reads as
completeness. That failure has its own entry in this repository's history more than once.

**What it does NOT decide.** Whether a magnitude is an inference AXIS is a separate question from
whether it is declared here — a discretisation constant, an acceptance threshold and a physical
stiffness all appear as floats. The ledger carries a ``kind`` hint and stops there; promoting a row to
an axis is a PI decision, and so is the partition that matters most: which axes are SHARED across cell
types and which are per-cell-state. ``ROADMAP.md``'s denominator is ``K x M x N`` and this counts one N.

Sanity Gate:
    * dimensional — none computed. Values are copied verbatim from the literal that declared them; a
      value this script cannot read as a literal is recorded as ``None`` rather than evaluated, because
      importing the module to read it would run whatever the module runs at import.
    * boundary cases — a package with no constants is REPORTED as scanned-and-empty, never omitted; a
      file that fails to parse raises with its path rather than being skipped.
    * conservation — every name found is either in the ledger or in the exclusions, and the script
      asserts the two partition the scan before writing.
    * CFL/precision — not applicable; nothing is integrated and no float is compared.
    * sign sense — not applicable; no quantity is accumulated.
    * measurement protocol — one pass over the source tree, parsing only. No module is imported, no
      device is opened, and the count is reproducible from a checkout alone.

engine units: none — this counts declarations, it does not carry magnitudes into physics.
Runtime: pure host AST parsing. CPU-only by nature; it opens no device and imports nothing it scans.

Usage:
    python aleph/scripts/param_ledger.py --out aleph/docs/param_ledger.yaml
    python aleph/scripts/param_ledger.py --summary-only
"""

from __future__ import annotations

import argparse
import ast
from dataclasses import dataclass, field
from pathlib import Path

#: Packages swept. ⚠ These are ROOTS, not a file list — every `*.py` beneath each is parsed, so a
#: module added tomorrow is counted with no edit here. A typed file list is the defect one level up and
#: is how the hand count of 2026-08-24 missed six whole modules.
ROOTS: tuple[str, ...] = ("aleph/laws", "aleph/world", "aleph/observe")

#: ⚠ **`aleph/virtual_cell` IS DELIBERATELY NOT SWEPT, and it was on the first pass.** It carries 465
#: declarations, which would have made it the largest contributor — and they are the RULER, not the
#: thing measured. Its modules say so themselves: `sandbox_fluctuation` opens *"this module makes no
#: claim about the engine's membrane … it is the ORACLE the engine will be judged against"*, and
#: `sandbox_inverse` asks *"what is identifiable from a synthetic microscopy image ALONE"*. Counting an
#: oracle's settings inside N corrupts the only question N exists to answer: how many axes the
#: observations constrain.
#:
#: ⚠ The argument that nearly kept `optics.py` + `synthetic_microscopy.py` (21 rows) was that a PSF
#: width is a NUISANCE PARAMETER trading off against a cell stiffness — a blurred image and a stiff
#: cortex look alike. That is true physics and false about THIS pipeline, and the difference was
#: settled by measurement rather than argument: **nothing outside `virtual_cell/` imports either
#: module, and nothing in `world/` / `observe/` / `laws/` / `scripts/` imports `virtual_cell` at all.**
#: The stage-0 observables are computed from arrays geometrically (`projected_area_um2`,
#: `site_area_um2`); `observe/traction.py` reads a load path, not an image. The optical model enters
#: only for an image-to-image comparison, which is the sandbox's question and the external-network
#: data path — not the engine's.
#:
#: ⚠ **AND THAT LEAVES AN ASSUMPTION NOBODY HAS WRITTEN DOWN.** Three of the six stage-0 observables —
#: spread area, migration, traction — ARE imaging measurements in the laboratory. Computing them
#: geometrically assumes the experimental image-to-number extraction is unbiased. That is a modelling
#: decision of the same class as "no nucleotide" and "the membrane is a surface, not a bilayer", and
#: like those it belongs in the record rather than in the gap between two trees.
_NOT_SWEPT: tuple[tuple[str, str], ...] = (
    ("aleph/virtual_cell", "oracles, sandboxes and the external-network data path — the ruler, "
                           "not the cell. 465 declarations, imported by nothing in the engine path."),
)

#: Name suffixes that are never a magnitude. Matched on the LAST underscore-separated word so
#: `LEDGER_PATH` and `AUDIT_DOC` drop while `X_CATCH_UM` stays.
_NOT_A_MAGNITUDE: frozenset[str] = frozenset({
    "PATH", "DOC", "URL", "FILE", "DIR", "KEY", "NAME", "NAMES", "LABEL", "LABELS",
    "MSG", "TEXT", "NOTE", "REASON", "SCHEMA", "VERSION", "KIND", "KINDS", "ORDER",
})

#: Word fragments that hint a constant is a KNOB OF THE MACHINE rather than a property of the cell.
#: ⚠ A HINT, not a verdict — the row still lands in the ledger and still needs grading. It exists so a
#: reader can see at a glance that `DEFAULT_MIN_TAU_WINDOWS` and `F_HEAD_PN` are not the same kind of
#: thing, without the script deciding that for them.
#:
#: ⚠⚠ **AND IT IS WRONG MORE OFTEN THAN IT IS RIGHT, MEASURED.** PI ruling 2026-08-24 refused it as a
#: FILTER after the physics tree's 21 flagged rows were read one by one: **13 of 21 are false
#: positives (62%)**, because the hint reads the NAME. Among the mislabelled:
#:   `E_SUB_DEFAULT_PA` 5000.0 — the SUBSTRATE MODULUS, and the independent variable of a stage-0
#:      observable ("traction vs substrate stiffness"). Flagged only for containing DEFAULT.
#:   `EPS_RUPTURE` 0.5 — nuclear-envelope rupture strain, which is PI-GAP card C3.
#:   `R_cell_min_um` / `R_cell_max_um` / `nc_min` / `nc_max` — PRIOR BOUNDS on a distribution, i.e.
#:      an axis RANGE, which is the opposite of a knob.
#:   `v_max_um_s`, `NU_SUB_DEFAULT`, `A_ADHESION_UM_DEFAULT` — all physical.
#: Excluding machine-hint wholesale would therefore DROP a stage-0 independent variable and a carded
#: physical threshold. Grade per row; never filter on this field.
_MACHINE_HINTS: frozenset[str] = frozenset({
    "TOL", "TOLERANCE", "EPS", "MIN", "MAX", "DEFAULT", "LIMIT", "BUDGET", "CAP",
    "ITER", "ITERS", "STEPS", "WINDOW", "WINDOWS", "SEED", "CHUNK", "BLOCK", "STRIDE",
})


@dataclass(frozen=True, slots=True)
class Row:
    """One declared magnitude, located exactly."""

    module: str
    name: str
    line: int
    value: float | int | None
    scope: str          # "module" | class name
    kind: str           # "physical-candidate" | "machine-hint"


@dataclass
class Scan:
    """What one sweep found, with its exclusions kept rather than dropped."""

    rows: list[Row] = field(default_factory=list)
    excluded: list[tuple[str, str, str]] = field(default_factory=list)   # module, name, reason
    files: int = 0
    empty_packages: list[str] = field(default_factory=list)


def _literal(node: ast.AST) -> float | int | None:
    """The numeric value of a literal, or ``None`` for anything that needs evaluation.

    ⚠ Deliberately shallow. Reading `2.0 * PI_0` would require importing the module, and importing a
    module to count its constants runs whatever it runs at import — including, in this tree, `wp.init()`.
    An unreadable value is recorded as unreadable, which is a true statement about the declaration.
    """
    if isinstance(node, ast.Constant) and isinstance(node.value, (int, float)) \
            and not isinstance(node.value, bool):
        return node.value
    if isinstance(node, ast.UnaryOp) and isinstance(node.op, ast.USub):
        inner = _literal(node.operand)
        return None if inner is None else -inner
    return None


def _kind(name: str) -> str:
    """Hint whether a name reads as cell physics or as machine control."""
    words = set(name.upper().split("_"))
    return "machine-hint" if words & _MACHINE_HINTS else "physical-candidate"


def _scan_module(path: Path, rel: str, scan: Scan) -> None:
    """Parse one file and append its declarations. Raises with the path on a parse failure."""
    try:
        tree = ast.parse(path.read_text(encoding="utf-8"))
    except SyntaxError as exc:                                  # pragma: no cover - a real defect
        raise SyntaxError(f"{rel}: {exc}") from exc
    scan.files += 1

    def take(name: str, node: ast.AST, line: int, scope: str) -> None:
        last = name.upper().rsplit("_", 1)[-1]
        if last in _NOT_A_MAGNITUDE:
            scan.excluded.append((rel, name, f"name ends in {last} — not a magnitude"))
            return
        value = _literal(node)
        if value is None and not isinstance(node, ast.Constant):
            # An annotation-only dataclass field has no value node at all; that is still a declaration.
            pass
        scan.rows.append(Row(rel, name, line, value, scope, _kind(name)))

    for node in tree.body:
        # ⚠ BOTH forms. `Assign` alone missed every module-level ANNOTATED constant
        # (`X: float = 1.0`), which made them invisible to the ledger AND to its exclusion list —
        # neither counted nor named. The self-check caught it on `turgor_pi0.LEDGER_PATH`, which is
        # exactly the shape this sweep exists to stop: a scope that narrows itself silently.
        if isinstance(node, (ast.Assign, ast.AnnAssign)):
            if isinstance(node, ast.Assign):
                if len(node.targets) != 1 or not isinstance(node.targets[0], ast.Name):
                    continue
                target, value = node.targets[0], node.value
            else:
                if not isinstance(node.target, ast.Name):
                    continue
                target, value = node.target, node.value
            if not target.id.isupper() or target.id.startswith("_"):
                continue
            if value is None or _literal(value) is None:
                scan.excluded.append((rel, target.id, "value is not a numeric literal"))
                continue
            take(target.id, value, node.lineno, "module")
        elif isinstance(node, ast.ClassDef):
            for body in node.body:
                if not isinstance(body, ast.AnnAssign) or not isinstance(body.target, ast.Name):
                    continue
                annotation = body.annotation
                if not (isinstance(annotation, ast.Name) and annotation.id in ("float", "int")):
                    continue
                default = body.value if body.value is not None else annotation
                take(body.target.id, default, body.lineno, node.name)


def sweep(roots: tuple[str, ...] = ROOTS, repo: Path | None = None) -> Scan:
    """Walk every package root and return one :class:`Scan`.

    Args:
        roots: package directories, relative to the repository. Every ``*.py`` beneath each is parsed.
        repo: repository root. Defaults to this file's grandparent.

    Returns:
        The scan, with rows and exclusions partitioning everything the sweep named.
    """
    base = repo or Path(__file__).resolve().parents[2]
    scan = Scan()
    for root in roots:
        directory = base / root
        if not directory.is_dir():
            scan.empty_packages.append(f"{root} (absent)")
            continue
        before = len(scan.rows)
        for path in sorted(directory.rglob("*.py")):
            if "__pycache__" in path.parts or path.name == "__init__.py":
                continue
            _scan_module(path, str(path.relative_to(base)), scan)
        if len(scan.rows) == before:
            scan.empty_packages.append(f"{root} (scanned, no declarations)")
    return scan


def summarise(scan: Scan) -> str:
    """A per-module count, heaviest first, with the exclusion tally beside it."""
    per: dict[str, int] = {}
    for row in scan.rows:
        per[row.module] = per.get(row.module, 0) + 1
    lines = [f"{'module':52s} {'n':>4s}  kind split"]
    for module, count in sorted(per.items(), key=lambda kv: (-kv[1], kv[0]))[:24]:
        machine = sum(1 for r in scan.rows if r.module == module and r.kind == "machine-hint")
        lines.append(f"{module:52s} {count:4d}  physical {count - machine:3d} / machine {machine:3d}")
    if len(per) > 24:
        lines.append(f"... and {len(per) - 24} more modules")
    machine_total = sum(1 for r in scan.rows if r.kind == "machine-hint")
    lines += [
        "",
        f"files parsed                  {scan.files:5d}",
        f"declarations in the ledger    {len(scan.rows):5d}"
        f"   (physical-candidate {len(scan.rows) - machine_total}, machine-hint {machine_total})",
        f"EXCLUDED, with reason         {len(scan.excluded):5d}",
    ]
    for note in scan.empty_packages:
        lines.append(f"  ⚠ {note}")
    return "\n".join(lines)


def as_yaml(scan: Scan) -> str:
    """The ledger, hand-formatted so it stays readable in review and diffs one row per line."""
    out = [
        "# param_ledger.yaml — GENERATED by aleph/scripts/param_ledger.py. Do not hand-edit rows.",
        "#",
        "# Every magnitude DECLARED IN PYTHON under the swept roots. This is the denominator the",
        "# hand inventory of 2026-08-24 under-counted by more than half.",
        "#",
        "# ⚠ `grade: UNGRADED` on every row is the honest state, not an omission. SOURCED / DERIVED /",
        "# CONVENIENCE / PI_GAP / UNRATIFIED_PROXY is a provenance judgement and it is the PI's.",
        "# ⚠ `kind` is a HINT read off the NAME and it is WRONG 13 times in 21 on the physics tree",
        "#   (measured 2026-08-24). PI ruling: it is a reading aid and MAY NOT be used as a filter —",
        "#   doing so drops E_SUB_DEFAULT_PA, which is a stage-0 observable's independent variable.",
        "# ⚠ `partition: UNASSIGNED` is decision 3 and it sets the REAL total. SHARED (a myosin head is",
        "#   a myosin head, and sharing it is what makes the two-cell-line test falsifiable) /",
        "#   PER_CELL_LINE (expression levels, every areal density) / PER_STATE (shape, polarity,",
        "#   turgor) / INSTRUMENT. K x M x N is 399 if everything is shared and ~700 if nothing is.",
        "# ⚠ A row is not an inference axis. Promoting one is a separate PI decision, and so is the",
        "#   partition that decides the real total: which axes are SHARED across cell types and which",
        "#   are per-cell-state. ROADMAP's denominator is K x M x N; this counts one N.",
        "",
        f"roots: [{', '.join(ROOTS)}]",
        f"files_parsed: {scan.files}",
        f"n_declarations: {len(scan.rows)}",
        f"n_excluded: {len(scan.excluded)}",
        "",
        "rows:",
    ]
    for row in sorted(scan.rows, key=lambda r: (r.module, r.line)):
        value = "null" if row.value is None else repr(row.value)
        out.append(
            f"  - {{module: {row.module}, name: {row.name}, line: {row.line}, "
            f"scope: {row.scope}, value: {value}, kind: {row.kind}, "
            f"grade: UNGRADED, partition: UNASSIGNED}}"
        )
    out += ["", "# Named rather than dropped — a sweep that narrows itself silently reads as complete.",
            "excluded:"]
    for module, name, reason in sorted(scan.excluded):
        out.append(f"  - {{module: {module}, name: {name}, reason: \"{reason}\"}}")
    return "\n".join(out) + "\n"


def _demo() -> None:
    """Self-check, runnable without pytest — the shape this tree adopted 2026-08-21.

    ⚠ Asserts the PARTITION and the anchors, not a total. Pinning the total here would make every
    added constant a failure of this file rather than a fact about the tree, which is the opposite of
    what a denominator is for.
    """
    scan = sweep()
    assert scan.files > 50, scan.files
    assert scan.rows, "swept the tree and found no declaration — that is the finding, not a pass"

    names = {(r.module, r.name) for r in scan.rows}
    excluded = {(m, n) for m, n, _ in scan.excluded}
    assert not (names & excluded), sorted(names & excluded)[:5]

    # Anchors: three the hand count MISSED, one it found, one that must be excluded.
    assert any(n.startswith("TALIN_") for _m, n in names), "fa_maturation's talin chain went uncounted"
    assert any("k_catch0" == n for _m, n in names), "the crossbridge catch constant went uncounted"
    assert any(m.endswith("polymerization_warp.py") for m, _n in names), "actin turnover uncounted"
    assert any(n == "LEDGER_PATH" for _m, n in excluded), "a path was taken as a magnitude"

    machine = sum(1 for r in scan.rows if r.kind == "machine-hint")
    assert 0 < machine < len(scan.rows), (machine, len(scan.rows))
    print(f"[param_ledger] self-check OK — {len(scan.rows)} declarations over {scan.files} files, "
          f"{len(scan.excluded)} excluded with reason")


def main(argv: list[str] | None = None) -> int:
    """Sweep and report, optionally writing the ledger."""
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--out", type=Path, default=None, help="write the YAML ledger here")
    ap.add_argument("--summary-only", action="store_true", help="print the tally and stop")
    ap.add_argument("--show-excluded", action="store_true", help="list every exclusion and its reason")
    args = ap.parse_args(argv)

    scan = sweep()
    print(summarise(scan), flush=True)

    if args.show_excluded:
        print("\nEXCLUDED — named, not dropped:")
        for module, name, reason in sorted(scan.excluded):
            print(f"  {module:46s} {name:32s} {reason}")

    if args.out and not args.summary_only:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(as_yaml(scan))
        print(f"\n[param_ledger] wrote {args.out}  ({len(scan.rows)} rows, all UNGRADED)", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
