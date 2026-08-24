#!/usr/bin/env python3
"""Manifest COVERAGE ratchet — the OMISSION side of the KB integrity gates.

WHY THIS EXISTS (2026-07-25, audit finding). `verify_params.py` and
`verify_runs.py` check every row a manifest *lists*. Neither could ever be wrong
about what a manifest **OMITS**. Measured consequence, both gates green:

  * `params_manifest.yaml` declared 42 constants — ALL of them under
    `aleph/validation/oracles/configs/phase1_unit*.yaml`, retired v1 oracle
    configs that no runtime path reads. The live runtime config tree
    `aleph/configs/**` carries 71 KU-tagged numeric constants, of which the
    manifest covered **zero**.
  * `results_manifest.yaml` declared 31 result-claims, none of them from `ff/`
    or `ac/`; 36 of the 37 committed `aleph/outputs/**/REPORT.md` files were
    referenced nowhere in it.

So both manifests froze in June while `make kb-check` stayed green. A gate that
cannot notice a whole engine is missing manufactures confidence — the exact
failure mode the audit named as worse than having no gate.

HOW IT WORKS — a RATCHET, not a cliff. 71 params + 29 REPORT headlines will not
be declared today, and a gate that instantly red-lines everything gets disabled
by the next session. So the CURRENT undeclared set is written down, item by item
with a one-line reason, in the checked-in `coverage_baseline.yaml`, and this
module gates that the set may only SHRINK:

  ok        item declared in the manifest                                (best)
  DEBT      item uncovered AND listed in coverage_baseline.yaml — counted and
            printed on every single run so it cannot be forgotten
  COVERAGE_GAP        uncovered and NOT in the baseline -> HARD FAIL
                      (a new KU-tagged constant / a new REPORT.md headline)
  FINGERPRINT_CHANGED an item in the baseline whose on-disk value (params) or
                      headline-number set (results) CHANGED -> HARD FAIL.
                      The baseline blesses a specific undeclared value, never a
                      licence for that value to drift silently.
  RATCHET_STALE       baseline entry that is now covered, or gone from disk ->
                      HARD FAIL until pruned. The baseline must equal the
                      uncovered set exactly; slack in it is a suppression list.
                      Fix with `--prune-baseline` (removes only; never adds).

Anti-erosion properties, deliberate:
  * a MISSING or unimportable coverage module makes the calling gate FAIL, never
    skip (see verify_params.py / verify_runs.py import blocks);
  * a MISSING baseline file makes every uncovered item a new gap -> loud fail;
  * an EMPTY config tree makes every baseline entry stale -> loud fail. You
    cannot go green by deleting the thing being measured;
  * `--prune-baseline` can only REMOVE entries. Widening the ratchet is
    `--add-baseline <id> --reason "<why>"`, which requires an explicit reason,
    stamps `added:` with today's date, and is reported separately from the
    legacy 2026-07-25 debt on every run.

WHAT COUNTS (mechanical definitions, so the inventory cannot be argued with):
  params  — a line in `aleph/configs/**/*.y{a}ml` that (a) is a `key: value`
            mapping entry whose value parses as a number and (b) mentions
            `KU-<digit>` anywhere on the line (i.e. the project's own
            "this constant is literature-anchored" marker).
            id = "<repo-relative config path>:<dotted.key>" — the exact pair
            params_manifest.yaml already uses, so declaring an item covers it.
            KU-tagged lines whose value is a list/bool/string (acceptance bands,
            on/off flags) are NOT coverage-eligible; they are counted and listed
            so the exclusion is visible rather than silent.
  results — a `aleph/outputs/**/REPORT.md` whose HEADLINE BLOCK (title line
            through the line before the first `##` section, capped at 80 lines)
            contains at least one numeric literal after masking inline-code
            spans, link targets, URLs, ISO dates, bare years and section refs
            (which strips commit hashes, filenames and §refs, leaving the prose
            claim numbers). id = the repo-relative REPORT.md path; covered when
            a results_manifest claim names it as `artifact:` or lists it under
            an optional `covers:` key.

Pure stdlib + PyYAML, sub-second, no network, no token, no duckdb — same
character as the gates it plugs into.

RUN:
  python kb_coverage.py                 # both inventories + the debt counts
  python kb_coverage.py --gate          # both ratchets; exit 1 on any hard fail
  python kb_coverage.py --prune-baseline           # remove stale entries only
  python kb_coverage.py --add-baseline <id> --reason "<why>"   # widen (loud)
"""
from __future__ import annotations

import hashlib
import re
import sys
from dataclasses import dataclass, field
from datetime import date
from pathlib import Path

import yaml

HERE = Path(__file__).parent
REPO_ROOT = HERE.parents[2]  # tag_kb -> outputs -> ffn_sim -> repo root
BASELINE = HERE / "coverage_baseline.yaml"
PARAMS_MANIFEST = HERE / "params_manifest.yaml"
RESULTS_MANIFEST = HERE / "results_manifest.yaml"

CONFIG_ROOT = REPO_ROOT / "aleph" / "configs"
OUTPUTS_ROOT = REPO_ROOT / "aleph" / "outputs"

RATCHET_DATE = "2026-07-25"  # the day the omission side was first written down

# --------------------------------------------------------------------------- #
# inventory: what EXISTS on disk that a manifest ought to declare
# --------------------------------------------------------------------------- #
_KEY_RE = re.compile(r"^(?P<ind>[ \t]*)(?P<key>[A-Za-z_][A-Za-z0-9_.\-]*)[ \t]*:(?P<rest>.*)$")
_LIST_RE = re.compile(r"^(?P<ind>[ \t]*)-[ \t]*(?P<rest>.*)$")
_KU_RE = re.compile(r"KU-\d")
_NOT_A_NUMBER = {"true", "false", "null", "~", "yes", "no", "on", "off"}


@dataclass(frozen=True)
class Item:
    """One thing on disk that a manifest ought to declare."""

    item_id: str
    fingerprint: str   # value (params) / headline-number sha (results)
    detail: str        # human-readable: unit-ish comment, or the number list


def _strip_comment(rest: str) -> str:
    """Drop a trailing YAML `#` comment, respecting simple quoting."""
    out: list[str] = []
    quote: str | None = None
    for i, ch in enumerate(rest):
        if quote:
            out.append(ch)
            if ch == quote:
                quote = None
            continue
        if ch in "\"'":
            quote = ch
            out.append(ch)
            continue
        if ch == "#" and (i == 0 or rest[i - 1] in " \t"):
            break
        out.append(ch)
    return "".join(out).strip()


def _as_number(value: str) -> float | None:
    if value == "" or value.lower() in _NOT_A_NUMBER:
        return None
    try:
        return float(value)
    except ValueError:
        return None


def _fmt_value(v: float) -> str:
    """Round-trippable float text. Quoted in the baseline YAML on purpose: PyYAML
    1.1 will not resolve `1.0e-6` (no exponent sign) as a float, so the baseline
    stores the literal and this module does the float() itself."""
    return repr(v)


def scan_config(path: Path) -> tuple[list[Item], list[Item], int]:
    """Scan one config YAML line-wise (comments carry the KU tags, so a parsed
    tree is useless here).

    Returns (ku_numeric_items, ku_non_numeric_items, ku_mentioning_lines).
    """
    rel = path.relative_to(REPO_ROOT).as_posix()
    numeric: list[Item] = []
    other: list[Item] = []
    mentions = 0
    stack: list[tuple[int, str]] = []
    for line in path.read_text().splitlines():
        if _KU_RE.search(line):
            mentions += 1
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        in_list = False
        m = _KEY_RE.match(line)
        if m:
            indent = len(m.group("ind").expandtabs(2))
        else:
            lm = _LIST_RE.match(line)
            if not lm:
                continue
            m = _KEY_RE.match(lm.group("rest"))
            if not m:
                continue
            in_list = True
            indent = len(lm.group("ind").expandtabs(2)) + 2
        key = m.group("key")
        value = _strip_comment(m.group("rest"))
        while stack and stack[-1][0] >= indent:
            stack.pop()
        dotted = ".".join([k for _, k in stack] + (["-"] if in_list else []) + [key])
        if value == "":                      # parent mapping key, not a constant
            stack.append((indent, key))
            continue
        if not _KU_RE.search(line):
            continue
        item_id = f"{rel}:{dotted}"
        num = _as_number(value)
        comment = line.split("#", 1)[1].strip() if "#" in line else ""
        if num is None:
            other.append(Item(item_id, value, comment or value))
        else:
            numeric.append(Item(item_id, _fmt_value(num), comment or value))
    return numeric, other, mentions


def params_inventory() -> tuple[list[Item], list[Item], int]:
    """Every KU-tagged numeric constant in the LIVE runtime config tree."""
    numeric: list[Item] = []
    other: list[Item] = []
    mentions = 0
    if not CONFIG_ROOT.is_dir():
        return numeric, other, mentions
    for path in sorted(CONFIG_ROOT.rglob("*.y*ml")):
        if path.suffix not in (".yaml", ".yml"):
            continue
        n, o, mm = scan_config(path)
        numeric += n
        other += o
        mentions += mm
    return sorted(numeric, key=lambda i: i.item_id), sorted(other, key=lambda i: i.item_id), mentions


_NUM_RE = re.compile(r"(?<![\w.])[-+]?(?:\d+\.\d+|\.\d+|\d+)(?:[eE][-+]?\d+)?")


def headline_block(text: str, max_lines: int = 80) -> str:
    """Title line through the line before the first `##` section heading."""
    out: list[str] = []
    for i, line in enumerate(text.splitlines()):
        if i and re.match(r"^#{2,}\s", line):
            break
        if i >= max_lines:
            break
        out.append(line)
    return "\n".join(out)


def headline_numbers(block: str) -> list[str]:
    """Numeric literals a reader would take as the report's claim.

    Masked out first (they are identifiers, not claims): inline-code spans
    (commit hashes / filenames / config keys), markdown link targets, autolinks
    and bare URLs, ISO dates, bare 19xx/20xx years, and `§n` section refs.
    """
    b = re.sub(r"`[^`]*`", " ", block)
    b = re.sub(r"\]\([^)]*\)", " ", b)
    b = re.sub(r"<[^>\s]*>", " ", b)
    b = re.sub(r"https?://\S+", " ", b)
    b = re.sub(r"\d{4}-\d{2}-\d{2}", " ", b)
    b = re.sub(r"(?<![\w.])(?:19|20)\d{2}(?![\d.])", " ", b)
    b = re.sub(r"§\s*\d+(?:\.\d+)*", " ", b)
    seen: list[str] = []
    for m in _NUM_RE.finditer(b):
        tok = m.group(0)
        if tok not in seen:
            seen.append(tok)
    return sorted(seen)


def results_inventory() -> tuple[list[Item], list[str]]:
    """Every committed REPORT.md whose headline block states a number.

    Returns (eligible_items, reports_with_no_headline_number).
    """
    items: list[Item] = []
    silent: list[str] = []
    if not OUTPUTS_ROOT.is_dir():
        return items, silent
    for path in sorted(OUTPUTS_ROOT.rglob("REPORT.md")):
        rel = path.relative_to(REPO_ROOT).as_posix()
        nums = headline_numbers(headline_block(path.read_text(errors="replace")))
        if not nums:
            silent.append(rel)
            continue
        sha = hashlib.sha1(",".join(nums).encode()).hexdigest()[:12]
        items.append(Item(rel, sha, ", ".join(nums)))
    return items, silent


# --------------------------------------------------------------------------- #
# what the manifests DECLARE
# --------------------------------------------------------------------------- #
def _claims(path: Path) -> list[dict]:
    if not path.exists():
        return []
    return (yaml.safe_load(path.read_text()) or {}).get("claims", []) or []


def declared_params() -> tuple[set[str], list[str]]:
    """(declared item ids, manifest rows whose config is NOT the runtime tree)."""
    ids: set[str] = set()
    off_tree: list[str] = []
    for c in _claims(PARAMS_MANIFEST):
        cfg = (c.get("config") or "").strip()
        key = (c.get("key") or "").strip()
        if not cfg or not key:
            continue
        ids.add(f"{cfg}:{key}")
        if not cfg.startswith("aleph/configs/"):
            off_tree.append(f"{c.get('id', '?')} -> {cfg}")
    return ids, off_tree


def declared_results() -> set[str]:
    ids: set[str] = set()
    for c in _claims(RESULTS_MANIFEST):
        art = (c.get("artifact") or "").strip()
        if art.endswith("REPORT.md"):
            ids.add(art)
        for extra in c.get("covers") or []:
            ids.add(str(extra).strip())
    return ids


# --------------------------------------------------------------------------- #
# the ratchet baseline
# --------------------------------------------------------------------------- #
def load_baseline() -> dict[str, dict[str, dict]]:
    """{'params': {id: entry}, 'results': {id: entry}}. A missing file is an
    EMPTY baseline on purpose: every uncovered item then reads as a new gap and
    the gate goes loudly red, instead of quietly green."""
    out: dict[str, dict[str, dict]] = {"params": {}, "results": {}}
    if not BASELINE.exists():
        return out
    raw = yaml.safe_load(BASELINE.read_text()) or {}
    for kind in ("params", "results"):
        for entry in raw.get(kind) or []:
            if isinstance(entry, dict) and entry.get("id"):
                out[kind][str(entry["id"])] = entry
    return out


_BASELINE_HEADER = f"""\
# coverage_baseline.yaml — the RATCHET DEBT LIST of the KB integrity gates.
#
# NOT a suppression file. This is the written-down, countable, embarrassing set
# of things that EXIST on disk and are NOT declared in a manifest. kb_coverage.py
# prints the remaining count on every `make kb-check`, and gates that the set may
# only SHRINK:
#
#   * a NEW KU-tagged constant in aleph/configs/** or a NEW REPORT.md headline
#     that is not declared and not listed here  -> HARD FAIL (COVERAGE_GAP)
#   * an entry here whose on-disk value / headline numbers CHANGED -> HARD FAIL
#     (the entry blesses one specific undeclared value, not future drift)
#   * an entry here that is now declared, or gone from disk -> HARD FAIL until
#     pruned (`python kb_coverage.py --prune-baseline`), so this list stays
#     exactly equal to the real uncovered set
#
# THE FIX for any line below is to DELETE it by declaring the item properly:
#   params  -> add a row to params_manifest.yaml (config + key + value + KU +
#              citation_key + honest `declared:`)
#   results -> add a claim to results_manifest.yaml naming the REPORT.md as
#              `artifact:` (or listing it under an existing claim's `covers:`)
#
# Fields: id | value (params: exact on-disk value, quoted so PyYAML 1.1 cannot
#         mis-resolve it) | headline_sha + numbers (results) | reason | added
#
# GENERATED once from disk on {RATCHET_DATE} (`--prune-baseline` rewrites it).
# `added: {RATCHET_DATE}` = pre-existing legacy debt. Any LATER `added:` date is
# debt accepted after the ratchet landed and is counted separately on every run.
"""


def write_baseline(baseline: dict[str, dict[str, dict]]) -> None:
    """Serialize the baseline deterministically (id-sorted, quoted scalars)."""
    lines = [_BASELINE_HEADER]
    for kind in ("params", "results"):
        lines.append(f"{kind}:")
        entries = baseline.get(kind) or {}
        if not entries:
            lines[-1] = f"{kind}: []"
            lines.append("")
            continue
        for item_id in sorted(entries):
            e = entries[item_id]
            lines.append(f'  - id: "{item_id}"')
            if kind == "params":
                lines.append(f'    value: "{e.get("value", "")}"')
            else:
                lines.append(f'    headline_sha: "{e.get("headline_sha", "")}"')
                lines.append(f'    numbers: "{e.get("numbers", "")}"')
            reason = str(e.get("reason", "")).replace('"', "'")
            lines.append(f'    reason: "{reason}"')
            lines.append(f'    added: "{e.get("added", RATCHET_DATE)}"')
        lines.append("")
    BASELINE.write_text("\n".join(lines))


def _entry_fingerprint(kind: str, entry: dict) -> str:
    return str(entry.get("value" if kind == "params" else "headline_sha", ""))


def _same_fingerprint(kind: str, item: Item, entry: dict) -> bool:
    recorded = _entry_fingerprint(kind, entry)
    if kind != "params":
        return recorded == item.fingerprint
    try:                                  # numeric compare: text form may differ
        return float(recorded) == float(item.fingerprint)
    except (TypeError, ValueError):
        return recorded == item.fingerprint


# --------------------------------------------------------------------------- #
# evaluation
# --------------------------------------------------------------------------- #
@dataclass
class Coverage:
    kind: str
    inventory: list[Item] = field(default_factory=list)
    declared_ids: set[str] = field(default_factory=set)
    covered: list[Item] = field(default_factory=list)
    debt_legacy: list[Item] = field(default_factory=list)
    debt_accepted: list[tuple[Item, str]] = field(default_factory=list)
    gaps: list[Item] = field(default_factory=list)            # HARD FAIL
    changed: list[tuple[Item, str]] = field(default_factory=list)  # HARD FAIL
    stale: list[str] = field(default_factory=list)            # HARD FAIL
    notes: list[str] = field(default_factory=list)            # informational

    @property
    def uncovered_n(self) -> int:
        return len(self.debt_legacy) + len(self.debt_accepted) + len(self.gaps) + len(self.changed)

    @property
    def failed(self) -> bool:
        return bool(self.gaps or self.changed or self.stale)


def evaluate(kind: str) -> Coverage:
    """Compare disk inventory vs manifest declarations vs the ratchet baseline."""
    cov = Coverage(kind=kind)
    baseline = load_baseline()[kind]
    if kind == "params":
        inv, non_numeric, mentions = params_inventory()
        declared, off_tree = declared_params()
        cov.notes.append(
            f"{len(non_numeric)} KU-tagged non-scalar lines (acceptance bands / on-off "
            f"flags) are not coverage-eligible")
        cov.notes.append(
            f"{mentions - len(inv) - len(non_numeric)} further KU mentions in "
            f"comments/prose are not `key: value` constants")
        if off_tree:
            cov.notes.append(
                f"{len(off_tree)}/{len(declared)} declared rows point OUTSIDE "
                f"aleph/configs/** (retired oracle configs, read by no runtime path)")
    else:
        inv, silent = results_inventory()
        declared = declared_results()
        cov.notes.append(
            f"{len(silent)} REPORT.md files state no headline number "
            f"(not coverage-eligible; adding one makes them a new gap)")
    cov.inventory, cov.declared_ids = inv, declared

    for item in inv:
        if item.item_id in declared:
            cov.covered.append(item)
            continue
        entry = baseline.get(item.item_id)
        if entry is None:
            cov.gaps.append(item)
        elif not _same_fingerprint(kind, item, entry):
            cov.changed.append((item, _entry_fingerprint(kind, entry)))
        elif str(entry.get("added", RATCHET_DATE)) == RATCHET_DATE:
            cov.debt_legacy.append(item)
        else:
            cov.debt_accepted.append((item, str(entry.get("added"))))

    live = {i.item_id for i in inv}
    for item_id in sorted(baseline):
        if item_id in declared:
            cov.stale.append(f"{item_id}  (now DECLARED — delete this debt line)")
        elif item_id not in live:
            cov.stale.append(f"{item_id}  (gone from disk / no longer KU-tagged)")
    return cov


# --------------------------------------------------------------------------- #
# reporting
# --------------------------------------------------------------------------- #
_WHAT = {"params": ("KU-tagged runtime constants in aleph/configs/**",
                    "params_manifest.yaml"),
         "results": ("REPORT.md headline claims in aleph/outputs/**",
                     "results_manifest.yaml")}


def summary_lines(cov: Coverage, tag: str) -> list[str]:
    what, manifest = _WHAT[cov.kind]
    out = [f"  [{tag}] coverage ratchet — {len(cov.inventory)} {what}",
           f"           declared in {manifest}: {len(cov.covered)}",
           f"           UNDECLARED (ratchet debt): {cov.uncovered_n}"
           f"   <-- may only SHRINK"]
    out.append(f"             legacy (added {RATCHET_DATE}): {len(cov.debt_legacy)}")
    out.append(f"             accepted after the ratchet: {len(cov.debt_accepted)}")
    for item, when in cov.debt_accepted:
        out.append(f"               + {item.item_id}  (added {when})")
    for note in cov.notes:
        out.append(f"           note: {note}")
    return out


def failure_lines(cov: Coverage, tag: str) -> list[str]:
    _, manifest = _WHAT[cov.kind]
    out: list[str] = []
    if cov.gaps:
        out.append(f"[{tag}] COVERAGE_GAP — {len(cov.gaps)} item(s) exist on disk, are "
                   f"NOT declared in {manifest}, and are NOT in the ratchet baseline:")
        for item in cov.gaps:
            out.append(f"   ✗ {item.item_id}")
            out.append(f"       on disk: {item.fingerprint}   ({item.detail})")
    if cov.changed:
        out.append(f"[{tag}] FINGERPRINT_CHANGED — {len(cov.changed)} undeclared item(s) "
                   f"changed value/headline while sitting in the ratchet baseline:")
        for item, was in cov.changed:
            out.append(f"   ✗ {item.item_id}")
            out.append(f"       baseline: {was}   ->   disk: {item.fingerprint}")
            out.append(f"       disk detail: {item.detail}")
    if cov.stale:
        out.append(f"[{tag}] RATCHET_STALE — {len(cov.stale)} baseline entry/entries no "
                   f"longer match the uncovered set (the list must stay exact):")
        for s in cov.stale:
            out.append(f"   ✗ {s}")
        out.append("   Fix: python kb_coverage.py --prune-baseline   (removes only)")
    if out:
        out.append(f"Fix a GAP or a CHANGE by DECLARING the item in {manifest} "
                   f"(preferred), or — if it genuinely cannot be declared yet — widen the "
                   f"debt list on purpose:")
        out.append('   python kb_coverage.py --add-baseline "<id>" --reason "<why>"')
        out.append("   (stamps today's date; reported as post-ratchet debt on every run)")
    return out


def gate(kind: str, tag: str) -> int:
    """Print the ratchet verdict for one manifest. Returns 0 (pass) / 1 (fail)."""
    cov = evaluate(kind)
    for line in summary_lines(cov, tag):
        print(line)
    if cov.failed:
        for line in failure_lines(cov, tag):
            print(line)
        return 1
    print(f"  [{tag}] coverage OK — no new undeclared items, no silent changes, "
          f"{cov.uncovered_n} debt item(s) remaining.")
    return 0


def check(kind: str, tag: str) -> None:
    """Non-destructive summary for refresh.sh / --check (never exits nonzero)."""
    cov = evaluate(kind)
    for line in summary_lines(cov, tag):
        print(line)
    for line in failure_lines(cov, tag):
        print(f"  ⚠️  {line}")


# --------------------------------------------------------------------------- #
# baseline maintenance
# --------------------------------------------------------------------------- #
def prune_baseline() -> int:
    """Remove ONLY stale entries (now declared, or gone from disk). Never adds."""
    baseline = load_baseline()
    removed: list[str] = []
    for kind in ("params", "results"):
        cov = evaluate(kind)
        live = {i.item_id for i in cov.inventory}
        for item_id in list(baseline[kind]):
            if item_id in cov.declared_ids or item_id not in live:
                del baseline[kind][item_id]
                removed.append(f"{kind}: {item_id}")
    if not removed:
        print("nothing to prune — baseline already equals the uncovered set.")
        return 0
    write_baseline(baseline)
    print(f"pruned {len(removed)} stale baseline entry/entries:")
    for r in removed:
        print(f"  - {r}")
    return 0


def add_baseline(item_id: str, reason: str) -> int:
    """Widen the ratchet on purpose. Requires a reason; stamps today's date."""
    if not reason.strip():
        print("refused: --add-baseline requires a non-empty --reason.")
        return 1
    for kind in ("params", "results"):
        cov = evaluate(kind)
        match = next((i for i in cov.inventory if i.item_id == item_id), None)
        if match is None:
            continue
        if item_id in cov.declared_ids:
            print(f"refused: {item_id} is already declared in its manifest.")
            return 1
        baseline = load_baseline()
        entry = {"id": item_id, "reason": reason.strip(), "added": date.today().isoformat()}
        if kind == "params":
            entry["value"] = match.fingerprint
        else:
            entry["headline_sha"] = match.fingerprint
            entry["numbers"] = match.detail
        baseline[kind][item_id] = entry
        write_baseline(baseline)
        print("!! RATCHET WIDENED — undeclared debt accepted on purpose:")
        print(f"     {kind}: {item_id}")
        print(f"     reason: {reason.strip()}")
        print(f"     added:  {entry['added']}  (printed as post-ratchet debt every run)")
        return 0
    print(f"refused: {item_id!r} is not in either disk inventory. "
          f"Run `python kb_coverage.py` to see the exact ids.")
    return 1


def rebuild_baseline(reason: str) -> int:
    """One-shot generator used to write down the state of the world at ratchet
    time. Refuses once a baseline exists — after that the only widening path is
    the explicit, per-item, reason-carrying --add-baseline."""
    if BASELINE.exists():
        print(f"refused: {BASELINE.name} already exists. Use --prune-baseline "
              f"(shrink) or --add-baseline (explicit, per item, with a reason).")
        return 1
    baseline: dict[str, dict[str, dict]] = {"params": {}, "results": {}}
    for kind in ("params", "results"):
        cov = evaluate(kind)
        for item in cov.inventory:
            if item.item_id in cov.declared_ids:
                continue
            entry = {"id": item.item_id, "reason": reason, "added": RATCHET_DATE}
            if kind == "params":
                entry["value"] = item.fingerprint
            else:
                entry["headline_sha"] = item.fingerprint
                entry["numbers"] = item.detail
            baseline[kind][item.item_id] = entry
    write_baseline(baseline)
    print(f"wrote {BASELINE} — params debt {len(baseline['params'])}, "
          f"results debt {len(baseline['results'])}")
    return 0


def main() -> None:
    argv = sys.argv[1:]
    if "--prune-baseline" in argv:
        sys.exit(prune_baseline())
    if "--rebuild-baseline" in argv:
        i = argv.index("--rebuild-baseline")
        reason = argv[i + 1] if len(argv) > i + 1 else f"pre-{RATCHET_DATE}, undeclared"
        sys.exit(rebuild_baseline(reason))
    if "--add-baseline" in argv:
        i = argv.index("--add-baseline")
        if len(argv) <= i + 1:
            print("usage: --add-baseline <id> --reason <why>")
            sys.exit(1)
        item_id = argv[i + 1]
        reason = argv[argv.index("--reason") + 1] if "--reason" in argv else ""
        sys.exit(add_baseline(item_id, reason))
    if "--gate" in argv:
        rc = gate("params", "params") | gate("results", "runs")
        sys.exit(rc)
    rc = 0
    for kind, tag in (("params", "params"), ("results", "runs")):
        check(kind, tag)
        cov = evaluate(kind)
        print(f"  [{tag}] inventory ({len(cov.inventory)} items):")
        for item in cov.inventory:
            state = ("declared" if item.item_id in cov.declared_ids else "UNDECLARED")
            print(f"     {state:10s} {item.item_id} = {item.fingerprint}")
        rc |= 1 if cov.failed else 0
    sys.exit(0)


if __name__ == "__main__":
    main()
