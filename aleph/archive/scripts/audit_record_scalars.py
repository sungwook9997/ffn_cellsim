"""Audit every DERIVED scalar in the committed run records against the raw fields beside it.

The question is not "is the arithmetic right" — it usually is.  It is the one today's session kept
answering the hard way: **does the field cover what its NAME promises?**  Four guards failed that test
in one afternoon (`fraction_of_native` sees only the cortex, `t0` is the first SAMPLED step rather than
the first accepted one, `balance_ok_d` tests adjoint wiring rather than convergence, and a cache check
of mine asked about directory names when the arch tag is a file component).  Each printed a pass.

So this walks the record corpus and, per record, recomputes what it can from the raw fields and reports
three separate things, because they are not the same defect:

  ARITHMETIC   the derived value disagrees with the raw fields it claims to summarise
  SCOPE        the value is right for a narrower quantity than the name states
  DECLARED     the value is an INPUT echoed back, not a measurement of what was built

Read-only.  Writes nothing, runs nothing on a device.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "aleph" / "outputs"

TOL = 1e-9


def close(a: float | None, b: float | None, tol: float = TOL) -> bool:
    """Return whether two optional floats agree to a relative tolerance."""
    if a is None or b is None:
        return False
    if not (math.isfinite(a) and math.isfinite(b)):
        return False
    scale = max(abs(a), abs(b), 1.0)
    return abs(a - b) / scale < tol


def audit(rec: dict, path: Path) -> list[tuple[str, str, str]]:
    """Return (kind, field, detail) findings for one record."""
    f: list[tuple[str, str, str]] = []
    cfg = rec.get("config") or {}
    cen = rec.get("census") or {}
    tim = rec.get("timing") or {}
    mea = rec.get("measurements") or {}
    traj = (mea.get("trajectory") or []) if isinstance(mea, dict) else []

    # ---- timing: pure arithmetic over fields in the same block ---------------------------------
    w, p, n, ni = (tim.get("wall_seconds"), tim.get("physical_time_s"),
                   tim.get("n_steps"), tim.get("n_inner_iterations"))
    if w is not None and p:
        if not close(tim.get("real_time_factor"), w / p):
            f.append(("ARITHMETIC", "timing.real_time_factor", f"{tim.get('real_time_factor')} != {w/p}"))
    if w is not None and n:
        if not close(tim.get("wall_seconds_per_step"), w / n):
            f.append(("ARITHMETIC", "timing.wall_seconds_per_step", f"{tim.get('wall_seconds_per_step')} != {w/n}"))
    if w is not None and ni:
        if not close(tim.get("wall_seconds_per_inner_iteration"), w / ni):
            f.append(("ARITHMETIC", "timing.wall_seconds_per_inner_iteration", "mismatch"))

    # ---- census: the scope question ------------------------------------------------------------
    if "fraction_of_native" in cen:
        fr = cen["fraction_of_native"]
        fil = cen.get("cortex_filaments")
        if fil is not None and not close(fr, fil / 70686.0):
            f.append(("ARITHMETIC", "census.fraction_of_native", f"{fr} != cortex_filaments/70686"))
        # scope: does it move when a non-cortex compartment is reduced?
        if cfg.get("membrane_subdiv") is not None:
            f.append(("SCOPE", "census.fraction_of_native",
                      f"= cortex_filaments/70686 only; this run has membrane_subdiv="
                      f"{cfg['membrane_subdiv']} and n_total_nodes={cen.get('n_total_nodes')}, "
                      f"neither of which can move it"))
    if "cortex_filaments" in cen and cfg.get("filaments") is not None:
        if cen["cortex_filaments"] == cfg["filaments"]:
            f.append(("DECLARED", "census.cortex_filaments",
                      "equals config.filaments exactly — it is the REQUESTED count echoed back, "
                      "not a count of what the builder produced"))

    # ---- t0: is it the first ACCEPTED step, as its note claims? --------------------------------
    t0 = rec.get("t0") or {}
    if t0 and traj:
        first = traj[0]
        step0 = first.get("step")
        ge = cfg.get("gamma_every")
        if step0 not in (0, None) or (ge and ge > 1):
            f.append(("SCOPE", "t0",
                      f"note says 'the first ACCEPTED step' but t0 = trajectory[0], which is step "
                      f"{step0} (t={first.get('t_s')}s) under gamma_every={ge}"))

    # ---- bound_fraction rows -------------------------------------------------------------------
    heads = cen.get("n_heads")
    if heads:
        bad = [r for r in traj
               if r.get("n_bound") is not None and r.get("bound_fraction") is not None
               and not close(r["bound_fraction"], r["n_bound"] / heads)]
        if bad:
            f.append(("ARITHMETIC", "trajectory.bound_fraction",
                      f"{len(bad)}/{len(traj)} rows disagree with n_bound/n_heads"))

    # ---- lifetimes_covered ---------------------------------------------------------------------
    lt, pt, cov = mea.get("bound_head_lifetime_s"), mea.get("physical_time_s"), mea.get("lifetimes_covered")
    if lt and pt is not None and cov is not None and not close(cov, pt / lt):
        f.append(("ARITHMETIC", "measurements.lifetimes_covered", f"{cov} != {pt}/{lt}"))

    # ---- gate ratio ----------------------------------------------------------------------------
    g = rec.get("gate") or {}
    if isinstance(g, dict) and g.get("signal"):
        want = g["residual"] / g["signal"]
        if not close(g.get("residual_over_signal"), want):
            f.append(("ARITHMETIC", "gate.residual_over_signal", "mismatch"))

    # ---- the honest guard: is comparable=False being carried? ----------------------------------
    if tim and tim.get("comparable") is False and tim.get("wall_seconds_per_step") is not None:
        f.append(("NOTE", "timing.comparable",
                  "False — this record's own s/step may not be quoted as engine speed"))
    return f


def main() -> int:
    """Walk the corpus, print findings by kind, and return 1 if any SCOPE/DECLARED finding remains."""
    records = []
    for p in sorted(OUT.rglob("*.json")):
        if "_historical" in p.parts or "obsidian" in str(p) or "tag_kb" in p.parts:
            continue
        try:
            r = json.loads(p.read_text())
        except Exception:
            continue
        if isinstance(r, dict) and (("census" in r and "timing" in r) or r.get("schema", "").startswith("run-record")):
            records.append((p, r))

    print(f"scanned {len(records)} run records under {OUT.relative_to(ROOT)}\n")
    kinds: dict[str, int] = {}
    per_field: dict[tuple[str, str], list[str]] = {}
    for p, r in records:
        for kind, field, detail in audit(r, p):
            kinds[kind] = kinds.get(kind, 0) + 1
            per_field.setdefault((kind, field), []).append(f"{p.relative_to(OUT)} :: {detail}")

    for kind in ("ARITHMETIC", "SCOPE", "DECLARED", "NOTE"):
        rows = {k: v for k, v in per_field.items() if k[0] == kind}
        if not rows:
            continue
        print("=" * 78)
        print(f"{kind}  ({kinds.get(kind, 0)} findings across {len(rows)} field(s))")
        print("=" * 78)
        for (_, field), examples in sorted(rows.items(), key=lambda kv: -len(kv[1])):
            print(f"\n  {field}   [{len(examples)} record(s)]")
            for e in examples[:2]:
                print(f"      {e}")
            if len(examples) > 2:
                print(f"      ... and {len(examples) - 2} more")
        print()

    hard = kinds.get("SCOPE", 0) + kinds.get("DECLARED", 0)
    if not records:
        print("NO RECORDS FOUND — refusing to report a pass over an empty corpus")
        return 2
    print(f"checks executed over {len(records)} records; "
          f"{kinds.get('ARITHMETIC', 0)} arithmetic, {hard} scope/declared finding(s)")
    return 1 if hard else 0


if __name__ == "__main__":
    raise SystemExit(main())
