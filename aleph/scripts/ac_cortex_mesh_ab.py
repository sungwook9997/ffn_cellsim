#!/usr/bin/env python
r"""Cortex mesh A/B orchestrator — does the PHYSIOLOGICAL finer cortex mesh change the physics?

The KEY verification the Lead runs the moment gbook (CUDA) returns.  This is a THIN wrapper: it does NO new
physics — it just drives the two EXISTING native drivers with two cortex configs and tabulates the result:

    COARSE (baseline, committed default):  cortex_seg_um=0.5   cortex_density_per_fil=20   →   494,802 actin nodes
    FINE   (physiological ~75 nm mesh):    cortex_seg_um=0.075 cortex_density_per_fil=40   → 2,898,126 actin nodes

For each of the two gates it runs COARSE then FINE and prints a side-by-side comparison:

  * GATE-A (convergence probe, ``ac_gate_a_fq_coarse_test.py``): converged max|PF| over the actin cortex
    (baseline OFF and fiber-quotient ON).  QUESTION: does the resting residual converge AS WELL on the fine
    mesh?  If the fine converged max|PF| is materially higher (or the ON-vs-OFF gap widens/shrinks), the mesh
    resolution CHANGES the resting-baseline physics → PI-relevant mesh sensitivity, not just a solver detail.
  * GATE-B (dynamic driver, ``ac_gate_b_cortex_motor_native.py``): steady-state bound% and emergent γ
    (source / network / total).  QUESTION: does the emergent cortical tension γ change with the mesh?  A finer
    mesh gives many more actin sites within a head's capture reach and a denser crosslink network to carry the
    tension — γ (and bound%) shifting materially between coarse and fine is the mesh-physics signal.

Bit-identical guard: run with the coarse config alone and the two drivers reproduce today's committed behaviour
exactly (the CLI defaults are 0.5 / 20 / 3.0 — the current CellConfig defaults), so this wrapper introduces no
drift; it only ADDS the fine arm.

────────────────────────────────────────────────────────────────────────────────────────────────────────────
RUN ON GBOOK (CUDA A5000; the native runs will NOT execute on the dev Mac — I0-A).  From ~/ffn_ac_native:

  # one-shot both gates, both meshes, with the comparison table:
  PYTHONPATH=~/ffn_ac_native:~/ffn_ac_native/ffn_sim \
    ~/miniconda3/envs/ffn_sim/bin/python \
    aleph/scripts/ac_cortex_mesh_ab.py --gate both --native --steps 40 --dt 0.01 --outer 40

  # or drive the underlying two drivers by hand (what the wrapper runs) — GATE-A convergence probe:
  #   COARSE:
  PYTHONPATH=~/ffn_ac_native:~/ffn_ac_native/ffn_sim ~/miniconda3/envs/ffn_sim/bin/python \
    aleph/scripts/ac_gate_a_fq_coarse_test.py --native --outer 40
  #   FINE:
  PYTHONPATH=~/ffn_ac_native:~/ffn_ac_native/ffn_sim ~/miniconda3/envs/ffn_sim/bin/python \
    aleph/scripts/ac_gate_a_fq_coarse_test.py --native --outer 40 --cortex-seg-um 0.075 --cortex-density 40

  # GATE-B dynamic driver:
  #   COARSE:
  PYTHONPATH=~/ffn_ac_native:~/ffn_ac_native/ffn_sim ~/miniconda3/envs/ffn_sim/bin/python \
    aleph/scripts/ac_gate_b_cortex_motor_native.py --steps 40 --dt 0.01 --outer 40
  #   FINE:
  PYTHONPATH=~/ffn_ac_native:~/ffn_ac_native/ffn_sim ~/miniconda3/envs/ffn_sim/bin/python \
    aleph/scripts/ac_gate_b_cortex_motor_native.py --steps 40 --dt 0.01 --outer 40 \
      --cortex-seg-um 0.075 --cortex-density 40

Long runs: launch under nohup and monitor the LOG FILE (ssh python is not on PATH; use the full env python).

────────────────────────────────────────────────────────────────────────────────────────────────────────────
MEMORY EXPECTATION (fine full stack on the 16 GB A5000):

  The fine cortex is 2,898,126 actin nodes (41 nodes/filament × 70,686) vs 494,802 coarse — a ~5.9× node
  blow-up, plus ~2× the crosslink count (density 20→40 ≈ 1.4M→2.8M crosslinks).  Node state is f64 vec3d
  (24 B/node ≈ 70 MB per whole-cell node array); each driver holds order 10-15 such node/solver work vectors
  → roughly 1-1.5 GB of node arrays, plus the crosslink/segment/bond SoA (small int/float, ~a few hundred MB
  at fine density), plus the fixed membrane(subdiv8)+nucleus+myosin blocks.  Estimated peak is order 3-6 GB
  for the fine full stack — comfortably inside 16 GB.  VERIFY empirically on the first gbook run (watch
  nvidia-smi); this is an estimate, not a measured peak.

  If the full stack turns out tight, reducing the membrane subdivisions for THIS A/B is acceptable: we are
  isolating the CORTEX-mesh effect, so trimming an unrelated compartment's resolution does not confound the
  measurement (pass --membrane-subdiv to GATE-A; GATE-B pins subdiv8 internally — drop it there in the driver
  if needed).  Never lower the cortex biological density to fit memory (CLAUDE.md) — use larger/multi-GPU
  hardware instead.

────────────────────────────────────────────────────────────────────────────────────────────────────────────
This wrapper is CPU-importable (stdlib only; it shells out to the drivers, which own the CUDA build).  On the
dev Mac use --dry-run to print the exact commands without executing them.
"""

from __future__ import annotations

import argparse
import os
import re
import subprocess
import sys
from pathlib import Path

# repo layout: this file is aleph/scripts/ac_cortex_mesh_ab.py
SCRIPTS_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPTS_DIR.parent.parent                        # …/ffn_cellsim
FFN_SIM = REPO_ROOT / "aleph"
GATE_A = SCRIPTS_DIR / "ac_gate_a_fq_coarse_test.py"
GATE_B = SCRIPTS_DIR / "ac_gate_b_cortex_motor_native.py"

# The two cortex configs under test (seg_um, density_per_fil). length stays at the 3.0 µm baseline.
COARSE = {"cortex_seg_um": 0.5, "cortex_density": 20.0}      # committed default → bit-identical baseline
FINE = {"cortex_seg_um": 0.075, "cortex_density": 40.0}      # geometry-verified ~75 nm mesh (2.9M nodes)


def _nodes(seg_um: float, length_um: float, n_fil: int) -> int:
    """Actin node count for a cortex config (nodes/filament = round(L/ℓ₀)+1), for the banner."""
    return (round(length_um / seg_um) + 1) * n_fil


def _gate_a_cmd(python: str, cfg: dict, args) -> list[str]:
    cmd = [python, str(GATE_A), "--outer", str(args.outer)]
    if args.native:
        cmd.append("--native")
    if args.myosin_fraction is not None:
        cmd += ["--myosin-fraction", str(args.myosin_fraction)]
    if args.membrane_subdiv is not None:
        cmd += ["--membrane-subdiv", str(args.membrane_subdiv)]
    cmd += ["--seed", str(args.seed),
            "--cortex-seg-um", str(cfg["cortex_seg_um"]),
            "--cortex-density", str(cfg["cortex_density"])]
    return cmd


def _gate_b_cmd(python: str, cfg: dict, args) -> list[str]:
    cmd = [python, str(GATE_B), "--steps", str(args.steps), "--dt", str(args.dt),
           "--outer", str(args.outer), "--filaments", str(args.filaments), "--seed", str(args.seed)]
    if args.catch_slip:
        cmd.append("--catch-slip")
    cmd += ["--cortex-seg-um", str(cfg["cortex_seg_um"]),
            "--cortex-density", str(cfg["cortex_density"])]
    return cmd


def _run(cmd: list[str]) -> str:
    """Run a driver, streaming its stdout to the console AND capturing it for parsing. Returns full stdout."""
    env = dict(os.environ)
    # make BOTH import styles work: GATE-A imports `ac.*` (needs ffn_sim on path), GATE-B imports `aleph.*`.
    extra = f"{REPO_ROOT}{os.pathsep}{FFN_SIM}"
    env["PYTHONPATH"] = f"{extra}{os.pathsep}{env['PYTHONPATH']}" if env.get("PYTHONPATH") else extra
    print(f"\n$ {' '.join(cmd)}\n", flush=True)
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                            text=True, bufsize=1, env=env)
    lines: list[str] = []
    assert proc.stdout is not None
    for line in proc.stdout:
        sys.stdout.write(line)
        sys.stdout.flush()
        lines.append(line)
    proc.wait()
    if proc.returncode != 0:
        print(f"  [warn] driver exited with code {proc.returncode}", flush=True)
    return "".join(lines)


def _parse_gate_a(out: str) -> dict:
    """Extract the converged actin max|PF| (baseline OFF + fiber-quotient ON) from the GATE-A verdict block."""
    off = re.search(r"baseline \(OFF\):.*?final\s+([\d.eE+-]+)", out)
    on = re.search(r"fiber-quotient\s+start.*?final\s+([\d.eE+-]+)", out)
    return {"off_final": float(off.group(1)) if off else None,
            "on_final": float(on.group(1)) if on else None}


def _parse_gate_b(out: str) -> dict:
    """Extract the steady-state (last-step) bound% and γ (source/network/total) from the GATE-B step table."""
    row = re.compile(
        r"^\s*(\d+)\s+(\d+)\s+([\d.]+)%\s+([\d.eE+-]+)\s+([\d.eE+-]+)\s+([\d.eE+-]+)\s+([\d.eE+-]+)\s+([\d.eE+-]+)\s*$",
        re.MULTILINE)
    last = None
    for m in row.finditer(out):
        last = m
    if last is None:
        return {"bound_pct": None, "g_source": None, "g_network": None, "g_total": None}
    return {"bound_pct": float(last.group(3)), "g_source": float(last.group(4)),
            "g_network": float(last.group(5)), "g_total": float(last.group(6))}


def _fmt(v) -> str:
    return "   n/a  " if v is None else f"{v:>9.4g}"


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Cortex mesh A/B: COARSE (0.5/20) vs FINE (0.075/40) on the GATE-A + GATE-B native drivers.",
        formatter_class=argparse.RawDescriptionHelpFormatter, epilog=__doc__)
    ap.add_argument("--gate", choices=["a", "b", "both"], default="both", help="which gate(s) to run")
    ap.add_argument("--dry-run", action="store_true",
                    help="print the exact coarse+fine commands and exit (CPU-safe; no CUDA needed)")
    ap.add_argument("--python", default=sys.executable, help="python interpreter for the driver subprocesses")
    # GATE-A knobs
    ap.add_argument("--native", action="store_true", help="GATE-A: full 70,686-filament cortex (else 1500)")
    ap.add_argument("--myosin-fraction", type=float, default=None, help="GATE-A: resting bound-myosin fraction")
    ap.add_argument("--membrane-subdiv", type=int, default=None,
                    help="GATE-A: membrane subdivisions (drop for the A/B if the fine stack is tight)")
    # GATE-B knobs
    ap.add_argument("--steps", type=int, default=40, help="GATE-B: accepted physical steps")
    ap.add_argument("--dt", type=float, default=0.01, help="GATE-B: physical timestep dt_phys [s]")
    ap.add_argument("--catch-slip", action="store_true", help="GATE-B: Kovacs catch-slip detach (proxy)")
    ap.add_argument("--filaments", type=int, default=70686, help="GATE-B: cortex F-actin count")
    # shared
    ap.add_argument("--outer", type=int, default=40, help="inner relax / descent iterations (both gates)")
    ap.add_argument("--seed", type=int, default=0, help="RNG seed (both gates)")
    args = ap.parse_args()

    n_fil = args.filaments
    print("=" * 100)
    print("  CORTEX MESH A/B — does the physiological finer cortex mesh CHANGE the physics?")
    print("-" * 100)
    print(f"    COARSE : seg_um={COARSE['cortex_seg_um']:<6} density={COARSE['cortex_density']:<5} "
          f"→ {_nodes(COARSE['cortex_seg_um'], 3.0, n_fil):>9,} actin nodes  (committed default → bit-identical)")
    print(f"    FINE   : seg_um={FINE['cortex_seg_um']:<6} density={FINE['cortex_density']:<5} "
          f"→ {_nodes(FINE['cortex_seg_um'], 3.0, n_fil):>9,} actin nodes  (~75 nm physiological mesh)")
    print("=" * 100, flush=True)

    gates = ["a", "b"] if args.gate == "both" else [args.gate]
    plan: list[tuple[str, str, list[str]]] = []
    for g in gates:
        for tag, cfg in (("COARSE", COARSE), ("FINE", FINE)):
            cmd = _gate_a_cmd(args.python, cfg, args) if g == "a" else _gate_b_cmd(args.python, cfg, args)
            plan.append((g, tag, cmd))

    if args.dry_run:
        print("\n[dry-run] commands that WOULD run (execute these on gbook / drop --dry-run):")
        for g, tag, cmd in plan:
            print(f"\n  GATE-{g.upper()} {tag}:\n    {' '.join(cmd)}")
        print("\n[dry-run] no CUDA touched. Remove --dry-run on gbook to run + tabulate.", flush=True)
        return

    results: dict[tuple[str, str], dict] = {}
    for g, tag, cmd in plan:
        print("\n" + "#" * 100)
        print(f"#  GATE-{g.upper()}  —  {tag}")
        print("#" * 100, flush=True)
        out = _run(cmd)
        results[(g, tag)] = _parse_gate_a(out) if g == "a" else _parse_gate_b(out)

    # ── comparison tables ──────────────────────────────────────────────────────────────────────────────────
    print("\n\n" + "=" * 100)
    print("  MESH-SENSITIVITY COMPARISON  (COARSE vs FINE)")
    print("=" * 100)
    if "a" in gates:
        ca, fa = results[("a", "COARSE")], results[("a", "FINE")]
        print("\n  GATE-A — converged max|PF| over the actin cortex (plateau to beat ≈ 1.16 pN; lower = better):")
        print(f"    {'':<20}{'COARSE':>12}{'FINE':>12}   verdict")
        print(f"    {'baseline (OFF)':<20}{_fmt(ca['off_final']):>12}{_fmt(fa['off_final']):>12}")
        print(f"    {'fiber-quotient ON':<20}{_fmt(ca['on_final']):>12}{_fmt(fa['on_final']):>12}")
        if ca["on_final"] is not None and fa["on_final"] is not None:
            worse = fa["on_final"] > ca["on_final"] + 1e-9
            print(f"    → fine converges {'WORSE (mesh sensitivity)' if worse else 'as well / better'} "
                  f"(Δ = {fa['on_final'] - ca['on_final']:+.4g} pN)")
    if "b" in gates:
        cb, fb = results[("b", "COARSE")], results[("b", "FINE")]
        print("\n  GATE-B — steady-state emergent tension (bound% + γ [pN/µm]):")
        print(f"    {'':<20}{'COARSE':>12}{'FINE':>12}")
        for key, lab in (("bound_pct", "bound %"), ("g_source", "γ_source"),
                         ("g_network", "γ_network"), ("g_total", "γ_total")):
            print(f"    {lab:<20}{_fmt(cb[key]):>12}{_fmt(fb[key]):>12}")
        if cb["g_total"] is not None and fb["g_total"] is not None:
            print(f"    → γ_total shifts {fb['g_total'] - cb['g_total']:+.4g} pN/µm coarse→fine "
                  f"({'MESH-SENSITIVE emergent tension' if abs(fb['g_total'] - cb['g_total']) > 1e-3 else 'mesh-insensitive'})")
    print("\n" + "=" * 100)
    print("  Interpretation: a materially higher fine max|PF| (GATE-A) or a shifted γ_total (GATE-B) means the")
    print("  cortex mesh RESOLUTION changes the physics → surface to PI. Matching values ⇒ mesh-converged.")
    print("=" * 100, flush=True)


if __name__ == "__main__":
    main()
