---
archived_on: 2026-07-28
superseded_by: aleph/docs/v2_audit/AC_EXECUTION_PLAN_2026-07-25.md
reason: >
  Written BEFORE the 2026-07-25 PI reframe, i.e. for a different objective — forward prediction
  and magnitude matching, rather than inferring per-cell-type parameters with the gate on
  mechanical connectedness. Archived, not deleted: its measurements and reasoning stand as a
  record of what was true then. Nothing in it may be quoted as current state; STATE.md is that.
  Selected mechanically: pre-reframe AND cited by no live file (code, STATE.md, CLAUDE.md,
  cell_engine/, gate_contracts/, tests, Makefile). Citations from run outputs and from other
  pre-reframe documents were not treated as protective.
---

# Execution plan — Warp parity spike (B) → de-cohesion make-or-break (A)

**PI decision 2026-06-19:** sequence **B → A**. Do the Warp parity spike first (get
re-platform DATA, not vibes), then the de-cohesion make-or-break.

## Why this is now defensible (the thing that changed today)
A day ago an agentic Warp effort would have produced **silent drift** (unverified code
masquerading as done). Today the **results-integrity system** we built — CI + bit-parity
tests + `run_audit`/`results_manifest.yaml` — makes it **honestly verifiable**: every Warp
piece is *labeled* verified / failing / unrun, never silently trusted. That neutralizes the
safety objection (PI's point, conceded). What the gate does NOT do: make the work complete or
correct. It is a **checker, not a doer** — you get an honestly-labeled *partial* port, not a
finished platform.

## Honest framing of the sequence
- **B is a SPIKE, not a full port.** It ports ~2 representative pieces + benchmarks them. It
  does **not** yield a usable Warp DCM. Therefore **A runs on the EXISTING HOOMD DCM**
  regardless of B — B does not block or enable A; it informs the *downstream* full-port call.
- **A is the platform-thesis decider** (does the shape-resolving DCM reach A/A0→7-10).
- **Full Warp re-platform = Phase C, GATED on B (parity feasible? pays?) + A (validated target /
  tissue-layer need).** Never before both.

## Guard-rails — apply to BOTH phases (the anti-drift contract)
1. Nothing is "done" without a **committed artifact + a passing automated test** (the gate).
2. Warp parity = bit-compare against the **committed HOOMD native reference outputs**, NOT a
   new agent-authored oracle. (Independent grading — the agent must not grade its own homework.)
3. GPU validation must actually **RUN on the gbook CUDA GPU** to render a verdict. Un-run = a
   `needs-regen`/unverified verdict, never `verified`. (Warp is CUDA-only; the Mac cannot validate.)
4. Every headline result enters `results_manifest.yaml`; A/A0 is **top-down silhouette only**.

---

## PHASE B — Warp parity spike (FIRST)  ·  goal: GO/NO-GO data on a full port

### B0 — prereqs (minute 0)
- `python -c "import warp"` on gbook; else `pip install warp-lang` into the `ffn_sim` env.
- **Pin the reference fixture**: commit a small HOOMD-native reference output (positions after
  K steps of overdamped BAOAB at fixed seed + kT, and the radial-shell force on a fixed config)
  — this is the bit-parity TARGET. Without a committed reference, B2/guard-rail-2 cannot grade.
- Decide the single-GPU schedule (see "Tonight" below).

### B1 — overdamped L-M BAOAB step in Warp  (gated)
- Implement the unconstrained overdamped BAOAB update as a Warp kernel (one thread/particle),
  same splitmix64 counter-RNG semantics as `native/.../baoab_kernel.cu`.
- **Gate:** bit-parity test vs the committed reference — kT=0 deterministic bit-for-bit
  (`< 1e-9`), kT>0 bond/position drift `< 1e-7`. Runs on gbook GPU. Done ONLY if it passes.

### B2 — one representative force in Warp  (gated)
- Port `FFNRadialShellForce` (turgor/membrane/nucleus radial shell) as a Warp kernel.
- **Gate:** parity vs the committed compartment-force reference (`1e-12..1e-15`-class), on GPU.

### B3 — micro-benchmark
- Warp vs native steps/s for the B1+B2 system at a fixed N, **same A5000**. Record the ratio.

### B4 — differentiability spot-check
- Take one force law (e.g. B2) and verify a clean reverse-mode gradient w.r.t. a parameter
  (e.g. shell stiffness). This is Warp's actual prize — confirm the autodiff hook is usable.

### B → DECISION (what the spike answers)
| signal | reading |
|---|---|
| parity feasibility (B1/B2) | bit/tol parity reached easily / with friction / infeasible |
| speed (B3) | Warp ≈ / faster / slower than native, same GPU |
| differentiability (B4) | gradient clean / awkward / broken |
| effort | LOC + friction per piece → extrapolated full-port cost |

→ **GO/NO-GO on Phase C (full re-platform).** GO only if parity is tractable AND
(speed-neutral-or-better OR differentiability genuinely pays).

**Risks:** Warp install/version friction; Warp's no-dynamic-topology constraint bites the DCM
remesh/pool later (not visible in the spike — flag it); self-authored-test weakness (mitigated
by guard-rail 2); GPU contention with the cleanball sweep.

**Effort:** ~1–2 days. B1/B2 code + tests are agent-authorable overnight; the GPU parity
verdicts are gated on gbook execution.

---

## PHASE A — de-cohesion make-or-break (SECOND)  ·  goal: does DCM reach A/A0→7-10

Runs on the **existing HOOMD DCM**. The machinery EXISTS and is validated — this is
**integration wiring, not physics invention** (verified by the 2026-06-19 mapping):
junction-switch (GPU-wired) + lamellipodium clutch/tether/ratchet (GPU-ready, single-cell
A/A0→2.03) are built but **unwired into the N=100 spheroid**, which currently uses a coarse
body-force proxy.

### A0 — wiring (CPU code, no GPU)
- `cell/dcm_gpu_build.py`: add `lamellipodium` + `disable_active_traction` flags to
  `build_gpu_dcm_simulation`; wire `detect_rim_cells` + `add_actin_pool_to_snapshot` +
  `ActinClutchAnchorGPU` + `LamellipodialTractionTetherGPU` + `GpuLamellipodiumAdvance`.
- `cell/dcm_active.py`: gate `ActiveRimTraction` (body-force proxy) behind a flag, default OFF.
- `scripts/dcm_two_stage_production.py`: `--lamellipodium` / `--disable-active-traction` CLI.

### A1 — CPU smoke (N=20)
Build + ~5000 steps → no NaN, no pool exhaustion, top-down A/A0 increases, V/V0 ≈ 1.

### A2 — GPU production (N=100, gbook)
1M agg + 50k spread (lamellipodium + junction-switch ON, proxy OFF). Measure **top-down
A/A0(t) + radial centroid drift + V/V0** together. Commit a **lightweight top-down summary
JSON** (not the bulk pkl — gate lesson C12) → add a `results_manifest.yaml` claim
(`declared: verified`) so the gate enforces it.

### A → DECISION (why make-or-break)
- **A/A0 2.5 monotone → 7-10 trend** = DCM validated, platform thesis lives. All other
  workstreams (speed, code-health, γ, gate) were infrastructure for this.
- **plateau 1.5–2** → tune levers (v_front, k_tether, W_cs) → else dynamic `junctional_actin`
  de-cohesion (⚠ needs PI sign-off on k_couple/catch-slip constants = magic-number) → else
  report the **measured** bound (not a center-based-borrowed estimate).
- **<1 compaction** = force-balance bug → debug.

**Guard-rail:** top-down A/A0 ALONE is forbidden (it hides z-flattening) — co-track radial
centroid drift + V/V0.

**Risks:** multi-cell de-cohesion + crawl have ZERO physics validation (parity tests only);
pool exhaustion; pressure-proxy calibration at the coarse scale.

**Effort:** wiring ~2–3 days; smoke + production ~1–2 days GPU.

---

## PHASE C — full Warp re-platform (GATED, NOT before B + A)
Only if B = GO **and** A's result warrants it (a validated DCM to port against, and/or the
need for a Warp differentiable tissue/continuum layer). Structure = SCALE-BRIDGE: keep the
fine-grained HOOMD reference as the parity oracle; Warp becomes the differentiable tissue
layer, ported **piece-by-piece, each gated by bit-parity vs the committed HOOMD reference**.

---

## Tonight — single-GPU scheduling (one A5000)
The cleanball sweep is running (slow; redundant re-derivation; the A/A0(R) figure already
exists). Options, in priority order per the B→A decision:
1. **B tonight:** install Warp on gbook → author B1/B2 + parity tests (agent, CPU) → run the
   parity verdicts on the GPU. This likely **preempts the sweep** (single GPU). The sweep's
   only loss is fresh committed pkls (the figure + old pkls remain).
2. **Sweep finishes, then B:** lower contention, B's GPU verdicts wait for the sweep to end.

Recommended: since B is now the priority and the sweep is largely redundant, **preempt the
sweep for the Warp spike** — but author B1/B2 + tests first (CPU/agent), and only take the GPU
for the parity verdicts + benchmark (short), so little GPU time is spent before the GO/NO-GO.

---

## One-line summary
B (gated Warp spike) buys a **data-backed GO/NO-GO** on re-platforming, safely, because the
gate makes it honest. A (de-cohesion, existing DCM) buys the **platform-thesis verdict**. The
full Warp port is downstream of both — never a blind bet before the cheap decisive tests.

---

## Orchestration revision (PI 2026-06-20) — 2 parallel workers + coordinator
**Decision:** run the **FULL Warp port (B-gated, piece-by-piece)** and the **de-cohesion
make-or-break** as two INDEPENDENT parallel tracks; integrate later. Feasible because
de-cohesion's value is the physics VERDICT + the validated-on-HOOMD wiring, which becomes the
**parity reference** the Warp de-cohesion is later validated against.

- **Track 1 — Warp port (background worker, worktree).** Full port, piece-by-piece, EACH gated
  by bit-parity vs the committed HOOMD reference. **Warp has a CPU backend → parity runs on
  CPU; gbook GPU only for the speed benchmark.** Milestone 1 = B0 (install + Mac-CPU
  feasibility) + B1 (BAOAB) + B2 (radial force) + GO/NO-GO data + effort signal. Comes to rest;
  coordinator dispatches milestone 2 (more pieces).
- **Track 2 — de-cohesion (background worker, worktree).** On the existing HOOMD DCM. A0 wiring
  + A1 CPU smoke now; A2 GPU N=100 is coordinator-gated (single A5000).
- **Track 3 — coordinator (the Lead session).** Wakes on each worker's completion notification;
  sequences the single A5000 (cleanball sweep → de-cohesion N=100 → Warp benchmark); on BOTH
  tracks done, runs the conditional integration.

### Conditional integration (the graft — answers "de-cohesion 나중에 올려도 되지?")
- **IF de-cohesion SUCCEEDED (A/A0→7-10) AND Warp port done → Track 4:** re-wire de-cohesion
  onto the Warp DCM, validated **bit-parity against the HOOMD de-cohesion run**. The graft is a
  **re-wire-on-Warp, NOT copy-paste** (the HOOMD wiring doesn't transfer verbatim; its RESULT
  is the oracle).
- **IF de-cohesion PLATEAUED → NO graft.** Report the measured bound; the Warp port stands as
  infrastructure (differentiability + GPU-residency) for whatever line comes next.

### GPU reality (single A5000)
Track-1 parity runs on CPU (Warp CPU backend), so it does NOT contend. Only the Warp speed
benchmark and de-cohesion N=100 need the GPU — coordinator-sequenced after the cleanball sweep.

### Honest scope caveat
The full Warp port is a multi-milestone program; a background worker advances ONE milestone per
run and comes to rest with honest gate-labeled status (verified / pending / unrun), then the
coordinator re-dispatches. This is "honest incremental progress", not "wake up to a finished
platform" — the gate guarantees the former, never the latter.
