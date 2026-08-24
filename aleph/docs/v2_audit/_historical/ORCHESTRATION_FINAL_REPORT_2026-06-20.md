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

# Overnight orchestration — final report (Track 3 coordinator)

**Date:** 2026-06-20
**Plan:** `WARP_SPIKE_THEN_DECOHESION_PLAN_2026-06-19.md` (§Orchestration revision, PI 2026-06-20)
**Coordinator log:** `/Users/sw1/ffn_orchestration/coordinator.log`
**Decision:** Track 2 verdict = **PLATEAU** → **NO graft** (Track 4 not executed). Warp port stands as infrastructure.

---

## ⚠️ CORRECTION / ADDENDUM — 2026-06-21 (supersedes the Track 1 section below)

The **Track 1 (Warp) section below reflects the 13:21 coordinator state and is now STALE.**
The Track-1 worker continued to **15:32** (after this report was written at 13:21) and
resolved both "honest remainders" it had flagged. The corrected record:

1. **REMESH is NOT an open gap — host-hybrid is the correct architecture, and it is built +
   validated.** The GPU side stays a *fixed-size node-pool* (dormant slots, SPLIT activates /
   COLLAPSE parks); the low-cadence remesh runs on the host. This mirrors the reference codes
   (SimuCell3D = CPU/OpenMP, CellSim3D = GPU-but-fixed-topology — **both avoid GPU dynamic
   remesh**). Resync validated single-cell 162→390 nodes (228 splits / 19 events) and multi-cell
   648→672; remesh is **<1 % of multi-cell runtime**. ⇒ "REMESH = untested Phase-C make-or-break"
   is **withdrawn**; no GPU dynamic-remesh kernel is needed.
2. **B3 GPU speed benchmark — RAN on the A5000, results committed** (JSON fixtures
   `dcm_{hybrid,multicell,grid}_bench_a5000.json`). Single-cell **30.3×** per-step vs the existing
   path (net 10–23.5× with remesh cadence); multi-cell hash-grid **10.9× (8 cells) → 36.8×
   (27 cells / 4374 nodes)** — speedup *grows* with cell count (grid O(N·k) vs brute O(N²)).
3. **End-to-end DIFFERENTIABLE loop built + verified** (autodiff vs FD 4e-12, grads w.r.t.
   dP0 / k_edge) — Warp's actual prize, now realized through the full hybrid loop.

Post-milestone commits: `8c22559` (hybrid + B3) · `b2e347f`/`2e9a65c` (multi-cell) · `3e97e68`
(hash-grid + remesh gate) · `075bb6d`/`f6b4c86`/`bef793f`/`5d775df` (differentiable loop, CUDA
graph, active-only) · `ef14380` (buffer-reuse + L-M loop) · `96afdfe` (persistent-grid honest
negative). 24 parity tests green.

**PI decision 2026-06-21:** the engine and physics decisions are **orthogonal**, so engine
adoption proceeds *regardless* of where the de-cohesion physics lands. Warp = GO and **adopted as
the Phase C DCM engine** going forward (differentiability + 30× GPU is exactly what the
fine-grained single-cell line, route a, wants). Phase C IS now triggered. De-cohesion is parked
for a later discussion. Honest residual adoption gates: GPU-*backend* parity (parity above was
measured on Warp's CPU backend; speed on GPU) and porting the substrate forces
(`DcmSubstrateForceGPU` z-well + `DcmSubstrateWettingGPU` in-plane wetting — the mechanistic
spreading driver).

> **⚠️ Track-2 (de-cohesion) PLATEAU verdict is ALSO superseded — see
> `DCM_LANDMINE_REGISTER_2026-06-21.md` (committed on `dcm/decohesion`).** A post-report landmine
> audit found the committed `A/A0 = 1.958/2.286 PLATEAU` was **driven by the `SettlingForce`
> body-force proxy, not the mechanistic lamellipodium it credited** (the `proxy_off:true` flag
> only zeroed `DcmActiveRimTraction`, not `SettlingForce`). The 2×2 controls show the headline
> matches **settle-only** (crawl_residual −1.72), while the only config isolating the
> lamellipodium (ctrl_C: settle-OFF) reaches **A/A0 3.0 — above the claimed ceiling — at genuine
> held-maxZ crawl**, but **diverges past 120k** (`ctrlC_..._RUNAWAY` figs). PI has since **deleted
> `SettlingForce` + added the `forces_manifest` recurrence gate** (commits `da5cff7`/`19b9141`),
> but the decisive clean re-run (settle=0 + `--substrate-wetting` + divergence guard, headline =
> `crawl_residual`) was **not executed**. ⇒ the de-cohesion make-or-break is **UNDECIDED, not
> PLATEAU**. This *reinforces* the engine/physics orthogonality (engine adoption stands either
> way) and is parked for the later de-cohesion discussion.

---

## One-line outcome

The de-cohesion make-or-break came back **PLATEAU (A/A0 ≈ 2.3, not 7–10)** — the platform thesis does **not** clear at the coarse N=100 scale, confirming the Layer-2/center-based structural limit. In parallel, the **Warp port milestone is GO** (15 pieces, fully parity-gated), so the differentiable GPU-resident substrate now exists as infrastructure for the next line — but the one Warp unknown that actually decides a full re-platform (dynamic-topology REMESH) was **not** touched by this spike and remains the true Phase-C make-or-break.

Because A = PLATEAU, **Phase C (full Warp re-platform) is NOT triggered** (it is gated on a validated A/A0→7–10 target, which did not materialize) and **no de-cohesion→Warp graft is performed**.

---

## Track 2 — de-cohesion make-or-break (the platform-thesis decider)

**Verdict: PLATEAU.** `aa0_peak = aa0_final = 2.286` (N=100 cleanball; lamellipodium + junction-switch ON, body-force proxy OFF).

Evidence chain (all guard-rails co-tracked — top-down A/A0 **with** maxZ + COM drift + V/V0, never A/A0 alone):

| stage | A/A0 | maxZ (µm) | V/V0 | COM drift | reading |
|---|---|---|---|---|---|
| A1 smoke N=20 | 1.0→1.448 monotone | 48→33 | ~1.0 | — | wiring functional, real flatten |
| A2 50k interim | 1.0→1.823 monotone | 93→28 (3.3×) | 0.92 | 0.50 µm | still rising |
| A2 250k asymptote | **1.958 flat 110k→250k** | 93→23 (4×) | 0.918 | 1.2 µm | **plateaued** |
| A2 lever band 1.5→3.0 | **2.286** flat from 80k | 93→16 | 0.90 | <1.6 µm | rim 2→11 cells, +17% for 5.5× cells |

- **z-flatten is real** (maxZ drops 4×), volume stable, no bulk translation → the plateau is a genuine spreading ceiling, not a projection or drift artifact.
- **Bottleneck diagnosed:** `detect_rim_cells` is one-shot at build → only ~2 bottom cells ever crawl. Forcing 5.5× more rim cells (contact_band lever) lifted A/A0 only +17% — **diminishing returns = structural ceiling**, not a tuning miss.
- **Interpretation:** consistent with the prior Layer-2 conclusion — the A/A0 **magnitude** belongs to the **fine-grained single-cell line (route a)**, not a center/coarse-CBM upgrade. FORM is reproducible; magnitude is structurally bounded here.
- **Methodology note:** used the PI-validated cleanball bypass (not 1M-agg) to isolate the spreading engine away from the over-compaction confound. Summary JSONs + `results_manifest.yaml` claim committed; commits `9169e84` → `eb4d38d` → `5cf0887`, branch `dcm/decohesion`, worktree `~/ffn_decohesion`.

**Untested escalation flagged FOR PI (not run — needs sign-off):** DYNAMIC rim re-detection (re-detect crawling cells over the run instead of one-shot). This is a **new mechanism** (k_couple / re-detect cadence are magic-numbers) → PI gate-contract/magic-number sign-off required before it can be tried.

---

## Track 1 — Warp parity port (GO/NO-GO on a full re-platform)

**Milestone: GO.** 15 pieces ported, **every one parity-gated vs the committed HOOMD reference on Warp's CPU backend** (guard-rail 2: independent oracle, agent did not grade its own homework). Gate: **21 claims / 15 VERIFIED / 0 drift, 22 pytest pass.** Branch `warp/port`, worktree `~/ffn_warp`, warp-lang 1.14.0.

Pieces verified: B0 install · B1 BAOAB (**bit-for-bit exact**, kT=0 and kT=1.5) · B2 radial-shell ×3 · B4 differentiability (**clean reverse-mode grad, 2.6e-15**) · M-SHAKE · Fixman (→ constrained-integrator stack complete) · compartment bond/angle/LJ (cortex network stack) · DCM node-face contact · DCM turgor · cohesion · lamellipodium. Force-law parities land at **machine-eps (1e-13…1e-16)**; the iterative solvers (M-SHAKE) hit deterministic tol < plan 1e-9.

**GO/NO-GO signals:** parity feasibility = **EASY** (bit/machine-eps); differentiability = **CLEAN** (the actual prize); effort = **~200 LOC/piece, low-friction**. Minor frictions logged (float64 literals, `wp.rint`≠`wp.round`, atomic ordering ~1e-12, PYTHONPATH for the parity script).

**Honest remainders (NOT done — labeled, never silently trusted):**
1. **REMESH (SWAP/SPLIT/COLLAPSE) = the dynamic-topology Phase-C make-or-break.** Every *fixed-mesh* force ported cleanly precisely because Warp's no-dynamic-topology constraint doesn't bite them. The remesh is exactly where it *does* bite — **this spike did not test it**, and it is the real hard part of any full port.
2. `DcmSubstrateForce` z-well — minor per-step force, not yet ported.
3. ~~**B3 GPU speed benchmark — coordinator-gated, UNRUN** (see below).~~ **[SUPERSEDED — see 2026-06-21 addendum at top: B3 RAN on A5000, 30.3× single-cell / 10.9–36.8× multi-cell, JSON committed.]**

Track 1's previous session **died silently** at ~02:48 (no DONE/BLOCKED, working tree clean — nothing lost); a fresh session took over at 12:35, re-verified the 16 committed parity tests, and completed the port. (Coordinator monitoring gap: ~8.5h spent believing a dead worker was active — fixed mid-run by adding status-file staleness detection.)

---

## GPU arbitration summary (single A5000)

- The cleanball sweep (auto-advancing N=200→N=300, **redundant** — the A/A0(R) figure already exists) held the A5000. Per the plan (sweep preemptible; de-cohesion N=100 = top priority) the coordinator **preempted it** to free the GPU for the decisive de-cohesion run. Only loss: fresh pkls; figure + prior pkls remain.
- de-cohesion N=100 owned the A5000 for its 50k + 250k + lever runs; released the lock on completion.
- Track-1 parity ran entirely on Warp's **CPU** backend → never contended.
- ~~**B3 Warp speed benchmark not run:**~~ **[SUPERSEDED 2026-06-21 — RAN; see top addendum.]** It was subsequently executed on the A5000 (30.3× single-cell per-step; 10.9–36.8× multi-cell with hash-grid; JSON fixtures committed).

---

## What this buys / what's next

- **Platform thesis at coarse scale:** does NOT clear (A/A0 bounded ~2.3). The decisive, honestly-measured **bound** is the deliverable — not a borrowed center-based estimate.
- **A/A0 magnitude → fine-grained single-cell line (route a).** This is where the missing magnitude lives; reaffirmed by the structural-ceiling lever test.
- **Warp infrastructure exists** (validated fixed-mesh parity + clean autodiff + GPU-residency path) and is ready to underpin a differentiable layer — **but** a full re-platform decision must first clear the **REMESH dynamic-topology** unknown, which this spike deliberately left as the flagged Phase-C make-or-break.
- **For PI (decisions):** (1) pursue DYNAMIC rim re-detection? (new-mechanism sign-off); (2) commit fine-grained single-cell as the A/A0-magnitude line; (3) want the B3 GPU speed number? (cheap, runnable now).

**No graft performed.** Both worker results are committed and gate-labeled; this report is the integration deliverable for the PLATEAU branch.
