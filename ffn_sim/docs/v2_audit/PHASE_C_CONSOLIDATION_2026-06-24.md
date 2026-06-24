# Phase C — consolidation + PI decision surface (2026-06-24)

Closeout of the autonomous-loop session (branch `h7/m1-ipc-contact`). The bounded autonomous work —
contact, spreading, aggregation, numerics — is characterized; what remains is PI-scoped. This document
consolidates the state so the next direction can be chosen cleanly. PI call this session: **consolidate
+ stop the autonomous loop** (cron `7d93bf4c` cancelled).

---

## 1. Spreading — EXHAUSTED and physiologically grounded (the headline scientific result)

**At the physiological MCF7 cortical tension (γ = 1e-2 N/m), NEITHER a single cell NOR the spheroid
passively spreads.** This is now proven on both scales and is the mechanistic reason the platform's
A/A0 = 7–10 does not emerge from any passive lever.

- **Single cell, physiological γ=1e-2:** A/A0 = **1.019** — stays round (non-wetting). `γ_basal = γ − w_cs
  = +7.15e-3 > 0`, the basal/apical differential is far too weak to overcome cortical tension. Douezan
  `S = w_cs − 2γ = −0.017 < 0` (non-wetting). [commits dee2409, 61b3139]
- **Single cell, γ=1e-4 (100× sub-physiological):** "spreads" — but this is TWO artifacts: (a) γ pushed
  100× below physiological into an artificial S>0 wetting regime, and (b) even there the raw top-down-hull
  A/A0 (~19) was inflated by **12/162 ejected perimeter nodes**; the genuine bulk spread is only **2.14×**.
- **Spheroid, every mechanism:** A/A0 ≈ 1.000 across proxy-free / IPC / deformable / full-stack +
  lamellipodium + polarization + well + penalty + stable-numerics. A collective structural limit.

**Conclusion (re-confirms the PI-ratified Layer-2 call, now at the force-kernel level on both scales):**
the FORM is reproduced; the A/A0 = 7–10 MAGNITUDE belongs to an **ACTIVE** mechanism (lamellipodium /
precursor-monolayer protrusion as the *driver*, Aslemarz/Gupta 2024), not passive wetting/polarization.
Passive differential tension is NOT the spreading mechanism at physiological values.

**→ PI decision (fork a):** build the active spreading mechanism (lamellipodium as DRIVER). Blocked
because its force scale was REFUTED/HALTED to PI in the SF/NMII force-scale work (the missing
motor-density datum). Choosing a force anchor is a PI call — auto-tuning it to make A/A0 rise is the
magic-number sin the PI explicitly forbade.

## 2. M1 contact — FULLY characterized (3 force methods + SimuCell3D projection)

The reconsider directive's lever#3 ("제약방식" / constraint method) is now closed. Full speed/accuracy
frontier (n100, full cad×40/ecm×167 bundle): figure `figs/m1_contact_method_frontier.png`.

| method | final pen | stable @ accel_dt 8e-4 (fast prod) | notes |
|---|---|---|---|
| capped penalty | 3.1 | ✅ | force-based; tunnels under the strong bundle |
| IPC barrier | 2.6 | ✅ | solver-coupled (Hessian in CG) + CCD |
| **IPC penalty-only** | **1.5** | ✅ | **the production recommendation** (good-enough-honest) |
| SimuCell3D projection | **0.0** | ❌ (diverges) | **TRUE penetration-free**, but only at smaller dt |

**Projection (lever#3) full verdict:**
- Achieves **pen → 0** (true geometric non-penetration) — the most accurate contact of any method:
  V/V0 1.000, cfl 0, zero drift. Figure: `figs/n100_projection_pen0_nucleus_montage.png`.
- **A 35-agent adversarial review was load-bearing:** it caught a RUN-BREAKING divergence the headline
  pen→0 hid (the one-sided projection injected momentum → V/V0 = 109 @ step 2527). Fixed with a
  momentum-conserving 50/50 barycentric split (self-test COM drift 0.0066 → 0.0000·me). Reviewer's
  critical #1 (face-normal direction) was NOT applied — judgment over deference (the closest-feature
  r_vec direction is correct for edge/vertex contacts; self-test-proven).
- **Stability ceiling = the timestep, not the gap.** A near-zero gap still diverged at accel_dt 8e-4;
  the cause is the large implicit timestep (Δx/dt blows up). Diagnostic bracket: **stable at 8e-5 (10×
  base, held past step 2400 — through the divergence point of the larger steps), diverges at 2e-4 (25×
  base, V/V0→7.5 by step 1826) and at 8e-4 (100× base, V/V0→1329).** So the stability ceiling sits in
  [8e-5, 2e-4): the projection is penetration-free at **accel_dt 8e-5 ≈ 10× the production wall-time**
  (not the 100× of the base dt, and not as low as 2e-4 which diverges).

**→ PI decision (fork, low-stakes):** accept **IPC penalty-only (pen 1.5)** as the production contact
(de-cohesion runs are A/A0-hull-robust at this pen), with the **projection at accel_dt ~2e-4** as the
penetration-free reference for decisive contact-fidelity checks. No further autonomous contact work.

## 3. Aggregation + numerics — working

Aggregation produces clean round all-touching spheroids (pen ~0.12, V/V0 ~1.0); visualized with the
nucleus-inclusive z-montage per the standing PI viz directive. Numerics: stable configs identified per
contact method + the projection-timestep incompatibility characterized.

## 4. Process notes (durable)

- **Sync:** Syncthing does NOT sync source code (only `outputs/`). New code must be rsync'd to gbook
  (`/home/sungwook/ffn_phase_c`) before launching — recorded in memory. gbook runs two syncthing
  instances (a half-finished upgrade; PI cleanup, irrelevant to code-sync).
- **Colab:** notebook prepared (`scripts/colab/ffn_phase_c_colab.ipynb`) but unused — there is **no
  autonomous Colab-execution tool** surfaced to Claude (only a notebook-EDIT tool). gbook A5000 was the
  only GPU driven. Colab parallelism needs (1) the branch pushed to GitHub (PI sign-off) and (2) the PI
  to drive the Colab UI or expose an execution interface.
- **Adversarial review pattern** (find → verify → judgment-triage) caught a divergence an aggregate
  metric hid — worth repeating before trusting any new-method headline.

## 5. PI decision surface (what the next direction depends on)

| fork | what it unblocks | blocker |
|---|---|---|
| **(a) active spreading** = lamellipodium as DRIVER | the A/A0 = 7–10 magnitude (the platform's goal) | force-scale anchor (REFUTED/HALTED to PI; choosing it is a PI call, not auto-tunable) |
| **(b) C6 contact-line substrate model** | a bounded wetting equilibrium (Young angle) | major build; payoff uncertain (spheroid is a collective limit, robust to substrate) |
| **(c) branch push** for Colab parallel | second GPU lane | `h7/m1-ipc-contact` (~100 commits) push needs PI sign-off (ffn/foundation rule) + an autonomous Colab exec path |

**Recommendation:** the platform's goal lives in **(a)** — but it is genuinely PI-gated on the force
anchor. (b) and (c) are secondary. The honest state is that further gbook cycling would be re-litigation
until the PI sets the (a) direction.
