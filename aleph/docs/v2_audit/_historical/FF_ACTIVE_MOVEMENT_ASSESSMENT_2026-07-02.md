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

# FF active cell movement (migration) — feasibility + gap analysis (2026-07-02)

**Date:** 2026-07-02  **Engine:** FF (Warp)  **Branch:** dcm/main
**PI directive (8h goal):** "ff엔진에서도 능동적 세포 움직임도 다 구현가능한지 보고 … 모든 것은 full compartment
기준으로." Assess whether directed migration/protrusion is implementable in FF, full-compartment.

## Verdict — MAJOR BUILD (not feasible-now, not a modest addition)

FF is a mature **cortical-mechanics / tension** engine with genuine *contractile* + *structural* actives and a
full passive whole-cell compartment stack — but it has **no engine of forward motion**. Net cell translocation
= a phase-scale build (~4–5 mechanistic modules + a new free-translation driver replacing the AFM-only
protocol), NOT a modest extension.

## What EXISTS (the contractile + structural half; file:line)

- Myosin contractile force `network_warp.py:44` (`myosin_kernel`) + engaged-only `:63` + Bell turnover KMC `:86`.
- Arp2/3 branch angle (harmonic 70°, mechanistic) `network_warp.py:135`.
- Hand-KMC kinetic layer `hand_kmc.py:178` (KINESIN/NMIIA/α-actinin/filamin + **INTEGRIN_A5B1** presets, lit).
- Unified builder `weave.py:141` — **lamellipodium** (Arp2/3 dendritic patch), **filopodium** (bundle), **stress
  fiber** (sarcomeric) all present as BUILDERS.
- Full passive compartments: nucleus (law 0), membrane (law 1), cytoplasm (η drag), turgor — wired in the
  whole-cell driver `network_warp.py:380` (but AFM-compression protocol only).

## What is MISSING (the propulsive half) — ranked

1. **Actin polymerization / treadmilling force — MISSING (decisive).** No kernel adds/removes monomers or
   applies a polymerization ratchet at the leading edge (grep `polymeriz|treadmill|elongat|retrograde` → only
   comments). Filaments are built once + only relax. **This is the engine of protrusion; without it there is no
   movement.** Lit anchor: Brownian/elastic ratchet (Mogilner-Oster 1996); force-velocity (Parekh 2005).
2. **Cell polarization / symmetry-breaking — MISSING.** No front-rear asymmetry field. Needed so myosin +
   polymerization sum to NET translation (current contraction is symmetric → zero net motion). Lit: Rac-front/
   RhoA-rear.
3. **Substrate + GPU-resident FA-clutch traction — PARTIAL/unwired.** `fa_anchor.py` (catch-slip clutch,
   INTEGRIN_A5B1 Kong 2009) exists but is a host-only numpy prototype, unwired to any GPU relax/driver, with no
   substrate plane (`relax_on_device` has no fixed-anchor support). Lit: motor-clutch (Chan-Odde 2008).
4. **Retrograde-flow ↔ clutch ↔ protrusion closed loop — MISSING.** The pieces exist separately; the
   Mitchison-Kirschner treadmill + molecular-clutch coupling is not assembled.
5. **LINC nucleus coupling + confining geometry — MISSING (for confined migration).** FF_STAGE6W:118–122
   explicitly omits LINC/perinuclear-actin/MT-cage; no confinement channel (only symmetric plates). Lit:
   LINC/nesprin-SUN (Lombardi 2011).

## Full-compartment coupling status

The whole-cell driver runs cortex+nucleus+membrane+turgor+cytoplasm-drag together (bit-parity validated) — but
only as an **AFM compression** protocol; NO driver permits net cell translation or applies a lateral/directional
force. Membrane could resist a leading edge but nothing pushes it out; nucleus feels cortex only via volume +
contact (LINC omitted).

## Recommendation (my call, per the 8h "remove limits" mandate)

Directed migration IS buildable but is a new phase. **Start with piece 1 (actin polymerization/treadmilling
force)** — it is the decisive gap and the foundation everything else couples to. Then polarization (2), then
wire the FA clutch to a substrate (3), then the closed loop (4). Each mechanistic + lit-anchored, built on the
existing full-compartment whole cell. This is surfaced to PI as a phase-scale decision; if greenlit I lead with
the polymerization ratchet kernel.

Sources: Mogilner-Oster 1996; Parekh 2005; Chan-Odde 2008; Lombardi 2011. Related: FF_STAGE6L (unified
architecture), FF_STAGE6W (LINC omission), FF_STAGE6Z (whole-cell driver).
