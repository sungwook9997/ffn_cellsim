# Composite-tension re-validation — DESIGN (KU-3.5 / KU-3.1 re-attribution)

> **STATUS: DESIGN doc — 2026-06-02. Prep, NOT a contract, NO code, NO gate change.**
> Parallel session (`AUTONOMOUS_LOG_2026-06-02.md`). Specifies the measurement protocol
> for the EXTEND memo's flagged "real scientific cost" (`CELL_MECHANICS_EXTEND_VS_REBUILD.md`
> §A.3 Issue 2 / open item): re-derive KU-3.5 (cortical tension) and KU-3.1 (rounding)
> as **composite** gates now that H.8 membrane + H.9 nucleus exist as co-players. This
> is a re-ATTRIBUTION protocol, **not** a band change — the actual gate-contract
> relabel is a PI decision (PI-gated below). No edits to gates/configs/Lead files.

## The problem (EXTEND Issue 2)

KU-3.5 and KU-3.1 currently attribute **all** whole-cell mechanics to the cortex,
because the cortex was the only mechanical element. The shipped H.8 module's own
docstring states it: attaching membrane surface tension makes the measured γ
**composite** — `γ_total = γ_membrane + γ_cortex` (`cell/membrane_surface.py:49-59`;
Sens & Plastino 2015; Fischer-Friedrich 2014). With a membrane and a nucleus sharing
the load, the cortex-only attribution is physically wrong:

- **KU-3.5 (tension):** measured γ = γ_cortex **+ γ_membrane** (composite Laplace balance).
- **KU-3.1 (rounding):** rounding resistance gains a **nuclear** term (H.9).

Until re-derived, current cortex-only results must be **labelled "cortex-attributed",
not "whole-cell"** (EXTEND §A.5 recommendation 3).

## The key insight — re-attribution, NOT a band change

The KU-3.5 / KU-3.1 **bands stay** (they are literature contracts). What changes is the
**decomposition** of the measured quantity into compartment contributions. Because every
EXTEND module is **additive + default-off** (bit-for-bit when off; the seven-subsystem
harness contract), toggling them gives a clean **superposition decomposition** — no new
measurement, no new mechanism.

## The protocol — a factorial toggle of the default-off modules

Run the existing KU-3.5 Laplace-balance / aspiration measurement at each toggle state
(everything else fixed, same seed):

| Run | H.8 membrane | H.9 nucleus | measures | isolates |
|---|---|---|---|---|
| A (current) | off | off | γ_A = γ_cortex | γ_cortex (the current KU-3.5 reading) |
| B | **on** | off | γ_B = γ_cortex + γ_mem | **γ_mem = γ_B − γ_A** |
| C | on | **on** | γ_C, rounding resistance | **nuclear term = C − B** (KU-3.1) |

- **γ_mem** isolated by the A→B difference; check it lands in the membrane band
  (KU-3.B1, T 0.03–0.3 mN/m).
- **γ_cortex** is the run-A reading; check it lands in the cortex band (KU-3.5) — and is
  now free to sit at its *correct* value rather than absorbing the whole-cell tension.
- **γ_total = γ_B** should match the **whole-cell apparent tension** (KU-3.B1 /
  electrodeformation ~10⁻² N/m). The win: cortex-alone over-weights γ_cortex to hit the
  whole-cell band; the composite lets γ_cortex sit lower while γ_total matches.
- **Nuclear rounding term** (C−B) feeds the KU-3.1 composite re-derivation.

## Superposition sanity check (the load-bearing assumption)

The toggle decomposition assumes the contributions **add linearly** at the operating
point (small-deformation Laplace regime). Verify: γ_B − γ_A (membrane increment with
cortex present) ≈ the membrane-only Laplace prediction `2γ_mem/R`, i.e. the membrane and
cortex tensions superpose on the shared shell. If they do NOT (e.g. the shell geometry
shifts enough that the cortex contribution changes between A and B), the decomposition
is coupled, not additive — flag it; the modules act on the same shell tags, so a
geometry shift is possible and must be checked, not assumed.

## Cancer relevance (ties to CANCER_CELLTYPE_PARAM_MAP)

For MCF7 / MDA the membrane + nucleus + cytoplasm are mechanically **co-dominant**
(EXTEND §A.3), so cortex-only attribution is *most* wrong for cancer cells. Composite
re-attribution is therefore a prerequisite for honest cell-type validation — the
whole-cell modulus a cell-type preset must reproduce is the *composite*, not γ_cortex.

## Sanity Gate (recorded now per CLAUDE.md)

1. **Dimensional.** γ [N/m]; ΔP = 2γ/R [Pa]; the decomposition sums tensions [N/m] ✓.
2. **Boundary / default-off.** Run A (both off) MUST equal the current KU-3.5 reading
   bit-for-bit (the additive-module contract) — the anchor of the whole decomposition.
3. **Superposition.** γ_B − γ_A ≈ 2γ_mem/R (membrane increment); else the decomposition
   is coupled — surface it.
4. **Band membership.** γ_mem ∈ KU-3.B1; γ_cortex ∈ KU-3.5; γ_total ∈ whole-cell band.
5. **Sign-sense.** Adding a positive γ_mem raises the measured γ (B > A); adding a stiff
   nucleus raises rounding resistance (C > B). Wrong sign = a module-sign bug.

## What is PI-gated (NOT done here)

- The **gate-contract relabel** itself: marking the current KU-3.5/3.1 results
  "cortex-attributed" and adding a **composite VG** (VG-H3-composite is already seeded
  `blocked` in Notion). That is a gate-contract change ⇒ **PI sign-off** (CLAUDE.md
  no-gate-loosening). This doc only specifies the *protocol*; running it + relabelling
  is the Lead/PI step.

## Open items for PI / Lead

- [ ] Ratify the **A/B/C toggle decomposition** as the composite re-validation protocol.
- [ ] Approve relabelling current KU-3.5/3.1 results "cortex-attributed" (interim) and
      unblocking **VG-H3-composite** once runs A/B/C land.
- [ ] Confirm the superposition check is part of the gate (catch coupled, non-additive
      cases rather than assuming linearity).
- [ ] Sequence after the H.8 composite-tension run + (for KU-3.1) the H.9 nucleus
      integration land — both Lead-owned.
