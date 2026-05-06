# Stage 1e — bounded outcomes (Sim A vs Sim B comparison)

This document records the **PI-pre-authorized bounded-outcome decision tree**
for the Stage 1e Sim A (1D radial ODE) vs Sim B (Stage 1d 3D pilot)
comparison. Written before runs execute.

The companion documents are:
- `docs/v1/stage1e_sanity.md` — six-check Sanity Gate (no new physics
  parameters; analysis layer only)
- `docs/05_radial_approximation.md` — pre-existing radial reduction
  framework
- `docs/06_radial_full_comparison.md` — pre-existing Sim A vs Sim B
  comparison protocol
- `docs/12_validation.md` — Pillar 3 (PI experimental data is cross-check)

## Mechanism question (one)

**How well does the simple 1D radial ODE reduction (Sim A) reproduce
the Stage 1d 3D pilot trajectory (Sim B), as quantified by the
A/A₀ = a + b/R + c/R² fit parameters and the trajectory RMS deviation?**

This is the project's **core methodological deliverable** per
`docs/00_project_vision.md` ("Quantify when/how the radial approximation
is valid"). The answer is a quantitative *validity domain* statement
about the radial-symmetric phenomenological model the PI uses for
experimental analysis.

## Bounded outcomes (4 buckets)

The buckets are defined on the *agreement* between Sim A and Sim B, as
characterised by (a, b, c) parameter relative deviations and A/A₀
trajectory RMS.

| Outcome | Definition | Action |
|---|---|---|
| **Bucket E1 — Sim A reproduces Sim B (radial reduction valid)** | All three of (a, b, c) match within 25% relative deviation AND A/A₀ trajectory RMS deviation ≤ 0.10 AND Pearson correlation ≥ 0.90 | STOP. Surface to PI as **the project's central methodological finding** (Stage 1 publication-strong): the PI's empirical `A/A₀ = a + b/R + c/R²` model is shown to be quantitatively recoverable from a simple 1D radial ODE derived from the same first-principles framework as the 3D simulation; the (a, b, c) parameters are not arbitrary fits but emerge from the layer-by-layer mechanism. PI decides whether to (i) advance to Stage 2 (Layer 6 chemistry/necrosis), (ii) extend to 3-phenotype Sim A vs Sim B (Stage 1e.b), (iii) consolidate the 16-simulation publication-strong narrative + Stage 1e validity finding. **No automatic Stage 2 entry.** |
| **Bucket E2 — partial agreement (validity regime quantified)** | At least one of (a, b, c) matches within 25%, but at least one deviates substantially (≥ 50%); OR A/A₀ trajectory RMS in [0.10, 0.30] | STOP. Surface to PI as **partial validity finding**: the radial reduction captures the qualitative trajectory but not the quantitative parameters; specific (a, b, c) deviations identify which 3D mechanisms (Marangoni, anisotropy, Layer 5) are NOT captured by the reduction. This is itself a publication-relevant finding about the validity boundary of phenomenological radial models. |
| **Bucket E3 — large difference (anisotropy / 3D-only physics important)** | A/A₀ trajectory RMS > 0.30 OR Pearson correlation < 0.50 OR all three (a, b, c) deviate by > 50% | STOP. Surface to PI. Radial reduction does NOT capture the Stage 1d 3D dynamics. The 3D mechanisms (substrate anchor, Marangoni, anisotropic spreading, Path C gravity, Layer 5 spatial heterogeneity) are essential and cannot be subsumed into a simple 1D R(t) ODE. Equally publication-relevant — quantitative *negation* of radial-reduction validity. PI decides among (i) extend Sim A with additional terms, (ii) report negative finding, (iii) advance to Stage 2 with the negative finding noted. |
| **Bucket E4 — Sim A or Sim B fails** | Sim A NaN/divergence/non-convergent fit, OR Sim B trajectory data missing/corrupt, OR scipy fit non-convergent | STOP **immediately**. Implementation issue. Diagnose before proceeding. |

## Per-run gates

Per Stage 1e sanity-md:
- All Stage 1d inherited gates carry over (Sim B already evaluated).
- **Sim A integration finite** (NEW GATE): no NaN/Inf in R(t).
- **Radial fit convergence** (Sim A AND Sim B, both gated separately).
- **(a, b, c), RMS deviation, Pearson** are diagnostics, not gates;
  feed into the Bucket E1/E2/E3 classification.

## Diagnostics (no gate, log-only)

- Per Stage 1e: Sim A trajectory R(t), A/A₀(t), φ̄(t), ζ_eff(t).
- Sim A vs Sim B fit comparison: a, b, c per simulation; relative
  deviations; trajectory RMS; Pearson correlation; residual analysis.
- Plot: A/A₀(t) curves overlaid (Sim A blue, Sim B red); fitted curves.

## Files of record

```
results/stage1e_comparison/{report.md, fits.json}
```

(No new HDF5 — Sim A is analytical, Sim B is the existing Stage 1d
pilot.)

## Auto-progression rules (PI 2026-04-29 full authorization, Stage 1e only)

**Auto-progression IS allowed within Stage 1e**:
- Sanity_md → outcomes → implementation → pytest → analysis script →
  Bucket E1/E2/E3/E4 classification → commit. All within one stage,
  automatic.

**Auto-STOP conditions**:
- Bucket E4 trigger.
- pytest fails before analysis.

**Auto-STOP at end of Stage 1e**:
- After comparison completes (or auto-STOP fires), STOP.
- Stage 2 (Layer 6 chemistry) auto-entry FORBIDDEN.
- Stage 1e.b (3-phenotype comparison) auto-entry FORBIDDEN.
- Inherited fixes (energy-monotone Cousin-Rule, anchor-force K_eff)
  auto-entry FORBIDDEN.

## Stop conditions (general; unchanged)

- Magic number 도입 금지 (Stage 1e introduces no new physical parameters).
- Gate semantic 변경 금지 beyond inherited Stage 1a++ contract changes.
- v13 anti-pattern (measurement-protocol consistency 위반) 시 즉시 STOP.
  Note: Sim A vs Sim B uses RECONCILED A/A₀ measurement protocol per
  sanity-md check 6 (R_eff_xy = √(A_hull/π) for Sim B); v13 risk
  managed.
- 한도 외 escalation 금지 — Stage 2 / 1e.b / inherited fixes 자동
  진입 금지.
