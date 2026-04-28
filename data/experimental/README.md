# Experimental Data — Read Only

These CSV files are the PI's experimental data from the 2026-03-13 spreading assay.

**FORBIDDEN USES**:
- Do NOT use these for parameter fitting / optimization
- Do NOT use these for model selection
- Do NOT load into the simulation as initial conditions

**ALLOWED USES**:
- Visualization overlay only (compare simulation A/A₀ curves alongside)
- Order-of-magnitude sanity check
- Independent qualitative comparison after simulation completes

See `docs/00_project_vision.md` and `docs/12_validation.md` for detailed rationale.

## Data summary

| Condition | N spheroids | Duration | Frame interval | A/A₀ final range | Initial R range |
|---|---|---|---|---|---|
| Bare | 8 | 82 hr | 60 min | 4.6 – 13.4 | 140–419 μm |
| Pre | 25 | 82 hr | 60 min | 4.8 – 26.3 | 103–386 μm |
| Lam4 | 26 | 83 hr | 60 min | 8.0 – 33.1 | 87–255 μm |

## CSV columns
Frame, Time_min, model_choice, used_thresh, review_mode, filename, Area_um2, Area_px, Perimeter_px, Solidity, Circularity, A0_um2, A_over_A0, EffectiveRadius_um, Series

## Mapping to simulation conditions (loose)
- Bare ≈ no protein on pV4D4 (low surface signaling)
- Pre ≈ pre-adsorbed laminin on pV4D4 (moderate signaling)
- Lam4 ≈ laminin in medium (4 μg/mL) (strong signaling)

These map LOOSELY to a φ-spectrum (low → moderate → high), but our 5-point sweep is in **mechanical-parameter space**, not in experimental conditions. The mapping is for visual overlay only.
