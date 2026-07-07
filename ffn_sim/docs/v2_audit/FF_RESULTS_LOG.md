# FF full-compartment cell — quantitative results log

Accumulating record of quantitative results from the FF single-cell engine, compared to the **measured MCF7
force–indentation law** (not vague bands): at every indentation you read model-vs-real directly.
Append a dated entry per result set. Figures: `outputs/ff/figs/`; machine-readable data: `outputs/ff/data/`.

Reality anchors (LITERATURE, KB-verified — NOT fits):

- **MCF7 whole-cell modulus, 10 µm colloidal bead** (matched large-contact geometry) **E ≈ 249 Pa** (Zbiral 2023; KB-6.1.1)
- **MCF7 sharp-tip adherent E ≈ 0.2–1.0 kPa** (Li 2008 BBRC 374:609; KB-6.1.1)
- MCF7 cytoplasm G′ ≈ 33 Pa, η ≈ 56 Pa·s (Dessard 2024; KB-3.B3.1) — the missing poroelastic piece
- resting cortical tension: interphase ~0.3–1 mN/m; MCF7 suspended ~10 mN/m (Moazzeni 2021)
- resting turgor ΔP ~40 Pa (Fischer-Friedrich 2014, HeLa) up to ~few-hundred Pa
- whole-cell parallel-plate compression force ~10–150 nN (Fischer-Friedrich)

Real force law: Hertz spherical F(δ) = (4/3)(E/(1−ν²))√R·δ^1.5, ν=0.5, R=7.5 µm; apparent stress σ = F/(πR²); linear σ=E·ε.

---

## 2026-07-07 — Full-compartment AFM indent (native, A5000)

**Config**: `simulate_whole_cell_compression_on_device` — cortex (Nc=494,802) + MT aster (40 tubes) + MTOC +
3000-bead nucleus + lumped membrane γ_mem; rigid parallel plates; strain sweep 0→27%.
**Figure**: `outputs/ff/figs/ff_fullcompartment_results.png` · **data**: `outputs/ff/data/ff_afm_results_2026-07-07.csv`.

| strain | F_plate [nN] | ΔP turgor [Pa] | apparent γ [mN/m] | V/V0 | σ_vm p95 [Pa] | areal strain p95 |
|---|---|---|---|---|---|---|
| 0.00 | 0.01 | 36 | 0.14 | 1.000 | 0 | 0.000 |
| 0.06 | 103 | 4842 | 18.2 | 0.994 | 58 | 0.164 |
| 0.12 | 801 | 18967 | 71.3 | 0.975 | 227 | 1.33 |
| 0.18 | 2806 | 42655 | 160.5 | 0.947 | 530 | 4.08 |
| 0.24 | 6975 | 77637 | 292.5 | 0.908 | 1024 | 9.08 |
| 0.27 | 10185 | 100530 | 379.0 | 0.886 | 1379 | 12.69 |

**Compartment coupling** (|ΔF_plate| at 10% strain when the compartment is removed): **nucleus +45.3%**
(V_cyto=V_hull−V_nuc hydrostatic — genuinely coupled), membrane +0.02% (γ_mem ≪ turgor — negligible),
MT aster +0.04% (L_mt<R, few arms reach the compression axis — negligible).

**Field patterns**: areal strain **localizes sharply at the indentation contact** (p95 → 12.7 at 27%); σ_vm
deviatoric median ≈ 0 (γ-floor), p95 concentrates at the contact rim (→ 1379 Pa).

### Press → stress, model vs measured MCF7 (Hertz E=249 Pa, colloidal, matched geometry)

| strain | depth δ [µm] | F model [nN] | F real [nN] | force ratio | σ model [Pa] | σ real [Pa] | stress ratio |
|---|---|---|---|---|---|---|---|
| 3%  | 0.22 | 13.9  | 0.13 | 108× | 79    | 7.5  | 11× |
| 6%  | 0.45 | 103   | 0.37 | 281× | 583   | 14.9 | 39× |
| 12% | 0.90 | 800   | 1.04 | 773× | 4530  | 29.9 | 152× |
| 18% | 1.35 | 2806  | 1.90 | 1475× | 15876 | 44.8 | 354× |
| 27% | 2.02 | 10185 | 3.49 | 2916× | 57636 | 67.2 | 857× |

### Reality check (model ÷ experiment) — quantitatively NOT matched, two opposite-direction gaps
| quantity | model | measured MCF7 | ratio | direction |
|---|---|---|---|---|
| resting cortical tension | 0.14 mN/m | ~10 mN/m suspended / 0.3–1 interphase | **0.014×** | ~70× TOO LOW |
| force @12% (δ≈0.9µm) | 800 nN | ~1 nN (E=249 Pa) / ~1–5 nN (sharp-tip) | **~770×** | TOO HIGH |
| apparent stress @12% | 4530 Pa | ~30 Pa (σ=E·ε, E=249) | **~150×** | TOO HIGH |
| turgor ΔP @10% strain | ~13 kPa | ~hundreds Pa | **~40×** | TOO HIGH |

The overshoot **grows with indentation** (108×→2916× in force) — the model stiffens super-linearly (incompressible
pressurized bag) where the real cell stays soft (poroelastic water efflux). This is the signature of the missing
cytoplasm piece, and it is why no single scale factor reconciles model and data.

- **Too LOW at rest** = the **γ-floor** (myosin active stress ~500× under the measured band from first principles) — long-standing open issue.
- **Too HIGH under compression** = **no poroelastic cytoplasm** (the cell is modeled as a too-incompressible pressurized bag; real cells shed water quasi-statically and are much softer). This is the missing cytoplasm piece.
- So no single scale factor reconciles it — two distinct physics deficits.

**Verdict**: geometry (R, nucleus 0.70R), nucleus-stiffening direction (+45%), and FEM field patterns are
qualitatively right; the quantitative mechanics are off by 1–2 orders in the two ways above. **No direct fit to a
measured MCF7 AFM curve has been done** — the comparison is order-of-magnitude vs the anchors above. Closing the
γ-floor (myosin generation) and adding poroelastic cytoplasm are the prerequisites for a quantitative reality claim.

**Crawl (mt_native_stress, Nc=483,000)**: does NOT translocate (disp = 0 vs real 10–100 nm/s) — traction piece +
GPU-implicit-path protrusion gap. Cortex coherent (V/V0=1.0), 19,545/19,545 FA integrin clutches bound.
