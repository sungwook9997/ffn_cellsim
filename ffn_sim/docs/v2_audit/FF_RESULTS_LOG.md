# FF full-compartment cell — quantitative results log

Accumulating record of quantitative results from the FF single-cell engine, with literature reality-bands.
Append a dated entry per result set. Figures: `outputs/ff/figs/`; machine-readable data: `outputs/ff/data/`.
Reality anchors (LITERATURE reference, not fits): MCF7 suspended cortical tension ~10 mN/m (Moazzeni 2021);
interphase cortical tension ~0.3–1 mN/m; resting turgor ΔP ~40 Pa (Fischer-Friedrich 2014, HeLa) up to ~few-hundred
Pa; whole-cell parallel-plate compression force ~10–150 nN (Fischer-Friedrich).

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

### Reality check (model ÷ experiment) — quantitatively NOT matched, two opposite-direction gaps
| quantity | model | literature | ratio | direction |
|---|---|---|---|---|
| resting cortical tension | 0.14 mN/m | ~10 mN/m (MCF7) | **0.014×** | ~70× TOO LOW |
| AFM force @20% strain | ~3900 nN | ~50–150 nN | **~30–80×** | TOO HIGH |
| turgor ΔP @10% strain | ~13 kPa | ~hundreds Pa | **~40×** | TOO HIGH |

- **Too LOW at rest** = the **γ-floor** (myosin active stress ~500× under the measured band from first principles) — long-standing open issue.
- **Too HIGH under compression** = **no poroelastic cytoplasm** (the cell is modeled as a too-incompressible pressurized bag; real cells shed water quasi-statically and are much softer). This is the missing cytoplasm piece.
- So no single scale factor reconciles it — two distinct physics deficits.

**Verdict**: geometry (R, nucleus 0.70R), nucleus-stiffening direction (+45%), and FEM field patterns are
qualitatively right; the quantitative mechanics are off by 1–2 orders in the two ways above. **No direct fit to a
measured MCF7 AFM curve has been done** — the comparison is order-of-magnitude vs the anchors above. Closing the
γ-floor (myosin generation) and adding poroelastic cytoplasm are the prerequisites for a quantitative reality claim.

**Crawl (mt_native_stress, Nc=483,000)**: does NOT translocate (disp = 0 vs real 10–100 nm/s) — traction piece +
GPU-implicit-path protrusion gap. Cortex coherent (V/V0=1.0), 19,545/19,545 FA integrin clutches bound.
