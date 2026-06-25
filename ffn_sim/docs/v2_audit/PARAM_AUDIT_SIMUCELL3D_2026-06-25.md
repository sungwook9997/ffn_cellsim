# Parameter audit vs literature + SimuCell3D — why our cells don't deform (2026-06-25)

**Trigger (PI):** the no-bowl aggregation runs produce ROUND-BALL packing (per-cell asphericity
1.000→0.983 — cells barely deform), not the confluent space-filling faceted tissue the PI wants.
PI: *"changing params to force deformation is itself magic-number tuning — audit the CURRENT params
via TAG first, then find the real literature, including how SimuCell3D sets these."*

## 1. TAG audit of the CURRENT parameters (KB ground truth)

| param | current value | TAG verdict | source |
|---|---|---|---|
| `K_vol` | **7.73e5** | ⚠️ **magic number** — no source row / no KU / no citation; tuned for stability | none |
| `turgor_dP0` | 133 Pa | ⚠️ tuned / band-implied (no sourced row) | none |
| `gamma_surf` (γ) | 1e-4 N/m | ⚠️ 5× below registered MCF7 default (0.5 mN/m); lit range low edge | KB-3.5 |
| `k_edge` | 0.001 N/m | ⚠️ not an independently-sourced row; should derive from γ·geometry | (KB-3.5) |
| cytoplasm η | 65.9 Pa·s | ✅ **source-derived** | Dessard 2024 |

→ 4 of 5 deformation-relevant params are NOT literature-anchored. The "round ball" result rests on
un-grounded constants — it cannot be presented as the physiological prediction.

## 2. Literature + SimuCell3D (deep-research, 104-agent verified, cited)

**SimuCell3D** (Runser, Vetter & Iber 2024, *Nature Computational Science* 4:299–309; bioRxiv
2023.03.28.534574; PMC11052725) — the flagship 3D deformable-cell model — uses LITERATURE-MEASURED
defaults (NOT tuned-to-outcome):

| SimuCell3D param | value | measured range |
|---|---|---|
| cortical surface tension **γ** | **1e-3 N/m** | 5e-4 – 2.5e-3 N/m |
| cytoplasmic **bulk modulus K** | **2500 Pa** | ~2250 Pa |
| max internal net pressure p_max | 2500 Pa | 300 – 2200 Pa |
| adhesion ω = repulsion ξ | 1e9 Pa/m | — |
| bending k_b | 2e-18 J | 1–2e-18 J |
| viscous damping ζ | 2.5e-10 kg/s | — |

- Cytoplasm = slightly-compressible fluid, **p = −K·ln(V/V₀)** — NO explicit Young's modulus, NO
  explicit viscosity.
- **Single-cell cortical tension (lit):** ~0.1–0.7 mN/m rounded/interphase (HeLa interphase 0.17;
  mouse ICM 0.18–0.67), rising ~10× in mitosis (1.3–1.6). Dominantly active (actomyosin); passive
  ≤0.2 mN/m. ⚠️ **No DIRECT MCF7 value exists** — must use the epithelial/cancer range.
- **Animal-cell turgor ≈ 0:** internal hydrostatic EXCESS only ~40 Pa (interphase) → ~400 Pa
  (mitosis), 3–4 orders below the osmolarity scale (Fischer-Friedrich 2014, HeLa).

**Faceting is set by surface tension, not adhesion.** SimuCell3D's dimensionless groups:
`γ̃ = γ/(K·l)`, `ω̃ = ω·l/K`, with `l = ⟨V₀⟩^(1/3)`. In their sweeps **raising γ̃ 0.02→0.10
stratifies/facets the tissue, while a 100× adhesion increase does NOT** — γ is the dominant lever.
(Cross-model: vertex/SPV shape index p₀=P/√A, rigidity transition p₀*≈3.81; deformable-polygon
asphericity confluence at A≈1.16.)

## 3. THE decisive number — our dimensionless surface tension

`l = V₀^(1/3) = ((4/3)π R³)^(1/3) = 12.1 µm` (R=7.5µm).

⚠️ **CORRECTION (verified from the run logs):** the no-bowl runs did NOT use the `K_vol=7.73e5`
*default* — they passed **`--k-vol 1e3` → K = 1000 Pa**, which is already in the SimuCell3D / measured
bulk-modulus ballpark (2250–2500 Pa). So the run's K was fine; the `7.73e5` magic-number default was
overridden. **The single remaining gap is γ.**

| | γ (N/m) | K (Pa) | **γ̃ = γ/(K·l)** | facets? |
|---|---|---|---|---|
| **SimuCell3D (lit)** | 1e-3 | 2500 | **0.033** | ✅ in band (0.02–0.10) |
| **OUR RUN** (`--k-vol 1e3`) | 1e-4 | 1000 | **0.0083** | ❌ ~4× below band |
| **fix: γ→1e-3** (K unchanged) | 1e-3 | 1000 | **0.083** | ✅ in band |
| driver DEFAULT (unused) | 1e-4 | 7.73e5 | 1.1e-5 | ❌ |

So the cells stayed round because **γ is 10× too low** (1e-4 vs the SimuCell3D 1e-3) — that ALONE puts
γ̃ at 0.0083, ~4× under the faceting band. **Raising γ to the literature value 1e-3 N/m → γ̃ ≈ 0.083
(in-band) → faceting, with NO change to K_vol** (the run's 1000 Pa is already grounded). This is the
minimal, lowest-risk, fully lit-anchored fix — and it is NOT outcome-tuning (γ=1e-3 is SimuCell3D's
measured-value default; we are correcting a too-low γ, not inventing one).

The `K_vol=7.73e5` *default* remains a magic number that should be re-anchored to ~2500 Pa for
correctness, but it did not affect these runs.

## 4. Literature-grounded fix (NO outcome-tuning — these are SimuCell3D's measured-value defaults)

| param | current | → grounded value | basis |
|---|---|---|---|
| `K_vol` | 7.73e5 | **→ ~2500 Pa** | cytoplasmic bulk modulus, Fischer-Friedrich 2014 / SimuCell3D |
| `gamma_surf` (γ) | 1e-4 | **→ ~1e-3 N/m** | SimuCell3D default, in 5e-4–2.5e-3 measured band |
| `turgor_dP0` | 133 | → ~40–400 Pa (keep ~100) | animal-cell hydrostatic excess (interphase~40) |
| `k_edge` | 0.001 | → derive from γ·mesh geometry | KB-3.5 |
| cytoplasm η | 65.9 Pa·s | keep | Dessard 2024 ✅ |

→ predicted γ̃ ≈ 0.033 → cells deform into confluent faceted tissue (SimuCell3D regime), **emergently**.

## 5. Implementation CAVEAT (must flag before re-running)

`K_vol=7.73e5` was originally chosen *for stability* (project memory, DCM spreading-foundation).
Dropping it ~300× to the real ~2500 Pa makes the cytoplasm genuinely soft/compressible → the cells
can deform (the goal) BUT the integrator may need the SimuCell3D-style compressible-fluid handling
(`p=−K·ln(V/V₀)`, smaller dt / the implicit solver) to avoid volume collapse. I.e. the stiff K_vol was
masking an integration-stability gap; fixing the physics exposes it. This is engineering, not a
parameter choice — the value 2500 Pa stands.

## Sources (verified)
- Runser, Vetter & Iber 2024, *Nat. Comput. Sci.* 4:299–309 (SimuCell3D) — PMC11052725; bioRxiv 2023.03.28.534574
- Fischer-Friedrich et al. 2014, *Sci. Rep.* 4:6213 (HeLa cortical tension + hydrostatic pressure)
- Boniface et al. 2024, *PRL* 132:248401 (mouse ICM tension, micropipette)
- Bi, Yang, Marchetti & Manning 2016, *PRX* 6:021011 (shape index p₀, rigidity transition)
- Boromand et al. 2018, *PRL* 121:248003 (deformable polygon asphericity)
- KB-3.5 (cortical tension); Dessard 2024 (MCF7 cytoplasm viscosity)
