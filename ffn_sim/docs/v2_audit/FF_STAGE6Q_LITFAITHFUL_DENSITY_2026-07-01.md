# FF Stage 6Q — Lit-faithful cortical densities + the clean γ_myo floor metric

**Date:** 2026-07-01  **Engine:** FF (Warp, GPU-native)  **Hardware:** gbook RTX A5000
**Branch:** dcm/main  **Supersedes (magnitude):** FF_STAGE6O native γ, FF_STAGE6P framing

This stage closes the user's "a, b, c — 다" tasks: bring every cortical density to its **direct measured
lit value at the MCF7 radius**, report the floor against the **MCF7-faithful active band**, and — the
substantive finding — replace the passive-contaminated `γ_active` sum with the **clean `γ_myo` channel**
as the floor metric. The γ-floor PERSISTS and DEEPENS to ~1700× vs the MCF7-active band; the deepening
is the honest consequence of the measured-sparse Nie myosin density plus the metric correction.

---

## (a) Per-motor stall — lit range, swept controlled variable

The NMIIA minifilament per-side stall is a **literature RANGE**, documented (not tuned to a number):

| anchor | heads/side | per-side stall | basis |
|---|---|---|---|
| Stam-Hocky (AFINES) | 10 | **5 pN** (10 × 0.5 pN, Kovacs) | coarse minifilament (conservative) |
| Stachowiak 2012 | ~17 | ~8–9 pN ensemble | ensemble force |
| Billington 2013 | ~29 (structural) | ~15 pN | cryo-EM NMIIA minifilament headcount |

`NMIIA_MINIFIL_STALL_PN = 5.0` is kept as the **conservative anchor**; f_myo is a **swept controlled
variable** (`gamma_floor` f_myo sweep) — NOT a magic number. The A5000 sweep confirmed the floor is robust
across the whole range: 216× (5 pN) → 72× (15 pN) → 64× (17 pN) → 36× (30 pN) vs the active band — the
floor PERSISTS even at the structural per-motor force. (These factors were computed pre-(c) on `γ_active`;
the post-(c) clean `γ_myo` factors are ~3–5× deeper, below.)

## (b) Band-reporting — MCF7-faithful active target

`gamma_floor_production` now reports the floor against the **active fraction (~70%) of the MCF7 suspended
IQR (Hosseini 2020, 180–400 → active 126–280 pN/µm)** as the PRIMARY target (`floor_vs_mcf7_active`),
plus `floor_vs_salbreux_active` (Chugh/Salbreux 70%) and the legacy `floor_factor_under_band` (generic
total) for continuity. FF γ is the ACTIVE channel → it must be compared to the active fraction, not the
total measured tension (which is ~70% myosin + passive floor + turgor partner; FF_STAGE6P).

## (c) Lit-faithful densities — actin (KB-3.18) and myosin (Nie), decoupled

Every cortical density is now its **direct measured areal value at MCF7 R = 7.5 µm** (area 4πR² ≈ 707 µm²),
per the physiological-baseline HARD rule — not a convenient ratio:

| species | areal density | source | N at R=7.5 |
|---|---|---|---|
| actin (F-actin) | **~100 /µm²** | KB-3.18 cortex composition | **70 686** |
| NMIIA minifilament | **0.625 /µm²** | Nie 2015 (only direct cortical datum) | **442** |
| α-actinin crosslinker | 1:1 with actin (γ-irrelevant, FF_STAGE6M) | — | 70 686 |

This **supersedes** the prior FF actin count (21375 ≈ 30/µm², itself the CLAUDE.md "~38000 native"
≈ 30/µm² scaled to R=7.5) — which was **~3× sparser** than the project's own KB-3.18 datum.

> ⚠️ **CLAUDE.md ↔ KB conflict surfaced to PI.** CLAUDE.md states "~38,000 native filaments per cell"
> (= 30/µm² at R=10); KB-3.18 states the cortical actin density is "~100/µm²" (its companion "N~1.3e4"
> is internally inconsistent with its own density — implies area ~130 µm², likely a stale/typo'd count).
> These disagree ~3×. The actin filament COUNT does **not** drive the γ-floor (γ is myosin-bound), so
> this is resolved here as a **citation/faithfulness fix** (use the KB-3.18 density), with the conflict
> flagged for a PI CLAUDE.md/KB reconciliation — not a floor lever.

**Decoupling note:** raising actin to 70686 via the old `n_myo = n_fil/10` ratio would have put myosin at
~10/µm² = 16× Nie. The myosin density is therefore DECOUPLED to the measured Nie 0.625/µm² (442) — the
physiological-baseline-correct operating point. Crosslinkers stay 1:1 with actin (γ-irrelevant).

---

## THE SUBSTANTIVE FINDING — γ_active is passive-contaminated; γ_myo is the clean metric

At the native lit-faithful density the `γ_active` SUM (Σ over actin + crosslink + myosin plane-crossing
tension) is **NOT a usable myosin metric** — it is dominated by a **structural passive actin-network
residual**:

```
gamma_active(f_myo=0) vs relaxation steps (A5000, N=70686):
  steps=  300 : ga(f0)=0.8388   ga(f5)=0.6513   gmyo(f5)=0.0993   gactin(f0)=0.8387
  steps= 1000 : ga(f0)=0.8575   ga(f5)=0.6701   gmyo(f5)=0.0993   gactin(f0)=0.8575
  steps= 3000 : ga(f0)=0.8579   ga(f5)=0.6705   gmyo(f5)=0.0993   gactin(f0)=0.8579
  steps= 8000 : ga(f0)=0.8566   ga(f5)=0.6692   gmyo(f5)=0.0993   gactin(f0)=0.8567
```

1. **`γ_active(f_myo=0) ≈ 0.84 pN/µm ≠ 0`, and it PLATEAUS** (300→8000 steps) — so it is **structural, not
   under-relaxation**. The carrier is `γ_actin` (the inextensibility constraint forces), NOT crosslinks
   (`γ_xl ≈ 0.001`). It is the frustration of a densely-crosslinked geodesic-arc network on the curved
   shell. It scales with actin density (invisible at the small-N / ×40 scale, dominant at native).
2. **Turning myosin ON DECREASES `γ_active`** (0.858 → 0.518 at f=10): myosin contraction RELAXES the
   passive residual. So `γ_active(f5) − γ_active(f0) < 0` — `γ_active` cannot measure myosin tension here.
3. **`γ_myo` is artifact-free**: linear in f_myo with zero intercept — `γ_myo = 0.0199 · f_myo`
   (0.0, 0.0199, 0.0496, 0.0993, 0.1986 at f = 0, 1, 2.5, 5, 10). This is the honest myosin-induced
   cortical tension. **The floor is now reported on `γ_myo`.**

**ADVERSARIAL CHECK — physical, not a measure_gamma bug.** A measurement bug would give a constant
offset independent of N; a physical density-frustration must vanish as density → 0. The residual scales
LINEARLY with filament areal density (A5000, settle 600 steps), carried by `γ_actin` not `γ_xl`:

```
 N_fil   areal   γ_active(f0)  γ_actin(f0)  γ_xl(f0)   γ_active/density
   300     0.4      0.0033       0.0035      0.0006      0.00781
  1000     1.4      0.0115       0.0115      0.0006      0.00810
  5000     7.1      0.0595       0.0594      0.0007      0.00841
 20000    28.3      0.2431       0.2431      0.0011      0.00859
 70686   100.0      0.8542       0.8541      0.0014      0.00854
```

`γ_active(f0)/density ≈ 0.0085 pN/µm per (fil/µm²)` is constant over a 235× density range → a genuine
density-linear network frustration (inextensible crosslinked filaments on the curved shell cannot relax
all force-free-formed crosslinks simultaneously), NOT a constant artifact. This is also WHY the ×40 /
small-N γ-floor work never saw it (at N=1000 it is 0.0115 pN/µm ≪ the myosin signal). At native density
it dominates `γ_active`, mandating the `γ_myo` metric.

**SOURCE — crosslink-network frustration, not arc geometry.** Toggling the crosslinks isolates it (A5000):

```
 N_fil   crosslinks   γ_actin(f0)
 20000      ON          0.2431
 20000      OFF         0.0000     ← bare great-circle arcs (geodesics) carry ZERO residual
 70686      ON          0.8541
 70686      OFF         0.0000
```

With crosslinks OFF the residual is EXACTLY zero — the bare filament arcs are geodesics (bending-free at
rest), so it is NOT a discretization / arc-construction artifact. The entire residual is the **crosslink
prestress**: force-free-at-formation crosslinks bridging different fibers become geometrically frustrated
against the inextensible segments as the network settles → a self-stressed (tensegrity-like) equilibrium.
This is a physically-reasonable passive cortex prestress (8.5e-4 mN/m at native density — well below the
~0.04 mN/m measured blebbistatin-insensitive passive floor, and far below the band), and it is correctly
EXCLUDED from the active tension by the `γ_myo` dipole channel (which is unaffected by the crosslink state).

This is a **metric correction, not a gate-loosening**: switching to `γ_myo` makes the floor DEEPER (more
conservative). `γ_myo` was already the established clean active channel (FF_STAGE6M, bit-identical under
the link_k correction). `γ_active` is retained in the output, **flagged** as the passive-contaminated sum.

---

## Production γ-floor at lit-faithful densities (A5000, 3 realizations, 27.7 s)

```json
{
  "gamma_myo_mean":   7.42e-05  (mN/m)   // CLEAN myosin channel — the floor metric
  "gamma_myo_std":    4.8e-07
  "gamma_active_mean":1.39e-04  (mN/m)   // passive-contaminated SUM (flagged, not the floor)
  "floor_vs_mcf7_active":     1699×       // PRIMARY — vs MCF7 Hosseini active (126 pN/µm)
  "floor_vs_salbreux_active": 3303×       // vs Chugh/Salbreux active fraction
  "floor_factor_under_band":  4719×       // vs generic total band (legacy)
  "n_fil": 70686, "n_xl": 70686, "n_myo": 442, "f_myo": 5.0
}
```

**Headline:** at the FULLY lit-faithful MCF7-R cortex — actin 100/µm² (KB-3.18), NMIIA 0.625/µm² (Nie 2015),
per-motor 5 pN (Stam-Hocky) — the clean active cortical tension is **γ_myo ≈ 7.4e-5 mN/m, floored ~1700×
under the MCF7-active band**. DEEPER than the prior reported ~416× (FF_STAGE6O/P), because (1) the
Nie-faithful myosin density is 4.8× sparser than the prior ratio, and (2) the clean γ_myo metric removes
the passive-residual inflation that had flattered γ_active.

**This CONFIRMS the γ-floor, it does not rescue it.** Using the measured densities lowers γ further (Nie is
sparse) — exactly the [[project-gamma-floor-likely-deficit]] conclusion. The one open lever is unchanged:
the missing experimental **MCF7-adherent load-engaged NMII density** datum (engaged density > Nie presence),
the SF/NMII twin. EXPERIMENTAL gap, not a model mechanism.

**Lit re-check (2026-07-01, focused PubMed):** "myosin minifilament density" → 4 hits, NONE a cortical/SF
areal density above Nie 2015 (van Loon 2021 *Mol Biol Cell* doi:10.1091/mbc.E21-05-0258 images individual
apical-cortex NMII minifilaments LIVE but reports no density — the closest methodology if a wet-lab
measurement is commissioned; the rest are qualitative / structural / RBC). The documented absence STANDS
(spot re-check, not exhaustive; Consensus quota was exhausted).

---

## Sanity gate

- **Dimensional:** γ_myo [pN/µm] = f_myo [pN] × (Σ plane-crossing geometric factor) / circumference [µm]. ✓
- **Boundary:** γ_myo(f_myo=0) = 0 exactly (no motor → no active tension). ✓
- **Monotonic / linear:** γ_myo = 0.0199·f_myo, zero-intercept, R²≈1 over f∈[0,10]. ✓ (controlled-variable, not tuned)
- **N-faithful:** densities = direct measured areal values (actin KB-3.18, myosin Nie) at MCF7 R; no count chosen to hit a band. ✓
- **Metric hygiene:** γ_active passive residual is structural (plateaus), carried by γ_actin; γ_myo is the artifact-free channel. ✓
- **No magic number:** the floor DEEPENED under the correction (conservative); no parameter moved to make a gate pass. ✓

## Files
- `ff/cortex_assembly.py` — CortexParams.n_filaments 21375 → 70686 (KB-3.18 100/µm²).
- `ff/gamma_floor.py` — PROD_N_FIL/XL/MYO = 70686/70686/442 (Nie myosin); production floor on clean γ_myo.
- `ff/architecture_spec.py` — CORTEX n_filaments 70686, motor density_per_fil 0.00625 (Nie/actin).
- `tests/ff/test_cortex_assembly.py`, `test_weave.py`, `test_gamma_floor.py` — updated for new densities + γ_myo metric.

Related: [[project-gamma-floor-likely-deficit]], [[project-ff-gpu-native]], [[feedback-no-param-tuning-to-outcome]],
[[feedback-physiological-baseline]], [[project-sf-nmii-forcescale-result]].
