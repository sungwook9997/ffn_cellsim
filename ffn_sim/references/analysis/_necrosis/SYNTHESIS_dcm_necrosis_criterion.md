---
id: SYNTHESIS_dcm_necrosis_criterion
topic: DCM 3-zone (proliferating / quiescent / necrotic) multi-factor necrosis & quiescence criterion
inputs: Q1_diffusion_critical_size.md, Q2_mcf7_livedead.md, Q3_solid_stress.md
purpose: a calibrated, mechanism-level per-cell state updater combining NUTRIENT/O2 access (rim-depth vs penetration depth) AND MECHANICAL pressure (local compressive stress) for the DCM spheroid layer
date: 2026-06-11
status: synthesis / criterion contract — written BEFORE the run; calibration targets are OVERLAY-ONLY (never fit), per project hard rules
---

# SYNTHESIS — DCM 3-zone necrosis / quiescence criterion (MCF-7 calibrated)

This file fuses the three research notes into a single per-cell state-update rule for the
DCM (multi-cell / spheroid composition) layer. A cell's fate is decided by **two
co-located drivers** that both peak adverse in the core:

- **(a) Nutrient/O2 access** — the cell's depth below the cluster surface vs the O2/nutrient
  diffusion **penetration depth** (the Greenspan/Thomlinson-Gray viable rim).
- **(b) Mechanical compression** — the cell's local **compressive solid stress** (proxied by
  neighbor/contact count or measured contact-normal force) vs a quiescence/death threshold.

Both are evaluated; the cell takes the **worst (most-arrested/most-dead) verdict** of the two
channels (a strict OR for death, OR for quiescence). This reproduces the canonical three-shell
MCTS architecture — proliferating rim → quiescent shell → necrotic core — with a viable rim
that stays ~constant in thickness while the necrotic core grows with diameter.

---

## 1. The calibrated band numbers (decision table)

These are the numbers the updater uses. Every value is anchored to a source in §5; the
single chosen value is in **bold**, the literature range in parentheses.

### 1a. Nutrient/O2 (rim-depth) channel — depth `d = R_cluster − r_cell` from the surface

| Band | Symbol | Value | Range (lit.) | Meaning |
|---|---|---|---|---|
| **Penetration depth** (O2/nutrient diffusion length, consumption-limited) | `L_pen` | **150 µm** | 100–200 µm (≈150–160 µm) | distance O2 keeps tissue *alive* from the source |
| **Proliferative shell depth** (outer, division-competent) | `d_prolif` | **40 µm** | 20–50 µm (~2–3 cell layers = 20 µm minimum) | `d < d_prolif` → division allowed |
| **Viable-rim thickness** (proliferating + quiescent, total live shell) | `L_viable` | **150 µm** | 100–200 µm (constant as core grows) | live cells live within `L_viable` of the surface |
| **Quiescent shell** | — | `d_prolif ≤ d < L_viable` = **40–150 µm** | — | viable but division-arrested |
| **Necrotic depth** (beyond diffusion-supported life) | — | `d ≥ L_viable` = **≥150 µm** | — | nutrient channel → necrotic |

Equivalent radial form (the "constant-rim, growing-core" invariant, used directly in code):
- Necrotic-core radius `R_nec = max(0, R_cluster − L_viable)` → cells with `r_cell < R_nec`
  fail the nutrient channel.
- Quiescent-onset radius `R_q = max(0, R_cluster − d_prolif)` → cells with `R_nec ≤ r_cell < R_q`
  are nutrient-quiescent.

### 1b. Mechanical (compressive-stress) channel — local stress σ per cell

| Band | Symbol | Value | Range (lit.) | Meaning |
|---|---|---|---|---|
| **Quiescence onset** (proliferation slows) | `σ_q_onset` | **1 kPa** | 0.5–1 kPa | smooth p27-like ramp begins |
| **Full proliferation arrest** (quiescence) | `σ_q_full` | **5 kPa** | 5 kPa (Delarue, Montel) | division rate → 0 by here |
| **Necrosis / mechanical death** (sustained) | `σ_death` | **10 kPa** | 6–16 kPa full inhibition; 10 kPa "growth stops" | sustained σ ≥ this over residence time → deactivate |
| **Volume-checkpoint proxy** (preferred read-out) | `ΔV/V₀` | **−5 % at 5 kPa** | ~5 % loss in ~30 min @ 5 kPa | arrest when volume loss exceeds checkpoint |
| Real-tumor residual-stress band (sanity) | — | 3.3–10.9 kPa | — | in-vitro thresholds are physiological |

Pressure proxy when σ is not directly measured (cheap DCM first pass — see §3):
contact/neighbor count `n_c`. For ~15 µm cells the bands map onto coordination number:
- `n_c ≲ 6` (surface/loose) → low stress (rim regime, σ < `σ_q_onset`).
- `n_c ≈ 7–11` (interior, partial confinement) → quiescence regime (σ between onset and full).
- `n_c ≳ 12` (fully caged, dense core packing) → high stress (σ ≥ `σ_q_full`, approaching death).
Prefer the **measured contact-normal force / HOOMD virial stress** over `n_c` whenever the
contact model is wired (no-abstraction rule); `n_c` is the fallback proxy only.

### 1c. Spheroid-scale geometry the bands must reproduce (MCF-7 calibration)

| Quantity | Value | Source |
|---|---|---|
| **Critical diameter — no dead core below** | **≈ 300–400 µm** (MCF-7 100 & 300 µm all-viable) | Q2 |
| **Critical diameter — necrotic-core onset** | **≈ 400–500 µm** (use **500 µm**) | Q1, Q2 |
| **Prominent necrotic core** | **≈ 700 µm / MCF-7 day 8** (none at day 6) | Q2 (Lee 2010) |
| Hypoxic-core onset (pre-necrotic) | ≈ 200 µm diameter | Q1 |
| Growth plateau (diffusion caps shell) | 600–800 µm | Q1 |
| Cell diameter (MCF-7) | ≈ 15 µm | Q1, Q2 |

Note the geometric self-consistency check: with `L_viable = 150 µm`, the necrotic core first
appears when `R_cluster > L_viable`, i.e. diameter > 300 µm, and is robust by ~400–500 µm —
matching the MCF-7 "all-viable at 300 µm, core onset ~500 µm" landmarks. (If a later overlay
shows MCF-7 onset closer to 500 µm than 300 µm, raise `L_viable` toward 180–200 µm within the
cited range — do NOT tune outside 100–200 µm and do NOT tune to fit a PI CSV.)

---

## 2. The multi-factor necrosis / quiescence criterion (the rule)

For each cell, compute two channel verdicts and take the worst.

Let:
- `d = R_cluster − r_cell` = rim-depth (distance from the cluster surface).
- `σ` = local compressive solid stress (measured contact-normal stress / virial; or the
  `n_c` proxy mapped to a stress band per §1b).
- `t_σ` = cumulative residence time the cell has spent at `σ ≥ σ_death`.

**Nutrient channel verdict `S_nut`:**
- `d < d_prolif` (40 µm) → `PROLIFERATING`
- `d_prolif ≤ d < L_viable` (40–150 µm) → `QUIESCENT`
- `d ≥ L_viable` (≥150 µm) → `NECROTIC`

**Mechanical channel verdict `S_mech`:**
- `σ < σ_q_onset` (1 kPa) → `PROLIFERATING`
- `σ_q_onset ≤ σ < σ_death` (1–10 kPa) → `QUIESCENT` (smooth p27-like ramp; division rate
  multiplier `g(σ)` falls from 1 at 1 kPa to 0 at `σ_q_full` = 5 kPa, then 0 up to death)
- `σ ≥ σ_death` (10 kPa) **for** `t_σ ≥ τ_death` → `NECROTIC` (else `QUIESCENT`, reversible)

**Combined verdict (worst-of-two, with PROLIFERATING < QUIESCENT < NECROTIC):**
```
S_cell = max_severity( S_nut , S_mech )
```
- A cell **PROLIFERATES** only if BOTH channels say PROLIFERATING (rim AND low-stress).
- A cell is **QUIESCENT** if either channel is quiescent and neither is necrotic — REVERSIBLE:
  if it later returns to the rim (cluster shrinks/reorganizes) or decompresses, it resumes.
- A cell **NECROSES** (deactivate / mark dead) if EITHER channel is necrotic:
  - nutrient: `d ≥ L_viable` (the diffusion-limit / media-contact-loss rule), OR
  - mechanical: sustained `σ ≥ σ_death` for `t_σ ≥ τ_death`.
  Necrosis is **irreversible** (once dead, the cell stays deactivated; it may still occupy
  volume / lose contact-bearing, matching core densification + later disintegration).

This is exactly the "compound" criterion the research notes converge on: the diffusion-limit
depth rule (Q1/Q2) OR-ed with the mechanical compression rule (Q3), both keyed to the core.

---

## 3. The concrete DCM updater (per-cell, per-decision-step)

Pseudo-code for the state updater, called every `Δt_DCM` (cell-fate timestep, ≫ the HOOMD
mechanical step). Phase-1 pass = **route A (geometric/contact proxy, no PDE)**; the mechanism
hook for route B (explicit O2 field) is noted where it slots in.

```python
# --- per-cluster setup (once per DCM step) ---------------------------------
com      = center_of_mass(cells)                 # cluster centroid
R_cluster = surface_radius(cells, com)           # outer radius (e.g. 95th-pct |r-com| or
                                                 #   convex-hull / alpha-shape surface)
# bands (µm, kPa) from §1
L_pen, d_prolif, L_viable = 150.0, 40.0, 150.0
sig_q_onset, sig_q_full, sig_death = 1.0e3, 5.0e3, 10.0e3   # Pa
tau_death = TAU_DEATH                            # residence time at >=sig_death (set @ Δt_DCM scale)

# --- per-cell update -------------------------------------------------------
for c in cells:
    if c.state == NECROTIC:
        continue                                 # irreversible; skip

    r = norm(c.pos - com)
    d = R_cluster - r                            # rim-depth (distance from surface)

    # (a) NUTRIENT / O2 channel  ------------------------------------------------
    # route A: geometric depth rule (Greenspan constant-rim invariant)
    if   d <  d_prolif:   S_nut = PROLIFERATING
    elif d <  L_viable:   S_nut = QUIESCENT
    else:                 S_nut = NECROTIC
    # route B hook: replace the depth test with local O2/glucose field lookup:
    #   c_O2 = O2_field(c.pos); c_glu = glucose_field(c.pos)
    #   S_nut = NECROTIC   if c_O2 < 0.02 mM or c_glu < 0.06 mM
    #         = QUIESCENT  if c_O2 < p_m (~0.5 mmHg / 4.6e-3 mM)
    #         = PROLIFERATING otherwise

    # (b) MECHANICAL / compression channel  -------------------------------------
    sigma = c.contact_normal_stress              # measured virial/contact stress (Pa), preferred
    # fallback proxy when contact stress is unavailable:
    #   n_c = len(c.contacts)
    #   sigma = stress_from_coordination(n_c)    # n_c<=6 ->~0; 7-11 -> 1-5 kPa; >=12 -> >=5 kPa
    if   sigma <  sig_q_onset: S_mech = PROLIFERATING
    elif sigma <  sig_death:   S_mech = QUIESCENT          # reversible p27-like arrest
    else:                                                   # sigma >= sig_death
        c.t_sigma += dt_DCM
        S_mech = NECROTIC if c.t_sigma >= tau_death else QUIESCENT
    if sigma < sig_death:
        c.t_sigma = 0.0                          # decompression resets the death clock

    # smooth division-rate multiplier from the mechanical channel (p27-like ramp)
    g_sigma = clip((sig_q_full - sigma) / (sig_q_full - sig_q_onset), 0.0, 1.0)

    # (c) COMBINE — worst-of-two  ----------------------------------------------
    S = worst(S_nut, S_mech)                      # severity: PROLIF < QUIESC < NECROTIC

    if S == NECROTIC:
        c.state = NECROTIC                        # irreversible deactivate (mark dead)
        c.divide_rate = 0.0
    elif S == QUIESCENT:
        c.state = QUIESCENT                       # reversible
        c.divide_rate = 0.0
    else:  # PROLIFERATING  (rim AND low-stress)
        c.state = PROLIFERATING
        c.divide_rate = base_rate * g_sigma       # nutrient-gated + stress-graded
```

Design notes (no-abstraction rule):
- **Do NOT impose σ(r).** Keep division rim-biased (it only fires when `S == PROLIFERATING`,
  i.e. in the outer `d_prolif` shell), enforce confinement (outer-shell/ECM boundary at
  physiological modulus), and let HOOMD contact mechanics propagate the load. The inward-rising
  core stress then **emerges** from radial cell anisotropy alone (Dolega 2017) — the
  `contact_normal_stress` read each step is the emergent quantity, not a hand-coded profile.
- **Volume-checkpoint variant (more mechanistic, preferred when cell volume is tracked):**
  replace the `sigma` test with a `ΔV/V₀` test — quiescence when `ΔV/V₀ < −5 %` (the
  Delarue 2014 cellular volume checkpoint at 5 kPa), death when osmoregulation relief
  (Na+/NHE1 efflux capacity) is exceeded. σ and ΔV are equivalent read-outs (osmotic ≈
  mechanical equivalence, Montel/Delarue); pick volume when the enclosed-volume Π term is live.
- **Physiological baseline (project hard rule):** start every cell at its resting setpoint
  (cytoplasm ~65 Pa·s, resting turgor ~40 Pa) — the kPa thresholds are ~100× resting turgor, so
  stress must genuinely build to matter. Never start from a null/floppy shell and expect
  emergent agreement; if route B is used, set the O2/glucose medium boundary to real culture
  values (medium pO2 ~100 mmHg / ~0.2 mM dissolved O2, glucose at DMEM concentration).

---

## 4. Calibration / sanity targets (OVERLAY ONLY — never fit)

The updater is correct if a growing MCF-7 cluster reproduces (do not tune to make these pass;
if they fail, surface to PI for a contract change, do not loosen bands):

1. **All-viable below ~300–400 µm diameter** (MCF-7 at 100 & 300 µm = no dead core, calcein+).
2. **Necrotic-core onset ~400–500 µm**, **prominent core ~700 µm** (≈ MCF-7 day 8; none day 6).
3. **Viable rim ~150 µm and ~constant** as the core grows (Greenspan/Mueller-Klieser invariant);
   only the necrotic core expands.
4. **Necrotic-core fraction-vs-diameter** tracks the BT-474 sibling curve
   (0.20 @ 0.50 mm → 0.29 @ 0.84 mm → 0.61 @ 0.86 mm) — the closest quantitative breast curve.
5. **Whole-spheroid viability time course** ≈ 95 % (d2) → 92 % (d7) → 75 % (d10) → 50 % (d14)
   for MCF-7 (Sajjadi 2025) — progressive once the core forms.
6. **Mechanical channel:** proliferation halves around σ ~5 kPa, near-zero by ~10 kPa, full
   inhibition 6–16 kPa; ~5 % volume loss at 5 kPa; inward-rising core stress profile (Dolega);
   quiescence channel fully reversible on decompression / return to rim.

---

## 5. Source anchors (cited)

Nutrient / diffusion (Q1, Q2):
- Thomlinson RH, Gray LH (1955) Br J Cancer 9(4):539–549, doi:10.1038/bjc.1955.55 — cord
  radius ~200 µm necrosis onset, none <100 µm, viable rim ≤180 µm, O2 diffusion ~150 µm.
- Greenspan HP (1972) Stud Appl Math 51(4):317–340 — three-radius moving boundary; constant
  viable rim, expanding necrotic core.
- Mueller-Klieser W (1997) Am J Physiol 273(4):C1109–23, doi:10.1152/ajpcell.1997.273.4.C1109,
  PMID 9357753 — viable rim ~100–200 µm; rim ~constant, core grows.
- Grimes DR et al. (2014) J R Soc Interface 11:20131124, doi:10.1098/rsif.2013.1124,
  PMID 24430128 — diffusion limit 232±22 µm, OCR 7.29±1.4×10⁻⁷ m³ kg⁻¹ s⁻¹, D_O2 (fit) ~3.8e-9,
  r_l 233 / r_n 155 µm (HCT116).
- Grimes/Bull (2016) PLOS One 11(4):e0153692, doi:10.1371/journal.pone.0153692 — D=2×10⁻⁹ m²/s,
  mitotic O2 cutoff p_m ~0.5 mmHg (std glucose) / ~5 mmHg (low glucose).
- Jiang Y et al. (2005) — necrosis thresholds O2 <0.02 mM, glucose <0.06 mM, lactate >8 mM.
- Lee SY et al. (2010) Oncol Rep 24(1):73–79, PMID 20514446 — MCF-7 ~700 µm/day 8 necrotic
  core, none at day 6; three-zone zonation.
- Sajjadi et al. (2025) PMC12638756 — MCF-7 viability 95→92→75→50 % (d2/7/10/14); necrotic
  layer ~day 6 (5 000-cell); calcein-AM/EthD-1 2/4 µM.
- Nieto C et al. (2026) Eur J Pharm Sci 216:107370, doi:10.1016/j.ejps.2025.107370 (KB paper 02)
  — BT-474 necrotic-fraction curve 0.20/0.29/0.61 at d5/6/7; glucose necrosis threshold ~0.08 mM.
- Front Oncol 6:105 (2016) doi:10.3389/fonc.2016.00105 — MCF-7 <0.5 mm no necrosis, 1.8 mm clear
  core. RSC Analyst (2020) doi:10.1039/d0an00979b — MCF-7 OCR per volume falls with radius.
- ffn KB: papers 02, 04 (viable rim ~100 µm constant), 06 (FUCCI zones), 15 (MCF-7 gaps 5–10 µm,
  cell ~15 µm), 16 (zonation >500 µm, plateau 600–800 µm, rim 20 µm), 21, 27.

Mechanical / solid stress (Q3):
- Helmlinger G, Netti PA, Lichtenbeld HC, Melder RJ, Jain RK (1997) Nat Biotechnol 15(8):778–783,
  doi:10.1038/nbt0897-778, PMID 9255794 — inhibiting stress 45–120 mmHg = 6–16 kPa, line-independent.
- Montel F et al. (2011) Phys Rev Lett 107(18):188102, doi:10.1103/PhysRevLett.107.188102,
  PMID 22107677 — ~10 kPa drastically reduces growth, proliferation loss mainly in the CORE; 5 kPa measurable.
- Delarue M et al. (2014) Biophys J 107(8):1821–1828, doi:10.1016/j.bpj.2014.08.031, PMID 25418163
  — compression → cell-volume drop (minutes) → p27^Kip1 (hours) → late-G1 arrest (days); ~5 %
  volume loss at 5 kPa in ~30 min; conserved across 5 lines; reversible (the volume checkpoint).
- Cheng G, Tse J, Jain RK, Munn LL (2009) PLoS ONE 4(2):e4632, doi:10.1371/journal.pone.0004632,
  PMID 19247489 — high-stress regions suppress proliferation + induce mitochondrial apoptosis (dormancy).
- Dolega ME et al. (2017) Nat Commun 8:14056, doi:10.1038/ncomms14056, PMID 28128198 — internal
  pressure RISES toward the core under compression, from radial cell anisotropy alone (let σ emerge).
- McGrail DJ … Dawson MR (2015) Biophys J 109(7):1334–1337, doi:10.1016/j.bpj.2015.07.046,
  PMID 26445434 — survival under solid stress = active Na+/NHE1 osmoregulation; breast lineage
  validated; death = osmoregulation failure.
- Hadjigeorgiou AG, Stylianopoulos T (2023) Biomech Model Mechanobiol 22:1625–1643,
  doi:10.1007/s10237-023-01716-3, PMID 37129689 — real-tumor residual stress 3.31–10.88 kPa.
- Stylianopoulos T et al. (2012) PNAS 109(38):15101–15108, doi:10.1073/pnas.1213353109,
  PMID 22932871 — solid stress collapses vessels → hypoxia (mechanics⇄nutrient coupling).
- Nia HT et al. (2017) Nat Biomed Eng 1:0004, doi:10.1038/s41551-016-0004, PMID 28966873 — solid
  stress increases with tumor size; confinement a major contributor.
- Boot RC, Koenderink GH, Boukany PE (2021) Adv Phys X 6(1):1978316, doi:10.1080/23746149.2021.1978316
  (KB paper 11) — viable rim ~100–200 µm; three-shell structure.

Caveats (carried from the source notes): MCF-7 necrotic-fraction-vs-size is only
semi-quantitative (the quantitative curve is BT-474, HER2+ sibling, same calcein/PI imaging);
no MCF-7-specific kPa threshold exists but the solid-stress thresholds are explicitly
cell-line-independent, so the 1/5/10 kPa band transfers; MCF-7's compact morphology
(5–10 µm gaps) suggests a rim near the lower end (justifies L_viable=150 over 200 µm);
MCF-7 spheroids loosen by ~7–10 days at some densities, confounding late viability.
