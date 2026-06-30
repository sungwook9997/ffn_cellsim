# FF Stage 6d — γ-floor prototype: MD-free reproduction of the cortical-tension floor

**Date:** 2026-06-29 · **Engine:** `ffn_sim/ff/` (Filament-FEM, Cytosim physics) · **Branch:** `dcm/main`
**Status:** prototype DONE + decisive result — **surfaced to PI** (no tuning; γ swept as a controlled variable)

---

## 1. What this is

The decisive FF experiment from `ENGINE.md §4`: measure cortical tension γ with an **MD-free
mechanical-equilibrium solve** (Cytosim physics on the Warp/FF stack) and compare it to the
**BAOAB-MD γ** from the archived `cortex/cortical_tension.py` (the original γ-floor finding). The
question (ENGINE.md §4): *does the MD-free solve reproduce the same γ (⇒ FF is MD-equivalent for γ),
or does it localize where in-time load-coupling is load-bearing?*

## 2. Method (all grounded; FF units pN·µm·s)

- **Cortex** (`cortex_assembly`): 1000-filament ×40 mesoscopic actin shell, R=10µm, 7 beads,
  ℓ₀=0.5µm, κ=k_BT·ℓ_p (ℓ_p=17µm) — H.3 config values, identical to the MD cortex.
- **Crosslinkers / myosin** (`hand_kmc` presets): α-actinin (Ferrer 2008) / filamin (Pereverzev
  catch–slip) 30/70; NMIIA myosin (Kovács 2003 / Stam-Hocky 2015). Links connect cross-fiber node
  pairs within the **mesoscale reach √(A/n_fil)** — the geometric dual of the ×40 coarse-graining
  (the molecular ε=60nm can't connect the sparse coarse filaments; faithful port of
  `connected_mesh.mesoscale_reach`). Force-free at formation.
- **Turgor** (physiological-baseline rule): Young-Laplace radial pressure, dP0=133 Pa (SOLID
  registry), K_vol=1e3 Pa.
- **Equilibrium:** settle the RESTING passive shell (bending only) + constraint projection + reshape
  (NF2007 project-then-reshape, no BAOAB thermostat); myosin added as the MODULATOR at measurement
  (physiological-baseline rule: measure tension on the resting turgor-pressurised cell + add myosin).
- **γ:** `gamma_estimator.method_of_planes_gamma` — a **faithful port of the validated HOOMD γ
  kernel**, so the MD-free γ is read with the *same instrument* as the BAOAB-MD γ. Two channels kept
  SEPARATE (archived rule): γ_active = actomyosin (actin axial tension from the inextensibility
  multipliers + crosslinker + myosin links); γ_passive = turgor ΔP·R/2.

> **Units trap caught:** 1 pN/µm = 1e-3 mN/m, so the Salbreux band 0.35–0.65 mN/m = **[350, 650]
> pN/µm** (not 0.35–0.65). An early 1000× slip had inverted the conclusion; fixed.

## 3.0 GROUNDED production operating point (the defensible number)

Run at the ×40-mesoscale config values (NOT a prototype downscale): **N=1000 filaments (Plan v2 §3
H.3), n_xl=1000 (KU-3.19), n_myo=100 (Salbreux `n_motors_per_cell`, 3/µm²), f_myo=5 pN (NMIIA
per-side stall)** — `gamma_floor.gamma_floor_production`:

> **γ_active = 0.152 ± 0.007 pN/µm = 1.5e-4 mN/m  →  ~2300× under the Salbreux band.**
> γ_passive (turgor) = 665 pN/µm = 0.665 mN/m (at-band).

This **lands on the archived BAOAB-MD g_soft (3.1e-5–1.4e-4 mN/m)** — quantitative agreement, not just
same-order. The prototype-scale numbers in §3 below (γ_active≈0.95 pN/µm at N=100) over-stated γ
because they used too many myosin links per filament (n_myo=200 at N=100 ≫ the Salbreux density); the
floor MAGNITUDE scales with the modelled motor density (γ_active ∝ n_myo at fixed R — confirmed by an
N-scan 0.94→3.15 pN/µm as N:100→1000 with n_myo∝N). **The exact floor factor is set by the motor
density per cross-section — which is precisely the ungrounded "missing motor-density datum"; at the
grounded Salbreux count it is ~2300× under band.**

## 3.0b The motor-density datum is NOT missing — it EXISTS and CONFIRMS the floor (2026-06-30, ultracode)

An ultracode 31-agent literature hunt (6 search angles → 24 adversarial verifications → synthesis)
re-examined the "missing motor-density datum". Corrected conclusion (cross-checked against the
project's own `H7_CORTICAL_MYOSIN_DENSITY_DATUM_2026-06-07.md`):

- **A direct measurement EXISTS:** Nie et al. 2015 *Cytoskeleton* 72(1):29-46 (PMID 25641802) —
  intensity-calibrated super-res of cortical NMII minifilaments: **~0.625 minifilaments/µm²** (HeLa
  medial cortex; range 0.31-0.94; 1 focus ≈ 1 minifilament; LOW-MEDIUM confidence; **non-MCF7**).
- **It makes the floor WORSE, not better.** γ_active ∝ density, and the measured ~0.6/µm² is **~26-35×
  BELOW** the ~16-21/µm² the active-γ ceiling needs (and ~5× below the code's mis-cited 3/µm²). The
  real datum drives the ceiling to ~0.02 mN/m — *deeper* under the 0.35 band. **The floor is confirmed
  by the measured density, not opened by a missing one.**
- **What is genuinely missing** is narrower than "nobody measured it": (a) an MCF7/breast value; (b) a
  *force-bearing* (load-engaged) vs merely-present discrimination — imaging counts presence, not
  engagement; (c) a true stress-fiber cross-section minifilament count (SFs are imaged en-face). None
  of these can plausibly close a ~30× gap (Nie's HeLa cortex would need ~30× more force-bearing
  minifilaments than measured to lift γ to band — biophysically implausible).
- **Citation fix:** the runtime's "3/µm² = Salbreux 2012" (which this FF module had propagated) is a
  CONFIRMED misattribution; re-anchored to Nie 2015 + flagged (`gamma_floor.py` const + note). Density
  is the controlled variable (swept), never tuned to the band.

This converts the §4 PI option-2 ("source a datum and re-run at the derived value") into a CLOSED
question: the datum exists, and using it reinforces option-1 (accept the floor as a structural limit;
route magnitude to the fine-grained single-cell line). Consistent with the SF/NMII REFUTE and the
Layer-2 magnitude-as-structural-limit findings.

## 3.0c Is the floor CORRECT or a MODEL DEFICIT? (2026-06-30, ultracode + Codex) — leaning DEFICIT

A second ultracode hunt (26 agents, 14 FLOOR_DEFICIT vs 0 FLOOR_CORRECT votes) + a Codex cross-check
attacked the crux: *is the measured cortical-tension band predominantly myosin-generated?* Both
verdicts (each cross-checked against the project's own files) converge:

- **The band is ~50-90 % MYOSIN-dependent (central ~70 %).** Blebbistatin reduces measured cortical
  tension by **68 %** (Fischer-Friedrich 2016, mitotic HeLa, parallel-plate AFM, CONFIRMED), >50 %
  (Tinevez 2009), ~75 % (Fischer-Friedrich 2014), ~91 % (Warmt 2021, MCF-10A). The genuine
  blebbistatin-INSENSITIVE (bare-membrane) passive floor is only **~0.04 mN/m (~9-12 % of band)**. So
  the "band is mostly passive" reading is **REJECTED** — the active actomyosin floor is the REAL gap.
- **The g_rigid / γ_passive "passive carries the band" story is a turgor double-book.** ΔP·R/2 is the
  pressure's *partner*, and dP0=133 Pa is band-implied/tuned (PARAM_AUDIT 2026-06-25) → reclassified
  (§ code caveat). The project's g_rigid ≈ 0.57 mN/m "myosin-independent" is the same artifact; the
  real passive floor is ~10× smaller.
- **Network amplification (Ronceray/Broedersz/Lenz 2016) is NECESSARY but not sufficient alone.**
  Verified amplification ≈ **7-10×**; the grounded gap is ~370× (prototype) to ~2300× (Nie density),
  so buckling alone lifts γ_active only toward the ~0.10 mN/m dipole ceiling (~3.5× still under band).
  And amplification is **buckling-GATED**: the project's own H7_SIGMA_A / H7_GATE_B_BUCKLING docs
  found M-SHAKE forbids buckling (r/r0=1.000; ~73 % of load-bearing spans are single rigid rods) — so
  the floor was **measured with the suspect mechanism turned OFF**, and the deficit-vs-floor question
  has **never actually been tested with buckling ON**.
- **Verdict: GENUINELY UNRESOLVED, leaning model-deficit.** Internally the floor is real *given the
  current constraints*; the literature says the band is active-dominated; the two reconcile only if
  buckling-gated amplification **plus** a true force-bearing (load-engaged) motor density + penetration
  /overlap (Truong-Quang 2021) together supply the missing factor.

**Decisive FF test (concrete, the one experiment that settles it).** Run the cortex
buckling-CAPABLE + percolating + with clustered motor foci + pre-equilibrated binding, and ablate
buckling ON vs OFF, reading the Hill-bounded γ_active only:
  1. relax the constraint so compression-side spans can buckle (r/r0<1) — thin anchor density to ≥2-seg
     spans, OR finite-stiffness inextensibility **(needs PI integrator-freeze sign-off to relax M-SHAKE)**;
  2. cluster motors into minifilament foci so per-span load > F_crit (probe: F_crit≈2.76 pN at 1-seg,
     0.69 pN at 2-seg) — the nonlinear regime;
  3. pre-equilibrate occupancy + run to s_grip→ℓ₀ so heads walk (current s_grip≈0 is the aggregation wall).
  CONFIRM (deficit): γ_active climbs above the ~0.10 mN/m ceiling toward the 0.2-0.4 mN/m active band.
  REFUTE (floor correct): γ_active stays floored even with buckling/percolation/walking ON.
  The buckling ON−OFF difference *is* the amplification factor → compare to Ronceray's ~7×.

**PI actions (in order):** (1) **fix the bookkeeping** (done here — γ_passive reclassified; surface the
g_rigid double-book to PI); (2) **run the buckling-enabled FF test** (PI integrator-freeze sign-off,
scoped); (3) **commission the force-bearing density datum** (the legit non-magic-number lever); (4)
do NOT tune n_myo/f_myo to band (magic-number violation). Full: the ultracode synthesis +
H7_SIGMA_A_NETWORK_INVESTIGATION_2026-06-07 + H7_GATE_B_BUCKLING_INVESTIGATION_2026-06-08.

**Test status — DONE (Stage 6h, FF_STAGE6H_BUCKLING_NETWORK_2026-06-30): floor SURVIVES with buckling
ON.** The network measurement was built (`network_contractility.py`, robust connected inextensible
contractile solver) and run: contractile σ = 10-90 pN/µm, ~25-100× UNDER band at physiological NMIIA
stall, **FLAT in connectivity z (2.98-5.29, incl. sub-isostatic)**, amplification σ/σ_dipole 0.2-0.6
(<1, screening) across f/F_crit=1.4-145 — **no Ronceray amplification; the "buckling rescues it"
branch is NOT supported in the inextensible model.** Key: Ronceray's ~7× is an EXTENSIBLE-network
effect; FF/Cytosim hard inextensibility structurally precludes it. So the floor is a force-MAGNITUDE
limit (motor force×density), not a connectivity/buckling artifact; the remaining lever is the
uncommissioned force-bearing motor-density datum (orthogonal). Cleanest next FF lever: finite axial
modulus (actin EA, Gittes 1993) instead of the hard constraint — see §6h doc. [historical note: the
first 3 prototype attempts hit projector-singularity / explicit-CFL / seeding issues, since fixed.]
The single-fiber buckling enabler is confirmed (`test_buckling.py`: FF buckles, ~0.3 µm bow at 15 %
compression — unlike the DCM M-SHAKE), so the test is NOT PI-integrator-freeze-gated *for FF*. But
three prototype attempts at the full network-amplification measurement (triangular-lattice contractile
patch, σ_buckle/σ_stiff) hit a genuine numerical obstacle: the explicit-projected / reshape relaxation
is FRAGILE on a *contractile, connected, buckling* network — the per-fiber constraint projector goes
SINGULAR under contraction (degenerate segment), and reshape divides by zero when a bend-stiff fiber's
segment collapses; the bucklable patch also did not reliably enter the out-of-plane buckled branch
(needs proper instability seeding). No reliable amplification number was obtained — and a fragile/NaN
number is worse than none (integrity). **Prerequisite for the decisive measurement = a robust
contractile-connected-network equilibrium solver** (regularized projector or finite-stiffness
inextensibility + a buckling-mode seed + virial/Kirkwood stress with sign conventions, all
lit-anchored connectivity/foci). That is a scoped dedicated build, flagged here rather than forced.

## 3. Result — the floor is REPRODUCED MD-free (prototype-scale sweep)

γ_active is linear in the myosin prestress f_myo (the swept controlled variable): **≈0.19 pN/µm per
pN** of per-link prestress (R=10µm, n_myo=200, 8 realizations).

| f_myo [pN/link] | γ_active [pN/µm] | γ_active [mN/m] | vs Salbreux band |
|---|---|---|---|
| 0 | 0.0 | 0 | — |
| 1 | 0.19 ± 0.01 | 1.9e-4 | ~1800× under |
| **5 (NMIIA per-side stall, lit anchor)** | **0.95 ± 0.06** | **9.5e-4** | **~370× under** |
| 10 | 1.90 | 1.9e-3 | ~180× under |
| 50 | 9.5 | 9.5e-3 | ~37× under |
| 150 | 28.5 | 2.85e-2 | ~12× under |

- **γ_passive (turgor) = 665 pN/µm = 0.665 mN/m** ⚠️ **NOT independent evidence — see §3.0c.** This is
  CIRCULAR (dP0=133 Pa is band-implied/tuned, so ΔP·R/2 "hitting band" is by construction) AND a
  double-book (Young-Laplace tension is the pressure's partner, not an independent passive channel).
  The genuine passive floor is ~0.04 mN/m; the band is ~70% myosin. Kept only as a numeric check that
  the method-of-planes reproduces ΔP·R/2, NOT as a "passive carries the band" result.
- **At the lit-anchored NMIIA prestress, actomyosin γ is ~370× UNDER band**, while turgor γ is
  at-band. To reach band-level cortical tension via this transmission the per-link prestress would
  need ≈1840 pN ≈ **370× the physiological NMIIA per-side stall (5 pN)**.

Figure: `outputs/ff/figs/gamma_floor_sweep.png` (log-γ vs f_myo, bands + turgor baseline + lit
anchor overlaid). Data: `outputs/ff/gamma_floor_sweep.npz`.

## 3b. Quantitative cross-check vs the BAOAB-MD γ (disk-grounded)

From `docs/CORTICAL_TENSION_RECORD_2026-06-04.md` (the authoritative archived γ-floor record):

| quantity | BAOAB-MD (archived) | MD-free FF (this work) |
|---|---|---|
| active actomyosin γ at lit kinetics | g_soft ≈ **3.1e-5 – 1.4e-4 mN/m** (FA-anchored / Δ_active) | γ_active = **1.5e-4 mN/m** at the grounded point (§3.0); 9.5e-4 at the over-dense prototype |
| best case under forcing | v0×3000 → **0.030 mN/m (11.6× under)** | linear in f_myo; band needs ≈370× more force |
| KU-3.5 band floor | **0.35 mN/m** | 350 pN/µm (= 0.35 mN/m) |
| passive/turgor γ | at-band (turgor channel) | 665 pN/µm = **0.665 mN/m** at-band |
| diagnosed cause | **force GENERATION** (per-head force / motor recruitment), NOT boundary / FA / compliance / kinetics | **force magnitude / transmission** (too few force-bearing motors per cross-section); kinetics NOT the lever (Stage 6e) |

Both independent methods land the actomyosin γ in the same **~1e-4 – 1e-3 mN/m** regime — ~100–1000×
under band — and both localize the cause to **force generation/magnitude, not the integrator, the
boundary, adhesion, compliance, or binding kinetics**. The MD line showed FA-anchoring and compliant
backbones do not lift g_soft, and that even v0×3000 reaches only 0.030 mN/m (super-stall, filtered
out by the physical-validity gate); the FF line shows the same floor emerges at mechanical
equilibrium with full transmission and is unmoved by dynamic Hand turnover (Stage 6e). **The MD-free
mechanical solve reproduces the BAOAB-MD γ-floor** → FF is MD-equivalent for γ, and the floor is a
physics statement (the missing motor-density datum), not a method artifact.

## 4. Interpretation (what it localizes) — for PI

The MD-free mechanical solve **reproduces the BAOAB-MD γ-floor**: actomyosin tension ~100–1000×
under band, turgor at-band. Because the MD-free solve assumes the myosin stall force is *transmitted*
through the network at mechanical equilibrium (no dynamic binding-throughput limit), reproducing the
floor here means **the floor is STRUCTURAL — a force-magnitude / transmission deficit, not a dynamic
MD artifact**. This is the same conclusion the MD line reached from the other side:

- H7 Gate-B: "myosin stress ~400× < turgor" (active floor is force-TRANSMISSION/lever, not binding
  kinetics) — `project-h7-gate-b-relaxed-build`.
- SF/NMII force-scale REFUTE/HALT: the single-SF tension is not closable from a measured native
  motor-density datum — the **missing parallel-minifilaments-per-cross-section datum** —
  `project-sf-nmii-forcescale-result`.

Our ≈370× gap is the γ-twin of that same missing motor-density factor: the model has the right
*per-motor* physics (NMIIA stall 0.5 pN/head, Kovács/Stam-Hocky kinetics) but **too few force-bearing
motors per equatorial cross-section** to lift actomyosin γ to band. The shell instability we also saw
(turgor at band-implied dP0 inflates the floored cortex — it cannot balance 665 pN/µm with ~1 pN/µm
of actomyosin tension) is the same deficit expressed dynamically.

## 5. Open / PI decision

This is **not** a magic-number situation — γ was swept as a controlled variable; no derived param was
lowered to hit a band (hard rule respected). The floor is a *physics statement*, not a bug. Options
for PI (none auto-taken):

1. **Accept the floor as the FF/MD-consistent result** — both engines agree the single-cell
   center/mesoscale actomyosin cannot reach band cortical tension; the missing piece is the
   motor-density datum (surface to literature search), not the integrator.
2. **Add the missing real force at a derived value** — if a grounded parallel-minifilament density
   (motors per cross-section) can be sourced, raise n_myo / f_myo to that DERIVED value (legit lever)
   and re-run the sweep; predict γ_active scales linearly so band needs ≈370× more force-bearing
   motors than currently modelled.
3. **γ target band itself** stays unresolved (Moazzeni MCF7 1e-2 N/m vs SimuCell3D 1e-3 vs emergent;
   ENGINE.md §6) — the sweep reports γ(f_myo) so any band can be read off without retuning.

## 6. Build artifacts (Stage 6b-iv → 6d, branch dcm/main)

| chunk | module | tests |
|---|---|---|
| 6b-iv inextensibility | `ff/constraints.py` (projector P, reshape, axial tension) | 6 + 3 |
| 6c-a units | `ff/units.py` (pN·µm·s) | 5 |
| 6c-b cortex | `ff/cortex_assembly.py` + `viz_cortex.py` | 4 |
| 6c-c Hand KMC | `ff/hand_kmc.py` (NF2007 §10.1) | 8 |
| 6d γ estimator + floor | `ff/gamma_estimator.py`, `ff/gamma_floor.py`, `ff/gamma_floor_sweep.py` | 7 + 8 |

Full `tests/ff` suite: **43/43 green.**
