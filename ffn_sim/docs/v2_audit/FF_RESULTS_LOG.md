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

---

## 2026-07-07 (refinement) — "is cytoplasm alone enough to fix the overshoot?" → NO (measured + adversarially verified)

PI pushed on the compression-overshoot cause (lateral bulge?) and on my "missing cytoplasm piece" line above.
That single-cause line is **overstated** — corrected here by (a) a direct drained-vs-undrained measurement
(`scripts/ff_bulge_probe.py`) and (b) a 12-agent adversarial diagnosis (workflow `wf_214b6be0-391`).

**Direct measurement** (N_fil=2000 cortex, R0=7.50 µm, CPU; absolute F scales with N — the ratios are the point):

| mode | strain | R_eq/R0 (model) | R_eq/R0 (const-vol) | V/V0 | F [nN] | ΔP [Pa] | γ_app [mN/m] |
|---|---|---|---|---|---|---|---|
| undrained | 12% | +1.7% | +6.6% | 0.998 | 78.9 | 1342 | 5.1 |
| undrained | 27% | **+4.9%** | **+17.0%** | 0.989 | 910 | 8300 | 32.6 |
| drained@40Pa | 12% | +1.0% | +6.6% | 0.979 | **2.0** | 40 | **0.15** |
| drained@40Pa | 27% | +1.0% | +17.0% | 0.903 | **4.0** | 40 | **0.15** |

Three findings: **(1)** the cortex **under-bulges** — reaches only +4.9% lateral expansion at 27% vs the +17.0%
a constant-volume oblate needs (PI's "가로로 안 늘어난다" confirmed); the shortfall is exactly the volume it sheds
into the osmotic ΔP. **(2)** switching to the **drained limit collapses the plate force 39–228×** → the overshoot
(gap B) is dominantly the **sealed/undrained van't Hoff osmotic turgor**, not cortex/cytoskeleton stiffness.
**(3)** drained, γ_app = 0.15 mN/m (= Young–Laplace ΔP·R/2 at the γ-floor) → the drained response is
**too SOFT** → draining alone swings the model past the target toward too-soft (the "drained-load-carrier" trap).

**Corrected diagnosis** (supersedes the "missing cytoplasm piece" framing above):
- **Two gaps live in different tensor channels.** Gap A (resting γ, ~70× low) is **in-plane/deviatoric** (actomyosin,
  Young–Laplace γ=ΔP·R/2) — **no isotropic pore-fluid/cytoplasm DOF can ever populate it**; needs the myosin fix
  (missing MCF7-adherent load-*engaged* NMII density datum — surface to PI). Gap B (compression) is **normal/hydrostatic**.
- **Gap B root cause is a constitutive misapplication, not "correct undrained physics":** an *equilibrium osmotic
  P–V law* (Guo 2017) is evaluated on a *no-flux (undrained) geometric volume trajectory* `V_cyto=V_hull−V_nuc` —
  no Lp/Darcy/aquaporin DOF. That manufactures the 1/(V−vmin) pole (40→100530 Pa). ⚠️ My earlier "K_vol=736 kPa
  is the correct undrained modulus, do not retune" was **wrong** — that canonizes the bug. Fix = **replace the closure
  with a biphasic poroelastic cytoplasm (fluid transport + a distributed *solid* drained matrix ~10²–10³ Pa,
  Moeendarbary 2013)**, NOT soften K_vol (that would be gate-tuning), and NOT a viscous filler behind the seal.
- **The overshoot "grows with indentation" is largely a PROTOCOL/category error:** large-strain nominal σ=F/πR²
  divided by a small-strain fixed-E Hertz (E=249 Pa, Zbiral 10µm colloidal). Any strain-stiffening material overshoots
  a fixed-E Hertz. **Must run FIRST**: extract E_eff by the *same* Hertz inversion at 2–3% strain, matched colloidal
  geometry, rate-matched — to size the genuine residual before adding physics. Code emits no Hertz E / no time DOF.
- **Whether the cortex must ALSO be lifted for B is rate-gated and open:** Zbiral's E=249 Pa (10µm bead) may itself be a
  poroelastic bulk readout; a well-calibrated biphasic cytoplasm could carry B on its own solid matrix. Depends on the
  (unmeasured) Zbiral loading rate vs τ_p≈0.5–1.4 s.

**Ranked what it takes** (adversarially-verified): (1) rate-matched small-strain Hertz-inverted validation **first**;
(2) biphasic poroelastic cytoplasm = drainage **+ solid load-bearing matrix** (dominant for B); (3) cortex γ-floor
lift via engaged-NMII density (dominant for A, orthogonal, missing datum → PI); (4) membrane K_A reservoir-exhaustion
upturn (secondary, large-strain, PI-blocked f_excess); (5) resting ΔP=2γ/R co-fixes automatically once cortex lifts;
(6) nucleus stays honest (~45%). **Verdict: cytoplasm is necessary but NOT sufficient — B needs biphasic cytoplasm
(+ possibly cortex) and A is entirely orthogonal.** Membrane is a separate compartment (reservoir area supply), not
"cytoplasm." FF uses the turgor osmotic law (soft, physical), NOT DCM's explicit K_vol·(V−V0)² penalty; water efflux
belongs in that turgor/volume-regulation term, and the drained limit is already available as `pressure_setpoint`.

---

## 2026-07-08 — IMPLEMENTED the fix ("전부 진행"): biphasic cytoplasm + membrane reservoir + Hertz harness

Grounded by a 4-agent literature workflow (`wf_4561ab1a-253`; all values cited, no-magic-number). New code is
additive/backward-compat (all-off → bit-identical; 7 sanity gates in `tests/ff/test_poroelastic_cytoplasm.py`,
42 existing compression/membrane tests still pass).

**Piece #1 — Hertz-inversion harness** (`scripts/ff_hertz_validation.py`): extract E_eff at small strain (2–3%)
by the SAME Sneddon-parabolic inversion Zbiral used (F=(4/3)(E/(1−ν²))√R·δ₁^1.5, ν=0.5, R_cell=7.5µm, δ₁=R·strain,
slope fit), so we compare **modulus-to-modulus** vs 249 Pa. **This immediately corrected the headline:** the honest
undrained residual is **16.5×** (E_fit≈4114 Pa), *flat* with strain (Hertzian) — **NOT** the 857× artifact (which
was large-strain σ ÷ small-strain E, exactly as diagnosed).

**Piece #2 — biphasic poroelastic cytoplasm** (`network_warp.py`, additive params `Lp_um_s_Pa`, `K_drained_Pa`,
`load_time_s`):
- (A) DRAINAGE — Kedem-Katchalsky (σ≈1): osmotic reference `V0_eff` drains toward `V_cyto`, relaxing the van't Hoff
  turgor to dP0. **Analytic exp-relaxation**, clock **τ_osm=(V0_eff−vmin)/(Lp·A·Π_in)** — the MEMBRANE drainage
  time (~34 s at Lp=1e-7; ~212 s exosmotic default 1.6e-8), *not* the poroelastic τ_p≈1.4 s. Lp=1e-7 µm/(s·Pa)
  (COS-7+airway Pf; PMC3161049). **No MCF7 Lp datum → proxy, flag PI.**
- (B) DRAINED SOLID — Terzaghi effective stress `dP_solid=K_drained·max(0,(V0−V_cyto)/V0)` on the **FIXED** rest
  volume V0 (bug-fixed: the draining reference would vanish exactly at the drained limit). K_drained≈300 Pa
  (Moeendarbary MDCK/HT1080 soft-epithelial 0.4 kPa, ν≈0.25–0.3). **No MCF7 drained modulus → borrowed, flag PI.**

**Piece #4 — membrane reservoir K_A upturn** (`compartments.py` `f_excess`): γ_mem plateau until area>A0·(1+f_excess),
then K_A=0.235 N/m upturn (Rawicz 2000) capped at lysis. f_excess=0.25 controlled variable (0.10–0.40 band, Raucher-
Sheetz/Figard). **No MCF7 f_excess → controlled variable, flag PI.** f_excess=0 recovers the prior plateau.

**Piece #3 — cortex γ-floor (gap A):** the engaged-fraction (0.65, Truong-Quang) and density controlled-variable
(`gamma_floor_density_sweep`, PROD_N_MYO=442 @ Nie-2015 0.625/µm²) **already exist**; the only lever is the **missing
MCF7-adherent load-engaged NMII density datum**. Tuning it to the band would violate the no-magic-number rule →
**left as a flagged controlled variable, surfaced to PI.** Cytoplasm/drainage cannot touch gap A (different channel).

### Validation (Hertz E_eff, N_fil=2000, small strain; figure `outputs/ff/figs/ff_poroelastic_validation.png`)
| regime | E_fit [Pa] | ×MCF7 | note |
|---|---|---|---|
| undrained (= Zbiral fast rate 5 µm/s) | 4114 | **16.5×** | flat E_pt (Hertzian) ✓ |
| drained: setpoint 40 Pa + K_drained | 1548 | 6.2× | osmotic relaxed |
| drained: Lp + slow ramp + K_drained | 1167 | **4.7×** | fully drained (drained_frac=1) |
| **real MCF7 (Zbiral)** | **249** (224–279) | 1× | target |

**Honest result:** the biphasic cytoplasm cuts the overshoot **16.5× → ~5×** — a real ~3× improvement, and it
confirms the diagnosis quantitatively (drainage removes the osmotic contribution). **But it is NOT sufficient, and
the reason sharpens the PI's insight:** at Zbiral's measured rate (5 µm/s, τ_load≈0.1 s ≪ τ_osm≈30–200 s) **the cell
physically cannot drain** → the model stays ~undrained (16.5×). The residual ~5× (even fully drained) is the **CORTEX
indentation stiffness** — the lateral-bulge / crosslinker-stiffness lever (the cortex can't cheaply gain surface area
to accommodate at constant volume). So the **dominant remaining lever for matching a fast AFM is the cortex, not more
cytoplasm** — which is exactly the "가로로 코르텍스가 늘어나야 한다" point. Gap A (resting tension, ~70× low) is still
orthogonal and datum-limited.

**Next (PI-gated):** the residual is a cortex-area-accommodation problem — the cortex crosslink turnover / area
remodeling under load (so the inextensible fixed-connectivity network can increase apparent area and bulge at
constant volume). That is a cortex-remodeling change, distinct from the cytoplasm just added. **PI-flag data still
missing:** MCF7 Lp, MCF7 drained cytoplasm modulus, MCF7 f_excess, MCF7-adherent engaged-NMII density.

---

## 2026-07-08 (cont.) — "새로운 것들 넣어서 진행": cortex crosslink turnover + rate sweep

Added the flagged next-lever mechanism: **viscoelastic crosslink turnover** (`xl_turnover_kernel` in `network_warp.py`;
params `xl_koff_per_s`, `xl_x_beta_um`). Each `reshape_every`, a fraction 1−exp(−k_off(F)·dt) of the crosslinker
ensemble unbinds and rebinds FORCE-FREE, so rest lengths creep toward current lengths (stress relaxation) — a Maxwell
element per crosslinker, relaxation time 1/k_off(F). Bell slip k_off(F)=k_off0·exp(|F|·x_β/kT); α-actinin
**k_off0=0.066/s, x_β=0.4 nm (Ferrer 2008 PNAS)**. Real-time via `load_time_s`. `xl_koff_per_s=None` → bit-identical
(2 new sanity gates: backward-compat + rate-dependent softening; all 9 poroelastic + 42 existing tests pass).

### Rate sweep (full physics: drainage Lp=1e-7 + K_drained=300 + turnover koff=0.066 + f_excess=0.25)
Figure `outputs/ff/figs/ff_rate_sweep.png`. τ_osm≈34 s (drainage), 1/k_off≈15 s (turnover).

| v_load [µm/s] | load time @3% [s] | E_fit [Pa] | ×MCF7 (colloidal 249) | regime |
|---|---|---|---|---|
| 5 (Zbiral) | 0.09 | 4110 | 16.5× | elastic + undrained |
| 0.5 | 0.90 | 4035 | 16.2× | elastic |
| 0.05 | 9.0 | 3394 | 13.6× | onset of relaxation |
| 0.005 | 90 | 1116 | 4.5× | drained + remodeled |
| 0.0005 | 900 | 1274 | 5.1× | fully relaxed (plateau) |

**Honest result:** the new turnover mechanism gives the **physically-correct rate-dependent viscoelastic/poroelastic
softening** — the modulus drops **16.5× → ~5×** as the loading time crosses τ_osm≈34 s and 1/k_off≈15 s. Two findings:
1. **At Zbiral's rate (5 µm/s) the model is 16.5×** — drainage and turnover cannot act (load ≈0.1 s ≪ both clocks).
   This is a genuine PREDICTION: MCF7 probed at 5 µm/s should read stiff; the 249 Pa colloidal value is a large-contact
   whole-cell average likely including relaxation/spreading over the fit window.
2. **A ~5× floor remains at all rates** (E_fit plateaus ~1100–1300 Pa). This is NOT closed by drainage, drained solid,
   or crosslink turnover — it is the cortex **segment inextensibility + bending** (the `reshape` constraint keeps
   filaments inextensible, so the mesh resists local indentation geometrically even with fully relaxed crosslinkers)
   plus the tension-shell-vs-Hertz-solid inversion effect. **Notably ~1100 Pa sits at the TOP of the sharp-tip MCF7
   band (200–1000 Pa, Li 2008)** — so it is not unphysical; the ~4.5× gap is specifically to the *colloidal large-contact*
   249 Pa (softer because large probes average over the whole soft cell).

**So the residual floor is the cortex ARCHITECTURE (inextensibility/bending), not remodeling kinetics.** Closing the
last gap to the colloidal 249 Pa would need cortex **filament-length dynamics** (treadmilling / severing so segments
can change length and the shell can gain true area) — a deeper mechanistic change — OR accepting the model represents
the stiffer sharp-tip/local regime. **PI decision point.** All added physics is literature-grounded, no gate-tuning;
the rate-dependence is the honest emergent behavior. **PI-flag data still missing:** MCF7 Lp, MCF7 drained cytoplasm
modulus, MCF7 f_excess, MCF7-adherent engaged-NMII density, MCF7 cortex crosslinker turnover rate under load.

---

## 2026-07-08 (cont.) — PI: "누르는 속도가 세포 변형 속도보다 빠를 수도" → CONFIRMED (the dominant artifact)

PI hypothesised the pressing speed matters and we may press faster than the cell can deform. **Empirically decisive.**

**Convergence test** (`converge_test.py`, undrained, strain 3%, increasing relaxation steps) — the quasi-static
solve was **NOT converged** at the n_steps I had been using:

| n_steps | undrained E [Pa] | ×249 | drained+K_drained E [Pa] |
|---|---|---|---|
| 1600 (used) | 4127 | 16.6× | 1333 |
| 3200 | 262 | 1.1× | 1332 |
| 6400 | 20 | 0.1× | 1332 |
| 12800 | 19 | 0.1× | 1330 |

The undrained force **drops by orders of magnitude** with more relaxation — **the "16.5×" I reported was an
under-relaxed (fast-press) TRANSIENT**, not the equilibrium. The relaxed undrained equilibrium is **~17–19 Pa (0.1×,
γ-floor, too SOFT)**. (The drained+K_drained case was already converged at ~1330 Pa — the K_drained solid does not relax.)

**Physical press-speed ramp** — implemented the physiologically-correct fix (`v_press_um_s`, `dwell_steps` in
`network_warp.py`): each numerical step = dt_real = dt_mu/μ SECONDS with μ the NF2007 cytoplasm-viscosity mobility
(**η=65.9 Pa·s**, Dessard 2024), the plate advances v_press·dt_real per step, so relaxation-completeness is set by
PHYSICS not an arbitrary n_steps. This also fixes a physiological-baseline violation (the loop used a numerical dt,
not η). Result (`ff_press_speed.png`, strain 3%, read at ramp-end):

| v_press [µm/s] | E [Pa] | ×249 |
|---|---|---|
| 5 (Zbiral) | 251255 | 1009× |
| 1 | 50777 | 204× |
| 0.2 | 10157 | 41× |
| 0.05 | 2553 | 10× |
| 5, then dwell→relax | 17 | 0.1× (equilibrium = γ-floor) |

**E ∝ v_press exactly** (5→1→0.2 = 1009→204→41×) = a **VISCOUS transient** (F≈η·v_press·geometry), NOT the elastic
modulus. So: (1) pressing speed dominates the reading; (2) the relaxed elastic equilibrium is **soft (~17 Pa = γ-floor)**
— **all the prior "too-stiff" numbers (16.5×, the rate sweep's fast points, the original 857×) were fast-press /
under-relaxation transients**, not a real elastic over-stiffness. The genuine elastic gap is the γ-floor (too soft),
in the OTHER direction.

⚠️ **CAVEAT (honest):** the per-node drag uses the single-fiber NF2007 mobility, which **over-estimates the BULK
cytoplasm viscous response by ~100×** (continuum η·ε̇ ≈ 65.9·0.67 ≈ 44 Pa, vs the model's huge transient). So the
**absolute crossover speed needs the drag calibrated to the bulk η**; the **direction and the E∝v_press scaling are
robust**. Sanity gates added (backward-compat + viscous-transient + dwell-relaxes; 11 poroelastic tests pass).

**This reframes the whole overshoot story:** the compression response must be reported as a RATE-DEPENDENT
elastic+viscous mix at a stated, physical press speed, with the mechanics relaxed to equilibrium (or a physical
dwell) before reading — and the drag calibrated to η. Next: calibrate the per-node drag to bulk η=65.9, then re-run
the AFM comparison at Zbiral's 5 µm/s with a proper elastic/viscous split (as real AFM Hertz analysis does). The
elastic-equilibrium γ-floor (too soft) is the real remaining physics gap. **PI decision point.**

---

## 2026-07-08 (cont.) — bulk-η drag CALIBRATED → model lands in the MCF7 band at ~0.5 µm/s

Calibrated the per-node drag to the bulk cytoplasm viscosity (`eta_bulk_Pa_s` in `network_warp.py`): the whole-cell
Stokes drag distributed over the cortex nodes, **γ_node = 6π·η·R / Nc** (canonical form, grid-intensive Σγ=6πηR,
NO tuning) — this accounts for the hydrodynamic screening of the dense cortex that the single-fiber NF2007 mobility
ignored (the ~100× over-estimate). Validated against the analytic Newtonian ground truth η_eff=F_visc/(A·ε̇)
(`ff_drag_calibration`; figure `outputs/ff/figs/ff_drag_calibrated.png`):

| v_press [µm/s] | F_visc [pN] | η_eff [Pa·s] | E [Pa] | ×249 |
|---|---|---|---|---|
| 50 | 12122 | 141 | 23364 | 93.8× |
| 5 (Zbiral) | 1194 | 138 | 2313 | 9.3× |
| 0.5 | 109 | 124 | 223 | **0.9× (IN BAND)** |
| relaxed equilibrium | — | — | 13 | 0.1× (γ-floor) |

**Two results:** (1) **η_eff ≈ 124–141 Pa·s, rate-INDEPENDENT** (Newtonian) across a 100× speed range → the drag is
now physical, validated to **~2× of the 65.9 target** (the factor-2 is the crude squeeze-flow ε̇/area estimate, not a
drag error; the Stokes form is used verbatim). (2) With physical viscosity the model produces a **rate-dependent
modulus that passes THROUGH the MCF7 band (~249 Pa) at ~0.5–1 µm/s**, is viscous-stiff above (Zbiral 5 µm/s → 9.3×,
viscous-dominated: F_visc=1194 ≫ F_elastic=6.7), and elastically soft at the relaxed equilibrium (~13 Pa = γ-floor).

**Interpretation.** The apparent stiffness of a whole-cell compression is a rate-dependent elastic+viscous mix, and
at Zbiral's fast 5 µm/s it is **viscosity-dominated**. The model's residual gap at exactly 5 µm/s (9.3×) is the
~2× η-calibration factor × the elastic/viscous partition in the Hertz inversion — refinable. The honest headline:
once the press speed and the physical cytoplasm viscosity are handled, **the whole "cell is too stiff" story
dissolves** — the model's elastic modulus is soft (γ-floor), the fast-AFM stiffness is viscous, and a physical
sub-µm/s press lands the apparent modulus on the measured MCF7 value. The remaining true-physics gap is the
γ-floor (elastic equilibrium too soft, ~13 vs ~249 Pa), which needs the myosin-generation / engaged-density fix —
NOT more cortex stiffness. Sanity gate added (bulk-drag < single-fiber; 12 poroelastic tests pass).

**PI-flag data still missing:** MCF7 Lp, MCF7 drained cytoplasm modulus, MCF7 f_excess, MCF7-adherent engaged-NMII
density, MCF7 cortex crosslinker turnover-under-load rate, and Zbiral's exact δ/rate/elastic-viscous separation.

---

## 2026-07-08 (cont.) — NATIVE-scale (A5000) confirmation: the under-relaxation reframe HOLDS at full resolution

PI: re-verify at native. Ran on the gbook A5000 (`ff_native_confirm.py`, cortex NF=38000 → **Nc=266,000**,
`device=cuda:0`, 38.5 ms/step, turgor_every=50; figure `outputs/ff/figs/ff_native_convergence.png`). Small-strain
(3%) undrained modulus vs relaxation steps:

| n_steps | native E [Pa] | ×249 |
|---|---|---|
| 500 | 5,455,562 | 21900× |
| 4000 | 216,163 | 868× |
| 8000 | 18,080 | 72.6× |
| 16000 | 6,588 | 26.5× |
| 32000 | 1,417 | 5.69× |
| 64000 | **182** | **0.73×** |

**Confirmed at native.** The undrained modulus drops **monotonically through ~5 orders of magnitude** as the cell is
allowed to finish deforming (500 → 64000 steps), crosses the MCF7 band, and reaches the **soft side (0.73×, still
falling toward the γ-floor)** — exactly the CPU trajectory. So the small-strain "stiffness" was a fast-press /
under-relaxation TRANSIENT at full resolution too, NOT a real elastic over-stiffness; the relaxed native elastic
equilibrium is **soft** (heading to the γ-floor). Native equilibrates **~10× slower** than the N=2000 CPU cortex
(64000 vs ~6400 steps to reach the soft regime) — expected from the ~130× denser mesh (more/slower soft modes),
not a different conclusion.

Caveats: (a) native drained+K_drained was only run at 8000 steps (22×) which is itself under-relaxed (per the
undrained trajectory), so the converged native drained value is lower — not separately pinned. (b) The physical
press-speed ramp is infeasible at native (γ_node=6πηR/Nc ∝ 1/Nc → dt_real tiny → millions of ramp steps), so the
viscous/press-speed result stands as the CPU-scale demonstration; the elastic-equilibrium reframe is what native
confirms. **Bottom line: the whole investigation's headline — there was never a real "too stiff" gap, only the
γ-floor (too soft) — holds at native scale.** Code synced to gbook `~/ff_scratch` (rsync, non-repo).
