---
archived_on: 2026-07-28
superseded_by: aleph/docs/v2_audit/AC_EXECUTION_PLAN_2026-07-25.md
reason: >
  Written BEFORE the 2026-07-25 PI reframe, i.e. for a different objective — forward prediction
  and magnitude matching, rather than inferring per-cell-type parameters with the gate on
  mechanical connectedness. Archived, not deleted: its measurements and reasoning stand as a
  record of what was true then. Nothing in it may be quoted as current state; STATE.md is that.
  Selected mechanically: pre-reframe AND cited by no live file (code, STATE.md, CLAUDE.md,
  cell_engine/, gate_contracts/, tests, Makefile). Citations from run outputs and from other
  pre-reframe documents were not treated as protective.
---

# Resting bound-myosin setpoint — sourcing dossier + force-budget reframing (2026-07-23)

Branch `codex/ff-ac-codex`. Fresh Lead session. PI directive this session: **do NOT close the resting gate
fast** ("절대 완료를 함부로 하지 마 … 너무 빠르게 닫아버리면 안돼"). This dossier is the rigorous sourcing +
force-budget analysis behind keeping the gate OPEN. It supersedes the framing in
`AC_ENGINE_COMPLETION_ROADMAP_2026-07-23.md §0` that the gate is blocked on "ONE PI-GAP = fraction + per-head
force." **Nothing here closes the gate.** It reframes what closing it actually requires and surfaces the
decisions to the PI.

## ⚠️ CORRECTIONS (2026-07-23, post Codex adversarial review + native A5000 run) — READ FIRST

An independent Codex read-only review + a native A5000 run overturned three of this dossier's conclusions. All
verified against code/data (not rubber-stamped). **Net: the gate is HARDER than the body of this doc says.**

- **RETRACTED: "physiological content is sufficient in total (1.2–4.8×)" (§5 table).** The force budget
  `N_mf·N_side·f_head·duty` is a whole-cortex scalar-capacity SUM; the physically-relevant quantity is the
  method-of-planes CUT force. Reproduced: in the oracle geometry the total scalar (168,960 pN) over-counts the
  cut force (6,362 pN) by **26.6×** — so the ratios in §5 do NOT establish sufficiency. Sufficiency is
  UNDETERMINED and leans deficient. Only "native decides conversion" survives.
- **CORRECTED decision E (§3a/§9): the modulus-vs-tension walk-back was WRONG.** The repo already carries a
  sourced MCF7 cortical tension `MCF7_IQR_PN_UM = (180, 400)` pN/µm (**Hosseini 2020**, `ff/gamma_estimator.py:41`)
  + an active-fraction framework (0.70× total band). The myosin target is ~0.70·(180–400) = **126–280 pN/µm**, and
  must be reconciled with the 40 Pa turgor-Laplace (140). The target is likely HIGHER than 140, not lower.
- **NATIVE-CONFIRMED cap (Codex finding 7): the seeding cannot reach physiological duty.** At native (70,686) with
  production `capture_radius = 0.05 µm`, only **1,506/8,840 = 17.0%** of heads are eligible (heads sit ~0.2 µm off
  the actin; the code's own CUDA test overrides capture to 0.6 µm). Requested duty 0.3 → realized **0.170**;
  `eligible_shortfall = True`. And the native outer step **rolled back (not stable)** — the seeded run did not even
  balance. So physiological duty is un-seedable until head-placement / capture is fixed.
- **Recon miss:** `ac/motor/params_i0b3.yaml` (not read in recon) already records: **NM2B-pure first baseline** (PI
  2026-07-21), N_side=10 to retire (~28 heads, Nagy 2013), the k_off0=0.35 ADP-release **mislabel** (real detach
  0.4/s), and `cortical_minifilament_density` as a listed GAP. Decisions A/C′/F below were partly pre-decided there.

The rest of this doc is retained as written for provenance; where it conflicts with the four points above, **the
corrections win.**

## TL;DR (honest gate status = OPEN)

1. The bound-myosin cortical-tension **mechanism is real** (built `9bb91e0a`, re-verified this session on the CPU
   oracle: membrane residual 40.74 → 2.17 pN, minimum exactly at the Laplace target γ_cortex = ΔP·R/2 − γ_mem =
   140 pN/µm).
2. The handoff's "single PI-GAP = fraction + per-head force" **understates the problem**. Closing the gate needs
   the whole **force-magnitude / density SET** at physiological values: heads-per-minifilament (code is ~3× under),
   per-head force (an 11× PI-gated spread), the **load-dependent** duty ratio, and the minifilament density.
3. The CURRENT placeholder content is robustly force-deficient (the historical ~500× floor). **[SUPERSEDED — see
   Corrections]** an earlier draft claimed physiological content is "sufficient in total (1.2–4.8×)"; that used a
   scalar-capacity sum that over-counts the method-of-planes cut force by ~26.6× (reproduced), so **sufficiency is
   UNDETERMINED and leans deficient**. Native decides.
4. **[SUPERSEDED — see Corrections]** decision E is a real open item: the repo carries a sourced MCF7 tension
   (Hosseini 2020, 180–400 pN/µm) + active-fraction 0.70 ⇒ myosin target ~126–280 pN/µm, likely **higher** than the
   140 pN/µm turgor-Laplace value. The earlier "modulus-confusion, target credible" walk-back was wrong.
5. The open question is **conversion** (does the force become γ prestress?) — and now also **seedability**: at
   native, production capture (0.05 µm) caps the seeded duty at 0.17 and the seeded step rolls back. Both are
   **resolvable only on the gbook A5000** (native run 2026-07-23: realized 0.170, `outer_rolled_back=True`).
6. A **deeper fidelity question** surfaced (§6): the code models NMII head–actin as a **pure Bell slip** bond
   (off-rate rises with load → duty falls), but the primary literature (Kovács 2007) shows resistive load
   **increases** duty for tension maintenance (slowed ADP release). Wrong-sign load-dependence may be the root
   cause of BOTH the emergent γ-floor AND why a seeded workaround is needed at all.

## 1. The mechanism (recap) and its CPU-oracle calibration

Bound NMII heads carry isometric tension → the cortex contracts → the radial (un-preloaded) ERM stretches →
pulls the membrane inward against turgor. Verified this session (`coupled_resting_balance`, subdiv-3):

| seeded per-dipole tension (pN) | emergent γ_cortex (pN/µm) | membrane residual mean (pN) |
|---:|---:|---:|
| 0 (OFF) | 0.00 | 40.74 (full turgor floor) |
| 80 | 124.8 | 4.86 |
| **88** | **137.3** | **2.17** (min) |
| 90 | 140.5 | 2.31 |
| 120 | 187.4 | 13.36 |

The residual bottoms out where the emergent γ_cortex reaches the Laplace target — the mechanism, reproduced. **But
the 88 pN/edge that achieves this is a TEST knob, not a physiological datum.** §5 shows what physiological myosin
actually supplies.

## 2. Native cortex myosin geometry (recon, verified against code)

| quantity | code value | file | physiological? |
|---|---|---|---|
| actin filaments | 70,686 | `assemble.py:167` | ✅ (100/µm²·4π·7.5², KB-3.18) |
| minifilament areal density | 0.625/µm² → **442 minifilaments** | `gamma_floor.py:116-120`, `weave.py:155` | Nie 2015 (PMID 25641802), **HeLa medial cortex** proxy, range 0.31–0.94. ("Salbreux 3/µm²" was a confirmed misattribution.) MCF7 value unknown. |
| heads per minifilament | **N_side=10 → 20 heads/mf** | `assemble.py:141` | ⚠️ **~3× UNDER.** Billington 2013 (PMID 24072716) EM: ~28–30 molecules/mf → ~56–60 heads. Flagged "AFINES claim_a; **N_side is a PI GAP**". The lumped `ff/myosin_linear` already uses `N_HEADS_MINIFIL=60`. |
| total native heads | **8,840** (442·20) | `assemble.py:359` | physiological would be ~26,520 |
| per-head force f_stall | **0.5 pN** | `assemble.py:148` | ⚠️ **low end**, "AFINES claim_a; GAP" — see §4 conflict |
| k_xb (crossbridge) | 1.0e3 pN/µm | `assemble.py:145/368` | MASTER knob, GAP (band 100–1000) |
| capture_radius | 0.05 µm | `assemble.py:154/371` | provisional |

## 3. γ estimator formula (method-of-planes)

`cortical_tension.measure_cortical_stress` → `ff.gamma_estimator.method_of_planes_gamma`:

    γ = mean over n̂ of | (1/2πR) · Σ_{elements crossing n̂} T · |û·n̂| |

Linear in the myosin force budget. Each seeded bound head is one dipole carrying `f_head`. This is the
physically-correct cut-force / stress that fixed the legacy `~1015×` radial-self-cancellation artifact — so the
force budget below is NOT a metric artifact.

### 3a. The target is a defensible resting-cell value; the MCF7 "conflict" is probably a quantity confusion (decision E)

According to PubMed, the 40 Pa / ~150 pN/µm target is a **self-consistent, sourced Laplace pair**, not an
artifact:
- **Fischer-Friedrich 2014** (PMID 25169063, DOI 10.1038/srep06213), the exact source of the code's 40 Pa:
  interphase HeLa has **≈40 Pa and 0.2 mN/m** (metaphase ≈400 Pa, 1.6 mN/m). Via Laplace γ=ΔP·R/2, 40 Pa at their
  R≈10 µm ⇒ 0.2 mN/m = 200 pN/µm; scaling to the code's R=7.5 µm ⇒ 0.15 mN/m = 150 pN/µm (140 after −γ_mem). So
  the target MATCHES the measured resting surface tension. The mN/m-range tensions are METAPHASE, not resting.
- **Stewart 2011** (PMID 21196934, DOI 10.1038/nature09642) validates the MECHANISM directly: "osmotic pressure
  is balanced by inwardly directed actomyosin cortex contraction … the actomyosin cortex is required to maintain
  this rounding force." This is exactly the turgor-vs-myosin-cortical-tension balance the resting_setpoint models.
- **Chugh & Paluch 2018** (PMID 30026344, DOI 10.1242/jcs.186254, review): "Myosin-generated forces create
  tension in the cortical network," modulated by cortex composition + organization.

⇒ The target is a **defensible resting-cell cortical tension**: self-consistent as a HeLa interphase proxy
(40 Pa ↔ 0.2 mN/m, Fischer-Friedrich), and 0.14–0.2 mN/m sits inside the general resting-cortex band 0.1–1 mN/m.
The mN/m-range values (1.6 mN/m) are METAPHASE, not resting.

**On the apparent conflict (calibrated this session — walking back an over-alarm I folded in from the workflow):**
the sourcing workflow, trusting the unverified project-memory note "Moazzeni MCF7 ~1e-2 N/m = 10 mN/m," flagged a
10–70× conflict. This session's own PubMed check finds MCF7 mechanics papers report an **elastic MODULUS, not a
surface tension** (Omidvar 2014, DOI 10.1016/j.jbiomech.2014.08.002 — MCF7 local stiffness + adhesion force; no
MCF7 cortical *surface tension* found). A 10 mN/m interphase cortical tension would exceed even HeLa metaphase
(1.6 mN/m), which is implausible for a resting cell. **So the Moazzeni note is most likely a modulus (or apparent
AFM stiffness), NOT a turgor-balanced surface tension, and the conflict is probably a quantity confusion, not
real.** Net: the 140 pN/µm target likely stands, and the §5 sufficiency conclusion is probably robust (it would
only flip if a DOI-verified MCF7 resting tension turned out ≫0.2 mN/m). **Genuine (but non-catastrophic) gap:** no
DOI-verified MCF7-specific resting cortical *surface tension* exists in what I found — PI to supply one (and
confirm it is a tension, not a modulus) or accept the HeLa-interphase proxy. Do NOT back-solve ΔP from 140 pN/µm.

## 4. Per-head force — the 11× PI-gated conflict (unchanged, still needs a PI call)

| value | source | context |
|---|---|---|
| 0.5 pN | AFINES claim_a (ac/ current), ff/hand_kmc NMIIA | placeholder |
| 0.5–1 pN | KB single-molecule IIA candidate | in vitro |
| 2 pN | KB-3.18 (ff/myosin_linear `F_HEAD_PN`) | per-head |
| ≈4 pN | AFINES/Cytosim (Freedman 2017 / Tam 2021) | ensemble-fit |
| 5.0–5.7 pN | Erdmann-Schwarz PCM | cross-check |
| 2 pN "per MOTOR" | KB `PARAM-F_stall_motor` (Bangasser 2013) | motor-clutch abstraction, NOT per-head |

Single-molecule NMII step ~5.4 nm (Norstrom 2010, PMID 20511646, DOI 10.1074/jbc.M110.123851). **PI must pick the
per-head force**; the budget in §5 scales linearly with it.

## 5. ⭐ FORCE BUDGET — the honest arbiter (geometry-independent)

Required hoop force a great circle must carry to hold turgor: **F_hoop = γ_target·2πR = 140·47.1 = 6,597 pN**
(= cortex share of ΔP·πR²). Total myosin contractile budget = N_mf·N_side·f_head·duty:

| case | total myosin (pN) | / F_hoop | verdict |
|---|---:|---:|---|
| ac/ current (N_side 10, f 0.5, duty 0.10) | 221 | **0.03×** | 33× short — robustly deficient (= historical "~500× floor") |
| ac/ current (duty 0.30) | 663 | 0.10× | 10× short — deficient |
| physio-ish (30, 2.0, 0.30) | 7,956 | **1.21×** | sufficient in total |
| physio-hi (30, 4.0, 0.30) | 15,912 | 2.41× | sufficient |
| physio-max (30, 5.7, 0.30) | 22,675 | 3.44× | sufficient |
| physio + loaded-duty (30, 4.0, 0.60) | 31,824 | 4.82× | comfortably sufficient |

**Reading:** current placeholder content is force-deficient by 10–33× regardless of geometry (this IS the floor
the project keeps hitting). Physiological content puts the total budget in the sufficient band. Whether the
sufficient force *converts* to the γ=140 hoop prestress depends on native actin-network transmission (crosslink
connectivity + tension accumulation along filaments vs. the 23c §3 "sparse crosslinks / force-free ERM transmit
little" concern). The coarse oracle can't resolve this (its short-edge dipoles under-predict; it disagrees with
this total-force anchor by ~26×). **The native gbook run measuring γ_network is the only arbiter.**

## 6. ⭐ Deeper fidelity question — is the duty-ratio load-dependence the WRONG SIGN?

According to PubMed:
- **Kovács 2003** (PMID 12847096, DOI 10.1074/jbc.M305453200): NMIIA "spends only a small fraction of its ATPase
  cycle time in strongly actin-bound states" → low unloaded duty (~0.1); NMIIB higher (slow ADP release).
- **Kovács 2007** (PMID 17548820, DOI 10.1073/pnas.0701181104): "**resistive load increases duty ratio to favor
  tension maintenance by two-headed attachment**"; resistive load slows ADP release 5-fold (2A) / 12-fold (2B).

So in a *tensioned* resting cortex the physiological duty is materially **above** the unloaded 0.1 — good for the
force budget (§5). BUT the code (`bell_kinetics_analytic.py:12`) models NMII head–actin as a **pure Bell SLIP**
bond: off-rate RISES with load → duty FALLS. That is the **opposite sign** of the Kovács-2007 tension-maintenance
behavior. Consequences, to be examined (not asserted):
- The emergent duty from pure-slip kinetics would *shed* heads under cortical load → under-tension → plausibly a
  root cause of the recurring ~500× γ-floor.
- The seeded resting_setpoint mechanism is then a **workaround** for a kinetic model that won't *maintain* the
  bound state on its own. The fine-grained-faithful fix (CLAUDE.md principle) would be to add the load-dependent
  ADP-release (catch-like tension maintenance, Kovács 2007) so the resting bound fraction EMERGES, rather than
  seeding it.
- Caveat: the Hill force–velocity term already slows stepping under load (partial tension-maintenance); whether
  the *net* emergent duty rises or falls under load needs the coupled Hill+Bell+attach analysis. **Surface as a
  fidelity question, do not act unilaterally.**

## 7. Emergent-vs-imposed, and the fraction/force coupling

The code's own design says the bound fraction should be EMERGENT (`φ_b(0)=k_on/(k_on+k_off0)`, a cross-check vs
Kovács ~0.1) — yet `50/(50+0.35)=0.993`, a 10× conflict, because the Bell mechanical off-rate (0.35/s) is NOT the
ATPase-cycle detachment rate that sets duty. Moreover fraction and per-head force are **not independent**:
Kovács 2007 couples them (higher resistive per-head load → slower ADP release → higher duty). Imposing them as two
independent PI-GAPs is a simplification; the self-consistent setpoint lies on the load-dependence curve. This
reinforces §6: the faithful path is emergent, kinetics-derived.

## 8. Actionable path to actually close the gate (in order)

1. **PI ratifies the physiological setpoint SET** (§9) — not two numbers: heads-per-minifilament, per-head force,
   duty (loaded), density. From verified literature (this dossier + sourcing workflow).
2. **Raise N_side to the Billington-sourced ~30** (a physiological correction the physiological-baseline HARD rule
   already requires — currently a "TEST/GAP" placeholder; the lumped model already uses 60 heads). Expose
   N_side / f_head / density as source-gated CLI (like the resting-bound-myosin flags) so no magic-number commit.
3. **gbook A5000 native run** with the seeded bound-myosin mechanism at physiological content → measure γ_network
   (method-of-planes) + membrane residual. `PYTHONPATH=~/ffn_ac_native:~/ffn_ac_native/ffn_sim`.
4. If residual < 0.21 at physiological content → **PI ratifies gate closure**. If γ_network < 140 despite a
   sufficient total budget → transmission/connectivity (crosslink density, ERM, or the §6 kinetics) is the deeper
   blocker.
5. Consider the §6 fix (load-dependent ADP-release) so the resting bound state is EMERGENT, not seeded.

## 9. PI decision list (reframed — the force-magnitude/density SET, not two numbers)

| # | decision | current | sourced options | note |
|---|---|---|---|---|
| A | heads per minifilament (N_side) | 10 (20 heads) | ~30 (~60 heads), Billington 2013 | 3× under; physiological correction |
| B | per-head force f_head | 0.5 pN | 0.5–4 pN (single-head weak; ~3.8 pN two-headed strain) | §12; PI: keep 0.5 (ensemble does the work) or raise to ~1–4 |
| C | resting duty ratio | unset (PI-GAP) | NM2A 0.05 / NM2B 0.20–0.25 unloaded → loaded ~0.2 / 0.75–0.85 (Kovács 2007) | §12; use the *loaded* value |
| C′ | cortex NM2A/NM2B isoform mix | not modeled | MCF7 expression/localization | NEW GAP — sets whether duty ~0.05 or ~0.25 |
| D | minifilament density | 0.625/µm² (Nie HeLa) | UNSOURCED for MCF7 (~0.1–10/µm²) | §12; ⭐ **dominant magnitude driver, historical γ-floor root — do NOT back-solve from 140** |
| E | ΔP / cortical-tension baseline (MCF7) | 40 Pa HeLa proxy → 140 pN/µm | HeLa-interphase proxy (0.2 mN/m) defensible; MCF7-specific tension not cleanly sourceable | §3a (calibrated): the workflow's 10–70× "conflict" is probably a modulus-vs-tension confusion (MCF7 AFM = modulus, Omidvar 2014). Target likely stands. PI: supply a DOI-verified MCF7 resting *tension* (confirm it's a tension, not a modulus) or accept the HeLa proxy; don't back-solve ΔP from 140 |
| F | §6/§7 fidelity: emergent load-dependent duty; k_off0 mis-sourcing | pure Bell slip + k_off0=0.35/s (= ADP-release rate misused) | impose t0 duty (§12) or implement full chemomechanical cycle | fine-grained-faithful fix |

## 10. Exact native command (for PI to run once values are ratified; NOT run this session)

```
# on gbook A5000 (serial, one GPU):
PYTHONPATH=~/ffn_ac_native:~/ffn_ac_native/ffn_sim \
  ~/miniconda3/envs/ffn_sim/bin/python -m aleph.components.incumbent.driver --from-resting \
  --n-filaments 70686 --membrane-subdivisions 6 \
  --inner-solver erm_gauss_seidel_tournament --n-inner 1600 \
  --resting-bound-myosin-fraction <PI: loaded duty> \
  --resting-bound-myosin-force <PI: per-head pN> \
  --resting-bound-myosin-source "<citation>" \
  --json outputs/ac/implicit/resting_bound_native.json
# NOTE: N_side (heads/mf) and f_head are NOT yet CLI-exposed — physiological head count needs the code change in §8.2.
```

## 11. Reproduce (dev Mac, CPU)

```
python -c "from aleph.components.incumbent.resting_balance_oracle import coupled_resting_balance as c; \
  print(c(dipole_fraction=1,per_dipole_tension_pn=88,subdivisions=3).membrane_residual_mean_pn)"
# force budget: F_hoop = 140*2*pi*7.5 = 6597 pN; total = 442*N_side*f_head*duty (see §5 table).
```

## 12. Sourcing-workflow verified values (folded in — `wfruvf0se`, 74 agents, 0 errors, adversarially verified)

A citation-integrity-gated literature workflow (6 dimensions → per-claim adversarial verify → synthesis)
confirmed/refined the setpoint values. **CONFIRMED, high confidence (the duty-ratio set):**

| quantity | verified value | source (DOI) |
|---|---|---|
| NM2A (MYH9) unloaded duty | **0.05** (→0.11 at 100 µM ADP) | Kovács 2003 (10.1074/jbc.M305453200); single-mol 0.05 Melli 2018 |
| NM2B (MYH10) unloaded duty | **0.20–0.25** (transient 0.23–0.40) | Nagy 2013 (10.1074/jbc.M112.424671); Wang 2003 (10.1074/jbc.M302510200); Rosenfeld 2003 |
| NM2C (MYH14) duty | 0.02–0.26 (splice/Mg-dependent) | Heissler & Manstein 2011 (10.1074/jbc.M110.212290) |
| loaded (isometric) duty | **derived** ~0.2 (NM2A) / **0.75–0.85** (NM2B) | Kovács 2007 (10.1073/pnas.0701181104): ADP release ×5 (2A)/×12 (2B) slower under resistive load |
| per-head tension-maintenance load | **~3–4 pN** (derived ln(12)·kT/d₀, d₀=2.7 nm) | Kovács 2007 + Veigel 2003 d₀ (⚠️ d₀ from SMOOTH muscle — borrowed) |
| heads per minifilament | ~56–60 (28–30 molecules) | Billington 2013 (DOI to confirm) — NOT in verified set, domain knowledge |

**Refinements to this dossier:**
- **Duty** (decision C): the cortex-relevant value depends on the NM2A/NM2B **isoform mix** (a new PI-GAP). Pure NM2A
  ≈0.05 unloaded; NM2B ≈0.2–0.25; mixed ~0.1. Under the resting cortex's resistive load these rise substantially
  (my §6 point, now quantified). Use the **loaded** value.
- **Per-head force** (decision B): preset 0.5 pN is ~7–8× below the ~3.8 pN derived two-headed strain. Two honest
  readings: (1) single NMII heads genuinely ARE weak and the cortex reaches tension by ENSEMBLE (heads × duty ×
  density) with each head sub-stall; or (2) 0.5 pN is under-set and should be ~1–4 pN. Literature does not cleanly
  give a single-head NMII stall → PI gap. The 3.8 pN is a two-headed cooperative strain, not a single-head stall.
- **Cortical minifilament density** (decision D): the workflow could **NOT independently source it** — "the
  parameter most likely to drive the cortical tension magnitude, and it CANNOT currently be sourced" (only
  order-of-magnitude ~0.1–10/µm²). The code's 0.625/µm² (Nie 2015 HeLa) is the sole anchor. ⭐ This is the
  **dominant magnitude driver** and the historical γ-floor root — and it **must NOT be back-solved from 140 pN/µm**.
- **k_off0 mis-sourcing (new fidelity flag)**: the preset `k_off0=0.35/s` used as a mechanical Bell slip rate is
  almost exactly the NM2B **ADP-release** rate (Wang 2003 0.35/s; Kovács 2007 0.31/s). The code likely inherited an
  ATPase-cycle rate and is using it as a mechanical-unbinding rate — a category error that compounds §6/§7.
- **Emergent-vs-imposed** (decision F): resolved — impose the t0 duty from ATPase-cycle literature (0.05 A /
  0.2–0.25 B / ~0.1 mixed), because the model carries only a mechanical Bell rate and **cannot self-consistently
  derive the biochemical duty** unless the full chemomechanical cycle (Pi 0.13/s rate-limiting, ADP load-dependent,
  ATP-detach) is implemented. The `φ_b(0)=0.993` assertion is a modeling-category error, not a data error.

**Citation-integrity honesty flags (from the adversarial pass, per the 3-hallucination rule):**
- Melli 2018 (10.7554/eLife.32871) characterizes NM2A/2B, **not** NM2C (task mislabel corrected); true NM2C source
  is Heissler & Manstein 2011.
- Isometric duties (0.75–0.85 / 0.2) and the 3.8 pN load are **DERIVED, not measured** (is_primary=false).
- d₀≈2.7 nm / κ≈0.45 pN/nm are **smooth-muscle** (Veigel 2003) — borrowed, must not be attributed to NMIIA/IIB.
- Heads-per-minifilament, cortical areal density, MCF7 cortical tension, MCF7 ΔP are project-memory/domain
  knowledge, **not** the adversarially-CONFIRMED set — low/medium confidence, routed to PI. No DOIs fabricated.
- **NO GATE-TUNING:** every value is from literature; nothing was chosen to hit 140 pN/µm.

**Net effect on the gate status:** unchanged (OPEN). The verified duty/force values (loaded duty 0.2–0.85,
per-head ~1–4 pN) sit within the §5 physiological band, so the force-budget "sufficient in total" conclusion holds
against the (defensible) 140 pN/µm target. The workflow's decision-E alarm (10–70× MCF7 conflict) is walked back
in §3a — it is probably a modulus-vs-tension confusion, so the target likely stands. The genuinely open items are
therefore: (1) the **cortical minifilament density** (decision D — the dominant magnitude driver, only a HeLa
proxy), (2) the **per-head force / N_side / isoform-mix** physiological set (decisions A–C, C′), and (3) the
**conversion** question — only the native gbook run measures whether the sufficient force becomes γ=140 prestress.

## 13. Native A5000 runs (2026-07-23) — the deepest finding: the seeded myosin is INERT at native

Two native `--from-resting` runs at 70,686 filaments + subdiv-6 membrane, TEST setpoint (duty 0.3, f_head 2 pN),
`erm_gauss_seidel_tournament`, n_inner 400. Artifacts: `outputs/ac/resting_native/`.

| run | capture | realized duty | residual_start (pN) | candidate | outer | note |
|---|---|---|---|---|---|---|
| #1 | 0.05 µm (production) | **0.170** (1506/8840; shortfall) | 2634.1191655118346 | 1483.1048741614643 | rolled back | capture cap |
| #2 | 0.6 µm (test) | 0.300 (8840 eligible) | 2634.1191655118346 | 1483.104875297183 | rolled back | full-seed |

Three native ground-truth findings:

1. **The capture cap is purely the 0.05 µm radius** (Codex finding 7): at 0.6 µm all 8,840 heads are eligible and
   duty 0.3 seeds fully. Heads sit ~0.2 µm off the cortex actin (the code's own CUDA test uses 0.6 µm). So
   production seeding is capped at 17% — a head-placement / capture geometry issue.
2. **The binding blocker is the inner convergence + step acceptance.** `inner_converged=False`; the inner solve
   finds a candidate that lowers the force residual 2634→1483 pN, but it is REJECTED (`outer_rolled_back=True`,
   `residual_end=residual_start`). So the cortex never actually deforms/contracts. This is the same unclosed
   resting solve of diagnoses 22c–23d (n_inner=400 here; the diagnoses used more iterations and still landed a
   normalized projected residual 0.64–0.71 vs the 0.21 gate).
3. **The seeded myosin cannot even be exercised.** `residual_start` is IDENTICAL to 16 digits and the candidate
   identical to ~1e-6 pN between 17% and 30% seeding. Two non-exclusive reasons: (i) at the *undeformed* seeded
   t0 the per-minifilament tension is internally balanced (Newton's 3rd law) so it does not perturb the start
   residual, and (ii) the step rolls back before the cortex can contract, so the ERM never stretches and the
   myosin's transmission is never engaged. ⟹ **whether the seeded mechanism transmits at native is UNTESTABLE
   until the solve accepts a step** — it is not disproven, but it is also not demonstrated at native. (Whether
   `residual_start` should capture the undeformed-seeded myosin at all is a code question worth checking.)

**Reframed conclusion (native-grounded):** the resting_setpoint "mechanism demonstrated" (`9bb91e0a`) is a
CPU-oracle result only (the oracle relaxes the cortex in a 20k-iter breathing mode so the ERM can stretch). The
native single-step solve does NOT reproduce that — it rolls back before the cortex contracts. So the resting
gate's binding blocker is **(a) inner convergence + step acceptance**, then **(b) cortex→ERM→membrane transmission
(untested until (a))**, and only then **(c) the myosin magnitude/target** (§Corrections). The sourcing work
(§1–§12) is necessary but not sufficient; native convergence must be fixed first. Native-first from here.
