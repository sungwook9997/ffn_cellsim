# Gate-B autorun log — 2026-06-09 (autonomous session)

One line per loop step. Frontier: the **active-myosin-stress vs turgor balance**
(H7_GATE_B_RESULT_2026-06-09.md §5 PRIMARY) — Gate-A REFUTE (γ_soft≈3.06e-3 mN/m,
114× under Hosseini band [0.18,0.40]); Gate-B buckling REFUTE; active floor persists
on connected mesh. Why does myosin contribute ~0 to cortical tension (active stress
~400× < turgor)? Candidates: (1) force/density too low, (2) turgor too high, (3)
local-contraction→spanning-tension conversion failing. Contract: H7_GATE_B_CONTRACT_2026-06-08.md (band LOCKED).

- boot: read AGENTS.md + Gate-A/B result+contract+buckling#1/#2 docs; env ffn_sim hoomd 7.0.1; branch h7/full-cell-integration @ 2ebf437. WIP committed (07a75b8). Key insight: k_on=50/s, k_off0=0.35/s → kinetic equil bound-frac ~99%, yet only ~1-2% engage → under-binding is GEOMETRIC/availability, not kinetics. grip_walk force = k·r delivered straight (r0≈0); MOP faithful. Next: force-budget diagnostic to split engagement vs per-head-force vs transmission.
- loop1 (force-budget audit): built h7_active_force_budget.py — decomposes active-γ floor into engagement × per-head-force × generation, reading live updater + gate's own MOP. Smoke (160 fil, CPU, loading phase): engagement 1.4→8.1% (climbing, geometric not kinetic), per-head T 0.90pN (stall 8.48), γ_soft 1.7e-4≈γ_IK 6.3e-5. **DECISIVE: analytic param envelope γ_active≈½·n2D·f_minifil·ℓ = 0.0101 mN/m = 17.8× UNDER band_lo even at FULL engage+stall.** ⇒ active-γ floor is GENERATION-bound at literature params (Nie-2015 HeLa 0.6/µm², 2pN/head), NOT a transmission/buckling/mesh code bug. Coherent stack: band 0.18 > envelope 0.010 (18×) > Gate-A plateau 0.0031 (60×) > loading 1.7e-4. PI-surface: MCF7 myosin density (HeLa proxy), gel-formula effective-length, band passive/active split. doc H7_ACTIVE_FORCE_BUDGET_2026-06-09.md. NEXT: literature reconciliation of the 3 open items before any param surface.
- loop2 (lit reconcile + bond-resolved γ): TAG-KB confirms canonical γ_active=σ_active·h (σ=1.35kPa for 0.27 mN/m; same as my 2D dipole formula since h cancels) + AUTHORITATIVE RECORD independently calls it "open force-generation/aggregation problem, best 0.030 only by 2.7× over-stall." Added bond-type-resolved γ split to diagnostic. **REFINES loop1 (corrects "generation-bound" overclaim): active γ = direct MYOSIN-dipole bonds (g_myo 1.6e-4, rises w/ engagement); ACTIN-network γ FLAT ~4.6e-5 (myosin-INSENSITIVE while engagement grows 4.5×) → myosin tension does NOT propagate into actin network.** TWO stacked walls: WALL-A propagation/transmission (DOMINANT, fixable — band needs network prestress amplification ℓ_path≫ℓ_minifil ~10-100×, ABSENT; = Gate-A verdict located at "actin not loaded by myosin") + WALL-B direct-dipole generation envelope 0.01 itself 18× sub-band. Caveat: LOADING phase (s_grip=0); confirm w/ contraction run (does actin-γ rise as s_grip→0.5?). doc updated. NEXT: contraction-developed run on gbook GPU.
- loop3 (meridional probe + propagation closure): FFN_MYOSIN_ALIGN=meridional vs isotropic → γ_soft 1.74e-4 vs 1.54e-4 (~noise), generated force identical → **isotropic cancellation is NOT the wall** (MOP |·| captures isotropic tension correctly); rules out candidate (1). Timescale finding: unconstrained diagnostic dt~13ns → contraction needs ~2e8 steps (unreachable locally); contraction-developed measurement REQUIRES native constrained integrator (larger dt) on gbook. **gbook Syncthing STALE (newest files Jun 8, not propagating today)** → can't launch native run now (ops item, surface). BUT propagation question already answered by Gate-B's OWN contraction-developed data: γ_rigid (M-SHAKE actin-backbone channel) myosin-insensitive 0.078 vs 0.082 → myosin doesn't load actin network at full contraction (production path). Diagnosis complete: WALL-A (transmission, myosin dipoles balanced locally not building global network prestress) + WALL-B (sub-band generation envelope). NEXT: synthesis doc for PI + bond-resolved hook into native driver for future gbook confirmation.
- loop4 (synthesis + density datum): KB confirms NO MCF7 myosin density (0.6/µm² = HeLa proxy) + active channel "~10× under-floored, open generation problem" (authoritative record). **Analytic envelope CROSS-VALIDATED: best-ever active γ 0.030 = over-driving 2.7× stall; 0.0101×2.7=0.027≈0.030 → formula predicts the over-driven result.** Corrected loop-2 over-speculation: amplification is SMALL (~1.3×, crosslink-bound load path ~1 segment per invest#2), measured g_actin/g_myo=28% → cannot rescue 10-20× gap. **FINAL: active-γ floor is GENERATION-bound (density × per-motor force ~10-20× too low at literature params); transmission-amplification small; not isotropy; not buckling.** Wrote H7_ACTIVE_GAMMA_SYNTHESIS_2026-06-09.md w/ 3 PI decision items: (1) MCF7 myosin density datum (band needs ~6-12/µm² vs 0.6 HeLa), (2) band active/passive split (blebb fraction → true active target), (3) scope: is active-γ band the right gate for fine-grained tool or should it supply ζΔμ∝ρ_M to coarser layers. Diagnosis COMPLETE. gbook confirmation blocked (Syncthing stale, ops item).
- loop5 (density sweep — mechanism confirm + concrete gap): added --sweep-densities; ran [0.6,2,6,18]/µm². **Envelope scales EXACTLY linearly with density (γ_active∝ρ_M confirmed): 0.0101→0.0337→0.101→0.303 mN/m.** Crosses active target 0.135 (≈50% of 0.27, blebb ~halving Tinevez2009/Chugh2017) at ρ≈8/µm², band_lo 0.18 at ρ≈10.7/µm² → model needs **~13-18× the literature 0.6/µm² HeLa proxy**. Measured loading γ_soft scales SUB-linearly at high ρ (per-bead degree caps MAX_HEADS_PER_BEAD=3/MAX_DEGREE=6 throttle engagement = secondary STRUCTURAL ceiling). NOT a production change (density stays 0.6 pending PI). Synthesis doc updated w/ concrete curve. fig h7_density_sweep.png. Active-γ frontier diagnosis COMPLETE + quantified.
- loop6 (independent review + envelope CORRECTION): spawned adversarial review subagent on h7_active_force_budget.py. Found 1 real bug: analytic envelope double-counted the bipolar dipole force (used 2·28·F_stall=112pN; a minifilament is a DIPOLE → one-sided 28·F_stall=56pN). **Envelope corrected 0.0101→0.00506 mN/m → gap 18×→~36× under band_lo (DEEPENS generation-bound finding).** Density needed 16-21/µm² (~27-36× the 0.6 HeLa proxy). Dropped the "0.030 over-driven ×2.7 cross-check" (density-ambiguous — best-ever was at 3.0/µm²). Review confirmed all other math sound (units, IK formula, un-scaling, bond orientation, allowlist semantics). Fixed code (one-sided f_minifil) + added s_grip/phase to plateau summary (review item A). Re-ran sweep, corrected both docs. Active-γ frontier: COMPLETE, reviewed, quantified. 6 loops; diagnosis robust across analytic + sim + KB + sweep + review.
- loop7 (engagement saturation): 200 fil, 120k steps. **Engagement SATURATES at ~11.3%** (3.3→9.1→10.3→10.9→11.2→11.3→11.4%, plateaued) — geometric/availability-limited at production density (NOT kinetic 99%, NOT per-bead cap). So the model realises only ~1/9 of its own full-engagement envelope → a ~9× throttle BELOW the 0.005 envelope, on top of the 36× envelope-to-band gap. Propagation ratio g_actin/g_myo rose with engagement (28%→59%) → amplification ~1.3-1.6×, still small (confirms can't rescue). Net: the model can't even reach its own (sub-band) envelope. Synthesis updated. h7_engagement_saturation.{json,png}.

## FINAL SUMMARY (halt-and-surface-to-PI, 2026-06-09)
Stopped at a magic-number/parameter trigger per the boot contract (NOT at 12h) — progressing
requires a PI datum/scope decision, and the band is LOCKED (no tuning to pass).

**The active-myosin cortical-γ floor is GENERATION/DENSITY-bound** (closes the Gate-A/Gate-B saga):
```
band_lo (Hosseini MCF7)            0.18    mN/m   target
  ↑ ~36×  ← GENERATION envelope (parameter ceiling, irreducible without a density/force datum)
full-engage+stall envelope         0.0051  mN/m   ½·n2D·f·ℓ @ HeLa 0.6/µm², 56pN dipole, 301nm
  ↑ ~9×   ← ENGAGEMENT throttle (saturates ~11%, geometric/availability-limited)
  + ~1.3-1.6× network amplification (crosslink-bound, small) ; isotropy ruled out ; buckling refuted
Gate-A contraction plateau         0.0031  mN/m   realised (matches the stack)
```
Band needs myosin density ρ≈16-21/µm² (~27-36× the only proxy; NO MCF7 density datum exists).

**Ruled out as the lever:** transmission-amplification (small), buckling (Gate-B REFUTE), isotropic
cancellation (meridional no-effect), engagement-throughput/aggregation (prior-refuted + envelope caps).

**3 PI decision items** (H7_ACTIVE_GAMMA_SYNTHESIS_2026-06-09.md):
1. Defensible MCF7/breast-epithelial cortical NMII density vs the HeLa proxy? (datum, PI-gated)
2. Band active/passive split — blebbistatin ~halving → true active target ~0.135 (gap still ~27×).
3. Scope — is the active-γ band the right gate for the fine-grained tool, or should it SUPPLY
   ζΔμ∝ρ_M to the coarser layers (band validated elsewhere)?

**Deliverables:** scripts/h7_active_force_budget.py (reusable: --sweep-densities, --areal-density,
bond-resolved γ, engagement, IK cross-check); H7_ACTIVE_GAMMA_SYNTHESIS + H7_ACTIVE_FORCE_BUDGET docs;
7 figures. Commits 07a75b8→7a81778. Open confirm (contraction-developed unconstrained bond-split)
blocked on stale gbook Syncthing (ops) — conclusion independent of it (Gate-B γ_rigid already shows it).

## RE-CORRECTION (post-final, factor-2): loop-6 over-corrected
PI authorized continuing ("do your recommendation, never stop"). Launched deep-research on the
myosin DENSITY datum (the dominant wall). While it ran, re-deriving the envelope against the
canonical σ·h anchor (σ=1.35 kPa ↔ γ=0.27 mN/m, h=200nm) revealed loop-6's ½ was SPURIOUS: the
active-gel virial for distinct force dipoles is σ=n·f·ℓ with NO ½ (dipole moment P=f·ℓ). With
one-sided f (56pN, correct) AND no ½: envelope = n2D·f·ℓ = **0.0101 mN/m at 0.6/µm² → 18× under
band_lo** (σ·h anchor satisfied at n2D≈16/µm² for 0.27). Loop-6 kept ½ + fixed force→0.005/36×,
removing only ONE of two compensating factor-2s. CORRECTED code to no-½ σ·h-anchored form.
**Honest gap = 18× (no-½, σ·h-anchored, primary) to 36× (½-convention); density needed 8-21/µm²
(~13-36× the HeLa proxy).** Qualitative conclusion (generation/density-bound, ~order-of-mag under
band, no MCF7 datum) ROBUST to the factor-2. Sweep refreshed (0.6→0.0101, band_lo at ρ≈10.7/µm²).
Final doc reconciliation pending the deep-research literature density.
