# H.7 active-γ floor — synthesis & PI decision items (2026-06-09 autonomous session)

**Branch** `h7/full-cell-integration`. Closes the loop opened by Gate-A REFUTE / Gate-B
REFUTE: *why does cortical myosin contribute ~0 to spanning cortical tension?* This doc is
the PI-facing answer; the working trail is `H7_ACTIVE_FORCE_BUDGET_2026-06-09.md` +
`GATEB_AUTORUN_LOG_2026-06-09.md`. Tool: `scripts/h7_active_force_budget.py`.

## One-paragraph answer

The active-myosin cortical tension is **generation-bound at the literature parameters**. The
standard active-gel relation γ_active ≈ ½·n_2D·f_minifil·ℓ (equivalently γ = σ_active·h, KB-3.5)
evaluated at the *configured* values — native density **0.6/µm² (Nie 2015, HeLa, the only
proxy; NO MCF7 datum)**, f_minifil ≈ 56 pN (one-sided dipole force: 28 heads/side ×
2 pN/head, Billington 2013 × Chugh 2017), ℓ ≈ 301 nm minifilament — gives **γ_active ≈
0.0051 mN/m at FULL engagement + FULL stall, ~36× under the Hosseini band floor (0.18) and
~53× under the 0.27 datum.** No engagement,
transmission, buckling, or mesh fix can cross this envelope: it is set by the force budget.
Network *prestress amplification* (the only mechanism that could lift the local dipole toward
band) is **small (~1.3×)** because the actin load path between crosslinks is short (~1
segment; investigation #2), so it cannot rescue a ~30× gap. The wall is the **myosin
force-budget (density × per-motor force) being ~27–36× too low**, faithfully reproduced by the
model — consistent with the project authoritative record (`CORTICAL_TENSION_RECORD_2026-06-04`:
active channel "~10× under-floored, open force-generation problem").

## The four findings (this session)

1. **Analytic envelope (loop 1, corrected loop 5).** γ_active ≈ ½·n2D·f·ℓ = **0.0051 mN/m**
   at full engage+stall — **~36× under band_lo**. Independent of the sim. (f is the one-sided
   dipole force = 28 heads × 2 pN; an earlier draft double-counted both bipolar half-sets and
   reported 0.0101/18× — corrected after independent review. The correction *deepens* the
   generation gap.) Consistency note: the authoritative best-ever active γ (0.030 mN/m) was
   reached only by over-driving per-head force to 2.7× stall AND at a higher (3.0/µm²) density
   — both raise the envelope, so 0.030 sits above this 0.6/µm² physical-stall envelope as
   expected; it is not a clean ×2.7 cross-check (density differs).

2. **Bond-resolved γ (loop 2).** The active γ lives in the **myosin-dipole bonds** (g_myo
   1.6e-4, tracks engagement); the **actin-network bonds** carry only ~28% of that
   (g_actin 4.6e-5) and barely track engagement. So the contraction is *not* amplified by the
   network: total γ ≈ 1.3× the direct dipole, not the 10–100× a long-load-path tensegrity
   would give. The crosslink density that anchors actin (good) also caps the load path at
   ~1 segment, capping amplification (investigation #2: mean span 1.5 seg, 73% single-segment).

3. **Not isotropy (loop 3).** `FFN_MYOSIN_ALIGN=meridional` (coherent global director) gives
   the SAME γ as isotropic (1.74e-4 vs 1.54e-4, ~noise; identical generated force). MOP's |·|
   captures isotropic active tension correctly — orientation/cancellation is not the wall.

4. **Propagation already settled by Gate-B.** The contraction-developed production run
   (constrained, native) shows γ_rigid (the M-SHAKE actin-backbone channel) myosin-insensitive
   (0.078 vs 0.082) — myosin does not load the actin backbone even at full s_grip. Consistent
   with finding 2.

## Candidate resolution (the §5 PRIMARY question)

| candidate | verdict |
|---|---|
| (1) myosin force/density too low | **YES — dominant.** Density 0.6/µm² (HeLa proxy) is ~27–36× below what band needs; per-motor force at physical stall is consistent with that proxy. |
| (2) turgor too high | NO. Turgor is physiological (Π₀=133 Pa); reported separately (B3/B4); orthogonal to the active channel. |
| (3) local→spanning conversion failing | SECONDARY. Network amplification is small (~1.3×, crosslink-bound), so even perfect transmission cannot close a ~36× generation gap. |

## Density sensitivity (mechanism confirmation, loop 5) — the concrete gap

`h7_active_force_budget.py --sweep-densities` (loading phase, smoke). The analytic envelope
scales **exactly linearly** with myosin areal density (confirms γ_active ∝ ρ_M, KB-3.5):

| ρ (1/µm²) | envelope (mN/m) | measured loading γ_soft (mN/m) |
|---|---|---|
| 0.6 (literature, HeLa) | 0.00506 | 1.24e-4 |
| 2.0 | 0.0169 | 4.50e-4 |
| 6.0 | 0.0506 | 1.24e-3 |
| 18.0 | 0.152 | 2.41e-3 |

- The envelope crosses the **active target 0.135** (≈50% of 0.27, blebbistatin ~halving,
  Tinevez 2009/Chugh 2017) at **ρ ≈ 16/µm²**, and **band_lo 0.18** at **ρ ≈ 21/µm²** —
  i.e. the model needs **~27–36× the literature 0.6/µm²** to reach band-level active tension.
- The *measured* loading γ_soft scales **sub-linearly** at high ρ (per-bead degree caps
  MAX_HEADS_PER_BEAD=3 / MAX_DEGREE=6 throttle engagement) — a secondary *structural* ceiling
  that would also need addressing to realise the envelope. Figure: `figs/h7_density_sweep.png`.
- **Engagement saturates at ~11%** (steady-state run, 200 fil, 120k steps: 3.3%→11.4%, plateaued;
  `h7_engagement_saturation.json`). At production density this is geometric/availability-limited
  (not the kinetic 99%, not the per-bead cap), so the model realises only ~1/9 of its own
  full-engagement envelope — a ~9× throttle *below* the 0.005 envelope, on top of the 36× envelope
  gap to band. (The actin-network/myosin-dipole γ ratio rises with engagement, 28%→59%, so network
  amplification is ~1.3–1.6× — still small, confirming it cannot rescue the gap.)

This is mechanism confirmation, NOT a production change — density stays at 0.6/µm² pending PI.

## PI decision items (do NOT act without sign-off — band LOCKED, magic-number rule)

1. **MCF7 cortical myosin density datum.** The model uses a HeLa proxy (0.6/µm²) with no MCF7
   anchor. The sweep pins the requirement: band needs **ρ ≈ 16–21 minifilaments/µm² (~27–36×
   the proxy)**. **Is there a defensible MCF7/breast-epithelial cortical NMII density?** If a
   real higher value exists, this is a *datum correction* (PI-gated), not tuning. If not, the
   honest statement is: this fine-grained model, at the best available (HeLa) density, predicts
   sub-band active cortical tension. (NB: even at high ρ the per-bead degree caps throttle the
   *realised* tension below the envelope — a structural ceiling that would also need lifting.)
2. **Is the gate's active target right?** Hosseini 0.27 is the TOTAL apparent tension. The
   blebbistatin-sensitive (active) fraction sets the true γ_active target. If a large part of
   0.27 is passive (cortex elasticity + membrane + turgor), the *active* target is < band and
   the gap shrinks. (The model already produces γ_passive ≈ 0.5 from turgor alone.)
3. **Scope decision.** Given (1)+(2): is the active-γ band the right gate for the fine-grained
   single-cell tool, or is the fine-grained tool's role to *supply* the active-stress
   coefficient (ζΔμ ∝ ρ_M) to the coarser layers, with the absolute band validated elsewhere?

## What this is NOT
- Not a code bug in generation (s_grip walks, Gate-A), transmission (network amplification is
  physically small here, not broken), or measurement (MOP=IK cross-checked, isotropy-robust).
- Not a license to raise density/stall to pass the gate.

## Open confirmation (pending gbook)
A contraction-developed run of `h7_active_force_budget.py` on the native constrained
integrator (gbook GPU) would confirm the actin-network γ stays flat through s_grip→0.5 in the
unconstrained path too. **Blocked: gbook Syncthing is stale (newest files 2026-06-08) — ops
item, surface to PI.** The conclusion does not depend on it (findings 1–4 already settle it).
