# H.7 active-γ — myosin-overlap mechanism redefinition (DESIGN, pre-execution)

**Date** 2026-06-09 · **Branch** `h7/full-cell-integration` · **Status** DESIGN
(Sanity Gate Protocol mandatory before first execution — this is that design).
**PI decision (2026-06-09):** active-γ closed → architecture rebuild → (B) faithful
mesh hits a gate wall (load-path is physiologically short) → **(iii) reframe the lever
as a MYOSIN-OVERLAP mechanism**: myosin transmits force along the actin-myosin overlap
contour, independent of crosslink spacing. This doc designs that mechanism.

## 1. Physical target (the literature the mean-field model must capture)

Chugh 2017 (Nat Cell Biol), Truong Quang 2021, Miyazaki 2015: cortical tension is set
substantially by **actin-filament length / actin-myosin OVERLAP**, *largely independent
of myosin number*. The mechanism is the sarcomere-like one: a bipolar myosin minifilament
in the overlap zone of (anti)parallel actin engages many cross-bridges along the overlap
and slides the filaments, so the **tension transmitted to the filament/bundle ends ∝
(overlap length) × (engaged cross-bridge density) × (per-bridge force)** and is borne by
the actin contour across the overlap — NOT terminated at the nearest crosslink.

This is why interphase→mitosis tension rises >3× despite LOWER myosin (Chugh): longer
filaments / more overlap, not more motors.

## 2. Current model — three compounding deficits (this session's force-budget diagnosis)

Measured at the production operating point (`h7_active_force_budget`, n=1000):
realized active γ_soft ≈ 1.5e-3 mN/m vs MCF7 active target 0.40 mN/m → ~260× under.
Decomposed:

**(D1) Per-head force is capped at k·r ≈ 0.74 pN ≪ F_stall 8.48 pN (~11× under).**
grip_walk delivers `F = k_head_actin·(r − r0_eff)`, `r0_eff = max(r − s_grip, 0)` →
`F = k_head_actin · min(s_grip, r)`. Once the grip walks past the bond length
(`s_grip ≥ r`), the force **saturates at k·r**, independent of further walking. With
`k_head_actin = 1 pN/µm` (×4.24 mesoscale = 4.24 pN/µm) and binding distance r ≈ 175 nm,
the cap is `k·r ≈ 0.74 pN`. To reach F_stall the head would need r ≈ F_stall/k ≈ 2 µm
(physically impossible). The Hill stall (`F_series = k_series·min(s,r)`,
`k_series = 1/(1/k_head_actin + 1/k_head_spring) = k/2`) is even softer — the **soft
perpendicular offset spring `k_head_spring = 1 pN/µm` sits in SERIES with the power
stroke and bottlenecks force transmission**. Real myosin cross-bridge stiffness is
~0.3–2 pN/**nm** = 300–2000 pN/µm (Veigel/Finer/Molloy) — the model's 1 pN/µm conflates
the structural offset spring with the load-bearing cross-bridge.

**(D2) No overlap accumulation — local shunting.** Each engaged head pulls its single
bound bead; the force dissipates to the nearest crosslink anchor (~1 segment; steps 1–3
of `H7_ARCHITECTURE_LOADPATH`). Heads engaged along the same filament do NOT add into a
filament tension that survives to the ends. So tension is NOT ∝ overlap — confirmed by
the bond-resolved γ (actin-network channel myosin-insensitive, Gate-A/B).

**(D3) Engagement ~11 % (geometric/availability).** Steady-state engaged fraction
plateaus ~11 % (not the kinetic 99 %), so only ~1/9 of even the capped envelope is
realized.

Net: γ ≈ envelope(18–36× under) × (1/9 engagement) × (1/11 per-head) → ~200–400× under
band — consistent with the project record ("active channel ~400× under turgor").

**The overlap mechanism must address D1 + D2** (D3 is geometric and partly downstream of
overlap geometry). D1 is the per-bridge force; D2 is the along-contour accumulation.
Both are needed for tension ∝ overlap × density × force.

## 3. The redefinition — tension generated and transmitted along the overlap

Replace "each head pulls its nearest bead with a soft k·r spring" by "a minifilament in
an actin overlap generates a contractile tension distributed along the overlap contour":

- **O2 (mechanism, primary):** when a minifilament's heads engage an (anti)parallel actin
  pair over an overlap, the per-head contractile forces ADD ALONG THE CONTOUR so the
  tension delivered to the filament backbone within the overlap is `Σ_heads F_head`. In a
  bead-spring filament this means the engaged heads on one actin filament must pull
  COHERENTLY toward the minifilament centre (already true for a bipolar dipole) AND the
  actin backbone within the overlap must carry the cumulative tension (so g_actin rises
  with overlap). This is the faithful sarcomere mechanism; it makes tension ∝ overlap
  length (number of engaged heads on the filament) rather than ∝ a single local bond.
- **O1 (per-bridge force, supporting):** the per-head force must be able to reach F_stall
  at physical sub-bead stretches — i.e. the load-bearing transmission stiffness must be
  the cross-bridge stiffness (~0.3–2 pN/nm), not the 1 pN/µm structural offset spring.
  This is a PARAMETER/structural-model correction (see §5 — magic-number/PI gate).

## 4. Implementation plan (faithful, incremental, sanity-gated)

Stage M1 — **diagnose & instrument** (cheap, no model change): add an overlap-resolved
readout to `h7_active_force_budget`: per minifilament, the number of engaged heads on
each bound actin filament and the cumulative tension delivered along that filament's
contour vs the single-bond tension. Confirms D2 quantitatively (tension does NOT scale
with engaged-heads-per-filament today). [in-lane: scripts/h7_, cortex/cortical_tension]

Stage M2 — **O2 overlap accumulation** in `cortex/myosin.py` grip_walk: ensure the
engaged heads of a minifilament on a given actin filament deliver forces that sum along
the filament contour (coherent minus-end-ward pull), and that the actin backbone bonds
within the overlap carry the cumulative tension. Verify g_actin rises ∝ engaged-heads.
[in-lane: cortex/myosin.py — opt-in `stepping_mode` extension, legacy byte-identical]

Stage M3 — **O1 cross-bridge stiffness** (ONLY if M2 insufficient and PI-approved): set
the head-actin transmission stiffness to the literature cross-bridge value, decoupled
from the structural offset spring. Re-check CFL (stiffer k → smaller τ; the current note
says k=1pN/µm is far above CFL — a 1000× stiffer bond may approach the membrane/ERM CFL).

Stage M4 — **validate**: contraction run (GPU, gbook) measuring active γ on the overlap
mechanism vs band; per the loop, this needs the deploy (PI remote-overwrite) — bank the
cheap M1/M2 CPU evidence first.

## 5. Sanity gates (CLAUDE.md Hard Rule — record before execution)

1. **Dimensional:** `Σ_heads F_head` [N] borne by backbone [N]; tension/(2πR) → γ [N/m]. ✓
2. **Boundary:** zero overlap (no engaged heads on a filament) → zero added filament
   tension; full overlap (all heads, full stall) → the n2D·f·ℓ envelope, NOT above it.
3. **Conservation:** Newton's third law — every added contractile force has its reaction
   on the minifilament backbone (no net momentum injection); total cortex force = 0 at
   equilibrium.
4. **Sign:** the accumulated overlap tension is CONTRACTILE (raises γ), never expandile;
   the bipolar antiparallel gate (`_bipolar_accepts`) must keep the two sides coherent.
5. **Measurement consistency:** g_actin (network channel) must RISE with engaged-heads-
   per-filament — the direct test that overlap now transmits (today it does not).
6. **CFL/precision:** if O1 stiffens k, re-derive τ = γ_b/k and confirm dt is safe;
   surface to PI if it forces an integrator-freeze change.
7. **No gate-chasing:** the band stays LOCKED; success = the MECHANISM produces
   tension ∝ overlap with literature per-bridge force, NOT a γ tuned to the band.

## 6. Open items → PI / literature (surfaced, not actioned)

- **(D1/O1) cross-bridge stiffness datum.** Is `k_head_actin = 1 pN/µm` (= the structural
  offset `k_head_spring`) correct for the LOAD-BEARING head-actin cross-bridge, or should
  it be the ~0.3–2 pN/nm cross-bridge stiffness? This is the single largest per-head force
  deficit (~11×). Changing it is a **magic-number/parameter trigger → PI sign-off +
  literature anchor** (the brief literal specified 1 pN/µm for the offset spring; the
  head-actin bond inheriting it may be the modeling conflation). DO NOT change inline.
- **(O2) is overlap-accumulation already physical or a new term?** Determine in M1/M2
  whether coherent multi-head pull on one filament ALREADY accumulates (and is just
  dissipated by crosslinks) or whether the grip_walk geometry actively prevents it.
- **(M4) validation needs a GPU contraction run** = gbook deploy (PI remote-overwrite).
- The honest fallback if O2+O1 still fall short: the active-γ floor is a genuine
  fine-grained-model prediction (sub-band) and the tool's role is to SUPPLY ζΔμ to the
  coarse layers, not to hit the absolute band (the long-standing scope question).

## 7. ⭐ ROOT-CAUSE FOUND (D1/O1): cortex head stiffness is 1000× softer than the canonical motor

The H.3 brief states the Stam-Hocky/AFINES motor kernel is SHARED via
`ffn_sim/bridge/motor.py` (H.4-owned). The two configs that instantiate that same motor
DISAGREE by exactly 1000×:

| param | bridge/motor (H.4, phase1_h4.yaml) | cortex/myosin (H.3, phase1_h3.yaml) | ratio |
|---|---|---|---|
| head spring k | `head_spring_k: 1.0e-3` N/m = **1 pN/nm** ("AFINES motor stiffness") | `k_head_spring: 1.0e-6` N/m = 1 pN/µm | **1000×** |
| head-actin k | (k_xb = head_spring_k = 1e-3) | `k_head_actin: 1.0e-6` ("= k_head_spring") | **1000×** |

The literature cross-bridge stiffness is ~0.3–2 pN/**nm** (Veigel 2002, Kaya-Higuchi 2010,
Finer 1994) = 300–2000 pN/µm. The canonical motor's 1e-3 N/m (1 pN/nm) is right; the
cortex 1e-6 N/m (1 pN/µm) is **1000× too soft — a pN/µm-vs-pN/nm unit slip in the H.3
brief literal**, inherited by `k_head_actin` ("= k_head_spring").

**Consequence (this IS deficit D1).** With k = 1e-6, the Hill stall needs
min(s,r) = F_stall/k_series ≈ 17 µm — unreachable, so the head walks to the 2ℓ₀ cap
delivering only k·r ≈ 0.74 pN. With the canonical k = 1e-3, the Hill stall binds at
min(s,r) = F_stall/k_series ≈ **4 nm**, so the head delivers **F_stall** — ~11× more force,
the dominant fixable deficit, AND it makes the grip_walk physically correct (the motor
stalls at a physiological sub-bead stretch instead of saturating on a soft spring).

**This is a brief-literal CONTRACT value** ("k_head_spring = 1 pN/μm = 1e-6 N/m brief
literal", phase1_h3.yaml:331). Per the hard rules, correcting it is a contract change →
**PI sign-off + literature anchor**, NOT an inline edit. Two coupled checks:
- **CFL (integrator-freeze adjacent):** at k = 1e-3, k_backbone = 10·k = 1e-2 N/m →
  τ_backbone = γ_b/k_backbone ≈ 3.91e-10/1e-2 ≈ 39 ns ≈ 3×dt(13 ns). The H.3 CFL note
  ("FAR ABOVE dt_CFL") no longer holds at the corrected stiffness — the step may need to
  shrink (cortex dt is set via cfl_safety_factor in the cortex config, NOT the
  PI-frozen integrator/, per this session's pattern), or the backbone treated rigid
  (M-SHAKE constrained mode already does this). MUST re-derive before running.
- **Mesoscale interaction:** the ×4.24 parallel-scaling multiplies k too; s_grip_max =
  F_stall/k stays invariant, so the correction composes cleanly with mesoscale.

This supersedes the framing that O1 is merely "a soft spring" — it is a concrete 1000×
divergence from the project's own canonical motor module. HALT → PI for the contract
sign-off (brief-literal change + CFL re-derivation).

## 8. ⚠️ EMPIRICAL: the stiffness fix is COUPLED to the grip_walk r0 convention (2026-06-09)

PI authorized the cross-bridge stiffness unit-slip correction (k_head_spring/k_head_actin
1e-6→1e-3). Applying it to the config ALONE and running h7_active_force_budget (n=1000, mesoscale)
revealed the coupling the design (§3 O1+O2) anticipated:

- **Per-head bond tension EXPLODED to ~322 pN** (vs F_stall 8.48 pN, ~38×) and g_myo (myosin
  dipole channel) jumped 135× (1.06e-4 → 1.43e-2), γ_soft 1.6e-3 → 1.3e-2 (114×→14× under band).
- ROOT: the grip_walk attach bond uses **r0 ≈ 0**, so a freshly-bound head at binding distance
  r (~76 nm, up to max_bind_dist 60 nm + offset) carries force = **k·(r−0) = k·r**. At soft k
  that was small (~0.7 pN); at the physical stiff k it is k·r ≈ 322 pN — UNPHYSICAL (a real
  cross-bridge binds force-free and generates ~F_stall via the nm-scale power stroke, not k·r
  over the whole binding distance).
- **The transmission channel stayed floored** (g_actin 1.68e-3, ~unchanged) even with the stiff
  cross-bridge — confirming PI's point that the **crosslink (fiber↔fiber, k_intra=k_attach=1e-7)
  is the SEPARATE soft transmission link**, downstream of generation.

**⇒ The unit-slip fix CANNOT be a config-only change.** The correct fix is ATOMIC:
1. **k_head_spring/k_head_actin = 1e-3** (the unit-slip correction, literature/canonical-anchored), AND
2. **grip_walk r0 = BOUND length** (store r_bind at binding; r0_eff = r_bind − s_grip; force =
   k·(r − r0_eff) = k·s_grip, the power-stroke force, capped at F_stall by the Hill stall) — so a
   freshly-bound head is force-free and the force is the nm-scale stroke, NOT k·r. (This is O2.)
3. **CFL re-derivation** (k_backbone=1e-2 → τ_backbone≈39 ns native; dt via cortex cfl_safety).
4. **Crosslink k_intra/k_attach** re-anchor (the fiber↔fiber transmission link; verify α-actinin
   literature — separate from the cross-bridge, but the SAME soft-coupling class).

Config reverted to 1e-6 pending the atomic fix (a config-only 1e-3 leaves an over-tensioned broken
state). The physical end-state: per-head → F_stall (not k·r), γ → the full-stall envelope
(18–36× under band = the density/overlap gap), with transmission gated by the crosslink fix.

## 9. ⛔ The bin-quantized attach bond CANNOT deliver the stiff-cross-bridge stroke (2026-06-09)

Attempted the atomic fix (stiffness 1e-3 + grip_walk r0 = max(r_bind − s_grip, 0) via re-binning).
Implemented all five edits, then found a HARD resolution wall in the bin architecture:

- The stiff-cross-bridge stroke is tiny: s_grip_stall = F_stall/k_series ≈ **4 nm** (the head
  walks only ~4 nm before the Hill stall caps it at F_stall).
- But the myosin binding range is LARGE: `head_actin_max_bind_dist = 3.3e-7` = **330 nm** (the
  head sits 200 nm off the backbone and binds actin up to ~330 nm away).
- The attach-bond r0 is QUANTIZED into `n_bins` types over [0, 330 nm]. Force granularity per
  bin = k·bin_width. To resolve the 4 nm stroke (granularity < F_stall = 8.48 pN) needs
  bin_width ≲ 2 nm → **n_bins ≳ 165** (165+ HOOMD bond types — impractical). At n_bins=60,
  bin_width=5.5 nm → granularity 23 pN ≈ **2.7× F_stall**: a single re-bin overshoots the stall
  force. The bins are fundamentally too coarse to deliver a controlled F_stall.

**⇒ The unit-slip fix needs a CONTINUOUS per-head force, not the bin scheme.** Three nested
findings now define the correct fix:
1. stiffness k_head 1e-6 → 1e-3 (the unit slip);
2. r0 = bound-length (force = k·power-stroke, not k·r) — required, else k·r explodes;
3. **continuous per-head force** — the bin quantization can't deliver the 4 nm stroke; replace
   the harmonic attach-bond-as-force with a CUSTOM FORCE delivering, per engaged head,
   F = min(k·s_grip, F_stall) along the head→bead unit vector (reaction on the head/backbone),
   with the attach bond retained ONLY for the Bell-Evans off-rate bookkeeping. (Plus: the 330 nm
   binding range is itself unphysical for a stiff cross-bridge — a stiff bridge binds at low
   strain ~nm; head_actin_max_bind_dist likely needs shrinking too, a coupled binding-geometry
   item.)

This is a deeper myosin-contractility redesign (custom force) than a config + bin change. All
edits reverted (no broken/under-resolved intermediate). Fix fully specified above; execute as a
focused, sanity-gated effort (md.force.Custom = per-step host-sync, the GPU-main-port territory —
so the cupy/native path is the production form). The CFL re-derivation (k_backbone=1e-2 →
τ≈39 ns) and the crosslink k re-anchor remain part of the atomic set.

## 10. ✅ EXECUTION (2026-06-09 loop14-17): continuous_stroke built + generation-fix verified

The §9 fix is built as opt-in `stepping_mode="continuous_stroke"` and the GENERATION
half (items 1-4) is verified; item 5 (crosslink) is surfaced to PI below.

- **Item 1+3 (continuous custom force) — DONE, hermetic (loop14, commit 8dc2680).**
  `MyosinHeadForce(md.force.Custom)`: per engaged head F=min(k·s_grip, F_stall) along
  head→bead (reaction on head, Newton-3). Attach bond demoted to k=0 (Bell-Evans
  off-rate + nlist exclusion only; Bell-Evans uses the delivered force). 13 sanity-gate
  tests PASS (dimensional, boundary 0@s=0 + cap@F_stall, Newton-3, contractile sign,
  delivered |F|=F_stall NOT k·r, legacy byte-identical). grip_walk/binned_r0 unchanged.
- **Item 2 (k 1e-6→1e-3) — DONE, MODE-COUPLED (loop16, commit ab376de).**
  `k_cross_bridge_continuous=1e-3` in phase1_h3.yaml (Magic-Number Block: canonical
  bridge/motor + Veigel 2002/Kaya 2010/Finer 1994, 0.3-2 pN/nm). Applied ONLY in
  continuous_stroke (the stiff k explodes the harmonic k·r in legacy modes = loop12).
- **Item 4 (CFL) — DONE (loop16).** `reconcile_dt` folds myosin k_backbone (=1e-2,
  τ_backbone=9.2ns meso) → integrator dt lowered 1.30e-8→9.22e-10 s (14.1×), NO edit
  to the frozen integrator/. Smoke STABLE (3300 steps, all finite, NO blowup);
  per-head delivered T=1.15pN tracking k·s_grip, capped << F_stall=8.48pN — NOT the
  k·r=930pN explosion. The generation fix works end-to-end.
- **Binding range 330nm — KEPT (not shrunk).** The §9 shrink rationale was bin-
  resolution-specific (r0 quantization can't resolve the 4nm stroke over 330nm). With
  the continuous custom force, r does NOT enter the force → no explosion → no shrink
  needed; the head reach is already gated by capture_perp=210nm≈head_rest_length.

### ⛔ Item 5 (crosslink k re-anchor) → PI DECISION (gate-contract / production-wide)

KB verification (tag_query.py, 2026-06-09):
- **KB-1.28** (verified, High): the shared crosslink bond Hamiltonian H=(k_xl/2)|r|²
  with **k_xl ≈ 1e-4–1e-2 N/m, default 1e-3 N/m (1 pN/nm)**.
- The cortex `dynamic_crosslinkers.k_intra = k_attach = 1.0e-7 N/m` is **~10⁴× below**
  the KB-1.28 default and 10³× below its range floor (1e-4).
- The config's anchor comment ("0.1 pN/µm, KU-3.19 (Furuike 2001)") does NOT match a
  stiffness datum: KB-3.19 specifies crosslinker KINETICS (off-rates) only — NO k; and
  Furuike 2001 is filamin UNFOLDING kinetics, not a 0.1 pN/µm spring constant. This
  looks like the same mis-attribution class the 2026-06-02/06-08 audits found.

⇒ The 1e-7 crosslink stiffness is very likely a soft-coupling slip (§8 named it the
"SEPARATE soft transmission link"). Re-anchoring to KB-1.28 (→ 1e-3 N/m) is the
literature-first call. BUT it is a **production-wide gate-contract change**: it affects
EVERY cortex build (grip_walk production + the other session's compartment platform +
the already-run Gate-A/B), and the molecular stiffness of a FLEXIBLE crosslinker vs the
simulation harmonic-bond stiffness are distinct concepts (the KB range spans 100×). CFL
is not tightened (τ_xl=391ns at 1e-3 ≫ the 9.2ns myosin dt). Per the mission rule
("crosslink k는 문헌검증 후, 불확실시 PI surface" + gate-contract→PI sign-off), NOT
changed inline. **PI: re-anchor k_intra/k_attach 1e-7→1e-3 (KB-1.28), or confirm the
1e-7 with a correct single-molecule anchor?** Note: the generation fix needs no crosslink
change to verify; crosslink stiffness is the downstream transmission lever (smoke WALL-A
propagation already 66.8% in the connected mesh).

### Remaining: M4 GPU magnitude validation (gbook)
The cap is verified; the MAGNITUDE (s_grip→F_stall over ~1e7 steps → γ → the full-stall
envelope, ~18× under band = the KNOWN density/overlap gap, band LOCKED) needs the long
GPU contraction run (Mac→gbook rsync, PI-approved remote-overwrite). Best run AFTER the
crosslink decision (the transmission lever co-determines the network γ).
