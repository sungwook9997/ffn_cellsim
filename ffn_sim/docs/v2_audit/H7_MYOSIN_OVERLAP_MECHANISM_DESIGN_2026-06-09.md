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
