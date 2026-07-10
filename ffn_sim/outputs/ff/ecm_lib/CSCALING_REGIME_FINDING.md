# Collagen G'(c) scaling — root-caused to the mechanical regime (2026-07-10)

The FF ECM library reproduces collagen's **absolute modulus at reference concentration** (G'(1.5 mg/mL)≈11 Pa,
Yang-Kaufman) but the concentration **exponent** was n≈1.07 vs the literature n≈2.0. This note root-causes
that gap — it is a **mechanical-regime + rigidity-percolation** effect, not a bug, now bracketed from both sides.

## The regime bracket (collagen, box 24 µm, 6 concentrations)

| axial mode | physics | measured exponent n | why |
|---|---|---|---|
| `spring` (finite EA, default) | fibers **stretch** to carry shear | **1.07** (too shallow) | stretch-dominated: σ_xz carried by fiber axial stretch (seg 1.59 / xl 0.009 / bending ~0 Pa @γ=0.1). At fixed geometry the stress ∝ density ⇒ ~c¹. |
| `reshape` (inextensible) | fibers **bend** (can't stretch) | **7.24** (too steep) | bending-dominated, but ⟨z⟩ sweeps the **rigidity-percolation threshold** across the c-range: c=1→⟨z⟩=2.18 (sub-critical, floppy, G≈0) → c=7→⟨z⟩=14.2 (super-critical, 5858 Pa). The steepness is the threshold crossing, not the asymptotic scaling. |

**Real collagen (n≈2.0) sits between these** — bending-dominated (like `reshape`) but measured **above** the
rigidity threshold in the developed regime (like a controlled `spring`), so it neither stretch-saturates
(c¹) nor crosses percolation (c⁷).

## What the correct c² needs (a modeling decision — PI-scoped, not auto-tuned)

The exponent is set by how connectivity ⟨z⟩ scales with concentration:
- At **fixed ⟨z⟩** (bending-dominated, above threshold): ℓc=L_f/⟨z⟩ fixed, ρ∝c ⇒ G∝c¹.
- Real collagen's c² comes from ⟨z⟩ **increasing modestly** with c (more fiber crossings) **while staying
  above the rigidity threshold** — the combined ρ∝c and ℓc↓ gives ~c².

So reproducing c² requires a **crosslink-density-vs-concentration law** that (a) keeps the network above the
sub-isostatic rigidity threshold at the lowest concentration (my current rule gives ⟨z⟩=2.18 at 1 mg/mL,
floppy), and (b) grows ⟨z⟩ modestly with c. This is a genuine modeling choice about collagen crosslink
biology (LOX density vs concentration) — it should be **PI-decided and literature-anchored**, not tuned to
hit n=2. Flagged for PI. The `reshape` mode's numerical stability with the affine-energy method would also
need fixing (it diverges — currently boundary-driven + virial is the stable reshape path used for this test).

## Physical-connectivity control added (`target_z`, KB-1.3 anchored)

The raw near-contact rule let ⟨z⟩ blow up unphysically with density (1 mg/mL → ⟨z⟩=2.18 floppy; 7 mg/mL →
⟨z⟩=14.2 over-connected). `build_fibrillar_ecm(target_z=3.2)` now subsamples fiber-pair crosslinks to the
**literature ⟨z⟩~3.2 (KB-1.3, roughly concentration-independent)** — collagen crosslinks only a LOX-set
subset of geometric contacts, not every one. This is anchored to an existing KB datum, NOT tuned to n.

Result at fixed physical ⟨z⟩=3.2 (box 26 µm):

| c (mg/mL) | ⟨z⟩ | G (Pa) |
|---|---|---|
| 1.0 | 1.95 (can't reach 3.2 — genuinely under-connected at low c) | 8.3 |
| 2.0 | 3.20 | 18.5 |
| 4.0 | 3.20 | 36.8 |
| 7.0 | 3.20 | 64.1 |

Exponent **n=1.04** — the clean fixed-connectivity prediction (G∝density∝c at fixed ⟨z⟩), no longer
confounded by the ⟨z⟩ blow-up. Reference G(1.5 mg/mL)=13.4 Pa (⟨z⟩=3.02) still matches Yang-Kaufman.

## Independent confirmation: the first normal-stress sign N1 (via the virial tensor)

The full virial Cauchy tensor (`ecm_material_stress`) gives a SECOND, independent diagnostic of the same
regime — the first normal-stress difference N1=σ_xx−σ_zz under shear (`shear_stress_curve`):
- **Model: N1 > 0** (collagen +0.6 Pa, fibrin +4.4 Pa @20% strain) — a STRETCH-dominated network pushes the
  plates apart (the Poynting effect, like rubber).
- **Real collagen/fibrin: N1 < 0** (Janmey 2007) — bending-dominated semiflexible networks pull the plates
  together (fibers buckle in compression, sustain tension → net inward pull), with |N1|~|σ_xz| at ~20% strain.

So the WRONG SIGN of N1 (and |N1|/|σ_xz|≈0.15-0.28 vs ~1) confirms — from a completely different observable
than the c-exponent — that the spring-mode model is stretch-dominated, not bending-dominated. `figs/normal_
stress_N1.png`. Two independent signatures (c-scaling exponent AND N1 sign) point to the same root, making the
regime diagnosis robust. Reaching the negative-N1, c² bending-dominated regime is the same PI-gated modeling
decision (⟨z⟩(c) growth law / stable inextensible dynamics).

## Bottom line

The reference-concentration Pa is correct (validated) and ⟨z⟩ is now physically controllable at the KB-1.3
value. The c-scaling gap is fully understood: **at physical fixed ⟨z⟩ the athermal model gives c¹** (density-
linear); the literature c²≈2 requires ⟨z⟩ to GROW with c (bending-dominated, above threshold) — a PI-decided,
literature-anchored crosslink-density-vs-concentration law (collagen LOX biology), NOT an n-fit. Bracketed
1.04 (fixed ⟨z⟩) / 1.07 (⟨z⟩ blow-up) ↔ 7.24 (reshape through threshold). Independently confirmed by the
positive N1 sign (stretch-dominated) vs the literature negative N1 (bending-dominated). No tuning was applied.
