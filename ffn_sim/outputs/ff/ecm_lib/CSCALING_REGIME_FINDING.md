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

## Bottom line

The reference-concentration Pa is correct (validated). The c-scaling gap is understood: it is the
stretch↔bending regime + rigidity-percolation, bracketed 1.07 (stretch) ↔ 7.24 (bending-through-threshold),
with real c²≈2 requiring a PI-decided ⟨z⟩(c) crosslink law kept above threshold. No tuning was applied.
