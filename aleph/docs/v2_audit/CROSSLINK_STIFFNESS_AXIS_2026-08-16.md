# The crosslink stiffness is not a constant to correct — it is an axis, and our point is off it

**Status:** FINDING + AXIS DRAFT. Decides nothing. Written 2026-08-16 by the `arena/solver` session at PI
request, after the PI observed that under standing ruling #1 (*every physiological value is a declared axis
with a band and the scope that band is true for*) this is an axis declaration, not a value correction.

**Why it matters here:** `k_xl` is what sets `kmax`, `kmax` sets `dt_mu = 0.1/kmax`, and `dt_mu` is what the
explicit relaxation steps by. It is the single constant standing between the resting solve and a tractable
timestep.

---

## 1. What the production cortex actually runs, and what it is

Measured on the committed full-native record `arena_prestudy/s1_pressure_on.json`:

| | value |
|---|---|
| median crosslink `k_xl` (running) | **8.2e5 pN/µm = 820 pN/nm** |
| `kmax_pN_per_um` | 2,828,099.29 — **crosslink-dominated** |
| `dt_mu = 0.1/kmax` | 3.536e-08 |

`kim_network.py:140` attributes the value: *"Lit anchors (Ferrer 2008 PNAS, AFM): α-actinin 455 pN/nm =
4.6e5 pN/µm, filamin 820 pN/nm = 8.2e5 pN/µm."*

### 1.1 The numbers ARE in Ferrer 2008 — and they are not a spring constant

Retrieved through PubMed (PMID 18591676 / PMC2453742, doi:10.1073/pnas.0706124105). The paper says:

> *"the α-actinin/actin interaction is slightly more flexible (k_B T κ_m = **455 ± 215 pN/nm**) than the
> filamin/actin interaction (k_B T κ_m = **820 ± 551 pN/nm**)"*

Reported in the same Results: **transition distances 1.94–2.75 Å**, **free-energy barriers 3.6–4.3 k_BT**,
rupture forces 40–80 pN, intrinsic off-rates 0.066–0.087 s⁻¹.

The harmonic curvature of a barrier of height ΔG and width x_β is `k = 2ΔG/x_β²`. With k_BT = 4.114 pN·nm:

| | ΔG | x_β | 2ΔG/x_β² | paper |
|---|---|---|---|---|
| filamin | 3.75 k_BT | 1.94 Å | **820 pN/nm** | 820 ± 551 |
| α-actinin | 4.30 k_BT | 2.75 Å | **468 pN/nm** | 455 ± 215 |

Both reproduce from the paper's own co-reported quantities. **455 and 820 pN/nm are the curvature of a
binding free-energy barrier over ~2 Å** — a Bell-Evans landscape parameter. The paper's own wording treats
them as a property of the *interaction*, not of the protein as a structural element.

⚠ **The engine uses them as the spring constant of a structural element spanning two filaments.** Filamin
is a ~160 nm flexible dimer; α-actinin a ~35 nm rod. The barrier is ~0.2 nm wide. The length scales differ
by **800× and 127×**.

⚠ **This is not the unit slip the code already documents.** `kim_network.py:141-143` warns of *"a confirmed
pN/µm-vs-pN/nm slip"*. The unit conversion is correct — 820 pN/nm IS 8.2e5 pN/µm. What is wrong is which
physical quantity the 820 describes. Fixing units changes nothing.

**Not claimed:** that `k = 2ΔG/x_β²` is the formula the paper derives. It is a back-calculation that
reproduces both values from co-reported quantities, which is strong but is not the paper's own derivation
read in place. The notation `k_BT κ_m` is also ambiguous as printed.

---

## 2. What the field uses for the same element

| source | crosslinker | value |
|---|---|---|
| AFINES | α-actinin (long, flexible) | **0.1 pN/µm** |
| AFINES | fascin (short, stiff) | 1.0 pN/µm |
| AFINES | strong-antiparallel / parallel | 2 – 5 pN/µm |
| aLENS | motor tether κ_xl | **100 pN/µm** |
| **this engine** | `k_xl` default | **10 pN/µm** |
| **this engine** | `hand_kmc` `link_k` | **0.1 pN/µm** |
| **this engine** | **production, running** | **8.2e5 pN/µm** |

**Field band ≈ 0.1 – 100 pN/µm.** Production sits **8,200× to 8,200,000×** above it.

⚠ **The engine's own defaults are already inside the band.** `kim_network.py` calls them *"UNSOURCED
magnitude knobs … an AFINES soft surrogate"* and names the Ferrer values as the *"lit anchors"* to correct
toward. That is backwards: the in-band values are what the field uses for network-element compliance; the
Ferrer values are barrier curvatures.

### 2.1 The independent literature says the physics is inverted too

* Kasza 2009 — networks are *"rigid filaments connected by multiple **flexible** linkers"*; individual
  crosslink loads **< 10 pN**.
* Gardel 2006 — filamin's **flexible hinge** is a molecular requirement for cell-like mechanics.
* Wang 2023 — network shear depends on filamin's bending and tensile stiffness and is
  *"almost insensitive to the tensile stiffness of actin filaments"*.

Our engine has the opposite ordering. Actin backbone `EA/L_seg` = 4.4e4 / 0.5 = **8.8e4 pN/µm**; crosslink
**8.2e5 pN/µm** — the crosslink is **9.3× stiffer than the filament it connects**, while `kim_network.py`'s
own docstring calls `k_xl` *"the soft deformable element"* and states the model's valid regime as
`k_xl ≪ EA/L_seg`. **We are an order of magnitude into the opposite regime.**

⚠ A second, compounding deviation is recorded in the same ledger:
`cortex_seg_um_deviation = "5–10× TOO COARSE (0.5 µm vs 0.05–0.10 µm)"`, provenance CONVENIENCE. At the
sourced segment length the backbone would be ~5.9e5 pN/µm and the ratio 1.4× rather than 9.3× — the coarse
segment makes the backbone artificially soft relative to the crosslink.

---

## 3. Axis draft — band and scope only, no value selected

* **Axis:** `k_xl`, crosslink junction stiffness [pN/µm].
* **Band:** ≈ 0.1 – 100 pN/µm.
* **Scope the band is true for:** reconstituted in-vitro F-actin networks with the crosslinker treated as a
  **compliant network element** spanning two filaments. Sources: AFINES parameterisation (0.1–5), aLENS
  tether (100).
* **Explicitly NOT on this axis:** Ferrer 2008's 455 / 820 pN/nm. Different physical quantity — the
  curvature of a ~2 Å binding barrier, not the compliance of a ~35–160 nm protein under network load.
* **Leverage:** `dt_mu = 0.1/kmax`, so moving from 8.2e5 into the band raises the timestep by **8,200× to
  8,200,000×** *if the crosslink still dominates kmax there* — which is exactly what is not yet known.

⚠ Declaring an axis is not free: standing ruling #1's own bound is that at most the SVD rank of the
log-coordinate Jacobian, computed BEFORE the sweep, may be reported as inferred.

---

## 4. Option ③ — treat crosslinks as constraints, not springs (aLENS)

**What aLENS does.** Their stated motivation is our problem verbatim: tether relaxation `1/λ ≈ 3e-6 s`, and
*"explicit timestepping schemes require Δt < C/λ"*. Their answer is to stop integrating stiff springs:
crosslinkers become **bilateral constraints** `K[Φb(C) − Φb0] = −γb`, steric contact becomes **unilateral
complementarity** `0 ≤ Φu ⊥ γu ≥ 0`, and both go into **one convex QP** solved by parallel Barzilai-Borwein
projected gradient descent. Claimed gain: *"timesteps two or more orders of magnitude larger than currently
available"*. Demonstrated to 1e6 filaments / 3e6 motors.

**Where this engine already stands.** `project_constraint_forces_kernel` implements the exact NF2007
projector `P = I − Jᵀ(JJᵀ)⁻¹J`, and it works because a filament is a **chain**: `JJᵀ` is symmetric
tridiagonal and one thread per fiber runs a Thomas solve. **Crosslinks break that** — the constraint graph
becomes a general sparse graph and `JJᵀ` is no longer tridiagonal, so the Thomas solve does not extend.
aLENS's answer is not to form the projector at all but to solve iteratively with matvecs only. This engine
already has matrix-free projected iterative machinery (`SFImplicitCG`, `ProjectedAnalyticCG`), so the work
is extending the projection from the per-fiber chain to the crosslinked graph rather than building a solver.

**⚠ Option ② does not remove the need for ③.** aLENS's own tether stiffness is **100 pN/µm — inside the
band we would correct into** — and they still found explicit intractable and went to constraints. Correcting
`k_xl` from 8.2e5 to ~100 lands us at aLENS's starting point; constraints are what they did *from there*.

② repairs the physics. ③ repairs the integration. They are different problems and neither substitutes.

---

## 5. Recommended order, and the one measurement that gates ③

1. **Declare the axis (② )** — PI decision. Band and scope above; Ferrer excluded as a different quantity.
2. **Read `kmax` at band endpoints.** `kmax_pN_per_um` is computed at BUILD, so this needs a build and no
   run and no accepted step. It answers the question nobody can currently answer: **when the crosslink stops
   dominating, what dominates?** Bending, steric, membrane and myosin terms all feed the same max.
3. **Decide ③ from step 2's answer.** If `kmax` stays large with `k_xl` in band, the dominant term is
   something else and constraint treatment must target *that*, not the crosslinks.

⚠ **Do not go straight to ③.** We know the crosslink dominates `kmax` today. We do not know what dominates
without it, and constraint-treating the wrong family is expensive and unfalsifiable.

---

## 6. What this document does not do

It selects no value, sets no default, edits no constant, and does not touch `kim_network.py` or
`hand_kmc.py` — the production constant is PI-gated by the code's own comment and stays that way. It also
does not claim the corrected axis makes the resting solve converge: that question is separately blocked on
(e) 1, where the acceptance criterion itself is undecidable as written
(`PI_DECISION_e1_STATIONARITY_2026-08-16.md`).

**Sources.** Ferrer et al. 2008 PNAS via PubMed, doi:10.1073/pnas.0706124105 · aLENS, Yan et al. 2022 eLife
74160 · AFINES, Freedman et al. 2017 Biophys J · Kasza et al. 2009 PRE · Gardel et al. 2006 PNAS · Wang et
al. 2023 Extreme Mech Lett.
