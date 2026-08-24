# H.10 — Cytoplasm viscoelasticity / poroelasticity — DESIGN (read-only prep)

> **STATUS: DESIGN doc — 2026-06-02. NOT a contract, NO code written.** Produced in a
> parallel prep session (autonomy log `AUTONOMOUS_LOG_2026-06-02.md`); no existing
> file modified, no `integrator/` touched. Feeds the Notion `MC-H10-cytoplasm`
> (`draft`) contract toward `draft→implemented`. Companion to the brief
> `briefs/H10_cytoplasm.md`, the EXTEND memo `CELL_MECHANICS_EXTEND_VS_REBUILD.md`
> §B.5, and the PI-exp map (Hu 2024 viscosity bands). **Every code block is a design
> artifact, NOT wired.**

## Why H.10 is the hard one (and why it is design-only here)

Membrane (H.8) and nucleus (H.9) are *forces over a tag range* (Template 1) — clean.
The cytoplasm is **not a force; it is the background medium** the whole cell moves
through, and a real cytoplasm has **memory** (it is viscoelastic: stress depends on
strain *history*, not just instantaneous strain) and is **biphasic / poroelastic**
(a solid cytoskeletal network bathed in interstitial fluid). The current runtime
models the medium as a *single* Markovian water-viscosity Stokes drag `γ_b = 6πηR`
per particle type, fed to the BAOAB integrator via `gamma_map`. Capturing
viscoelasticity faithfully reaches into the **drag model / integrator**, which is
**frozen** (CLAUDE.md; PI-gate). So this unit is scoped separately and only the
*design* is produced autonomously; the build is PI-gated.

## What the cytoplasm must reproduce (KU-3.B3, literature-anchored)

All bands are **acceptance oracles**, never fitting targets. Probe-length matters:
microrheology viscosity `η ∝ L²` — state the wire length `L`.

| Quantity | Band (lit) | Source | Conf |
|---|---|---|---|
| Cytoplasm **viscosity η** (MRS, L=3±1 µm wires) | MCF10A 41.6 / MCF7 56.4 / MDA 10.7 Pa·s | H.10 brief / KU-3.B3.1 | — |
| Cytoplasm **viscosity η** (MRS, refined) | MCF-7 **65.9±11.4** vs MDA **12.0±5.7** Pa·s (~5×; band **12–66**) | Hu 2024 *Nanoscale Adv* PMC10929591 | high (3-0) |
| Cytoplasm **elastic shear modulus G** | MCF10A 79.3 / MCF-7 32.9 / MDA 38.6 Pa (band **33–79**) | Hu 2024 PMC10929591 | high (3-0) |
| **Poroelastic diffusion Dp** | 40–60 µm²/s (4–6×10⁻¹¹ m²/s); pore ξ≈14 nm; poroelastic regime for events faster than ~0.5 s | H.10 brief / KU-3.B3.2 | — |
| (contrast anchor, distinct quantity) Whole-cell AFM E₀(1 s) | MCF10A 1.14 / MCF-7 0.26 / MDA 0.46 kPa; fluidity 0.186/0.234/0.147 | Yubero 2020 *Commun Biol* | high (3-0) |

**The discriminator:** cytoplasm **viscosity** (MCF7 ≈ 5× MDA) is the single most
robust metastatic separator in the dataset — far cleaner than whole-cell stiffness
*ordering*, which is method-dependent and disputed (AFM → MCF7 softer;
electrodeformation → MCF7 stiffer; PI-exp map §Contested). So H.10's primary
validation target is the **viscosity contrast**, not an absolute stiffness.

**Three DISTINCT moduli — do not conflate** (a recurring error this doc guards):
1. cytoplasm **viscosity η** (Pa·s) — the loss/dissipative part, dominates this unit;
2. cytoplasm **shear modulus G** (Pa, ~30–80) — intracellular elastic storage;
3. whole-cell **AFM modulus E₀** (kPa, ~0.3–1) — the composite cortex+membrane+nucleus
   response, ≈10–20× above G, an H.3/H.8/H.9 composite quantity, NOT H.10.

## The current medium model (code-grounded, what exists)

- BAOAB reads a `gamma_map` keyed by particle type; each type gets a single scalar
  `γ_b = 6π η R` (Stokes), with η = **water** (~1e-3 Pa·s). BAOAB asserts every present
  type has an entry and is the sole position integrator (`methods=[]`).
- Fluctuation–dissipation is built in: BAOAB's per-type thermal noise amplitude is set
  from that same per-type `γ_b` and `kT`, so the equilibrium diffusion is `D = kT/γ_b`.
  **This is the load-bearing fact for the cheapest option below.**
- The medium is therefore **Markovian (memoryless) + monophasic + water-viscous**.
  Missing: memory (viscoelastic G(t)), the elastic storage modulus, biphasic pore flow.

## Three attach options (ranked by fidelity vs integrator-invasiveness)

### Tier 1 — per-type **effective viscosity** `η_eff` in `gamma_map` (cheapest; gamma_map-only)

**Idea.** Replace the water η with a cytoplasm η for cytoplasm-immersed bead types:
`γ_b = 6π η_eff R`, η_eff ∈ {33–66 Pa·s}. This is **~5×10⁴× the water drag** — and it
*is* the metastatic discriminator (set η_eff per cell-type → MCF7 ≈5× MDA emerges in
the diffusive dynamics).

**Why it is *almost* free.** It uses the **existing `gamma_map` contract** (Template-2
step 2) — a *value* change, not an integrator-code change. FDT is preserved
automatically: BAOAB scales its noise from the same per-type `γ_b`, so `D = kT/γ_b`
stays correct as η_eff rises (slower diffusion, the right sign). Default η_eff = water
⇒ pre-H.10 runs bit-for-bit identical.

**What it captures / misses.** Captures the **viscous (loss) contrast** — the
discriminator. Does **NOT** capture the **elastic storage G** (no memory: stress
relaxes instantly) nor poroelastic pore-flow. So Tier-1 is a *viscosity-contrast*
model, honest about being loss-only.

**The PI-gate (why it is not auto-built).** A 5×10⁴× drag jump changes the velocity
relaxation time `τ_v = m/γ_b` by 5×10⁴ and the diffusive exploration per step. In the
overdamped BAOAB limit this is *intended* (the cell is deeply overdamped in cytoplasm),
but it (a) is a large, physics-altering constant that "touches the drag model" (brief
§B.5 flags a PI-gate), and (b) interacts with every other CFL/timescale gate already
tuned at η_water. **Surface the η_eff magnitude + a re-run of the existing
timestep/CFL gates to PI before enabling.** Not a magic number (η_eff is the measured
Hu-2024 value), but a frozen-medium-policy change.

```python
# DESIGN ONLY — NOT a module. The change is a gamma_map value, default = water.
# In the (Lead-owned) gamma block of build_cortex_full_simulation, gated default-off:
#   eta_eff = p_cytoplasm.eta_cyto if p_cytoplasm else ETA_WATER
#   for t in cytoplasm_immersed_types:
#       gamma_map[t] = 6.0 * math.pi * eta_eff * radius_of[t]
# p_cytoplasm = None  ->  eta_eff = ETA_WATER  ->  bit-for-bit identical.
```

### Tier 2 — auxiliary-variable **Markovian GLE** for the elastic storage `G` (integrator-adjacent; PI-GATE)

**Idea.** True viscoelasticity = a Generalized Langevin Equation with a memory kernel
`K(t)` whose Laplace transform gives the complex modulus `G*(ω) = G' + iG''`. BAOAB is
Markovian and cannot carry memory directly. The standard fix is a **Markovian
embedding**: add per-bead auxiliary DOFs `s_k` obeying
`ṡ_k = -ν_k s_k - c_k v + noise`, coupling back a force `Σ c_k s_k` on the bead, so the
*marginal* dynamics reproduce an exponential-sum memory kernel `K(t)=Σ c_k² e^{-ν_k t}`
(Prony series fit to the measured `G(t)`; Ceriotti 2010 "colored-noise"; Baczewski &
Bond 2013). One or two modes already give a storage modulus + finite relaxation time.

**Why PI-gated.** This is a **new integration scheme** — it adds DOFs and modifies the
per-step propagation. It is *additive + default-off* in principle (zero auxiliary modes
⇒ recovers BAOAB exactly), but it **edits the frozen `integrator/`** (or adds a
sibling integrator), which is a CLAUDE.md PI-only gate. **Design only here.** Maps to
`G ∈ 33–79 Pa` (the storage modulus the Prony series must reproduce) + the relaxation
time set by `η/G` (~0.5–2 s, consistent with the poroelastic ~0.5 s crossover).

### Tier 3 — explicit **poroelastic filler particles** (most fine-grained; CLAUDE.md-ideal; heaviest)

**Idea.** The CLAUDE.md "full-fidelity, fine-grained, mechanistic" option: model the
biphasic cytoplasm explicitly — a solid cytoskeletal mesh (already partly present:
cortex + crosslinkers) plus **interstitial cytosol filler particles** that carry the
fluid phase, with a drag/permeability coupling reproducing Darcy flow through the pore
network (ξ≈14 nm) and the poroelastic diffusion `Dp = 40–60 µm²/s`. Template-2 (new
`cytosol_bead` type + `gamma_map` entry + redistribution Action). The poroelastic
timescale `τ_p ~ L²/Dp` (L≈3 µm ⇒ τ_p ≈ 0.18 s) then *emerges* rather than being
imposed, matching "poroelasticity dominates for events faster than ~0.5 s."

**Cost.** Adds a large particle population (densest of the three) and a pore-flow
coupling — the heaviest to build and run, but the only option that gives both storage
*and* pore-flow mechanistically. Best deferred until a poroelastic experiment (AFM
indentation relaxation, optical-tweezer creep) actually needs the pore-flow regime.
Template-2 contract; still default-off (no filler ⇒ identical), but real new physics.

## Recommended staging

1. **Tier-1 (η_eff drag contrast)** — the discriminator (MCF7 ≈5× MDA), gamma_map-only,
   FDT-safe by construction, default-off. **PI-gate = the η_eff magnitude + a re-run of
   the existing CFL/timestep gates at the new drag.** This alone validates KU-3.B3.1
   (the viscosity contrast) — the highest-value, lowest-risk first step.
2. **Tier-2 (Markovian-GLE storage)** — adds the elastic modulus `G` + finite
   relaxation. **PI-gate = frozen-`integrator/`** (new scheme / sibling integrator).
3. **Tier-3 (poroelastic filler)** — full biphasic mechanism for the pore-flow regime
   (`Dp`, ξ). Defer to a poroelastic experiment.

## Sanity Gate (to be satisfied at build time — recorded now per CLAUDE.md)

1. **Dimensional.** `γ_b = 6π η R` → [Pa·s · m] = [N·s/m] ✓. `D = kT/γ_b` → [m²/s] ✓.
   GLE: `c_k²/ν_k` has units of a force-correlation ⇒ a modulus after the area factor.
2. **Boundary.** `η_eff → η_water` (or `p_cytoplasm=None`, zero GLE modes, no filler)
   **must recover the current run bit-for-bit** — the default-off contract.
3. **Conservation / FDT (the critical one).** Tier-1: BAOAB's noise amplitude is set
   from the *same* per-type `γ_b`, so raising η_eff keeps `⟨v²⟩ = kT/m` and `D=kT/γ_b`
   exact — **verify BAOAB reads γ per-type from `gamma_map` for the noise term**, not a
   global γ (if it ever used a global γ for noise, Tier-1 would silently break FDT).
   Tier-2: the auxiliary-DOF noise must satisfy the embedded FDT (`⟨ξ_k ξ_k⟩ ∝ ν_k kT`)
   or temperature drifts — the core correctness test for the GLE.
4. **Numerical / CFL.** `τ_v = m/γ_b` shrinks ∝ 1/η_eff; confirm the overdamped BAOAB
   limit still holds and that the *existing* CFL/stiffness gates (ERM, enclosed-volume,
   membrane) re-pass at the new drag (they are τ = γ_b/k_eff gates — a larger γ_b
   *relaxes* them, so likely safe, but re-assert, do not assume).
5. **Sign-sense.** Higher η_eff ⇒ slower diffusion / more sluggish relaxation
   (MCF7 more viscous than MDA ⇒ MCF7 beads diffuse slower). Wrong sign = bug.
6. **Measurement-protocol consistency.** Report η at a stated probe length L (η∝L²);
   compare G to *intracellular* microrheology G, never to whole-cell AFM E₀.

## Validation acceptance (candidate gates — Notion VG-H10)

| Gate | Criterion | KU |
|---|---|---|
| Cytoplasm viscosity contrast | per-type η_eff ⇒ MCF7-tagged bead diffusion ≈ 5× slower than MDA-tagged (D = kT/γ_b ⇒ ~5× lower D) | KU-3.B3.1 |
| Cytoplasm modulus band (Tier-2) | emergent G ∈ 33–79 Pa, ≪ whole-cell AFM E₀ | KU-3.B3.1 |
| Poroelastic timescale (Tier-3) | relaxation faster than ~0.5 s; emergent Dp ∈ 40–60 µm²/s | KU-3.B3.2 |

## Open items for PI

- [ ] Ratify the **staging** (Tier-1 first; Tier-2/3 deferred).
- [ ] **PI-gate Tier-1**: approve the η_eff magnitude (33–66 Pa·s, cell-type-set) +
      a re-run of the existing CFL/timestep gates at the raised drag.
- [ ] Confirm Sanity Gate #3 against the BAOAB source: noise term reads γ **per-type**
      from `gamma_map` (FDT preserved under per-type η_eff). If not, Tier-1 needs a
      noise-term change → escalates to a frozen-`integrator/` PI-gate too.
- [ ] **PI-gate Tier-2/3**: frozen-`integrator/` (GLE scheme) / Template-2 (filler).
- [ ] Which bead types are "cytoplasm-immersed" (interior actin/myosin/xlink) vs
      membrane/cortex-shell (water-exposed) — sets which `gamma_map` entries change.
