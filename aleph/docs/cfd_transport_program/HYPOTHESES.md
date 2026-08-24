# HYPOTHESES — Intracellular Transport & Field Actuation of the Cytoskeleton

**Program:** ffn_cellsim CFD-engine pivot (fluid-first Biot poroelastic substrate + transported G-actin monomer field + emergent actomyosin network)
**Date:** 2026-07-16 · **Branch context:** `ff/mech-hierarchy` → `ac/foundation`
**Scope:** Each hypothesis is grounded in physics established this session (low-Re poroelastic CFD, membrane electrical shielding, Péclet transport analysis, momentum/scallop constraints, formin reaction-kinetics) and is stated to be falsifiable either by an in-silico native full-compartment run or by a literature/analytic contradiction. Magnitude verdicts inherit the density floor and the I0 magic-number gate; where a value is a target rather than a claim it is marked.

**Baseline physics constants used below** (session-established):
- Actin-monomer diffusivity `D_G ≈ 3–6 µm²/s`; cell length `L_cell ≈ 10 µm`, FA/protrusion length `L_FA ≈ 1 µm`.
- Cytosol Reynolds number `Re ≈ 1.5×10⁻¹⁰` (physiological), best-case fantasy `~10⁻²`; turbulence onset `Re ≳ 2000`.
- Membrane/interior conductivity ratio `σ_m/σ_i ≈ 10⁻⁷`; Schwan transmembrane response `V_m = 1.5 R E cosθ`.
- Poroelastic fluid: `∂p/∂t = c_v ∇²p`, Darcy `u = −(k/µ)∇p`; pressure-diffusivity `D_p ≈ 40–60 µm²/s`, `τ_p ~ 1–3 s`; drained-steady spatial pressure deviation `~1 Pa`; `K_drained ≈ 300 Pa`.
- Monomer transport field: `∂c/∂t = D_G∇²c − u·∇c + src − sink`; current engine G-actin = fixed scalar `20 µM` (infinite reservoir — a known fidelity gap).
- Single stress-fiber tension target: Kumar 2006, `10–30 nN`.

---

## H1 — Stress-fiber growth at a resting focal adhesion is REACTION-limited (formin/tension), not monomer-transport-limited

**Statement.** At a resting, adherent focal adhesion the elongation rate of a ventral stress fiber is set by formin processivity and mechanical (tension) gating of the barbed end, not by the rate at which G-actin monomer arrives. Diffusive resupply is far faster than consumption, so growth is insensitive to the local monomer flux.

**Rationale.** A single FA barbed cluster consumes only `~6%` of the diffusive monomer supply (`J ≈ 232 s⁻¹` consumed vs `I_max ≈ 3800 s⁻¹` diffusive arrival); diffusion refills the depleted zone `~260×` faster than it is drained. Formin FH1 domains locally concentrate monomer to `~10 mM`, bypassing the bulk pool, and thymosin-β4 sequesters `>80%` of monomer into a buffered reservoir that clamps free `[G-actin]`. The rate-controlling steps are therefore FH2 gating/processivity and the pN-scale tension dependence of barbed-end off-rate.

**Quantitative prediction.** Flux-elasticity of growth `d ln(growth)/d ln(flux) ≲ 0.1–0.2` across a physiological ±5–10× swing in local monomer flux. Equivalently, a `2×` change in delivered flux moves elongation rate by `<15%`. By contrast, a `2×` change in applied fiber tension or formin processivity moves it by order-1.

**Competing / null hypothesis (H1₀).** Growth is diffusion/transport-limited: `d ln(growth)/d ln(flux) → 1`, and starving the FA of monomer stalls the fiber proportionally.

**Falsification.** In the CFD engine, sweep the FA-local monomer sink/source (or `D_G`, or an imposed `u`) over ≥1 decade at fixed formin/tension and measure elongation elasticity. Elasticity `> 0.3` refutes H1. Wet-lab analogue: latrunculin/thymosin titration that lowers free `[G-actin]` should stall growth *before* any transport manipulation does — if instead local flow manipulation dominates, H1 fails.

---

## H2 — Advective enhancement of monomer delivery is significant only at Péclet ≳ 1 (v ≳ D/L ≈ 0.5 µm/s), i.e. in fast extended protrusions, not at FAs

**Statement.** Advection changes monomer delivery relative to diffusion only when the flow speed exceeds the diffusive velocity scale `D_G/L`. At the cell scale and at FAs this threshold is not met; it is met only in fast, geometrically extended protrusions.

**Rationale.** `Pe = vL/D_G`. With `D_G ≈ 3–6 µm²/s`: cell-scale threshold `v* = D_G/L_cell ≈ 0.5 µm/s`; FA-scale threshold `≈ 3–6 µm/s`. Measured animal-cell cytoplasmic velocities give `Pe ≈ 0.002–0.3 ≪ 1`. The lever exists only where both `v` is high and `L` is long: fast lamellipodia/extended protrusions reach `Pe ≈ 3` (2023 BiophysJ report of `2.3–3.2×` slower recovery half-times when advection is present), a regime where advective resupply is genuinely required.

**Quantitative prediction.** Delivery enhancement `≈ 1 + O(Pe)`: `<10%` at FAs (`Pe ≲ 0.1`), but `≈ 2–4×` in a `Pe ≈ 3` protrusion. Crossover at `v ≈ 0.5 µm/s` (cell) / `3 µm/s` (FA).

**Competing / null hypothesis (H2₀).** Advection matters everywhere flow exists (delivery scales with `v` at all `Pe`), including at FAs at physiological streaming speeds.

**Falsification.** Engine sweep of imposed `u` from `0.05 → 5 µm/s` at fixed geometry; delivery-enhancement curve must be flat below `v*` and rise only above it. A monotone, threshold-free rise from `v = 0` refutes H2. If FAs at `Pe ≈ 0.1` show `>25%` advective enhancement, H2 fails.

---

## H3 — An external DC electric field cannot drive intracellular cytosol flow or electrophorese monomers; galvanotaxis is signaling-mediated

**Statement.** A DC field applied outside the cell is shielded by the plasma membrane by `~10³–10⁴×`, leaving an interior field far too small to drive cytosol flow or monomer electrophoresis. Directional galvanotaxis is therefore a *signaling* response (membrane-receptor electrophoresis + PI3Kγ/PTEN polarity), not a direct cytoskeletal body force.

**Rationale.** Membrane/interior conductivity ratio `σ_m/σ_i ≈ 10⁻⁷`; a `1000 V/m` external DC field yields only `~0.1–1 V/m` inside (Schwan `V_m = 1.5 R E cosθ` drops the field almost entirely across the bilayer). Physiological ionic strength (150 mM, Debye length `~0.8 nm`) collapses the field-alignment mechanism that acts on the F-actin `~−4 e/nm` charge (Tang & Janmey 1996; alignment observed only at `500–1900 V/cm` in low-salt buffers). Documented galvanotaxis (fields `0.1–6 V/cm`) re-aims the cell via surface-receptor electrophoresis and PI3K/PTEN (Zhao 2006; Allen 2013) — the *same* membrane physics, acting at the surface, not the interior.

**Quantitative prediction.** Interior monomer drift Péclet `Pe ≈ 4×10⁻⁴ – 4×10⁻²` for physiological external DC fields — 1.5–3.5 decades below the `Pe ≈ 1` needed to matter. `Pe = 1` would require an interior field only reachable *after* electroporation (`~1 V` transmembrane, order `kV/cm` external), i.e. after the cell is permeabilized.

**Competing / null hypothesis (H3₀).** The external field penetrates enough to directly electrophorese monomers / drive cytoplasmic flow, and galvanotaxis is a direct electromechanical body force on the cytoskeleton.

**Falsification.** Couple a Schwan/Laplace membrane-shielding model to the CFD interior; predicted interior drift must stay `Pe ≪ 1` below the electroporation field. Wet-lab: galvanotaxis must survive interior-field nulling but die when receptor electrophoresis / PI3K–PTEN is blocked. If pharmacological block of PI3Kγ/PTEN leaves directional response intact while it scales with the *interior* field, H3 fails.

---

## H4 — Magnetic fields penetrate the cell but require embedded magnetic particles to act; rotating-bead actuation produces LOCAL stirring, not directed net flux, and cannot net-translate a cell

**Statement.** Because tissue is essentially non-magnetic (`μ_r ≈ 1`), a magnetic field passes through unimpeded but exerts negligible direct force on the diamagnetic cytoskeleton; actuation requires embedded magnetic particles. A rotating bead then produces *local* stirring/mixing of the cytosol, not a directed net monomer flux, and by momentum conservation cannot net-translate the cell.

**Rationale.** Direct diamagnetic force on a filament would need `~10–27 T` field gradients (physically absurd). Magnetic actuation is real only via tagged beads (magnetic twisting cytometry; Wang & Ingber 1993). A steadily rotating bead drives a closed, recirculating (zero-net-displacement) flow in the surrounding fluid; the famous "rotating actin filament" is myosin-driven (Nishizaka 1993), not field-driven. Any purely internal actuator conserves total momentum: recirculation ⇒ no center-of-mass translation.

**Quantitative prediction.** Around a bead rotating at frequency `f`, time-averaged net monomer transport `⟨∮ c·u dA⟩ ≈ 0` over a closed interior surface (net flux `< 1%` of the instantaneous circulating flux); local mixing (enhanced effective `D_eff`) can rise by order-1 in a `~1 bead-radius` shell but decays as `~1/r³`. Cell COM displacement from any internal bead protocol `≈ 0` absent external substrate traction.

**Competing / null hypothesis (H4₀).** A rotating internal bead produces a directed net monomer flux and/or can roll/translate the whole cell through internal flow alone.

**Falsification.** Engine: place a rotating-source in the Biot/advection field; integrate net monomer flux through interior control surfaces and track COM. A sustained directed net flux or nonzero COM drift (without an external traction term) refutes H4. Wet-lab: intracellular bead rotation should mix but not translate an untethered cell.

---

## H5 — Turbulent cytosol is impossible at cell scale (Re ≈ 10⁻¹⁰); internal flow cannot net-translate or roll a cell without external traction (scallop theorem)

**Statement.** At cellular scale the cytosol is deep in the Stokes regime; turbulence is categorically excluded. Consequently no internal, time-reversible-or-not flow can produce net self-propulsion or net rotation without external traction (Purcell's scallop theorem + momentum conservation).

**Rationale.** `Re ≈ 1.5×10⁻¹⁰` physiologically; even a fantastical `10⁻²` is `5–13` orders below the `Re ≳ 2000` turbulence threshold. The engine's fluid is Biot/Darcy diffusion (`∂p/∂t = c_v∇²p`) — no inertia, no advective velocity field with a nonlinear `(u·∇)u` term — so turbulent structure is not even representable. In this regime, inertialess mechanics forbids net locomotion from internal cycling alone; directed motion needs anchored external traction (adhesion clutch on ECM).

**Quantitative prediction.** No inertial range, no energy cascade, no vortex shedding at any achievable actuation. Any claimed "internal-flow-driven" net cell translation/rotation predicts COM/orientation drift `= 0` in a closed momentum budget; measurable net motion appears *only* when an external traction term (FA→ECM) is present.

**Competing / null hypothesis (H5₀).** Sufficiently strong internal stirring can create turbulence and/or net "swimming/rolling" of the cell.

**Falsification.** Enforce a closed-momentum diagnostic in the engine: internal-force-only configurations must give zero net COM/angular displacement to numerical precision. Any run showing net locomotion without an external traction channel is either a momentum-leak bug or refutes H5. Introducing a Navier–Stokes solver and finding `Re` never exceeds `~10⁻²` at physiological forcing further supports H5.

---

## H6 — Field-imposed filament MISALIGNMENT degrades organized retrograde flow and thus DIRECTED migration (hinder-by-disruption) — the defensible actuation direction

**Statement.** While fields cannot *drive* a cell (H3–H5), a field that partially disorients cytoskeletal filaments (e.g. via bead torques or low-salt electro-orientation) *degrades* the organized, polarized retrograde actin flow that migration depends on — a disruption effect that slows directed motion without supplying any propulsive force. This is the physically defensible (`d1`) actuation direction.

**Rationale.** Directed migration requires a coherently polarized actin network feeding retrograde flow and clutch engagement. Randomizing filament orientation reduces the net vectorial component of retrograde flow (an order-parameter reduction), lowering the traction that can be rectified into forward motion. This is subtractive (destroy order) not additive (inject momentum), so it evades the momentum/scallop constraints.

**Quantitative prediction.** Migration speed scales roughly with the nematic order parameter of the retrograde-flow field: imposing a disorientation that drops `S` from `~0.6 → ~0.3` reduces directed speed by `≈ 40–60%`, with net traction anisotropy falling commensurately. Effect is monotone in disruption strength and saturates when `S → 0` (undirected random walk), never reversing sign into enhanced migration.

**Competing / null hypothesis (H6₀).** Field-imposed alignment *enhances* directed migration (fields as a positive steering actuator), or filament orientation has no effect on migration speed.

**Falsification.** Engine: impose a controlled misalignment torque field, measure `S(t)` of retrograde flow and directed speed. If speed rises with disruption, or is insensitive to `S`, H6 fails. Wet-lab analogue: bead-torque disorientation should slow, never accelerate, persistent migration.

---

## H7 — Cell migration is driven by actin polymerization + actomyosin + adhesion-clutch traction; cytosol flow is a consequence/coupling medium, except in bleb/amoeboid/plant regimes

**Statement.** The primary motor of mesenchymal migration is protrusive actin polymerization plus actomyosin contraction transmitted to ECM through adhesion clutches. Cytosol flow is a *downstream consequence* of that mechanics and a coupling medium (it redistributes pressure and transports solutes) — it is not the migration driver, with the explicit exceptions of bleb-based/amoeboid motility and (non-animal) plant streaming, where hydraulic pressure/flow becomes a primary actor.

**Rationale.** Established mechanism: barbed-end polymerization pushes the membrane, NMII contracts the network, and catch-slip clutches (native pN-validated in the engine) rectify retrograde flow into ECM traction. The Biot fluid carries `~1 Pa` spatial pressure deviation at drained steady state — too small to be the mechanical driver — but it is essential for bleb inflation (fluid pressure detaching membrane) and hydraulic long-range coupling, which is where the exception regimes live.

**Quantitative prediction.** In the mesenchymal regime, ablating polymerization or clutch traction abolishes directed migration (speed `→ 0`), while turning the fluid coupling off (`FSI → OFF`) changes steady migration speed by `< few %`. In a bleb/amoeboid regime, the ordering flips: fluid pressure ablation abolishes bleb-driven motion.

**Competing / null hypothesis (H7₀).** Cytosol flow is the primary migration driver in the mesenchymal regime (flow-first migration everywhere).

**Falsification.** Knockout matrix in the engine: {polymerization, actomyosin, clutch, fluid-coupling} on/off × {mesenchymal, bleb}. If disabling fluid coupling kills mesenchymal migration, or if disabling polymerization/clutch does not, H7 fails for that regime.

---

## H8 — In the CFD engine, imposing (or myosin-driving) a directed u-field modulates SF growth ONLY where transport is limiting — an in-silico mirror of the wet-lab

**Statement.** Coupling a directed velocity field `u` (either externally imposed or generated by myosin contraction acting as the Biot source) to the transported monomer field will measurably change stress-fiber growth *only* in configurations that are transport-limited (high `Pe`, high consumption), and will leave reaction-limited FAs essentially unchanged. This is the direct in-silico test of H1+H2.

**Rationale.** Combines H1 (FAs are reaction-limited) and H2 (advection matters only at `Pe ≳ 1`). The engine solves `∂c/∂t = D_G∇²c − u·∇c + src − sink` with FA barbed-end sinks; growth reads local `c`. Where the sink barely dents `c` (reaction-limited), advecting `c` around does nothing; where the sink outpaces diffusion (transport-limited protrusion), advection changes delivered `c` and hence growth.

**Quantitative prediction.** Growth response to imposed `u`: at a resting FA (`Pe ≲ 0.1`), `Δgrowth < 10%` across a decade of `u`; in a `Pe ≈ 3` protrusion, `Δgrowth ≈ 2–4×`-scaled by the H2 enhancement. The response map of `Δgrowth` vs `Pe` collapses onto a single threshold curve crossing at `Pe ≈ 1`.

**Competing / null hypothesis (H8₀).** Imposed `u` modulates SF growth uniformly regardless of local transport regime (flow is a universal growth knob).

**Falsification.** Run the `u`-sweep across FA-type and protrusion-type sites; if reaction-limited FAs show `>25%` growth modulation, or if the `Δgrowth`–`Pe` curve has no `Pe ≈ 1` threshold, H8 fails. Cross-check: the in-silico threshold must reproduce the wet-lab observation (advective sensitivity only in fast protrusions).

---

## H9 — The cytosol fluid's mechanical role is rate-dependence and poroelastic τ_p transients, NOT a shifted static baseline

**Statement.** Turning the poroelastic fluid ON (`FSI-ON`) changes the *transient/rate-dependent* response and the pore-pressure structure, but does not materially change the drained static force balance or resting cell shape. "Fluid-first" is an architectural claim (an always-co-solved medium), not a claim that the fluid dominates statics.

**Rationale.** At drained steady state the fluid's spatial pressure deviation is `~1 Pa` (verified this session), negligible against `K_drained ≈ 300 Pa` and cortical tension. The fluid's real signature is rate-dependence: fast loading (`t_ramp ≪ τ_p`) traps undrained pressure (`p_max ≈ 8 Pa` in the CPU smoke), slow loading drains it (`p_max ≈ 0.03 Pa`). `τ_p ~ 1–3 s` from `D_p ≈ 40–60 µm²/s`.

**Quantitative prediction.** `FSI-ON` vs `FSI-OFF` static cortex `γ`, resting radius, and static AFM force differ by `< few %`; but the fast/slow AFM force ratio and the poroelastic relaxation time differ measurably (`τ_p ≈ 1–3 s`, undrained stiffening at fast rates). The discredited "emergent 65.9 Pa·s viscosity" is dropped (dimensional estimate `340–1100 Pa·s`, `5–17×` over) — no changed resting baseline is claimed.

**Competing / null hypothesis (H9₀).** The fluid substantially changes the static equilibrium (resting shape/tension) — fluid dominates statics.

**Falsification.** Native `--from-resting` full-compartment run, `FSI-ON` vs `FSI-OFF`, matched inner/outer steps to remove relaxation confounds. A `>5%` change in static `γ`/shape refutes H9 (fluid is statically dominant); no measurable `τ_p` rate-dependence would instead refute the "rate-dependence is the deliverable" half.

---

## H10 — A mass-conserving finite monomer pool shifts growth control to global concentration + thymosin buffering; the current fixed-scalar reservoir over-predicts sustained localized growth

**Statement.** Replacing the fixed-scalar `G_actin = 20 µM` (infinite reservoir) with a mass-conserving transported pool + thymosin-β4 buffering will show that sustained, spatially localized high-demand polymerization is limited by the *global* free-monomer budget and buffer exchange, not by local delivery — and that the current infinite-reservoir model over-predicts sustained growth under concentrated demand.

**Rationale.** The engine currently treats G-actin as an unlimited scalar (a known fidelity gap): polymerization draws from it without depletion. Physiologically the total actin pool is conserved and `>80%` of monomer is thymosin-sequestered, so free `[G-actin]` is buffered and finite. Under concentrated demand the reservoir depletes, and it is the buffer-release kinetics and global concentration that cap growth.

**Quantitative prediction.** With conservation on, a burst of localized polymerization draws down the free pool with a global relaxation set by buffer exchange; sustained growth rate saturates at a level below the infinite-reservoir prediction (over-prediction `≳ 20–50%` for high-demand configs). Growth becomes sensitive to *total* pool concentration (order-1 elasticity) while remaining insensitive to local flux (H1).

**Competing / null hypothesis (H10₀).** A finite conserved pool behaves identically to the fixed scalar under physiological demand (reservoir never limits) — the fidelity gap is inconsequential.

**Falsification.** Implement mass-conservation + buffer; compare sustained localized growth to the fixed-scalar baseline. If they agree within a few % under physiological demand, H10₀ holds and the fixed-scalar is adequate. If global-concentration elasticity is `< 0.3`, the "control-by-global-pool" claim fails.

---

## H11 — AC fields above the membrane-charging frequency penetrate the interior but yield zero net monomer drift; electroporation precedes DC electrophoretic Pe = 1

**Statement.** Above the membrane charging frequency (`~1.6 MHz`) an AC field penetrates the interior (the membrane capacitance shorts out), but because it oscillates the time-averaged monomer drift is zero. For DC, the interior field needed for `Pe = 1` is only reached after electroporation. So there is no field protocol that produces a sustained directed interior monomer drift below the damage threshold.

**Rationale.** Membrane shielding is frequency-dependent: DC/low-frequency is blocked (`10³–10⁴×`, H3), but above `~1.6 MHz` the bilayer capacitor passes the field. However, an oscillatory interior field gives `⟨u⟩ = 0` over a cycle — no net electrophoretic transport. DC penetration sufficient for `Pe ≈ 1` requires transmembrane `~1 V` (order `kV/cm` external), which is the electroporation regime.

**Quantitative prediction.** AC (`>1.6 MHz`): interior instantaneous `Pe` can approach 1 but net time-averaged `Pe ≈ 0` (`|⟨drift⟩| < 1%` of instantaneous). DC: crossing `Pe = 1` requires `E_interior` at the electroporation threshold (`~1 V` transmembrane), so any "useful" DC drift co-occurs with membrane permeabilization.

**Competing / null hypothesis (H11₀).** A sub-damage AC or DC protocol produces a sustained directed net interior monomer drift with `Pe ≳ 1`.

**Falsification.** Frequency-resolved shielding model coupled to the advection field: compute time-averaged net drift vs frequency and vs DC amplitude, with the electroporation threshold overlaid. Any sub-threshold protocol delivering net `Pe ≳ 1` refutes H11.

---

## H12 — Myosin contraction as the Biot fluid source is the only mechanistically-correct internal actuator that can locally reach Pe ≳ 1; external fields cannot substitute

**Statement.** The physically legitimate way to raise intracellular advective flux to the regime where it matters (`Pe ≳ 1`, H2/H8) is myosin-driven cytoplasmic streaming entering the CFD as the poroelastic source (`∂p/∂t = c_v∇²p + source(∇·v_solid)`), not an external field. Only actomyosin contraction can locally drive `u` fast enough, and only in specific geometries (fast/extended protrusions).

**Rationale.** External fields fail on shielding (H3, H11), require embedded particles that only stir (H4), and cannot self-propel (H5). Actomyosin contraction, in contrast, is the native pressure source: contracting the solid skeleton pressurizes the pore fluid and drives Darcy flow `u = −(k/µ)∇p`. This is the same primitive that generates traction (H7) and is the sanctioned "actuator" for the in-silico transport test (H8).

**Quantitative prediction.** Myosin contraction sourcing the Biot field produces local `u` reaching `Pe ≈ 1–3` only in extended high-flow geometries (matching the H2 lamellipodium regime), and `Pe ≲ 0.1` in the bulk/at FAs. External-field actuation produces interior `Pe ≪ 1` everywhere (H3/H11). Thus the myosin-source channel is the *only* one whose `Pe` map ever crosses 1.

**Competing / null hypothesis (H12₀).** An external field (E or B) can substitute for myosin as the actuator and drive interior `Pe ≳ 1`, or myosin streaming never reaches `Pe ≳ 1` in any regime.

**Falsification.** In the engine, compare the interior `Pe` maps produced by (a) myosin contraction as the Biot source and (b) any external-field protocol. If an external field reaches `Pe ≳ 1` interior, H12 fails; if myosin streaming never reaches `Pe ≳ 1` in any protrusion geometry, the "only actuator that matters" claim fails and advection is irrelevant everywhere (which would also collapse H8's high-`Pe` branch).

---

## Summary map

| ID | Claim | Key quantitative prediction | Refuted if |
|---|---|---|---|
| H1 | SF growth is reaction-limited | flux-elasticity `≲ 0.1–0.2` | elasticity `> 0.3` |
| H2 | Advection matters only at `Pe ≳ 1` | crossover `v* ≈ 0.5 µm/s`; protrusion `Pe ≈ 3` | FA `Pe≈0.1` shows `>25%` enhancement |
| H3 | DC field shielded; galvanotaxis = signaling | interior drift `Pe ≈ 4e-4–4e-2` | directional response scales with interior field, survives PI3K/PTEN block |
| H4 | B needs beads; rotation = local stir, no net flux/translation | net flux `<1%`, COM drift `≈0` | sustained net flux or COM drift without external traction |
| H5 | No turbulence; no self-translation | `Re ≪ 2000`; closed-momentum COM `=0` | net locomotion from internal forces only |
| H6 | Field misalignment degrades directed migration | `S: 0.6→0.3` ⇒ speed `−40–60%` | disruption speeds up / no `S`-dependence |
| H7 | Migration = polymerization+actomyosin+clutch; flow is medium | fluid-off Δspeed `<few %` (mesenchymal) | fluid-off kills mesenchymal migration |
| H8 | Imposed `u` modulates SF growth only where transport-limited | FA Δgrowth `<10%`; protrusion `2–4×` | reaction-limited FA modulates `>25%` |
| H9 | Fluid = rate-dependence, not static baseline | FSI on/off static Δ `<few %`; `τ_p ≈ 1–3 s` | `>5%` static shift |
| H10 | Conserved monomer pool ⇒ global-concentration control | over-prediction `≳ 20–50%` under high demand | finite pool ≈ fixed scalar physiologically |
| H11 | AC penetrates but net drift 0; DC needs electroporation | AC `⟨Pe⟩≈0`; DC `Pe=1` at `~1 V` transmembrane | sub-damage protocol gives net `Pe ≳ 1` |
| H12 | Myosin-as-Biot-source is the only interior `Pe≳1` actuator | myosin `Pe≈1–3` in protrusions; fields `Pe≪1` | external field reaches interior `Pe ≳ 1` |

**Global falsifier / honesty note.** All magnitude verdicts (SF `10–30 nN`, nucleus flatten, `γ`) inherit the `~530×` active density floor and are gated on the I0 magic-number ratification; these hypotheses are about *mechanism, direction, and scaling* (elasticities, Péclet thresholds, shielding decades, momentum budgets), which are robust to the density floor. Every in-silico test must be run at NATIVE full-compartment scale (`--from-resting`, all compartments ON) before its verdict is reported.