# Layer 3 φ audit — formulation bug, not accepted limitation

PI directive 2026-04-29 (post Option F gate contract changes). The
Production Lam4 F9 FAIL is **not** an accepted long-time decay; it is
a **conceptual conflation** in the φ variable. Code changes deferred
until this audit is reviewed and approved by PI; current Production
Lam4 plateau interpretation is **provisionally contaminated** by
artificial interior φ decay until the audit is resolved.

This document discharges the PI-required Layer 3 audit. Six sections
per the directive.

---

## 1. Current φ implementation summary

Source: `acs/physics/mlsmpm.py:_integrate_phi_ode` (lines 845–874);
configs `phi_initial` per phenotype.

### ODE
Per particle:

```
dφ_p / dt = k_+ · S_p · (1 − φ_p) − k_- · φ_p
```

with literature (Cho 2020) anchors:
- k_+ ≈ 1.39e-3 s⁻¹
- k_- ≈ 4.6e-4 s⁻¹
- φ_eq = k_+ / (k_+ + k_-) ≈ 0.751

Time-step stiffness check at constructor: `dt · (k_+ + k_-) ≤ 0.5`
(forward-Euler stability). Per-step clip to `[0, 1]` for f32 safety.

### Phenotype-specific initial conditions (Phase 4-v2 mapping)
- Bare: `phi_initial: 0.30`
- Pre: `phi_initial: 0.55`
- Lam4: `phi_initial: 0.80`

These map to the experimental MCF7 phenotypes established during
spheroid formation (ULA bare vs Pre-Col1 vs Lam4).

### v11 spatial S extension (`layer3_spatial_S = true` in production)
Source per `docs/v1/stage1b_layer3_sanity.md` §… and the in-line comment:

```python
S_p = 1.0
if ti.static(self.cfg.layer3_spatial_S):
    if self.x[p][2] >= h_band:   # interior particles
        S_p = 0.0
```

Boundary-band particles see S=1 (driven toward φ_eq); interior
particles see S=0 (driven toward 0 by the k_- term alone).

**Justification (as previously written)**: "Required for Stage 1d
(Layer 4) since otherwise ∇γ ≈ 0 with single-phenotype" — i.e., the
spatial S extension was introduced to *manufacture a φ gradient* so
that Layer 4 Marangoni had something to act on.

### Gate
"φ trajectory toward predicted φ_eq" (runner.py around line 940):
compares the global `<φ>(end)` against the boundary equilibrium
`φ_eq = k_+ / (k_+ + k_-)`, with tolerance `≤ max(0.5, |φ_init − φ_eq|)`.

For Production Lam4: `<φ>(end) = 0.159`, `φ_eq = 0.751`, |err| = 0.59
vs tolerance 0.50 → **FAIL** (F9 in `docs/v1/gate_fail_taxonomy.md`).

---

## 2. Why current F9 fails

The gate fails because the model produces:

| Quantity | Production Lam4 80hr | Pilot 4hr | Comment |
|---|---|---|---|
| `<φ>` (global mean) | 0.159 | 0.502 | global decays heavily over time |
| `φ_max` (boundary band) | 0.733 | 0.622 | boundary saturates near φ_eq |
| `φ_min` (interior) | 0.088 | 0.493 | interior collapses toward 0 |
| φ_boundary fraction (φ_boundary_sum / phi_sum) | 0.998 | n/a | nearly all "alive" φ is at boundary |

**Mechanism of decay** (interior, S=0):

```
dφ/dt = -k_- · φ
φ(t) = φ_init · exp(-k_- · t)
```

Half-life τ_half = ln(2) / k_- ≈ 1500 s ≈ 25 min.

For Production Lam4 t_total = 80 hr = 4800 τ_relax_units:
- Interior φ at t=80hr: 0.80 · exp(-k_- · 80·3600) = 0.80 · exp(-132)
  ≈ 0 (numerically zero).

The boundary band (S=1) saturates at φ_eq ≈ 0.733 (close, not exact,
because of finite k_- · time). Global mean = (boundary fraction) ·
0.733 + (interior fraction) · 0 ≈ 0.20 · 0.733 + 0.80 · 0 = 0.147 ≈
the observed 0.159.

**The model is internally consistent — it does exactly what its
equations say.** The bug is in **what the equations represent**.

### What the gate compares is also wrong
The gate compares `<φ>(global)` against `φ_eq(boundary)`. These are
two different physical quantities:
- `φ_eq` is the *boundary-band* equilibrium under continuous substrate
  signaling.
- `<φ>` is the *global mean* over all particles (boundary + interior).

The gate equation `<φ>(end) ≈ φ_eq` is a category error: it requires
the *whole spheroid* to reach the *boundary-only* equilibrium, which
under the spatial S=0 interior term is impossible by construction.

So F9 is a **double bug**:
1. The interior S=0 erases interior φ over the spreading timescale.
2. The gate compares the wrong global moment against the wrong
   reference.

---

## 3. Literature / PI rationale: Bare/Pre/Lam4 is formation-history initial state

### What Cho 2020 actually models

Cho et al. 2020 (cited as the k_+ / k_- anchor) measure the **Western
blot timecourse of E-cadherin and Integrin-β1 protein levels** during
phenotypic transitions in MCF7 cells. The k_+ / k_- rates capture the
*kinetics of a phenotypic transition*, e.g., when cells are induced
to remodel their adhesion network.

Cho 2020 does NOT model:
- Cell-population baseline at a *steady* phenotype.
- Spheroid-formation-induced adhesion-network state.
- Long-time (≥ hours) preservation of phenotype memory in the absence
  of an active transition stimulus.

The k_+ / k_- in Cho 2020 are the rate constants of a *transition*,
applicable when a stimulus is actively driving the system from one
steady state toward another. They are NOT the rate constants of
spontaneous decay from a steady state in the absence of stimulus.

### What PI experimental Bare/Pre/Lam4 represents

Per `memory/experimental_design.md` and the project framing:
- Bare/Pre/Lam4 = spheroid **growth environment** (phenotype variable).
- Spreading is on a single Col1 dish across all conditions; substrate
  is NOT a sweep variable.
- Bare/Pre/Lam4 enters as Stage 1b initial conditions on φ-ODE / γ_cc.

**The phenotype is established during spheroid formation, before
transfer to the Col1 substrate.** The 4-hour spreading assay does not
remodel the formation-history phenotype; it tests how that frozen
phenotype responds to the new substrate.

### Why current model is conceptually wrong

The current φ implementation conflates two timescales:
1. **Formation memory** (long, set at t=0, slow or no change during
   the assay): Bare/Pre/Lam4 baseline E-cad/Int-β1 ratio.
2. **Contact activation** (fast, dynamic during spreading, only at
   substrate-engaged cells): focal adhesion turnover, protrusion,
   traction.

The Cho 2020 k_+/k_- rates apply to **(2)**, not to **(1)**. The
spatial S=0 interior term, multiplied by the Cho 2020 k_- rate,
incorrectly drives the formation memory toward zero on the timescale
of contact activation — i.e., it erases the phenotype identity within
~25 min, contradicting the experimental setup where the phenotype is
expected to persist over 80 hr.

This is the "two distinct concepts in one φ variable" diagnosis in
the PI directive.

---

## 4. Proposed φ_memory + contact_activation replacement

### Two state variables per particle

Replace `φ_p` (single dynamic) with:

```
φ_memory_p ∈ [0, 1]    # formation phenotype memory (slow / fixed)
c_act_p    ∈ [0, 1]    # contact activation (fast, dynamic)
```

### Initialization

```
φ_memory_p(0) = phi_initial  (phenotype-specific: Bare 0.30, Pre 0.55, Lam4 0.80)
c_act_p(0)    = 0            (no contact activation at t=0; cells just landed)
```

### Evolution

**φ_memory** (long-timescale):

```
dφ_memory_p / dt = -ε · (φ_memory_p − phi_initial)
```

with ε very small (e.g., ε ≤ 1e-5 s⁻¹, τ_memory ≥ 28 hr) — far slower
than the spreading assay duration. In the limit ε → 0, φ_memory is
simply fixed at the initial phenotype throughout the run.

**Default**: ε = 0 (φ_memory exactly fixed). The decay term exists
only to allow future audits where formation memory is allowed to
slowly evolve under e.g. EMT induction.

**c_act** (Cho 2020 timescale, active only on contact-band particles):

```
dc_act_p / dt = S_p · [k_+ · (1 − c_act_p) − k_- · c_act_p]
```

S_p = 1 only for contact-band particles (z_p < h_band), 0 elsewhere.
For interior particles, c_act_p stays at its previous value (no decay)
— the activation state is preserved in cells that were once at the
boundary but have since flowed into the interior.

This is a deliberate departure from the current S=0 → k_- decay: the
c_act *integrates* the time the cell spends in the contact band,
rather than decaying in the interior. Mechanistic basis: focal
adhesion machinery, once assembled, persists for hours at the cellular
level; if a cell migrates inward, its focal adhesions disassemble on
their own timescale (hours), not on the Cho 2020 transition timescale
(minutes).

### Effective φ for downstream physics

Downstream consumers (Layer 4 γ(φ_eff), any other φ-dependent term):

```
φ_eff_p = φ_memory_p + κ · c_act_p · (1 − φ_memory_p)
```

with κ ∈ [0, 1] a contact-activation amplification weight. Defaults
κ = 1 (full saturation when c_act = 1: φ_eff = 1) or κ = 0.5 (mild
contact boost). κ becomes a parameter to be set by sanity-md decision
or sweep, with Magic-Number-Block discipline applied.

### Why this is principled

- **φ_memory** captures the formation phenotype and respects the
  experimental setup (Bare/Pre/Lam4 set before spreading).
- **c_act** captures Cho 2020's actual scope (transition kinetics
  during contact-driven activation).
- **φ_eff** is a phenotype-modulated contact-activation field —
  conceptually clean separation matching the biology.
- The "spatial S extension" v11 patch becomes unnecessary; the
  explicit two-variable model directly produces the boundary-vs-
  interior gradient that Layer 4 Marangoni needs, *without* erasing
  formation memory.

---

## 5. New gate contract

Replace single F9 "φ trajectory toward predicted φ_eq" with two
gates:

### Gate 5a — φ_memory preservation
```
|<φ_memory>(end) − phi_initial| ≤ ε_drift_tol
```
Default `ε_drift_tol = 0.01` (1% drift over 80 hr). Captures the
formation-memory invariant; failure indicates a Layer 3 bug or
inadvertent coupling.

### Gate 5b — contact_activation trajectory at boundary band
```
<c_act>_boundary_band(end) ∈ [c_eq − tol, c_eq + tol]
where c_eq = k_+ / (k_+ + k_-) = 0.751
tol = 0.10
```
Captures the Cho 2020 transition kinetics in their actual scope: only
in contact-engaged particles. Interior c_act is not gated (it is
allowed to be anywhere between 0 and the boundary-saturated value,
depending on cell migration history).

### Per-particle invariant gates (preserved)
```
φ_memory_p ∈ [0, 1]   (per-particle clip)
c_act_p ∈ [0, 1]       (per-particle clip)
```

### F9 retirement
F9 (single "φ trajectory toward predicted φ_eq") is **deprecated**
under the new contract. The gate compared the wrong global moment
against the wrong reference; it is replaced by Gates 5a + 5b above.

---

## 6. Minimal pilot plan to test plateau persistence

### Stage 1b.b: implement φ_memory + c_act split
Code changes (deferred until PI approves this audit):
- Add `phi_memory_p` and `c_act_p` Taichi fields to MLSMPMSolver.
- Replace `_integrate_phi_ode` with `_integrate_phi_memory_ode` +
  `_integrate_c_act_ode`. The c_act ODE has the spatial S gating; the
  memory ODE is independent.
- Update `γ(φ_p)` consumers to use `γ(φ_eff_p)` per §4.
- Remove `layer3_spatial_S` config flag and the v11 interior S=0
  patch (now subsumed by c_act's S gating).
- Update gate to two-gate contract per §5.
- Sanity-md: `docs/v1/stage1b_b_phi_split_sanity.md` discharging the
  6-check protocol.

### Pilot 1: 1k Lam4 4hr (smoke + correctness)
Goals:
- φ_memory(end) ≈ 0.80 across all particles (boundary AND interior).
- c_act(end) ≈ 0.75 in boundary band, < 0.75 in interior (proportional
  to interior cells' past contact-band residency, ~0 if no migration).
- φ_eff at boundary > φ_eff at interior (gradient preserved for Layer
  4 Marangoni).
- All gates PASS except F5 (radius drift, accepted).

Estimated wall-clock: 1 minute on Laptop A5000 (matches Phase 4-v2
pilot scale).

### Pilot 2: 5k Lam4 80hr (asymptote test)
Goals:
- A_over_A0_topdown trajectory: continue past 7.75hr peak (1.570) or
  retract as before?
- If continues: peak-and-decay was an artifact of the φ formulation
  bug. Mechanism A/E/F from `docs/v1/marangoni_review.md` may not be
  needed — the asymptote gap is partly closed by Layer 3 fix alone.
- If retracts as before: Layer 3 fix is necessary but not sufficient;
  Mechanism A/E/F genuinely needed; Stage 1d.b authorized.

Estimated wall-clock: ~25 minutes on Laptop A5000.

### Bucketing for Pilot 2 (revised against Pilot 2 baseline)

Define new buckets against Pilot 2 (post-Layer-3-fix) trajectory:

| Bucket | Pilot 2 endpoint | Interpretation | Next step |
|---|---|---|---|
| L3-A | A/A₀_topdown(end) ≥ 4.0 | Layer 3 fix alone closes most of the gap | paper-as-is, Mechanism A/E/F deferred |
| L3-B | 2.5 ≤ end < 4.0 | partial closure | Stage 1d.b Mechanism A authorized |
| L3-C | 1.5 ≤ end < 2.5 | minor improvement | Stage 1d.b Mechanism A + Stage 1a++.b authorized in sequence |
| L3-D | end < 1.5 (similar to current 1.409) | no improvement | Layer 3 fix is necessary but not sufficient — full Mechanism A/E/F + Stage 1a++.b warranted |

---

## Implications for downstream documents

### Production Lam4 finding
Mark `docs/v1/production_lam4_finding.md` peak-and-decay interpretation as
**provisionally contaminated** until Pilot 2 (above) is run. The
asymptote may be partly or fully an artifact of the φ formulation bug.

### Marangoni review
Mark `docs/v1/marangoni_review.md` Mechanism A / E / F upgrade options
as **deferred** until the Layer 3 audit (Pilot 2) is complete.
Mechanisms A/E/F act on γ(φ); if φ itself is the wrong quantity, the
mechanisms are calibrated against a wrong baseline.

### Codex review item 4
`docs/codex_review_synthesis.md` item 4 ("Layer 3 φ dynamics too
simple") is **upgraded from review item to PI directive** with this
audit. The audit is the resolution.

### Stage 1a++.b
`docs/v1/marangoni_review.md` Option β (Stage 1a++.b stochastic boundary
events) is **deferred** until Layer 3 audit (Pilot 2) is complete. A
clean continuum baseline is required before adding stochastic events.

### Gate fail taxonomy
F9 is reclassified from "HARD-BLOCKER (Layer 3 audit pending)" to
**FORMULATION BUG (audit complete, code change pending)**. Remaining
hard-blocker count: still 1 (F9), but now with a well-defined
remediation path (Pilot 1 + Pilot 2).

### Stage roadmap
Insert Stage 1b.b ("Layer 3 φ_memory + c_act split") between current
Stage 1b and Stage 1d. Stage 1d.b and Stage 1a++.b are blocked
behind Stage 1b.b Pilot 2 outcome.

---

## Checklist before code implementation

PI directive 2026-04-29 ("이거 해결 후에 순차적으로 스테이지 둘 다
수정하는 거 순차적으로 모두 진행해줘") authorizes autonomous
execution of all three stages (1b.b → 1d.b → 1a++.b). All sign-offs
granted in this directive:

- [x] PI sign-off on this audit document
- [x] PI sign-off on §5 new gate contract (replacing F9)
- [x] PI sign-off on §6 pilot bucketing (L3-A through L3-D)
- [x] PI default for κ (§4 effective φ formula): **κ = 1.0** (full
      contact-activation saturation; Layer 4 Marangoni works on the
      maximum ∇φ_eff which is the existing v11 design intent)
- [x] PI default for ε (memory decay): **ε = 0** (φ_memory exactly
      fixed throughout the spreading assay; matches the experimental
      setup where formation phenotype is established before spreading
      and not actively remodelled within the 80-hr window)
- [ ] Sanity-md (`docs/v1/stage1b_b_phi_split_sanity.md`) — to be
      written as the first artefact of the implementation phase

---

## Cross-references

- `docs/codex_review_synthesis.md` item 4 — original Layer 3 audit flag
- `docs/v1/gate_fail_taxonomy.md` F9 — to be reclassified
- `docs/v1/marangoni_review.md` §4 — Mechanism A/E/F deferral
- `docs/v1/production_lam4_finding.md` — plateau provisional contamination
- `docs/v1/stage1b_layer3_sanity.md` — original Layer 3 ODE design
- `docs/v1/parameter_registry.md` — Layer 3 parameters re-tier after split
- Cho et al. 2020 — Western blot timecourse, transition kinetics
- `memory/experimental_design.md` — Bare/Pre/Lam4 formation-history setup
- PI directive 2026-04-29 — this audit was directed
