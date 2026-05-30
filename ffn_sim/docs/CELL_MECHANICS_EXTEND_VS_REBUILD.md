# Cell mechanics — Extend vs Rebuild (membrane / nucleus / cytoplasm)

**Status**: DRAFT for PI review — 2026-05-30. Produced in a READ-ONLY debug/search
session; **no existing code modified**. Proposed code blocks below are design
artefacts, NOT wired into the build.

**Trigger (PI, 2026-05-30)**: realised that whole-cell mechanics has been
implicitly equated with **cortex mechanics**; plasma-membrane surface mechanics,
nucleus, and cytoplasm viscoelasticity are absent. With cortex now carrying the
entire mechanical weight, the question raised was: *"do we have to tear the whole
thing down (갈아엎다)?"*

**Verdict (this memo)**: **EXTEND, do not rebuild.** The runtime is a HOOMD
fine-grained particle/bond simulator; membrane, nucleus, and cytoplasm are
*additional particle types and custom forces*, not a different paradigm. Tearing
down a particle-based engine to add more particles is self-defeating. What is
genuinely new is (1) the **physics of each module** and (2) **re-validating** the
KU gates that currently attribute all whole-cell mechanics to cortex. Neither
requires discarding the existing runtime, BAOAB, Cell composition skeleton,
cortex, or ECM.

This doc has two parts:
- **Part A** — the extend-vs-rebuild decision, with the reusable-vs-new map.
- **Part B** — a PoC-level *design verification* that membrane-surface mechanics
  and a nucleus module attach cleanly through the existing additive pattern,
  with the exact attach points (file:line) they hook into.

---

## Part A — Decision: Extend, not Rebuild

### A.1 Why this is not a teardown

The architectural principle (CLAUDE.md) is that the **runtime is HOOMD
particle/bond dynamics** and a `Cell` is *"bookkeeping over particle groups +
force computes + Updaters"* ([cell.py:13-17](cell/cell.py#L13-L17)). That is
exactly the substrate you want for adding membrane/nucleus/cytoplasm.

Hard evidence that the codebase *adds subsystems cleanly* — `build_cortex_full_simulation`
already bolts on **seven** subsystems using one repeated **additive +
default-off** pattern, each *"bit-for-bit identical when off"*:

| Subsystem | Pattern entry | Attach mechanism |
|---|---|---|
| crosslinkers | [cell.py:589-607](cell/cell.py#L589-L607) | snapshot extend + tag range + Updater |
| cortical myosin | [cell.py:609-628](cell/cell.py#L609-L628) | snapshot extend + tag range + Updater |
| lamellipodium (H.5) | [cell.py:630-652](cell/cell.py#L630-L652) | snapshot extend + tag range + Updaters |
| FA / clutch (H.4) | [cell.py:654-705](cell/cell.py#L654-L705) | snapshot extend (appended LAST) + Updaters + pin |
| enclosed-volume pressure | [cell.py:886-894](cell/cell.py#L886-L894) | **custom force over tag range** |
| actin turnover | [cell.py:972-981](cell/cell.py#L972-L981) | bond-topology Updater |
| membrane *load* | [cell.py:1029-1044](cell/cell.py#L1029-L1044) | load Action feeding lamellipodium |

Membrane-surface mechanics, nucleus, and cytoplasm are the **8th, 9th, 10th**
entries in this same table. Part B verifies the attach contracts concretely.

### A.2 What "membrane" is today (fact-check)

The existing [membrane.py](cell/membrane.py) is **not** plasma-membrane surface
mechanics. Per its own docstring ([membrane.py:1-67](cell/membrane.py#L1-L67))
it is a **leading-edge tension reservoir**: a single mobile plane DOF supplying
the per-barbed-end reaction load `F = γ_mem · s_fil` that the lamellipodium's
Bell-Evans elongation/capping laws consume. It is diagnostic-only on the plane
position and **adds no particles** (Phase 1 forbade adding a membrane sheet).

So the PI's instinct is correct: **whole-cell membrane mechanics (area
elasticity + Helfrich bending + membrane–cortex composite tension) does not
exist yet.** What exists is (a) this leading-edge load, and (b) the **ERM
tether** ([erm.py](cortex/erm.py)) — a radial harmonic holding cortex beads at
`R_cell`, which is the *seed* of membrane–cortex adhesion.

### A.3 The two real issues (neither is "rebuild")

**Issue 1 — Roadmap gap, not architecture flaw.** Phase 1 was *deliberately*
cortex-first (Plan v2 sequence H.1 ECM → H.2 filament → H.3 cortex → H.5 → H.7).
Membrane-surface / nucleus / cytoplasm were never first-class roadmap items.
"Cortex is too heavy" is a **phasing artefact**, not lock-in: in a real cell,
membrane tension and cortex tension are a *coupled composite* — cortex is not
oversized, it is currently **alone**. The weighting self-corrects by **adding**
the co-players, not removing cortex. This is a PI roadmap-amendment decision,
not a code disaster.

**Issue 2 — Re-attribution / re-validation (the real scientific risk).** The
KU-3.5 cortical-tension and KU-3.1 rounding gates currently attribute **all**
whole-cell mechanics to cortex. Once a membrane surface and nucleus share the
mechanical load, those gates may be **mis-attributing** — the measured tension is
γ_total = γ_membrane + γ_cortex, and rounding resistance gains a nuclear term.
This is a **gate re-derivation** task (composite tension), not a teardown. It is
*especially* important for cancer cells, where nucleus + cytoplasm are
mechanically co-dominant (see `CELL_TYPE`/cancer discussion, session 2026-05-30).

### A.4 Reusable vs new (the actual change surface)

| Component | Grade | Note |
|---|---|---|
| HOOMD runtime, BAOAB integrator, `methods=[]` contract | **reuse as-is** | no change |
| `Cell` / `build_cortex_full_simulation` composition skeleton | **reuse + extend** | add 3 optional blocks, default-off |
| cortex (H.3), ECM (H.1), crosslinkers, myosin, FA, lamellipodium | **reuse as-is** | untouched |
| ERM tether | **reuse + promote** | becomes the membrane–cortex linker |
| Plasma-membrane **surface** mechanics (tension + bending + area elasticity) | **NEW additive module** | custom force over shell tags (+ optional sheet particles) |
| Nucleus (lamin-set stiffness, radial confinement, deformability) | **NEW additive module** | new particle group + confinement force |
| Cytoplasm viscoelasticity | **NEW (more invasive)** | background medium / per-type γ + memory; current model has water-drag only |
| Cell–cell cadherin junctions / multi-cell | **NEW + composition layer** | single-cell → multi-cell |
| KU-3.5 / KU-3.1 gates | **RE-VALIDATE** (code kept) | re-derive as composite |

Estimated reuse of the existing runtime + plumbing: **~85–90%.** A teardown would
*discard* the BAOAB freeze, the seven-subsystem additive harness, and the
validated cortex/ECM — and then **re-implement the same particle/bond runtime** to
hold the new particles. That is strictly more work for zero architectural gain.

### A.5 Recommendation

1. **Do not rebuild.** Adopt the additive-module path verified in Part B.
2. **Amend the roadmap** to make membrane-surface mechanics, nucleus, and
   cytoplasm first-class units (propose: H.8 membrane surface, H.9 nucleus,
   H.10 cytoplasm; cadherin/multi-cell later).
3. **Open a re-validation item**: re-derive KU-3.5 (tension) and KU-3.1
   (rounding) as **composite** (membrane + cortex + nucleus) gates. Until then,
   label current cortex-only results as "cortex-attributed", not "whole-cell".
4. **Promote the whole-cell mechanics altitude**: the model's claim should become
   "membrane–cortex composite + nucleus", not "cortex ≈ cell".

---

## Part B — PoC design verification: do membrane + nucleus attach cleanly?

**Question**: can plasma-membrane surface mechanics and a nucleus be added
through the existing additive pattern, default-off, without touching the frozen
integrator or the existing subsystems?

**Method**: trace the two reusable attach templates already in the code, then map
each new module onto them and flag friction. (Templates extracted read-only,
file:line cited.)

### B.1 Template 1 — custom force over a tag range (ERM / enclosed-volume)

Both [erm.py](cortex/erm.py) and [enclosed_volume.py](cortex/enclosed_volume.py)
implement the *"a custom force acting on a contiguous tag range, read by BAOAB
via net_force"* contract:

- Inherit **`hoomd.md.force.Custom`**; compute in `set_forces(timestep)`.
- Read `self._state.cpu_local_snapshot` → mask `(tag >= start) & (tag < end)` →
  compute per-bead force → write `self.cpu_local_force_arrays`
  (ERM [erm.py:201-234](cortex/erm.py#L201-L234);
  enclosed-volume [enclosed_volume.py:437-470](cortex/enclosed_volume.py#L437-L470)).
- **CFL gate at attach**: `τ = γ_b / k_eff`, require `dt ≤ cfl_safety_factor · τ`,
  else `RuntimeError`
  (ERM [erm.py:283-298](cortex/erm.py#L283-L298);
  enclosed-volume k_eff = `K_vol·S²/(N²·V0)`
  [enclosed_volume.py:419-435](cortex/enclosed_volume.py#L419-L435)).
- The force **participates in HOOMD `net_force`** that the BAOAB Updater reads
  each step (ERM docstring [erm.py:10-12](cortex/erm.py#L10-L12)) — no extra
  registration beyond `ig.forces.append(...)`.
- Attach functions `attach_erm_to_simulation` / `attach_enclosed_volume_to_simulation`
  take `(sim, resolved_params, tag_range, gamma_b, cfl_safety_factor)` and return
  the force object.

→ **Membrane surface tension/bending** and **nucleus radial confinement** are
*exactly* this shape (forces over a tag range). They reuse Template 1 verbatim.

### B.2 Template 2 — new particle subsystem + BAOAB contract

For modules that add *particles* (a fully-dynamical membrane sheet, or nucleus
beads), the contract from [baoab.py](integrator/baoab.py) +
[cell.py](cell/cell.py):

1. **Register the type pre-snapshot** (add to `snap.particles.types`, set
   `typeid`/positions/mass) — the snapshot-extension pattern used by every
   subsystem (e.g. lamellipodium [cell.py:643-652](cell/cell.py#L643-L652)).
2. **Supply a `gamma_map` entry** keyed by type name, finite-positive
   `γ_b = 6π η R`, **before** BAOAB attach
   ([cell.py:899-918](cell/cell.py#L899-L918)). BAOAB **asserts every present
   type is covered** and raises `RuntimeError` on a missing entry
   ([baoab.py:283-294](integrator/baoab.py#L283-L294)) — so a forgotten type
   fails loudly at attach, never silently.
3. **`methods=[]`** stays empty; BAOAB is the sole position integrator
   ([baoab.py:276-281](integrator/baoab.py#L276-L281)).
4. **Quasi-static elements** (a fixed membrane plane, an anchored nucleus
   centroid) use the **post-BAOAB position-reset Action** pattern — BAOAB has
   **no integration-group exclusion**, so a subsystem that must not move is held
   immobile by an Action appended *after* the BAOAB updater that resets its tags
   each step (`SubstrateLigandPin` [cell.py:177-237](cell/cell.py#L177-L237),
   instantiation [cell.py:1096-1115](cell/cell.py#L1096-L1115)).

### B.3 Verification verdict per module

**Membrane surface mechanics — CLEAN via Template 1 (+ optional Template 2).**
- Minimal mechanistic form = a **surface-tension + curvature custom force over
  the cortex shell tags `[0, n_cortex_actin)`** (γ_mem contributing to
  γ_total = γ_membrane + γ_cortex; Helfrich bending κ_m as the curvature term).
  This is the enclosed-volume template with a different `k_eff`. **No new
  particles, no integrator change.**
- Membrane–cortex adhesion = **promote ERM** ([erm.py](cortex/erm.py)) from
  "radial tether to center" to "spring to a membrane reference surface" — same
  custom-force shape.
- A **fully-dynamical lipid sheet** (separate membrane particles outside the
  cortex shell) is the heavier option — Template 2 (new `membrane_bead` type +
  gamma_map + bonds for area/bending elasticity). Also clean, but adds particles;
  matches the refinement [membrane.py:64-67](cell/membrane.py#L64-L67) already
  flagged as "a later refinement … would add particles".

**Nucleus — CLEAN via Template 2 + Template 1.**
- New `nucleus_bead` particle group near the cell centroid (Template 2: type
  registration + gamma_map entry).
- Mechanics via a **radial confinement / shell-stiffness custom force**
  (Template 1) parameterised by lamin-A/C-set stiffness `k_nuc`.
- If treated as a quasi-rigid body initially, anchor the centroid with the
  post-BAOAB pin (Template 2 step 4) — exactly the `SubstrateLigandPin` pattern.

### B.4 Proposed PoC skeletons (design artefacts — NOT wired)

Following the verified Template 1 contract verbatim:

```python
# PROPOSED ffn_sim/cell/membrane_surface.py  (NOT IMPLEMENTED — design only)
class MembraneSurfaceTension(hoomd.md.force.Custom):
    """γ_mem surface tension + κ_m Helfrich bending over the cortex shell tags.
    Mirrors EnclosedVolumePressure: read snapshot, mask shell tags, write force."""
    def __init__(self, p_mem, shell_tag_range):
        super().__init__()
        self.p, self.t0, self.t1 = p_mem, *shell_tag_range
    def set_forces(self, timestep):
        with self._state.cpu_local_snapshot as snap:
            tag = np.asarray(snap.particles.tag)
            pos = np.asarray(snap.particles.position)
        mask = (tag >= self.t0) & (tag < self.t1)
        # F_i = γ_mem-driven inward/tangential surface force on shell bead i
        #       + κ_m curvature term (Helfrich). [physics TBD — the module's real work]
        with self.cpu_local_force_arrays as f:
            f.force[:] = ...     # per-bead, zero outside mask
            f.potential_energy[:] = ...

def attach_membrane_surface(sim, p_mem, shell_tag_range, gamma_b, cfl_safety_factor=0.1):
    ig = sim.operations.integrator
    # CFL: tau = gamma_b / k_eff(gamma_mem); assert dt <= cfl_safety_factor*tau  (Template 1)
    force = MembraneSurfaceTension(p_mem, shell_tag_range)
    ig.forces.append(force)
    return force
```

```python
# PROPOSED ffn_sim/cell/nucleus.py  (NOT IMPLEMENTED — design only)
class NucleusConfinement(hoomd.md.force.Custom):
    """Radial harmonic shell stiffness (lamin-A/C-set k_nuc) on nucleus tags."""
    def __init__(self, p_nuc, nucleus_tag_range):
        super().__init__()
        self.p, self.t0, self.t1 = p_nuc, *nucleus_tag_range
    def set_forces(self, timestep):
        with self._state.cpu_local_snapshot as snap:
            tag = np.asarray(snap.particles.tag); pos = np.asarray(snap.particles.position)
        mask = (tag >= self.t0) & (tag < self.t1)
        d = pos - np.asarray(self.p.center); r = np.linalg.norm(d, axis=1)
        rhat = d / np.where(r > 0, r, 1.0)[:, None]
        F = (-self.p.k_nuc * (r - self.p.R_nuc))[:, None] * rhat
        F[~mask] = 0.0
        with self.cpu_local_force_arrays as f:
            f.force[:] = F
            f.potential_energy[:] = np.where(mask, 0.5*self.p.k_nuc*(r-self.p.R_nuc)**2, 0.0)
```

Hook points in `build_cortex_full_simulation` (additive + default-off, mirroring
the enclosed-volume block at [cell.py:886-894](cell/cell.py#L886-L894)):
- nucleus particles: extend the snapshot **before** `create_state_from_snapshot`
  ([cell.py:711](cell/cell.py#L711)), append `gamma_map["nucleus_bead"] = γ_nuc`
  in the gamma block ([cell.py:899-918](cell/cell.py#L899-L918)).
- forces: `if p_membrane_surface: attach_membrane_surface(...)` and
  `if p_nucleus: attach_nucleus_confinement(...)` after the integrator is set
  ([cell.py:875](cell/cell.py#L875)), parallel to enclosed-volume.

### B.5 Friction / risks surfaced by the verification

1. **CFL stiffness gate (magic-number / PI-gate item).** Membrane tension and
   nuclear stiffness are *stiff* springs. If `k_eff` makes `τ = γ_b/k_eff < dt`,
   the Template-1 gate raises. Precedent: `k_ERM` was **softened 1000×**
   (0.1 → 1e-4 N/m, PI-ratified 2026-05-26) rather than shrinking dt 33×
   ([phase1_h3.yaml:198-222](configs/phase1_h3.yaml#L198-L222)). Same decision
   will recur for membrane/nucleus stiffness — surface to PI, do **not** soften
   silently.
2. **No integration-group exclusion in BAOAB.** A rigid-ish nucleus or fixed
   membrane plane needs the post-BAOAB pin Action; there is no native "freeze
   these tags" hook ([cell.py:185-194](cell/cell.py#L185-L194)). Known, solved
   pattern — just must be remembered.
3. **gamma_map coverage is load-bearing.** Any new particle type *must* get a
   gamma entry or BAOAB halts at attach ([baoab.py:283-294](integrator/baoab.py#L283-L294)).
   This is a feature (fails loud), but it's the one step easy to forget.
4. **Re-validation, not plumbing, is the real cost.** Adding γ_mem changes what
   KU-3.5 measures (composite tension); adding the nucleus changes KU-3.1
   rounding resistance. The *attach* is mechanical; the *acceptance gates* must
   be re-derived (Part A.3 Issue 2).
5. **Cytoplasm viscoelasticity is the genuinely harder one.** Tension/nucleus
   are forces-over-tags (clean). A viscoelastic cytoplasm is a *background
   medium with memory* — beyond a per-type `γ`. It likely needs either a
   memory-kernel drag or filler particles, and touches the BAOAB drag model more
   deeply than Templates 1–2. Flag as the one module that may need integrator-
   adjacent design (still additive, but not "just another custom force").

### B.6 PoC verdict

**Both membrane-surface mechanics and the nucleus attach cleanly through the
existing additive pattern** — Template 1 (custom force over tag range) and
Template 2 (new particle subsystem + gamma_map + optional pin) — with **zero
changes to the frozen integrator or existing subsystems**, default-off so every
current run stays bit-for-bit identical. The teardown hypothesis is **rejected**:
the only genuinely new work is each module's physics and the gate re-validation,
both of which a rebuild would also require (plus the cost of rebuilding the
runtime). Cytoplasm viscoelasticity is the lone module that may need
integrator-adjacent design and should be scoped separately.

---

## Open items for PI

- [ ] Ratify **EXTEND** (this memo) vs rebuild.
- [ ] Roadmap amendment: H.8 membrane surface / H.9 nucleus / H.10 cytoplasm as
      first-class units; cadherin + multi-cell later.
- [ ] Re-validation item: re-derive KU-3.5 / KU-3.1 as composite gates;
      re-label current results as "cortex-attributed".
- [ ] Cancer-cell KU anchor table (MCF10A / MCF7 / MDA-MB-231 moduli, cytoplasm
      viscosity, cortical tension) — prerequisite for cell-type parameterisation
      (session 2026-05-30 cancer-mechanics discussion).
- [ ] Decide membrane fidelity: tension-only custom force (cheap, Template 1) vs
      full lipid sheet (Template 2, adds particles).
