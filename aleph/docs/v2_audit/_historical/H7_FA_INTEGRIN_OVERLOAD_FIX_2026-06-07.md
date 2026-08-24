---
archived_on: 2026-07-28
superseded_by: aleph/docs/v2_audit/AC_EXECUTION_PLAN_2026-07-25.md
reason: >
  Written BEFORE the 2026-07-25 PI reframe, i.e. for a different objective — forward prediction
  and magnitude matching, rather than inferring per-cell-type parameters with the gate on
  mechanical connectedness. Archived, not deleted: its measurements and reasoning stand as a
  record of what was true then. Nothing in it may be quoted as current state; STATE.md is that.
  Selected mechanically: pre-reframe AND cited by no live file (code, STATE.md, CLAUDE.md,
  cell_engine/, gate_contracts/, tests, Makefile). Citations from run outputs and from other
  pre-reframe documents were not treated as protective.
---

# H.7 GATE-B — Constrained FA build integrin–ligand overload diagnosis + fix

*Author: diagnostic subagent · 2026-06-07 · branch `h7/full-cell-integration`*
*Scope: analysis + proposed patch ONLY. No code edited, no git, no run.*

## TL;DR

The crash is **(a) over-stretched construction of the dynamic integrin–ligand
bond**, NOT a unit error and NOT (primarily) the rigid-backbone load. The
`integrin_ligand` Pereverzev catch bond is wired with **`r0 = 0`**
(`bridge/fa.py:660`, `cell/cell.py:1028`, from `configs/phase1_h4.yaml:56`
`integrin.r0: 0.0`). When `IntegrinBondUpdater` binds an integrin to its ligand
it does so for *any* integrin within `capture_radius_R_FA = 1.5 µm`, and with
`r0 = 0` the harmonic force is `F = k_int·r` — **already 50 pN (7.1× F_c) at the
nominal 50 nm contact gap, and up to 1500 pN (214× F_c) for a freshly-bound
far integrin.** This is the OPPOSITE of what the sibling S2 `fa_actin_clutch`
bond does — that one is deliberately placed force-free (`r0 = exact realised
separation`, `cell/cell.py:431–435`) precisely so the warm-up is stable. The
rigid backbone then prevents the cortex from relieving the stretch, so the
integrin is dragged to ~6.6 µm before the next batch tick, where
`pereverzev_k_off(F)` overflows its guard (`|F|/F_c ≈ 940`,
`validation/pereverzev.py:186`, called from `bridge/integrin_bonds.py:272`).

**Recommended fix:** make the dynamic integrin–ligand bond **force-free at the
moment of binding** — give it the same exact-separation `r0` treatment the S2
clutch already uses, applied per-bond inside `IntegrinBondUpdater.act()` when a
bond is *created* (candidate iii adapted to the dynamic case). This is the
physically honest fix: the catch bond should sit at its catch-peak operating
force (F* ≈ 7 pN) and ramp from there, not start at 50–1500 pN. See §3.

---

## 1. The integrin–ligand bond rest length vs the actual binding separation

### 1.1 Where r0 is set

The `integrin_ligand` bond type carries a **single global** `r0`:

* Config: `configs/phase1_h4.yaml:56` → `integrin.r0: 0.0` (commented
  "idealised attachment").
* Resolve: `bridge/fa.py:246` → `integrin_r0 = float(integrin["r0"])` = 0.0.
* Wired (cell path): `cell/cell.py:1027–1029`
  `bond.params[FA_BOND_INTEGRIN_LIGAND] = dict(k=p_fa.k_int_bare, r0=p_fa.integrin_r0)`.
* Wired (standalone path): `bridge/fa.py:660`
  `bond.params[BOND_TYPE_INTEGRIN] = dict(k=p.k_int_bare, r0=p.integrin_r0)`.

There are **no bonds at construction** (the comment at `cell/cell.py:1024–1026`
and `bridge/fa.py:577` is correct — the static `potential_energy ≈ 0` sanity
check holds). The bonds are created at RUNTIME by `IntegrinBondUpdater.act()`
(`bridge/integrin_bonds.py:311–321`) for every unbound integrin within
`capture_radius_R_FA` (`bridge/integrin_bonds.py:300–308`). Every such bond
inherits the type-level `r0 = 0`.

### 1.2 The actual binding separation

Contact-footprint geometry (the cell-integration path, `bridge/fa.py:456–501`):

* integrin sits `h_integrin_above_substrate` below its anchor cap bead:
  `pos_int.z = centres_z − h` (`fa.py:463`), `h = 5.0e-8 m = 50 nm`
  (`phase1_h4.yaml:43`).
* its ligand sits a further `h` below: `lig_z = centres_z − 2h` (`fa.py:492–493`).
* So the **nominal integrin↔ligand gap at construction = h = 50 nm.**

But the binding *eligibility* window is the full capture radius:
`capture_radius_R_FA = 1.5e-6 m = 1.5 µm` (`phase1_h4.yaml:51`), and binding
fires for any `d ≤ capture_radius` (`bridge/integrin_bonds.py:303`). Integrins
are scattered over a disk of radius `√(A/π)` around the anchor and the cortex
shell thermalises during the BAOAB drain, so realised binding separations span
**50 nm up to 1.5 µm.**

### 1.3 |F| / F_s at construction/binding (k_int_bare, geometry, F_s/F_c)

`k_int_bare = 1.0e-3 N/m` (`phase1_h4.yaml:55`), `F_s = 30 pN`, `F_c = 7 pN`
(`phase1_h4.yaml:82–84`). Pereverzev guard fires when `|F|/min(F_s,F_c) =
|F|/F_c > 700` (`validation/pereverzev.py:118,184–190`), i.e. `|F| > 700·F_c =
4.9 nN`.

With `r0 = 0`, `F = k_int·r`:

| Binding separation r | F = k_int·r | |F|/F_s | **|F|/F_c** | guard (700)? |
|---|---|---|---|---|
| 50 nm (nominal gap, h) | **50 pN** | 1.67 | **7.14** | no |
| 1.5 µm (capture edge, allowed bind) | **1500 pN** | 50.0 | **214** | no (yet) |
| **6.58 µm** (dragged under rigid load) | **6580 pN** | 219 | **≈ 940** | **YES → crash** |

Reverse-engineering the reported `|F|/F_c ≈ 940`: `F_crash = 940·F_c = 6.58 nN`
⇒ `r = F_crash/k_int = 6.58 µm` — that is **4.4× the capture radius.** A bond
cannot be *created* that far out (capture caps creation at 1.5 µm), so the
integrin must have been **dragged from a within-capture bind to 6.6 µm between
two batch ticks** before any `k_off` break could fire. The bond starts
over-stretched (≥ 7× F_c the instant it binds) and the rigid backbone removes
the only relief path.

**The catch bond is tuned for F* ≈ 6.99 pN** (analytic catch-peak,
`validation/pereverzev.py:25–31`), i.e. a clutch stretch of **~7 nm** at
`k_int = 1e-3 N/m`. The bond is being created at 50 nm (50 pN) to 1500 pN —
**7× to 214× past its design operating force** before dynamics even begins.

---

## 2. Which of (a)/(b)/(c) is it?

* **(a) over-stretched construction (placement) — YES, root cause.** `r0 = 0`
  on a bond that is *dynamically created* at separations of 50 nm–1.5 µm means
  the bond is born at 50–1500 pN, far above its F* ≈ 7 pN operating point. The
  S2 clutch sibling explicitly avoids exactly this by using
  `r0 = exact separation` (`cell/cell.py:431–435,463`) — the S1 integrin–ligand
  bond was left at the legacy `r0 = 0` idealisation. This is an internal
  inconsistency, not a physics requirement.
* **(b) unit error (N vs pN) — NO.** All forces are in N; `k_int_bare = 1e-3 N/m`
  with metre positions gives N. The 50 pN figure is dimensionally correct, not a
  1e12 slip. (The guard message even warns about the pN/N confusion, but that is
  not what is happening here — the magnitudes are genuinely physical-N and
  genuinely too large because of r0=0.)
* **(c) rigid-backbone load with no cushion — CONTRIBUTING, not primary.** The
  rigid backbone is why the *unconstrained* build survives: the soft cortex
  stretches, capping integrin displacement so even an r0=0 bond stays ≲ a few
  hundred pN and the next batch tick can break it before overflow. With the
  rigid backbone the integrin is dragged to 6.6 µm (940× F_c) before relief.
  But the rigid backbone only *exposes* the latent bug — if the bond were
  force-free at birth (F starts at 0, rises through F* ≈ 7 pN), the rigid build
  would not see a 940× spike on the very first loaded step. The cushion masks a
  bad initial condition; removing the cushion is not the disease.

Verdict: **(a) primary, (c) amplifier.** Fix (a) and the constrained build
inherits the same stability the unconstrained build has.

---

## 3. Recommended fix — force-free-at-binding integrin–ligand bond

### Candidate ranking

| Candidate | Verdict |
|---|---|
| **(i)/(iii) merged: per-bond exact-separation r0 at bind time** | **RECOMMENDED.** Honest, mirrors the S2 clutch, no magic number, respects Pereverzev. |
| (ii) gentle FA-load ramp during softstart | Palliative; does not fix the r0=0 birth force, only delays it. Rejected as primary. |
| (iii alone) seat integrins exactly at ligand plane | Only fixes the *nominal* 50 nm bond; a thermalised integrin that binds at 1.5 µm still births at 1500 pN. Insufficient on its own. |

The cleanest physically-honest fix unifies (i) and (iii): **a freshly-created
integrin–ligand catch bond must have its rest length equal to the exact
integrin↔ligand separation at the moment of binding**, so `F = 0` at birth and
the catch bond *earns* its load as the cell pulls, ramping naturally through its
F* ≈ 7 pN catch-peak. This is the identical principle the S2 clutch already
documents (`cell/cell.py:429–435`): "to keep the clutch FORCE-FREE at
construction … each clutch bond gets r0 = its EXACT initial separation … ½k·0²
= 0 J per clutch at construction."

### Physics honesty check (does it falsify the adhesion model?)

No. A nascent integrin–ligand bond forms at whatever separation the receptor and
ligand happen to be at — there is no physical reason the bond should be
pre-loaded to 50–1500 pN at the instant of catch. Setting `r0 = separation at
bind` is the *correct* zero-stress reference for a freshly nucleated bond; the
Pereverzev catch-slip kinetics (`k_off(F)`) then act on the load that develops
as the actomyosin/clutch pulls — which is exactly the regime Pereverzev 2005
(Biophys J 89:1446-54, DOI 10.1529/biophysj.105.062158) and the motor-clutch
operating point (Bangasser 2013, Biophys J 105:581-92, DOI
10.1016/j.bpj.2013.06.027) describe (O(pN), peaking near F* ≈ 7 pN here). The
KU-2.5 catch-peak validation gate is **unaffected**: that test
(`tests/validation/test_ku25_catch_peak.py`, per `integrin_bonds.py:82–89`)
holds an integrin at a *fixed prescribed* distance `F_target/k_int` and counts
breaks — it never relies on the runtime bind-time r0, so it keeps measuring
`k_off` against the oracle on the same force grid.

> Caveat — measurement-protocol consistency (`bridge/fa.py:62–75`): the S1
> docstring currently states "construction r0 = 0 means the gate readout is the
> *engaged* clutch force only, no baseline subtraction needed." Per-bond r0
> moves the zero-stress reference to the bind separation, so any **virial/per-bond
> force readout that assumes r0=0 must subtract the per-bond r0** (the force is
> `k·(|Δr|−r0)`, which the existing `F_mag` line already computes correctly).
> No measurement code reads "r raw"; the analytic `F_mag` at
> `integrin_bonds.py:271` is already `k·clip(r−r0,0,None)`, so it stays correct
> as long as `self.p.integrin_r0` is replaced by the per-bond r0 (see below).
> Update the docstring at `fa.py:62,72` to reflect the per-bond reference.

### The implementation requires a per-bond r0 (the dynamic analogue of S2)

The harmonic bond force in HOOMD is per **bond-type**, so a per-bond r0 means
one bond *type* per realised bind separation — the exact pattern the S2 clutch
uses (`cell/cell.py:431–435,509–517`, `fa_actin_clutch[_b{i}]`). For the
*dynamic* S1 bond this must be done inside `IntegrinBondUpdater.act()` at the
moment of binding, because bonds are created at runtime. Two equivalent
realisations:

**Option A (preferred, minimal, matches S2): per-bind r0-binned bond types,
registered up front.** Pre-register a family of `integrin_ligand_b{i}` bond
types covering r0 ∈ [0, capture_radius] (e.g. a fixed grid of `N_r0_bins`
bins), all with `k = k_int_bare`; at bind time pick the bin whose r0 is closest
to (or just below) the realised separation, so `F = k·(r − r0_bin) ≈ 0`. The
break loop computes `F_mag` per-bond from that bond's bin r0.

**Option B (exact, one type per live bond): not recommended** — unbounded type
growth as bonds churn; HOOMD type count would grow without bound over a long
run. The S2 clutch can do one-type-per-bond because those bonds are static;
dynamic S1 bonds churn, so a *binned* r0 (Option A) is the correct dynamic
analogue.

### Exact before/after sketch (Option A)

The change is localised to `IntegrinBondUpdater` and its construction wiring.
The key edits:

**(1) `bridge/integrin_bonds.py` — store a per-bond r0 and use it in `F_mag`.**

The break-loop currently uses the single `self.p.integrin_r0`:

```python
# bridge/integrin_bonds.py:271  (BEFORE)
F_mag = self.p.k_int_bare * np.clip(r - self.p.integrin_r0, 0.0, None)
```

```python
# (AFTER) per-bond r0 from the bond's TYPE (binned at bind time).
# r0_for_bond maps each surviving integrin_ligand bond -> its bin r0.
r0_bond = self._r0_for_bond(bt_int_bonds)   # shape (n_int_bonds,)
F_mag = self.p.k_int_bare * np.clip(r - r0_bond, 0.0, None)
```

where `bt_int_bonds` are the typeids of the current integrin bonds and
`_r0_for_bond` looks up the registered bin r0 (mirroring how `clutch_bin_r0`
is keyed by type name in `cell/cell.py:1035–1038`).

**(2) `bridge/integrin_bonds.py` — at bind (currently lines 311–321), assign the
new bond the r0-bin TYPE nearest the realised separation `d`:**

```python
# (AFTER) bind with a force-free r0 bin instead of the global integrin_ligand
d_bind = d[in_range][bind_mask]                      # realised separations
bin_typeids = self._nearest_r0_bin_typeid(d_bind)    # type per new bond, F≈0
# ... carry bin_typeids into the rebuilt bonds.typeid instead of a flat
#     `integrin_type` fill (replaces the constant fill at line 324).
```

**(3) Registration — register the `integrin_ligand_b{i}` bond-type family with
`k = k_int_bare` and the bin r0 grid.** In the cell path this slots in beside
the S2 clutch family at `cell/cell.py:1027–1029`:

```python
# cell/cell.py:1027 (BEFORE)
bond.params[FA_BOND_INTEGRIN_LIGAND] = dict(k=p_fa.k_int_bare, r0=p_fa.integrin_r0)
```

```python
# (AFTER) the canonical type stays (r0=0, used only when a bond happens to bind
# at ~0 separation), PLUS an r0-binned family spanning [0, capture_radius].
bond.params[FA_BOND_INTEGRIN_LIGAND] = dict(k=p_fa.k_int_bare, r0=0.0)
for name, r0_bin in integrin_ligand_r0_bins(
    n_bins=p_fa.n_integrin_r0_bins, max_sep=capture_radius
):
    bond.params[name] = dict(k=p_fa.k_int_bare, r0=float(r0_bin))
```

The bin types must also be added to `snap.bonds.types` in
`_extend_snapshot_with_fa` (alongside `FA_BOND_INTEGRIN_LIGAND` at
`cell/cell.py:498–499`) so the updater can swap typeids at runtime without a
state rebuild, and the updater's `act()` partition at `integrin_bonds.py:244–256`
must treat *all* `integrin_ligand_b{i}` typeids as integrin bonds (not just the
single `integrin_type`).

`n_integrin_r0_bins` is the only new knob; choose it so the worst-case residual
birth force `k_int · (capture_radius / n_bins)` stays well under F* ≈ 7 pN. With
`capture_radius = 1.5 µm`, `k_int = 1e-3 N/m`, requiring residual < F_c = 7 pN
gives `n_bins > k_int·capture/F_c = 1500 pN / 7 pN ≈ 215`; round up to a safe
**256 bins** (residual ≤ ~5.9 pN < F_c, and most binds at the 50 nm nominal gap
land in bin 0 ⇒ F ≈ 0). This is a *grid-resolution* parameter derived from the
F_c ceiling, not an empirical fit — it satisfies the Magic-Number Block exactly
as the S2 `n_bins` does. (Alternatively, snap the integrin onto the F*-stretch
plane at seeding so the nominal bind is at ~7 nm and r0 bin 0 covers it — but
binning is needed regardless for the thermalised far-binders.)

### Why this fixes the constrained build

A force-free birth means the integrin–ligand bond starts at F ≈ 0 and rises
through its F* ≈ 7 pN catch-peak as the rigid backbone loads it. The break-loop
then sees physically-bounded O(pN) forces every batch tick, k_off stays in
range, no overflow. The rigid backbone no longer matters because the bad initial
condition (50–1500 pN birth force) is gone — the same reason the unconstrained
build was stable, now provided by the bond's own rest length instead of the
soft cortex's compliance.

### Secondary hardening (defence-in-depth, not the fix)

Independently of the r0 fix, the `act()` break-loop should **break-before-it-
overflows**: clamp/saturate `k_off` (or treat `|F| > F_break_max` as a forced
unbind with probability 1) rather than letting `pereverzev_k_off` raise. A bond
stretched past the slip regime is physically already broken; raising a
`FloatingPointError` is a guard meant to catch *unit* bugs, not legitimate
high-load breaks. This makes the updater robust to any future over-stretch
without masking it. This is a guard-handling change, NOT gate-loosening: the
KU-2.5 oracle band is untouched; we only convert an overflow into the physically
correct "this bond breaks" outcome above the slip range.

---

## 4. Risks / things to verify before/after the patch

1. **Type-count growth.** Use the *binned* r0 family (Option A, fixed
   `n_bins`), never one-type-per-live-bond — dynamic churn would otherwise grow
   the HOOMD bond-type table without bound. 256 bins is fixed and small.
2. **Standalone H.4 path.** `build_h4_simulation` (`bridge/fa.py:651–742`) and
   its tests assume the single `integrin_ligand` type. Keep the canonical type
   registered (r0=0) and make the binned family additive so
   `tests/test_h4_topology.py` / `test_ku25_catch_peak.py` stay green; the
   KU-2.5 fixed-force harness binds at a prescribed distance and can keep using
   the flat type.
3. **`act()` bond partition.** `integrin_bonds.py:244–256` currently splits
   bonds by the single `integrin_type`; it must recognise the whole
   `integrin_ligand_b*` family as integrin bonds, else binned bonds leak into
   `other_bg` and never break.
4. **Docstring/measurement note.** Update `bridge/fa.py:62,72` — the "r0=0 ⇒ no
   baseline subtraction" claim no longer holds; the per-bond r0 is the
   zero-stress reference. The analytic `F_mag` already subtracts r0 correctly.
5. **Citation-integrity flag (surface to PI, unrelated to the patch).**
   `validation/pereverzev.py:85–87` and the H.4 brief cite KU-2.18 as
   *"Bangasser 2013 … Determinism and stochasticity during maturation of the
   zyxin focal adhesion"*. The real Bangasser 2013 at the cited locus (Biophys J
   105(3):581-92, DOI 10.1016/j.bpj.2013.06.027, PMID 23931306) is titled
   **"Determinants of maximal force transmission in a motor-clutch model of cell
   traction in a compliant microenvironment."** Volume/issue/pages/year match
   exactly — this is a wrong-title metadata drift, not a fabricated source. The
   Pereverzev 2005 anchor (Biophys J 89(3):1446-54, DOI
   10.1529/biophysj.105.062158, PMID 15951391) is verified correct. Recommend a
   SourceEvidence title fix (matches the 2026-06-02 metadata-drift audit
   pattern); does not affect the parameter values used.

## References (verified via PubMed, 2026-06-07)

- Pereverzev YV, Prezhdo OV, Forero M, Sokurenko EV, Thomas WE (2005). "The
  two-pathway model for the catch-slip transition in biological adhesion."
  *Biophys J* 89(3):1446-54. DOI 10.1529/biophysj.105.062158. PMID 15951391.
  — the `k_s·exp(F/F_s) + k_c·exp(−F/F_c)` form + the O(pN) operating range;
  catch-peak F* ≈ 7 pN with KU-2.18 params.
- Bangasser BL, Rosenfeld SS, Odde DJ (2013). "Determinants of maximal force
  transmission in a motor-clutch model of cell traction in a compliant
  microenvironment." *Biophys J* 105(3):581-92. DOI 10.1016/j.bpj.2013.06.027.
  PMID 23931306. — KU-2.18 illustrative clutch parameters; motor-clutch
  load-and-fail operating regime (O(pN)). [NOTE: codebase docstrings carry the
  wrong title for this DOI — see §4.5.]
