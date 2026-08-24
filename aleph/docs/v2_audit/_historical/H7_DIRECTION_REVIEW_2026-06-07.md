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

# H.7 Direction Review — emergent γ (GATE-B) + lamellipodium spreading geometries

**Date**: 2026-06-07 · **Reviewer**: critical review subagent · **Scope**: direction
soundness of the H.7 reboot (full physiological MCF7 cell → emergent cortical tension γ via
3 channels + two single-cell lamellipodium spreading geometries). **No code modified, no git
run.** This is a direction critique, not a sign-off.

Files audited: `CLAUDE.md`; `aleph/docs/briefs/H7_LAMELLIPODIUM_SPHERICAL_INTEGRATION.md`;
`aleph/configs/mcf7_baseline.yaml`; `aleph/cell/manifest.py`;
`aleph/cell/cell.py` (`build_cortex_full_simulation`, `Cell.build`);
`aleph/cortex/cortical_tension.py`; `aleph/cell/equilibration.py`;
`aleph/cell/lamellipodium.py` (+ `lamellipodium_basal_ring.py`,
`lamellipodium_polarized_patch.py`); `aleph/scripts/h7_gate_b.py`;
`aleph/scripts/h7_spreading_compare.py`; `aleph/cortex/enclosed_volume.py`.

---

## 0. Verdict

**SOUND in framing, with two substantive confounds and one concrete bug that must be fixed
or explicitly bounded before the GATE-B number or the spreading comparison is treated as
authoritative.** The 3-channel separated-γ design is the *correct* way to test KU-3.5 and is a
real methodological advance over the prior single-total approach — it directly implements the
physiological-baseline rule (measure on a turgor-pressurised cell, myosin as modulator) and
keeps the passive term honest. But:

1. The **passive (turgor) channel is band-circular by construction** — `turgor_dP0 = 133 Pa`
   was *derived from* the band-centre γ = 0.50 mN/m, and `gamma_passive = ΔP·R/2` recovers
   exactly 0.50 mN/m back. The passive channel currently cannot independently confirm or
   refute the band; it is an identity. This is the single most important thing to surface.
2. The **band itself is the wrong target for an adherent MCF7** (rounded/de-adhered
   HeLa/L929/neutrophil proxy; no MCF7 datum) — already known in MEMORY, but the GATE-B
   harness still prints it as *the* overlay with no adherent-MCF7 alternative.
3. The **`CappingUpdater` hardcodes the +ŷ membrane normal** (`lamellipodium.py:894`,
   `cos_theta = tangent[1]`), which is physically wrong for the `basal_ring` geometry whose
   tangents lie in the (x,y) plane — a real (if step-1-tolerable) bug.

Detail and citations below.

---

## 1. Is the GATE-B framing correct?

### 1.1 The 3-separate-channels design — YES, this is right

Reporting active (soft method-of-planes), rigid (M-SHAKE Lagrange), and passive (turgor
Young-Laplace) **separately and never folding turgor into a total** is the correct test of
KU-3.5, and the code enforces it as a hard contract: `gamma_structural = nansum(soft, rigid)`
explicitly excludes the passive term (`cortical_tension.py:552-556`), the docstring states the
B3 separation rule (`cortical_tension.py:25-33, 495-499`), and `h7_gate_b.py` prints "turgor
NOT folded into a total" (`h7_gate_b.py:126`). This is the mechanistically honest decomposition
the literature supports: cortical tension is myosin-generated network stress (Chugh 2017 [1];
Chugh 2018 [4]; Kelkar 2020 [3]), distinct from the osmotic pressure balance that pre-tensions
the shell. Keeping them separate is exactly what prevents the prior "in-band total mistaken for
active closure" failure mode.

The B4 bond-type denylist (excluding `integrin_ligand` / `fa_actin_clutch*` from the soft
channel, `cortical_tension.py:108-147, 346-352`) is also correct and well-reasoned: the FA load
path is a closed, stable set, the actin set is open/growing, so a denylist is the robust choice
(`cortical_tension.py:78-99`). The `|û·n̂|` absolute-value fix in the method-of-planes
(`cortical_tension.py:255-257`) is a genuine correctness fix (√N vs N cancellation) and is
documented with a synthetic-shell cross-check to Irving-Kirkwood <0.5%.

### 1.2 CONFOUND A — the passive channel is band-circular (HIGH severity)

`mcf7_baseline.yaml:42-51` derives the turgor from the band:

> `Young-Laplace Π₀ = 2·γ/R with band-centre γ = 0.50 mN/m at R = 7.5 µm → 133 Pa`

and `_gamma_passive` then computes `gamma = dP·R/2 = 133·7.5e-6/2 = 0.499e-3 N/m ≈ 0.50 mN/m`
(`cortical_tension.py:473-474`). So `gamma_passive` is **algebraically pinned to the band
centre** — it is `(2γ_band/R)·R/2 = γ_band` by construction. The passive channel therefore
*cannot* be read as evidence that the cell sits in-band; it is an identity that returns the
input. This directly contradicts the GATE-B claim that γ is "emergent at the natural operating
point." The passive number is a *设定* (a setpoint), not a measurement.

This matters because the realistic interphase turgor anchor in the same codebase is
**Stewart 2011 ~40 Pa interphase → ~400 Pa metaphase** (`enclosed_volume.py:25, 71-82,
182-184, 216`), which the manifest explicitly overrides ("Stewart 2011 interphase 40 Pa is the
lower alternative, superseded for the MCF7 baseline by this PI call", `mcf7_baseline.yaml:50-51`).
40 Pa → `gamma_passive = 0.15 mN/m`; 133 Pa → 0.50 mN/m. The choice of 133 over 40 is a
~3.3× swing in the passive channel that comes *from the band the run is supposed to be
testing against*. **Recommendation**: either (i) anchor turgor to an independent osmotic
measurement (Stewart 2011, or a measured MCF7 osmolarity), and let the band comparison be a
genuine prediction, or (ii) keep 133 Pa but state in every GATE-B report that the passive
channel is a definition, not a measurement, and that only `gamma_soft + gamma_rigid` is a
prediction. The current code does the honest thing for the *structural* total but the report
still prints `gamma_passive` next to the band as if it were an independent line.

### 1.3 CONFOUND B — the band is a rounded/de-adhered non-MCF7 proxy (HIGH severity, known)

MEMORY already records this; the review confirms it from primary literature. The
`CORTICAL_TENSION_BAND_N_PER_M = [0.35, 0.65] mN/m` is Salbreux/Charras/Paluch 2012
(`cortical_tension.py:182-185`), a review band built from **rounded, non-adherent** cells
(mitotic HeLa, L929, neutrophils). Primary numbers: passive human neutrophils ~30–37 pN/µm =
**0.030–0.037 mN/m** (Hochmuth 2000 [7]; Evans 1989 — *an order of magnitude below the band*);
S180 active tension γ ≈ 0.4 nN/µm = **0.4 mN/m** (Smeets 2018 [5], a fit not a direct
adherent-cortex measurement). An **adherent, spread** epithelial cell is not a liquid-drop
cortical-tension object at all — its surface mechanics are dominated by stress fibres and
substrate traction (Kumar 2019 [6]: "at high spread areas the cortex is thinner... sensitive
to myosin/traction"), so a single scalar "cortical tension" band is the wrong observable for
the spread state the H.7 cell is being built into (the manifest scope is `active_spreading`,
`mcf7_baseline.yaml:24`). **Recommendation**: the GATE-B target for the *spread/adhered*
operating point should be a traction/strain-energy or an apical-cortex-tension-on-a-rounded-cell
measurement, not the spread-cell-incompatible droplet band. If the only available MCF7-relevant
γ is on a rounded cell, then GATE-B should measure on a *rounded* (low-adhesion) build and the
spread build should be scored by traction, not γ. Surfacing this to PI is warranted because it
is a gate-contract (target-observable) question, not an inline tweak.

### 1.4 The demo-scale soft ≈ 0 finding — a real risk for any smoke conclusion

The `h7_gate_b.py` smoke path forces `demo_mode=True` at `n_filaments` ≈ 160
(`h7_gate_b.py:57-58`). Chugh 2017 [1] and Truong Quang 2021 [7] both show cortical tension is
set by **network architecture** — filament length, mesh size, and how far myosin minifilaments
penetrate the mesh — with a tension *maximum at intermediate filament length* (a non-monotone
dependence, [1]). A ×40-coarsened, demo-scale mesh has the wrong mesh size and connectivity, so
a near-zero `gamma_soft` at smoke scale is expected and uninformative — it is a pipeline check,
which the script correctly labels (`h7_gate_b.py:14-18, 146-148`). The risk is only if anyone
reads a smoke `gamma_soft ≈ 0` as physics. The labelling is good; keep it. But note the
corollary: even the **full-scale** soft channel will be sensitive to the connected-mesh build
parameters (`cm_z_struct`, `cm_bundle_mult`, `cm_reach`, `cell.py:719-722`), so γ_soft is only
meaningful once the mesh architecture itself is validated against Chugh's thickness/length data,
not just the bead count.

### 1.5 Is 2000 warmup steps enough for myosin contraction? — LIKELY NO; check the timescale

The cortex remodels and builds tension on **timescales of tens of seconds** via myosin
contraction + turnover (Kelkar 2020 [3]: "remodelling on timescales of tens of seconds";
Chugh 2018 [4]). The GATE-B `--warmup` default is 300 and the doc mentions ~2000
(`h7_gate_b.py:157`; `h7_spreading_compare` warmup 120, `h7_spreading_compare.py:149`). The
integrator dt is `dt_cfl = cfl_safety·τ_min` with `τ_min = min(τ_xl, τ_stretch, τ_bend)`
(`cortex.py:450-457`) — i.e. a *bending/stretch* CFL timestep, typically ~ns–sub-µs. Over
2000 steps that is microseconds of physical time, **many orders of magnitude short of the
tens-of-seconds myosin-contraction timescale**. So the "settled operating point" the warmup
produces is a *force-balance / overlap-drained* state (which is what `equilibration.py` is
designed for — it drains construction overlaps and thermalises, `equilibration.py:28-37`), **not
a state where myosin contraction has developed cortical tension**. Reading γ_soft after a
µs-scale warmup measures the construction stress, not the steady contractile tension.
**Recommendation**: separate two questions — (a) is the cell mechanically equilibrated (overlaps
drained, BAOAB stable)? (the current warmup answers this); (b) has myosin contraction reached
steady state? (needs a run length set by the *myosin/turnover* timescale, not the CFL timestep).
Report the physical time elapsed (`n_steps · dt_used`, both available at `h7_gate_b.py:78`) next
to γ so the reader can see whether contraction had time to develop. If full physical-time
equilibration is infeasible at the CFL dt, that is itself a finding to surface (it motivates the
constrained/rigid backbone, which removes the stiff stretch CFL — see §4).

---

## 2. Physiological-baseline completeness — what is still missing, ranked by γ impact

The manifest is genuinely thorough for the *mechanical* baseline (cortex+myosin+xlinks,
enclosed-volume turgor, nucleus, cytoplasm viscosity, membrane-surface; `mcf7_baseline.yaml:34-63`),
and the loader *raises* if any of these is off (`manifest.py:132-169`) — good enforcement of the
physiological-baseline rule. Missing structures, ranked by likely impact on the **cortical
tension γ** (not on general cell mechanics):

| Rank | Missing structure | Status in repo | Likely Δγ impact | Evidence |
|---|---|---|---|---|
| **1** | **Actin turnover** (cofilin severing / re-annealing) | declared, OFF (`mcf7_baseline.yaml:104-105`; `turnover.py` resolvable, `manifest.py:234-238`) | **HIGH** | Turnover sets filament length, and tension is *maximized at intermediate length* (Chugh 2017 [1]); without turnover the mesh length distribution is frozen at construction → γ_soft is pinned to an arbitrary architecture, not the regulated one. Kelkar 2020 [3]: turnover is half of what makes the cortex contractile on the relevant timescale. |
| **2** | **Stress fibres / contractile actin bundles** | not modelled (lamellipodium is branched Arp2/3 only) | **HIGH for the spread state** | In a *spread/adherent* cell, surface mechanics are dominated by stress fibres + traction, not the thin apical cortex (Kumar 2019 [6]; Laly 2021 [5-IF] shows keratin↔stress-fibre coupling). The GATE-B cell is `active_spreading` scope but has no stress-fibre load path, so the spread-state γ is structurally incomplete. |
| **3** | **Adherens / tight junctions, contractile belt** | not modelled (single cell) | MEDIUM (single-cell N/A; HIGH for Layer-2 hand-off) | Irrelevant for an isolated cell, but the γ is being fed to a *collective* Layer-2 model where junctional/belt tension is a large fraction of the effective surface tension. Flag at the seam, not here. |
| **4** | **Microtubules** | not modelled | **LOW in this model, despite biology** | Biologically MT depolymerization *raises* cortical tension — but **indirectly**, via GEF-H1 release → RhoA → ROCK → myosin-LC20 phosphorylation (Chang 2008 [MT-1]; Kolodney 1995 [MT-7]; Wang 2023 [MT-8]; Rape 2011 [MT-2]). This model has **no Rho/ROCK signaling layer**, so adding explicit MT struts would only contribute a small *compressive/load-bearing* mechanical term (Brangwynne/Pogoda: MT networks *soften* in compression), NOT the tension-raising biology. So microtubules are low-priority for γ *in this purely-mechanical model* — and adding them naively (as struts) could even bias γ the wrong way. This is a subtle but important point: the biology that makes MTs matter for tension is signaling the model does not have. |
| **5** | **Intermediate filaments (keratin; vimentin only on EMT)** | not modelled | **LOW for γ, MEDIUM for cytoplasm** | Guo 2013 [IF-10] is explicit: "VIFs contribute *little to cortical stiffness* but are critical for *intracellular* mechanics" (doubles cytoplasmic shear modulus to ~10 Pa). MCF7 natively expresses *keratin*, not vimentin, and gains vimentin only on EMT induction (Sivagurunathan 2022 [IF-7]). So IFs mostly belong in the *cytoplasm* compartment (already lumped as effective viscosity 65.9 Pa·s, `mcf7_baseline.yaml:36-38`), not in the γ budget. Low priority for GATE-B. |

**Top-line**: the two structures most likely to change the GATE-B γ are **actin turnover (#1)**
and **stress fibres (#2)** — both already-known gaps, both directly evidenced by Chugh 2017 and
Kumar 2019. Microtubules and IFs are *low* priority for γ in this model specifically because the
model lacks the signaling (Rho) layer through which they act on tension. That is a defensible
scoping call, but it should be stated explicitly so "we added all compartments" is not
over-claimed: the cell is mechanically complete but **signaling-free**, and tension regulation in
real cells is substantially signaling-driven.

---

## 3. The two lamellipodium geometries — is A/A0(isotropic-rim vs polarized-patch) the right
discriminator?

### 3.1 Partly. A/A0 distinguishes the *shapes* but not against MCF7 the way intended

`h7_spreading_compare.py` builds `basal_ring` (isotropic) vs `polarized_patch` (migrating) and
compares footprint A/A0 + a polarization index (`_polarization_index`, circular resultant length,
`h7_spreading_compare.py:72-79`). The metric is well-chosen — circular resultant is branch-cut
robust, unlike a max-min span — and isotropic-ring → ~0, polarized → ~1 is the correct signature.

But A/A0 alone is a **weak discriminator against MCF7 data**, for two reasons:

1. **The geometry is a preset, not an emergent outcome.** Whether a cell spreads isotropically or
   migrates polarized is governed by emergent Rac-Rho antagonism + ECM feedback, with cells
   switching among persistent-polarized / random / oscillatory modes depending on substrate and
   signaling (Holmes 2016 [LM-5]; Sadhu 2023 [LM-1]; Krause 2014 [LM-4]). MCF7 itself does
   *both*: a symmetric circular lamellipodium during spreading and an IGF-I-induced *polarized*
   leading edge during chemotaxis (Mañes 1999 [LM-10]). So presetting `geometry=basal_ring`
   vs `polarized_patch` and asking "which matches MCF7" is mis-posed — the answer is
   "both, depending on state." The right framing is: which geometry matches *the specific PI
   assay condition* (a sub-confluent MCF7 spreading on col-I/glass → most likely **isotropic
   circular spreading** early, Mañes 1999, so `basal_ring` is the default), and the polarized
   case is a separate migration assay, not a competing hypothesis for the same data.

2. **A/A0 conflates protrusion with adhesion-limited spreading.** Footprint area in a real
   spreading cell is set by the protrusion–adhesion–contraction balance, and is strongly
   spread-area-/traction-dependent (Kumar 2019 [6]; Holmes 2016 [LM-5] couples lamellipodium
   growth to integrin engagement). The model's A/A0 here is computed from *lamellipodium actin
   bead positions only* (`_basal_actin_xy` filters `lamel*actin` types,
   `h7_spreading_compare.py:45-57`) — i.e. it measures how far the dendritic array advanced, not
   the membrane footprint, and with the **membrane reaction force OFF** (step-1) there is no
   load opposing advance, so A/A0 is closer to "free barbed-end reach" than "spreading area."

**What would actually distinguish them against MCF7 data**: (i) spreading *kinetics* A(t) shape
— MCF7 isotropic spreading follows a characteristic fast-then-saturating curve (the PI's
`A/A0 = a + b/R + c/R²` form lives in the Layer-2 line); a polarized patch gives a directional,
non-saturating extension. (ii) The **traction/adhesion distribution** (ring of FA at the leading
edge for isotropic, front-loaded for polarized) — which the FA clutch already produces and is a
stronger, more mechanistic discriminator than hull area. **Recommendation**: score the geometries
by A(t) *kinetics* + FA/traction azimuthal distribution, with the membrane load ON, against a
*specified* MCF7 spreading assay — not by membrane-OFF lamellipodial-bead hull area, and not as a
binary "which is MCF7."

### 3.2 BUG — `CappingUpdater` hardcodes the +ŷ membrane normal (MEDIUM severity, real)

`CappingUpdater.act` computes the load angle as `cos_theta = float(tangent[1])  # tangent·ŷ`
(`lamellipodium.py:894`) and `sin_theta = √(1 − cos²θ)`, feeding the Bell-Evans cap rate
`k_cap = k_cap0·exp(−F·δ_cap·sinθ/kT)` (`lamellipodium.py:896-897`). This assumes the membrane
normal is **+ŷ** — correct for the flat-plane reconstitution geometry, but **wrong for
`basal_ring`**, whose mother tangents are outward-radial in the *(x,y)* plane
(`lamellipodium_basal_ring.py:25-26, 80-85`) and for `polarized_patch` (forward along p̂). For a
basal-ring tangent lying in (x,y), `tangent[1]` is just its y-component — an arbitrary projection
that has nothing to do with the actual local membrane normal — so the per-barbed-end capping rate
is computed against the wrong angle. The membrane *reaction force* is OFF in step-1
(`h7_spreading_compare.py` header; `mcf7_baseline.yaml:99-101` membrane_load OFF), which limits
the damage, but **capping is always on** (it is an Arp2/3 D1 updater, not gated by membrane), so
the dendritic-density / abortive-branching observable is already biased for the non-flat
geometries even in step-1. This contradicts the brief's claim that only the *geometry frame*
changes while the mechanism is "physically unchanged" (`H7...md:135`). **Recommendation**: make
the capping membrane normal per-WAVE (read from `tangent_of` / a stored per-WAVE normal) before
the basal_ring/polarized_patch A/A0 numbers are trusted; add a Sanity-Gate sign/sense test that
the cap-angle term uses the *local* edge normal. This is the same per-WAVE-normal generalization
the brief already plans for `membrane.py:633-790` (`H7...md:152-158, 291-295`) — capping needs it
too, and the brief's geometry section does not call out the CappingUpdater.

### 3.3 Membrane-OFF step-1 is otherwise a reasonable staging choice

Running geometry/advance dry (load-free) first, with a Sanity-Gate sign/sense test on the
standalone path before generalizing the membrane contact test, is exactly the staged order the
brief recommends (`H7...md:297-307`) and is sound — *except* that A/A0 measured load-free is not
a spreading-area observable (§3.1.2) and the capping bug (§3.2) is live even with membrane off.

---

## 4. Method risks

- **M-SHAKE rigid backbone shunting active force (KNOWN, real).** GATE-B runs `constrained=True`
  (`h7_gate_b.py:69`), zeroing the harmonic backbone force (`k=0 if constrained`, `cell.py:991-992`)
  and carrying backbone tension as a Lagrange multiplier recovered as `T = λ·r₀/Δt`
  (`cortical_tension.py:437`). The rigid channel is necessary (a soft-bond-only sum under-reports
  the backbone ~200×, `cortical_tension.py:16-18`). The risk the task names — a *rigid* segment
  locally absorbing active force instead of transmitting it — is real: a perfectly rigid bond
  cannot store the strain that a real semiflexible filament would, so myosin force applied across
  a rigid segment appears instantly as λ rather than as a propagating tension wave. Whether this
  biases the *isotropic average* γ is plausibly small (the method-of-planes sums |projected
  tension| over all crossing bonds, and λ does carry the constraint force), but it should be
  validated: run the SAME settled state **both** constrained and unconstrained (soft backbone)
  and check `gamma_rigid(constrained) ≈ gamma_soft_backbone(unconstrained)` on the backbone bonds
  alone. If they disagree beyond ~few %, the rigid shunt is biasing γ. This cross-check is cheap
  and is the decisive test; I did not find it in `h7_gate_b.py`.

- **Mesoscale ×40 coarse-graining and γ (real, partly unavoidable).** The ×40 filament-scale CG
  is the one sanctioned coarse-graining (CLAUDE.md), but γ depends on *mesh size and filament
  length* (Chugh 2017 [1], non-monotone tension-vs-length). A ×40 mesh has √40 ≈ 6.3× coarser
  spacing; the connected-mesh build tries to recover physiological connectivity (z≈3.3, giant
  99%, per MEMORY) but the *absolute* γ_soft from a coarsened mesh is not guaranteed to match the
  fine-grained value even after connectivity is fixed, because the tension-maximum-at-intermediate-
  length physics is length-scale-sensitive. This is the reason MEMORY already routes *magnitude*
  to a fine-grained line and keeps Layer-2 FORM-only. GATE-B should likewise treat γ *magnitude*
  as scale-bridged, and report γ as "ratio to a self-consistent reference" or with an explicit
  CG-bias caveat, not as an absolute mN/m to be matched to the (already-wrong) band.

- **CFL / equilibration (see §1.5).** The warmup drains overlaps and thermalises but does not
  reach the myosin-contraction timescale. The constrained run *removes* the stiff stretch CFL
  (backbone is rigid), which is the main lever to reach longer physical time per step — so the
  constrained channel is also the one most able to actually develop contraction. Report physical
  time elapsed alongside γ.

- **FA-adhesion equilibration prelude (sound).** The softstart→BAOAB prelude
  (`equilibration.py`; `manifest.py:280-283`; `h7_gate_b.py:64-74`) is the correct fix for the
  construction-time substrate↔cortex spatial disjointness (`equilibration.py:7-26`) — a raw
  FA-adhered free run trips the int32 image guard. One risk: the prelude settles the cell onto
  the substrate under softstart (capped-displacement steepest descent), and the *settled contact
  footprint* — which sets the FA cap and (for basal_ring) the lamella ring radius
  (`lamellipodium_basal_ring.py:53-61`) — therefore depends on the prelude budget. A too-short
  prelude gives a too-small footprint → wrong ring radius → wrong A0. Verify the footprint /
  clutch count has *converged* in the prelude (the diagnostics dict exists,
  `cell.py:1593, 1605`) before reading A/A0 or γ.

---

## 5. Top-5 prioritized recommendations

1. **De-circularize the passive channel (or label it a setpoint).** `gamma_passive` is currently
   `≡ γ_band` because `turgor_dP0 = 2γ_band/R` (§1.2). Either anchor turgor to an *independent*
   osmotic datum (Stewart 2011 ~40 Pa interphase, already in `enclosed_volume.py`, or a measured
   MCF7 osmolarity) so the band comparison becomes a real prediction, or have every GATE-B report
   state that the passive line is a definition and only `gamma_soft + gamma_rigid` is a
   prediction. *Rationale*: without this, GATE-B's headline "emergent γ at the operating point"
   is partly an identity, and a reader could mistake the recovered 0.50 mN/m for a result.

2. **Fix the GATE-B target observable for the SPREAD state (surface to PI — gate contract).**
   The [0.35,0.65] mN/m band is a rounded/de-adhered non-MCF7 droplet band (Hochmuth 2000 [7]
   neutrophils ~0.03 mN/m; Smeets 2018 [5] S180 ~0.4 mN/m fit); an adherent spread MCF7's surface
   mechanics are stress-fibre/traction-dominated (Kumar 2019 [6]). Decide: measure γ on a
   *rounded/low-adhesion* build (band-appropriate) and score the *spread* build by
   traction/strain-energy instead. *Rationale*: matching a spread cell to a rounded-cell band is
   a category error; this is a contract-level target question, not an inline tweak.

3. **Fix the `CappingUpdater` +ŷ-hardcoded normal before trusting basal_ring/polarized_patch
   A/A0.** Make `cos_theta` use the per-WAVE local edge normal (from `tangent_of`), not
   `tangent[1]` (`lamellipodium.py:894`), and add a Sanity-Gate sign/sense test. *Rationale*: it
   is a concrete, live bug (capping runs even with membrane load OFF) that biases the dendritic
   density observable for exactly the two non-flat geometries the H.7 spreading comparison is
   built on; the brief's geometry plan omits it.

4. **Separate "mechanically equilibrated" from "myosin contraction developed," and report
   physical time.** The µs-scale warmup (300–2000 steps at a bending/stretch CFL dt) cannot reach
   the tens-of-seconds myosin/turnover timescale (Kelkar 2020 [3]); γ_soft after it measures
   construction stress, not steady contractile tension. Print `n_steps·dt_used` next to γ, and use
   the constrained (rigid-backbone) channel — which removes the stiff stretch CFL — as the one
   that can actually advance physical time toward contraction steady state. *Rationale*: the
   "emergent operating point" must be a contracted state, not a drained-overlap state.

5. **Turn on actin turnover (#1) and add a stress-fibre load path (#2) before the full-scale γ
   is authoritative; cross-check the rigid channel against an unconstrained soft-backbone run.**
   Turnover sets the filament-length distribution that the tension-maximum physics (Chugh 2017 [1])
   depends on; stress fibres dominate spread-cell surface mechanics (Kumar 2019 [6]). And run the
   same settled state constrained vs unconstrained to confirm `gamma_rigid ≈ soft-backbone γ`
   (§4) so the M-SHAKE shunt is bounded. *Rationale*: these are the two compartments most likely
   to move γ, plus the one cheap test that bounds the dominant method risk. (Microtubules and IFs
   are correctly *low* priority for γ in this signaling-free mechanical model — state that
   scoping call explicitly rather than implying full completeness.)

---

## References (verified via Consensus/PubMed; existence confirmed — no invented citations)

Cortical tension / cortex architecture:
- [1] Chugh et al. 2017, *Nat Cell Biol* — Actin cortex architecture regulates cell surface
  tension (tension max at intermediate filament length).
  https://consensus.app/papers/details/3edb1475cd1e5469b21c32e88c94c184/
- [3] Kelkar, Bohec, Charras 2020, *Curr Opin Cell Biol* — Mechanics of the cellular actin
  cortex (remodelling on tens-of-seconds timescale).
  https://consensus.app/papers/details/5a92bb77c8d3586a9866a2b4d3412a61/
- [4] Chugh & Paluch 2018, *J Cell Sci* — The actin cortex at a glance.
  https://consensus.app/papers/details/312b4ce369a45835b9c514e0d57f70ff/
- [5] Smeets et al. 2018, *Biophys J* — Cortical elasticity & active tension (S180 γ≈0.4 nN/µm).
  https://consensus.app/papers/details/35f5a2f1e74a5cdebc7442bdda572a68/
- [6] Kumar et al. 2019, *BBA Mol Cell Res* — Spread area & traction determine cortex thickness.
  https://consensus.app/papers/details/a96453e6ff675e8bafccd977075dd9d8/
- [7] Hochmuth 2000, *J Biomech* — Micropipette aspiration (neutrophil cortical tension ~30 pN/µm).
  https://consensus.app/papers/details/23ddf245a08c50daa5f9d1ef0a98198b/
- Truong Quang et al. 2021, *Nat Commun* — Myosin penetration depth regulates surface mechanics.
  https://consensus.app/papers/details/2966ecdddab55c39a1be3e3a196c2a33/

Microtubule → contractility (indirect, via GEF-H1/RhoA — no direct cortical-tension mechanism):
- [MT-1] Chang et al. 2008, *Mol Biol Cell* — GEF-H1 couples nocodazole-induced MT disassembly to
  contractility via RhoA. https://consensus.app/papers/details/2ccab94e302757c5933e83777af8195b/
- [MT-2] Rape et al. 2011, *J Cell Sci* — MT depolymerization increases traction (two pathways).
  https://consensus.app/papers/details/bd473f56666a517c803764f5e54aa778/
- [MT-7] Kolodney & Elson 1995, *PNAS* — MT disruption → MLC phosphorylation → contraction.
  https://consensus.app/papers/details/9a13994125eb5a7ab97eaf5939e3e42f/
- [MT-8] Wang et al. 2023, *Curr Biol* — Mitotic cortical tension drug screen → Rho pathway;
  MT depolymerization raises tension via GEF-H1/RhoA.
  https://consensus.app/papers/details/e64b872691a25c70aeabbbeb5fab4131/

Intermediate filaments (little effect on cortical tension; matter for cytoplasm):
- [IF-10] Guo et al. 2013, *Biophys J* — VIFs contribute little to cortical stiffness, key to
  intracellular mechanics. https://consensus.app/papers/details/0ebcf79e229e5cd29307b207a0e71e56/
- [IF-7] Sivagurunathan et al. 2022, *Front Cell Dev Biol* — Vimentin in MCF7 (EMT-induced;
  native MCF7 is keratin). https://consensus.app/papers/details/d2e0ceb62ab65d12baa3ce4f068c3435/

Lamellipodium / spreading geometry:
- [LM-1] Sadhu et al. 2023, *J Cell Sci* — Minimal lamellipodia-based migration model.
  https://consensus.app/papers/details/24b83a6ffc9054a1b46dc4e61f85124e/
- [LM-4] Krause & Gautreau 2014, *Nat Rev Mol Cell Biol* — Lamellipodium dynamics / persistence.
  https://consensus.app/papers/details/ce29fa2add135da8bee3606c644bf888/
- [LM-5] Holmes et al. 2016, *PLoS Comput Biol* — Polarity-adhesion coupling explains
  persistent/random/oscillatory lamellipodial modes.
  https://consensus.app/papers/details/ef8de77ae9cb512296abb5c52462a8e9/
- [LM-10] Mañes et al. 1999, *EMBO J* — IGF-I induces polarized migration in MCF7 (vs symmetric
  spreading). https://consensus.app/papers/details/550352af96b25ea5a4d7a819bcff29a6/

In-repo provenance for the turgor/band/Stewart anchors: `enclosed_volume.py:25,71-82,182-184,216`
(Stewart 2011 40→400 Pa); `cortical_tension.py:182-185` (Salbreux/Charras/Paluch 2012 band);
`mcf7_baseline.yaml:42-51` (133 Pa band-derived turgor).
