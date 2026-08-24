# Active Cell engine — parallel session launch pack (2026-07-16)

Verified plan from the 16-agent parallelization workflow (`wf_d1755894-d6e`), adversarially revised.
Each `### Session` block below is a **self-contained boot prompt** — paste one into a fresh Claude Code
session. Base branch for everything new-engine: **`ac/new-engine`** (I1a analytic foundation landed, Biot
`alpha=1.0` ratified, I0-B1 closed).

## 0. Quick-start

**Concurrency = 5** (front-loaded). Launch these 5 sessions AT ONCE (all have conflict-free layer-0/1 work):

| Session | Track / branch | Increment(s) now | Worktree |
|---|---|---|---|
| **A** | fluid-spine `ac/fluid-spine` | I1a Warp solver (then I1b→I1c) | yes |
| **B** | nucleus `ac/nucleus` | I2 deformable-mesh nucleus (I7 later) | yes |
| **C** | excluded-volume `ac/solid-ev` | I2b steric WCA | yes |
| **D** | motor-nmii `ac/motor-nmii` | I3 head-resolved NMII | yes |
| **E** | emergence `ac/emergence` | I5-detector + synthetic oracles only | yes |

Gated (launch when the dep lands): **F** unified-weave `ac/weave` (I4, after D+A modules) · **G** adhesion
`ac/adhesion` (I6, after F) · **H** nucleus-loadpath (I7, continues session B after D+F) · **I** mt-coupling
`ac/mt-coupling` (I8 oracle-core after D; module after B+F) · production `ac/production` (I9, terminal).

**⭐ Single-GPU native reality (from the adversarial verify — internalize):** the native `--from-resting`
gates form a **strict data-dependent serial chain I1a→I1b→…→I9 on the ONE gbook A5000**. Adding GPUs does
NOT shorten it (I_n's cell literally contains I_{n−1}'s live state). There is no 2nd GPU, so even per-gate
**ensemble width** (I5's ~5-control ablation, I7's MCF7-vs-MDA, I8's ON/OFF, I9's seeds) also serializes.
The 5-way parallelism buys exactly one thing: **every increment's Warp kernel + closed-form oracles are
authored-and-CPU-green before its GPU slot opens, so the A5000 never idles.** It cannot compress the spine.

---

## 1. Coordination protocol (BINDING — every session reads this first)

1. **Runtime = NVIDIA Warp on CUDA GPU only (I0-A).** HOOMD is never imported/executed. The dev Mac cannot
   run Warp-CUDA at all — you author kernel SOURCE + run only **host-side numpy/scipy oracles** locally.
2. **Your worktree, your net-new `ac/<pkg>/` subpackage ONLY.** Do NOT edit any `ff/` file. Do NOT edit
   another track's files. Parallel worktrees stay conflict-free because each creates disjoint net-new files.
3. **Shared `ff/` edits are DEFERRED to the lead-owned serial native-integration step.** Where your increment
   must retire/delete/extend an `ff/` mechanism (e.g. I1a retiring the `6πηR→γ_solid` clock + scalar turgor +
   held-face FSI in `network_warp.py`; I3 deleting lumped `f_myo`/`myosin_kernel`; I6/I7 adding `cell_type.py`
   presets; any `hand_kmc.py`/`polymerization_warp.py` change), you do NOT touch `ff/` in your worktree —
   instead write the change as a documented **patch-note + adapter** inside your `ac/` package
   (`ac/<pkg>/INTEGRATION.md`: exact file, symbol, before→after, the double-count guard it satisfies). The
   lead applies these in spine order.
4. **Frozen interface contracts** (consume read-only; if not present yet, code against the signature below +
   a local test-double and flag it):
   - `aleph/components/fluid/domain.py :: Domain` — the live membrane-minus-nucleus domain on the conservative field grid.
     **fluid-spine OWNS it.** It exposes `set_nucleus_boundary(provider)` where `provider` yields the inner
     relative-no-flux mask; I1a ships a STATIC-sphere provider placeholder, the nucleus track ships a
     `DeformableNucleusMaskProvider` conforming to the same signature. Neither co-edits the other's file.
   - **Inner force-assembly loop** — each force primitive is `accumulate(state, out_force)` adding its
     per-node contribution. `MyosinForce`(I3, `ac/motor/`), `StericForce`(I2b, `ac/solid/`),
     `PressureCoupling`(I1a, `ac/fluid/`) each own their module; the lead-owned integrator sums them.
   - **Per-head hand/KMC device API** — `aleph/components/motor/hand.py` (I3 OWNS). Downstream tracks (I4/I6/I8) add
     presets in their OWN `ac/<pkg>/presets.py` that CONSTRUCT hands via the I3 API — never edit the KMC core.
5. **analytic-first** (the I1a pattern, already committed under `ac/fluid/`): write the closed-form numpy
   oracle(s) + pytest FIRST (pure numpy/scipy, ZERO Warp import — Warp-only-contract exempt), then the Warp
   kernel source. Your local green = the analytic gates. **You do NOT run native/CUDA gates** — you hand the
   lead a CPU-green module + a written native-gate spec.
6. **I0-B\<n\> params**: any value whose literature grounding is not closed → do NOT choose. Write it as
   `value: null` + `evidence_status: "GAP — PI"` in `ac/<pkg>/params_i0bN.yaml` (mirror `ac/fluid/params_i0b1.yaml`)
   and surface to PI. Provisional PI-authored values must be labeled with owner + expiry-trigger, never tuned
   to pass a gate.
7. **HARD rules still bind**: no magic numbers (Magic-Number Block), no gate-loosening, physiological baseline,
   validate at full native population, Sanity Gate before first execution, no param-tuning-to-outcome.
8. **Closeout**: commit on your `ac/<track>` branch (do NOT push, do NOT merge). Hand the lead: green analytic
   suite + native-gate spec + `INTEGRATION.md` patch-notes + I0-B\<n\> gaps. Answer in Korean; code/commits/docs
   by convention.

### 1.9 Visualization (BINDING — the professor's HARD viz rules)

Text-only npz/json/log hides the physics; the visual check catches outliers / wrong-sign / sampling
artefacts that pass/fail bands miss (PI 2026-05-21). **Every track visualizes at each milestone.**
- **WHEN**: your module's analytic suite goes green · a (lead-run) native gate lands · a sanity-gate sweep
  flags a non-trivial finding.
- **WHERE**: `aleph/outputs/ac/<track>/figs/*.png` + a `REPORT.md` in that dir with a `## Figures` section
  (each figure + a one-line caption). Keep ONE regeneration entry-point script per track
  (`aleph/scripts/ac_<track>_vis.py`); commit figures alongside data. The **fluid-spine pattern is already
  seeded** (`scripts/ac_fluid_vis.py` + `outputs/ac/fluid-spine/`) — copy its structure.
- **ANALYTIC figures** (Mac, matplotlib): for EACH analytic gate, plot the numeric/Warp result WITH the
  closed-form oracle or band OVERLAID — never the measurement alone.
- **INTEGRITY (hard)**: no axis truncation; no log/linear flip without a note; ALWAYS overlay the
  reference/band; annotate units (SI default); ensemble plots = per-realisation thin lines + ensemble-mean
  overlay.
- **NATIVE cell morphology** (the lead runs it on gbook, but YOUR native-gate spec must NAME the field/
  geometry to render): interactive 3-D CELL-MORPHOLOGY **HTML** — real geometry, peel/slab/cut views, NOT a
  matplotlib png (PI 2026-07-02); **full-res, no downsampling** (PI 2026-07-07); **browser-verified** via
  `browser_check.py` — a WebGL grep is meaningless, the render must be eyeballed (PI 2026-07-07/08).
- **ANOMALY**: a run that looks wrong → STOP, render the figures, surface to PI; never self-correct or wave
  it off as "explosion" (PI 2026-06-12).

Per-track figure spec (what each session must produce):

| Track | Analytic figures (Mac, oracle/band OVERLAID) | Native-render field (HTML, lead) |
|---|---|---|
| fluid-spine | Terzaghi U(T_v)+field vs oracle; Green spread ⟨r²⟩=2d·c_v·t; MMS eigenmode decay + moving-interval residual; Darcy flux linearity; total-actin conservation; FRAP | spatial pore-pressure + v_f/flux field |
| nucleus | bending vs 8πκ (FD-grad); volume conservation; lamin-knee modulus crossover; rupture onset; (I7 capstan) | deformable oblate nucleus mesh |
| solid-ev | WCA potential + FD-grad sign; zero-overlap→zero-force; compressed-pair resistance; k_EV grid-invariance | no-interpenetration / crowding |
| motor-nmii | Hill FV vs Hill-1938; Bell off-rate vs load; ensemble-stall emergence; per-head Newton residual | engaged-head map / minifilament geometry |
| emergence | S~0 null vs S=1 aligned; planted-bundle localization; rotational+label invariance; effect-size | (bundle-condensation overlay, deferred) |
| weave (I4) | branch-angle θ₀±σ_θ dist; weave_cell([CORTEX]) parity | full region set / dendritic net |
| adhesion (I6) | Pereverzev off-rate + F* peak; two-sided Newton; FA-growth fixed point | traction map + adhesion classes |
| mt-coupling (I8) | off-diagonal Jacobian; cross-response vs k_am; dynein single-hand FV | actin-MT net + MTOC→nucleus |

---

## 2. Worktree setup (lead runs once, before launching)

```bash
cd /Users/sw1/ffn_cellsim
for t in fluid-spine nucleus solid-ev motor-nmii emergence; do
  git branch ac/$t ac/new-engine
  git worktree add ../ffn_ac-$t ac/$t
done
# each session opens its own worktree dir: ../ffn_ac-<track>
```
Gated tracks (`weave adhesion mt-coupling production`) get their branch+worktree when their layer opens, off
the then-current `ac/new-engine` (so they inherit landed upstream modules).

---

## 3. Launch-now session prompts

### Session A — fluid-spine (I1a → I1b → I1c)

```
You are the FLUID-SPINE track of the new Active Cell (ac) engine, working in git worktree ../ffn_ac-fluid-spine
on branch ac/fluid-spine. Read (binding): aleph/docs/v2_audit/AC_PARALLEL_SESSIONS_2026-07-16.md §1
(coordination protocol), AGENTS.md, CLAUDE.md, and NEW_ENGINE_BUILD_PLAN_2026-07-16.md §1(P1)/§3(I1a,I1b,I1c)/
§4. Companion: FSI_WIRING_DESIGN_2026-07-16.md, cfd_transport_program/ENGINE_ARCHITECTURE_PLAN.md.

State: ac/fluid/ already holds the I0-B1 ledger (params_i0b1.yaml, alpha=1.0 ratified, native unblocked) and
the committed analytic oracle suite (consolidation_analytic/greens_analytic/manufactured, 25 tests).

Build, in order, NET-NEW under ac/fluid/ (do NOT edit ff/; defer network_warp.py edits to INTEGRATION.md):
- I1a: the conservative moving-domain Biot p/mass Warp solver — aleph/components/fluid/domain.py (OWN the Domain +
  set_nucleus_boundary(provider) contract from §1.4; static-sphere no-flux provider placeholder),
  field_grid.py (conservative FV device arrays + separate Warp HashGrid), biot_substrate.py (S*p_dot +
  alpha*div(v_s) + div(q) = s_water; p=p_ext+p_bar+p_excess), boundary.py (membrane hydraulic-flux BC
  s_water=L_p(dPi-dP); nucleus relative no-flux), scheduler.py (outer physical clock + inner mechanical solve;
  the 6πηR→γ_solid retirement is a network_warp.py patch-note in INTEGRATION.md, NOT an edit here).
- I1b: q=-(k/mu)(grad p - rho b), reconstruct v_f=v_s+q/phi. New oracle: Darcy slab flux linear in dp and k/mu;
  v_f-v_s=q/phi. (phi is I0-B1b — GAP if unsourced, surface to PI.)
- I1c: conservative RAD transport d(phi c)/dt + div(phi c v_f - phi D_c grad c)=R on the same moving domain;
  barbed-end consumption + pointed-end release; exact monomer conservation. New oracle: total-actin
  conservation to machine precision, FRAP recovery↔D_c, advection front follows v_f not q. G_actin
  scalar→monomer-field swap in polymerization_warp.py is an INTEGRATION.md patch-note.

Gates you own (CPU, local): extend the committed oracle suite — constant-state preservation, impermeable-limit
mass conservation, global fluid-content==integrated membrane flux, pressure-work sign, adjoint IBM, net-force
projection, refinement/convergence; plus the I1b/I1c oracles above. Native-gate spec (lead runs on gbook, in
spine order): FSI-OFF regression==Warp-FF bit-identical; resting Π₀=40 Pa pre-tension; drained↔undrained
rate dependence + spatial τ_p≈1 s; device-residency zero-roundtrip; wall-time; CFD-ON HTML.

I0-B params: I1b needs phi, permeability k/mu already provisional (see params_i0b1.yaml); I1c needs D_c
(draft ~2-6 µm²/s, KB-DRAFT-7-03) + c_0 (MCF7 GAP). Surface GAPs, do not choose. Deliver green analytic
modules + native-gate spec + INTEGRATION.md. Answer in Korean.
```

### Session B — nucleus (I2 deformable-mesh)

```
You are the NUCLEUS track of the new Active Cell (ac) engine, in worktree ../ffn_ac-nucleus on branch
ac/nucleus. Read (binding): AC_PARALLEL_SESSIONS_2026-07-16.md §1, AGENTS.md, CLAUDE.md,
NEW_ENGINE_BUILD_PLAN_2026-07-16.md §1(P4)/§3(I2), and companions _historical/MEMBRANE_HELFRICH_DESIGN_2026-07-16.md
(the Helfrich machinery to REUSE) + _historical/CELL_MECHANICS_FRAMEWORK_2026-07-15.md #6 (lamin-A/C vs lamin-B).

Build NET-NEW under ac/nucleus/ (do NOT edit aleph/laws/nucleus_envelope.py — reuse its geometry read-only; stage any
change as INTEGRATION.md): a triangulated nuclear-envelope (lamina) shell using the membrane Helfrich pattern
— dihedral bending + area elasticity + nucleoplasm volume (ν→½) + LINC tether stubs + chromatin polymer net.
Framework-#6 depth: lamin-A/C vs lamin-B split (small-strain→chromatin, large-strain→lamin-A/C
strain-stiffening), nucleoplasm viscosity, envelope-rupture tension threshold.

CRITICAL INTERFACE (§1.4): the fluid track OWNS aleph/components/fluid/domain.py. You provide a
DeformableNucleusMaskProvider conforming to set_nucleus_boundary(provider) (the moving oblate mesh is I1a's
inner relative-no-flux boundary). Code against that signature + a local test-double; do NOT edit ac/fluid/.

Gates you own (CPU, host numpy — mirror the I1a oracle style): sphere→8πκ bending vs FD-gradient (like the
membrane gate); volume conservation; Young-Laplace; small→large-strain modulus crossover at the lamin knee
(report-not-tune); rupture emerges past a sourced envelope-tension threshold (do NOT tune). Native-gate spec
(lead, gbook): resting nucleus stable at physiological E_nuc; oblate flatten at the resting adherent baseline
(volume conserved); reconcile the no-flux mask against I1a's live domain.

I0-B2 params (write params_i0b2.yaml, GAP the unknowns to PI): lamin-A/C vs lamin-B moduli, strain-stiffening
knee, nucleoplasm viscosity, envelope-rupture tension, E_nuc 1-10 kPa (method-dependent), oblate aspect 1.5-3.
Do NOT choose ungrounded values. (I7 cap→LINC→flatten is this track's NEXT increment, gated behind I3+I4 —
pre-author only its standalone oracles now: capstan ∮T·κ ds, oblate volume conservation, LINC/IF
force-extension.) Deliver green module + native spec + INTEGRATION.md. Answer in Korean.
```

### Session C — excluded-volume (I2b steric)

```
You are the EXCLUDED-VOLUME track of the ac engine, in worktree ../ffn_ac-solid-ev on branch ac/solid-ev.
Read (binding): AC_PARALLEL_SESSIONS_2026-07-16.md §1, AGENTS.md, CLAUDE.md, NEW_ENGINE_BUILD_PLAN §2(ADD-row)/
§3(I2b), cfd_transport_program/ENGINE_ARCHITECTURE_PLAN.md §3 (living-cortex steric gap). This implements the
CLAUDE.md HARD worked-example "LJ repulsive excluded volume ON from Phase 1".

Build NET-NEW under ac/solid/: a soft repulsive (WCA/LJ) filament-filament pair potential over filament nodes
via a device wp.HashGrid (port the neighbour-search pattern from dcm/dcm_neighbor_warp.py, read-only). Expose
StericForce.accumulate(state, out_force) per §1.4. Fully field-independent — reads no other increment's live
state.

Gates you own (CPU, host numpy — new aleph/components/solid/wca_analytic.py + tests): pair-potential FD-gradient sign
arbiter; zero-overlap==zero-force (OFF==bit-identical regression); k_EV/CFL bookkeeping (add k_EV to kmax;
assert grid-invariant per the Magic-Number Block — k_EV is a numerical repulsion scale, NOT tuned to a
crowding outcome); compressed-pair limit resists collapse (finite steric volume). Native-gate spec (lead,
gbook): a compressed patch resists volumetric collapse; NO fiber interpenetration on the native cell.
IMPORTANT (verify finding): your pre-I3 native crowding gate certifies a PASSIVE (no active myosin
contraction) cell — scope it to structural/regression only; the load-bearing no-interpenetration-under-
contractile-load verdict DEFERS to a re-run after I3/at I9. Also flag in INTEGRATION.md: the all-fiber EV must
SUBSUME (replace, not add to) the ff/ MT-tip↔cortex soft_contact_kernel, and reconcile with I7's LINC/plate
soft-contact (single-channel invariant — hard-truth #6).

I0-B params: k_EV (numerical, your call via Magic-Number Block, grid-invariance gate — NOT a lit value).
Deliver green module + native spec + INTEGRATION.md. Answer in Korean.
```

### Session D — motor-nmii (I3 head-resolved NMII)

```
You are the MOTOR-NMII track of the ac engine, in worktree ../ffn_ac-motor-nmii on branch ac/motor-nmii. Read
(binding): AC_PARALLEL_SESSIONS_2026-07-16.md §1, AGENTS.md, CLAUDE.md, NEW_ENGINE_BUILD_PLAN §1(P2)/§3(I3)/§4,
companions AFINES_ALGORITHM_NOTES.md and archive/hoomd_legacy/cortex/myosin.py + bridge/motor.py (READ-ONLY
port spec — HOOMD is NEVER imported/executed; it is geometry/kinetics reference only).

Build NET-NEW under ac/motor/: the head-resolved Stam-Hocky NMII minifilament ENTIRELY in Warp — explicit
backbone beads + explicit heads on both sides + per-head actin attachment, Hill force-velocity stepping,
compliance, per-head Bell detachment, all CUDA-resident with no per-step host state. OWN aleph/components/motor/hand.py (the
per-head hand/KMC device API from §1.4, ported from aleph/laws/hand_kmc.py read-only) and MyosinForce.accumulate.
The lumped f_myo / myosin_kernel / gamma_floor._myosin_force DELETE in aleph/laws/network_warp.py is an INTEGRATION.md
patch-note (same-commit-as-fine-motor guard), NOT an edit here.

Gates you own (CPU, host numpy — new ac/motor/*_analytic.py + tests): topology/count invariants; single-head
Hill FV vs Hill-1938 closed form; Bell slip monotonicity; per-head Newton closure; ensemble stall EMERGES
from bound heads (analytic reference — NOT imposed N_side·F_head); ATP/step work sign; force-free at v0.
Native-gate spec (lead, gbook): ensemble stall emerges on the native cell; GPU-residency zero-roundtrip; native
γ-floor (report-not-tune, density-floored FINDING — never add heads/density to close a floor).

I0-B3 params (write params_i0b3.yaml, GAP to PI — do NOT choose): F_stall_head (0.5 vs 2.0 pN by isoform/assay),
N_side (10 AFINES vs 29-30 Billington), v0 (0.12 vs 0.2 µm/s by assay), k_xb (⚠ MASTER force knob, physical
~100-1000 pN/µm — 1 pN/µm breaks stall), r0_head, L_bb/n_bb, k_on. Source by isoform/assay, never to lift γ.
Deliver green module + native spec + INTEGRATION.md. Answer in Korean.
```

### Session E — emergence-detector (I5 detector-core only)

```
You are the EMERGENCE track of the ac engine, in worktree ../ffn_ac-emergence on branch ac/emergence. Read
(binding): AC_PARALLEL_SESSIONS_2026-07-16.md §1, AGENTS.md, CLAUDE.md, NEW_ENGINE_BUILD_PLAN §1(P3)/§3(I5)/§5,
companion _historical/STRESS_FIBER_TARGET_2026-07-16.md.

SCOPE THIS SESSION = the DETECTOR-CORE + synthetic-config oracles ONLY (verify finding: the ablation-toggle
harness needs I3/I4 hooks that do not exist yet — defer it). Build NET-NEW under ac/emergence/: a label-blind
nematic/bundle-condensation detector coded against the FF fiber-array contract (pos / fiber_offsets / n_fibers
— PREDATES I4), reusing the validated Q-tensor S. It must NOT read region/type LABELS (anti-coupling firewall).

Gates you own (CPU, host numpy on SYNTHETIC configs — new ac/emergence/*_analytic.py + tests): isotropic→S~0
null band; aligned→S=1; planted-bundle localization; rotational + permutation/label invariance; monotonicity;
the pre-registered effect-size rule. All run on the dev Mac.

DEFER (write as a spec in ac/emergence/INTEGRATION.md, do NOT build now): the 4 ablation-control toggle hooks
(myosin-off needs I3, unanchored needs FA seed, KMC-off + angle-off need I4) and the native emergence PROOF
(needs I1c + I3 + I4 live). Note the falsifiability contract: a flat result → FINDING to PI + SEEDED-labeled
scaffold fallback, NEVER re-tuned to force condensation.

Deliver green detector + oracle suite + the deferred-harness spec in INTEGRATION.md. Answer in Korean.
```

---

## 4. Gated session prompts (launch when the dep lands)

### Session F — unified-weave (I4) — start after Sessions D (I3 hand API) + A (I1c monomer field) land their module APIs
```
You are the UNIFIED-WEAVE track (ac/weave worktree). Read AC_PARALLEL_SESSIONS §1, NEW_ENGINE_BUILD_PLAN
§1(P3)/§3(I4)/§3A, companions _historical/FILAMENT_SUBSYSTEM_PLAN_2026-07-16.md + _historical/STRESS_FIBER_TARGET_2026-07-16.md. Build
NET-NEW ac/weave/: weave_cell() + WovenCell + topology-reforming crosslink KMC (hash-grid reattach) + the FULL
region set (cortex + ventral/dorsal SF + transverse arc + perinuclear cap + filopodium fascin bundle +
lamellipodium dendritic Arp2/3 REBUILD, §3A.b angle-harmonic branch θ0≈70°±σ_θ). Consume I3's aleph/components/motor/hand.py
API + I1c's monomer field read-only; register KMC presets in aleph/components/weave/presets.py (do NOT edit aleph/laws/hand_kmc.py —
INTEGRATION.md). Pick ONE crosslinker-relaxation channel (double-count guard). Gates: Gate-1
weave_cell([CORTEX]) all-OFF→bit-identical γ; branch-angle dist matches θ0±σ_θ. I0-B4 params: θ0+k_θ, Arp2/3
branch/capping rate, NPF density, fascin bundle stiffness/spacing/n-per-bundle, tip-nucleator density — GAP to
PI. Answer in Korean.
```

### Session G — adhesion (I6) — start after Session F (I4 region API) lands
```
You are the ADHESION track (ac/adhesion worktree). Read AC_PARALLEL_SESSIONS §1, NEW_ENGINE_BUILD_PLAN §3(I6)/
§3A.a/d, companion FILAMENT_SUBSYSTEM_PLAN. Build NET-NEW ac/adhesion/: ONE α2β1-collagen catch-slip clutch
(reuse ff/fa_clutch_warp/fa_ecm/fa_maturation/fa_anchor read-only) on the FULL architecture (ventral SF→FA→ECM
+ lamellipodial nascent + filopodial tip); FA maturation; z-low basal areal-density sampler. The cell_type.py
α2β1 preset is an INTEGRATION.md patch-note (verify: cell_type.py is co-edited with I7 — DO NOT edit it in the
worktree; the lead sequences it). Gates: Pereverzev off-rate + F* oracle; two-sided ECM Newton Σf=0; both
traction channels present; label-blind tip/shaft/basal classes emerge. I0-B6: α2β1-collagen off-rate/F*
(needs its OWN audited evidence — NOT relabeled α5β1-FN), integrin density — GAP to PI. Answer in Korean.
```

### Session H — nucleus load-path (I7) — CONTINUE Session B (nucleus worktree) after D (I3 cap tension) + F (I4 cap region)
```
Continue the NUCLEUS track (ac/nucleus worktree) into I7. Read NEW_ENGINE_BUILD_PLAN §3(I7)/§3A.c, companion
FILAMENT_SUBSYSTEM_PLAN §3A.c. Build ac/nucleus/ cap→LINC→nucleus: contractile perinuclear actin cap draped
over the I2 deformable-mesh nucleus (tangential-tension→normal pressure capstan, NOT a radial strut) + cell-
type keratin(MCF7)/vimentin(MDA) IF secondary passive net. Consume I3 motor + I4 cap region read-only. IF
preset in aleph/components/weave/presets.py; the cell_type.py IF preset is an INTEGRATION.md note (co-edited with I6 — lead
sequences). Gates: analytic capstan ∮T·κ ds (pre-authored in session B); native flatten at rest (oblate,
volume conserved); cap-OFF stays spherical; COMBINED cap+IF baseline (double-load-path guard, same-commit
LINC+soft_contact). Force-scale: K_nuc≈25 nN/µm → 5-6 nN cap ≈0.2 µm ≈10× short — report honestly, do NOT tune.
I0-B7: k_linc (8 pN is a TENSION not a stiffness — GAP), n_cap_fil, IF persistence/strain-stiffening/density,
k_anchor — GAP to PI. Answer in Korean.
```

### Session I — mt-coupling (I8) — oracle-core after D (I3 KMC); module after B (I2 mesh) + F (I4 net)
```
You are the MT-COUPLING track (ac/mt-coupling worktree). Read AC_PARALLEL_SESSIONS §1, NEW_ENGINE_BUILD_PLAN
§1(P5)/§3(I8), companion _historical/CYTOSKELETON_CYTOSOL_REAUDIT_2026-07-16.md §2. SPLIT (verify finding): build the
standalone ORACLE-CORE now (two-node off-diagonal Jacobian; cross-response monotone in k_am; dynein single-hand
Hill/Bell reusing the I3 KMC oracle pattern; MTOC anchor spring) — host numpy, after I3's hand API exists. DEFER
the MODULE wiring (plectin/MACF actin↔MT crosslinker on link_spring reads I4's live actin net; MTOC→nucleus
anchor lands on I2's moving envelope node) until I2+I4 modules land — write it as a spec in INTEGRATION.md.
Dynein = POPULATION of individual hands (engaged count emerges), NOT lumped N×f. Native gate is confounded by
the I7 IF cage → attribution via cross-block Frobenius-norm ON/OFF ablation (lead, gbook). I0-B8: k_am
(re-anchor to plectin/MACF single-molecule — NOT the IF-failure mode), k_mtoc_nuc, N_dyn/k_on_dyn/f0_dyn —
GAP to PI. Answer in Korean.
```

### production (I9) — terminal, lead-owned; not module-parallelizable
Top-level ac/ driver + full-config YAML composing I1a→I8 at physiological setpoints. Only inert scaffolding
(argparse/ledger/viz template) is draftable ahead; the combined `--from-resting` FSI-ON acceptance run is the
heaviest single GPU job at the end of the serial chain (wall-time UNMEASURED — hard-truth #5). All magnitudes
density-floored → report-not-tune (§6.2).

## Change log
- 2026-07-16: created from the verified 16-agent parallelization workflow (wf_d1755894-d6e), adversarially
  revised (shared-file fencing, I5/I8 split, single-A5000 serial-native reality).
