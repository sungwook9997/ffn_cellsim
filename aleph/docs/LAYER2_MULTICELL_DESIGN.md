# Layer-2 multicellular spheroid model — design brief + scale-bridge contract

> **STATUS: PREP brief — 2026-06-02. NOT yet a ratified contract.** Specifies a NEW,
> isolated, parallel simulation line (multicellular tumor spheroid) that complements the
> existing fine-grained single-cell line. This doc is an **additive, untracked NEW file** —
> the Lead's working tree, the live single-cell runtime, the `integrator/` freeze, and the
> existing test suite are **bit-for-bit unaffected** (same isolation discipline as
> `PARALLEL_PREP_INDEX_2026-06-02.md`). Two background surveys back this brief: a literature
> model-class survey (7 spheroid papers) and a codebase reuse inventory of `aleph/`.
> Companion: `PI_EXP_VALIDATION_MAP.md` (Layer-1/Layer-2 split), `MIKADO_3D_DESIGN.md`
> (the shared 3D-ECM asset). PI sign-off needed before any code lands (touches the
> anti-coarse-graining hard rule + echoes retired v1 territory → see §9).

---

## 0. One-paragraph summary

Build a lightweight **center-based multicellular spheroid simulator (CBM: 1 particle per
cell)** on the **same HOOMD-blue + BAOAB stack** as the single-cell line, reusing the
substrate, ligand-kinetics, integrator, and I/O machinery already in `aleph/`. It reaches
the **collective / multi-cell** terms of the PI experiment's spreading law
`A/A₀ = a + b/R + c/R²` (cohesion `c`, `A/A₀` itself) that the single-cell line can never
reach in wall-time. The two lines are **isolated as code, coupled as physics**: the
single-cell line is the *ground-truth calibration source* for the spheroid's per-cell
parameters (cortical tension, traction, adhesion energy) via an explicit **scale-bridge
contract** (§6). The discipline that separates this from a re-tread of retired v1
(Taichi MLS-MPM spheroid continuum) is: parameters are **calibrated bottom-up / literature-
anchored, never free-fit**, and `A/A₀` is an **overlay-only validation target, never a
fitting target**.

---

## 1. Why now, why parallel

- **Reaches the experiment's actual body.** PI's KSME study is a *collective* MCF7-spheroid
  phenomenon. `A/A₀` and its `c/R²` cohesion term are multi-cell; the single-cell line
  reaches only the `b` (traction) term (`PI_EXP_VALIDATION_MAP.md` §"scale bridge"). 38k
  filaments × hundreds of cells is wall-time-impossible → a coarse multicell line is not
  optional if we want to touch `A/A₀`.
- **De-risks the slow bottom-up build.** Gives an early, end-to-end, experiment-touching
  vertical slice while the fine line matures.
- **Cheap on shared infra.** CBM = 1 HOOMD particle/cell; a 10²–10³-cell spheroid ≈ one
  ×40-meso single cell in particle count. Reuses BAOAB, substrate, ligand kinetics, GSD.
- **Fits the gbook concurrent-CPU model.** Light enough for many parallel CPU jobs;
  R-sweeps (§7 phase L2.3) are embarrassingly parallel.
- **The literature just landed.** Today's uploads are exactly this corpus (spheroid
  mechanics, directed invasion, multiscale tumor) → §5 model-class survey.

---

## 2. Architecture — isolated as code, coupled as physics

```
  PI experiment (validation north star · OVERLAY-only, never a fitting target)
  MCF7 spheroid · A/A₀ = a + b/R + c/R²    [a baseline · b traction/curvature · c cohesion]
  ─────────────────────────────────────────────────────────────────────────────
        a, b/R ──┐                                  ┌── c/R²,  A/A₀ itself
                 ▼                                  ▼
  ┌───────────────────────────┐      ┌────────────────────────────────────┐
  │ ① SINGLE-CELL LINE         │      │ ② MULTICELL LINE  (THIS BRIEF)      │
  │   (main · existing · slow) │      │   (new · isolated · light)          │
  │                           │      │                                    │
  │  fine-grained mechanism    │      │  CBM: 1 particle / cell            │
  │  filament/motor/clutch     │      │  cell–cell adhesion = pair force   │
  │  cortex/FA (×40 meso)      │      │  cell–ECM traction = active force  │
  │  → traction b, β1 pattern  │      │  → cohesion c, A/A₀                │
  └─────────────┬─────────────┘      └──────────────▲─────────────────────┘
                │      ⑥ SCALE-BRIDGE (the only      │
                └───────  intended coupling) ────────┘
                   per-cell params (bottom-up):
                   cortical tension · adhesion energy · traction · drag
        (both on the same HOOMD + BAOAB stack, same SI units, same box scale)
```

- **Code isolation:** Layer-2 is a NEW package `aleph/spheroid/` + NEW oracle subpkg +
  NEW config tree on a **separate branch** (`layer2/spheroid-cbm`). It *imports* reuse-ready
  modules read-only; it does **not** edit any Lead-owned file, the `integrator/` freeze, or
  any existing gate. New files ⇒ existing line bit-for-bit unaffected.
- **Physics coupling:** exactly one channel — the scale-bridge (§6). Single-cell observables
  parameterize Layer-2 cells. No other coupling.

---

## 3. Shared / reuse-ready (no new work) — from the codebase inventory

| Module | Provides | Reuse | Note |
|---|---|---|---|
| `integrator/baoab.py` | L-M BAOAB overdamped Langevin `Action`; per-type γ | **REUSE** | works for 1-particle/cell (no bonds needed). `make_baoab_updater(kT,gamma,dt,seed)`. **CPU-only** (per-step `cpu_local_snapshot`) — fine for 10²–10³ cells; GPU port deferred (template = `constrained_baoab.py`). |
| `bridge/ligand_species.py` | col-I→α2β1, laminin-111→α6β1, FN→α5β1 kinetics; `DEFAULT_LIGAND_FOR_CONDITION` (Bare/Pre/Lam4) | **REUSE verbatim** | pure-python, device-agnostic, import-isolated. Source-of-truth for cell–ECM adhesion identity. |
| `validation/pereverzev.py` | Pereverzev/Bell-Evans catch-slip oracle | **REUSE** | calibration + the cadherin-catch template. |
| `ecm/substrate.py` | `SubstrateAnchorSpring` `F=−k_sub(r−a)`, tunable rigid↔soft | **REUSE** | body-agnostic → spheroid basal cells anchor the same way. **This is the immediately-shared ECM (§4).** |
| `ecm/equilibrate.py` | clipped-Brownian + BAOAB settle (anti-blowup warmup) | **REUSE** | warmup pattern for the initial aggregate. |
| `common/gsd_traj.py` | GSD trajectory writer (live + offline, scalar fields, bonds) | **REUSE** | viz unlock; hook `attach_gsd_writer`. |
| `common/checkpoint.py`, `common/integrity.py` | fingerprinted restart; oracle-import guard | **REUSE** | adopt the discipline. |
| **Architectural patterns** | `resolve_derived(cfg)→Resolved* dataclass`; dual-config (`oracles/configs/*` + runtime `configs/*`); `acceptance:` band block; additive default-off module | **ADOPT** | Layer-2 follows the same conventions (Magic-Number Block, gates-before-run). |

## 3b. Must build NEW

1. **Center-based cell body** — 1 particle/cell with soft adhesive+repulsive pair (no such
   body type exists; everything today is sub-cell filaments).
2. **Cell–cell adhesion** — `junction/` is **empty** (0-byte `__init__.py`). No cadherin
   runtime exists anywhere. MVP = static pair well; mechanistic target = KU-4.2 catch-bond
   (§5, §7-L2.5). The Young contact-angle oracle (`oracles/junction/contact_angle.py`)
   already exists as the acceptance gate.
3. **Spheroid-scale observables** — spread area (convex-hull / α-shape of cell centers),
   radial density profile, detached-cell fraction (point-cloud metrics, ref [3]); the
   `A/A₀(R)` fit harness.
4. **Layer-2 oracle subpkg + config tree** — `validation/oracles/spheroid/` +
   `configs/layer2_*.yaml` with `acceptance:` bands; runtime-import-forbidden.
5. **Scale-bridge** — single-cell → per-cell-param coupling (§6); none exists.
6. **(Shared, later) 3D Mikado ECM + cell↔fiber binding** — for spheroid-in-collagen
   invasion; already designed in `MIKADO_3D_DESIGN.md`. See §4.

---

## 4. The shared-ECM finding (the user's headline question)

**Surprise from the inventory: the single-cell line never actually embeds in the Mikado
fiber network.** The H.1 Mikado builder and the Cell builder are never combined; the cell's
"ECM" is a synthetic **z=0 ligand plane** held immobile by a position-reset pin
(`cell/cell.py::SubstrateLigandPin`), optionally swapped for the finite-stiffness
`SubstrateAnchorSpring`. `fa.py` even documents "ligand plane replaced by H.1 ECM particles
in week 3" — **never done.**

Consequences for "can both lines share one ECM?":

- **Immediately shared = the substrate layer.** `ecm/substrate.py` + the z=0 ligand/pin
  pattern is body-agnostic. A spheroid's basal cells anchor to the same tunable-compliance
  substrate the single cell uses. Units are already identical (kT=4.28e-21 J, T=310 K,
  η=6.913e-4 Pa·s). **Use this for the 2D-spreading conditions (Bare/Pre/Lam4 on the dish).**
- **Box scale already favors Layer-2.** ECM box = 200 µm cube; cell boxes = 20–30 µm. A
  100–300 µm spheroid fits the 200 µm box *better* than a single cell does → Layer-2 adopts
  the ECM `L_box` scale.
- **A true shared 3D fiber ECM = build-new, but already designed and wanted by both lines.**
  Current Mikado is z=0-only (2D-extruded); no body↔fiber binding exists. `MIKADO_3D_DESIGN.md`
  already specifies the 2D→3D generator (and closes the single-cell line's logged 2D-⟨z⟩
  gap). So the 3D fiber ECM is a **shared deliverable**: the single-cell line gets its 3D
  gap closed, Layer-2 gets spheroid-in-3D-collagen invasion (ref [4], `Directed invasion of
  cancer cell spheroids inside 3D collagen matrices`). Deferred to phase L2.6 — not on the
  MVP critical path.

**Net:** share the **substrate** now (free); share the **3D Mikado fiber ECM** later (one
build, two beneficiaries).

---

## 5. Model class — decision + literature basis

**Decision: Center-Based Model (CBM, 1 particle/cell) first cut → graduate to a few-element
Subcellular-Element Model (SEM) only if cell shape / deformation-under-confinement /
division asymmetry become first-order to a result.**

Why (from the 7-paper survey):
- The only papers that run whole spheroids of 10²–10³ cells fast on CPU **as particles with
  explicit pairwise forces** are **Chen & Zou 2018 [7]** (off-lattice ellipsoids, overdamped
  force balance) and the **Palsson–Othmer** family it builds on — i.e. exactly the HOOMD
  primitive set (pair potential + Langevin). Lowest-friction port.
- The CPM papers (**Herold 2023 [3]**, **Rosenbauer 2023 [6]**) are scientifically adjacent
  (same observables, cell sorting, invasion) but are **lattice Monte-Carlo** — they violate
  the HOOMD-MD constraint and their energies aren't SI-calibratable. **Use them for
  observables + validation, not as the runtime.**
- The deformable-mesh models (**Odenthal/Smeets 2013 [2]**, **Fang & Lai 2016 [5]**) are the
  *right* level for the single-cell line (where cortical tension / traction / adhesion energy
  are physically defined) but far too heavy for 10³-cell spheroids — confirming the two-layer
  split.
- CBM is the natural **bottom-up calibration target**: every parameter is a physical
  mechanical quantity the single-cell line (or [2,5]) hands up.

### Interaction forms → HOOMD primitives

| Physics | MVP form (lit-anchored, runnable now) | Mechanistic target (hard-rule-compliant) |
|---|---|---|
| **Cell–cell adhesion + excluded volume** | **Morse** `U=D_e[e^(−2β(r−r₀))−2e^(−β(r−r₀))]` → `md.pair.Morse`; `D_e`=adhesion energy/contact, `r₀`=cell rest sep, `β`=inverse range | **KU-4.2 cadherin catch-bond** (force-dependent bind/unbind) replacing the static well — Bell-Evans/Pereverzev, the project's ratified cadherin mechanism. Optional higher-fidelity contact: **Palsson** (Chen & Zou Eq.12) or **JKR/Hertz** via `md.pair.Table`. |
| **Cell–ECM traction** | active force toward substrate / along polarity → `md.force.Active`/custom; magnitude from [5] Table I, [7] T_ij | coarse integrin clutch reusing `ligand_species.py` kinetics (Pereverzev/Bell-Evans) on the z=0 substrate |
| **Active motility** | self-propulsion, constant magnitude, persistence τ, reorient every τ → `md.force.Active` | fiber-guided bias in 3D ECM (k∥ along / k⊥ across, [4]) |
| **Integrator / drag** | overdamped BAOAB (`integrator/baoab.py`); per-cell γ from [7] μ_f/μ_cell/μ_s | (same) |
| **(SEM upgrade)** intra-cell bonds | FENE/harmonic beads + soft volume restraint; cortical constants from [5] WLC+POW (μ₀↔k_POW,L_p) | reuses the *same* cortical-tension numbers the single-cell line produces |

> **Hard-rule note.** A cell = 1 particle makes cell–cell adhesion *intrinsically* a coarse
> interaction — this is the sanctioned ×N scale (same logic as the ratified ×40 meso). The
> Morse well is an explicit **MVP scaffold** (interpretable `D_e`, lit-anchored), and the
> contract commits to the **KU-4.2 catch-bond** mechanistic upgrade from day 1 (additive
> default-off pattern, like `cell.py`). That upgrade path is what keeps this on the
> mechanistic side of the hard rule rather than a lumped continuum.

---

## 6. SCALE-BRIDGE CONTRACT (the crux)

For each Layer-2 cell parameter: the single-cell (Layer-1) observable that calibrates it,
and the literature anchor used as **fallback** (until Layer-1 produces the number) **and**
as an independent cross-check. **No parameter is free-fit; `A/A₀` is never an input.**

| Layer-2 param | Symbol | Calibrated from (Layer-1) | Literature anchor (fallback + cross-check) |
|---|---|---|---|
| cell–cell adhesion energy | `D_e` | single-cell cadherin catch-bond × contact area (junction, KU-4.2) — *not built yet → starts lit* | W=h₀σ₀ Odenthal [2]; σ(adh,cortical) Boot [1]; E_b Chen&Zou [7] |
| cell stiffness / repulsion | `β`,`Ê` | whole-cell modulus from cortex (H.3) | cortical shear μ₀ Fang&Lai [5]; AFM E₀ Yubero (PI_EXP map) |
| cell rest separation | `r₀` | single-cell diameter (cell geometry) | cell radius (lit) |
| **cell–ECM traction** | `F_tr` | **single-cell traction = KU-3.5 grip-walk = the live γ-production** | traction→area Fang&Lai [5]; T_ij Chen&Zou [7] |
| per-cell drag | `γ` | single-cell Stokes drag (η, cell size) | μ_f/μ_cell/μ_s Chen&Zou [7] |
| motility persistence | `τ` | single-cell polarity timescale (lamellipodium) | Geiger [4]; Rosenbauer [6] |
| ECM adhesion kinetics | `k_off` | **(shared)** `ligand_species.py` col-I/laminin/FN | already lit-anchored |

> **The beautiful coupling:** `F_tr` (Layer-2's traction, the `b` term) comes directly from
> the single-cell **KU-3.5 grip-walk** — i.e. *the main line's current γ-production output
> is literally Layer-2's traction input.* The `c` (cohesion) term waits on the single-cell
> cadherin junction (KU-4.2), which `junction/` doesn't have yet → starts literature-anchored
> (matches the EXTEND deferral), upgrades when Layer-1 produces it.

---

## 7. Phased build plan (each phase = runnable + gated)

- **L2.0 — Isolated line setup.** Branch `layer2/spheroid-cbm`; package `aleph/spheroid/`;
  `resolve_layer2(cfg)→ResolvedL2` (all-SI, Magic-Number Block); dual config
  `configs/layer2_cbm.yaml` + `oracles/configs/layer2_cbm.yaml` (`acceptance:` bands written
  before any run); `oracles/spheroid/` subpkg. Reuse baoab + Morse + gsd_traj.
- **L2.1 — CBM MVP (hello-spheroid).** N cells = N particles, Morse pair (lit `D_e,β,r₀`),
  BAOAB overdamped, no motility. **Gate G1 (written now):** a loose aggregate neither
  explodes nor disperses — equilibrium nearest-neighbor spacing within [0.9, 1.2]·r₀, RDF
  first peak at r₀, packing fraction stable over 10⁵ steps. (= excluded-volume + adhesion
  balance correct.)
- **L2.2 — Motility + substrate.** Add active force + persistence + `SubstrateAnchorSpring`
  (z=0). **Gate G2:** spheroid on substrate spreads monotonically; spread rate ↑ with
  `F_tr`, ↓ with `D_e` (sign-sense), matching [5] traction→area sign and [2] a∼√(2W/γ)·√t.
- **L2.3 — A/A₀(R) sweep.** R₀ sweep (parallel CPU jobs) → spread-area observable → fit
  `A/A₀ = a + b/R + c/R²`. **Gate G3:** the 1/R, 1/R² terms emerge from the edge/bulk cell
  ratio (surface/volume ∼1/R) with no hard-coded term; **PI poster overlaid, NOT fit.**
- **L2.4 — Ligand-identity axis.** Bare/Pre/Lam4 via `ligand_species.py`
  `DEFAULT_LIGAND_FOR_CONDITION` → 3-condition contrast; β1-pattern analog (diffuse/
  peripheral/uniform). **Gate G4:** condition ordering matches the experiment's *qualitative*
  spreading-threshold ranking (overlay, not fit).
- **L2.5 — Cadherin catch-bond (mechanistic upgrade).** Replace the static Morse well with
  the KU-4.2 force-dependent bond updater (Pereverzev template). **Gate G5:** Young
  contact-angle oracle; catch-bond lifetime-vs-force peak.
- **L2.6 — 3D Mikado ECM (shared asset).** Build the 3D fiber network (`MIKADO_3D_DESIGN.md`)
  + cell↔fiber binding → spheroid-in-collagen directed invasion. **Gate G6:** invasion
  distance ∼√t, fiber-orientation-dependent (Geiger [4]).
- **Later — SEM upgrade** only if shape/division/confinement become first-order.

---

## 8. Sanity gates & discipline

- Every numeric constant satisfies the **Magic-Number Block** (literature-anchored per §6 or
  a ratified numerical-policy choice); no value chosen to pass a gate.
- Gates G1–G6 are **written before the run** (contract), in `configs/layer2_*.yaml`
  `acceptance:` blocks. No inline gate-loosening — a wrong gate goes to PI.
- **No fitting to PI data.** `A/A₀` poster values are overlay-only at comparison time.
- Auto-viz at every phase closeout (`spheroid/layer2_vis.py`, baked into the driver per the
  production-driver-auto-viz rule).
- CFL / numerical-sanity check at the BAOAB timestep for the new force magnitudes.

## 9. Governance

This touches (a) the anti-coarse-graining hard rule's *interpretation* (resolved by §5–6:
mechanistic upgrade path + bottom-up calibration), (b) retired v1 territory (resolved: CBM
discrete + calibrated ≠ v1 lumped continuum + fit), (c) a parallel workstream. → **PI
sign-off before L2.0 code lands.** Branch stays off `ffn/foundation`. The `integrator/`
freeze is untouched (Layer-2 reuses baoab read-only).

## 10. References (today's uploads)

- [1] Boot, Koenderink & Boukany 2021, *Adv. Phys. X* 6:1, 1978316 — spheroid mechanics review (σ = f(adhesion, cortical tension); DAH/DITH/HIT). *Calibration backbone + validation.*
- [2] Odenthal, Smeets, Van Liedekerke et al. 2013, *PLoS Comput Biol* 9(10):e1003267 — deformable-cell spreading; **W = h₀σ₀**, JKR/Maugis-Dugdale; a∼√(2W/γ)√t. *Adhesion-energy anchor.*
- [3] Herold et al. 2023, *PLoS Comput Biol* 19(3):e1010471 — spheroid scoring; point-cloud observables + Wasserstein. *Observables (reject CPM runtime).*
- [4] Geiger et al. 2022, *PLoS ONE* 17(3):e0264571 — directed spheroid invasion in 3D collagen; fiber-guided random walk (k∥/k⊥); invasion ∼√t. *L2.6 validation + motility recipe.*
- [5] Fang & Lai 2016, *Phys. Rev. E* 93, 042404 — traction-driven spreading; cortical μ₀↔WLC+POW; traction→area table. *Traction + cortical anchors.*
- [6] Rosenbauer et al. 2023, *J. Phys. Chem. B* 127, 3607 — multiscale spheroid (CPM); adhesion ∝ contact area; sorting. *Concept + validation (reject CPM runtime).*
- [7] Chen & Zou 2018, *Math. Biosci. Eng.* 15(2):361–392 — off-lattice center-based spheroid (Palsson–Othmer Eq.12); SI-tabulated E_b, α, μ, growth. **Primary portable model + ready CBM param set.**

The load-bearing form to port: **Chen & Zou Eq. 12** (Palsson adhesion+compression). The
load-bearing adhesion identity: **W = h₀σ₀** with JKR/Maugis-Dugdale forms (Odenthal [2]).
