# Remaining Compartment Activations — Design Decisions for PI (2026-06-09)

Under the PI ownership grant, **3 compartments were activated EXPERIMENTAL→LIVE**
(osmotic_regulation, microtubules, intermediate_filaments) + the registry-driven
γ-denylist infra + the `internal_live` integration capstone (all 3 compose in one
cell, no contamination). 308 compartment/cell/cortex tests pass.

Each REMAINING compartment hinges on a genuine **physics/fidelity decision** (not
mechanical wiring). They are NOT rushed unilaterally because the choice changes the
physics interpretation. Each below is a crisp decision + my recommendation; with a
1-line PI answer the wiring is the same proven snapshot-extension pattern.

## 1. LINC — mesoscale acceptor topology + r0 (the make-or-break)
**Problem.** At the ×40 mesoscale the nucleus (R_nuc) ↔ cortex (R_cell) gap is
~µm, but a real nesprin spans ~50 nm. An explicit 50 nm nesprin is SUB-GRID — it
cannot be represented faithfully against the bare cortex (bonds would form at ~µm,
not 50 nm → either zero bonds at r0=50 nm, or massive construction stretch).
**Options.**
- **(A) Bond nucleus → perinuclear IF cage (recommended).** intermediate_filaments
  is now LIVE; its cage sits at [R_nuc, R_nuc+thickness] — perinuclear cytoskeleton.
  LINC bonds nucleus_bead → nearest if_bead at a small, in-grid r0. Physically
  faithful (LINC couples the nucleus to perinuclear IF, a real load path). Needs:
  LINC recipe requires intermediate_filaments; cell.py passes if_bead tags as
  cytoskeleton_tags; per-r0-bin linc_nesprin (module currently single-r0 → small
  addition). f_rest≈8 pN oracle becomes meaningful.
- **(B) Coarse-grained nucleus↔cortex tether.** Bond nucleus → cortex with
  capture~R_cell−R_nuc and r0=construction separation (per-r0 binned). Activatable
  standalone but r0≠50 nm → f_rest oracle N/A; it is an EFFECTIVE coupling, not a
  faithful nesprin. Honest but less faithful.
- **(C) Seed a dedicated perinuclear acceptor cap** (new particles at R_nuc+50 nm).
  Most faithful to the 50 nm span but adds a new particle layer (more wiring).
**Decision needed:** A / B / C. **Recommend A** (faithful + reuses LIVE IF; the
single→per-r0-bin linc_nesprin is a small module change).

## 2. ventral_stress_fibers — FA operating point + sf_myosin_ prefix (BLOCKER #2)
**Problem 1 (prefix).** SF NMII reuses `cortex_myosin_*` bond types = the active-γ
signal. A LIVE SF build with NMII would contaminate cortical γ. Needs a distinct
`sf_myosin_*` prefix (parameterize cortex/myosin.py's builder) so the SF motor load
is denylisted SEPARATELY (the registry-driven denylist already excludes `sf_`; it
must NOT exclude `cortex_myosin_`). **PI sign-off: parameterize cortex/myosin.py.**
**Problem 2 (operating point).** SF anchors on FA clutches → needs an FA-adhered,
equilibrated cell (FA is LIVE but the adhered build needs the equilibration prelude
— the same raw-run BAOAB-guard issue). **Decision:** activate SF (a) passive backbone
only now (NMII off, like the smoke — no prefix needed, no FA-adhered run), or (b)
full active SF (needs the prefix split + the FA-adhered equilibrated build).
**Recommend (a) passive backbone LIVE now** (the bundle structure + tension measure;
NMII deferred to the prefix-split + FA-adhered gate). N_filaments=20 PI-ratify.

## 3. cadherin_junction — the two-cell doublet builder
**Problem.** Cell.build is single-shell. GATE-J needs a `build_cell_doublet()`
(two cortex baselines in one box + interface-cap cadherin seeding + box sizing +
the equilibration prelude). The catch-bond binder + γ-denylist (`cadherin_`) are
ready (the smoke proved the binder forms A↔B trans-dimers). **Decision:** authorize
the new two-cell assembler (the single biggest new build piece, ~1.5-2 sessions).
**Recommend:** yes — it unlocks the entire multicell line; build it as a separate
`build_cell_doublet` (not entangled with the single-cell path).

## 4. membrane_reservoir — own mem_node layer + bleb σ_crit
**Problem 1 (BLOCKER-1).** membrane_surface rides the cortex shell tags, so a tether
mesh would be degenerate self-pairs. Needs an own radially-offset mem_node bead
layer (the smoke fabricated one standalone). **Problem 2.** σ_crit_bleb + f_excess
are None (PI-pending); the rupture/reservoir paths stay blocked. **Decision:**
activate the STATIC tether mesh (own mem_node layer + builtin mem_tether on the
shared bond — needs the offset-layer seeding in cell.py), bleb/reservoir deferred
until σ_crit_bleb + f_excess ratified. **Recommend:** static mesh LIVE; bleb deferred.

## 5. junctional_actin — STUB build path (depends on #3)
The enabled+anchored HOOMD build raises NotImplementedError (reserved); the scalar
laws work (smoke: F*≈6 pN, engaged f≈0.91). Depends on cadherin_junction (#3) first.
**Decision:** ratify the α-catenin catch set (k_catch0=1.0/k_slip0=0.02 ORDER) +
implement the reserved build path AFTER #3. Lowest priority (third leg).

## 6. (engineering, no physics decision) MTOC multi-bead core
microtubules is LIVE but n_mt is CAPPED at 7 (single-hub MTOC degree vs the nlist
exclusion cap). The PI-candidate dense aster (n_mt=20-250) needs a multi-bead MTOC
core in build_mt_topology (distribute arms across ceil(n_mt/5) core beads). Pure
engineering (no physics decision) but churns the 35 MT tests + 3 smokes. **Do when
a dense aster is needed; n_mt=7 sparse aster is LIVE meanwhile.**

---
*Companions: PLATFORM_PI_QUEUE.md, COMPARTMENT_SMOKE_REPORT_2026-06-09.md,
COMPARTMENT_ACTIVATION_DETAILED_PLANS_2026-06-09.md, PLATFORM_AUTORUN_LOG_2026-06-09.md.*
