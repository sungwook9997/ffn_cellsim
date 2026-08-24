# Wiring the detectors: what tether and manifest now read, and what the other three would need

**Status:** two wired (`tether`, `manifest`), three analysed and **left for PI judgement**
(`fluctuation`, `kinematics`, `traction`), one status-checked only (`stationarity`).
Follows `PARALLEL_AUDIT_PARAMS_DETECTORS_2026-08-20.md` §1. Session C, 2026-08-20.

Every claim here is a code fact checked on this tree today, or a test that runs. Nothing was
measured on a device and no number below is a physics result.

---

## 0. Correction to the audit's framing — read this first

The audit says *"`tether` is the forward operator that would read that instrument"*, meaning the
~3-minute AFM protocol the 6,514× throughput target is sized from. **That is the wrong pairing, and
the real situation is worse than the audit states.**

| | what the audit implies | what is on disk |
|---|---|---|
| the AFM harness that IS wired | a tether pull | a **spherical indenter** — `outputs/ac/afm_cortical_tension_readiness/`, `engine/afm_indenter.py`, `scripts/ac_afm_*` |
| what `tether.py` reads | that instrument | a **bead-pull tether plateau** — a different experiment |
| indentation forward operator | — | **does not exist**, anywhere in `aleph/observe/` (grep, 2026-08-20: `indent`/`Hertz` appear only in docstring examples) |

So the throughput-budget protocol does not have an operator with zero call sites. It has **no
operator at all**. `tether` was still the right first thing to wire — see §1 for why — but closing
the audit's own stated gap needs an indentation operator authored, which is a new-module decision
and is surfaced here rather than taken.

---

## 1. `tether` — WIRED

### What it reads
The **membrane** population, and only that: its per-frame total apparent area [µm²], the quantity
`membrane_area_reduce_kernel` (`laws/membrane_surface.py:102`) already sums on device over `faces_d`.
In arena terms, the FACE block of the membrane strand-set. Everything else is a resolved constant.

    membrane.tension_pn_per_um            <- area_tension(A(t), A0, MembraneAreaCard)
    membrane.bending_rigidity_pn_um       <- ResolvedMembrane.kappa       (0.0828 pN·µm, KB-3.B1.2)
    cortex.membrane_attachment_energy…    <- GAMMA_MCA_PN_UM              (10.0 pN/µm, KB-3.B1.4)

### What it puts out
Per-frame plateau force [pN], with the tether radius [µm] and the apparent tension [pN/µm] on the
support, reduced to one mean with a **correlated** standard error (`σ√(2τ/N)`, never `σ/√N`).

### Is the output comparable to experiment?
Yes in kind, with two named limits the operator states itself: `f² ` is *linear* in σ_app, so a
tether force is a **square-root probe** — a 2× tension change moves the force 41%; and one force
cannot split σ_app into σ_bilayer + W. The radius is the observable that breaks the (κ, σ) scaling
degeneracy, which is why it is on the record.

### ⚠ Is the setting right? Partly — and the wiring is what makes that visible

**The forward map was never missing from the runtime. It was there stripped of its protocol.**
`laws/membrane_surface.py:230 erm_rupture_force` is the *same closed form* the operator implements:

    erm_rupture_force(κ, γ_mem, γ_MCA)  ==  equilibrium_tether_force_pn(κ, γ_mem + γ_MCA)

`components/incumbent/compartments.py:644` evaluates it on every native run and hands it to
`erm_tether_kernel` as a **rupture threshold** — no manifest, no error bar, no identifiability
statement, no record. `test_the_runtime_already_computes_this_force_as_a_rupture_threshold` pins the
two together, so a silent divergence fails a test instead of surfacing in a run six weeks later.

**On today's runtime the number is an echo of three constants.** γ_mem is the constant
Raucher–Sheetz plateau (`compartments.py:579`, the `K_A` upturn default-OFF), κ is resolved
constant, `GAMMA_MCA_PN_UM` is a literal. Zero scatter. The operator computed
`standard_error = 0.0` for that case — a claim of *infinite precision* on a number that measured
nothing. Two fixes, both root-cause rather than per-caller:

* `tether.apparent_quantity` reports `standard_error = None` with the reason on the record when the
  plateau has no scatter. `ApparentQuantity.standard_error` was already documented as *"never zero
  as a stand-in for 'not computed'"*; nothing enforced it.
* `ApparentQuantity.__post_init__` now **refuses** `standard_error == 0.0`. One guard, all seven
  operators, because every one of them reaches it by the same route — a constant input series.

The driver stamps `INPUT_ECHO_NOT_A_MEASUREMENT` for that case. **That verdict is the point of
wiring it**: the detector was always going to say this, and until something called it, nothing did.

The one input that can make the force move is the membrane's own **area**, through
`engine/membrane_area.area_tension` — which had **zero call sites of its own**. The driver calls it,
so a run whose area varies yields a varying force and a real error bar, with no edit to the driver.
Both branches are pinned by tests.

### PI-GAPs this opened
1. **`GAMMA_MCA_PN_UM = 10.0` is a literal, not something the ERM population produces.**
   `erm_density_per_um2` is `ABSENT_FROM_CONTRACT_GRAPH` (audit §3). While W is a constant, `tether`
   is permanently an input echo. Source it, derive it from the ERM population, or declare it an
   inference output — all three are PI calls.
2. **`MembraneAreaCard.capacity` has no sourced value** for any cell type here; the retired
   `RESERVOIR_STRAIN = 0.60` says in its own comment that it is unsourced for MCF7. Kept required,
   so the driver refuses rather than inventing it.
3. **No indentation forward operator** — §0.

### Files
`aleph/scripts/ac_observe_tether_native.py`, `aleph/tests/ac/engine/test_observe_tether_wiring.py`
(7 checks, all passing), and two edits in `aleph/observe/{tether,operator}.py`.

---

## 2. `manifest` — WIRED, in the same diff

`ProtocolManifest.is_complete()` and `not_recorded_fields()` are called for the first time in this
tree. Every one of the 34 fields is required of the declaration, with no defaults; a missing field
is a refusal that names it.

**It is deliberately NOT a hard gate.** `is_complete()` is false while *any* field is
`NOT_RECORDED`, and refusing to write on that would push a caller to assert `NOT_APPLICABLE`
instead — which is the exact substitution the two-member enum exists to prevent. A confession must
not be more expensive than a false assertion. So the confession rides on the record
(`manifest_complete`, `manifest_not_recorded_fields`) and **downgrades the verdict** to
`MEASURED_UNDER_INCOMPLETE_PROTOCOL`.

Why the downgrade has teeth: `evidence/dataset.py`'s `ExternalDatasetManifest` already carries
`protocol_manifest_digest` and joins an external record to a simulated one **on the manifest hash**.
Until now only the external side had a manifest, so the join had one end. Two protocols that differ
in any field — including `NOT_APPLICABLE` vs `NOT_RECORDED` — get different digests, pinned by
`test_declaring_the_last_eleven_fields_lifts_the_downgrade_and_moves_the_hash`.

---

## 3. The three left for PI judgement

Analysed, not wired, as instructed. Each answers the four questions: what it reads, what it puts
out, whether that is comparable to experiment, and whether the setting is right.

### 3.1 `fluctuation` (904 lines) — the cheapest to wire, and the setting is *wrong today*

**Reads** `membrane.normal_displacement_um`, shape `(T, V)` — the per-vertex normal displacement of
the membrane over T accepted frames — plus a `ModeBasis` supplying `Φ (V,M)`, per-mode `q` [1/µm],
per-vertex areas and mean vertex spacing.

**Puts out** a mode power spectrum `P(q)` [µm²] over the mesh-measured trusted band
(`q_max = π/a` Nyquist, and a patch fundamental below), with per-mode `n_samples` **and**
`n_effective`, and the standard error taken from `n_effective`.

**Comparable?** Yes, and it is the *most* comparable of the seven: no probe geometry, no contact
model, no regularisation — which is exactly why `ALEPH-DQ-102` picked it first. Flicker-spectroscopy
`P(q)` is directly published.

**⚠ Is the setting right? No — and this is the substantive finding of §3.**
The operator measures **thermal** fluctuation. There is no thermal noise term on the membrane in
this runtime, verified three ways on this tree: the membrane's `accumulate()`
(`compartments.py:579-590`) launches only Helfrich bending, area tension and the ERM tether — all
deterministic; `membrane_surface.py` and `compartments.py` contain **no** `rand`/`randn`/gaussian
call; and `BAOAB` appears **only under `aleph/dcm/`**, which `STATE.md` (a) records as parked and
reached solely through `laws/relax.py::implicit_overdamped_step`. `kBT` does appear in the engine,
but as an *energy scale* for WLC and Bell kinetics, never as a fluctuating force.

So wiring `fluctuation` today reproduces `tether`'s failure exactly one level up: a spectrum from a
deterministic relaxation is not a thermal spectrum, and it would have the right units and a
plausible shape. **The prerequisite is a PI decision on whether the membrane gets a fluctuating
force term**, not an operator call. Recommend: do **not** wire until that is answered.

*Cheapest real prerequisite:* the driver would also need `normal_displacement` recorded per accepted
frame — the membrane block records positions, not displacements from a reference surface, so a
reference and a projection have to be declared. Both are cheap; the noise term is not.

### 3.2 `kinematics` (728 lines) — partly wireable now, and one of its three is missing

**Reads** three different things through three operators:

| operator | array | shape | does the runtime produce it? |
|---|---|---|---|
| `spreading_area_rate_operator` | `cell.projected_area_um2` | `(T,)` | ❌ **nothing computes a projected/silhouette area** in the current tree |
| `volume_flux_operator` | `cell.volume_um3` | `(T,)` | ✅ computed every step for turgor — `engine/fluid_core.py:946` |
| `piv_field_speed_operator` | `field.velocity_um_per_s` | `(T,N,3)` | ~ `surface_velocity_d` exists (`fluid_core.py:290`); a *field* PIV analogue would have to be declared |

**Puts out** `d(A)/dt` [µm²/s], **signed** `d(V)/dt` [µm³/s] and mean PIV speed [µm/s].

**Comparable?** Only in a restricted sense, and the module is unusually careful about it. The
external contrasts (eLife 72381: `initial_spreading_area_rate = 34.18`, `initial_volume_flux =
−16.23`) **state no units**. So magnitudes are refused by name and `direction_agreement()` is
offered instead — a sign survives an unknown positive scale, a magnitude does not. And
`DirectionAgreement` carries `is_acceptance_criterion = False` as a *field*, so a direction check
cannot quietly harden into a gate.

**Is the setting right?** For volume flux, yes — the quantity exists, is signed, and the sign is the
comparable part. This is the single cheapest genuine detector wiring available: it needs one scalar
recorded per accepted frame from a value the turgor law already computes. For projected area,
**no** — the observable does not exist and would have to be defined (top-down silhouette is the
project's established convention; it is not implemented). PIV needs a declared field.

**Recommendation to PI:** wire `volume_flux` only, and only against a *direction* contrast. It is
the one detector here whose input exists, varies, and whose comparison discipline is already
settled.

### 3.3 `traction` (619 lines) — ⚠ refusing on a citation whose code no longer exists

**Reads** `adhesion.traction_force_pn` and `adhesion.site_area_um2`, **plus** a separate audit owner
`load_path_audit` carrying `booked_force_pn`, `owner_displacement_um`, `deliver_call_count` — one
entry per connector row.

**Puts out** traction force per cell [pN] and stress [Pa].

**Comparable?** Magnitudes: **no**. The external TFM tables arrive tagged `source_unit_unknown` and
the builder "never infers missing units", so the external traction values are directions and ratios,
not pN. The operator says so itself.

**⚠ Is the setting right? The refusal machinery is sound; its stated reason is stale.**
The two open defects it names —

* `aleph/scenarios/whole_cell_couplings.py:735` (the discarded `ProtrusionActinTarget`)
* `aleph/scenarios/whole_cell.py:1643` (`_probe_delivery` scoring from the ledger, not from motion)

— **cannot be checked: `aleph/scenarios/` does not exist anywhere in this tree.** `ProtrusionActinTarget`,
`deliver_forces_pn` and `_probe_delivery` now appear **only inside `traction.py`'s own docstring and
`observe/__init__.py`**. So the module's promise — *"when the owning lanes fix the engine, this
operator begins reporting without a line of it being edited"* — points at a lane that was deleted.

This does **not** break the guard, and that is the design working: the precondition is a runtime
check against state-visible arrays (*a row that books force and moves nobody*), not a file path, and
an **absent** audit is a `MISSING_STATE_OWNER` refusal rather than a default to "sound". Fail-closed
is correct here.

But the practical consequence is that traction is unblockable by construction: **nothing in the tree
produces a `load_path_audit` block.** `GlobalCellLedger` (`engine/ledger.py:394`) has an
`add_cell_traction` slot, so the booked side has a home; `owner_displacement_um` per connector row
does not exist.

**Recommendation to PI, two separable decisions:**
1. **Re-anchor the citation.** The defect text must be re-measured against the current connector
   runtime (`STATE.md`: 2 of 38 connectors device-run, 28 `NOT_BOUND`) or retired. A refusal whose
   evidence was deleted is a claim nobody can check either way — the same shape as the entries in
   `STATE.md` (c).
2. **Decide whether `load_path_audit` is produced.** It is three arrays per connector row and it
   would answer the *booked-vs-moved* question for **every** connector, not only traction's. That
   makes it the highest-leverage item in this document, and it is a connector-runtime change, not an
   observation one.

---

## 4. `stationarity` — status only, not touched

As instructed. Verified on this tree, 2026-08-20:

* Two implementations stand: `aleph/observe/stationarity.py` (1,030 lines, the observation stack —
  `fluctuation`, `kinematics`, `tether` all import it) and `aleph/engine/observe/stationarity.py`
  (the production GATE-B driver).
* `tests/architecture/test_stationarity_modules_agree.py` — **7 checks, all passing**.
* The **statistic is identical** (`tau_int` agrees to a ratio of 1.00000 on white noise, AR(1)
  ρ=0.9, ρ=0.99 and a drifting series). The **verdict contracts** disagree on **24.8%** of cases,
  down from 30.7% after the PI-adopted `min_tau_windows = 50.0`.
* `drift_sigma` (3.0 vs 1.0) is **HELD** and open to the PI in
  `ENGINE_FORWARD_ACCEPTANCE_2026-08-20.md` — applying it made agreement *worse* (74 → 89).

Consolidation is a separate task and nothing here was folded into it.

---

## 5. Where the line count now stands

The audit's "call sites" column and the column below are **not the same measurement** — the audit
counted call sites, this counts *files* under `aleph/scripts/**` + `aleph/engine/**` that import the
module — so the two are reported separately rather than differenced.

| module | lines | audit's call sites (2026-08-20) | importing files, now | note |
|---|---:|---:|---:|---|
| `stationarity` | 1,030 | 3 | 1 | untouched; dual implementation still pinned |
| `operator` | 493 | 3 | 1 | + the tether driver; gained the shared zero-`sem` guard |
| `tether` | 682 | **0** | **1** | wired, and reports `INPUT_ECHO` until W or the area moves |
| `manifest` | 432 | **0** | **1** | 34 fields now required, `is_complete()` called for the first time |
| `fluctuation` | 904 | 0 | 0 | **blocked on a PI decision**: no thermal term on the membrane |
| `kinematics` | 728 | 0 | 0 | `volume_flux` is wireable now; projected area does not exist |
| `traction` | 619 | 0 | 0 | **blocked**: cited defect files deleted; nothing produces `load_path_audit` |

Of the 3,365 lines the audit found unwired, **1,114 (`tether` + `manifest`) now have a caller**.
The remaining 2,251 are blocked on decisions, not on typing.
