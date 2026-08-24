# I2 nucleus — INTEGRATION.md (deferred ff/ + fluid patch-notes, native-gate spec)

Per AC_PARALLEL_SESSIONS §1.3: where an increment must retire/extend an `ff/` mechanism or reconcile with
another track's owned file, the worktree does NOT edit it — the change is staged here as a documented
patch-note + adapter, and the **lead applies it in spine order**. This track edited ZERO `ff/` files and
ZERO `ac/fluid/` files; everything below is a hand-off.

---

## 1. `ff/nucleus_envelope.py` — REUSE-WITH-FIX → SUPERSEDED by `ac/nucleus/` (patch-note, do NOT edit in worktree)

`ff/nucleus_envelope.py` is the partial scaffold (build-plan §2b "REUSE-WITH-FIX"): it has a SINGLE global
bilinear areal tension + rupture, reuses `build_membrane_mesh`, but has **no** Helfrich bending kernel
wired, **no** volume/nucleoplasm, **no** chromatin, **no** LINC, and **no** lamin-A/C-vs-lamin-B split.

**Before → after (the lead applies at the I2 native-integration step):**

| Symbol / mechanism | Before (`ff/nucleus_envelope.py`) | After (`ac/nucleus/`) | Double-count guard |
|---|---|---|---|
| bending | declared `KAPPA_NE_PN_UM` but **no kernel** | reuse `ff.membrane_surface.helfrich_bending_kernel` + `calibrate_kappa_tilde` (κ̃=8πκ/Σ_ref) | bending is a NEW channel (was absent) — no overlap |
| areal tension | ONE global σ, single bilinear `RATIO_LAMIN=3.0` | PER-FACE local strain → framework-#6 3-regime split (`lamina_areal_tension_kernel`) | replaces the global-σ scaffold; do not run both |
| rupture | per-face `ruptured` flag on global σ | per-face local-strain rupture (`rupture_update_kernel`) | same flag semantics; the strain source changes global→local |
| volume/nucleoplasm | **absent** (bead ball had radial spring) | incompressible ν→½ penalty (`nucleoplasm_volume_*_kernel`) | NEW; replaces the retired `nucleus_shell_kernel` radial spring |
| chromatin | **absent** | WLC internal net (`chromatin_wlc_kernel`) | NEW |
| LINC | **absent** | nonlinear tether (`linc_tether_kernel`) | NEW; the I7 load path |
| nucleoplasm viscosity | **absent** | `nucleoplasm_viscosity_kernel` | NEW |

**RETIRE (build-plan §2b "CANNOT USE"):** `network_warp.py::nucleus_shell_kernel` (radial bead ball — "no
internal connectivity, one deformation mode") is DELETED in the same commit that wires the `ac/nucleus/`
deformable mesh into the assembly (same-commit guard, so no run ever has both the bead ball and the mesh).
The `E_NUC_PA=399` / `RATIO_LAMIN=3.0` / `EPS_RUPTURE=0.50` constants in `ff/nucleus_envelope.py` are
**scaffold proxies** — the authoritative values now live in `params_i0b2.yaml` (E_nuc anchored, the split
+ rupture GAP-PI). Do not import the ff constants as if sourced.

**Population-ledger (build-plan §A.1):** the nucleus registers unique IDs — envelope nodes/edges/faces
(subdiv=3 → 642 nodes / 1280 faces / 1920 hinges) + explicit chromatin-polymer link IDs — NOT hidden in a
radial spring or a bulk modulus. Report N_nodes/N_faces/N_hinges + chromatin-link count at the native gate.

---

## 2. `ac/fluid/domain.py` — the §1.4 no-flux reconciliation (fluid-spine OWNS it; adapter, do NOT edit)

`ac/fluid/domain.py :: Domain.set_nucleus_boundary(provider)` does not exist yet (fluid-spine is at its I1a
analytic stage). This track ships the conforming provider and codes against the frozen signature:

- **We ship:** `mask_provider.DeformableNucleusMaskProvider(verts, faces)` with `inside_mask(grid_points)
  → bool[N]` (winding-number occupancy of the LIVE mesh) + `boundary() → (centroids, outward_normals)` +
  `update(verts)` (refresh each outer step). It satisfies the `NucleusBoundaryProvider` Protocol.
- **I1a ships (placeholder):** a STATIC-sphere provider. Our `StaticSphereMaskProvider` mirrors it so the
  contract is testable now; when fluid-spine lands the real `Domain`, it calls `set_nucleus_boundary(
  DeformableNucleusMaskProvider(...))` — a drop-in swap, no co-edit.
- **The reconciliation the lead must enforce (build-plan I2 row):** the I1a Biot no-flux mask must track
  the **LIVE deformable oblate mesh about `centre_nuc`, NOT a static sphere R_nuc**. The provider is
  `update`-d from the mechanical solve's envelope nodes every outer step so the mask moves with the
  surface. The **relative** no-flux BC (q·n̂=0 in the moving-boundary frame) needs the envelope node
  velocities — the fluid track reads them from the shared state at integration; the provider supplies the
  geometry (faces + outward normals), the fluid solver supplies the velocity coupling. Single source of
  the boundary: the `ac/nucleus/` mesh, consumed read-only by `ac/fluid/`. Neither edits the other.

---

## 3. Inner force-assembly (the §1.4 `accumulate` contract; lead-owned integrator sums)

Each nucleus force primitive is a Warp kernel adding its per-node contribution to a shared `out_force`
(mirrors `PressureCoupling`/`MyosinForce`/`StericForce`). The lead-owned integrator launches, in the inner
mechanical solve, on the nucleus node block: `helfrich_bending_kernel` (reused) + `lamina_areal_tension_
kernel` + `nucleoplasm_volume_force_kernel` (after a `nucleoplasm_volume_reduce_kernel` → host p_vol) +
`nucleoplasm_viscosity_kernel` + `linc_tether_kernel` + `chromatin_wlc_kernel`, then `rupture_update_
kernel` once per OUTER step (irreversible tear). Add κ̃/ℓ³, K_areal, K_vol, γ_nuc/dt, k_linc, K_chrom to
`kmax` for the CFL (Sanity Gate). All float64, atomic_add, one-thread-per-element — no per-step host state.

**Single-channel invariants (hard-truth #6, cross-track):**
- The nucleus↔cytoskeleton load path is EITHER the LINC tether (I7 actin cap → LINC → nucleus) — there is
  no second radial-strut nucleus anchor. (The retired IF-cage polar spokes do NOT return as the driver;
  IF re-enters only as the §3A.c cell-type secondary passive net, gated with I6.)
- Reconcile with **I2b (excluded-volume)**: the solid-EV track flags that its all-fiber EV must reconcile
  with I7's LINC/plate soft-contact. The nucleus envelope's steric interaction with cortex filaments is
  the I2b `StericForce` (one steric channel), NOT a second nucleus-specific contact kernel.

---

## 4. Native-gate spec (lead runs on gbook, in spine order — this track authored + CPU-greened it)

The dev Mac cannot run Warp-CUDA (CPU-only build). The lead runs these `--from-resting` native gates on the
A5000 after I1a's live domain lands (the nucleus mask consumes it):

1. **Resting stability** — the deformable nucleus is stable (no drift/blow-up) at the physiological E_nuc
   (399 Pa MCF7 anchor) under the resting Π₀=40 Pa turgor; envelope bending energy ≈ 1.0×8πκ at rest
   (κ̃-calibrated), volume conserved to machine precision.
2. **Oblate flatten at the resting adherent baseline** (I7, after cap/LINC) — the nucleus flattens to the
   MCF7 oblate aspect (GAP-PI, 1.5–3 band) with **volume conserved** (a∝A^⅓); cap-OFF stays spherical.
   ⚠ Force-scale is a REPORTED finding (K_nuc≈25 nN/µm → 5–6 nN cap ≈0.2 µm ≈10× short) — do NOT tune.
3. **Reconcile the no-flux mask against I1a's live domain** — the fluid mask == the live oblate mesh
   occupancy (not a static R_nuc sphere); global fluid-content balance still closes with the moving
   nucleus inclusion.
4. **Rupture emergence** (magnitude INVALID until eps_rupture closes) — rupture appears only past the
   sourced envelope threshold; structural on/off is CPU-greened, the native magnitude is a PI-blocked
   finding.
5. **Zero GPU↔CPU roundtrip** inside the physical-time loop (the I0-A device-residency profiler gate).
6. **Viz** — the interactive 3-D HTML in `outputs/ac/nucleus/REPORT.md` §Native-gate viz spec.

**Magnitude gates are INVALID until `params_i0b2.yaml` closes** its 5 GAP-PI unknowns (lamin-A/C-vs-B
split, eps_rupture, eta_nucleoplasm, k_linc, MCF7 oblate aspect). The STRUCTURAL gates (§ REPORT table)
are magnitude-independent and pass now.

---

## 5. I0-B2 GAP surface (PI decisions needed before the native magnitude verdicts)

1. **lamin-A/C vs lamin-B modulus split** — the 3-regime STRUCTURE is verified (KB-3.B2.2); the numerical
   partition of E_nuc(399 Pa) into K_chrom / K_lamin_b / K_lamin_ac is unsourced. Do NOT reuse the ff
   `RATIO_LAMIN=3.0` proxy as if sourced.
2. **E_nuc method-dependence** — 399 Pa (MCF7 in-situ, Fischer2020, low-strain/chromatin regime) is the
   anchor; 1–10 kPa (KB-3.B2.1) is the isolated/large-strain band. Confirm the low-strain anchor + whether
   the lamin-A/C tangent should be pinned to the kPa isolated value.
3. **knee strain-measure** — KB-3.B2.2 gives the knee as an EXTENSION (~3 µm); the model uses areal strain.
   Ratify the conversion (do not silently keep 0.10).
4. **envelope-rupture criterion** — tension-governed vs strain-governed (framework asks tension; KB-3.B2.3
   gives a ~10% cross-section CONFINEMENT/compression criterion). Ratify criterion + MCF7 value.
5. **nucleoplasm viscosity** — 1–100 Pa·s band (KB-DRAFT-3.B-03); ⚠ keep SEPARATE from the 65.9 Pa·s
   cytoplasm bulk drag (the same double-count trap the fluid track guards for mu_pore). Pin the MCF7 value.
6. **k_linc** — HALT→PI (KB-DRAFT-3.B-09: nesprin nonlinear, no Hookean k; ~8 pN is a TENSION not a
   stiffness). Needs the nesprin force-extension law + magnitude before the I7 flatten magnitude gate.
7. **MCF7 oblate aspect** — KB-6.2.2 has MDA (2.2–4.7) but not MCF7; the flatten is REPORT-not-tune (let it
   emerge from the cap load), not a set input.
