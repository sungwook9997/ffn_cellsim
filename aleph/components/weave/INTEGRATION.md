# ac/weave — INTEGRATION notes (I4 unified weave -> lead-owned serial native integration)

Track **unified-weave** (`ac/weave`), increment **I4**. This file is the contract between this worktree's
NET-NEW `ac/weave/` package and the shared `ff/` runtime the **lead** edits in spine order (AC_PARALLEL
SESSIONS 1.3). Nothing here is edited into `ff/` from this worktree.

Contents: (1) the interfaces this track OWNS + the upstream APIs it CONSUMES read-only; (2) the `ff/`
patch-notes (exact file . symbol . before->after . double-count guard) the lead applies WHEN I4 lands native;
(3) the native-gate spec (gbook, lead-run); (4) honest-scope + forbidden-shortcut notes.

---

## 1. Interfaces — OWNED by this track + CONSUMED read-only

### 1a. OWNED (net-new `ac/weave/`)
- `woven_cell.py :: weave_cell()` / `WovenCell` — the ONE unified network builder. Regions concat leaves-out
  into GLOBAL indices; cross-region bonds are first-class (`xbond_*`). `to_crosslinked_cortex()` hands the
  (lead-owned) gamma/relax harness a `CrosslinkedCortex`. `population_ledger()` reports unique-active /
  dormant / node counts (build-plan A.1; no cortex/SF double-count).
- `branch_angle.py` (+ `branch_angle_warp.py`) — the 3-body angle-harmonic Arp2/3 junction
  (U = 1/2 k_theta (theta-theta0)^2, thermal sigma_theta = sqrt(kT/k_theta)). The host oracle matches the
  device kernel bit-for-formula.
- `lamellipodium.py` — the dendritic Arp2/3 REBUILD (mother/daughter tree, N-fixed dormant-daughter pool).
- `nucleation.py` — flux-limited Arp2/3 branch nucleation + capping (consumes the I1c monomer field read-only).
- `crosslink_kmc.py` (+ `crosslink_kmc_warp.py`) — the topology-reforming crosslink KMC + the single-channel
  double-count guard.
- `walk_dir.py` + `presets.py` — the I3 power-stroke `walk_dir` hand-off (host mirror + device kernel).
- `stress_fiber.py` — the SF prestress + traction-loop representation (FILAMENT_SUBSYSTEM_PLAN §3 G1/G2, the
  "biggest structural gap": stress fibers were seeded but not IN the whole-cell loop). `wire_stress_fibers(cell)`
  discovers each stress fiber LABEL-BLIND (a motorized + FA/LINC-anchored connected bundle — never `region_id`)
  and computes its EMERGENT axial prestress from the discrete NMII head reactions accumulated along the bundle
  backbone (NOT a lumped `k_SF`), plus the traction DIPOLE the SF hands its anchors — the SF's entry into the
  substrate loop (FA end) and the nucleus loop (LINC end, I7). Magnitude is GAP-guarded: `f_head_pN=None`
  returns the magnitude-independent SHAPE (axis, contractility sign, closed-vs-open dipole) with `tension=None`;
  a provisional per-head force yields a density-floored FINDING, never tuned to a traction band.

### 1b. CONSUMED read-only (frozen upstream interfaces, section 1.4)
- **I3 `ac/motor/hand.py`** — `NMIIHandParams` + `attach_prob`/`detach_prob`/`bell_off_rate` (`@wp.func`) +
  `allocate_hand_state`. `presets.py` CONSTRUCTS hands via this API; the crosslink KMC kernel reuses the same
  `@wp.func` closed forms. **This track never edits the KMC core.**
  - **walk_dir hand-off (HARD, from `ac/motor/INTEGRATION.md` 1).** `allocate_hand_state` carries a per-head
    `walk_dir` (vec3d), zero = passive (OFF default). **I4 fills it from the actin barbed-end polarity when a
    head (re)binds.** OWNED here: `presets.fill_walk_dir_kernel` (device) + `walk_dir.fill_walk_dir` (host
    mirror). The lead launches `fill_walk_dir_kernel` right AFTER the I3 `attach_kernel` each KMC tick so the
    power stroke walks along the ACTUAL actin grabbed (not the frozen minifilament-axis default). Provide
    `node_fiber` (N,) and `barbed_node` (F,) from `WovenCell` to that kernel.
- **I1c `ac/fluid/transport.py :: MonomerField`** — the G-actin monomer density `u = phi c`. `nucleation.py`
  consumes a sampled scalar/array `c` (the flux limit); it never runs the Warp field. `D_c` (~2-6 um^2/s draft)
  + `c0` (MCF7 GAP) block only the NATIVE nucleation gate, not the analytic build.
- **I2 `ac/nucleus`** — the deformable envelope geometry. The perinuclear cap's `linc_sites` are the cap-side
  tether endpoints; the nucleus-side endpoints are supplied by `ac.nucleus` at **I7** (read-only; the tether
  wire-in is I7's, gated with the cap driver + IF secondary net — the cap+IF double-load-path guard).

---

## 2. `ff/` patch-notes — lead applies these WHEN I4 lands native (do NOT apply from this worktree)

### 2a. `ff/weave.py` + `ff/architecture_spec.py` — DEMOTE the scripted bundle (build-plan 2b)

| Symbol | File | Before | After (I4) |
|---|---|---|---|
| `ff.weave._build_bundle` | `ff/weave.py:59-85` | pre-made parallel/graded sarcomeric bundle used as the FINAL SF/filopodium structure | **DEMOTE to a non-authoritative control.** `ac.weave.regions._build_emergent_bundle` seeds isotropic-conditional on the manifold; the bundle CONDENSES at I5. `_build_bundle` is reachable only via `RegionSpec(emergent=False)` (a labeled control). |
| `ff.weave._build_lamellipodium_patch` | `ff/weave.py:88-138` | scripted +-35 deg two-mode array (pre-aligned) | **REPLACE** with `ac.weave.lamellipodium.build_lamellipodium`: mother orientations isotropic-conditional on the NPF field; the +-35 deg two-mode EMERGES from theta0=70 deg branching, not scripted. |
| `ARP23_BRANCH_ANGLE_RAD` / `ARP23_BRANCH_SIGMA_DEG` / `ARP23_BRANCH_K` | `ff/architecture_spec.py:24-26` | Faessler-2020 branch anchors (kept) | **REUSE** — mirrored in `ac.weave.branch_angle` (theta0=70 deg, sigma_theta=9 deg, k_theta=0.173 pN.um/rad^2 = kT/sigma^2). No change; the ac oracle round-trips the equipartition relation the ff constant was derived from. |

### 2b. `ff/network_warp.py :: branch_angle_kernel` — port to the ac engine

`ac.weave.branch_angle_warp.branch_angle_kernel` is a NET-NEW port of `ff/network_warp.py:291
branch_angle_kernel` (same standard harmonic-angle force). It adds an `active` mask so the N-fixed dormant
daughter pool is force-neutral until a nucleation event activates it. The lead wires the ac kernel into the
new-engine inner force-assembly loop; the ff kernel stays in `ff/` untouched (read-only reference).

### 2c. `ff/motility_warp.py :: xl_turnover_kernel` — DISABLE (double-count guard, HARD)

**Double-count guard (build-plan 6, "I2-crosslinker-relaxation"):** there are TWO crosslinker-relaxation
channels and running both counts crosslinker fluidization twice.

| Symbol | File | Before | After (I4) |
|---|---|---|---|
| `xl_turnover_kernel` | `ff/motility_warp.py:357` | rest-length Maxwell creep `rest += frac*(L-rest)` (topology FIXED, stress relaxes) — the `r0_creep` channel | **DISABLE in the native driver** when the ac `kmc_reattach` channel is on. `ac.weave.crosslink_kmc.assert_single_channel` enforces EXACTLY ONE; `weave_cell(..., xl_reattach=True, r0_creep=False)` selects the topology-reforming KMC (which the r0-creep cannot do — needed for bundle condensation). The lead must not launch `xl_turnover_kernel` in the same run. |

`ff/ecm_mechanics.py:731-800` also drives `xl_turnover_kernel` on the ECM crosslinks — that is a SEPARATE
(ECM Maxwell) use and is NOT the actomyosin crosslinker channel; leave it (scope: cell actomyosin crosslinkers
only).

### 2d. `ff/clutch_*` / `ff/ecm_mechanics.py` — SF traction-dipole hand-off (build-plan G1/G2, lead wires at I6)

The `stress_fiber.wire_stress_fibers` output exposes each SF's FA-anchor nodes (`anchor_kind == "fa"`) as the
substrate-clutch endpoints and each cap SF's LINC nodes (`anchor_kind == "linc"`) as the I7 nucleus tether
endpoints. The **lead** wires the FA end into the existing `ff` α2β1–collagen clutch/ECM loop so the SF contractile
dipole (currently faked by the cortex shell gripping the substrate, §2.2 of the plan) becomes a real SF→FA→ECM→FA
loop — the I6 traction gate then requires BOTH channels (lamellipodial ratchet + SF sarcomeric, §5.3). This track
supplies only the anchor endpoints + the emergent axial force; it never edits `ff/clutch_*` from this worktree.
The SF axial force is assembled natively from the I3 `ac.motor` head kernels (no lumped `k_SF`).

### 2e. `ff/polymerization_warp.py :: G_actin` — monomer-field swap (shared with I1c)

The dendritic nucleation flux-limit consumes the I1c transported monomer field, not the fixed scalar
`G_actin`. This is the SAME patch-note I1c already records (`ac/fluid` INTEGRATION); named here so the lead
sequences the lamellipodium nucleation consumer after the field swap lands.

---

## 3. Native-gate spec (lead runs on the gbook A5000, in spine order I4 after I3)

The dev Mac cannot run Warp-CUDA (CPU-only build) — these gates are the lead's. Each maps to a CPU-green
analytic oracle in this package.

- **3a. Gate-1 full-population cortex parity (native, bit-identical).** `weave_cell([CORTEX])` all-OFF at the
  FULL 70,686 filaments -> bit-identical gamma to the current ff cortex. The CPU gate proved delegation
  bit-identity on a small cortex (population-independent); the native gate confirms it at full population +
  logs the population ledger (unique-active / dormant / nodes / peak GPU bytes). Oracle: `test_cortex_parity`.
- **3b. Branch-angle distribution (native).** On the native lamellipodium the device `branch_angle_kernel`
  relaxes to a branch-angle distribution matching theta0 +- sigma_theta (Faessler). The angle FLUCTUATES (soft
  potential); do NOT clamp. Oracle: `branch_angle` + `test_lamellipodium_oracle`. Device force ==
  `angle_forces` host oracle to tolerance (NG-0-style kernel parity).
- **3c. Dendritic protrusion + retrograde flow (native, report-not-tune).** Barbed ends push the membrane
  (Mogilner-Oster ratchet, reused); protrusion velocity + retrograde-flow band are FINDINGS (density-floored,
  build-plan 6.2) — NEVER tune `k_arp0`/`k_cap`/NPF to hit a band. Blocked on I0-B4 GAPs + the I1c monomer
  field live.
- **3d. Topology-reforming crosslink KMC (native).** Device bound fraction -> k_on/(k_on+k_off) unloaded
  (matches host); the bound-pair set reforms while the count is steady; crosslinker count conserved; OFF
  (k_on=k_off0=0) bit-identical to a static-crosslink regression. Oracle: `crosslink_kmc` +
  `test_crosslink_kmc_oracle`. **Confirm `xl_turnover_kernel` is OFF (2c guard).**
- **3e. walk_dir hand-off (native, HARD).** After each KMC tick the device `fill_walk_dir_kernel` sets bound
  heads' walk_dir from actin barbed-end polarity; a passive (unbound) head stays zero; GPU-residency
  zero-roundtrip (walk_dir never round-trips to host). Bit-identical to `walk_dir.fill_walk_dir`. Directed
  contraction emerges (the ac.motor powerstroke gate certifies the sign).
- **3f. Emergence PROOF (native, I5 — FALSIFIABLE, NOT this increment's claim).** With I1c + I3 + I4 live +
  myosin/KMC ON, the label-blind detector (`ac.emergence`) shows SF/cap/arc/filopodium bundles CONDENSE from
  the isotropic-conditional seed while cortex stays isotropic; controls (myosin-off / KMC-off / angle-off /
  unanchored) stay flat. **Flat result -> FINDING to PI + SEEDED scaffold fallback, NEVER re-tuned** (build-plan
  6.1). The region-id labels are DIAGNOSTIC and must not be read by the detector (firewall).
- **3g. Stress-fiber prestress + traction loop (native — SHAPE gate CPU-green, MAGNITUDE is GAP).** The SF is
  now IN the whole-cell loop (FILAMENT_SUBSYSTEM_PLAN G1/G2, the biggest structural gap). `wire_stress_fibers`
  discovers each SF label-blind (motorized + FA/LINC-anchored connected bundle) and produces its axial-prestress
  load path + traction dipole. CPU-green SHAPE gate (`test_stress_fiber_oracle`): an organized antiparallel
  FA<->FA bundle is a CLOSED contractile dipole (residual->0, uniform positive tension); the isotropic seed is
  OPEN (uncondensed) — the closed dipole is what CONDENSES at I5 (3f), not claimed by the seed. The native gate:
  with I3 head-resolved NMII live, the SF axial tension emerges from the discrete head reactions (reuse the
  `ac.motor` Hill kernels — NOT a lumped `k_SF`) and the FA-end dipole drives the substrate clutch (§2e); the SF
  MAGNITUDE (single SF ~5-6 nN) is a **density-floored FINDING**, report-not-tune, blocked on the I0-B3 per-head
  stall + I0-B6 engaged-head density (single SF magnitude is geometrically not back-solvable — build-plan §4C,
  [[project-sf-nmii-forcescale-result]]). Oracle: `stress_fiber` + `test_stress_fiber_oracle`.

### Native viz (lead, HTML — name the fields):
Interactive 3-D cell-morphology **HTML**: (i) **full region set** — cortex + ventral/dorsal SF + arc + cap +
filopodium + dendritic lamellipodium colored by DIAGNOSTIC region-id (post-hoc, firewall-safe); (ii)
**dendritic net** — mother/daughter tree with the angle-harmonic branch junctions + live crossbridges. Full-res,
no downsampling (PI 2026-07-07); browser-verified via `browser_check.py`. ANOMALY -> STOP, render, surface to
PI (never self-correct / never "explosion").

---

## 4. Honest-scope + forbidden-shortcut notes (state these; do not overclaim)

- **FULL fidelity vs scaffold (be precise).** FULL: cortex (bit-identical), lamellipodium dendritic Arp2/3
  (angle-harmonic branch + N-fixed activation + flux-limited nucleation), topology-reforming crosslink KMC,
  walk_dir hand-off. HONEST SEEDS (region builders, not yet emergence-proven): ventral SF, dorsal SF,
  transverse arc, perinuclear cap, filopodium — correct manifold + isotropic-conditional seed + anchors; the
  bundle CONDENSATION is the I5 native proof. DEFERRED: the myosin-turnover -> local-density SF condensation
  TRIGGER (build-plan 6.1), full growing-N nucleation (I9+), the cap -> nucleus tether nucleus-side wire-in
  (I7, ac.nucleus).
- **SF wiring (`stress_fiber.py`) — REPRESENTATION + hand-off delivered, MAGNITUDE is GAP.** The SF is now in the
  whole-cell loop (G1/G2): each SF is discovered label-blind and produces an emergent axial-prestress load path
  (from the discrete NMII head reactions, no lumped `k_SF`) + the FA/LINC traction dipole (§2d hand-off). The
  emergent-prestress SHAPE (closed dipole condenses, seed is open) is CPU-green; the single-SF force MAGNITUDE
  (~5-6 nN, density-floored) is a FINDING blocked on I0-B3/I0-B6 — never tuned to a traction band. The clean
  closed contractile dipole is what CONDENSES at I5, not claimed by the isotropic seed.
- **Regions are SEEDED, never scripted (P3 / 5-confirm-2/3).** Local orientations are isotropic conditional on
  the declared manifold + nucleator + NPF/polarity field; `ff.weave._build_bundle` is DEMOTED to a labeled
  control (`emergent=False`). The `+-35 deg` two-mode lamellipodium is NOT hard-set — it emerges from theta0=70
  deg branching. region_id is a DIAGNOSTIC label; the emergence detector must not read it (anti-coupling
  firewall).
- **N-FIXED (5-confirm-1).** Branch nucleation ACTIVATES pre-allocated dormant daughters (active mask) — it
  does NOT allocate new nodes. True growing-N nucleation/severing stays I9+.
- **No magnitude tuned to an outcome.** theta0/sigma_theta are Faessler-2020 sourced; k_theta is DERIVED by
  equipartition (k_theta = kT/sigma^2), never independently fit. Every nucleation/capping/NPF/fascin magnitude
  is an I0-B4 GAP surfaced to PI (`params_i0b4.yaml`); the native protrusion/bundle gate is INVALID until they
  close. Density-floored magnitudes are FINDINGS (report-not-tune, build-plan 6.2).
- **Single crosslinker channel (double-count guard).** `kmc_reattach` XOR `r0_creep`; enforced at construction
  (`assert_single_channel`). The lead disables `ff.motility_warp.xl_turnover_kernel` in the native driver.
- **Warp is SOURCE only on the Mac.** `branch_angle_warp` / `crosslink_kmc_warp` / `presets.fill_walk_dir_kernel`
  are authored + codegen-clean + CUDA-device-guarded but NEVER launched here (I0-A: Warp runs on CUDA GPU
  only). HOOMD is never imported (the warp-only-contract test enforces it in `ac/`).
