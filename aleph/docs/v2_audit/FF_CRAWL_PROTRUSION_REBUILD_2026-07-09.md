# FF crawl — mechanistic protrusion rebuild (retire the body-force proxy)

**Date:** 2026-07-09 · **Owner:** Lead session · **Status:** PLAN (PI-approved to proceed)
**Supersedes the crawl-drive mechanism in** `_historical/FF_MOTILITY_CRAWL_2026-07-06.md` (the `leading_edge_push`/`protrusion_reaction` body-force pair).

---

## 1. Problem — the crawl was a tearing artifact, not translocation

The PI flagged the crawl as a failure. Diagnosis (proven on the coarse cell, host implicit,
fil 250, 1500 steps, dt 50 ms):

| condition | front-rear EXTENT | COM "drift" | verdict |
|---|---|---|---|
| **protrude OFF** (traction/contraction only) | 14.9 → **15.8 µm (Δ+0.8)** | +521 nm | ✅ COMPACT |
| **protrude ON** (current body-force proxy) | 15.0 → **115.5 µm (Δ+100)** | +18 883 nm | ❌ TEARS |

The cell does not crawl — the front flies forward (+71 µm) and the rear is pushed back (−6.5 µm);
the "COM displacement" is that asymmetric **stretch**, not the cell moving. With protrusion OFF the
cortex stays perfectly compact, which **proves the cortex is built correctly** (reshape + crosslinks
hold it) — the tearing is 100 % the protrusion.

### Root cause (two layers)
1. **The protrusion is a lumped BODY-FORCE proxy**, not real polymerization. `leading_edge_push_kernel`
   adds `+f_pro·phat` to the front-cap nodes; `protrusion_reaction_kernel` subtracts the total
   **uniformly over all cortex nodes** (`F_i −= total/Nc·phat`, `motility_warp.py:87`). Net force = 0
   (Newton pair), but it is a **stretching dipole**: front +, everything else −. A fixed elastic
   network under this dipole elongates without bound. This **violates the mechanistic-first principle**
   (a proxy in place of the real barbed-end polymerization we already built).
2. **The implicit crawl loop also dropped the NF2007 §5.3 `reshape` segment-inextensibility projection**
   (device arrays `foff_d/soff_d/sr_d` were set up but never launched). Wiring it back helps
   (136 → 115 µm) but is **necessary, not sufficient** — the dominant tear is the body-force dipole.

## 2. Blast radius — verified ISOLATED from AFM / resting / DCM

6-agent audit (workflow `wf_9ed77723-aac`). The rebuild edits exactly two things:
**`aleph/laws/motility_warp.py` protrusion kernels + the `ff_crawl_on_substrate.py` crawl loop.**

| component | uses protrusion kernels? | uses crawl loop? | verdict |
|---|---|---|---|
| **AFM** (ff_hertz_validation, ff_press_speed_probe → `simulate_whole_cell_compression_on_device`) | no | no | **ISOLATED** |
| **resting / adherent** (gamma_floor, `--from-resting`, resting full-compartment) | no | no | **ISOLATED** (crawl→resting is one-way) |
| **DCM** (separate engine, imports nothing from `ff/`) | no | no | **ISOLATED** |
| motility reverse-deps | crawl-ONLY | crawl-ONLY | PARTIAL — *file* hygiene only |
| tests | none assert protrusion | — | PARTIAL — *import* hygiene only |

**Answer to the PI's question — this does NOT touch AFM (or resting/adherent/DCM).** AFM uses a
completely separate path: `build_crosslinked_cortex` (the shared cortex *builder*, unchanged) and
`simulate_whole_cell_compression_on_device` (its own compression sim). Neither the AFM scripts nor
gamma_floor/network_warp ever import the protrusion kernels or the crawl loop, and no shared module
carries mutable state (all module globals are immutable constants / frozen dataclasses / stateless
kernels). The protrusion kernels are consumed by **`ff_crawl_on_substrate.py` alone**.

**Two care-points (file hygiene, not physics):**
- `motility_warp.py` also exports generic helpers used by 2 non-crawl demos
  (`axpy_physical_kernel` → ecm_remodel; `cortex_volume_kernel`/`physical_node_gammas`/`volume_gradient`
  → membrane_cell). Do NOT break those or `__all__`.
- If a kernel symbol is deleted/renamed, update the crawl driver import line (`ff_crawl_on_substrate.py:39-42`)
  + `__all__` (`motility_warp.py:295`) in lockstep, or 2 crawl tests fail at *collection* (ImportError).
  No test asserts protrusion physics.

## 3. Mechanistic design — real polymerization, emergent reaction

Replace the body-force dipole with the actual actin-treadmill physics we already have the pieces for:

```
FRONT: directed barbed-end polymerization  →  reshape advances the barbed node  →  pushes the front membrane out
         (grow front-cap, forward-oriented fibers' tip rest length; Mogilner-Oster ratchet throttles by membrane load)
REACTION: emergent, NOT a body force  →  link_spring propagates the tip reaction down the filament
         →  into the basal network  →  bound clutches resist retrograde flow against the substrate  →  TRACTION
REAR: pointed-end depolymerization (treadmill)  →  removes subunits at the back  →  mass-conserving  →  cell stays compact
BODY: pulled forward by the net forward traction; held compact by cortical tension + hard V=V0
```

This is mechanistically honest: the only external force on the cell is substrate traction through the
clutches (F*≈7 pN Pereverzev catch-slip, as-recorded per PI); polymerization is internal; the reaction
is carried by the network, not smeared as a body force. The existing `--audit` gate (clutches OFF ⇒ net
COM drift ≈ 0) is exactly the correctness check for "polymerization alone cannot translocate."

### Reuse (all validated, GPU-resident) vs build-new

**REUSE:** `resolve_polymerization` (Pollard/Mogilner-Oster constants) · `polymerization_kernel`
(`polymerization_warp.py` — already takes a **curated** barbed_node/prev_node/barbed_seg index list +
per-filament `load_share`, exactly what directed front growth needs — preferred over editing
`barbed_end_growth_kernel`) · `reshape_kernel` (grown rest length → forward advance) · `link_spring_kernel`
(tip reaction → basal network) · `clutch_spring_kernel` + `clutch_catchslip_kmc_kernel` + fa_maturation
stack (emergent load-gated traction) · `xl_turnover_kernel` (keeps cortex fluid).

**BUILD NEW:** (a) host **front-fiber selector** — filter `S['front']` cap fibers by
`tip_tangent·phat > cos θ_max` → emit the polymerization_kernel index arrays (directed front growth);
(b) **rear pointed-end depolymerization** kernel — mirror of barbed growth, shrinks rear tip-segment
rest length (mass-conserving treadmill); (c) **bead-insertion** topology op to sustain front growth past
`seg_max = 2·seg0`; (d) **DELETE** `protrusion_reaction_kernel` (+ retire `leading_edge_push_kernel`).

## 4. Staged implementation + validation gates

Each step is gated; do not advance on a red gate.

**Step 0 — keep the reshape wiring** (already done, working-tree): reshape in both implicit paths.
Necessary infra; retained.

**Step 1 — CORE: retire body-force, directed front polymerization + emergent reaction.**
Delete `protrusion_reaction_kernel` smear; replace the `elif protrude:` launch with directed front
barbed-end polymerization (front-cap forward-oriented fibers) + reshape advance. Front growth is capped
(`seg_max`), so the front advances a bounded amount while traction pulls the body.
- **Gate 1a (compact):** extent Δ ≤ ~1 µm over the run (vs +100 µm now).
- **Gate 1b (translocation):** front AND rear COM both advance forward together (not front-only).
- **Gate 1c (traction-driven):** `--audit` clutches-OFF ⇒ net COM drift ≈ 0.

**Step 2 — rear treadmill** (sustained crawl): rear pointed-end depolymerization → mass conservation.
- **Gate 2:** the cell crawls > 1 cell-length while extent stays ~const (no monotonic elongation).

**Step 3 — bead insertion** (unbounded): topology op past `seg_max`.
- **Gate 3:** sustained multi-length crawl, seed-robust (seeds 7/11/17/23), physiological speed band.

**Step 4 — NATIVE re-test** (the deliverable): Nc≈266k on A5000. Compact + physiological.

## 5. Out of scope (separate open items)
- **Native crawl SPEED / whole-cell Stokes drag** (Σγ ∝ Nc): the `--com-drag` modal-drag fix is a
  *solver-conditioning* item, orthogonal to this polymerization rebuild. Track separately; do not
  conflate with the tearing fix.
- `build_crosslinked_cortex`, `simulate_whole_cell_compression_on_device`, gamma_floor, DCM — untouched.

## 6. Commit discipline
One coherent commit per gated step (reshape + Step 1 together as the "retire body-force proxy" fix),
kb-check clean, do not push `ffn/foundation` without PI sign-off.
