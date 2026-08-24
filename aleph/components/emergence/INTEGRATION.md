# emergence (I5) — INTEGRATION notes + DEFERRED harness spec

Per AC_PARALLEL_SESSIONS §1.3, this worktree edits **no `ff/` file**. The one `ff/` change this increment
implies is written below as a **patch-note + adapter** for the lead to apply in spine order. This file also
carries the DEFERRED work this session was scoped OUT of (the 4 ablation-control toggle hooks + the native
emergence proof) as an implementation spec, plus the **pre-registered decision rule** that must be LOCKED
before any native run (the anti-tuning record).

Scope reminder (Session E): detector-core + synthetic-config oracles ONLY. The ablation harness needs I3/I4
hooks that do not exist yet; the native proof needs I1c + I3 + I4 live. Both are deferred here, not built.

---

## 1. `ff/` patch-note — retire the VACUOUS `bundle_count`, route to the I5 detector

- **File:** `aleph/laws/architecture_metrics.py`
- **Symbol:** `bundle_count(net) -> int`  (currently `return int(net.n_fibers)`)
- **Before:** `bundle_count` "detects" a bundle by returning the number of fibers it was HANDED under the
  label "bundle" — it is `n_fibers` verbatim. It cannot tell a condensed bundle from an isotropic gas of the
  same fiber count, and it reads the construction label implicitly (you only call it on the thing you already
  built as a bundle). This is the `architecture_metrics.py` REBUILD flagged in NEW_ENGINE_BUILD_PLAN §2b
  ("`bundle_count`=n_fibers is VACUOUS — build a real condensation/nematic detector — I5").
- **After (adapter, no physics moved into `ff/`):** keep `parallel_order_parameter` (it is the validated S
  the detector reuses). **Deprecate `bundle_count` as an emergence metric.** Any caller that wants a bundle
  COUNT must call the label-blind detector:

  ```python
  from aleph.components.emergence.detector import EmergenceDetector
  report = EmergenceDetector().detect(net.pos, net.fiber_offsets)
  report.n_emerged_bundles      # honest O(1) count of CONDENSED patches (not n_fibers)
  report.condensed_fraction     # fraction of fibers inside any recovered bundle
  report.is_ordered_global      # global S beat the finite-N isotropic null at the pre-registered z
  ```

  Recommended concrete edit for the lead: leave `bundle_count` in place but add a one-line docstring
  deprecation ("DIAGNOSTIC ONLY — vacuous for emergence; use `ac.emergence.detector`") so no gate reads it as
  evidence of condensation. Do NOT move the detector into `ff/` — it lives in `ac/emergence/` and consumes the
  `ff` fiber-array contract read-only.
- **Double-count guard it satisfies:** the detector is the SINGLE source of the "did a bundle emerge?" verdict.
  `bundle_count`, `weave(bundle)`, and any region label must never be a second, parallel emergence claim
  (§5 confirm-(3): `weave(bundle)` is DEMOTED to a non-authoritative control). One verdict channel, label-blind.

No other `ff/` symbol is touched by I5. The detector imports **only** numpy/scipy and reads `pos` /
`fiber_offsets` — it never imports Warp or HOOMD (Warp-only contract, host acceptance layer).

---

## 2. Detector native-consumption interface (how it plugs into the live run — lead)

The detector is a **post-step host analysis**, not a hot-loop kernel. After the native `--from-resting`
mechanical solve settles at an outer-clock checkpoint, the lead reads the live fiber geometry off the GPU ONCE
(a single authorized device→host copy of `pos` + `fiber_offsets` for measurement — this is diagnostics, not
authoritative per-step state, so it does not violate the zero-roundtrip hot-loop rule) and calls:

```python
report = EmergenceDetector().detect(pos_host, fiber_offsets_host)
```

- **Input contract:** `pos (N,3) float`, `fiber_offsets (F+1,) int`, contiguous node ownership, ≥2 nodes/fiber.
  The native unified network (I4) exposes exactly this (it PREDATES the I4 region labels by construction).
- **Label firewall:** the lead MUST pass geometry only. Do NOT pass region/type arrays — the detector has no
  parameter for them, and `test_detector_contract.py` fails if one is ever added.
- **Scale note (native 70,686+ cortical fibers):** `compute_local_fields` loops per fiber in pure Python
  (one 3×3 eigen per fiber + a KD-tree ball query). At ~70k fibers this is a seconds-scale post-processing
  pass, acceptable for a checkpoint analysis. If it becomes a bottleneck at full whole-cell population, port
  the neighbor reduction to the shared device HashGrid (the same pattern I2b uses) — but the **analytic core
  here stays the reference** the GPU version is gated against (analytic-first).

---

## 3. PRE-REGISTERED decision rule (LOCK before the native run — anti-tuning record)

These are declared NOW, on synthetic ground truth, so they cannot be adjusted post-hoc to make a native run
pass (HARD: no gate-loosening, no param-tuning-to-outcome). Any change is a gate-contract change → PI sign-off.

| Knob | Value | Justification (Magic-Number Block) | Where |
|---|---|---|---|
| `EFFECT_SIZE_Z_CRIT` | **5.0** | 5σ discovery convention; dimensionless z, N-invariant; a bar the physics clears, not a fit | `null_model.py` |
| `NEIGHBORHOOD_SCALE` | **3.0** × median-NN spacing | neighbor-count sufficiency (O(10) neighbors); grid-invariant (scales with the config); recovery verified across the k∈{2.5,3,3.5,4} plateau | `condensation.py` |
| `ALIGN_COS_MIN` | **1/√2** (45° half-cone) | a member is aligned WITH the patch, not merely adjacent; geometric, dimensionless | `condensation.py` |
| `MIN_BUNDLE_NEIGHBORS` | **4** | a bundle needs several co-aligned fibers; structural sufficiency floor | `condensation.py` |
| null band | MC `(μ,σ)` at the run's N, cross-checked to `3/(2N)` | closed-form-anchored finite-N null | `null_model.py` |

**Pre-registered native decision (declared here, before the run):**
- **Ensemble:** ≥5 seeds per condition (§I5 "pre-registered ensemble seed count"); report per-realisation
  traces + mean (visualization-integrity), never a single seed.
- **Persistence:** a bundle counts as EMERGED only if it clears `Z_CRIT` **and persists** over a settling
  window (not a single-checkpoint transient). Window = the outer-clock interval over which the mechanical
  solve is converged; recorded with the run, not chosen after seeing it.
- **Positive verdict:** the treatment (myosin + KMC + anchoring ON) shows `n_emerged_bundles ≥ 1`
  (or `is_ordered_global`) **persistently and across seeds**, AND every one of the 4 controls (§4) stays flat.

---

## 4. DEFERRED — the 4 ablation-control toggle hooks (need I3/I4 — NOT built this session)

The detector is the READOUT; these controls are the FALSIFIERS. Each toggles OFF one candidate driver and
must make the emergence signal COLLAPSE toward the isotropic null — if a bundle "emerges" with the driver off,
the positive was an artifact. Verify finding: the hooks do not exist yet, so this is a spec.

| Control | Toggle | Needs | Expected detector signature (vs treatment) |
|---|---|---|---|
| **myosin-off** | set the head-resolved NMII active force to 0 (no contractile alignment) | **I3** `ac/motor` MyosinForce on/off | S / `n_emerged_bundles` collapse to the null band; STRESS_FIBER_TARGET says myosin-driven alignment is the driver |
| **unanchored** | remove the FA anchors (no boundary tension to organize against) | **I6** FA seed (or an I4/I1 anchor stub) | bundles fail to localize / drift; condensed_fraction → ~0 |
| **KMC-off** | freeze the crosslink KMC (no topology reattachment) | **I4** `ac/weave` crosslink KMC | no persistent bundle (order cannot be locked in); S transient only |
| **angle-off** | disable the angle-dependent crosslinker lifetime (§4D: emergence must pass with it OFF) | **I4** weave presets | emergence must STILL pass — this control proves the angle-choice is not secretly doing the work |

Harness shape (to build when I3/I4 land): one driver runs treatment + the 4 controls at matched seeds, calls
`EmergenceDetector().detect(...)` on each at the same checkpoints, and asserts
`z(treatment) > Z_CRIT ≥ z(control)` for all controls, per seed, with per-realisation traces rendered. Wire it
as `ac/emergence/ablation.py` + `tests/ac/emergence/test_ablation_native_spec.py` (native, lead-run on gbook).

---

## 5. DEFERRED — native emergence PROOF (needs I1c + I3 + I4 live — the falsifiable gate)

The actual P3 claim: seed an isotropic actomyosin network (filaments + NMII + α-actinin + FA, orientations
sampled isotropically conditional on the manifold/nucleator fields — NO pre-made bundle), run the coupled
solve, and show ordered+dense bundles CONDENSE where the SF/cap/filopodium candidate regions are, while the
cortex stays isotropic — read out by THIS detector, blind to the construction labels.

- **Prereqs:** I1c live monomer field (density feedback), I3 head-resolved NMII (the alignment driver), I4
  unified weave + crosslink KMC (topology reformation). Until all three are live, only synthetic proof exists.
- **Honest risk (hard-truth #1):** the `myosin-turnover → local-density-increase` trigger that
  STRESS_FIBER_TARGET names as the PRIMARY SF driver is not yet in the pre-allocated dormant-pool encoding, so
  emergence may NOT condense in the overdamped quasi-static solve. This is why the claim is genuinely
  falsifiable, not guaranteed by the architecture.

---

## 6. Falsifiability contract (HARD — §5 confirm-(5), hard-truth #1)

If the native run comes back FLAT — global S inside the null band AND `n_emerged_bundles == 0` persistently
across seeds with the treatment fully ON — that is a **genuine NEGATIVE, not a bug to tune away**:

1. **Surface it to PI as a FINDING** (with the per-realisation traces + the ablation table). Do not touch
   `Z_CRIT`, `NEIGHBORHOOD_SCALE`, the seed count, or any physical parameter to manufacture a positive.
2. **Fallback = a SEEDED scaffold, explicitly labeled `SEEDED`.** It validates force mechanics ONLY and must
   NEVER be reported as self-organization (§5 confirm-(3); the seeded `weave(bundle)` is a non-authoritative
   control, not evidence).
3. The detector staying willing to return "nothing emerged" is the whole point — a detector that always finds
   a bundle would be the vacuous `bundle_count` again.

---

## Change log
- 2026-07-16: created with the I5 detector-core. Patch-note (architecture_metrics.bundle_count retirement),
  native-consumption interface, pre-registered decision rule, and the deferred 4-toggle ablation + native-proof
  specs. Detector + synthetic oracle suite are built and green (44 host tests); the deferred items await
  I1c/I3/I4.
