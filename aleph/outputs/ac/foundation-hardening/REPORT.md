# FF-AC Codex foundation hardening — live FSI, transactional mechanics, and device retry

Status: **🚧 conservation, FSI transfer, transaction rollback, retry-partition invariance, exact-memory
instrumentation, NF2007 constraint projection, unilateral ERM mechanics, and source-gated Bell on/off CUDA
kinetics pass; sourced MCF7 ERM density/kinetic values, a single-ERM preload-force basis, physiological
equilibration, and the admin-disabled NG-9
measurement remain open**
Branch: `codex/ff-ac-codex`
Start commit: `765083d`
Native device: NVIDIA RTX A5000 Laptop GPU, Warp 1.14.0 CUDA

> **2026-07-21 scientific correction.** The earlier `13.102 Pa` / `2,297 ERM` capacity calculation below
> multiplied a continuum membrane tube-extraction force by a molecular linker count. That force is not a
> measured single-ezrin rupture force, so the calculation is retained only as historical diagnostic arithmetic
> and can no longer satisfy NG-3. The corrected gate requires a source-grounded single-ERM force basis plus
> MCF7 density and complete Bell kinetics. No previously reported runtime trajectory is relabelled.

## Outcome

This milestone replaces the main false-closure paths found in the composed Active Cell foundation:

1. static/host-remapped boundaries → persistent CUDA `wp.Mesh` membrane and nucleus classifiers;
2. voxel/L1 membrane area → live-triangle Kedem−Katchalsky surface quadrature with conservative Peskin
   deposition to the pressure grid;
3. uniform pressure stored only as a bulk field → live spatial membrane traction with a mask-aware affine
   interior pressure trace;
4. fixed inner-iteration count called stable → device-latched displacement, filament-constraint, finite-state,
   and force diagnostics;
5. failed mechanics left partially committed → device-only rollback of geometry, both pressure buffers, mask,
   `div(v_s)`, membrane source channels, moving-face content, and accepted physical time;
6. hidden per-step host diagnostics → device reports plus a directional memcpy profiler gate;
7. non-finite fluid state and invalid physical time could evade a transaction verdict → device finite latches and
   pre-mutation `dt_phys` validation;
8. irreversible ERM/nuclear rupture during candidate iterations → force-only candidate evaluation plus a
   device-predicated, accepted-transaction commit. A super-threshold ERM contributes zero candidate force
   without changing `bound`, so the committed topology cannot invalidate the converged force balance;
9. independent conservation checks with no scheduler-wide closure → one CUDA ledger integrating the true
   membrane surface flux, deposited membrane source, interior source, consumed `-alpha div(v_s)` source, and
   moving-face transport across every Biot subcycle of one accepted transaction;
10. unverified solid-to-fluid map → CUDA rigid-translation, affine-dilation, and Peskin partition-of-unity
    mass/momentum gates for `SolidDilatationCoupling`;
11. boundary-clipped bulk pressure gradient → affine-complete masked reconstruction plus one FLUID-masked,
    per-node-normalized Peskin stencil shared by fluid→solid interpolation and solid→fluid spread;
12. analytic sphere volume mixed with the actual triangular pressure surface → the live membrane's exact
    discrete enclosed volume is now the actin pressure/FSI control volume, closing internal force without an
    empirical correction or post-hoc COM projection.
13. single-attempt scheduler with no retry policy → same-time inner mechanics retry chunks plus a persistent
    CUDA run-active latch. Convergence disables later fixed-launch chunks, exhaustion rolls back once, and the
    first rejected outer transaction device-disables every remaining pre-scheduled physical step.
14. partial Warp mempool/endpoint memory readings with no exact route → a package-free NVML binding keyed by
    Warp's GPU UUID. It accepts only the driver accounting statistic `maxMemoryUsage` as exact, records bytes,
    PID, UUID, and mode, and keeps disabled/unsupported/missing values explicitly non-exact.
15. an unsourced one-ERM-per-membrane-node topology silently treated as physiological → an explicit areal-density
    build path whose count is derived from the represented triangular area, plus an upper-bound capacity contract
    that keeps NG-3 open when density/provenance is absent or mechanically insufficient;
16. NF2007 finite-step `reshape` used without its required force-side partner → a variable-length Warp-CUDA
    tridiagonal implementation of `P = I - Jᵀ(JJᵀ)⁻¹J`, launched on every inner iteration and validated against
    the independent dense NumPy oracle;
17. a generic “preload not initialized” diagnosis → force-family decomposition of the rejected full-native CUDA
    candidate, including a D2D pre-rollback snapshot and separate bending, crosslink, membrane, ERM, pressure,
    steric, and bulk-pressure measurements;
18. a two-sided Hookean ERM that became an unphysical compressive strut → a unilateral tether whose energy and
    force vanish at and below rest length while preserving the accepted-transaction rupture predicate;
19. any nonempty ERM source label could satisfy the NG-3 provenance condition → an explicit
    `erm_density_mcf7_production` latch that cross-cell diagnostics must leave false;
20. every retry chunk ended with an extra strong reshape operator → intermediate chunks now partition only the
    iteration budget, while the final attempt owns closeout reshape. A CUDA gate bounds the 3-iteration
    continuous-vs-partitioned position difference at `3.25e-22 µm`, far below the existing
    `7.45e-9 µm` position tolerance.

The moving-domain remap now reads an immutable old pressure field, interior and membrane water sources have
separate ownership, the pressure CFL responds to the live device pressure jump, and accepted solid motion feeds
`div(v_s)` on the refreshed domain. The CUDA operator also has native Terzaghi and Green acceptance gates.
The canonical accepted-step oracle now closes all five nonzero conservation terms in one transaction. It
keeps the staggered time index explicit: Biot consumes `div(v_s)^n`; the post-remap spread constructs
`div(v_s)^(n+1)` for the next step and is never back-dated into the current balance.

The result is intentionally **not** labelled a physiological resting baseline. The required population and
component laws execute on the A5000, but the 40 Pa initial state does not satisfy the registered mechanics
convergence criterion and is therefore rejected and rolled back.

## Gate verdicts

| Gate | Result | Evidence |
|---|---:|---|
| Live membrane−nucleus geometry | **PASS** | Rest and prescribed volume-preserving deformation: 0 query failures, 0 non-surface mask mismatches; mesh IDs and point allocations unchanged. Twelve resting surface ties are reported under the closed-set policy. |
| NG-2 constant state | **PASS** | Uniform 40 Pa drift = exactly 0 Pa. |
| NG-2 impermeable content | **PASS** | Relative content drift `4.91e-15` (fixed limit `1e-11`). |
| NG-2 source + `div(v_s)` balance | **PASS** | Relative error `2.58e-11` (fixed limit `1e-10`). |
| NG-2 live-triangle membrane flux | **PASS** | True triangle area `703.4902 µm²`; surface/grid deposition error `7.22e-15`; content-balance error `2.64e-11`; 0 unresolved faces. |
| NG-2 moving-domain remap | **PASS** | Relative closure error `5.56e-13`; remapped fluid cells remain exactly 40 Pa. |
| NG-2 Terzaghi | **PASS** | At `T_v=0.3`, degree-of-consolidation error `3.10e-5` and profile RMSE `7.28e-6`. |
| NG-2 Green spread | **PASS** | CUDA `d⟨r²⟩/dt = 299.999999999988 µm²/s` vs `2dc_v = 300`; relative error `4.11e-14`. |
| NG-2 component suite | **PASS** | All rows above pass on CUDA. |
| **NG-2 canonical coupled outer step** | **PASS** | One accepted 2-subcycle CUDA transaction consumes an affine `div(v_s)^n` made by the production FSI spread and has all five terms nonzero: content `Δ=-0.117144`, membrane surface/grid `4.57518e-7`, interior `0.0122625`, solid dilatation `-0.451032`, moving face `0.321625`; physical relative closure error `9.34e-14`, deposition error `7.16e-21` (fixed limit `1e-10`). |
| FSI dilatation transfer | **PASS** | Rigid-translation max divergence `3.89e-16 s⁻¹`; affine core divergence/integral errors exactly 0; Peskin spread weight error 0 and momentum relative error `7.40e-16`. |
| NG-3 uniform pressure surface law | **PASS** | CUDA/reference error `2.66e-15`; zero jump, sign reversal, closed-force, virial/work, and unresolved-face gates pass. |
| NG-3 linear pressure trace | **PASS** | Actual-surface affine pressure trace error `3.66e-15`; nodal force error `1.27e-15`; 0 unresolved faces. |
| NG-3 physiological t0 | **FAIL / OPEN** | The molecular capacity basis is absent: the historical 11.43 pN value is a continuum membrane tube-extraction scale, not a single-ERM rupture datum. No MCF7 active density or complete kinetic parameter set is registered. The 600 µm⁻² cross-cell proxy has its MCF7 provenance latch false and remains unconverged after both 6,000 and 30,000 iterations, so NG-3 correctly stays open. |
| NG-5 whole closed-system force/work | **PASS** | At the full 70,686-filament cortex, the affine bulk force matches `-α∇p V` to `3.42e-12`, the live membrane surface force matches `∮p n dA` to `1.17e-15`, and their net-force residual is `3.42e-12` (fixed limit `1e-10`). The actual CUDA Peskin pair gives nodal/grid work `2299.418011047/2299.418011027 pN·µm/s`, relative error `8.86e-12`, with positive dissipation. |
| Rejected outer-step transaction | **PASS** | Unconverged state restores `pos`, `p`, `p_new`, `mask`, `div(v_s)`, `s_membrane`, and `s_total` bit-exactly; non-finite fluid state rejects and restores pressure; invalid `dt_phys` raises before mutation. A forced ERM at `2.0000000000000013×` its derived rupture threshold contributes zero candidate force, leaves `bound` bit-exact through rejection, and ruptures exactly one tether only under the explicit accepted device predicate. |
| Device retry and multi-step fail-stop | **PASS** | Three retry chunks were pre-scheduled on CUDA. A first-chunk-converged solve executed 1 attempt and matched the one-chunk control bit-exactly; an unconverged solve executed all 3 attempts, restored all coupled fields bit-exactly, and produced outer enable states `[true, false, false]`. The accepted control committed all 3 steps at `0.001, 0.002, 0.003 s`. Partitioning 3 candidate iterations as 1×3 changes position by only `3.25e-22 µm`, versus the existing `7.45e-9 µm` tolerance. |
| NG-6 zero authoritative DtoH | **PASS** | Small and full-native hot loops: 0 `memcpy DtoH`; intentional post-loop `.numpy()` controls: 1. |
| NG-9 exact peak GPU bytes | **INCOMPLETE — ADMIN SWITCH** | The exact NVML process-lifetime counter is implemented and ABI/device-tested, but the A5000 reports `ACCOUNTING_DISABLED`. Population and isolated timing pass; Warp mempool high-water remains partial and is not relabelled exact. |

## Full-native population, timing, and partial GPU ledger

The authoritative attempt used the mandated cortical population; no biological density was reduced.

| Quantity | Value |
|---|---:|
| Unique active cortical F-actin | **70,686 fibers** |
| Actin nodes | **494,802** |
| Whole composed nodes | **511,114** |
| NMII particles / explicit heads | 15,028 / 8,840 |
| Nucleus / membrane nodes | 642 / 642 |
| Conservative field cells | 79,507 |
| Pressure/FSI control volume (live discrete membrane) | 1,751.9375 µm³ |
| Analytic seed-sphere volume (not used for discrete coupling) | 1,767.1459 µm³ |
| Warp mempool high-water (**partial, not exact whole-device peak**) | 258,522,474 bytes |
| Whole-device bytes observed after the attempt | 487,653,376 of 16,760,569,856 bytes |
| Exact-accounting backend / probe | `NVML_PROCESS_LIFETIME_ACCOUNTING` / **`ACCOUNTING_DISABLED`** |
| Accounted device / process | `GPU-ede1b117-f5f3-1bfc-f96b-a77f8fa87c2d` / PID 3,889,451 |
| Isolated attempted outer step / inner budget | **1.131 s / 1 attempt × 60 iterations** |
| Build / precheck / postcheck / total runner wall time | 7.166 / 0.067 / 0.067 / 8.433 s |

The earlier wording “exact Warp mempool high-water” was incorrect: the mempool counter is exact only for
allocations routed through that pool. `gpu_memory_accounting.py` now queries NVML's process-lifetime maximum
total allocation in bytes, which includes pooled and non-pooled allocations. It does not enable accounting or
request administrator authority. The hard exact value therefore remains open while the driver mode is disabled.

## Full-native mechanics finding

The rejected 60-iteration **candidate** remained finite and changed its maximum force residual from
`2634.119` to `1864.635 pN`. Its segment constraint residual was `1.44e-15 µm`, below the derived
`7.45e-9 µm` tolerance. However, the last projected displacement was `9.946e-4 µm`—about `1.34e5` times the
tolerance—so `inner_converged=False`.

That candidate is no longer mistaken for committed state. The scheduler rolled it back on CUDA:

- committed residual returned to `2634.119 pN`;
- committed COM drift and moving-face content were exactly zero;
- membrane/nucleus query failures and unresolved pressure/flux faces were zero;
- accepted physical time remained `0.0 s`;
- `outer_accepted=False`, `outer_rolled_back=True`, and `STABLE=False`.

This is the scientific result of the hardening pass: the original fixed-iteration/finite verdict hid a
non-equilibrated physiological initial condition. No tolerance, population, or physical magnitude was changed
to manufacture a pass. The native command explicitly used `max_inner_retries=0`, so its one-attempt result is
directly comparable to the prior baseline; retry count is an operator-visible numerical budget, not a hidden
parameter tuned to force convergence.

## Preload diagnostic correction and NF2007 projection finding (2026-07-21)

The historical acceptance algebra produced the following arithmetic, now explicitly quarantined from the
production gate because its force basis is a continuum membrane tether/tube scale rather than a single ERM:

| Quantity | Default diagnostic topology |
|---|---:|
| Discrete membrane area | 703.490216 µm² |
| Membrane-only Laplace support at 10 pN/µm | 2.666667 Pa |
| ERM count / represented density | 642 / 0.912593 µm⁻² |
| Continuum membrane tether/tube force scale | 11.434707 pN |
| Diagnostic multiplied pressure value | 13.101921 Pa |
| Diagnostic ratio to required 40 Pa | 0.327548× |
| Historical inferred density/count (**not molecularly authorized**) | 3.264914 µm⁻² / 2,297 |

These values do **not** establish a necessary molecular capacity bound. The Notion Contract-Graph/TAG query
found no registered MCF7 ERM areal-density, single-linker force basis, or complete kinetic datum, so the
production defaults remain `None` and require an explicit density, source label, MCF7-production provenance
latch, and Bell contract. A cross-cell source string alone can no longer close NG-3.
The density builder creates one independent ERM state per derived molecule, unique cortex endpoint, rest length,
and rupture flag; it does not multiply the force of 642 aggregate springs.

For diagnosis only, the primary-source estimate from
[Diz-Muñoz et al., PLOS Biology 2010](https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.1000544)
was tested at 600 linkers/µm². The paper reports that value for control zebrafish mesendoderm cells, not MCF7,
so the artifact is labelled `cross-cell-proxy_NOT-MCF7` and is not installed as a production default. Counting
against the actual triangular surface produced 422,094 explicit states (599.999816 µm⁻²). The former
`171.587×` capacity ratio is not molecular evidence and is quarantined. Its MCF7-production provenance latch
remains false. The run retains the PI-ratified single-linker
stiffness of 4,600 pN/µm, not the retired provisional 100 pN/µm. The provenance annotation was corrected:
[Braunger et al. 2014](https://pmc.ncbi.nlm.nih.gov/articles/PMC3975028/) reports a most-probable system
stiffness of 2.2 pN/nm, while the 4.6 pN/nm single-bond estimate is derived in Braunger's dissertation; the
value remains PI-ratified but its Contract-Graph registration is still pending. After 6,000 no-WCA projected
iterations on the A5000:

- all 422,094 ERM states remained bound; 317 compressed links now contributed exactly zero force, while the
  median/max tensile load was 0.1585/0.6681 pN versus the 11.4347 pN rupture force, with zero threshold crossings;
- force-family maxima were 50.814 pN pressure, 47.202 pN ERM, 3.397 pN membrane area tension, 1.146 pN
  crosslinks, and 1.546 pN after cancellation;
- the last NF2007 reshape moved any node by only 3.47e-15 µm, confirming that the new force projector removed
  the former projection/reshape split artifact;
- the strict projected displacement remained 4.33e-8 µm, `5.8109×` the unchanged 7.45e-9 µm tolerance; the
  exact projected force maximum was 0.3549 pN on an actin node, so
  the candidate was rejected. No tolerance or physical parameter was changed to manufacture a pass.

A follow-up 30,000-iteration diagnostic did **not** validate extrapolation from the early decay. Its strict
projected displacement was 4.13e-8 µm (`5.5466×` tolerance), raw residual rose to 2.409 pN, and exact projected
force maximum remained 0.3388 pN on an actin node. It accumulated 8,015 compressed zero-force ERM links;
the remaining tensile links stayed below rupture (maximum 1.4461 pN), with zero rupture events. The late
explicit trajectory is nonmonotonic and effectively stalled, so neither artifact closes NG-3. The 30,000-step
artifact predates the retry-boundary operator correction and is retained only as a stall diagnostic; the current
operator is separately covered by the partition-invariance gate above.

The unilateral force law and full Bell on/off **mechanism** are now implemented without inventing parameters:
ERM transmits tension but cannot push the membrane and cortex apart; accepted outer steps apply stochastic
Bell detachment and capture-gated rebind, while rejected steps are bit-exact no-ops. A successful rebind records
its formation length and therefore introduces no hidden prestress. The production parameter values remain
intentionally absent.
[Korkmazhan and Dunn 2022](https://pmc.ncbi.nlm.nih.gov/articles/PMC9355349/) provides single-molecule
slip-bond lifetime/distance evidence, and
[Fritzsche et al. 2014](https://pmc.ncbi.nlm.nih.gov/articles/PMC3907236/) measures cellular ezrin turnover,
but the assays do not supply a registered MCF7 areal density and complete on/off parameter set for this runtime.

The retained projector is not a parameter change. The per-fiber Warp kernel solves the exact symmetric
tridiagonal `JJᵀ` system for variable-length fibers, conserves total force, and matches the dense oracle within
2e-13 while giving `J(PF)=0`. The remaining blocker is a genuinely convergent GPU-resident implicit or
principled preconditioned solve plus PI-ratified MCF7 density, Bell values, and a single-ERM capacity
contract—not a looser gate or more
iterations of the stalled explicit mode.

## Remaining blockers

- Register or PI-ratify an MCF7 ERM areal density and complete Bell on/off contract. The 600 µm⁻² zebrafish
  mesendoderm estimate is a useful cross-cell proxy but cannot silently become the MCF7 production value;
  the explicit provenance latch now enforces that distinction.
- Complete a Warp-GPU implicit/preconditioned physiological preload solve that balances 40 Pa across membrane,
  ERM, constraint-projected cortex, and the remaining compartments at `t0`. The 30,000-iteration diagnostic
  falsified simple explicit-budget extrapolation, and the convergence tolerance remains unchanged.
- Enable the A5000's NVML accounting mode before launching a fresh native process (admin action:
  `sudo nvidia-smi -am 1`), then rerun the existing driver. The exact per-process lifetime counter is already
  wired; the Warp mempool and endpoint samples remain correctly labelled partial until the switch is enabled.
- The first whole-system production baseline still requires separately counted lamellipodium, SF/arcs/cap,
  filopodia, MT, IF, collagen ECM, α2β1−collagen clutches, and their explicit states. This milestone does not
  claim I5+ scope.

## Verification

- Dev Mac: **325 passed, 21 CUDA-only skipped** (`346` collected); no Warp-CPU simulation path is launched.
- RTX A5000: **346/346 `ffn_sim/tests/ac` tests pass** in 44.02 s.
- RTX A5000 focused ERM/preload/retry scope: **33/33 pass**; the retry-partition device gate passes.
- Focused Ruff scope for all changed foundation modules and tests: **PASS**. The repository still contains older,
  out-of-scope lint debt in untouched modules.
- Static contracts reject hot-loop `.numpy()`, Warp-CPU test launches, and hard-coded CUDA ordinals.
- The persistent assembled-state sidecar was regenerated with the transactional dumper: it records
  `inner_converged=false`, `outer_accepted=false`, `outer_rolled_back=true`, committed residual restored to
  `2634.119 pN`, and `stable=false`. The 2026-07-17 boot handoff carrying the former false `STABLE` verdict is
  prominently superseded and cannot serve as current evidence.
- Native artifacts:
  - `native_gates.json` — live-domain, transactional rollback, device retry/multi-step, NG-2/3/6, and NG-9 evidence;
  - `native_gates_arrays.npz` — figure source arrays;
  - `native_resting_70686.json` — full config, population, convergence, timing, memory, and profiler ledger;
  - `native_resting_70686.log` / `native_gates.log` — complete A5000 console records.
  - `preload_force_probe_proxy600_tension_only.json` — 6,000-iteration full-native, no-WCA cross-cell-proxy
    force-family decomposition, unilateral ERM distribution, device convergence history, capacity ledger, and
    exact projected-force diagnostics;
  - `preload_force_probe_proxy600_tension_only_30k.json` — rejected 30,000-iteration explicit-stall diagnostic,
    retained as evidence rather than production ratification.

## Figures

- **`figs/foundation_ng2_ng3_ng6.png`** — live-boundary parity; six NG-2 component residuals plus the coupled
  ledger, FSI spread momentum, and NG-5 force/work residuals against fixed limits; uniform and linear
  CUDA/reference pressure-force identity; full-native convergence residuals; directional DtoH trace with
  rollback, device-retry, and explicit NVML probe status; and the 150-vs-10 pN/µm preload gap. Axes start at
  zero except the explicitly labelled logarithmic convergence panel.
- **`figs/preload_capacity_and_relaxation.png`** — default-vs-proxy pressure-capacity ratios; 60-, 6,000-, and
  30,000-iteration displacement ratios against the unchanged gate; the late explicit convergence trajectory;
  and unilateral ERM tension quantiles with compressed zero-force counts against the derived rupture force.
  Log axes are explicitly labelled and the source-specific proxy warning is printed on the figure.

Regenerate the figure:

```bash
PYTHONPATH=. python ffn_sim/scripts/ac_foundation_vis.py
```
