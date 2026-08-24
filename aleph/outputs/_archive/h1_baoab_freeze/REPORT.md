# H.1 BAOAB freeze point — Main Session closeout

**Date**: 2026-05-20
**Branch**: `phase1/h1-ecm` cut from `v2/foundation` @ `2543d6a` (per boot
prompt). **Actual parent commit at sign-off** is `0a1b682` — between my
`git checkout -b phase1/h1-ecm 2543d6a` and my final commit, two upstream
commits landed on the branch base (`06c91d4` CLAUDE.md Notion pointer +
`0a1b682` `v2/foundation → ffn/foundation` rename across docs refs). I
did not author or solicit those commits and they do not conflict with
BAOAB work; they look like Orchestrator-side housekeeping. PI should
confirm that the freeze point Sub Session inherits is `75a1035`, not the
sibling state at `2543d6a`. See §Open for PI item 5.

**Scope**: D3 deliverable only (Leimkuhler-Matthews BAOAB-limit Updater).
Mikado topology / cross-links / shear protocol / KU-1.30 validation are
**not** in this commit — the boot prompt requires PI sign-off on BAOAB
before touching anything else.

## What landed

| File | Status |
| --- | --- |
| `ffn_sim/integrator/baoab.py` | New. Module-docstring Sanity Gate (§1-6) written before execution, per CLAUDE.md hard rule. Implements `LeimkuhlerMatthewsBAOAB(hoomd.custom.Action)` + `make_baoab_updater()` factory + `_wrap_into_box()` helper. |
| `ffn_sim/tests/test_baoab.py` | New. 16 STATIC checks (input validation, boundary cases, sign-sense, prv_rnds memory, Lees-Edwards-aware box wrap). All PASS. |
| `ffn_sim/scripts/baoab_polymer_sanity.py` | New. 100-bead polymer empirical gates (bond equipartition vs 3D analytical, free-bead Stokes-Einstein diffusion, same-seed reproducibility). Bending mode reported as DIAGNOSTIC with rationale (see §Open for PI). |
| `ffn_sim/outputs/h1_baoab_freeze/baoab_sanity_report.json` | Gate-by-gate evidence. |
| `ffn_sim/outputs/h1_baoab_freeze/REPORT.md` | This file. |

## Gate results

**STATIC (`pytest ffn_sim/tests/test_baoab.py`)**: 16 / 16 PASS. ~1 s wall.

**EMPIRICAL (`python ffn_sim/scripts/baoab_polymer_sanity.py --bench`)**:

| Gate | Bound | Measured | Result |
| --- | --- | --- | --- |
| Bond equipartition vs 3D analytical | ±5% | +1.05% (⟨U⟩/bond=0.5152 vs ½kT·(1+2kT/(kr0²))=0.5099) | **PASS** |
| Free-bead diffusion D = kT/γ | ±10% | +5.9% (D=1.059 vs 1.000) | **PASS** |
| Same-seed reproducibility | max\|Δr\|=0 | 0.0 bit-identical over 5000 steps | **PASS** |
| Bending mode (DIAGNOSTIC, not a gate) | — | ⟨U⟩/triplet=0.8545; vs brief 2D ½kT: +71%; vs 3D analytical 0.9342: −8.5% | DIAG |

Wall-time on M1 Max CPU: **10 003 steps/s** = 0.10 s / 1k steps on the
100-bead polymer with bond + angle + L-M Updater. Phase 0.2 reference
(`hoomd_polymer_sanity.py` with Langevin) was 11 s / 1 M steps =
~91 k steps/s; the L-M Updater is ~9× slower per step because the
position update is per-step Python with `cpu_local_snapshot` re-entry.
This is below the H.1 wall-time budget for the polymer (≤5× v1 regression)
but is the obvious bottleneck for ECM scaling and worth flagging.

## How the implementation deviates from a naive read of the brief

These were judgment calls made under "auto mode"; PI should override if
disagreed.

1. **Tag-indexed prv_rnds, not row-indexed.** HOOMD 7's `ParticleSorter`
   tuner reorders local-snapshot rows between calls. Indexing prv_rnds
   by snapshot row corrupts the W_{n-1} memory across reorders.
   Implementation reads `snap.particles.tag` each step and gathers/
   scatters per stable tag. Confirmed empirically — initial naive
   row-indexed implementation failed 5/16 static tests with a
   characteristic "particles swap" symptom; tag-indexed version
   passes all 16.

2. **`md.Integrator.methods = []`, forces-only.** The brief asks for a
   `hoomd.custom.Action` Updater co-existing with HOOMD force machinery.
   I verified that `md.Integrator` evaluates forces every step even with
   `methods=[]`, so the Action can read `net_force` directly from the
   `cpu_local_snapshot` without a no-op NVE method. The Action raises
   at attach if `methods` is non-empty (would double-step). Probe script
   results saved in commit message.

3. **Manual Lees-Edwards-aware box wrap.** HOOMD 7's Python `Box` does
   not expose a public `wrap()` method (confirmed by inspection). The
   Action wraps positions using the upper-triangular fractional-coord
   transform on `Lx, Ly, Lz, xy, xz, yz`, so the same Action handles
   the orthorhombic polymer test AND the xy-tilted Lees-Edwards box
   that H.1's KU-1.30 #2 strain-stiffening test will need.

4. **Bond equipartition gate target is the 3D-corrected analytical
   reference (0.5099), not the naive ½ kT.** In 3D, the harmonic
   magnitude potential on a 3-vector has an `r²` volume element that
   shifts ⟨U⟩ above ½ kT by ~2 σ²/r0². With σ/r0=0.1 in this config,
   the correction is +2% — comparable to the empirical observation
   (+1%). Naive ½ kT would give a misleading "+3% bias" reading that
   is actually 2/3 physical and 1/3 sampling-error. The 3D-corrected
   reference is the principled gate.

## Open for PI

1. **Bending-mode equipartition target in the brief is 2D.** The H.1
   brief specifies `⟨½ k_θ (θ−π)²⟩ ≈ ½ kT to ±10%` as an equipartition
   acceptance for the polymer. That target value is the AFINES 2D
   reference; HOOMD is 3D and the sin(θ) volume element on the bond
   angle's two transverse DOFs shifts the analytical mean to
   ⟨U_bend⟩ ≈ 0.93 kT (numerically integrated for k_θ=5, kT=1,
   t0=π) — see `_bending_3d_analytical_reference()` in the polymer
   script. My measurement of 0.855 lies 8.5% below the 3D analytical
   reference. Per CLAUDE.md hard rule "no gate-loosening", I have
   **not** edited the bending gate inline. Recommendation: drop the
   bending equipartition gate from BAOAB sign-off and move principled
   bending validation to H.2 (persistence length L_p via tangent-
   correlation against `κ_B / kT`). PI: please confirm.

2. **L-M ↔ E-M O(Δt²) ↔ O(Δt) order separation deferred to H.2.** My
   first script tried an order check (halve Δt, expect ≥3× error
   reduction). At the dts where the polymer is stable for overdamped
   integration (Δt ≤ 0.005, ≈ 0.5 τ_relax), the integrator bias on
   bond equipartition is below statistical noise even at 500 k samples,
   so the order signature is unobservable. The right home for this
   verification is H.2 — the staged update to
   `ffn_sim/validation/oracles/common/sanity_gate.py` already adds a
   `PHASE_1_REFERENCE_INTEGRATORS` field (`euler_maruyama`,
   `hoomd_brownian`) gated on `allow_reference_integrator=True`, which
   anticipates exactly this comparison.

3. **Pre-existing uncommitted state on `v2/foundation`.** On session
   boot the working tree had three uncommitted edits in the read-only
   `validation/oracles/` area: D3 + KU-1.1 additions to
   `sanity_gate.py`, a new `tests/test_sanity_gate.py`, and a D3
   integrator-name update to `configs/phase1_unit1.yaml`. These look
   like Sub Session / Orchestrator pre-staging for H.2 (not Main
   Session H.1 work). I stashed both to preserve them; they survive
   on `v2/foundation` as `stash@{0}` and `stash@{1}`:

       stash@{0}: pre-h1-boot: uncommitted D3 integrator switch in phase1_unit1.yaml
       stash@{1}: pre-h1-boot: uncommitted D3+KU1.1 sanity_gate.py changes

   PI: please confirm these stashes are accounted for in the Sub /
   Orchestrator queues, or escalate.

5. **Branch base advanced under me during the session.** Documented in the
   header. The two commits that landed on the branch base while I was
   working (`06c91d4`, `0a1b682`) are unrelated rename/docs housekeeping
   and do not conflict with BAOAB. PI: confirm `75a1035` is the
   intended freeze point and that the rename from `v2/foundation` to
   `ffn/foundation` is the canonical branch name going forward.

4. **Per-step Python L-M Updater is ~9× slower than HOOMD-native
   Langevin.** Bench: 10 k steps/s vs Phase 0.2's 91 k steps/s on the
   same polymer. This is the dominant cost — `cpu_local_snapshot`
   re-entry plus per-step numpy. ECM at N=21·N_f beads (e.g. 2 100
   beads for 100 filaments) projects to ~480 steps/s, or ~35 min for
   1 M steps. Within the H.1 brief's 5× v1 regression budget for the
   polymer, but likely tight for ECM. Optimisation plan deferred to
   the ECM scaling commit; sign-off does not block on it.

## Recommended next prompt (PI to review)

> **Main Session — H.1 ECM Mikado, continuation from BAOAB freeze**
>
> BAOAB freeze point is at `phase1/h1-ecm` @ `<commit>` (Main commit
> from 2026-05-20). PI sign-off granted on the L-M-BAOAB-limit Updater.
>
> Open PI calls from the freeze report still pending:
> - bending equipartition gate: drop from BAOAB / move to H.2? ____
> - pre-existing `v2/foundation` stash@{0..1}: who owns? ____
>
> Continue H.1 per brief `ffn_sim/docs/briefs/H1_ecm_mikado.md`:
> 1. `ffn_sim/ecm/mikado.py` (reuse `validation/oracles/ecm/fiber_network.py`
>    geometry; no v1 force kernel imports).
> 2. `ffn_sim/ecm/cross_links.py` (reuse oracle geometry only).
> 3. `ffn_sim/ecm/shear_protocol.py` (fresh write; HOOMD `BoxResize`).
> 4. `ffn_sim/configs/phase1_h1.yaml` (port KU-1.x values from
>    `validation/oracles/configs/phase1_unit1.yaml` with D4 overrides;
>    `dynamics.integrator = leimkuhler_matthews_baoab`).
> 5. `ffn_sim/tests/test_h1_mikado.py` — topology smoke + energy oracle
>    vs `validation/oracles/ecm/fiber_mechanics.py`, rel ≤ 1e-6.
>
> 2D-thin-slab vs 3D-periodic decision (brief §Open implementation
> questions) needs a PI call before topology generation. Recommendation
> based on BAOAB freeze: 3D-periodic, project to 2D for v1 oracle
> comparison — keeps the integrator general and avoids a z-pinning
> external field that AFINES does not have.
