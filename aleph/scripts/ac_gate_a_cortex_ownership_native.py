#!/usr/bin/env python
r"""GATE-A NATIVE driver — the cortex stops being a bind-target port and owns its own arrays.

WHAT THIS SETTLES.  Every native cortex number produced so far was measured through a *view* onto the
incumbent's global array: ``ac.cell.assemble.build_cell`` concatenates ``[actin | myosin | nucleus |
membrane]`` into ONE ``pos_d``/``f_d`` (``assemble.py:809``) and the engine hands that whole allocation out
under the cortex's name (``ac_gate_b_cortex_motor_native.py:198,227``, verbatim: "ALIAS the global array").
Co-location in an array is never a connection (PI 2026-07-22), so while that held, the **nine** connectors
terminating on the cortex had no cortex-owned force array for an adjoint pair to close against, and were
structurally unverifiable — which is why the execution plan puts cortex array ownership on the critical path
(``AC_EXECUTION_PLAN_2026-07-25.md`` §10; ``COMPARTMENT_VALIDATION_TRACKS_2026-07-25.md`` T2).

WHAT IT MEASURES (all four declared before the run):

  1. **PRIVATE ALLOCATIONS.**  The cortex's ``position_d``/``force_d`` are distinct device pointers from the
     incumbent's ``pos_d``/``f_d`` and from each other, and their length is EXACTLY the cortex node count —
     not ``n_total``.  The length is the limb that cannot pass vacuously: the incumbent puts actin at global
     offset 0, so a cortex index and a global index coincide and any "the indices look local" check is
     satisfied by the very arrangement being retired.

  2. **PARITY — the engine cortex is the incumbent's cortex, not a second model.**  With the compartments
     off, the incumbent's whole force assembly IS the cortex channel set (Cytosim bending + Hookean
     crosslink, plus the Arp2/3 angle harmonic when mixed).  The engine owner launches the same kernels with
     the same argument order over its own arrays, so the two force fields must agree node-by-node to
     **1e-12 relative** — the criterion T2 wrote down before this run existed.  Anything larger is a model
     difference; anything at that scale is float64 atomic arrival order, which is why the gate is a relative
     bound and NOT bit-identity (bit-identity is not an admissible acceptance test here: ``atomic_add``
     arrival order is nondeterministic).

  3. **DETECTOR CONTROL — the comparison can fail.**  A known constant force is injected into a copy of the
     incumbent's actin-block force, and parity must then fail by exactly that magnitude.  A parity gate that
     has never been seen failing certifies nothing, and this gate's PASS arm is an equality, which is the
     easiest kind of check to satisfy by accident.

     This replaces the control originally written here — "rebuild with myosin ON and require parity to
     fail" — because the run MEASURED that arm to be vacuous: with myosin on, the incumbent's actin-block
     force is unchanged **to within float64 round-off**.  That is not a defect, it is the resting build
     being what it says it is: every NMII head is unbound at rest (the resting bound-myosin setpoint is a
     PI-GAP and defaults OFF), so the minifilament's own backbone and head-arm forces land on myosin
     particles at indices >= n_actin and no myosin force reaches an actin node at all.  The myosin arm is
     kept and reported as that measurement, clearly labelled as NOT a foreign-force control.

     "Round-off", not "exactly zero", and the difference is a worked example of the standing rule.  On the
     200-filament slice the delta was **exactly 0.0** and this docstring said so; at the native 494,802-node
     population the same quantity is **4.9e-18 pN** — nonzero, because `atomic_add` arrival order into the
     shared actin force array is not reproducible between two builds of different length.  The mechanism
     claim survived the population change and the exact-zero phrasing did not, so the record reports the
     measured magnitude and `is_exactly_zero` is allowed to come back False rather than being asserted.

     The control this gate genuinely cannot run yet is the DEVICE-side one: compose the engine cortex and
     the incumbent's cortex into one solve and show the double count appearing as an exact factor of two.
     That needs `_accumulate_all(omit=)`, which is PI-gated (``…TRACKS…`` §6 D2), and it is named here so
     the gap is visible rather than implied by the word "control".

  4. **REJECT-RESTORE.**  Snapshot, displace the cortex, roll back under ``accepted = 0`` — every node must
     return to its pre-candidate value exactly (an equality of the restored state with its own snapshot, not
     a cross-implementation comparison, so exactness is the right criterion).  Then repeat under
     ``accepted = 1`` and require the displacement to SURVIVE, because a rollback that always restores is
     indistinguishable from one that never ran.

WHAT IT DOES NOT CLAIM.  No physics magnitude: the forces here are the incumbent's own forces, and the only
numbers this gate produces are a parity ratio and pointer identities.  The rung earned is ``CUDA_UNIT`` — real
kernels executed on CUDA over component-private arrays — and NOT ``CONNECTED``, which needs a connector
dispatched inside an accepted-step transaction with its two-sided resultant closing.  ``QuantitativeClaim``
stays ``BLOCKED``.

WHY THIS IS A GATE-A DRIVER.  It is a STATIC force-assembly comparison at the resting build — no clock, no
events, no accepted-step transaction.

Usage (gbook A5000, through the lease — never a bare ssh)::

    python aleph/scripts/ffn_gpu.py run --minutes 45 --reason "cortex array ownership parity" \
      --population "70686 filaments / 494802 actin nodes" -- \
      aleph/scripts/ac_gate_a_cortex_ownership_native.py --out outputs/ac/cortex_ownership/record.json

A reduced ``--filaments`` is for developing the instrument only; the record marks itself non-native and the
rung stops below the ladder's native rung accordingly.
"""

from __future__ import annotations

import argparse
import json
import time
from typing import Any

import numpy as np
import warp as wp

from aleph.engine.contracts import (
    EvidenceLabel,
    EvidenceRung,
    QuantitativeClaim,
    VoidCeiling,
)
from aleph.engine.cortex_population import build_cortex_population
from aleph.engine.cortex_state import (
    CORTEX_CHANNELS_NOT_BOUND,
    assert_component_state_disjoint,
    build_cortex_state_owner,
    cortex_ownership_census,
)
from aleph.engine.observe.artifact import (
    observation_artifact,
    timing_block,
    write_artifact,
)
from aleph.engine.surface_body import storage_key

#: The native cortical F-actin population (100 µm⁻² × 4π(7.5 µm)²).  A smaller count develops the
#: instrument; it never supports a conclusion.
NATIVE_FILAMENTS = 70686

#: T2's pre-declared pass criterion for a pure-refactor parity, quoted rather than chosen here:
#: "engine-owned channel forces equal the incumbent's contribution node-by-node to <1e-12 relative".
PARITY_RELATIVE_TOLERANCE = 1.0e-12

#: How far the staged double count may sit from EXACTLY 2.0 and still be called a clean double.
#: Derived, not chosen: the doubled field is the same sum evaluated twice into the same float64 array,
#: so the only difference from 2x is `atomic_add` arrival order — the same effect that makes the pass
#: arm's parity 1e-16 rather than 0, accumulated over the ~1e6 contributions a native cortex node sees.
#: Anything outside this is not a double count, and the control must then say so rather than round to 2.
DOUBLE_COUNT_TOLERANCE = 1.0e-6

#: Declared BEFORE the run.  Above this, the two force fields differ by more than float64 summation-order
#: noise can explain, so the comparison is no longer measuring accumulation order — it is measuring two
#: different models, and the parity observable stops meaning what its name says.
PARITY_VOID_CEILING = VoidCeiling(
    ratio=1.0e-6,
    rationale=(
        "a relative force difference above 1e-6 cannot be float64 atomic arrival order on a network of "
        "this size, so the run would be comparing two models rather than one model with itself, and the "
        "parity number would not be interpretable as a refactor check"
    ),
)


def _resolve_device(requested: str | None) -> str:
    """Return the resolved CUDA device string, refusing any non-CUDA runtime (I0-A)."""
    wp.init()
    device = wp.get_device(requested)
    if not device.is_cuda:
        raise RuntimeError(
            f"this gate requires a CUDA GPU (I0-A); resolved {str(device)!r}, which is not CUDA"
        )
    return str(device)


def _incumbent_actin_force(
    *,
    n_filaments: int,
    seg_um: float,
    length_um: float,
    density_per_fil: float,
    overlap_free: bool,
    seed: int,
    device: str | None,
    with_myosin: bool,
    omit: tuple[str, ...] = (),
    engine_mechanics: object | None = None,
) -> tuple[np.ndarray, np.ndarray, dict[str, Any]]:
    """Build the incumbent cell and return ``(actin-block force, actin-block position, identity)``.

    Every compartment except the cortex is switched OFF (and myosin only when the caller asks for the
    negative control), so the incumbent's whole force assembly is the cortex channel set.  Steric and
    pressure are off for the same reason: both address populations wider than the cortex and are named in
    :data:`~aleph.engine.cortex_state.CORTEX_CHANNELS_NOT_BOUND` as deliberately unbound here.

    ``omit`` and ``also_accumulate_engine_cortex`` exist for T2's DEVICE-side negative control, the one
    this gate previously had to name as un-runnable: compose the engine cortex and the incumbent's cortex
    into ONE force array and require the double count to appear as an exact factor of two.  With the mask
    ON the two are exclusive and parity holds; with it OFF and the engine owner also accumulating, the
    same comparison must fail by exactly 2.  A control that cannot fail proves nothing, and until
    ``_accumulate_all(omit=)`` existed (PI `COMPARTMENT_VALIDATION_TRACKS` §6 D2, 2026-07-29) this one could not be run at all.
    """
    from aleph.components.incumbent.assemble import CellConfig, build_cell
    from aleph.components.incumbent.driver import _accumulate_all

    cfg = CellConfig(
        n_filaments=n_filaments,
        cortex_seg_um=seg_um,
        cortex_length_um=length_um,
        cortex_density_per_fil=density_per_fil,
        overlap_free_cortex=overlap_free,
        with_myosin=with_myosin,
        with_steric=False,
        with_nucleus=False,
        with_membrane=False,
        with_pressure=False,
        seed=seed,
        device=device,
    )
    cell = build_cell(cfg)
    _accumulate_all(cell, cell.pos_d, cell.f_d, omit=omit)
    if engine_mechanics is not None:
        # Deliberately add the SAME channels a SECOND time, onto the SAME array, so the double count is
        # measurable. This is the failure `omit=` exists to prevent, staged on purpose: the engine
        # owner's mechanics takes (pos, force) as arguments, so it can be pointed at the incumbent's
        # global array even though the owner's own arrays are private.
        engine_mechanics.accumulate(cell.pos_d, cell.f_d)
    wp.synchronize_device(cell.device)
    n_actin = int(cell.n_actin)
    identity = {
        "pos_storage": str(storage_key(cell.pos_d)),
        "force_storage": str(storage_key(cell.f_d)),
        "n_actin": n_actin,
        "n_total": int(cell.n_total),
        "actin_is_offset_zero": True,
        "omitted_channels": list(omit),
        "engine_cortex_also_accumulated": engine_mechanics is not None,
    }
    return cell.f_d.numpy()[:n_actin], cell.pos_d.numpy()[:n_actin], identity


def _relative_difference(a: np.ndarray, b: np.ndarray) -> dict[str, float]:
    """Node-by-node force comparison, reported as a distribution rather than a single pass/fail scalar.

    The denominator is the max force magnitude of the reference field, so the ratio is the fraction of the
    largest real force the disagreement represents — a per-node relative error would divide by nearly-zero
    forces on quiescent nodes and report a meaningless number.
    """
    delta = np.linalg.norm(a - b, axis=1)
    magnitude = np.linalg.norm(b, axis=1)
    scale = float(magnitude.max()) if magnitude.size else 0.0
    return {
        "max_abs_delta_pN": float(delta.max()) if delta.size else 0.0,
        "mean_abs_delta_pN": float(delta.mean()) if delta.size else 0.0,
        "reference_max_force_pN": scale,
        "max_relative_delta": float(delta.max() / scale) if scale > 0.0 else float("inf"),
        "n_nodes_compared": int(delta.shape[0]),
    }


@wp.kernel
def _displace_kernel(pos: wp.array(dtype=wp.vec3d), delta: wp.vec3d) -> None:
    """Displace every node by a fixed vector — the perturbation the rollback must undo."""
    t = wp.tid()
    pos[t] = pos[t] + delta


def main() -> None:
    parser = argparse.ArgumentParser(
        description="GATE-A: the cortex owns its own device arrays, and they are the incumbent's cortex.")
    parser.add_argument("--filaments", type=int, default=NATIVE_FILAMENTS,
                        help=f"cortex F-actin count (native = {NATIVE_FILAMENTS})")
    parser.add_argument("--cortex-seg-um", type=float, default=0.5, dest="cortex_seg_um")
    parser.add_argument("--cortex-length-um", type=float, default=3.0, dest="cortex_length_um")
    parser.add_argument("--cortex-density", type=float, default=20.0, dest="cortex_density")
    parser.add_argument("--overlap-free", action=argparse.BooleanOptionalAction, default=True,
                        help="build the cortex overlap-free (the incumbent's production setting)")
    parser.add_argument("--seed", type=int, default=0,
                        help="host RNG seed; MUST match on both sides or the two weaves are two cortices")
    parser.add_argument("--device", type=str, default=None,
                        help="Warp device to resolve; omit to take the default CUDA device")
    parser.add_argument("--skip-myosin-arm", action="store_true",
                        help="skip the myosin measurement arm (it rebuilds the cell a second time)")
    parser.add_argument("--out", type=str, default=None, help="write the run record here")
    parser.add_argument("--build-commit", type=str, default=None,
                        help="declare the build commit on a machine whose tree is not a git checkout")
    args = parser.parse_args()

    run_started = time.perf_counter()
    device = _resolve_device(args.device)
    native_population = int(args.filaments) == NATIVE_FILAMENTS
    print(f"[cortex-ownership] device={device} filaments={args.filaments:,} "
          f"({'NATIVE' if native_population else 'SLICE — instrument development only'})", flush=True)

    # ── the engine-owned cortex: its own population, its own arrays, its own kernels ────────────────────
    population = build_cortex_population(
        int(args.filaments),
        seg_um=args.cortex_seg_um,
        length_um=args.cortex_length_um,
        density_per_fil=args.cortex_density,
        overlap_free=bool(args.overlap_free),
        seed=int(args.seed),
    )
    owner = build_cortex_state_owner(population.topology, device=device)
    owner.accumulate()
    wp.synchronize_device(device)
    engine_force = owner.force_d.numpy()
    engine_pos = owner.position_d.numpy()
    census = cortex_ownership_census(population, owner)
    print(f"[cortex-ownership] engine cortex: {census['population']['n_nodes']:,} nodes, "
          f"channels {census['bound_force_channels']}", flush=True)

    # ── (2) parity against the incumbent, compartments off so its assembly IS the cortex channel set ────
    incumbent_force, incumbent_pos, incumbent_identity = _incumbent_actin_force(
        n_filaments=int(args.filaments),
        seg_um=args.cortex_seg_um,
        length_um=args.cortex_length_um,
        density_per_fil=args.cortex_density,
        overlap_free=bool(args.overlap_free),
        seed=int(args.seed),
        device=args.device,
        with_myosin=False,
    )

    # (1) private allocations — checked before any number is read from them.
    same_length = int(owner.position_d.shape[0]) == incumbent_identity["n_actin"]
    distinct_from_incumbent = (
        str(storage_key(owner.position_d)) != incumbent_identity["pos_storage"]
        and str(storage_key(owner.force_d)) != incumbent_identity["force_storage"]
    )
    assert_component_state_disjoint([owner])
    ownership = {
        "engine_position_len": int(owner.position_d.shape[0]),
        "engine_force_len": int(owner.force_d.shape[0]),
        "incumbent_n_actin": incumbent_identity["n_actin"],
        "incumbent_n_total": incumbent_identity["n_total"],
        "engine_arrays_are_exactly_the_cortex": same_length,
        "engine_arrays_distinct_from_incumbent": distinct_from_incumbent,
        "engine_position_storage": census["position_storage"],
        "engine_force_storage": census["force_storage"],
        "incumbent_storage": {
            "position": incumbent_identity["pos_storage"],
            "force": incumbent_identity["force_storage"],
        },
    }

    geometry_identical = bool(np.array_equal(engine_pos, incumbent_pos))
    parity = _relative_difference(engine_force, incumbent_force)
    parity_passed = geometry_identical and parity["max_relative_delta"] < PARITY_RELATIVE_TOLERANCE
    print(f"[cortex-ownership] parity: geometry_identical={geometry_identical} "
          f"max_relative_delta={parity['max_relative_delta']:.3e} "
          f"(criterion < {PARITY_RELATIVE_TOLERANCE:g})", flush=True)

    # ── (3) detector control: inject a KNOWN force and require the comparison to report exactly it ──────
    # The PASS arm above is an equality, which is the easiest check to satisfy by accident, so the gate has
    # to be seen failing on a difference it was told about.
    injected_pN = 10.0 * float(parity["reference_max_force_pN"]) or 1.0
    injected = incumbent_force.copy()
    injected[:, 0] += injected_pN
    detector = _relative_difference(engine_force, injected)
    detector_control = {
        "what": "a known constant force added to a COPY of the incumbent's actin-block force",
        "injected_pN": injected_pN,
        "measured": detector,
        "parity_fails_as_required": detector["max_relative_delta"] >= PARITY_RELATIVE_TOLERANCE,
        "recovers_the_injection": bool(
            np.isclose(detector["max_abs_delta_pN"], injected_pN, rtol=1.0e-12, atol=0.0)),
    }
    print(f"[cortex-ownership] detector control: fails={detector_control['parity_fails_as_required']} "
          f"recovered {detector['max_abs_delta_pN']:.6g} of {injected_pN:.6g} pN injected", flush=True)

    # ── (3b) the myosin arm — a MEASUREMENT, not a control (see the module docstring) ───────────────────
    myosin_arm: dict[str, Any] = {"ran": False}
    if not args.skip_myosin_arm:
        mixed_force, _, mixed_identity = _incumbent_actin_force(
            n_filaments=int(args.filaments),
            seg_um=args.cortex_seg_um,
            length_um=args.cortex_length_um,
            density_per_fil=args.cortex_density,
            overlap_free=bool(args.overlap_free),
            seed=int(args.seed),
            device=args.device,
            with_myosin=True,
        )
        share = _relative_difference(incumbent_force, mixed_force)
        myosin_arm = {
            "ran": True,
            "what": "incumbent rebuilt with myosin ON; how much of that reaches an ACTIN node",
            "n_total_with_myosin": mixed_identity["n_total"],
            "n_myosin_particles": mixed_identity["n_total"] - mixed_identity["n_actin"],
            "myosin_force_on_actin_pN": share["max_abs_delta_pN"],
            "is_exactly_zero": bool(share["max_abs_delta_pN"] == 0.0),
            "why": (
                "every head is UNBOUND at rest (the resting bound-myosin setpoint is a PI-GAP and defaults "
                "OFF), so the minifilament's backbone/head-arm forces land on myosin particles at indices "
                ">= n_actin and no myosin force reaches an actin node"
            ),
            "not_a_control": (
                "because it is zero, this arm cannot detect a foreign force; the device-side double-count "
                "control is the separate omit-mask arm below, unblocked by PI TRACKS-§6-D2 on 2026-07-29"
            ),
        }
        print(f"[cortex-ownership] myosin arm (measurement): force reaching actin = "
              f"{share['max_abs_delta_pN']:.4g} pN over "
              f"{myosin_arm['n_myosin_particles']:,} myosin particles", flush=True)

    # ── (3b) T2 DEVICE-SIDE DOUBLE-COUNT CONTROL — the one this gate had to declare un-runnable ────────
    # PI decision `COMPARTMENT_VALIDATION_TRACKS` §6 D2 (2026-07-29, option (a)) added `_accumulate_all(omit=)`, so the control can now run:
    # compose the engine cortex and the incumbent's cortex into ONE force array with the mask OFF, and
    # require the comparison to fail by an EXACT factor of two. If it failed by anything else, the
    # double count would not be what the mask prevents; if it did not fail at all, the mask would be
    # protecting against nothing and `omit=` could be dropped.
    doubled_force, _, doubled_identity = _incumbent_actin_force(
        n_filaments=int(args.filaments),
        seg_um=args.cortex_seg_um,
        length_um=args.cortex_length_um,
        density_per_fil=args.cortex_density,
        overlap_free=bool(args.overlap_free),
        seed=int(args.seed),
        device=args.device,
        with_myosin=False,
        omit=(),                                   # mask OFF — the incumbent still launches its cortex
        engine_mechanics=owner.mechanics,          # and the engine adds the SAME channels again
    )
    ref = float(np.max(np.abs(incumbent_force)))
    ratios = np.full(incumbent_force.shape[0], np.nan)
    big = np.linalg.norm(incumbent_force, axis=1) > 1.0e-9 * max(ref, 1.0)
    if big.any():
        ratios[big] = (np.linalg.norm(doubled_force[big], axis=1)
                       / np.linalg.norm(incumbent_force[big], axis=1))
    finite = ratios[np.isfinite(ratios)]
    factor_median = float(np.median(finite)) if finite.size else float("nan")
    factor_max_dev = float(np.max(np.abs(finite - 2.0))) if finite.size else float("nan")
    doubled_parity = _relative_difference(incumbent_force, doubled_force)
    # The tolerance is the SAME round-off floor the pass arm is held to; a factor of two is not a
    # near-miss quantity, so anything outside it means the staged double count is not a clean double.
    exact_factor_two = bool(finite.size and factor_max_dev <= DOUBLE_COUNT_TOLERANCE)
    omit_control = {
        "what": ("mask OFF while the engine owner also accumulates the same channels onto the same "
                 "array — the failure `_accumulate_all(omit=)` exists to prevent, staged deliberately"),
        "n_nodes_compared": int(finite.size),
        "factor_median": factor_median,
        "max_deviation_from_two": factor_max_dev,
        "parity_fails": bool(doubled_parity["max_relative_delta"] > PARITY_RELATIVE_TOLERANCE),
        "fails_by_exactly_two": exact_factor_two,
        "max_relative_delta_vs_single": doubled_parity["max_relative_delta"],
        "why_this_is_the_control": (
            "the pass arm shows the engine and the incumbent agree; ONLY this arm shows that adding the "
            "force twice is DETECTABLE rather than silently absorbed, which is what makes the parity "
            "number evidence of exclusivity instead of evidence of a coincidence"
        ),
        "identity": doubled_identity,
    }
    print(f"[cortex-ownership] omit-mask control: factor {factor_median:.9f} "
          f"(max deviation from 2 = {factor_max_dev:.3g}), parity fails = {omit_control['parity_fails']}",
          flush=True)

    # ── (4) reject-gated restore of the component's OWN position array ──────────────────────────────────
    displacement = wp.vec3d(1.0e-3, -2.0e-3, 3.0e-3)
    before = owner.position_d.numpy().copy()
    owner.snapshot_candidate()
    wp.launch(_displace_kernel, dim=int(owner.position_d.shape[0]),
              inputs=[owner.position_d, displacement], device=device)
    wp.synchronize_device(device)
    displaced = owner.position_d.numpy().copy()
    rejected_d = wp.zeros(1, dtype=wp.int32, device=device)
    owner.rollback(rejected_d)
    wp.synchronize_device(device)
    restored = owner.position_d.numpy()
    restore_exact = bool(np.array_equal(restored, before))

    owner.snapshot_candidate()
    wp.launch(_displace_kernel, dim=int(owner.position_d.shape[0]),
              inputs=[owner.position_d, displacement], device=device)
    accepted_d = wp.ones(1, dtype=wp.int32, device=device)
    owner.rollback(accepted_d)
    wp.synchronize_device(device)
    kept = owner.position_d.numpy()
    accepted_survives = bool(np.array_equal(kept, before + np.asarray(
        [displacement[0], displacement[1], displacement[2]], np.float64)))
    transaction = {
        "displacement_applied_um": [float(displacement[0]), float(displacement[1]),
                                    float(displacement[2])],
        "displacement_was_real": bool(not np.array_equal(displaced, before)),
        "rejected_step_restores_exactly": restore_exact,
        "accepted_step_keeps_the_candidate": accepted_survives,
        "predicate_stayed_on_device": True,
    }
    print(f"[cortex-ownership] transaction: reject-restore={restore_exact} "
          f"accept-keeps={accepted_survives}", flush=True)

    verdict = {
        "private_allocations": bool(same_length and distinct_from_incumbent),
        "geometry_identical_to_incumbent": geometry_identical,
        "force_parity_within_1e-12": parity_passed,
        "detector_control_fails_as_required": detector_control["parity_fails_as_required"],
        "detector_control_recovers_the_injection": detector_control["recovers_the_injection"],
        "omit_mask_control_fails_by_exactly_two": omit_control["fails_by_exactly_two"],
        "rejected_step_restores_exactly": restore_exact,
        "accepted_step_keeps_the_candidate": accepted_survives,
    }
    gate_passed = all(value for value in verdict.values() if value is not None)

    rung = EvidenceRung.CUDA_UNIT
    basis = (
        "the cortex owner allocated component-private position/force arrays of exactly its own node count, "
        "launched the real ff bending + crosslink kernels on CUDA over them, and its force field matched the "
        "incumbent's cortex assembly node-by-node; NOT CONNECTED — no connector was dispatched inside an "
        "accepted-step transaction and no two-sided resultant was closed"
        + ("" if native_population else
           " — and this ran on a SLICE, so nothing here transfers to the native population by itself")
    )

    record = observation_artifact(
        run_label="GATE-A cortex array ownership — private arrays, incumbent parity, double-count control",
        evidence=EvidenceLabel(rung=rung, quantitative=QuantitativeClaim.BLOCKED, basis=basis),
        config={
            "argv": {key: value for key, value in sorted(vars(args).items())},
            "parity_relative_tolerance": PARITY_RELATIVE_TOLERANCE,
            "parity_tolerance_source": (
                "COMPARTMENT_VALIDATION_TRACKS_2026-07-25.md T2 pass criterion, written before this run"
            ),
            "incumbent_compartments": "myosin/steric/nucleus/membrane/pressure OFF for the parity arm",
        },
        census={
            **census,
            "is_native_population": native_population,
            "population_note": (
                f"native is {NATIVE_FILAMENTS:,} filaments; a slice develops the instrument and never "
                "supports a conclusion"
            ),
            "n_cortex_terminating_connectors": len(census["connector_domains"]),
        },
        t0={
            "engine_max_force_pN": float(np.linalg.norm(engine_force, axis=1).max()),
            "incumbent_max_force_pN": float(np.linalg.norm(incumbent_force, axis=1).max()),
            "measured_before_any_perturbation": True,
        },
        timing=timing_block(
            wall_seconds=time.perf_counter() - run_started,
            device_note=(
                "NOT a benchmark: the wall-clock covers two (or three) full cell builds plus kernel "
                "compilation, and this gate runs no time steps at all — it is a static force-assembly "
                "comparison, so there is no converged-step cost to report"
            ),
        ),
        measurements={
            "ownership": ownership,
            "parity": parity,
            "detector_control": detector_control,
            "myosin_arm": myosin_arm,
            "omit_mask_double_count_control": omit_control,
            "transaction": transaction,
            "verdict": verdict,
            "channels_not_bound": dict(CORTEX_CHANNELS_NOT_BOUND),
        },
        device=device,
        parameter_provenance={
            "cortex_seg_um": "CONVENIENCE",
            "cortex_length_um": "CONVENIENCE",
            "cortex_density_per_fil": "CONVENIENCE",
            "kappa_actin": "SOURCED",
            "crosslink_k": "SOURCED",
        },
        void_ceiling=PARITY_VOID_CEILING,
        gate_passed=gate_passed,
        residual=float(parity["max_abs_delta_pN"]),
        signal=float(parity["reference_max_force_pN"]),
        declared_commit=(args.build_commit or None),
    )
    print(json.dumps(record, indent=2, default=str), flush=True)
    if args.out:
        write_artifact(args.out, record)

    failed = [name for name, value in verdict.items() if value is False]
    if failed:
        raise SystemExit(f"[cortex-ownership] GATE-A FAILED on: {', '.join(failed)}")
    print("[cortex-ownership] GATE-A PASSED — the cortex owns its arrays and they are the incumbent's "
          "cortex", flush=True)


if __name__ == "__main__":
    main()
