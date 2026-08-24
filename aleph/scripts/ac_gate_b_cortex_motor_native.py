#!/usr/bin/env python
r"""GATE-B NATIVE driver — emergent cortical NMII tension on the FULL 70,686-filament cortex.

Demonstrates the GATE-B MECHANISM on the real composed cell (CUDA / gbook A5000):

    myosin heads bind by k_on EVENTS  →  the split crossbridge scatters force onto the cortex actin
    →  cortical tension γ EMERGES over accepted physical steps, starting from ZERO bound heads.

There is NO static prestress and NO resting-bound-myosin seed: at t0 every straddle-placed NMII head is
UNBOUND, and the tension is grown purely by the accepted-step event runtime
(:class:`aleph.engine.cortex_motor_slice.CortexMotorSlice` driving
:class:`aleph.engine.transaction.CellTransaction`).  Each accepted step: propose attach events (live
point-to-segment query) → inner mechanical relax (the cortex responds to the crossbridge pull) → per-head
Hill/Bell loads → commit the k_on-Poisson attach + catch-slip/Bell detach KMC → advance the device clock.

────────────────────────────────────────────────────────────────────────────────────────────────────────
RUN ON GBOOK (this needs a CUDA GPU; it will NOT run on the dev Mac):

    ssh gbook
    cd ~/ffn_ac_native
    PYTHONPATH=~/ffn_ac_native:~/ffn_ac_native/ffn_sim \
      ~/miniconda3/envs/ffn_sim/bin/python \
      aleph/scripts/ac_gate_b_cortex_motor_native.py --steps 40 --dt 0.01 --outer 40

    # physiological catch-slip detach (Kovacs 2007) instead of the Bell-slip smoke path:
    PYTHONPATH=~/ffn_ac_native:~/ffn_ac_native/ffn_sim \
      ~/miniconda3/envs/ffn_sim/bin/python \
      aleph/scripts/ac_gate_b_cortex_motor_native.py --steps 40 --dt 0.01 --outer 40 --catch-slip

    # elastic-artifact signature — sweep k_xb; the r0_bind-FIXED residual must NOT scale ∝ k_xb (old one did):
    for K in 250 1000 4000; do PYTHONPATH=~/ffn_ac_native:~/ffn_ac_native/ffn_sim \
      ~/miniconda3/envs/ffn_sim/bin/python \
      aleph/scripts/ac_gate_b_cortex_motor_native.py --steps 20 --dt 0.01 --outer 40 --k-xb $K ; done

Long run: launch under nohup and monitor the LOG FILE (ssh python is not on PATH; use the full env python).
────────────────────────────────────────────────────────────────────────────────────────────────────────

WHAT TO EXPECT (with the PROVISIONAL params below — a MECHANISM demonstration, NOT quantitative production):
  * `bound` climbs from 0 toward ~50-65% of the heads over the first ~5-15 steps (k_on≈50/s, dt=0.01 →
    P_attach≈1-exp(-0.5)≈0.39/step for a head within capture; straddle places ~60% of heads within the
    0.05 µm capture reach, so the reachable subset saturates and the rest stay free — an EMERGENT duty, not
    an imposed fraction);
  * γ_source (the crossbridge dipoles) rises from 0 with the bound count; γ_network (actin-constraint +
    crosslink) carries it onward; γ_total is the observable. With f_stall≈0.5 pN/head (PROVISIONAL, GAP) and
    a partially-bound cortex the magnitude is LOW (order 1e-1–1e1 pN/µm) — this demo shows the SIGN and the
    RISE (γ grows monotonically with the bound population), NOT the physiological ~100s pN/µm magnitude,
    which is a sourced-parameter (f_stall / N_side / duty) PI decision, not a runtime outcome to tune;
  * max|PF| (the projected residual) stays finite and the cortex tracks the crossbridge pull each step.

FORCE-ACCEPT: the demo force-accepts every step (`accepted_d = ones`) so binding accumulates and the
mechanism is visible; a resting-baseline PRODUCTION run would instead gate on the force-balance predicate
(reject a non-converged candidate) once γ is high enough to hold turgor.  Pass --balance-gate to switch.

CONFIG CHOICES (said out loud, per the brief):
  * `resting_bound_myosin_fraction=None` (default) is the exact knob that BUILDS the straddle-placed
    minifilament STRUCTURE with ALL heads UNBOUND — the minifilaments are placed together with the cortex
    regardless of the resting seed; the seed only *pre-binds* a fraction.  None ⇒ purely dynamic, which is
    what GATE-B needs.  (Confirmed in `ac/cell/assemble.py::build_cell` + `_build_myosin`: heads use
    `allocate_hand_state` default bound=0; the resting seed is a separate, source-gated, default-OFF path.)
  * `nmii_backbone_lp_um` is set to the DIAGNOSTIC convergence-fixture value (a PI/KB GAP) with an explicit
    DIAGNOSTIC source, PURELY so `_build_myosin` materialises the F6 backbone-angle topology the slice's
    `BackboneArmMechanics` binds.  It does NOT affect the unbound-at-rest condition, and the backbone is
    straight at its rest so the bending force is ≈0 at t0.  (Left as None the resting build omits the
    backbone-angle device array, which `BackboneArmMechanics` requires.)
  * INTEGRATION: the assembled cell composes [actin | myosin | …] into ONE global `pos_d`/`f_d`.  The slice
    was written for split ownership; here the actuator and the cortex port both ALIAS that one global array
    with GLOBAL head/segment indices, so the "split" is logical (nmii owns the myosin index block, cortex
    owns the actin block) — exactly how `assemble.py` already composes the motor and connects it via the
    crossbridge kernel.  Newton's 3rd law still holds within the global force array (head +f, two cortex
    nodes −(1−t)f/−t·f, at distinct indices).  `cell.myosin` is replaced by a thin adapter so the driver's
    own inner relax (`make_inner_solve` → `_accumulate_all`) drives the slice's internal mechanics + split
    crossbridge each iteration (one binding SoA, no double count).

Runtime: NVIDIA Warp on CUDA only (I0-A).  Authored on the dev Mac (no CUDA) — the Lead runs it on gbook.
"""

from __future__ import annotations

import argparse
import dataclasses
import math
import time
from pathlib import Path

import numpy as np
import warp as wp

from aleph.components.incumbent.assemble import (
    NMII_BACKBONE_LP_DIAGNOSTIC_UM,
    NMII_CAPTURE_UM,
    NMII_F0,
    NMII_F_STALL_TEST,
    NMII_HEAD_OFFSET_UM,
    NMII_KAPPA_TEST,
    NMII_KOFF0,
    NMII_KON_TEST,
    NMII_K_XB_TEST,
    NMII_N_BB,
    NMII_N_SIDE,
    NMII_V0_TEST,
    R_CORTEX_UM,
    CellConfig,
    build_cell,
)
from aleph.components.incumbent.assembled_cortical_stress import actin_axial_tension, crosslink_tension
from aleph.components.incumbent.cortical_tension import measure_cortical_stress
from aleph.components.incumbent.driver import _accumulate_all, make_inner_solve
from aleph.engine.contracts import EvidenceLabel, EvidenceRung, QuantitativeClaim
from aleph.engine.observe.artifact import (
    observation_artifact,
    timing_block,
    write_artifact,
)
from aleph.engine.observe.stationarity import (
    StationarityVerdict,
    assess_observables,
)
from aleph.engine.cortex_motor_slice import (
    CortexMotorParams,
    build_cortex_motor_slice,
    build_nmii_actuator_state,
)
from aleph.engine.nmii_actuator import NMIIBendingStiffness
from aleph.components.motor.bell_kinetics_analytic import bound_state_lifetime
from aleph.components.motor.segment_motor import SegmentDetachKinetics

#: Observables judged for stationarity, in the order the summary prints them.  ``gamma_total`` is the
#: one a run is steered by; the rest are recorded so a reader can see WHICH part is still moving —
#: a settled network carrying a drifting source is a different finding from neither having settled.
STATIONARITY_OBSERVABLES: tuple[str, ...] = (
    "gamma_total_pn_per_um", "gamma_network_pn_per_um", "gamma_source_pn_per_um", "bound_fraction",
)

# ── PROVISIONAL catch-slip constants (Pereverzev two-pathway) — NOT SOURCED, PI/KB GAP ────────────────
# There is NO audited NMII head-actin catch-slip parameterisation in the Contract-Graph.  These proxy values
# only give the physiological SHAPE (off-rate falls with load up to F*, then rises — Kovacs 2007) so the
# --catch-slip mechanism path runs; they are constrained in spirit by the 5×/12× ADP-slowdown + NM2B duty
# 0.2-0.3 but are a MECHANISM-DEMO PROXY, never a sourced magnitude.  PI must source them before any
# quantitative claim.  (k in 1/s, x in µm.)
CATCH_SLIP_PROVISIONAL = dict(k_catch0=0.35, x_catch=1.0e-3, k_slip0=0.35, x_slip=0.6e-3)


def _peak_gpu_bytes(device: object) -> int | str:
    """Return the memory-pool high-water mark in bytes, or a STRING saying why there is no number.

    The hard runtime contract requires "exact peak GPU bytes" on every native run, and no record in this
    repo carried one — which is why a ``--membrane-subdiv`` reduction taken as "an A5000 memory safety
    valve" could never be checked against the memory it actually saved.

    **It returns a reason, never a bare ``None``.**  The first version took a ``Device``, was handed the
    driver's ``dev`` — which is a ``str`` (line ~591) — raised ``AttributeError`` inside a blanket
    ``except``, and wrote ``None``.  The record then read "not available" when the truth was "this
    telemetry crashed", and the two are not the same fact.  That is the exact defect class this session
    spent the day auditing, reproduced inside its own remedy, so the signature now resolves device-likes
    and every failure path names itself in the artifact.

    ``wp.get_device`` accepts a ``Device`` as readily as a string, so there is deliberately NO isinstance
    guard here: the first attempt had one keyed on ``wp.context.Device``, and ``wp.context`` does not
    exist in Warp 1.14 (it moved to ``warp._src.context``).  The guard raised ``AttributeError`` on every
    call and the reason string then blamed the device for not resolving.  Narrowing an accepted input to
    check its type is how that happened; the function now just asks Warp.

    Args:
        device: Anything Warp accepts as a device — ``"cuda:0"``, a ``Device``, or ``None``.

    Returns:
        Peak pool bytes as an ``int``, or a short string explaining why no number exists.
    """
    try:
        resolved = wp.get_device(device)
    except Exception as exc:                             # noqa: BLE001 — the reason is the payload
        return f"unavailable: device {device!r} did not resolve ({type(exc).__name__}: {exc})"
    if not resolved.is_cuda:
        return f"unavailable: {resolved} is not a CUDA device"
    if not resolved.is_mempool_supported:
        return f"unavailable: {resolved} has no CUDA memory pool"
    try:
        return int(wp.get_mempool_used_mem_high(resolved))
    except Exception as exc:                             # a telemetry read must not fail a run...
        return f"unavailable: mempool high-water read failed ({type(exc).__name__}: {exc})"  # ...but it must SAY so


def _banner(params: CortexMotorParams, mode: str) -> None:
    """Print the loud PROVISIONAL-PARAMS banner listing every value + its PI-GAP status."""
    gap = "GAP — PI (params_i0b3.yaml value: null / still_gap_under_nm2b)"
    prov = "PROVISIONAL (draft/Medium; params_i0b3.yaml)"
    phys = "physical constant"
    rows = [
        ("k_on   [1/s]", params.k_on, gap + "  (provisional 50/s = ff/hand_kmc NMIIA preset, NOT audited)"),
        ("f_stall[pN] ", params.f_stall, gap + "  (single-molecule NM2 stall NOT FOUND; 0.5 = AFINES claim_a)"),
        ("v0     [µm/s]", params.v0, gap),
        ("kappa  [-]  ", params.kappa, "Kovács 2003 Hill curvature a/F0=0.5 (dimensionless; NOT a force)"),
        ("k_xb   [pN/µm]", params.k_xb, gap + "  (MASTER knob, mid physical band 100-1000)"),
        ("r0_head[µm] ", params.r0_head, gap + "  (= head↔backbone arm offset)"),
        ("r0_xb  [µm] ", params.r0_xb, "≈0 — a bound head sits on the actin site"),
        ("k_off0 [1/s]", params.k_off0, prov + "  (Bell slip prefactor; Stam-Hocky/Tam)"),
        ("f0     [pN] ", params.f0, phys + "  (= kBT/x_beta ≈ 7.13 pN; Veigel 2002 x_beta)"),
        ("capture[µm] ", params.capture_radius, prov),
    ]
    if params.detach_kinetics is SegmentDetachKinetics.CATCH_SLIP:
        for name in ("k_catch0", "x_catch", "k_slip0", "x_slip"):
            rows.append((f"{name:11s}", getattr(params, name),
                         "PROVISIONAL PROXY — NOT SOURCED, PI/KB GAP (mechanism-shape only)"))
    print("=" * 100)
    print(f"  GATE-B PROVISIONAL PARAMETERS  —  detach mode: {mode}")
    print("  MECHANISM DEMONSTRATION, NOT quantitative production — every magnitude below is a PI decision.")
    print("-" * 100)
    for label, value, status in rows:
        print(f"    {label} = {float(value):>10.5g}   {status}")
    print("=" * 100, flush=True)


def _params(catch_slip: bool, k_xb: float) -> CortexMotorParams:
    """Assemble the slice params from the params_i0b3 provisional magnitudes (+ the catch-slip proxy).

    ``k_xb`` (the crossbridge MASTER knob) is the CLI-overridable scaling used to confirm the elastic-artifact
    signature: BEFORE the r0_bind fix the residual scaled ∝ k_xb (a purely-elastic placement force); with the
    attach-unstrained fix the passive force is 0 at bind, so the residual should NOT scale with k_xb.
    """
    kw = dict(
        v0=NMII_V0_TEST, f_stall=NMII_F_STALL_TEST, kappa=NMII_KAPPA_TEST, k_xb=float(k_xb),
        r0_head=NMII_HEAD_OFFSET_UM, r0_xb=0.0, capture_radius=NMII_CAPTURE_UM,
        k_on=NMII_KON_TEST, k_off0=NMII_KOFF0, f0=NMII_F0,
    )
    if catch_slip:
        return CortexMotorParams(detach_kinetics=SegmentDetachKinetics.CATCH_SLIP,
                                 **CATCH_SLIP_PROVISIONAL, **kw)
    return CortexMotorParams(detach_kinetics=SegmentDetachKinetics.SLIP, **kw)


def _build_actuator(cell, device, k_xb: float):
    """Extract the straddle-placed minifilament actuator from the built cell → an NMII state-owner.

    The slice's ``BackboneArmMechanics`` binds the SAME device topology (backbone/head bonds + F6 angle
    triples) and constants the assembled ``cell.myosin`` was built with, over the global ``pos_d`` (GLOBAL
    particle indices).  The per-minifilament/per-head metadata is synthesised (not consumed by the mechanics/
    crossbridge/KMC kernels), so no biological count is invented.  ``k_xb`` feeds the F6 head-arm bending
    stiffness for consistency with the CLI crossbridge knob (the arm is ≈force-free at the straight rest).
    """
    myo = cell.myosin
    if myo is None or myo.segment_runtime is None:
        raise RuntimeError("cell.myosin (with segment runtime) is required; build with with_myosin=True")
    if myo.backbone_angles is None:
        raise RuntimeError("F6 backbone-angle topology absent — build with nmii_backbone_lp_um set (see header)")
    n_heads = int(myo.head_node.shape[0])
    n_part = NMII_N_BB + 2 * NMII_N_SIDE
    n_min = int(cell.ledger.get("myosin_n_minifilaments", n_heads // (2 * NMII_N_SIDE)))
    n_actin = int(cell.n_actin)
    r0_backbone = float(myo.r0_backbone)
    bending = NMIIBendingStiffness(
        persistence_length_um=NMII_BACKBONE_LP_DIAGNOSTIC_UM, segment_len_um=r0_backbone,
        k_xb_pn_per_um=float(k_xb), r0_head_um=NMII_HEAD_OFFSET_UM,
    )
    head_side = np.tile(np.array([0] * NMII_N_SIDE + [1] * NMII_N_SIDE, np.int32), n_min)[:n_heads]
    return build_nmii_actuator_state(
        n_minifilaments=n_min, n_heads=n_heads,
        position_d=cell.pos_d,                                                    # ALIAS the global array
        minifilament_offset_d=wp.array((n_actin + np.arange(n_min + 1) * n_part).astype(np.int32),
                                       dtype=wp.int32, device=device),
        active_minifilament_d=wp.array(np.ones(n_min, np.int32), dtype=wp.int32, device=device),
        minifilament_id_d=wp.array(np.arange(n_min, dtype=np.int64), dtype=wp.int64, device=device),
        particle_role_d=wp.array(np.zeros(int(cell.n_total), np.int32), dtype=wp.int32, device=device),
        head_node_d=myo.head_node,                                               # GLOBAL head indices
        head_id_d=wp.array(np.arange(n_heads, dtype=np.int64), dtype=wp.int64, device=device),
        head_side_d=wp.array(head_side, dtype=wp.int32, device=device),
        backbone_bonds_d=myo.backbone_bonds, head_bonds_d=myo.head_bonds,
        backbone_angles_d=myo.backbone_angles, head_arm_angles_d=myo.head_arm_angles,
        k_backbone_pn_per_um=float(myo.k_backbone), r0_backbone_um=r0_backbone,
        k_head_spring_pn_per_um=float(myo.k_head_spring), r0_head_um=float(myo.params.r0_head),
        bending=bending, device=device,
    )


def _build_port(cell, device):
    """Build the cortex actin bind-target port over the LIVE segment topology (aliasing the global array)."""
    from aleph.engine.nmii_actuator import CORTEX_COMPONENT, FilamentMotorPortView

    sr = cell.myosin.segment_runtime
    n_seg = int(sr.seg_node_a.shape[0])
    seg_a_host = sr.seg_node_a.numpy()
    fiber_of_node = cell.node_fiber_d.numpy() if cell.node_fiber_d is not None else np.zeros(cell.n_actin, np.int64)
    # persistent identity: the fiber each segment belongs to; material arc-coords are structural placeholders
    # (the KMC/crossbridge kernels do not consume them — they carry persistent identity for a future remap).
    return FilamentMotorPortView(
        component=CORTEX_COMPONENT,
        position_d=cell.pos_d, force_d=cell.f_d,                                  # ALIAS the global array
        segment_node_a_d=sr.seg_node_a, segment_node_b_d=sr.seg_node_b, segment_polarity_d=sr.seg_polarity,
        persistent_filament_id_d=wp.array(fiber_of_node[seg_a_host].astype(np.int64),
                                          dtype=wp.int64, device=device),
        material_s0_d=wp.array(np.arange(n_seg, dtype=np.float64), dtype=wp.float64, device=device),
        material_s1_d=wp.array(np.arange(1, n_seg + 1, dtype=np.float64), dtype=wp.float64, device=device),
        topology_epoch_d=wp.array(np.zeros(1, np.int32), dtype=wp.int32, device=device),
    )


class _SliceMyosinAdapter:
    """Replaces ``cell.myosin`` so the driver's inner relax drives the SLICE's motor (one binding SoA).

    ``_accumulate_all`` calls ``cell.myosin.accumulate(pos, f)`` each inner iteration; this adapter adds the
    minifilament's own internal backbone/head-arm mechanics + the slice's split crossbridge over WHATEVER
    ``pos``/``f`` the inner solve hands it (array-agnostic — it builds a view/port over those arrays).  The
    crossbridge is zero until the KMC binds heads, so the cortex only feels the motor once binding EMERGES.
    """

    segment_runtime = None  # the slice's transaction owns the KMC commit; the driver's must not also run it.

    def __init__(self, slc, original) -> None:
        self._mechanics = slc.actuator_state.mechanics
        self._connector = slc.connector
        self._view0 = slc.actuator_state.geometry()
        self._port0 = slc.port
        self._original = original                    # defensive: delegate any other MyosinForce attr the driver reads

    def __getattr__(self, name: str):
        # only reached when normal lookup fails (accumulate/segment_runtime are found first).
        return getattr(object.__getattribute__(self, "_original"), name)

    def accumulate(self, pos: wp.array, f: wp.array) -> None:
        self._mechanics.accumulate_internal(pos, f)                              # backbone rod + head arm + F6
        view = dataclasses.replace(self._view0, position_d=pos, force_d=f)
        port = dataclasses.replace(self._port0, position_d=pos, force_d=f)
        self._connector.accumulate_candidate(view, port)                         # split crossbridge (+f / −f)


def _crossbridge_family(slc, pos):
    """Bound-head crossbridge dipoles → (head_pos, attach_pos, tangential tension) for γ_source, else None.

    The tension uses the PER-HEAD zero-strain reference ``r0_bind`` (attach-unstrained), the SAME quantity the
    device crossbridge scatters onto the cortex — so γ_source is the real ACTIVE power-stroke tension, not the
    removed coarse-mesh placement artifact.
    """
    st = slc.connector.state
    idx = np.nonzero(st.bound_d.numpy() == 1)[0]
    if idx.size == 0:
        return None
    hn = slc.connector.head_node_d.numpy()[idx]
    sa, sb = st.seg_a_d.numpy()[idx], st.seg_b_d.numpy()[idx]
    t = st.bary_t_d.numpy()[idx][:, None]
    absc = st.abscissa_d.numpy()[idx][:, None]
    wd = st.walk_dir_d.numpy()[idx]
    r0_bind = slc.connector._r0_bind_d.numpy()[idx]           # per-head attach-unstrained reference
    head = pos[hn]
    attach = (1.0 - t) * pos[sa] + t * pos[sb] + absc * wd
    k_xb = float(slc.connector.params.k_xb)
    tension = k_xb * (np.einsum("ij,ij->i", attach - head, wd) - r0_bind)
    return head, attach, tension


def _measure_gamma(cell, slc, n_planes: int = 64, omit: frozenset[str] | tuple[str, ...] = ()):
    """Method-of-planes γ (source / network / total) on the live cortex — an out-of-hot-loop device readback.

    ``omit`` MUST be the mask the inner solve runs with.  This is not a tidiness argument: the measurement
    assembles the force field independently of the solve, so a mask applied to one and not the other makes
    the estimator read a DIFFERENT force field from the one being integrated.  Measured on 2026-07-29 —
    with ``--engine-cortex`` the A/B arms' dynamics were bit-identical (``max_f_cortex_pn`` and ``n_bound``
    agreed to the last digit at every sample) while γ differed by a CONSTANT 0.41526 pN/µm, because the
    engine's cortex channels arrive through the ``cell.myosin`` hook — which is not omittable — and the
    unmasked measurement therefore added the incumbent's bending/crosslink on top of the engine's.
    Dynamics right, measurement double-counted.  A constant offset with identical trajectories is the
    signature; if the relocation itself had been wrong the trajectories would have diverged.
    """
    dev = wp.get_device(cell.device)
    f = wp.zeros(cell.n_total, dtype=wp.vec3d, device=cell.device)
    _accumulate_all(cell, cell.pos_d, f, omit=omit)      # includes the slice crossbridge (via the adapter)
    wp.synchronize_device(dev)
    pos, f_ext = cell.pos_d.numpy(), f.numpy()
    n_actin, foff = int(cell.n_actin), cell.foff_d.numpy()
    sa, sb, t_actin = actin_axial_tension(f_ext[:n_actin], pos[:n_actin], foff)
    families: dict[str, tuple] = {"actin": (pos[sa], pos[sb], t_actin)}
    if int(getattr(cell, "n_xl", 0)):
        xa, xb, t_xl = crosslink_tension(pos, cell.xl_d.numpy(), cell.kxl_d.numpy(), cell.r0xl_d.numpy())
        families["crosslink"] = (pos[xa], pos[xb], t_xl)
    cb = _crossbridge_family(slc, pos)
    if cb is not None:
        families["crossbridge"] = cb
    return measure_cortical_stress(
        R_CORTEX_UM, families, pos[:n_actin], f_ext[:n_actin], n_planes=n_planes,
        source_families=("crossbridge",), network_families=("actin", "crosslink"))


def _bound_head_lifetime(slc, params: CortexMotorParams) -> dict[str, float | int | str]:
    """Return this run's driving correlation time by asking the SHARED kinetics closed form.

    The arithmetic lives in :func:`~aleph.components.motor.bell_kinetics_analytic.bound_state_lifetime`, not
    here: a driver-local copy is how the literal ``1/0.4`` — wrong in all three of its numbers for the
    runs it stamped — survived, and any other lane needing a driving time should get the fixed version
    for free rather than re-deriving it.  This wrapper only reads the bound heads' loads off the device.
    """
    state = slc.connector.state
    bound = state.bound_d.numpy() == 1
    return bound_state_lifetime(
        state.loads_bell_d.numpy()[bound],
        catch_slip=params.detach_kinetics is SegmentDetachKinetics.CATCH_SLIP,
        k_off0=params.k_off0, f0=params.f0,
        k_catch0=params.k_catch0, x_catch=params.x_catch,
        k_slip0=params.k_slip0, x_slip=params.x_slip, kT=params.kT,
    )


def _assess_series(trajectory: list[dict], sample_dt_s: float, *, driving_tau_s: float,
                   min_windows: float, drift_sigma: float) -> dict[str, dict]:
    """Judge every recorded observable — a thin call into the shared engine helper."""
    return assess_observables(
        trajectory, sample_dt_s, observables=STATIONARITY_OBSERVABLES,
        driving_correlation_time_s=driving_tau_s,
        min_windows=min_windows, drift_sigma=drift_sigma)


def _dump_viz(cell, slc, *, mode: str, n_heads: int, args, out_path: str,
              omit: frozenset[str] | tuple[str, ...] = ()) -> str:
    """Extract the dynamic steady-state cortex + bound-NMII into ONE v2 ``.npz`` for the dev-Mac viewer.

    Host-side only, OUT of the hot loop (I0-A): one device readback (``.numpy()`` after the step loop) of
    arrays the slice/cell already own on device.  Encodes, through the SINGLE v2 writer
    (:func:`aleph.scripts.dump_state.write_v2_dump`), so ``ac_cell_assembled_viz.py`` reads it with no
    special-casing:

      (a) **cortex actin** as a ``filaments`` component whose per-node scalar ``fmag`` is the composed per-node
          |F| [pN] — the SAME encoding as ``ac/cell/dump_state``'s ``f_total_post`` (a fresh ``_accumulate_all``
          → per-node norm), so the emergent-tension field colours exactly like the assembled dump.  It also
          carries the per-node axial tension (the method-of-planes γ source) as an extra ``axial_tension`` field.
      (b) the **bound NMII minifilaments** as their own ``lines_raw`` component (backbone line verts + head-arm
          verts) carrying which heads are BOUND (``head_bound``), the head node positions (``head_pos``), and the
          per-head Bell crossbridge load (``head_load``).
      (c) the **bound crossbridges** as an explicit ``MOTOR`` connector (nmii→cortex; endpoints head→actin,
          per-head crossbridge tension as the load) — the force PATH that produced the tension.

    Plus enough metadata (bound %, γ_source/network/total) to title the scene.  Reuses the driver's own
    :func:`_measure_gamma` / :func:`_crossbridge_family` / :func:`actin_axial_tension`.
    """
    from pathlib import Path

    from aleph.scripts.ac_viz_common import seg_pairs_from_offsets
    from aleph.scripts.dump_state import _c, _k, write_v2_dump

    dev = wp.get_device(cell.device)

    # (0) one out-of-loop composed-force readback (matches ac/cell/dump_state f_total_post) + geometry.
    f = wp.zeros(cell.n_total, dtype=wp.vec3d, device=cell.device)
    _accumulate_all(cell, cell.pos_d, f, omit=omit)            # crossbridge + crosslink + bending + steric + …
    wp.synchronize_device(dev)
    pos = cell.pos_d.numpy()
    f_ext = f.numpy()
    n_actin = int(cell.n_actin)
    foff = cell.foff_d.numpy().astype(np.int64)

    # (a) cortex actin — per-node composed |F| (the emergent-tension field) + per-node axial tension (γ source).
    pos_actin = pos[:n_actin]
    seg = seg_pairs_from_offsets(foff, n_actin)
    fmag = np.linalg.norm(f_ext[:n_actin], axis=1)
    sa_ax, sb_ax, t_ax = actin_axial_tension(f_ext[:n_actin], pos_actin, foff)
    node_tension = np.zeros(n_actin, np.float64)
    if t_ax.size:
        np.maximum.at(node_tension, sa_ax, np.abs(t_ax))       # each node ← max |axial tension| of its segments
        np.maximum.at(node_tension, sb_ax, np.abs(t_ax))
    fib_of_node = np.zeros(n_actin, np.int64)
    fib_of_node[foff[:-1]] = 1
    fib_of_node = np.cumsum(fib_of_node) - 1                    # global-unique actor id = owning fiber

    n_bound = int(slc.bound_head_count())
    bound_pct = 100.0 * n_bound / max(n_heads, 1)
    rep = _measure_gamma(cell, slc, omit=omit)                 # method-of-planes γ for the scene title

    cortex = _c("cortex", f"cortex F-actin · emergent γ ({n_actin:,} nodes)", "filaments", pos_actin,
                seg=seg, fmag=fmag, actor_id=fib_of_node, order=0,
                array_meta={"axial_tension": node_tension},
                meta={"n_actors": int(cell.n_fibers), "field": "per-node composed |F| [pN]",
                      "n_bound_heads": n_bound})

    # (b) NMII minifilaments — backbone line verts + head arms + which heads BOUND + per-head crossbridge load.
    mech = slc.actuator_state.mechanics
    bb = mech.backbone_bonds_d.numpy().astype(np.int64).reshape(-1, 2)     # GLOBAL node pairs into pos
    hb = mech.head_bonds_d.numpy().astype(np.int64).reshape(-1, 2)
    head_node = slc.actuator_state.head_node_d.numpy().astype(np.int64)
    bound = slc.connector.state.bound_d.numpy().astype(np.int32)
    head_load = np.abs(slc.connector.state.loads_bell_d.numpy())          # per-head Bell |F| load [pN]
    nmii = _c("nmii", f"NMII minifilaments · {n_bound:,}/{n_heads:,} heads bound ({bound_pct:.0f}%)",
              "lines_raw", pos[bb].reshape(-1, 3), order=1,
              array_meta={"head_verts": pos[hb].reshape(-1, 3), "head_pos": pos[head_node],
                          "head_bound": bound, "head_load": head_load},
              meta={"n_actors": int(cell.ledger.get("myosin_n_minifilaments", 0)), "n_heads": int(n_heads),
                    "n_bound": n_bound})

    # (c) bound crossbridges → an explicit MOTOR connector (nmii→cortex), load = per-head crossbridge tension.
    cons: list[dict] = []
    cb = _crossbridge_family(slc, pos)
    if cb is not None:
        head_cb, attach_cb, tension_cb = cb
        ep = np.stack([head_cb, attach_cb], axis=1).reshape(-1, 3)         # [head0,attach0, head1,attach1, …]
        cons.append(_k("nmii_cortex_motor", "MOTOR", f"NMII→cortex crossbridge · {tension_cb.shape[0]:,} bound",
                       "nmii", "cortex", ep, load=np.abs(tension_cb)))

    # Two-axis evidence (PI D1-B, 2026-07-28). The rung is STRUCTURAL and this run earns CONNECTED:
    # `nmii_cortex_motor` is genuinely dispatched and the head binding it reports is generated by k_on
    # Poisson events rather than seeded. The MAGNITUDES are a separate axis and stay BLOCKED — every
    # NMII parameter here is a PI-GAP, the crossbridge term is ~0.5% of the reported gamma, and the
    # steps are force-accepted rather than converged. Emitting one bare "CONNECTED" string is what let
    # this artifact previously read as though its gamma were quotable.
    evidence = EvidenceLabel(
        rung=EvidenceRung.CONNECTED,
        quantitative=QuantitativeClaim.BLOCKED,
        basis=(
            f"nmii_cortex_motor dispatched; {n_bound:,}/{int(n_heads):,} heads bound from k_on events "
            f"over {int(args.steps)} steps at dt_phys={float(args.dt)}"
        ),
    )
    stage = {"name": f"GATE-B dynamic cortex-motor (N={n_actin:,})", **evidence.as_artifact_fields(),
             "bound_pct": bound_pct, "n_bound_heads": n_bound, "n_heads": int(n_heads),
             "gamma_total_pn_per_um": float(rep.gamma_total_pn_per_um),
             "gamma_source_pn_per_um": float(rep.gamma_source_pn_per_um),
             "gamma_network_pn_per_um": float(rep.gamma_network_pn_per_um),
             "detach": mode, "steps": int(args.steps), "dt_phys": float(args.dt), "k_xb": float(args.k_xb)}
    # slim, JSON-clean ledger (write_v2_dump's manifest json.dumps has no default=float coercion).
    ledger = {"myosin_n_minifilaments": int(cell.ledger.get("myosin_n_minifilaments", 0)),
              "myosin_n_heads": int(n_heads), "myosin_n_bound": n_bound}
    meta = {"stage": stage, "report": dict(stage), "ledger": ledger, "R_cell": float(R_CORTEX_UM)}
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    out = write_v2_dump(out_path, [cortex, nmii], cons, meta)
    print(f"[gate-b] dump-viz wrote {out}\n"
          f"          cortex |F| field ({n_actin:,} nodes) + {n_bound:,}/{n_heads:,} bound NMII heads "
          f"({bound_pct:.1f}%) + {len(cons)} MOTOR connector · γ_total={rep.gamma_total_pn_per_um:.4g} pN/µm",
          flush=True)
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description="GATE-B: emergent cortical NMII tension on the native cortex.")
    ap.add_argument("--steps", type=int, default=40, help="number of ACCEPTED physical steps to run")
    ap.add_argument("--dt", type=float, default=0.01, help="physical timestep dt_phys [s] (attach/detach tick)")
    ap.add_argument("--outer", type=int, default=40, help="inner mechanical relax iterations per step (n_inner)")
    ap.add_argument("--catch-slip", action="store_true",
                    help="use the physiological Pereverzev catch-slip detach (proxy constants) instead of Bell slip")
    ap.add_argument("--balance-gate", action="store_true",
                    help="gate acceptance on the force-balance predicate instead of force-accepting every step")
    ap.add_argument("--k-xb", type=float, default=NMII_K_XB_TEST, dest="k_xb",
                    help="crossbridge stiffness [pN/µm] override (default = params_i0b3 provisional 1000); "
                         "sweep to confirm the OLD residual scales ∝ k_xb but the r0_bind-FIXED one does not")
    ap.add_argument("--filaments", type=int, default=70686, help="cortex F-actin count (native = 70686)")
    ap.add_argument("--cortex-seg-um", type=float, default=0.5, dest="cortex_seg_um",
                    help="cortex segment rest length ℓ₀ [µm] → mesh/pore size (baseline 0.5; fine mesh 0.075)")
    ap.add_argument("--cortex-density", type=float, default=20.0, dest="cortex_density",
                    help="crosslinks / filament (baseline 20; fine mesh 40)")
    ap.add_argument("--cortex-length-um", type=float, default=3.0, dest="cortex_length_um",
                    help="representative cortical filament contour length L [µm] (baseline 3.0)")
    ap.add_argument("--cortex-arp23-fraction", type=float, default=0.0, dest="cortex_arp23_fraction",
                    help="Arp2/3 fraction of cortical actin BY MASS (~0.33 Bovellan; 0.0 = formin-only, "
                         "bit-identical baseline). >0 splits the fixed --filaments budget into formin + short "
                         "branched Arp2/3 (areal density preserved). NOTE: _accumulate_all does NOT launch "
                         "branch_angle_kernel, so the 70° junctions carry only the anchor spring, not the "
                         "angle-harmonic restoring force (see report).")
    ap.add_argument("--membrane-subdiv", type=int, default=6, dest="membrane_subdiv",
                    help="plasma-membrane icosphere subdivision level. DEFAULT 6 = 40,962 verts (~131 nm "
                         "spacing, ~4 nodes across a 0.5 µm bleb neck) — the PI-RATIFIED production baseline "
                         "(assemble.py:267, PI 2026-07-22, reaffirmed 2026-07-29); 7 = 163,842 (~66 nm) is the "
                         "ratified VALIDATION resolution and costs 1.24x. This default was 8 (655,362 verts) "
                         "until 2026-07-29 — a level in no ratification, costing 2.20x for resolution nothing "
                         "asked for. The memory argument that came with it does not survive measurement: peak "
                         # `%%` because argparse runs every help string through `% params`, and a bare
                         # `% f` is a conversion specifier that wants a number and gets a dict. This is
                         # what made `--help` raise TypeError and `test_gate_b_help_exposes_flags` fail.
                         "pool at subdiv 8, full native, is 1.222 GB against the A5000's 16.76 GB (92.7%% free), "
                         "so lower levels were never an 'A5000 memory safety valve' and are not offered as one")
    ap.add_argument("--seed", type=int, default=0, help="base RNG seed")
    ap.add_argument("--engine-cortex", action="store_true", dest="engine_cortex",
                    help="compute the cortex's own force channels through the ENGINE component "
                         "(ac/engine/cortex_state) instead of the incumbent's inline launches: the "
                         "incumbent omits them via _accumulate_all(omit=) and the component's mechanics "
                         "supplies them. Default OFF and byte-identical. This is what makes ac/engine "
                         "fire a kernel in a production run at all — until now every kernel in a native "
                         "step came from ac/cell, ac/motor, ac/fluid or ff/")

    # ── stationarity: the run judges its OWN settling, and the contract is declared here ─────────────
    # These four are the PRE-DECLARED contract (CLAUDE.md's no-gate-loosening rule): they are recorded
    # verbatim in the artifact next to the verdict, so a later reader judges the verdict against the
    # criterion that produced it rather than one chosen after seeing the series.
    ap.add_argument("--min-windows", type=float, default=5.0, dest="min_windows",
                    help="CONTRACT: shortest analysable span, in bound-head lifetimes, that may return "
                         "anything but TOO_SHORT (default 5). Declare it before the run")
    ap.add_argument("--drift-sigma", type=float, default=1.0, dest="drift_sigma",
                    help="CONTRACT: stationary requires |drift across the window| < this many sigma "
                         "of the series' own fluctuation (default 1.0)")
    ap.add_argument("--driving-tau-s", type=float, default=0.0, dest="driving_tau_s",
                    help="override the driving correlation time [s]. Default 0 = derive it from THIS "
                         "run's own detach kinetics (see _bound_head_lifetime), taking the larger of "
                         "the zero-load and the mean-load value")
    ap.add_argument("--gamma-every", type=int, default=1, dest="gamma_every",
                    help="record the gamma telemetry every Nth step (default 1 = every step). The "
                         "per-step readback is ~2/3 of this driver's wall-clock at native, so a stride "
                         "buys physical time at fixed cost; the sample spacing is N*dt and the "
                         "correlation-time estimator is given THAT, not dt")
    ap.add_argument("--until-stationary", action="store_true", dest="until_stationary",
                    help="ignore --steps as a target and keep stepping until gamma_total is judged "
                         "STATIONARY, bounded by --max-steps and --max-wall-minutes. This is the "
                         "protocol a swept parameter that CHANGES the relaxation time requires: a "
                         "fixed step count compares points at different times, which is how the "
                         "2026-07-28 blebbistatin sweep came out non-monotonic")
    ap.add_argument("--max-steps", type=int, default=0, dest="max_steps",
                    help="hard step budget for --until-stationary (default 0 = 20x --steps)")
    ap.add_argument("--max-wall-minutes", type=float, default=0.0, dest="max_wall_minutes",
                    help="hard wall-clock budget for --until-stationary [min] (default 0 = none). The "
                         "run stops and records NOT_SETTLED rather than being killed by the lease")
    ap.add_argument("--out", type=str, default="",
                    help="run-record@2 artifact path; the trajectory is written there")
    ap.add_argument("--build-commit", default=None,
                    help="build the run is measured on, supplied by ffn_gpu.py run")
    ap.add_argument("--dump-viz", type=str, default=None, metavar="PATH", dest="dump_viz",
                    help="after the loop reaches steady state, write ONE v2 .npz (cortex |F| field + BOUND NMII "
                         "heads/crossbridges) for the dev-Mac viewer ac_cell_assembled_viz.py; out-of-loop "
                         "host readback only (I0-A safe)")
    args = ap.parse_args()

    wp.init()
    dev = str(wp.get_device())
    if not wp.get_device().is_cuda:
        raise RuntimeError(f"GATE-B native driver needs a CUDA GPU (I0-A). Resolved {dev!r} is not CUDA.")
    mode = "CATCH_SLIP (Kovacs 2007, PROVISIONAL proxy)" if args.catch_slip else "SLIP (Bell, provisional-sourced)"
    mesh_label = "FINE" if args.cortex_seg_um < 0.5 or args.cortex_density > 20.0 else "COARSE baseline"
    print(f"[gate-b] device={dev}  steps={args.steps}  dt={args.dt}s  inner={args.outer}  k_xb={args.k_xb}  "
          f"detach={mode}  (crossbridge attach-UNSTRAINED: r0_bind removes the coarse-mesh placement artifact)",
          flush=True)
    print(f"[gate-b] cortex mesh: seg_um={args.cortex_seg_um}  density_per_fil={args.cortex_density}  "
          f"length_um={args.cortex_length_um}  ({mesh_label})", flush=True)

    # 1. Build the FULL native cell — myosin STRUCTURE present, ALL heads UNBOUND (resting seed OFF).
    t0 = time.time()
    cfg = CellConfig(
        n_filaments=args.filaments, with_myosin=True, overlap_free_cortex=True,
        membrane_subdivisions=args.membrane_subdiv, nucleus_subdivisions=3, erm_radial_pairing=True,
        cortex_seg_um=args.cortex_seg_um,                         # ← mesh resolution (baseline 0.5; fine 0.075)
        cortex_density_per_fil=args.cortex_density,               # ← crosslinks/filament (baseline 20; fine 40)
        cortex_length_um=args.cortex_length_um,
        cortex_arp23_fraction=args.cortex_arp23_fraction,         # ← 0.0 formin-only; >0 mixed formin+Arp2/3
        resting_bound_myosin_fraction=None,                       # ← heads UNBOUND at t0 (purely dynamic)
        nmii_straddle_placement=True,                             # heads placed ON actin (~nm crossbridge)
        nmii_backbone_lp_um=NMII_BACKBONE_LP_DIAGNOSTIC_UM,       # materialise F6 angle topology (see header)
        nmii_backbone_lp_source="MECHANISM_DEMO_DIAGNOSTIC_L_p_FIXTURE_NOT_PRODUCTION_PI_GAP",
        seed=args.seed,
    )
    cell = build_cell(cfg)
    n_heads = int(cell.myosin.head_node.shape[0])
    print(f"[gate-b] built cell in {time.time() - t0:.1f}s: n_actin={cell.n_actin} n_total={cell.n_total} "
          f"n_heads={n_heads} n_minifilaments={cell.ledger.get('myosin_n_minifilaments')} "
          f"straddle_frac={cell.ledger.get('myosin_straddle_placed_frac', 0.0):.2f} "
          f"resting_bound={cell.ledger.get('resting_bound_myosin_n_bound')}", flush=True)

    # 2. Extract the actuator + cortex port, build the slice (aliasing the global pos_d/f_d).
    params = _params(args.catch_slip, args.k_xb)
    _banner(params, mode)
    actuator = _build_actuator(cell, dev, args.k_xb)
    port = _build_port(cell, dev)
    max_seg_um = float(cell.srest_d.numpy().max())
    slc = build_cortex_motor_slice(
        actuator_state=actuator, port=port, params=params, base_seed=args.seed, device=dev,
        max_segment_length_um=max_seg_um)

    # 3. Route the slice's motor through the driver's own inner relax (one binding SoA, no double count).
    cell.myosin = _SliceMyosinAdapter(slc, original=cell.myosin)
    # ── optional: hand the cortex's force channels to the ENGINE component ─────────────────────────
    # The incumbent stops launching them (`omit=`) and `ac/engine/cortex_state` supplies them from its
    # own mechanics, over the SAME topology built from the SAME seed — which is why the gate that
    # measured node-by-node parity at 2.80e-16 licenses this swap rather than merely accompanying it.
    #
    # HONEST SCOPE. The arrays are still shared: the component computes the cortex channels but writes
    # into the incumbent's `f`, because a genuinely private force array needs the inner solve to relax
    # TWO arrays together, and that is `SOLVE_COUPLED` — an empty phase (STATE.md (c) 4). So this closes
    # "ac/engine fires no kernel" and gives the cortex an owning object; it does NOT end co-location.
    cortex_owner = None
    omit_channels: frozenset[str] = frozenset()
    if args.engine_cortex:
        from aleph.components.incumbent.driver import OMITTABLE_CHANNELS
        from aleph.engine.cortex_population import build_cortex_population
        from aleph.engine.cortex_state import build_cortex_state_owner, incumbent_channels_to_omit

        # EVERY parameter the incumbent's cortex is built from must be forwarded, or the two are two
        # different cortices and omitting the incumbent's channels DROPS force instead of relocating it.
        # `arp23_fraction` in particular splits the filament budget into more, shorter filaments — the node
        # count check below would catch that, but a guard that only fires on a count is not a reason to
        # leave the parameter unforwarded.
        population = build_cortex_population(
            int(args.filaments), seg_um=args.cortex_seg_um, length_um=args.cortex_length_um,
            density_per_fil=args.cortex_density, arp23_fraction=args.cortex_arp23_fraction,
            arp23_length_um=cfg.cortex_arp23_length_um, arp23_seg_um=cfg.cortex_arp23_seg_um,
            arp23_mother_fraction=cfg.cortex_arp23_mother_fraction,
            overlap_free=cfg.overlap_free_cortex, overlap_mode=cfg.cortex_overlap_mode,
            overlap_span=cfg.cortex_overlap_span, seed=int(args.seed))
        cortex_owner = build_cortex_state_owner(population.topology, device=dev)
        if int(cortex_owner.position_d.shape[0]) != int(cell.n_actin):
            raise RuntimeError(
                f"engine cortex owns {int(cortex_owner.position_d.shape[0]):,} nodes but the incumbent's "
                f"actin block is {int(cell.n_actin):,} — two different cortices, so omitting the "
                f"incumbent's channels would drop force rather than relocate it")
        # DERIVED from what the owner actually bound, never typed here. A hand-written pair is right for the
        # formin-only default and silently DOUBLES the Arp2/3 branch-angle force for a mixed cortex
        # (`--cortex-arp23-fraction > 0`), because the owner binds that channel while the incumbent is never
        # told to stop launching it. Deriving it makes the two sets one set.
        omit_channels = incumbent_channels_to_omit(cortex_owner.mechanics)
        assert omit_channels <= OMITTABLE_CHANNELS
        if len(omit_channels) != len(cortex_owner.mechanics.channels):
            raise RuntimeError(
                f"engine cortex binds {cortex_owner.mechanics.channels} but only {sorted(omit_channels)} "
                "would be omitted — two engine channels collapsed onto one incumbent name, so some force "
                "would be counted twice")
        print(f"[gate-b] ENGINE CORTEX: {int(cortex_owner.position_d.shape[0]):,} nodes, channels "
              f"{cortex_owner.mechanics.channels}; incumbent omits {sorted(omit_channels)}", flush=True)

    inner_solve = make_inner_solve(cell, n_inner=args.outer, reshape_every=20, inner_solver="explicit",
                                   omit=omit_channels)

    if cortex_owner is not None:
        # The adapter already routes the slice's motor through the inner relax; chain the engine cortex
        # onto the same hook so its channels are summed at exactly the point the incumbent's were.
        _slice_adapter = cell.myosin
        _engine_mechanics = cortex_owner.mechanics

        class _CortexPlusMotor:
            """Supply the ENGINE cortex channels alongside the motor, on whatever arrays solve hands us."""

            segment_runtime = None

            def __getattr__(self, name: str):
                return getattr(_slice_adapter, name)

            def accumulate(self, pos: wp.array, f: wp.array) -> None:
                _engine_mechanics.accumulate(pos, f)      # ac/engine computes the cortex, in production
                _slice_adapter.accumulate(pos, f)

        cell.myosin = _CortexPlusMotor()

    def cortex_relax() -> None:
        # The driver's inner mechanical relax: bending + inextensibility + crosslink + steric + turgor +
        # nucleus + membrane/ERM + (via the adapter) the slice's internal mechanics + split crossbridge.
        inner_solve(float(args.dt))

    if args.balance_gate:
        raise SystemExit(
            "--balance-gate requires a GlobalCellLedger + tol_sq_d wired by the Lead (see "
            "build_cortex_motor_slice ledger/tol_sq_d args + the slice's ledger contributors); the "
            "MECHANISM demo path is force-accept, which lets binding accumulate so γ can be seen to rise.")
    accepted_ones = wp.array(np.ones(1, np.int32), dtype=wp.int32, device=dev)

    # 4. Accepted-step event loop — γ EMERGES from the binding events, and the run judges its own settling.
    sample_dt = float(args.dt) * max(1, int(args.gamma_every))
    max_steps = int(args.max_steps) if args.max_steps > 0 else (20 * int(args.steps)
                                                               if args.until_stationary else int(args.steps))
    wall_budget_s = float(args.max_wall_minutes) * 60.0 if args.max_wall_minutes > 0.0 else float("inf")
    lifetime = _bound_head_lifetime(slc, params)          # pre-run: nothing bound yet, so zero-load only
    driving_tau = (float(args.driving_tau_s) if args.driving_tau_s > 0.0
                   else float(lifetime["conservative_s"]))
    print(f"\n[gate-b] stationarity contract: min_windows={args.min_windows:g} x driving tau, "
          f"drift_sigma={args.drift_sigma:g}; driving tau = {driving_tau:.4g} s "
          f"({'CLI override' if args.driving_tau_s > 0.0 else lifetime['derivation']})")
    print(f"[gate-b] required analysable span = {args.min_windows * driving_tau:.4g} s of physical time "
          f"AFTER the transient = {args.min_windows * driving_tau / float(args.dt):.0f} steps at dt="
          f"{args.dt}; sampling every {max(1, int(args.gamma_every))} step(s) (sample dt {sample_dt:g} s)"
          + (f"; budget {max_steps} steps" if args.until_stationary else ""), flush=True)
    print(f"\n{'step':>6} {'t_s':>7} {'bound':>7} {'bound%':>7} {'g_source':>10} {'g_network':>10} "
          f"{'g_total':>10} {'maxPF_all':>10} {'maxF_ctx':>10}", flush=True)
    # THE TRAJECTORY IS THE RESULT, so it is kept rather than printed and dropped. Until 2026-07-28 this
    # driver wrote NO artifact at all — the "tension emerges 0 -> 3.72 pN/um" headline exists only as
    # stdout, which is why STATE.md records "None in a clean checkout" for it. A per-step row costs
    # nothing next to a 7.7 s step and is the difference between a run that can be re-judged and one
    # that has to be re-run.
    trajectory: list[dict] = []
    #: Assessments taken WHILE the loop ran, kept so the stopping decision is auditable rather than
    #: inferable only from the final verdict.
    assessment_log: list[dict] = []
    stop_reason = "step target reached"
    loop_t0 = time.perf_counter()
    step = 0
    while step < max_steps:
        slc.step(cortex_relax, dt_phys=float(args.dt), accepted_d=accepted_ones)  # force-accept: watch γ rise
        step += 1
        # `step == 1` is sampled UNCONDITIONALLY, whatever the stride. Without it `trajectory[0]` is the
        # first SAMPLED step, not the first accepted one — and `t0` below reads `trajectory[0]` while its
        # own note says "the state after one dt". Measured 2026-07-29 on one run: at `--gamma-every 5`
        # that mislabels step 4 as t0, where n_bound is 2.31x and gamma_total 17.9x their step-0 values.
        if step != 1 and step % max(1, int(args.gamma_every)) and step != max_steps:
            continue

        # out-of-hot-loop telemetry (device→host readback BETWEEN steps only; I0-A honoured inside the loop)
        n_bound = slc.bound_head_count()
        rep = _measure_gamma(cell, slc, omit=omit_channels)
        max_pf = float(inner_solve.convergence_force_d.numpy()[0])   # projected residual (whole-cell max|PF|)
        fh = cell.f_d.numpy()[: cell.n_actin]
        max_f_ctx = float(np.nanmax(np.linalg.norm(fh, axis=1))) if fh.size else 0.0
        trajectory.append({
            "step": int(step - 1), "t_s": float(step * args.dt),
            "n_bound": int(n_bound), "bound_fraction": float(n_bound) / float(max(n_heads, 1)),
            "gamma_source_pn_per_um": float(rep.gamma_source_pn_per_um),
            "gamma_network_pn_per_um": float(rep.gamma_network_pn_per_um),
            "gamma_total_pn_per_um": float(rep.gamma_total_pn_per_um),
            "max_pf_pn": max_pf, "max_f_cortex_pn": max_f_ctx,
        })
        print(f"{step - 1:>6d} {step * args.dt:>7.2f} {n_bound:>7d} "
              f"{100.0 * n_bound / max(n_heads, 1):>6.1f}% "
              f"{rep.gamma_source_pn_per_um:>10.4g} {rep.gamma_network_pn_per_um:>10.4g} "
              f"{rep.gamma_total_pn_per_um:>10.4g} {max_pf:>10.4g} {max_f_ctx:>10.4g}", flush=True)

        if not args.until_stationary:
            if step >= int(args.steps):
                break
            continue

        # ── the stopping decision, taken on the series rather than on a step count ────────────────
        # Re-derive the driving time each assessment: with heads now bound the catch-bond lifetime is
        # measurable and is LONGER than the zero-load value, so the requirement can only tighten.
        if len(trajectory) >= 8 and len(trajectory) % 8 == 0:
            lifetime = _bound_head_lifetime(slc, params)
            if args.driving_tau_s <= 0.0:
                driving_tau = max(driving_tau, float(lifetime["conservative_s"]))
            reports = _assess_series(trajectory, sample_dt, driving_tau_s=driving_tau,
                                     min_windows=float(args.min_windows),
                                     drift_sigma=float(args.drift_sigma))
            gamma = reports["gamma_total_pn_per_um"]
            assessment_log.append({"step": int(step), "t_s": float(step * args.dt),
                                   "driving_tau_s": driving_tau,
                                   "verdict": gamma["verdict"], "n_samples": gamma["n_samples"],
                                   "drift_over_std": gamma.get("drift_over_std")})
            print(f"       ├─ stationarity[gamma_total] = {gamma['verdict']}  "
                  f"(tau_drive {driving_tau:.3g} s, {gamma['n_samples']} samples, "
                  f"drift/std {gamma.get('drift_over_std', float('nan')):.3g})", flush=True)
            if gamma["verdict"] == StationarityVerdict.STATIONARY.value:
                stop_reason = "gamma_total judged STATIONARY"
                break
        if time.perf_counter() - loop_t0 > wall_budget_s:
            stop_reason = f"wall budget {args.max_wall_minutes:g} min exhausted before settling"
            break
    else:
        stop_reason = f"step budget {max_steps} exhausted before settling"

    loop_s = time.perf_counter() - loop_t0
    n_steps_run = step

    # Final judgement, on the conservative driving time — measured now that heads are bound.
    lifetime = _bound_head_lifetime(slc, params)
    if args.driving_tau_s <= 0.0:
        driving_tau = max(driving_tau, float(lifetime["conservative_s"]))
    stationarity = _assess_series(trajectory, sample_dt, driving_tau_s=driving_tau,
                                 min_windows=float(args.min_windows),
                                 drift_sigma=float(args.drift_sigma))

    print(f"\n[gate-b] DONE: {n_steps_run} accepted steps ({n_steps_run * float(args.dt):.3g} s physical, "
          f"{n_steps_run * float(args.dt) / driving_tau:.2f} bound-head lifetimes). "
          f"Bound {slc.bound_head_count()}/{n_heads} heads "
          f"({100.0 * slc.bound_head_count() / max(n_heads, 1):.1f}%).  stop: {stop_reason}", flush=True)
    print("[gate-b] stationarity verdicts (RECORDED, not yet gate-authoritative — the amendment that "
          "would make settling the acceptance criterion is PROPOSED, PI ratification pending):", flush=True)
    for name in STATIONARITY_OBSERVABLES:
        rec = stationarity[name]
        mean = rec.get("mean")
        shown = (f"mean {mean:.6g} +/- {rec['sem']:.3g}" if mean is not None
                 else "no mean returned — a mean here would be the average of a transient")
        print(f"    {name:<26} {rec['verdict']:<11} {shown}", flush=True)
    record = observation_artifact(
        run_label="gate_b_cortex_motor — cortical tension trajectory from NMII binding events",
        evidence=EvidenceLabel(
            rung=EvidenceRung.CUDA_UNIT,
            quantitative=QuantitativeClaim.BLOCKED,
            basis=("recorded the per-sample cortical tension, its source/network split and the "
                   "bound-head count over accepted physical steps at the stated population, and judged "
                   f"each series for stationarity (gamma_total: "
                   f"{stationarity['gamma_total_pn_per_um']['verdict']}). BLOCKED because every NMII "
                   "magnitude here is a PI-GAP (k_xb, f_stall, kappa, N_side, L_bb) and every step is "
                   "force-accepted rather than admitted on the balance predicate"),
        ),
        config={
            "steps_requested": int(args.steps), "steps_run": int(n_steps_run),
            "dt_phys": float(args.dt), "n_inner": int(args.outer),
            "detach": mode, "k_xb": float(args.k_xb), "filaments": int(args.filaments),
            "cortex_seg_um": float(args.cortex_seg_um), "cortex_density": float(args.cortex_density),
            "cortex_length_um": float(args.cortex_length_um),
            "membrane_subdiv": int(args.membrane_subdiv), "seed": int(args.seed),
            "acceptance": "force-accept (accepted_d = ones); --balance-gate is not implemented",
            "engine_cortex": bool(args.engine_cortex),
            "omitted_channels": sorted(omit_channels),
            "gamma_every": int(args.gamma_every), "sample_dt_s": float(sample_dt),
            "until_stationary": bool(args.until_stationary),
            "max_steps": int(max_steps), "max_wall_minutes": float(args.max_wall_minutes),
            "stationarity_contract": {"min_windows": float(args.min_windows),
                                      "drift_sigma": float(args.drift_sigma),
                                      "driving_tau_s_override": float(args.driving_tau_s)},
        },
        census={
            # MEASURED from the built cell wherever the cell can be asked. `cortex_filaments` was
            # `args.filaments` — the REQUESTED count echoed back — in all 25 records audited 2026-07-29,
            # so the population stamp STATE.md requires on every result contained no measurement at all.
            "cortex_filaments": int(getattr(cell, "n_fibers", args.filaments)),
            "cortex_filaments_requested": int(args.filaments),
            "n_actin_nodes": int(cell.n_actin),
            "n_total_nodes": int(cell.n_total), "n_heads": int(n_heads),
            # The membrane is the compartment `fraction_of_native` could not see, and the one the
            # architectural principle names as the cost. Record what the subdivision actually produced.
            "membrane_vertices": int(10 * 4 ** int(args.membrane_subdiv) + 2),
            "membrane_subdiv": int(args.membrane_subdiv),
            # NAMED for what it measures. There is no whole-cell fraction here because the membrane's
            # physiological subdivision is sourced nowhere, so its denominator does not exist yet.
            "fraction_of_native_cortex": float(args.filaments) / 70686.0,
            "exact_peak_gpu_bytes": _peak_gpu_bytes(dev),
        },
        t0={
            "n_bound": int(trajectory[0]["n_bound"]) if trajectory else 0,
            "gamma_total_pn_per_um": float(trajectory[0]["gamma_total_pn_per_um"]) if trajectory else 0.0,
            "note": "the first ACCEPTED step; heads start unbound, so t0 is the state after one dt",
        },
        timing=timing_block(
            wall_seconds=loop_s, physical_time_s=float(n_steps_run) * float(args.dt),
            n_steps=int(n_steps_run), n_inner_iterations=int(n_steps_run) * int(args.outer),
            device_note=(f"telemetry (gamma readback + host force norm) runs BETWEEN steps every "
                         f"{max(1, int(args.gamma_every))} step(s) and is included in this wall-clock; "
                         "it is not free — at native it dominated the per-step cost until the stride "
                         "existed"),
        ),
        measurements={"trajectory": trajectory,
                      "physical_time_s": float(n_steps_run) * float(args.dt),
                      "sample_dt_s": float(sample_dt),
                      # DERIVED from this run's own kinetics, not the 1/0.4 literal this driver used to
                      # stamp — see _bound_head_lifetime for why all three numbers in that were wrong.
                      "bound_head_lifetime": lifetime,
                      "bound_head_lifetime_s": float(driving_tau),
                      "lifetimes_covered": float(n_steps_run) * float(args.dt) / driving_tau,
                      "stationarity": stationarity,
                      "stationarity_during_run": assessment_log,
                      "stop_reason": stop_reason},
        device=str(dev),
        # The LABEL and the VALUE, together.  Recording only "PI_GAP" makes the record un-re-analysable: on
        # 2026-07-29 the sweep's working-stroke validity (`f_stall / k_xb`) had to be recomputed from the
        # CURRENT build's constant, because no record carried the `f_stall` the run actually used.  A run made
        # under a different value would have been judged against today's and nothing would have said so.
        parameter_provenance={
            "k_xb": {"provenance": "PI_GAP", "value": float(args.k_xb), "units": "pN/um"},
            "f_stall": {"provenance": "PI_GAP", "value": float(NMII_F_STALL_TEST), "units": "pN"},
            "kappa_hill": {"provenance": "PI_GAP", "value": float(NMII_KAPPA_TEST), "units": "1"},
            "N_side": {"provenance": "PI_GAP", "value": int(NMII_N_SIDE), "units": "heads/side"},
            "L_bb": {"provenance": "PI_GAP", "value": int(NMII_N_BB), "units": "backbone beads"},
            "k_on": {"provenance": "PI_GAP", "value": float(NMII_KON_TEST), "units": "1/s"},
            "k_off0": {"provenance": "SOURCED", "value": float(NMII_KOFF0), "units": "1/s"},
            "_schema": "each entry is {provenance, value, units}; the value is what THIS run used",
        },
        notes={"not_a_claim": ("no steady state is asserted by the run label. Whether this trajectory "
                               "reached one is answered by measurements.stationarity, and the series is "
                               "here so it can be re-judged without re-running"),
               "verdict_is_not_yet_a_gate": (
                   "the stationarity verdicts are RECORDED, not acceptance criteria. The amendment that "
                   "would make settling what stage 1's outer gate scores is PROPOSED and PI ratification "
                   "is pending, so nothing here rejects a run; a DRIFTING or TOO_SHORT verdict says the "
                   "mean must not be quoted, which is a statement about the number and not about the run")},
        declared_commit=getattr(args, "build_commit", None),
    )
    out_path = Path(args.out) if getattr(args, "out", None) else Path(
        "outputs/ac/cell_assembled/gate_b_cortex_motor_trajectory.json")
    write_artifact(out_path, record)
    print(f"[gate-b] wrote {out_path}  ({len(trajectory)} sample rows over {n_steps_run} steps, "
          f"{record['measurements']['lifetimes_covered']:.2f} bound-head lifetimes)", flush=True)

    if args.dump_viz:
        # out-of-loop viz dump: host readback of the steady-state cortex + bound NMII (I0-A honoured — the loop
        # is done; no authoritative GPU→CPU roundtrip inside the physical-time loop).
        _dump_viz(cell, slc, mode=mode, n_heads=n_heads, args=args, out_path=args.dump_viz,
                  omit=omit_channels)


if __name__ == "__main__":
    main()
