#!/usr/bin/env python
r"""GATE-B NATIVE driver — a stress fiber's ACTIVE tension EMERGES from NMII binding events.

The SF counterpart of ``ac_gate_b_cortex_motor_native.py``, driving the ``nmii_sf_motor`` MOTOR edge that the
reference graph has declared since the architecture landed and that nothing had a runtime for::

    myosin heads bind by k_on EVENTS  →  the split crossbridge pulls the two anti-parallel SF filaments inward
    →  SF AXIAL TENSION and the FA-anchor traction EMERGE over accepted physical steps, starting from ZERO

There is NO lumped ``k_SF`` bundle spring and NO prestress seed: at t0 every straddle-placed head is UNBOUND and
the SF is at its force-free rest length, so any tension measured later was produced by the accepted-step event
runtime (:class:`~aleph.engine.sf_motor_slice.SFMotorSlice` driving
:class:`~aleph.engine.transaction.CellTransaction`).

WHAT IT SHOWS (all on-device; host reads are OUT-OF-LOOP diagnostics between steps only):
  1. **REST** — heads unbound, SF at rest: the composed candidate force ≈ 0 on BOTH owners (emergent-not-lumped).
  2. **EMERGENCE** — over accepted steps the bound-head count climbs 0 → N by ``k_on`` Poisson events, and the
     SF axial tension ``T = k_axial·(L−r0)`` plus the FA-anchor reaction climb WITH it (reported per step, so a
     flat or non-monotone trace is visible rather than hidden).
  3. **SIGN** — the ventral fibers' FA anchors are pinned, so the emergent traction pulls INWARD (contractile),
     never outward; the per-head loads sit at or below stall.
  4. **TRANSACTION** — one accepted step = snapshot(all three participants) → propose events → inner relax →
     per-head loads on the converged geometry → one device predicate to every rollback/commit → clock +1, with
     the ``sf_arc``/``nmii`` populations re-asserted disjoint every step.
  5. **ROLLBACK** — a REJECTED step advances no binding state and no clock (checked explicitly).
  6. **ACCEPTANCE IS PHYSICAL** — the predicate is no longer ``accepted_d = ones``.  Each body reduces its OWN
     force array into one half of the :class:`~aleph.engine.ledger.GlobalCellLedger` balance gate, whose
     device flag decides the step; the tolerance is DERIVED on-device from the ledger's own resultants (PI D8),
     so nothing is compared against a supplied force magnitude.  A **positive control** then breaks the adjoint
     pair by a known amount and shows the gate rejecting and the step bit-restoring — without it a gate that
     always passes is indistinguishable from a gate that cannot fail.
  7. **THE ENERGY LEDGER CLOSES** — ``ΔU + dissipated + event_jump = W_active`` per step, every term formed from
     the launches the step already performs (see :mod:`aleph.engine.observe.energy`), with the residual's
     own MOBILITY scaling saying whether it is the explicit integrator's first-order consistency error or a
     genuine leak.  No energy function is written anywhere: a second implementation of the forces would be free
     to drift from the forces it claims to differentiate.

WHY THE HEADS REACH THE ACTIN.  The sarcomere is built with the lateral separation and axial overlap DERIVED
from the minifilament that must sit in it (``lateral = 2·head_offset``, ``overlap = backbone_length``;
:func:`~aleph.engine.sf_motor_slice.sf_sarcomere_geometry_for`).  Host measurement for the default
population: head→filament-line residual median 0.000 µm (max 0.055) vs 0.200 µm for the legacy co-located
sarcomere, and the ``+``/``−`` head groups land on the two DIFFERENT filaments (bipolar dot ≈ −0.92).

SPLIT OWNERSHIP IS PHYSICAL HERE (unlike the cortex lane, where the actuator and the port alias one global
``cell.pos_d``): the minifilament particles and the SF nodes are two separate device arrays, so the crossbridge
is a genuine two-array adjoint transfer and both arrays are relaxed by the inner solve.

ONLY STRAIGHT SARCOMERES CARRY A MOTOR.  The perinuclear cap is an ARCH, so its two filaments are not
anti-parallel over the shared span and a rigid bipolar minifilament cannot straddle them (measured bipolar dot
−0.600 instead of −1).  Curved bundles are therefore excluded from the motor stations BY CONSTRUCTION and the
exclusion is reported by class in the run record — never silently approximated.

PARAMS ARE PI-GAPs (report-not-tune).  ``k_on``, ``f_stall``, ``v0``, ``kappa``, ``k_xb``, the catch-slip
constants, the SF axial ``k_axial``, the NMII internal stiffnesses AND the minifilament layout (``n_bb`` /
``n_side`` / ``backbone_len`` / ``head_offset``) are unresolved KB/PI GAPs (cards N1–N9, S1).  Every one of them
must be passed on the command line — none is defaulted.  The only two defaults are ``--k-off0``
(provisional-SOURCED) and ``--f0`` (the physical constant ``kBT/x_beta``).  This gate demonstrates the MECHANISM
and closes NO quantitative band.

CONVERGENCE IS REPORTED, NOT ASSUMED.  The α-actinin dorsal↔arc crosslinks (4.6e5 pN/µm, several meeting at one
arc apex) dominate the Gershgorin bound, so the CFL-stable explicit step is ~1e-7 µm/pN and the inner relax
needs THOUSANDS of iterations per outer step.  A tension read from an unconverged relax would be a transient,
not a force the fiber holds, so every step reports the free-node residual and its ratio to the tension, and the
verdict REQUIRES that ratio below 1%.

────────────────────────────────────────────────────────────────────────────────────────────────────────
RUN ON GBOOK (needs a CUDA GPU; will NOT run on the dev Mac):

    ssh gbook
    PYTHONPATH=~/ffn_ac_native:~/ffn_ac_native/ffn_sim \
      ~/miniconda3/envs/ffn_sim/bin/python \
      aleph/scripts/ac_gate_b_sf_motor_native.py \
        --k-axial 1000 --k-xb 1000 --k-backbone 1000 --k-head-arm 100 --backbone-lp 1.0 \
        --k-on 50 --f-stall 0.5 --v0 0.12 --kappa 0.5 --capture 0.210 \
        --n-bb 14 --n-side 10 --backbone-len 0.301 --head-offset 0.200 \
        --catch-slip --k-catch0 0.35 --x-catch 1.0e-3 --k-slip0 0.35 --x-slip 0.6e-3 \
        --steps 40 --dt 0.01 --relax 4000

``--relax 4000`` is the value at which the reported residual actually falls below 1% of the tension for the
default population; a smaller value prints a loud residual and FAILS the verdict rather than quietly reporting a
transient.  Long run: launch under nohup and monitor the LOG FILE (ssh python is not on PATH; use the full env
python).
────────────────────────────────────────────────────────────────────────────────────────────────────────

Runtime: NVIDIA Warp on CUDA only (I0-A).  Authored on the dev Mac (no CUDA) — the Lead runs it on gbook.
"""

from __future__ import annotations

import argparse
import json
import time

import numpy as np
import warp as wp

from aleph.engine.ledger import assemble_balance_tolerance_ratio, make_global_cell_ledger
from aleph.engine.observe.artifact import (
    observation_artifact,
    timing_block,
    write_artifact,
)
from aleph.engine.observe.energy import (
    accumulate_squared_displacement_kernel,
    balance_convergence,
    step_energy_balance,
)
from aleph.engine.sf_mechanics import build_sf_mechanics_topology
from aleph.engine.sf_motor_slice import build_sf_motor_slice, sf_sarcomere_geometry_for
from aleph.engine.sf_population import SFClass, build_sf_arc_population
from aleph.components.motor.backbone_warp import (
    arm_bending_cfl_stiffness,
    arm_orientation_k_theta,
    backbone_bending_cfl_stiffness,
    backbone_bending_k_theta,
    backbone_bending_kappa,
)
from aleph.components.motor.minifilament_topology import MinifilamentTopology
from aleph.engine.contracts import (
    EvidenceLabel,
    EvidenceRung,
    QuantitativeClaim,
    VoidCeiling,
    classify_gate,
)
from aleph.components.motor.segment_motor import SegmentDetachKinetics
from aleph.laws.network_warp import axpy_kernel

#: The residual/signal ceiling above which this gate's tension observable is not interpretable, DECLARED
#: HERE — at module scope, before any run — because a ceiling chosen after seeing the residual is
#: gate-loosening (PI decision D7, 2026-07-28).
#:
#: This gate is exactly the case ``GateVerdict.VOID`` was created for.  Its first landed run reported
#: ``MECHANISM PASS`` alongside a free-node residual at **15.5% of the tension it was reporting**, and an
#: earlier SF result — the withdrawn **1.92 pN** — was published as a PASS from precisely that state.  A
#: run in that condition has not measured a tension at all, so it must not be able to report PASS or FAIL
#: on a tension criterion; VOID is the honest third outcome.
#:
#: The ratio is the same 1% the driver's own convergence criterion already uses, so this introduces no
#: new constant: what D7 adds is that exceeding it now VOIDS the verdict instead of merely flipping one
#: boolean inside a dict of eight others.
SF_MOTOR_VOID_CEILING = VoidCeiling(
    0.01,
    "a free-node residual above 1% of the reported tension makes the tension a relaxation TRANSIENT "
    "rather than a force the fiber holds — measured on this very slice, a 10x longer relax moved T_max "
    "by 3.4x at an identical bound population, so neither PASS nor FAIL of a tension criterion is "
    "interpretable above this ratio",
)

#: Versioned gate contract this driver is measured against
#: (``aleph/docs/v2_audit/gate_contracts/GATE_B_sf_motor.yaml``). The contract pins the
#: thresholds, THE OBSERVABLE definitions, and the configuration; its content hash is
#: stamped into the run record so a later threshold/observable/configuration change cannot
#: turn a FAIL into a PASS unnoticed.
GATE_CONTRACT_ID = "GATE-B.sf_motor"


def _load_gate_contract_stamp(gate_id: str) -> dict:
    """Compute the gate-contract stamp for this run.

    ``aleph/outputs/tag_kb`` is not an importable package, so the reference
    implementation is loaded by path relative to the repo root (which also works on gbook,
    where the tree lives at ``~/ffn_ac_native``). HARD FAIL if it cannot be loaded: an
    unstamped gate artifact is precisely the hole this closes.
    """
    import importlib.util
    from pathlib import Path

    module_path = Path(__file__).resolve().parents[2] / "aleph" / "outputs" / "tag_kb" / "gate_contract.py"
    if not module_path.exists():
        raise SystemExit(f"GATE-CONTRACT MISSING: {module_path} not found — a gate may not run "
                         f"without a versioned contract to stamp its artifact with.")
    spec = importlib.util.spec_from_file_location("gate_contract", module_path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)  # type: ignore[union-attr]
    return module.stamp_from_gate_id(gate_id, provenance="RUNTIME")


#: Explicit-Euler (overdamped) Courant target for the inner relax — the same dimensionless CFL safety margin
#: the sf_arc native gate uses (hard limit C < 2; target 0.5 leaves a ≥2× margin that also absorbs the axial
#: geometric stiffening).  Stated and stability-derived, never fitted.
COURANT_TARGET = 0.5


@wp.kernel
def _pin_force_kernel(force: wp.array(dtype=wp.vec3d), pinned: wp.array(dtype=wp.int32)) -> None:
    """Zero the force on pinned (Dirichlet) nodes so the relax cannot move them (device-only)."""
    i = wp.tid()
    if pinned[i] != wp.int32(0):
        force[i] = wp.vec3d(wp.float64(0.0), wp.float64(0.0), wp.float64(0.0))


@wp.kernel
def _inject_unbalanced_force_kernel(force: wp.array(dtype=wp.vec3d), amount: wp.float64) -> None:
    """POSITIVE CONTROL: add a one-sided force to a single particle, breaking the Newton pair exactly.

    One thread (``dim=1``).  The balance gate asserts that the two bodies' force resultants cancel; this
    adds ``amount`` to ONE body and nothing to the other, so the resultant sum becomes exactly that
    vector and the gate must reject.  Without such a step in the same run, a gate that accepted every
    step is indistinguishable from a gate that cannot fail.

    The magnitude is the run's own ``f_stall`` — one head's worth of unbalanced force — so the control is
    scaled by the physics under test rather than by a chosen number.
    """
    force[0] = force[0] + wp.vec3d(amount, wp.float64(0.0), wp.float64(0.0))


class _TrajectoryRecorder:
    """Record ONE inner solve's trajectory: a subsampled polyline plus the EXACT ``Σ|Δx|²``.

    The two halves are measured differently on purpose.  A line integral survives subsampling to within
    its quadrature error, so the energy ledger's work terms may be taken on a coarse polyline; a sum of
    squared displacements does not (a chord is shorter than the path it replaces, so subsampling would
    understate the dissipation by an arbitrary factor).  The squared displacement is therefore
    accumulated on EVERY inner iteration by a device kernel and read once, out of the loop.

    Solver-agnostic: it reads the displacement the update actually produced rather than assuming
    ``Δx = μF``, so the explicit and implicit inner paths are measured by the same instrument.
    """

    def __init__(self, sf_owner: object, actuator: object, *, device: str, stride: int) -> None:
        self.sf_owner = sf_owner
        self.actuator = actuator
        self.device = device
        self.stride = max(1, int(stride))
        self.prev_sf_d = wp.zeros_like(sf_owner.position_d)
        self.prev_nmii_d = wp.zeros_like(actuator.position_d)
        self.squared_d = wp.zeros(1, dtype=wp.float64, device=device)
        self.samples: list[np.ndarray] = []
        self.enabled = False

    def _concatenated(self) -> np.ndarray:
        return np.concatenate(
            (np.asarray(self.sf_owner.position_d.numpy(), np.float64),
             np.asarray(self.actuator.position_d.numpy(), np.float64)), axis=0)

    def begin(self) -> None:
        """Start a fresh recording at the current configuration."""
        self.squared_d.zero_()
        self.samples = [self._concatenated()]

    def before_update(self) -> None:
        """D2D-copy both arrays so the next update's displacement can be differenced against them."""
        wp.copy(self.prev_sf_d, self.sf_owner.position_d)
        wp.copy(self.prev_nmii_d, self.actuator.position_d)

    def after_update(self, index: int) -> None:
        """Accumulate this iteration's ``|Δx|²`` on device, and sample the configuration on the stride."""
        for position, previous in ((self.sf_owner.position_d, self.prev_sf_d),
                                   (self.actuator.position_d, self.prev_nmii_d)):
            wp.launch(accumulate_squared_displacement_kernel, dim=int(position.shape[0]),
                      inputs=[position, previous, self.squared_d], device=self.device)
        if (index + 1) % self.stride == 0:
            self.samples.append(self._concatenated())

    def finish(self) -> tuple[np.ndarray, float]:
        """Close the recording at the final configuration and return ``(path, Σ|Δx|²)``."""
        final = self._concatenated()
        if not self.samples or not np.array_equal(self.samples[-1], final):
            self.samples.append(final)
        wp.synchronize_device(self.device)
        return np.asarray(self.samples, np.float64), float(self.squared_d.numpy()[0])


def _max_force_norm(force_d: wp.array, free_mask: np.ndarray | None = None) -> float:
    """OUT-OF-LOOP host diagnostic: max per-node force magnitude (device→host read between phases only).

    With ``free_mask`` the maximum is taken over the FREE (non-Dirichlet) nodes only — that is the inner solve's
    RESIDUAL.  A pinned node legitimately keeps a non-zero reaction after convergence (it is the traction), so
    including it would make the residual look permanently large and hide non-convergence.
    """
    f = force_d.numpy()
    if f.shape[0] == 0:
        return 0.0
    norms = np.linalg.norm(f, axis=1)
    if free_mask is not None:
        norms = norms[free_mask]
    return float(np.max(norms)) if norms.size else 0.0


def _sf_axial_tension(pos_d: wp.array, links: np.ndarray, link_k: np.ndarray, link_r0: np.ndarray) -> dict:
    """OUT-OF-LOOP: SF axial tension per segment ``T = k·(L−r0)`` [pN] (positive = stretched/under tension)."""
    pos = pos_d.numpy()
    a, b = links[:, 0].astype(np.int64), links[:, 1].astype(np.int64)
    length = np.linalg.norm(pos[b] - pos[a], axis=1)
    tension = link_k * (length - link_r0)
    return {
        "T_max_pN": float(np.max(tension)) if tension.size else 0.0,
        "T_mean_pN": float(np.mean(tension)) if tension.size else 0.0,
        "T_abs_max_pN": float(np.max(np.abs(tension))) if tension.size else 0.0,
    }


def _fa_traction(force_d: wp.array, pos_d: wp.array, fa_sites: np.ndarray, centroid: np.ndarray) -> dict:
    """OUT-OF-LOOP: the reaction each pinned FA anchor carries, and whether it points INWARD (contractile).

    A ventral SF is a closed dipole: both outer ends are FA-anchored, so a contracting fiber pulls each anchor
    TOWARD the fiber's interior.  The radial projection ``⟨F, r̂⟩`` with ``r̂`` pointing from the population
    centroid to the anchor is therefore NEGATIVE for contraction and positive for extension — the sign check
    that separates a real contractile dipole from a sign-flipped one.
    """
    if fa_sites.size == 0:
        return {"n_fa": 0, "F_total_pN": 0.0, "radial_mean_pN": 0.0, "inward_fraction": 0.0}
    force = force_d.numpy()[fa_sites]
    pos = pos_d.numpy()[fa_sites]
    radial = pos - centroid
    norm = np.linalg.norm(radial, axis=1, keepdims=True)
    radial = np.divide(radial, norm, out=np.zeros_like(radial), where=norm > 1e-12)
    projection = np.einsum("ij,ij->i", force, radial)
    return {
        "n_fa": int(fa_sites.size),
        "F_total_pN": float(np.sum(np.linalg.norm(force, axis=1))),
        "radial_mean_pN": float(np.mean(projection)),
        "inward_fraction": float(np.mean(projection < 0.0)),
    }


def _cfl_step(topology, *, k_xb: float, k_backbone: float, k_head_arm: float, backbone_lp: float,
              r0_backbone: float, r0_head: float, n_sf_nodes: int, n_nmii_particles: int) -> tuple[float, dict]:
    """Return the CFL-stable explicit ``dt/γ`` for BOTH relaxed arrays, from the stiffnesses actually built.

    Gershgorin row-sum bound on the material tangent of each array (the same construction the ``sf_arc`` native
    gate documents), then the step is ``C_target / λ_max`` over the WORSE of the two arrays:

      * SF nodes: ``2k`` per Hookean link/joint endpoint, ``4α``/``8α`` per bending end/centre node, plus
        ``2·k_xb`` for a bound crossbridge pulling on that node;
      * NMII particles: ``2·k_backbone`` (rod) + ``2·k_head_arm`` (lever) + the angle-harmonic contributions
        (via the landed ``*_bending_cfl_stiffness`` maps) + ``2·k_xb`` for the crossbridge.

    No magic number: every term is a stiffness this run actually assembled.
    """
    g_sf = np.zeros(int(n_sf_nodes), dtype=np.float64)
    for (i, j), k in zip(topology.links, topology.link_k):
        g_sf[int(i)] += 2.0 * float(k)
        g_sf[int(j)] += 2.0 * float(k)
    for (i, j), k in zip(topology.arc_joints, topology.arc_k):
        g_sf[int(i)] += 2.0 * float(k)
        g_sf[int(j)] += 2.0 * float(k)
    for (a, b, c), alpha in zip(topology.bend_triples, topology.bend_alpha):
        g_sf[int(a)] += 4.0 * float(alpha)
        g_sf[int(b)] += 8.0 * float(alpha)
        g_sf[int(c)] += 4.0 * float(alpha)
    g_sf += 2.0 * float(k_xb)                        # a bound head can pull on any SF node
    lambda_sf = float(g_sf.max()) if g_sf.size else 0.0

    # The angle-harmonic CFL contributions come from the landed derivation maps — the SAME ones the actuator's
    # NMIIBendingStiffness uses — so the relax step and the assembled bending stiffness cannot drift apart.
    k_theta_bb = backbone_bending_k_theta(backbone_bending_kappa(float(backbone_lp)), float(r0_backbone))
    k_theta_arm = arm_orientation_k_theta(float(k_xb), float(r0_head))
    lambda_nmii = (2.0 * float(k_backbone) + 2.0 * float(k_head_arm) + 2.0 * float(k_xb)
                   + float(backbone_bending_cfl_stiffness(k_theta_bb, float(r0_backbone)))
                   + float(arm_bending_cfl_stiffness(k_theta_arm, float(r0_head))))

    lambda_max = max(lambda_sf, lambda_nmii)
    if not np.isfinite(lambda_max) or lambda_max <= 0.0:
        raise SystemExit("[sf_motor GATE-B] FAIL: assembled stiffness is non-positive — nothing to relax")
    detail = {
        "lambda_sf_pN_per_um": lambda_sf,
        "lambda_nmii_pN_per_um": lambda_nmii,
        "lambda_max_pN_per_um": lambda_max,
        "courant_target": COURANT_TARGET,
        "n_sf_nodes": int(n_sf_nodes),
        "n_nmii_particles": int(n_nmii_particles),
    }
    return COURANT_TARGET / lambda_max, detail


def main() -> None:
    ap = argparse.ArgumentParser(
        description="GATE-B native: SF active tension emerges from nmii_sf_motor binding events.")
    # SF + NMII mechanical GAPs (all REQUIRED — no default anywhere).
    ap.add_argument("--k-axial", type=float, required=True, help="SF actin axial stiffness [pN/µm] (GAP)")
    ap.add_argument("--k-xb", type=float, required=True, help="crossbridge stiffness [pN/µm] (GAP)")
    ap.add_argument("--k-backbone", type=float, required=True, help="NMII backbone rod [pN/µm] (GAP)")
    ap.add_argument("--k-head-arm", type=float, required=True, help="NMII head-arm lever [pN/µm] (GAP)")
    ap.add_argument("--backbone-lp", type=float, required=True, help="NMII backbone L_p [µm] (GAP)")
    # head kinetics GAPs.
    ap.add_argument("--k-on", type=float, required=True, help="per-head attach rate [1/s] (GAP)")
    ap.add_argument("--f-stall", type=float, required=True, help="per-head stall [pN] (GAP)")
    ap.add_argument("--v0", type=float, required=True, help="unloaded head velocity [µm/s] (GAP)")
    ap.add_argument("--kappa", type=float, required=True, help="Hill curvature a/F0 [-] (dimensionless)")
    ap.add_argument("--capture", type=float, required=True, help="head→actin capture reach [µm]")
    # NOT GAPs, so these two keep documented defaults: k_off0 is provisional-SOURCED (Stam-Hocky/Tam) and f0 is
    # the physical constant kBT/x_beta (Veigel 2002 x_beta).  Everything else must be supplied.
    ap.add_argument("--k-off0", type=float, default=0.35, help="Bell slip prefactor [1/s] (provisional-sourced)")
    ap.add_argument("--f0", type=float, default=4.28e-3 / 0.6e-3, help="Bell f0 = kBT/x_beta [pN] (physical)")
    ap.add_argument("--catch-slip", action="store_true",
                    help="physiological Pereverzev catch-slip detach (PROXY constants) instead of Bell slip")
    # catch-slip constants are KB/PI GAPs: no default, and REQUIRED once --catch-slip selects that hazard.
    ap.add_argument("--k-catch0", type=float, default=None)
    ap.add_argument("--x-catch", type=float, default=None)
    ap.add_argument("--k-slip0", type=float, default=None)
    ap.add_argument("--x-slip", type=float, default=None)
    # minifilament reference topology — every one of these is a KB/PI GAP (I0-B3 layout), so none is defaulted.
    ap.add_argument("--n-bb", type=int, required=True, help="backbone beads per minifilament (GAP)")
    ap.add_argument("--n-side", type=int, required=True, help="heads per anti-parallel side (GAP)")
    ap.add_argument("--backbone-len", type=float, required=True, help="backbone contour [µm] (GAP)")
    ap.add_argument("--head-offset", type=float, required=True, help="head↔backbone arm offset [µm] (GAP)")
    # SF population.
    ap.add_argument("--n-ventral", type=int, default=8)
    ap.add_argument("--n-dorsal", type=int, default=4)
    ap.add_argument("--n-arc", type=int, default=4)
    ap.add_argument("--n-cap", type=int, default=4)
    ap.add_argument("--n-per-fiber", type=int, default=9)
    # run control.
    ap.add_argument("--steps", type=int, default=40, help="accepted physical steps")
    ap.add_argument("--dt", type=float, default=0.01, help="outer physical timestep [s]")
    ap.add_argument("--relax", type=int, default=200, help="inner overdamped relax iterations per step")
    # implicit inner solve (the fix for the non-convergence that BLOCKED the quantitative claim). The explicit
    # path stays the A/B control: same gate, same diagnostics, so the two can be compared directly.
    ap.add_argument("--implicit", action="store_true",
                    help="use the implicit preconditioned-CG inner solve (ac/engine/sf_implicit.py) instead of "
                         "the explicit overdamped relax; the explicit path remains available as the control")
    ap.add_argument("--cg-iterations", type=int, default=200,
                    help="PCG launch budget per Newton iteration (a compute budget, never physical time)")
    ap.add_argument("--newton", type=int, default=4,
                    help="implicit (Newton-like) iterations per outer step; the geometry is nonlinear so the "
                         "linear solve is repeated, with the KMC binding state frozen across them")
    ap.add_argument("--cg-check-every", type=int, default=0,
                    help="host early-exit cadence for the PCG loop (0 = full fixed budget). A control-flow "
                         "readback of the device convergence latch only; dx is bit-identical either way")
    ap.add_argument("--implicit-step-scale", type=float, default=1.0,
                    help="mobility (dt/γ) multiplier vs the explicit CFL bound. The IMEX step is "
                         "unconditionally stable so this may exceed 1; it is a CONVERGENCE-RATE knob and the "
                         "equilibrium must not depend on it (that invariance is the convergence evidence)")
    # acceptance predicate.  The balance gate is the DEFAULT (plan §10 step 2); the constant predicate is
    # kept only as the A/B control, because "the gate accepted every step" means nothing unless the same
    # run also shows a step the gate rejects.
    ap.add_argument("--constant-accept", action="store_true",
                    help="CONTROL: accept every step with a constant device predicate (the pre-2026-07-28 "
                         "behaviour) instead of letting the force-balance gate decide")
    ap.add_argument("--no-positive-control", action="store_true",
                    help="skip the balance-gate positive control (a deliberately unbalanced step that the "
                         "gate must reject); the control is what makes a passing gate non-vacuous")
    # energy ledger.
    ap.add_argument("--energy-every", type=int, default=1,
                    help="measure the step energy ledger every k accepted steps (0 disables it)")
    ap.add_argument("--energy-samples", type=int, default=16,
                    help="trajectory segments per measured step; the work integrals use this polyline "
                         "while the dissipation is accumulated on EVERY inner iteration")
    ap.add_argument("--no-energy-mobility-scaling", action="store_true",
                    help="skip the post-loop mobility-halving pass that says whether the balance residual "
                         "is the integrator's first-order consistency error or a genuine leak")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--device", type=str, default="cuda:0")
    ap.add_argument("--out", type=str, default="")
    ap.add_argument("--native-population", action="store_true",
                    help="ASSERT that this run's SF/NMII inventory is the physiological one derived from "
                         "density x geometry, not a slice. It is the ONLY thing that lets the record "
                         "stamp CONNECTED rather than CUDA_UNIT, and it is a flag rather than a node-count "
                         "threshold on purpose: 'native' is a statement about where the numbers came from, "
                         "which no count can check. DEVELOP ON A SLICE, CONCLUDE AT NATIVE (PI 2026-07-28)")
    ap.add_argument("--build-commit", type=str, default="",
                    help="the commit this run is measured on. REQUIRED on the native machine, whose tree is "
                         "synced rather than checked out: without it the build stamp degrades silently to "
                         "'unknown', which is the missing-build defect the record exists to close")
    args = ap.parse_args()

    if args.catch_slip:
        missing = [name for name in ("k_catch0", "x_catch", "k_slip0", "x_slip")
                   if getattr(args, name) is None]
        if missing:
            raise SystemExit(
                f"REQUIRED-PARAM: --catch-slip selects the Pereverzev catch-slip hazard, whose constants "
                f"{missing} are KB/PI GAPs with no audited NMII head-actin parameterisation. Pass them "
                f"explicitly (they are then reported as PROXY in the run record) or use the Bell-slip path."
            )

    # GATE CONTRACT — resolved BEFORE any GPU work, so a missing/renamed contract costs
    # milliseconds instead of a 46-minute native run. The hash covers the threshold, the
    # OBSERVABLE definition, and the configuration this gate is measured at; it is written
    # into the run record below so `verify_gate_contracts.py` can refuse a PASS whose
    # contract moved since the FAIL it supersedes (audit failure mode M2).
    contract_stamp = _load_gate_contract_stamp(GATE_CONTRACT_ID)
    print(f"[sf_motor GATE-B] gate contract {contract_stamp['gate_id']} "
          f"{contract_stamp['contract_hash']} ({contract_stamp['contract_path']})")

    run_started = time.perf_counter()
    wp.init()
    if not wp.get_device(args.device).is_cuda:
        raise SystemExit("GATE-B SF-motor native gate requires CUDA (I0-A); no CPU simulation path exists.")
    device = args.device
    from aleph.engine.cortex_motor_slice import NMIIMotorParams

    topology_nmii = MinifilamentTopology(
        n_bb=int(args.n_bb), n_heads_per_side=int(args.n_side),
        backbone_length_um=float(args.backbone_len), head_offset_um=float(args.head_offset))
    geometry = sf_sarcomere_geometry_for(topology_nmii)
    print(f"[sf_motor GATE-B] device={device}")
    print(f"[sf_motor GATE-B] DERIVED sarcomere geometry from the minifilament: "
          f"lateral={geometry['sarcomere_lateral_um']:.4g} µm (=2·head_offset)  "
          f"overlap={geometry['sarcomere_overlap_um']:.4g} µm (=backbone contour)")
    print(f"[sf_motor GATE-B] PI-GAP params (report-not-tune): k_axial={args.k_axial:g} k_xb={args.k_xb:g} "
          f"k_on={args.k_on:g} f_stall={args.f_stall:g} v0={args.v0:g} capture={args.capture:g} "
          f"detach={'CATCH_SLIP(proxy)' if args.catch_slip else 'SLIP'}")

    # 1) the DISJOINT sf_arc population, built with the straddle geometry the motor requires.
    population = build_sf_arc_population(
        n_ventral=args.n_ventral, n_dorsal=args.n_dorsal, n_arc=args.n_arc, n_cap=args.n_cap,
        n_per_fiber=args.n_per_fiber, **geometry)
    population.assert_partitioned()
    topology_sf = build_sf_mechanics_topology(population, k_axial_pn_per_um=float(args.k_axial))

    params = NMIIMotorParams(
        v0=float(args.v0), f_stall=float(args.f_stall), kappa=float(args.kappa), k_xb=float(args.k_xb),
        r0_head=float(args.head_offset), r0_xb=0.0, capture_radius=float(args.capture),
        k_on=float(args.k_on), k_off0=float(args.k_off0), f0=float(args.f0),
        # NOT `float(...)`: these four are `float | None` and default to None on purpose — they are
        # PI-GAPs that `--catch-slip` makes mandatory and the SLIP path never reads. `argparse` already
        # yields a float when one is supplied, so the cast added nothing and turned "not supplied" into
        # `TypeError: float() argument must be ... not 'NoneType'` — which made the DEFAULT slip path
        # unrunnable. The guard that belongs here already exists one layer down:
        # `CortexMotorParams.__post_init__` requires all four IFF `detach_kinetics is CATCH_SLIP`.
        k_catch0=args.k_catch0, x_catch=args.x_catch,
        k_slip0=args.k_slip0, x_slip=args.x_slip,
        detach_kinetics=(SegmentDetachKinetics.CATCH_SLIP if args.catch_slip else SegmentDetachKinetics.SLIP),
    )

    # 2) COMPOSE the slice: sf_arc owner + nmii owner + the nmii_sf_motor connector, ONE CellTransaction.
    #    The ledger comes with it, because the acceptance predicate is now the balance gate rather than a
    #    constant: each body reduces its OWN force array into one half of the gate, and the tolerance is
    #    derived on-device from those same resultants (PI D8) so nothing supplies a force magnitude.
    ledger = make_global_cell_ledger(n_topology_slots=1, device=device)
    tol_sq_d = wp.zeros(1, dtype=wp.float64, device=device)
    slice_ = build_sf_motor_slice(
        sf_population=population, minifilament_topology=topology_nmii, params=params,
        sf_k_axial_pn_per_um=float(args.k_axial),
        nmii_k_backbone_pn_per_um=float(args.k_backbone),
        nmii_k_head_spring_pn_per_um=float(args.k_head_arm),
        nmii_backbone_persistence_length_um=float(args.backbone_lp),
        base_seed=int(args.seed), device=device,
        ledger=ledger, tol_sq_d=tol_sq_d)
    sf_owner, actuator, connector = slice_.sf_owner, slice_.actuator_state, slice_.connector
    n_heads = int(actuator.n_heads)
    n_nmii_particles = int(actuator.position_d.shape[0])

    # D8's DEVICE side, which had no caller until here.  ``n_terms`` is the number of accumulated force
    # contributions across BOTH channels — one per node of each body — and it is a COUNT known at build
    # time, never a magnitude.  The float64 summation bound for that count is uploaded once; the tolerance
    # itself is derived per step, on the device, from the ledger's own |reaction| + |traction|.
    n_balance_terms = int(sf_owner.n_nodes) + n_nmii_particles
    gamma_n = assemble_balance_tolerance_ratio(n_balance_terms)
    gamma_n_d = wp.array(np.array([gamma_n], np.float64), dtype=wp.float64, device=device)
    slice_.balance_gamma_n_d = gamma_n_d
    print(f"[sf_motor GATE-B] BALANCE GATE: {n_balance_terms} accumulated contributions ⇒ float64 "
          f"summation ratio γ_n = {gamma_n:.4e} (D8; the tolerance is γ_n·(|reaction|+|traction|), "
          f"derived on-device each step from the ledger's own resultants)")
    station_census = population.motor_station_census()
    print(f"[sf_motor GATE-B] bound {len(slice_.participants)} participants: sf_arc(n_nodes={sf_owner.n_nodes}, "
          f"n_links={topology_sf.n_links}) + nmii(n_minifilaments={actuator.n_minifilaments}, "
          f"n_heads={n_heads}) + {connector.name}")
    print(f"[sf_motor GATE-B] motor stations: {station_census['N_stations']} on "
          f"{station_census['ready_by_class']}; EXCLUDED (curved sarcomere, cannot be straddled): "
          f"{station_census['excluded_by_class'] or 'none'}")
    print(f"[sf_motor GATE-B] SPLIT OWNERSHIP: sf pos ptr={sf_owner.position_d.ptr} vs nmii pos ptr="
          f"{actuator.position_d.ptr} (two never-merged arrays)")

    # 3) pin the FA anchors of the ventral fibers (a closed contractile dipole needs its ends held).
    fa_sites = np.asarray(population.fa_sites, np.int64)
    pinned_host = np.zeros(int(sf_owner.n_nodes), np.int32)
    pinned_host[fa_sites] = 1
    pinned_d = wp.array(pinned_host, dtype=wp.int32, device=device)
    free_mask = pinned_host == 0
    centroid = np.asarray(population.pos, np.float64).mean(axis=0)
    print(f"[sf_motor GATE-B] pinned {int(pinned_host.sum())} FA anchor nodes (Dirichlet) of "
          f"{sum(1 for b in population.bundles if b.sf_class == SFClass.VENTRAL)} ventral + dorsal fibers")

    # 4) the CFL-stable inner relax step, DERIVED from the assembled stiffnesses of both arrays.
    step_val, cfl_detail = _cfl_step(
        topology_sf, k_xb=float(args.k_xb), k_backbone=float(args.k_backbone),
        k_head_arm=float(args.k_head_arm), backbone_lp=float(args.backbone_lp),
        r0_backbone=float(topology_nmii.segment_length_um), r0_head=float(args.head_offset),
        n_sf_nodes=int(sf_owner.n_nodes), n_nmii_particles=int(actuator.position_d.shape[0]))
    step = wp.float64(step_val)
    print(f"[sf_motor GATE-B] CFL: λ_max ≤ {cfl_detail['lambda_max_pN_per_um']:.4g} pN/µm "
          f"(sf {cfl_detail['lambda_sf_pN_per_um']:.4g} / nmii {cfl_detail['lambda_nmii_pN_per_um']:.4g}) "
          f"⇒ dt/γ = {COURANT_TARGET:g}/λ_max = {step_val:.3e} µm/pN")

    # 4b) OPTIONAL implicit inner solve (--implicit).  The explicit path below is kept as the A/B CONTROL:
    #     the convergence claim is only credible if the two paths agree on the equilibrium tension, and if
    #     the implicit answer is invariant to its own iteration budget.  Nothing is removed.
    implicit_solver = None
    implicit_detail: dict = {}
    if args.implicit:
        from aleph.engine.sf_implicit import SFImplicitCG, SFImplicitTopology

        implicit_topology = SFImplicitTopology(
            device=device, n_sf=int(sf_owner.n_nodes), n_nmii=int(actuator.position_d.shape[0]),
            sf_topology=topology_sf, actuator=actuator, connector_state=connector.state,
            k_xb_pn_per_um=float(args.k_xb), r0_xb_um=float(params.r0_xb),
            pinned_sf_nodes=fa_sites)
        implicit_solver = SFImplicitCG(
            implicit_topology, max_iterations=int(args.cg_iterations),
            cg_check_every=int(args.cg_check_every))
        # The implicit step is the IMEX overdamped update (a = 1/(dt/γ)); it is unconditionally stable, so the
        # mobility may exceed the explicit CFL bound.  The scale is a CONVERGENCE-RATE knob only — the
        # equilibrium it converges to must not depend on it, which is exactly what the A/B check below tests.
        implicit_step = float(step_val) * float(args.implicit_step_scale)
        implicit_detail = {
            "cg_iterations": int(args.cg_iterations),
            "newton_iterations": int(args.newton),
            "mobility_step_um_per_pN": implicit_step,
            "mobility_step_scale_vs_cfl": float(args.implicit_step_scale),
            "n_pinned": int(implicit_topology.n_pinned),
            "regularizer_a_pN_per_um": 1.0 / implicit_step,
        }
        print(f"[sf_motor GATE-B] IMPLICIT PCG: {args.cg_iterations} CG iters x {args.newton} Newton, "
              f"mobility dt/γ = {implicit_step:.3e} µm/pN ({args.implicit_step_scale:g}x CFL) ⇒ "
              f"regularizer a = {1.0 / implicit_step:.4g} pN/µm; Dirichlet mask on "
              f"{implicit_topology.n_pinned} FA nodes (in the OPERATOR, not just the RHS)")

    # The energy ledger's dissipation channel uses the mobility the inner solve is ACTUALLY integrating
    # with, so the two paths report their own: the explicit CFL step, or the IMEX step when --implicit.
    if implicit_solver is not None:
        mobility_for_energy = float(implicit_step)
        energy_stride = 1                       # only `newton` iterations, so every one is a sample
    else:
        mobility_for_energy = float(step_val)
        energy_stride = max(1, int(args.relax) // max(1, int(args.energy_samples)))
    recorder = _TrajectoryRecorder(sf_owner, actuator, device=device, stride=energy_stride)
    # COUNTED, not computed from (steps x relax): the rollback probe, the positive control and the
    # mobility-scaling passes all run inner iterations too, and a formula written beside the loop drifts
    # from it. Cost per iteration is the term that dominates this lane, so it has to be the real count.
    inner_iterations = {"n": 0}

    def inner_solve_explicit(iterations: int | None = None, mobility: float | None = None) -> None:
        """Overdamped relax of BOTH owner arrays under the composed candidate force (pinned FA held).

        ``iterations``/``mobility`` default to the run's own; they are parameters only so the post-loop
        mobility-halving pass can walk the SAME relaxation at a smaller step, which is what separates an
        integrator-consistency residual from a genuine leak in the energy balance.
        """
        count = int(args.relax) if iterations is None else int(iterations)
        move = step if mobility is None else wp.float64(float(mobility))
        for index in range(count):
            inner_iterations["n"] += 1
            if recorder.enabled:
                recorder.before_update()
            slice_.zero_forces()
            slice_.accumulate()
            wp.launch(_pin_force_kernel, dim=int(sf_owner.n_nodes),
                      inputs=[sf_owner.force_d, pinned_d], device=device)
            wp.launch(axpy_kernel, dim=int(sf_owner.n_nodes),
                      inputs=[sf_owner.position_d, move, sf_owner.force_d], device=device)
            wp.launch(axpy_kernel, dim=int(actuator.position_d.shape[0]),
                      inputs=[actuator.position_d, move, actuator.force_d], device=device)
            if recorder.enabled:
                recorder.after_update(index)

    def inner_solve_implicit() -> None:
        """Newton-like implicit iterations: rebuild the force, solve (aI + MKM + a(I−M)) dx = M F, apply dx.

        The geometry is nonlinear (the pair tangent and the crossbridge walk direction depend on position), so
        the linear solve is repeated a few times per outer physical step.  The KMC binding state is FROZEN
        across these iterations — it commits only on the accepted step — which is what makes the tangent the
        exact derivative of the force the candidate is accumulating.
        """
        for index in range(int(args.newton)):
            inner_iterations["n"] += 1
            if recorder.enabled:
                recorder.before_update()
            slice_.zero_forces()
            slice_.accumulate()
            implicit_solver.step(
                sf_position_d=sf_owner.position_d, sf_force_d=sf_owner.force_d,
                nmii_position_d=actuator.position_d, nmii_force_d=actuator.force_d,
                mobility_step=implicit_step)
            if recorder.enabled:
                recorder.after_update(index)

    def inner_solve() -> None:
        """Run the selected inner solve, then leave the CONVERGED force in place for the diagnostics.

        The final pass deliberately OMITS the pin mask: the force left on a pinned FA node IS the Dirichlet
        reaction that anchor carries (the traction). Masking it here would zero the very quantity the traction
        diagnostic reads, and the solve above is already finished, so leaving the reaction in place moves nothing.
        """
        if implicit_solver is not None:
            inner_solve_implicit()
        else:
            inner_solve_explicit()
        slice_.zero_forces()
        slice_.accumulate()

    links = np.ascontiguousarray(topology_sf.links.reshape(-1, 2))
    link_k = np.ascontiguousarray(topology_sf.link_k)
    link_r0 = np.ascontiguousarray(topology_sf.link_r0)

    # ── the energy ledger's two force closures ────────────────────────────────────────────────────────
    # Both write positions into the two owner arrays and accumulate the composed candidate force; neither
    # fires a KMC event or commits anything, so the discrete state is frozen and the integrals are
    # properties of the FIELD.  The only difference between them is the binding SoA, which is what makes
    # their difference the active (crossbridge) channel and nothing else.
    n_sf_nodes = int(sf_owner.n_nodes)
    binding_arrays = connector.binding_field_arrays()

    def _force_at(positions: np.ndarray) -> np.ndarray:
        array = np.ascontiguousarray(np.asarray(positions, np.float64))
        sf_owner.position_d.assign(np.ascontiguousarray(array[:n_sf_nodes]))
        actuator.position_d.assign(np.ascontiguousarray(array[n_sf_nodes:]))
        slice_.zero_forces()
        slice_.accumulate()
        wp.synchronize_device(device)
        return np.concatenate(
            (np.asarray(sf_owner.force_d.numpy(), np.float64),
             np.asarray(actuator.force_d.numpy(), np.float64)), axis=0)

    def _measure_step_energy(
        path: np.ndarray, squared: float, mobility: float, pre_binding: list[np.ndarray],
    ) -> tuple[object, bool]:
        """Form the step's balance from the field the step ACTUALLY integrated, then restore everything.

        Two states have to be put back afterwards and both are checked rather than trusted: the binding
        SoA (reinstated to its pre-step values, because the accepted step's commit has already moved it
        and the trajectory was walked under the OLD binding) and the two position arrays (the path
        integral moves them along the polyline).  A restore that is verified is the same discipline the
        heads-detached control uses; an unverified one would leave the next step measuring a different
        configuration than the one it reported.
        """
        post_binding = [np.asarray(array.numpy()).copy() for array in binding_arrays]
        post_sf = np.asarray(sf_owner.position_d.numpy(), np.float64).copy()
        post_nmii = np.asarray(actuator.position_d.numpy(), np.float64).copy()
        # The FORCE arrays are part of what this observation destroys, and they are read afterwards as
        # the step's solver residual — so they are snapshotted and restored bit-exactly rather than
        # recomputed. Recomputing is subtly wrong and was measured to be: an `accumulate` after the step
        # sees the POST-commit binding (heads that just attached, walk references the Hill step just
        # moved), while every unmeasured step reports the PRE-commit force its solve left behind. That
        # difference is the KMC's own discontinuity, not a residual, and it made the measured steps look
        # 8x less converged than their neighbours.
        post_sf_force = np.asarray(sf_owner.force_d.numpy(), np.float64).copy()
        post_nmii_force = np.asarray(actuator.force_d.numpy(), np.float64).copy()
        for array, values in zip(binding_arrays, pre_binding, strict=True):
            array.assign(values)
        bound_pre = np.asarray(connector.state.bound_d.numpy()).copy()
        detached = np.zeros_like(bound_pre)

        def force_passive(positions: np.ndarray) -> np.ndarray:
            connector.state.bound_d.assign(detached)
            value = _force_at(positions)
            connector.state.bound_d.assign(bound_pre)
            return value

        ledger_record = step_energy_balance(
            path=path, force_total=_force_at, force_passive=force_passive,
            squared_displacement_um2=squared, mobility_um_per_pN=mobility,
            # event_jump is ZERO BY DERIVATION here: the only binding-dependent force in this slice is the
            # crossbridge, which lives entirely in the active channel, so an accepted KMC commit at fixed
            # configuration moves no passive potential.  It is passed explicitly rather than omitted.
            event_jump=0.0,
        )
        for array, values in zip(binding_arrays, post_binding, strict=True):
            array.assign(values)
        sf_owner.position_d.assign(post_sf)
        actuator.position_d.assign(post_nmii)
        sf_owner.force_d.assign(post_sf_force)
        actuator.force_d.assign(post_nmii_force)
        restored = (
            all(np.array_equal(np.asarray(array.numpy()), values)
                for array, values in zip(binding_arrays, post_binding, strict=True))
            and np.array_equal(np.asarray(sf_owner.position_d.numpy(), np.float64), post_sf)
            and np.array_equal(np.asarray(actuator.position_d.numpy(), np.float64), post_nmii)
            and np.array_equal(np.asarray(sf_owner.force_d.numpy(), np.float64), post_sf_force)
            and np.array_equal(np.asarray(actuator.force_d.numpy(), np.float64), post_nmii_force)
        )
        return ledger_record, restored

    # --- (1) REST: unbound heads, force-free rest length ⇒ composed force ≈ 0 ---------------------------
    slice_.zero_forces()
    slice_.accumulate()
    wp.synchronize_device(device)
    rest = {
        "n_bound": slice_.bound_head_count(),
        "sf_maxF_pN": _max_force_norm(sf_owner.force_d),
        "nmii_maxF_pN": _max_force_norm(actuator.force_d),
        **_sf_axial_tension(sf_owner.position_d, links, link_k, link_r0),
    }
    print(f"[sf_motor GATE-B] (1) REST  bound={rest['n_bound']}/{n_heads}  sf max|F|={rest['sf_maxF_pN']:.3e} "
          f"nmii max|F|={rest['nmii_maxF_pN']:.3e} pN  T_abs_max={rest['T_abs_max_pN']:.3e} pN  (≈0 expected)")
    if rest["n_bound"] != 0:
        raise SystemExit("[sf_motor GATE-B] FAIL: heads are pre-bound at t0 — the emergence claim would be void")

    # --- (2) EMERGENCE over accepted steps -------------------------------------------------------------
    # THE PREDICATE.  ``accepted_d=None`` hands acceptance to the device balance flag the ledger writes:
    # the two bodies' force resultants must cancel to within the tolerance derived on-device from those
    # same resultants.  The constant predicate is retained only as the A/B control (--constant-accept).
    constant_accept_d = (
        wp.array(np.ones(1, np.int32), dtype=wp.int32, device=device) if args.constant_accept else None)
    predicate_source = "constant_ones_CONTROL" if args.constant_accept else "ledger.balance_ok_d"
    print(f"[sf_motor GATE-B] acceptance predicate = {predicate_source}")
    energy_every = int(args.energy_every)
    trace: list[dict] = []
    energy_rows: list[dict] = []
    clock_before = int(slice_.clock.accepted_step_index_d.numpy()[0])
    for s in range(int(args.steps)):
        measure_energy = energy_every > 0 and (s + 1) % energy_every == 0
        pre_binding = ([np.asarray(array.numpy()).copy() for array in binding_arrays]
                       if measure_energy else [])
        recorder.enabled = measure_energy
        if measure_energy:
            recorder.begin()
        slice_.step(inner_solve, dt_phys=float(args.dt), accepted_d=constant_accept_d)
        wp.synchronize_device(device)
        recorder.enabled = False
        if measure_energy:
            path, squared = recorder.finish()
            energy, energy_restored = _measure_step_energy(
                path, squared, mobility_for_energy, pre_binding)
            energy_rows.append({
                "step": s + 1, "state_restored_exactly": bool(energy_restored),
                # A REJECTED step is rolled back, so its final configuration is not where the recorded
                # trajectory ended and the balance is not a balance of anything. Recorded rather than
                # filtered, so the artifact shows which rows are interpretable.
                "step_was_accepted": int(ledger.balance_ok_d.numpy()[0]) == 1,
                **energy.as_artifact_fields(),
            })
        tension = _sf_axial_tension(sf_owner.position_d, links, link_k, link_r0)
        traction = _fa_traction(sf_owner.force_d, sf_owner.position_d, fa_sites, centroid)
        loads = connector.state.loads_bell_d.numpy()
        # RESIDUAL of the inner solve on the FREE nodes: without this, a tension that is merely a
        # still-relaxing transient is indistinguishable from one the bound heads are actually holding.
        residual_sf = _max_force_norm(sf_owner.force_d, free_mask)
        residual_nmii = _max_force_norm(actuator.force_d)
        row = {
            "step": s + 1,
            "n_bound": slice_.bound_head_count(),
            **tension,
            "fa_radial_mean_pN": traction["radial_mean_pN"],
            "fa_inward_fraction": traction["inward_fraction"],
            "head_load_max_pN": float(np.max(loads)) if loads.size else 0.0,
            "residual_sf_free_pN": residual_sf,
            "residual_nmii_pN": residual_nmii,
            "residual_over_tension": (residual_sf / tension["T_abs_max_pN"]
                                      if tension["T_abs_max_pN"] > 0.0 else float("inf")),
            # the acceptance predicate as it was evaluated for THIS step (out-of-loop host read of the
            # flag the device wrote; the step itself never touched it on the host)
            "balance_ok": int(ledger.balance_ok_d.numpy()[0]),
            "balance_residual_pN": float(np.sqrt(float(ledger.balance_residual_sq_d.numpy()[0]))),
            "balance_tolerance_pN": float(np.sqrt(float(tol_sq_d.numpy()[0]))),
            "reaction_resultant_pN": float(np.linalg.norm(ledger.reaction_resultant_d.numpy()[0])),
            "traction_resultant_pN": float(np.linalg.norm(ledger.traction_resultant_d.numpy()[0])),
        }
        trace.append(row)
        if s % max(1, int(args.steps) // 10) == 0 or s == int(args.steps) - 1:
            print(f"[sf_motor GATE-B] step {row['step']:3d}  bound={row['n_bound']:4d}/{n_heads}  "
                  f"T_max={row['T_max_pN']:+.4g} pN  T_mean={row['T_mean_pN']:+.4g}  "
                  f"FA_radial={row['fa_radial_mean_pN']:+.4g} pN (inward {100*row['fa_inward_fraction']:.0f}%)  "
                  f"head_load_max={row['head_load_max_pN']:.4g} pN  "
                  f"resid_free={residual_sf:.3e} ({100*row['residual_over_tension']:.2f}% of T)  "
                  f"balance_ok={row['balance_ok']} "
                  f"(|R+T|={row['balance_residual_pN']:.3e} vs tol {row['balance_tolerance_pN']:.3e} pN)",
                  flush=True)
        if measure_energy:
            last = energy_rows[-1]
            print(f"[sf_motor GATE-B]   energy: ΔU={last['delta_potential_pN_um']:+.4e} "
                  f"D={last['dissipated_pN_um']:.4e} W_act={last['active_input_pN_um']:+.4e} "
                  f"resid={last['residual_pN_um']:+.3e} pN·µm (closure {last['closure_ratio']:.3e}); "
                  f"path-independence of ΔU {last['path_independence_ratio']:.3e}", flush=True)
    clock_after = int(slice_.clock.accepted_step_index_d.numpy()[0])

    # --- (3) ROLLBACK: a rejected step must change no binding state and no clock ------------------------
    rejected_d = wp.array(np.zeros(1, np.int32), dtype=wp.int32, device=device)
    bound_before = connector.state.bound_d.numpy().copy()
    seg_before = connector.state.seg_id_d.numpy().copy()
    abscissa_before = connector.state.abscissa_d.numpy().copy()
    sf_pos_before = sf_owner.position_d.numpy().copy()
    nmii_pos_before = actuator.position_d.numpy().copy()
    clock_pre_reject = int(slice_.clock.accepted_step_index_d.numpy()[0])
    slice_.step(inner_solve, dt_phys=float(args.dt), accepted_d=rejected_d)
    wp.synchronize_device(device)
    rollback_detail = {
        "binding_soa_restored": bool(np.array_equal(bound_before, connector.state.bound_d.numpy())
                                     and np.array_equal(seg_before, connector.state.seg_id_d.numpy())
                                     and np.array_equal(abscissa_before,
                                                        connector.state.abscissa_d.numpy())),
        # the slice docstring claims the SF/NMII POSITIONS are bit-restored on rejection too — check it, do not
        # take the claim on trust (the inner solve moves both arrays, so this is the load-bearing half).
        "sf_position_bit_restored": bool(np.array_equal(sf_pos_before, sf_owner.position_d.numpy())),
        "nmii_position_bit_restored": bool(np.array_equal(nmii_pos_before, actuator.position_d.numpy())),
        "clock_did_not_advance": int(slice_.clock.accepted_step_index_d.numpy()[0]) == clock_pre_reject,
    }
    rollback_ok = all(rollback_detail.values())

    # --- (4) POSITIVE CONTROL: the balance gate must be able to REJECT --------------------------------
    # A gate that accepted every step is indistinguishable, in an artifact, from a gate that cannot fail.
    # So break the Newton pair by exactly one head's worth of stall force on ONE body and require the gate
    # to reject it — with the same run, the same tolerance, and the same predicate path.  This is the
    # ng6-style positive control the residency gate already uses, applied to acceptance.
    positive_control: dict = {"ran": False}
    if not args.no_positive_control and not args.constant_accept:
        kick = float(args.f_stall)
        pos_sf_before = sf_owner.position_d.numpy().copy()
        pos_nmii_before = actuator.position_d.numpy().copy()
        bound_before_control = connector.state.bound_d.numpy().copy()
        clock_pre_control = int(slice_.clock.accepted_step_index_d.numpy()[0])

        def inner_solve_unbalanced() -> None:
            """The identical inner solve, then one one-sided force so the two channels cannot cancel."""
            inner_solve()
            wp.launch(_inject_unbalanced_force_kernel, dim=1,
                      inputs=[actuator.force_d, wp.float64(kick)], device=device)

        slice_.step(inner_solve_unbalanced, dt_phys=float(args.dt), accepted_d=None)
        wp.synchronize_device(device)
        control_residual = float(np.sqrt(float(ledger.balance_residual_sq_d.numpy()[0])))
        positive_control = {
            "ran": True,
            "injected_unbalanced_force_pN": kick,
            "injected_as": "one head's stall force on a single nmii particle — the run's own f_stall, so "
                           "the control is scaled by the physics under test rather than by a chosen number",
            "balance_ok": int(ledger.balance_ok_d.numpy()[0]),
            "balance_residual_pN": control_residual,
            "balance_tolerance_pN": float(np.sqrt(float(tol_sq_d.numpy()[0]))),
            "residual_matches_injection": bool(
                abs(control_residual - kick) <= 1.0e-6 * max(kick, 1.0e-30)),
            "gate_rejected": int(ledger.balance_ok_d.numpy()[0]) == 0,
            "sf_position_bit_restored": bool(
                np.array_equal(pos_sf_before, sf_owner.position_d.numpy())),
            "nmii_position_bit_restored": bool(
                np.array_equal(pos_nmii_before, actuator.position_d.numpy())),
            "binding_soa_restored": bool(
                np.array_equal(bound_before_control, connector.state.bound_d.numpy())),
            "clock_did_not_advance": (
                int(slice_.clock.accepted_step_index_d.numpy()[0]) == clock_pre_control),
        }
        print(f"[sf_motor GATE-B] (4) POSITIVE CONTROL: injected {kick:g} pN one-sided ⇒ "
              f"balance_ok={positive_control['balance_ok']} "
              f"(|R+T|={control_residual:.6g} pN, expected {kick:g}); "
              f"rolled back exactly: "
              f"{positive_control['sf_position_bit_restored'] and positive_control['nmii_position_bit_restored']}",
              flush=True)

    # --- (5) MOBILITY SCALING: is the balance residual the integrator, or a leak? ----------------------
    # Pure observation, outside the transaction: relax from the current configuration at μ and again at
    # μ/2 with twice the iterations (the SAME relaxation), restoring the state bit-exactly between the
    # two.  A first-order integrator's consistency error halves; a genuine leak does not move.
    mobility_scaling: dict = {"ran": False}
    if energy_every > 0 and not args.no_energy_mobility_scaling and implicit_solver is None:
        anchor_sf = sf_owner.position_d.numpy().copy()
        anchor_nmii = actuator.position_d.numpy().copy()
        anchor_binding = [np.asarray(array.numpy()).copy() for array in binding_arrays]

        def _relax_pass(iterations: int, mobility: float) -> tuple[object, bool]:
            recorder.stride = max(1, iterations // max(1, int(args.energy_samples)))
            recorder.enabled = True
            recorder.begin()
            inner_solve_explicit(iterations=iterations, mobility=mobility)
            recorder.enabled = False
            path, squared = recorder.finish()
            return _measure_step_energy(path, squared, mobility, anchor_binding)

        def _restore_anchor() -> bool:
            for array, values in zip(binding_arrays, anchor_binding, strict=True):
                array.assign(values)
            sf_owner.position_d.assign(anchor_sf)
            actuator.position_d.assign(anchor_nmii)
            return bool(np.array_equal(sf_owner.position_d.numpy(), anchor_sf)
                        and np.array_equal(actuator.position_d.numpy(), anchor_nmii))

        base_ledger, base_restored = _relax_pass(int(args.relax), mobility_for_energy)
        restored_a = _restore_anchor()
        fine_ledger, fine_restored = _relax_pass(2 * int(args.relax), 0.5 * mobility_for_energy)
        restored_b = _restore_anchor()
        convergence = balance_convergence(base=base_ledger, refined=fine_ledger)
        mobility_scaling = {
            "ran": True,
            "base": base_ledger.as_artifact_fields(),
            "halved_mobility": fine_ledger.as_artifact_fields(),
            "convergence": convergence.as_artifact_fields(),
            "state_restored_exactly": bool(
                base_restored and fine_restored and restored_a and restored_b),
            "protocol": "two observation-only relaxations from the SAME configuration — μ with N "
                        "iterations and μ/2 with 2N — with every position and binding array restored "
                        "bit-exactly between them and afterwards; no event fires and nothing commits",
        }
        print(f"[sf_motor GATE-B] (5) MOBILITY SCALING: residual "
              f"{base_ledger.balance.residual:+.4e} → {fine_ledger.balance.residual:+.4e} pN·µm at half "
              f"mobility ⇒ exponent {convergence.residual_exponent} "
              f"(1 = integrator consistency, 0 = leak) → {convergence.verdict}", flush=True)
        recorder.stride = energy_stride

    # --- verdict ---------------------------------------------------------------------------------------
    bound_series = np.array([row["n_bound"] for row in trace], dtype=np.int64)
    tension_series = np.array([row["T_abs_max_pN"] for row in trace], dtype=np.float64)
    load_series = np.array([row["head_load_max_pN"] for row in trace], dtype=np.float64)
    residual_ratio = np.array([row["residual_over_tension"] for row in trace], dtype=np.float64)
    final = trace[-1] if trace else {}
    ledgers = {led.component: {"block": list(led.block), "active": led.active_count}
               for led in slice_.transaction.population_ledgers}

    # Causality, not correlation.  Both the bound count and the tension rise monotonically with the step index,
    # so a positive correlation between them is structurally guaranteed and proves nothing.  What DOES separate
    # "the heads produced this tension" from "a still-relaxing transient" is:
    #   (a) while NO head is bound the tension must be zero (nothing else in the slice can make tension);
    #   (b) the tension must EXCEED the largest single-head load — a single crossbridge cannot account for it,
    #       so the load has to have been transmitted and accumulated along the fiber;
    #   (c) the inner solve's free-node residual must be small COMPARED TO the tension being reported, else the
    #       number is a transient of the relax rather than a force the fiber is holding.
    unbound_steps = bound_series == 0
    tension_zero_while_unbound = (bool(np.all(tension_series[unbound_steps] < 1.0e-9))
                                  if unbound_steps.any() else True)
    # MECHANISM claims — structural, sign, and causal.  None of them depends on the inner solve having reached
    # equilibrium: they say the motor edge exists, fires on events, transmits load in the contractile direction,
    # and obeys the accepted-step transaction.
    mechanism = {
        "heads_unbound_at_t0": rest["n_bound"] == 0,
        "rest_force_free": rest["sf_maxF_pN"] < 1.0e-6 and rest["T_abs_max_pN"] < 1.0e-6,
        "heads_bound_by_events": int(bound_series.max()) > 0 if bound_series.size else False,
        "tension_zero_while_no_head_is_bound": tension_zero_while_unbound,
        "tension_exceeds_single_head_load": (
            float(tension_series[-1]) > float(load_series[-1]) > 0.0 if tension_series.size else False),
        "fa_traction_is_inward": float(final.get("fa_inward_fraction", 0.0)) > 0.5,
        "clock_advanced_once_per_accepted_step": clock_after - clock_before == int(args.steps),
        "rejected_step_changed_nothing": rollback_ok,
        # ACCEPTANCE.  Every step of the trace was admitted by the device balance flag, not by a constant.
        "every_step_accepted_by_the_balance_gate": (
            all(int(row["balance_ok"]) == 1 for row in trace) if trace else False),
        # …and the gate can say no: the positive control breaks the Newton pair by one head's stall force
        # and the same predicate path rejects it and bit-restores the step.  Without this, "the gate
        # accepted every step" is not evidence of anything.
        "balance_gate_rejects_a_broken_adjoint_pair": (
            bool(positive_control.get("gate_rejected")
                 and positive_control.get("residual_matches_injection")
                 and positive_control.get("sf_position_bit_restored")
                 and positive_control.get("nmii_position_bit_restored")
                 and positive_control.get("clock_did_not_advance"))
            if positive_control.get("ran") else True),
        # ENERGY.  The ledger has to be measurable without disturbing the run it measures.
        "energy_ledger_restored_state_exactly": (
            all(bool(row["state_restored_exactly"]) for row in energy_rows) if energy_rows else True),
        # …and its residual is the integrator's first-order consistency error, not a leak.  This is the
        # only claim here that could catch a double count, so it is a mechanism criterion, not a note.
        "energy_balance_residual_is_not_a_leak": (
            mobility_scaling["convergence"]["verdict"] != "leak"
            if mobility_scaling.get("ran") else True),
    }
    if args.constant_accept:
        mechanism["every_step_accepted_by_the_balance_gate"] = True  # the control deliberately bypasses it
    # QUANTITATIVE claim — separate, and currently the binding constraint.  Any tension MAGNITUDE is only
    # meaningful once the inner solve has converged; below that the number is a relaxation transient (measured:
    # a 10x longer relax moved T_max by 3.4x while the bound population was identical).
    quantitative = {
        "inner_solve_converged_relative_to_tension": (
            float(residual_ratio[-1]) < 0.01 if residual_ratio.size else False),
    }
    verdict = {**mechanism, **quantitative}

    # D7 classification. The observable this gate reports is the SF axial tension, so the residual it is
    # judged against is the inner solve's free-node residual on the SAME configuration. VOID overrides
    # both PASS and FAIL: a run whose residual swamps its signal measured nothing, whatever its
    # criterion booleans say.
    tension_signal = float(tension_series[-1]) if tension_series.size else 0.0
    residual_signal = (
        float(residual_ratio[-1]) * tension_signal if residual_ratio.size else 0.0)
    gate_verdict = classify_gate(
        passed=all(verdict.values()), residual=abs(residual_signal), signal=abs(tension_signal),
        ceiling=SF_MOTOR_VOID_CEILING)

    report = {
        "gate": "GATE-B native: nmii_sf_motor — SF active tension emerges from myosin binding events",
        # the contract this run was measured under (threshold + observable + configuration)
        "gate_contract": contract_stamp,
        "device": device,
        "sarcomere_geometry_derived": geometry,
        "population": {
            "sf": population.census(),
            "motor_station_census": station_census,
            "n_sf_nodes": int(sf_owner.n_nodes), "n_sf_links": int(topology_sf.n_links),
            "n_minifilaments": int(actuator.n_minifilaments), "n_heads": n_heads,
        },
        "population_ledgers_disjoint": ledgers,
        "cfl": cfl_detail,
        "solver": {
            "mode": "implicit_pcg" if implicit_solver is not None else "explicit_overdamped_relax",
            "explicit_relax_iterations": int(args.relax),
            **implicit_detail,
        },
        "params_pi_gap": {
            "k_axial": args.k_axial, "k_xb": args.k_xb, "k_backbone": args.k_backbone,
            "k_head_arm": args.k_head_arm, "backbone_lp": args.backbone_lp, "k_on": args.k_on,
            "f_stall": args.f_stall, "v0": args.v0, "kappa": args.kappa, "capture": args.capture,
            "detach": "CATCH_SLIP_PROXY" if args.catch_slip else "SLIP",
        },
        "rest": rest,
        "trace": trace,
        "final": final,
        "clock": {"before": clock_before, "after": clock_after, "steps": int(args.steps)},
        "rollback_detail": rollback_detail,
        "acceptance": {
            "predicate": predicate_source,
            "balance_gate": {
                "channels": "sf_arc body force resultant (reaction) vs nmii body force resultant "
                            "(traction) — the two never-merged arrays of this cut, accumulated "
                            "independently, so their cancellation tests the crossbridge's adjoint "
                            "scatter and not the convergence",
                "n_accumulated_contributions": n_balance_terms,
                "gamma_n": float(gamma_n),
                "tolerance_derivation": "PI D8: γ_n = (n−1)u/(1−(n−1)u) for float64, scaled ON DEVICE by "
                                        "the ledger's own Σ|f_i|; the gate takes no supplied magnitude",
            },
            "positive_control": positive_control,
        },
        "energy_ledger": {
            "per_step": energy_rows,
            "mobility_scaling": mobility_scaling,
            "terms": "ΔU = −∫F_passive·dx along the step's own trajectory; W_active = the same walk with "
                     "the force evaluated with heads bound and heads detached, differenced; dissipated = "
                     "Σ|Δx|²/μ accumulated on device every inner iteration; event_jump = 0 BY DERIVATION "
                     "(the only binding-dependent force here is the crossbridge, which is the active "
                     "channel, so an accepted commit at fixed configuration moves no passive potential)",
            "no_energy_function": "no potential is implemented anywhere — a second implementation of the "
                                  "forces would be free to drift from the forces it claims to "
                                  "differentiate, which is the 'gate a broken build cannot fail' pattern",
        },
        "verdict": verdict,
        "verdict_mechanism": mechanism,
        "verdict_quantitative": quantitative,
        "gate_verdict": {
            "verdict": gate_verdict.value,
            "criterion_met": all(verdict.values()),
            "observable": "SF axial tension max|T| [pN] on the final accepted step",
            "residual_pN": abs(residual_signal),
            "signal_pN": abs(tension_signal),
            "residual_over_signal": (
                float(residual_ratio[-1]) if residual_ratio.size else None),
            "void_ceiling_ratio": SF_MOTOR_VOID_CEILING.ratio,
            "void_ceiling_rationale": SF_MOTOR_VOID_CEILING.rationale,
        },
        "quantitative_claim_status": (
            "CONVERGED — tension magnitudes in this record are equilibrium values"
            if all(quantitative.values()) else
            "BLOCKED — the explicit overdamped relax has NOT converged (free-node residual "
            f"{100.0 * float(residual_ratio[-1]) if residual_ratio.size else float('nan'):.1f}% of the reported "
            "tension). Every tension/traction MAGNITUDE here is a relaxation TRANSIENT and must not be quoted: "
            "measured, a 10x longer relax changes T_max by ~3.4x at an identical bound population. The MECHANISM "
            "claims above are unaffected (they are structural/sign/causal). Closing the quantitative claim needs "
            "an implicit/CG inner solve for this slice — the alpha-actinin dorsal-arc crosslink (4.6e5 pN/um, "
            "several meeting at one arc apex) dominates the Gershgorin bound, so the CFL-stable explicit step is "
            "~1e-7 um/pN while the SF axial mode relaxes orders of magnitude slower. This is the same wall the "
            "cortex hit at GATE A, where ProjectedAnalyticCG + the fiber-arclength multigrid resolved it."),
        "note": "MECHANISM demo: tension rises with the bound-head population from ZERO, with no lumped k_SF "
                "and no prestress seed. Magnitudes are PROVISIONAL — every kinetic/mechanical constant is a "
                "KB/PI GAP (cards N1-N9, S1). No quantitative band is closed here.",
    }
    # Wrap the gate report in the ONE self-stamping run record (D1-B two axes, D7 verdict, build stamp).
    # The rung is what this run MEASURED, never a hand-typed label: acceptance now runs on a physical
    # device predicate with a positive control, which is the CONNECTED criterion — but only when the run
    # is at native population.  A slice run is a mechanism check, so it stamps CUDA_UNIT and says why.
    # DEVELOP ON A SLICE, CONCLUDE AT NATIVE (PI 2026-07-28).
    native_population = bool(args.native_population)
    rung = EvidenceRung.CONNECTED if native_population else EvidenceRung.CUDA_UNIT
    rung_basis = (
        "the accepted-step predicate is the device force-balance flag (not a constant), its tolerance is "
        "derived on-device from the ledger's own accumulated magnitudes, a positive control shows the gate "
        "rejecting a deliberately broken adjoint pair and bit-restoring the step, and the energy balance is "
        "recorded per step with its residual identified by its own mobility scaling"
        + ("" if native_population else
           " — measured on a SLICE, so the rung stops at CUDA_UNIT: nothing here is at the native "
           "population, and this repo has already had a slice result fail to transfer (the fiber-quotient "
           "coarse operator)"))
    record = observation_artifact(
        run_label="GATE-B nmii_sf_motor — physical acceptance predicate + closed energy ledger",
        evidence=EvidenceLabel(
            rung=rung, quantitative=QuantitativeClaim.BLOCKED, basis=rung_basis),
        config={
            "argv": {key: value for key, value in sorted(vars(args).items())},
            "gate_contract": contract_stamp,
            "sarcomere_geometry_derived": geometry,
        },
        census={
            "n_sf_nodes": int(sf_owner.n_nodes),
            "n_nmii_particles": n_nmii_particles,
            "n_nodes_total": int(sf_owner.n_nodes) + n_nmii_particles,
            "n_dof": 3 * (int(sf_owner.n_nodes) + n_nmii_particles),
            "n_sf_links": int(topology_sf.n_links),
            "n_arc_joints": int(topology_sf.n_arc_joints),
            "n_minifilaments": int(actuator.n_minifilaments),
            "n_heads": n_heads,
            "is_native_population": native_population,
            "population_note": "the native cell is 511,114 nodes; a slice is for developing the "
                               "instrument, never for concluding with it",
            "sf": population.census(),
            "motor_station_census": station_census,
            "population_ledgers_disjoint": ledgers,
        },
        t0=rest,
        timing=timing_block(
            wall_seconds=time.perf_counter() - run_started,
            physical_time_s=float(args.dt) * int(args.steps),
            n_steps=int(args.steps),
            n_inner_iterations=int(inner_iterations["n"]),
            device_note=(
                "NOT a benchmark: the energy ledger instruments "
                f"{len(energy_rows)} of {args.steps} steps (2 extra device copies + 2 launches per inner "
                "iteration on those steps, plus host path integrals), and the mobility-scaling pass adds "
                "3x one step's relax. Build/compile time is included."),
        ),
        measurements=report,
        device=device,
        parameter_provenance={
            name: "PI_GAP" for name in (
                "k_axial", "k_xb", "k_backbone", "k_head_arm", "backbone_lp", "k_on", "f_stall",
                "v0", "kappa", "capture", "n_bb", "n_side", "backbone_len", "head_offset")
        } | {"k_off0": "SOURCED_PROVISIONAL", "f0": "DERIVED"},
        void_ceiling=SF_MOTOR_VOID_CEILING,
        gate_passed=all(verdict.values()),
        residual=abs(residual_signal),
        signal=abs(tension_signal),
        declared_commit=(args.build_commit or None),
    )
    text = json.dumps(record, indent=2, default=str)
    print(text, flush=True)
    if args.out:
        write_artifact(args.out, record)

    print(f"[sf_motor GATE-B] D7 gate verdict = {gate_verdict.value} "
          f"(residual/signal {residual_ratio[-1] if residual_ratio.size else float('nan'):.4g} vs "
          f"void ceiling {SF_MOTOR_VOID_CEILING.ratio:g})", flush=True)

    mechanism_failed = [name for name, ok in mechanism.items() if not ok]
    quantitative_failed = [name for name, ok in quantitative.items() if not ok]
    if mechanism_failed:
        raise SystemExit(f"[sf_motor GATE-B] MECHANISM FAILED: {mechanism_failed} (see report)")
    print("[sf_motor GATE-B] MECHANISM PASS: nmii_sf_motor is a real MOTOR runtime — heads bound by k_on "
          "EVENTS from zero, the tension is zero while no head is bound, it exceeds any single head's load "
          "(so load is transmitted along the fiber), the FA reaction is INWARD (contractile), one predicate "
          "reached every rollback/commit, and a rejected step restored the binding SoA and both position "
          "arrays bit-exactly.", flush=True)
    print(f"[sf_motor GATE-B] EVIDENCE rung = {rung.value} "
          f"({'native' if native_population else 'SLICE — develop here, conclude at native'})", flush=True)
    if quantitative_failed:
        raise SystemExit(
            "[sf_motor GATE-B] QUANTITATIVE CLAIM BLOCKED — the inner solve has not converged "
            f"(residual {100.0 * float(residual_ratio[-1]):.1f}% of the reported tension). The MECHANISM is "
            "established; the tension MAGNITUDES in this record are relaxation transients and must not be "
            "quoted. Fix = an implicit/CG inner solve for this slice, not a longer explicit relax "
            "(see quantitative_claim_status in the report).")
    print("[sf_motor GATE-B] QUANTITATIVE PASS: the inner solve converged, so the reported tension/traction are "
          "equilibrium magnitudes (still PROVISIONAL in absolute terms — the params are PI-GAPs).", flush=True)


if __name__ == "__main__":
    main()
