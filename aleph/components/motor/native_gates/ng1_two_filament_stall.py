r"""NG-1 (I3 motor) — two-filament sarcomere ISOMETRIC-STALL gate (run on the gbook A5000 by the lead).

The cleanest bug-vs-density DECISION on real hardware: ONE production-topology Stam-Hocky NMII minifilament
(``ac.cell.assemble`` NMII_* values: n_bb=14, N_side=10, L_bb=0.301 µm, offset=0.200 µm) wedged between TWO
straight ANTIPARALLEL actin filaments — outward barbed ends, the actin held rigid (every actin node clamped =
the strict isometric limit; a far-end clamp of a rigid track bears the same net reaction). The minifilament is
centred so all 20 heads overlap actin at zero crossbridge strain and can bind.

It runs the REAL device KMC + inner mechanical relaxation (NOT a closed form): each tick relaxes the
minifilament mechanics (backbone rod + head arms + power-stroke crossbridge, :meth:`MyosinForce.accumulate`)
with the actin pinned, refreshes the per-head Hill load (:meth:`MyosinForce.compute_loads`), then advances one
attach → step+Bell-detach KMC tick (:meth:`MyosinForce.step_kinetics`). The heads bind, walk (abscissa grows),
their crossbridge tension rises toward ``f_stall`` so the Hill velocity self-limits, and the Bell slip sets an
EMERGENT engaged fraction — nothing is imposed. At settled isometric stall it measures:

  NG-1.m1  the two clamp reactions  (net crossbridge force delivered to filament A and to filament B) — should
           be equal-and-opposite ≈ ±F_side, the emergent per-side stall from ``two_filament_reference``.
  NG-1.m2  the TRANSMISSION RATIO   |reaction_B| / |reaction_A|  (crossbridge → arm → backbone → opposite
           filament). ≈ 1 ⇒ the motor + series compliance are faithful; ≪ 1 ⇒ a real series-compliance /
           geometry BUG (the analytic reference proves the ideal collinear ratio is exactly 1).
  NG-1.m3  per-head Hill load       → ``f_stall``  (isometric: the walk stalled).
  NG-1.m4  emergent bound fraction  → φ_b = k_on/(k_on + k_off0 e^{f_stall/f0})  (Bell, not an imposed duty).

VERDICT (CLAUDE.md §6.2 decision): transmission ≪ 1 ⇒ fix the motor (density-independent). transmission ≈ 1
AND per-side ≈ analytic F_side ⇒ the motor is faithful and the LOW per-minifilament force is real — carry it to
the full-cell DENSITY question (``two_filament_reference.sweep_gap_claims``): the cortical tension is floored by
the sparse *measured* minifilament areal density (Nie 2015 ~0.625 µm⁻²), NOT a broken motor. NEVER add
heads/density to lift the reaction — that is the whole point of separating the two failure modes here.

Tolerances are physically motivated (N_side=10 ⇒ Binomial KMC shot noise on the engaged count), NOT tuned to
pass: a genuine NG-1 failure is a real on-device transmission/stall bug the CPU analytic reference cannot show.
Do NOT loosen a tolerance or edit a kernel to force a pass — surface to PI.

I0-A: launching the Warp kernels REQUIRES a CUDA device, so this runner raises on any non-CUDA device and
only executes on the gbook A5000.

Run (from the repo root on gbook):
    PYTHONPATH=. ~/miniconda3/envs/ffn_sim/bin/python -m aleph.components.motor.native_gates.ng1_two_filament_stall
"""

from __future__ import annotations

import argparse
import sys

import numpy as np
import warp as wp

# production topology + mechanical/kinetic magnitudes — the EXACT ac.cell.assemble t0 values (no drift).
from aleph.components.incumbent.assemble import (
    NMII_BACKBONE_LP_DIAGNOSTIC_UM,
    NMII_CAPTURE_UM,
    NMII_F0,
    NMII_F_STALL_TEST,
    NMII_HEAD_OFFSET_UM,
    NMII_K_BACKBONE_TEST,
    NMII_K_HEAD_ARM_TEST,
    NMII_K_XB_TEST,
    NMII_KAPPA_TEST,
    NMII_KOFF0,
    NMII_KON_TEST,
    NMII_L_BB_UM,
    NMII_N_BB,
    NMII_N_SIDE,
    NMII_V0_TEST,
)
from aleph.components.motor import two_filament_reference as ref
from aleph.components.motor.backbone_warp import (
    arm_orientation_k_theta,
    backbone_bending_k_theta,
    backbone_bending_kappa,
)
from aleph.components.motor.hand import NMIIHandParams
from aleph.components.motor.minifilament_topology import MinifilamentTopology
from aleph.components.motor.minifilament_warp import (
    MyosinForce,
    build_minifilament_nodes,
    crossbridge_kernel,
)

DEVICE: str | None = None   # resolve the current CUDA device; the hardware contract forbids fixed ordinals
SEED = 20260717

# actin track discretization — a discrete NODE clamp track (NOT the production cortex actin, which is a 0.5 µm
# segmented filament bound by a point-to-segment query; this nearest-NODE model is the NG-1 stand-in and a
# known simplification — see the module note). NOT physics magnitudes: numerical placement of the immovable
# clamp. Post-F6 the head-arm is RIGID, so a head that binds a node offset in x carries that offset as a real
# crossbridge strain the arm can no longer swing away (the pre-fix floppy arm silently absorbed it).
# ⚠ HONEST NOTE (Codex audit): the worst-case binding offset spacing/2 = 2 nm is ~4× the stall stretch
# f_stall/k_xb = 0.5 nm — NOT "well below" it. A freshly-bound head therefore starts OVER-strained and relaxes
# onto its node (abscissa = 0) before it walks; finer spacing shrinks that transient but does not remove it.
# The production point-to-segment binding avoids the node-offset artifact entirely.
ACTIN_SPACING_UM = 0.004           # ≪ 2·NMII_CAPTURE_UM (0.210); worst-case head→node offset 0.002 µm (= 4× stall stretch)
ACTIN_HALF_SPAN_UM = 0.45          # covers the head x-range (±L_bb/2 ≈ 0.15 µm) with margin

# PASS tolerances (physically motivated; documented, NOT tuned to a band)
TRANSMISSION_LO, TRANSMISSION_HI = 0.80, 1.25   # |reaction_B|/|reaction_A|; ≪1 is a series/geometry bug
REACTION_RTOL = 0.35               # per-side reaction vs analytic F_side; Binomial(10, φ_b) shot noise ⇒ ~±3 heads
LOAD_LO_FRAC, LOAD_HI_FRAC = 0.50, 1.15         # mean per-head load / f_stall (REPORTED; GAP-sensitive — see sweep)
# ⭐ STALL criterion (Codex audit + the L_p sweep): the physically-correct isometric-stall test is that the walk
# has STOPPED, i.e. the Hill velocity ratio v/v0 → 0 — NOT that the absolute load hit a band (which the sweep
# showed is GAP-sensitive). v/v0 = (1 − f/f_s)/(1 + (f/f_s)/κ) is 0 exactly at stall; a settled v/v0 ≤ 0.15 means
# the head walks at ≤15% of its unloaded speed ⇒ essentially stalled, robustly across the backbone/arm GAP.
V_RATIO_STALL = 0.15
BOUND_FRAC_LO = 0.70               # emergent φ_b ≈ 0.99; allow KMC / mid-cycle spread


@wp.kernel
def _relax_free_nodes_kernel(
    pos: wp.array(dtype=wp.vec3d),
    force: wp.array(dtype=wp.vec3d),
    dt_mu: wp.float64,
    n_pinned: wp.int32,
):
    """Overdamped inner step ``x += dt_mu·F`` for the FREE (minifilament) nodes only; actin nodes are pinned.

    Nodes ``[0, n_pinned)`` are the clamped actin tracks (never moved); ``[n_pinned, n_total)`` are the
    minifilament backbone + heads that relax under the composed force.
    """
    t = wp.tid()
    if t < n_pinned:
        return
    pos[t] = pos[t] + dt_mu * force[t]


@wp.kernel
def _nearest_actin_kernel(
    pos: wp.array(dtype=wp.vec3d),
    head_node: wp.array(dtype=wp.int32),
    n_actin: wp.int32,
    out_dist: wp.array(dtype=wp.float64),
    out_idx: wp.array(dtype=wp.int32),
):
    """Device-resident nearest-actin neighbor query — one thread per head, brute-force over the clamped
    actin block ``pos[0:n_actin)`` (I0-A: the neighbor query stays on the GPU; NO per-tick D2H/H2D roundtrip).

    The actin nodes ``[0, n_actin)`` never move (``_relax_free_nodes_kernel`` pins them), so scanning the live
    ``pos`` is identical to scanning a fixed snapshot. For the tiny two-filament system (n_heads·n_actin ≈
    20·452) this dense scan is trivial; the production cortex uses a device point-to-segment query instead.
    """
    h = wp.tid()
    hp = pos[head_node[h]]
    best = wp.float64(1.0e30)
    bi = wp.int32(-1)
    for a in range(n_actin):
        d = wp.length(pos[a] - hp)
        if d < best:
            best = d
            bi = a
    out_dist[h] = best
    out_idx[h] = bi


def _build_two_filaments(offset_um: float, spacing_um: float, half_span_um: float):
    """Two straight ANTIPARALLEL actin tracks flanking the minifilament along the backbone (x) axis.

    Filament A at ``y = +offset`` (barbed end at −x), filament B at ``y = −offset`` (barbed end at +x): the +
    heads (walk_dir = −x) bind A, the − heads (walk_dir = +x) bind B, so the two barbed ends point OUTWARD
    (antiparallel). Returns ``(posA, posB, barbedA_hat, barbedB_hat)`` — the two node blocks + unit barbed
    directions for reporting.
    """
    xs = np.arange(-half_span_um, half_span_um + 0.5 * spacing_um, spacing_um)
    posA = np.column_stack([xs, np.full_like(xs, +offset_um), np.zeros_like(xs)])
    posB = np.column_stack([xs, np.full_like(xs, -offset_um), np.zeros_like(xs)])
    barbedA = np.array([-1.0, 0.0, 0.0])   # A's barbed end is at min-x (the + heads walk −x toward it)
    barbedB = np.array([+1.0, 0.0, 0.0])   # B's barbed end is at max-x
    return posA, posB, barbedA, barbedB



def _build_params() -> NMIIHandParams:
    """The production NMII hand params (ac.cell.assemble t0 values; the GAP magnitudes are provisional)."""
    p = NMIIHandParams()
    p.k_on = wp.float64(NMII_KON_TEST)
    p.k_off0 = wp.float64(NMII_KOFF0)
    p.f0 = wp.float64(NMII_F0)
    p.v0 = wp.float64(NMII_V0_TEST)
    p.f_stall = wp.float64(NMII_F_STALL_TEST)
    p.kappa = wp.float64(NMII_KAPPA_TEST)
    p.k_xb = wp.float64(NMII_K_XB_TEST)
    p.r0_head = wp.float64(NMII_HEAD_OFFSET_UM)
    p.r0_xb = wp.float64(0.0)
    p.capture_radius = wp.float64(NMII_CAPTURE_UM)
    return p


def run(ticks: int = 600, tau: float = 1.0e-3, inner_iters: int = 3000, avg_last: int = 80,
        device: str | None = DEVICE, lp_um: float | None = None, arm_mult: float = 1.0) -> dict:
    """Build the two-filament sarcomere, run the isometric-stall KMC, and return the measured decision bundle."""
    wp.init()
    # I0-A hard guard — never launch on a non-CUDA device (dev Mac). Check the device inventory BEFORE
    # resolving `device` (default None → the current CUDA device; no fixed ordinal per the hardware contract),
    # so a machine with no CUDA gives this clear message rather than an opaque device-identifier error.
    if not any(d.is_cuda for d in wp.get_devices()):
        raise RuntimeError(
            "ng1_two_filament_stall requires a CUDA GPU (I0-A); no CUDA device is present. "
            "This runner is authored source — run it on the gbook A5000.")
    dev = wp.get_device(device)
    if not dev.is_cuda:
        raise RuntimeError(
            f"ng1_two_filament_stall requires a CUDA GPU (I0-A). Resolved device {str(dev)!r} is not CUDA.")

    topo = MinifilamentTopology(n_bb=NMII_N_BB, n_heads_per_side=NMII_N_SIDE,
                                backbone_length_um=NMII_L_BB_UM, head_offset_um=NMII_HEAD_OFFSET_UM)

    # ── geometry: two antiparallel actin tracks + one centred minifilament (heads overlap actin) ──────────
    posA, posB, barbedA, barbedB = _build_two_filaments(NMII_HEAD_OFFSET_UM, ACTIN_SPACING_UM,
                                                        ACTIN_HALF_SPAN_UM)
    nA, nB = posA.shape[0], posB.shape[0]
    n_actin = nA + nB
    idxA = np.arange(0, nA)
    idxB = np.arange(nA, nA + nB)

    d = build_minifilament_nodes(topo, np.zeros(3), np.array([1.0, 0.0, 0.0]))
    mini_pos = d["positions"]
    backbone_bonds = (d["backbone_bonds"] + n_actin).astype(np.int32)
    head_bonds = (d["head_bonds"] + n_actin).astype(np.int32)
    head_node = (d["head_node"] + n_actin).astype(np.int32)
    backbone_angles = (d["backbone_angles"] + n_actin).astype(np.int32)   # F6 backbone bending triples
    head_arm_angles = (d["head_arm_angles"] + n_actin).astype(np.int32)   # F6 head-arm orientation triples
    walk_dir = d["walk_dir"]                 # +heads → −x (A barbed), −heads → +x (B barbed): the actin polarity
    n_hands = head_node.shape[0]

    pos_np = np.ascontiguousarray(np.concatenate([posA, posB, mini_pos], axis=0), np.float64)
    n_total = pos_np.shape[0]

    pos = wp.array(pos_np, dtype=wp.vec3d, device=device)
    bb_d = wp.array(backbone_bonds, dtype=wp.int32, device=device)
    hb_d = wp.array(head_bonds, dtype=wp.int32, device=device)
    hn_d = wp.array(head_node, dtype=wp.int32, device=device)
    ba_d = wp.array(backbone_angles, dtype=wp.int32, device=device) if backbone_angles.shape[0] else None
    ha_d = wp.array(head_arm_angles, dtype=wp.int32, device=device)
    n_seg = max(NMII_N_BB - 1, 1)
    r0_backbone = NMII_L_BB_UM / n_seg
    # Backbone bond stiffness = the params_i0b3 PER-SEGMENT value (applied per bond by harmonic_bond_kernel).
    # The ×(n_bb−1) "rigid-rod / end-to-end" reinterpretation is WITHDRAWN (Codex cross-audit + PI: the contract
    # is per-bond). The chain is therefore SOFT (end-to-end = k/(n_bb−1)); the analytic reference is corrected to
    # match (k_eff uses (n_bb−1)/k_backbone). ⚠ grid-invariance of a fixed per-bond k is an OPEN PI/KB contract.
    k_backbone_seg = NMII_K_BACKBONE_TEST
    # F6 fix: derived (Magic-Number-Block) bending stiffnesses — rigid rod + rigid head-arm so the crossbridge
    # tension transmits collinearly (transmission → 1). NOT tuned to a band (backbone_warp).
    # lp_um / arm_mult are CONVERGENCE-SWEEP diagnostics (default = the sourced values): they vary the GAP
    # k_θ to PROVE the 6/6 PASS is insensitive to the L_p GAP (transmission is symmetry-driven, load is
    # Hill-driven — the RESULT is k_θ-independent; k_θ only sets the convergence rate, and the rigid-rod
    # backbone (13000) dominates cfl so dt is ~unchanged across the sweep). NEVER a knob to make a gate pass.
    lp = float(lp_um) if lp_um is not None else NMII_BACKBONE_LP_DIAGNOSTIC_UM
    k_theta_bb = backbone_bending_k_theta(backbone_bending_kappa(lp), r0_backbone)
    k_theta_arm = arm_orientation_k_theta(NMII_K_XB_TEST, NMII_HEAD_OFFSET_UM) * float(arm_mult)
    myo = MyosinForce(bb_d, hb_d, hn_d, _build_params(),
                      k_backbone=k_backbone_seg, r0_backbone=r0_backbone,
                      k_head_spring=NMII_K_HEAD_ARM_TEST, n_hands=n_hands, walk_dir=walk_dir, device=device,
                      backbone_angles=ba_d, head_arm_angles=ha_d,
                      k_theta_backbone=k_theta_bb, k_theta_arm=k_theta_arm, segment_len_um=r0_backbone)

    f_d = wp.zeros(n_total, dtype=wp.vec3d, device=device)
    xb_only = wp.zeros(n_total, dtype=wp.vec3d, device=device)
    # CFL: fold the F6 bending + rigid-rod effective linear stiffness into kmax (MyosinForce.cfl_stiffness). A
    # backbone bead carries several stiff bonds (2 backbone + head arms + 2 backbone-angle + head-arm angles),
    # so its true diagonal stiffness exceeds the max single-bond kmax — use a 0.01 safety factor (not 0.1) so
    # the overdamped inner relaxation is stable (the soft-backbone 0.02 was a limit-cycle), compensated by more
    # inner_iters. The load then builds to f_stall over the PHYSICAL series-compliance chain stretch
    # (F_side/k_eff, matching the analytic) — hence the large tick budget (Hill self-limit is asymptotic).
    kmax = myo.cfl_stiffness
    dt_mu = 0.01 / kmax
    # Device-resident neighbor-query buffers (I0-A: the nearest-actin assignment is computed on the GPU each
    # tick by _nearest_actin_kernel — NO pos.numpy()/H2D roundtrip in the KMC hot loop).
    nd = wp.zeros(n_hands, dtype=wp.float64, device=device)
    ni = wp.full(n_hands, -1, dtype=wp.int32, device=device)

    def _xb_force() -> np.ndarray:
        """Fresh crossbridge-only nodal force [pN] (actin feels ONLY the crossbridge; backbone/arm are internal)."""
        xb_only.zero_()
        wp.launch(crossbridge_kernel, dim=n_hands,
                  inputs=[pos, hn_d, myo.state["bound"], myo.state["anchor"], myo.state["abscissa"],
                          myo.state["walk_dir"], myo.params.k_xb, myo.params.r0_xb],
                  outputs=[xb_only], device=device)
        wp.synchronize_device(device)
        return xb_only.numpy()

    # ── isometric-stall KMC loop coupled to the inner mechanical relaxation ───────────────────────────────
    reac_A_hist, reac_B_hist, load_hist, bound_hist, vratio_hist = [], [], [], [], []
    fstall_v = float(NMII_F_STALL_TEST)
    kappa_v = float(NMII_KAPPA_TEST)
    for tk in range(ticks):
        # nearest-actin neighbor query ON DEVICE (no D2H): fills nd/ni for the KMC attach in-place.
        wp.launch(_nearest_actin_kernel, dim=n_hands,
                  inputs=[pos, hn_d, wp.int32(n_actin)], outputs=[nd, ni], device=device)
        with wp.ScopedDevice(dev):
            for _ in range(inner_iters):
                f_d.zero_()
                myo.accumulate(pos, f_d)                    # backbone rod + head arms + power-stroke crossbridge
                wp.launch(_relax_free_nodes_kernel, dim=n_total,
                          inputs=[pos, f_d, wp.float64(dt_mu), wp.int32(n_actin)], device=dev)
            wp.synchronize_device(dev)
        myo.compute_loads(pos)                              # per-head tangential crossbridge tension (Hill/Bell)
        myo.step_kinetics(nd, ni, tau, SEED + tk)           # attach → step(Hill) → detach(Bell), device RNG

        if tk >= ticks - avg_last:                          # sample the settled isometric window
            fxb = _xb_force()
            reac_A_hist.append(float(fxb[idxA, 0].sum()))   # net crossbridge x-force delivered to filament A
            reac_B_hist.append(float(fxb[idxB, 0].sum()))   # ... to filament B (equal & opposite ⇒ ratio 1)
            loads = myo.loads.numpy()
            bound = myo.state["bound"].numpy().astype(bool)
            load_hist.append(float(loads[bound].mean()) if bound.any() else 0.0)
            bound_hist.append(float(bound.mean()))
            if bound.any():                                 # per-head Hill v/v0 = (1−fr)/(1+fr/κ), clamped ≥ 0
                fr = loads[bound] / fstall_v
                vr = np.maximum(0.0, (1.0 - fr) / (1.0 + fr / kappa_v))
                vratio_hist.append(float(vr.mean()))
            else:
                vratio_hist.append(1.0)

    reac_A = float(np.mean(reac_A_hist)) if reac_A_hist else 0.0
    reac_B = float(np.mean(reac_B_hist)) if reac_B_hist else 0.0
    mean_load = float(np.mean(load_hist)) if load_hist else 0.0
    bound_frac = float(np.mean(bound_hist)) if bound_hist else 0.0
    v_over_v0 = float(np.mean(vratio_hist)) if vratio_hist else 1.0
    transmission = abs(reac_B) / abs(reac_A) if reac_A != 0.0 else 0.0

    analytic = ref.two_filament_isometric(NMII_N_SIDE, NMII_F_STALL_TEST)
    return {
        "device": str(dev), "n_actin": n_actin, "n_heads": n_hands, "ticks": ticks, "tau": tau,
        "reaction_A_pN": reac_A, "reaction_B_pN": reac_B, "reaction_mag_pN": 0.5 * (abs(reac_A) + abs(reac_B)),
        "transmission_ratio": transmission, "mean_head_load_pN": mean_load, "bound_fraction": bound_frac,
        "v_over_v0": v_over_v0,
        "F_side_analytic_pN": analytic["f_side"], "phi_b_analytic": analytic["phi_b"],
        "f_stall_pN": NMII_F_STALL_TEST, "k_eff_pN_per_um": analytic["k_eff_pN_per_um"],
        "ideal_transmission": analytic["ideal_transmission"],
        "barbed_A": barbedA.tolist(), "barbed_B": barbedB.tolist(),
        "lp_um": lp, "lp_status": "DIAGNOSTIC_UNSOURCED_SWEEP", "arm_mult": float(arm_mult),
        "k_theta_bb": float(k_theta_bb),
        "k_theta_arm": float(k_theta_arm), "cfl_stiffness": float(myo.cfl_stiffness),
    }


def _verdict(r: dict) -> bool:
    """Physically-motivated PASS/FAIL on the measured decision bundle (tolerances documented, not tuned)."""
    fside = r["F_side_analytic_pN"]
    fstall = r["f_stall_pN"]
    checks = {
        "transmission≈1": TRANSMISSION_LO <= r["transmission_ratio"] <= TRANSMISSION_HI,
        "reaction_A≈F_side": abs(abs(r["reaction_A_pN"]) - fside) <= REACTION_RTOL * fside,
        "reaction_B≈F_side": abs(abs(r["reaction_B_pN"]) - fside) <= REACTION_RTOL * fside,
        "reactions_opposite": (r["reaction_A_pN"] * r["reaction_B_pN"]) < 0.0,
        # STALL = the walk has stopped (v/v0 → 0), the GAP-robust criterion — NOT the GAP-sensitive load band.
        "v/v0→0 (stalled)": r["v_over_v0"] <= V_RATIO_STALL,
        "φ_b_emergent": r["bound_fraction"] >= BOUND_FRAC_LO,
    }
    print(f"NG-1 two-filament sarcomere isometric-stall gate | device={r['device']} | "
          f"warp={wp.config.version} | actin_nodes={r['n_actin']} heads={r['n_heads']} "
          f"ticks={r['ticks']} seed={SEED}")
    print(f"  production topology: n_bb={NMII_N_BB} N_side={NMII_N_SIDE} L_bb={NMII_L_BB_UM}µm "
          f"offset={NMII_HEAD_OFFSET_UM}µm | k_xb={NMII_K_XB_TEST} k_arm={NMII_K_HEAD_ARM_TEST} "
          f"k_bb={NMII_K_BACKBONE_TEST} pN/µm | f_stall={fstall}pN")
    print(f"  antiparallel barbed ends: A={r['barbed_A']}  B={r['barbed_B']}  (outward)")
    print(f"  MEASURED clamp reactions:  A={r['reaction_A_pN']:+.4f} pN   B={r['reaction_B_pN']:+.4f} pN   "
          f"(mag {r['reaction_mag_pN']:.4f} pN)")
    print(f"  ANALYTIC  F_side = N_side·φ_b·f_stall = {fside:.4f} pN  (φ_b={r['phi_b_analytic']:.4f})")
    print(f"  TRANSMISSION ratio |B|/|A| = {r['transmission_ratio']:.4f}   "
          f"(ideal collinear = {r['ideal_transmission']:.1f}; k_eff = {r['k_eff_pN_per_um']:.2f} pN/µm)")
    print(f"  STALL: v/v0 = {r['v_over_v0']:.4f} (≤ {V_RATIO_STALL} ⇒ stalled)   "
          f"per-head load = {r['mean_head_load_pN']:.4f} pN (reported; GAP-sensitive)   "
          f"bound fraction = {r['bound_fraction']:.4f}")
    for name, ok in checks.items():
        print(f"    [{'PASS' if ok else 'FAIL'}] {name}")
    passed = all(checks.values())
    if passed and TRANSMISSION_LO <= r["transmission_ratio"] <= TRANSMISSION_HI:
        print("  DECISION: transmission ≈ 1 and per-side ≈ analytic ⇒ motor faithful → carry the low "
              "per-minifilament force to the full-cell DENSITY question (sweep_gap_claims / Nie-2015 floor).")
    elif not (TRANSMISSION_LO <= r["transmission_ratio"] <= TRANSMISSION_HI):
        print("  DECISION: transmission ≪ 1 ⇒ a motor / series-compliance BUG, INDEPENDENT of density → "
              "fix the motor before any density claim. Surface to PI.")
    print(f"NG-1 OVERALL: {'PASS' if passed else 'FAIL'} ({sum(checks.values())}/{len(checks)} checks)")
    return passed


def main() -> int:
    p = argparse.ArgumentParser(description="NG-1 two-filament sarcomere isometric-stall decision gate (A5000).")
    p.add_argument("--ticks", type=int, default=600)
    p.add_argument("--tau", type=float, default=1.0e-3, help="KMC tick [s]")
    p.add_argument("--inner-iters", type=int, default=3000,
                   help="inner overdamped relaxation iters/tick (raised for the stiffer F6 bending + stable dt)")
    p.add_argument("--avg-last", type=int, default=80, help="settled-window ticks to average the reactions over")
    p.add_argument("--device", type=str, default=DEVICE)
    p.add_argument("--lp-um", type=float, default=None,
                   help="backbone persistence-length GAP override [µm] (convergence-sweep diagnostic; "
                        "default = an explicitly non-production diagnostic sweep fixture)")
    p.add_argument("--arm-mult", type=float, default=1.0,
                   help="head-arm k_theta multiplier (convergence-sweep diagnostic; default 1.0)")
    a = p.parse_args()
    r = run(ticks=a.ticks, tau=a.tau, inner_iters=a.inner_iters, avg_last=a.avg_last, device=a.device,
            lp_um=a.lp_um, arm_mult=a.arm_mult)
    return 0 if _verdict(r) else 1


if __name__ == "__main__":
    sys.exit(main())
