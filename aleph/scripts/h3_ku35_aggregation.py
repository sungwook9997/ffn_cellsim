#!/usr/bin/env python
"""KU-3.5-ACTIVE force-AGGREGATION estimator (2026-06-04, authoritative gate).

The KU-3.5 gate splits (PI 2026-06-04, ``docs/KU_3_5_GATE_STRUCTURE_2026-06-04.md``)
into ACTIVE (``g_soft`` = myosin), PASSIVE (``g_rigid`` = structural), and a
Layer-2 BRIDGE anchor. The ACTIVE gate is still OPEN and the blocker is force
GENERATION / AGGREGATION: per-head myosin force IS generated (sub-stall,
delivered) but does NOT aggregate into a sustained shell tension. Connectivity,
percolation, turnover and coherence were all falsified (2026-06-03); the wall is
upstream — *why does the sum of per-head forces not become hoop tension?*

This module is the quantitative answer. For a live (or fixture) cortex snapshot
it computes a FORCE-AGGREGATION LEDGER that decomposes the chain

    per-head generation  →  coherent ceiling  →  realized shell tension (g_soft)

into orthogonal, individually-bounded efficiencies so the loss can be localised
to ONE branch:

  (a) GENERATION   Σ|F_head| over the bound myosin head-actin attach bonds
      (= exactly the force HOOMD applies: F = k_head_actin·(r − r0), r0 ≈ 0 in
      grip_walk). Per-head distribution + Σ|F| vs the band-required cut force.
  (b) HILL VALIDITY F/F_stall distribution. The ACTIVE gate is only valid at
      F/F_stall ≤ 1 (Hill-valid / sub-stall); super-stall rows are rejected.
      Reports BOTH the bond-frame load k_ha·(r−r0) and — when the live
      ``MyosinStepUpdater`` is supplied — the SERIES load k_series·min(s,r)
      that the Hill kernel actually stalls on (myosin.py:1280-1283).
  (c) REALIZED g_soft via the method-of-planes (the gate readout), decomposed
      per bond channel (attach / myosin-internal / xlink).
  (d) AGGREGATION EFFICIENCY = realized / expected-from-Σ|F_head|, split into
      a radial-projection loss (η_agg: are the pulls tangential, or do they go
      into local in/out bead motion?) and a medium-cancellation loss (η_medium:
      do the soft xlinks k=1e-7 / internal springs compress and cancel the
      attach tension?).

BRANCH DIAGNOSIS the ledger settles:
  * γ_ceiling < band  → GENERATION-LIMITED (per-head × N too small even if
    perfectly aggregated). Lever = recruitment / per-head force.
  * η_agg ≪ 1         → AGGREGATION-LIMITED (radial): myosin pulls are radial,
    force funnels into local bead displacement not hoop tension.
  * η_medium ≪ 1      → AGGREGATION-LIMITED (medium): the soft crosslink gauze /
    minifilament springs absorb/oppose (k=1e-7 dipole-direction cancellation).

CALIBRATION. The Irving-Kirkwood whole-shell hoop tension γ_IK = Σ T·L·sin²ψ /
(8πR²) is the rigorous "shell tension these bond forces produce"; the project
already validated (``h3_ku35_estimator_audit.py``) that the corrected
method-of-planes ≈ γ_IK to <0.5% on a synthetic isotropic shell, so the two
estimators are on the same scale and the efficiency ratios are meaningful. This
module reuses that exact 8πR² convention.

This is a MEASUREMENT module — it changes no physics. It is fixture-validated
(``tests/test_ku35_aggregation.py`` + ``--self-test`` here): known forces →
known aggregation verdict (tangential / radial / cancellation).

Run:
    conda run -n ffn_sim python -m aleph.scripts.h3_ku35_aggregation --self-test
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

PKG = Path(__file__).resolve().parents[1]

# KU-3.5 active band (N/m) and single-head stall (N), repo-pinned constants.
BAND = (0.35e-3, 0.65e-3)          # N/m, KU-3.5 cortical-tension band
F_STALL_BASE = 0.5e-12             # N, single-head NMIIA stall (config literal)

# Channel labels for the soft (harmonic) bond set. The myosin head→actin ATTACH
# bonds are the ONLY force injected into the cortex actin (the generators); the
# minifilament-internal springs and the crosslink gauze form the transmission
# medium that can cancel.
CH_ATTACH = "attach"               # cortex_myosin_attach_b* (generators)
CH_MYO_INTERNAL = "myosin_internal"  # cortex_myosin_backbone + _head_backbone
CH_XLINK = "xlink"                 # xlink_intra + xlink_attach_b*
CH_OTHER = "other"

# Verdict-label cutoffs (REPORTING heuristics, not physics gates — the raw
# numbers are always printed so the verdict is transparent on a continuum).
_ETA_AGG_LOW = 0.30                # below → radial-dominated aggregation loss
_ETA_MEDIUM_LOW = 0.50             # below → medium-cancellation dominates


def classify_channel(type_name: str) -> str:
    """Map a HOOMD bond type name to its KU-3.5 aggregation channel."""
    if type_name.startswith("cortex_myosin_attach"):
        return CH_ATTACH
    if type_name in ("cortex_myosin_backbone", "cortex_myosin_head_backbone"):
        return CH_MYO_INTERNAL
    if type_name.startswith("xlink"):
        return CH_XLINK
    return CH_OTHER


# --------------------------------------------------------------------------- #
# Pure per-bond geometry / tension (no HOOMD; fixture-testable)
# --------------------------------------------------------------------------- #
def bond_geometry(pos: np.ndarray, bg: np.ndarray, center: np.ndarray | None = None):
    """Per-bond direction, length, midpoint-radial and tangential fraction.

    Args:
        pos: ``(N, 3)`` positions (a coherent frame; NOT re-centered here).
        bg: ``(M, 2)`` bond endpoint indices into ``pos``.
        center: ``(3,)`` shell centre for the radial direction; defaults to the
            centroid of ``pos`` (the cell COM ≈ origin by construction).

    Returns:
        Tuple ``(u, L, sin2)`` where ``u`` is the ``(M, 3)`` unit bond vector,
        ``L`` the ``(M,)`` bond length and ``sin2`` the ``(M,)`` tangential
        fraction ``1 − (r̂·u)²`` (1 = bond purely tangential to the shell, 0 =
        purely radial).
    """
    if center is None:
        center = pos.mean(axis=0)
    rA = pos[bg[:, 0]]
    rB = pos[bg[:, 1]]
    d = rB - rA
    L = np.linalg.norm(d, axis=1)
    L_safe = np.where(L > 0, L, 1.0)
    u = d / L_safe[:, None]
    rmid = 0.5 * (rA + rB) - center
    rn = np.linalg.norm(rmid, axis=1)
    rn_safe = np.where(rn > 0, rn, 1.0)
    rhat = rmid / rn_safe[:, None]
    cos2 = np.sum(rhat * u, axis=1) ** 2
    sin2 = np.clip(1.0 - cos2, 0.0, 1.0)
    return u, L, sin2


def mop_tension(pos: np.ndarray, bg: np.ndarray, T: np.ndarray, R: float,
                center: np.ndarray | None = None, n_planes: int = 12) -> float:
    """Corrected method-of-planes cortical tension (matches the live gate).

    ``γ = mean_planes |Σ_{crossing} T·|u·n̂|| / (2πR)``. The ``|u·n̂|`` (abs on
    the projection, signed T retained) is the 2026-06-02 estimator fix that
    removes the spurious √N cancellation; identical algorithm to
    ``h3_ku35_tension._tension_method_of_planes``. Returns 0.0 for an empty
    bond set.
    """
    if bg.shape[0] == 0:
        return 0.0
    if center is None:
        center = pos.mean(axis=0)
    posc = pos - center
    rA = posc[bg[:, 0]]
    rB = posc[bg[:, 1]]
    d = rB - rA
    L = np.linalg.norm(d, axis=1)
    L_safe = np.where(L > 0, L, 1.0)
    u = d / L_safe[:, None]
    phi = (1 + 5 ** 0.5) / 2
    ii = np.arange(n_planes, dtype=np.float64)
    z = 1 - 2 * (ii + 0.5) / n_planes
    rxy = np.sqrt(np.clip(1 - z * z, 0.0, None))
    th = 2 * np.pi * ii / phi
    normals = np.stack([rxy * np.cos(th), rxy * np.sin(th), z], axis=1)
    gammas = []
    for nh in normals:
        a_side = rA @ nh
        b_side = rB @ nh
        cross = (a_side * b_side) < 0
        if not cross.any():
            gammas.append(0.0)
            continue
        f_cut = T[cross] * np.abs(u[cross] @ nh)
        gammas.append(float(np.sum(f_cut)) / (2.0 * np.pi * R))
    return float(np.mean(np.abs(np.asarray(gammas))))


def ik_hoop(L: np.ndarray, sin2: np.ndarray, T: np.ndarray, R: float) -> float:
    """Irving-Kirkwood whole-shell hoop tension: Σ T·L·sin²ψ / (8πR²) [N/m].

    Signed (T keeps its tension/compression sign). This is the rigorous
    mechanical surface tension the bond set produces; the project's estimator
    audit validated it ≈ corrected MOP on a synthetic isotropic shell.
    """
    if L.size == 0:
        return 0.0
    return float(np.sum(T * L * sin2) / (8.0 * np.pi * R ** 2))


def ik_ceiling(L: np.ndarray, T: np.ndarray, R: float) -> float:
    """Coherent tangential ceiling: Σ |T|·L / (8πR²) [N/m].

    The shell tension this bond set WOULD produce if every bond were purely
    tangential (sin²ψ → 1) and tensile (|T|, no sign cancellation) — i.e. the
    upper bound on aggregation for the given generated force magnitudes.
    """
    if L.size == 0:
        return 0.0
    return float(np.sum(np.abs(T) * L) / (8.0 * np.pi * R ** 2))


# --------------------------------------------------------------------------- #
# Core ledger (operates on extracted arrays; fixture-testable)
# --------------------------------------------------------------------------- #
def aggregation_ledger_from_bonds(
    *,
    pos: np.ndarray,
    bg: np.ndarray,
    channel: np.ndarray,
    T: np.ndarray,
    R_cell: float,
    F_stall: float,
    F_series: np.ndarray | None = None,
    g_soft_gate: float | None = None,
    center: np.ndarray | None = None,
    n_planes: int = 12,
    band: tuple[float, float] = BAND,
) -> dict:
    """Compute the KU-3.5-active force-aggregation ledger from bond arrays.

    Args:
        pos: ``(N, 3)`` positions (coherent frame; the cell COM is used as the
            shell centre unless ``center`` is given).
        bg: ``(M, 2)`` bond endpoint indices.
        channel: ``(M,)`` channel label per bond (see ``classify_channel``).
        T: ``(M,)`` signed scalar bond tension (+ tension, − compression).
        R_cell: shell radius [m] (cut circumference 2πR, shell area 4πR²).
        F_stall: single-head stall force [N] (mesoscale-scaled if applicable).
        F_series: optional ``(n_attach,)`` series load k_series·min(s,r) the
            Hill kernel stalls on, aligned to the attach bonds in ``bg`` order;
            when given the Hill-validity stats use it instead of the bond load.
        g_soft_gate: optional authoritative g_soft (mN/m → pass in N/m) from the
            live gate; reported as the official readout alongside the internal
            cross-check MOP.
        center: optional shell centre; defaults to ``pos.mean(0)``.
        n_planes: method-of-planes cut count.
        band: KU-3.5 active band (N/m).

    Returns:
        A nested dict: ``generation``, ``hill``, ``tension``, ``efficiency``,
        ``verdict``.
    """
    if center is None:
        center = pos.mean(axis=0)
    u, L, sin2 = bond_geometry(pos, bg, center=center)

    is_attach = channel == CH_ATTACH
    is_myo = is_attach | (channel == CH_MYO_INTERNAL)
    n_attach = int(is_attach.sum())

    # ---- (a) GENERATION (attach channel = force injected into actin) -------
    T_attach = T[is_attach]
    F_head = np.abs(T_attach)                 # |k·(r−r0)| per bound head
    sumF = float(F_head.sum())
    F_band = band[0] * 2.0 * np.pi * R_cell   # cut force required for the band
    gen = dict(
        n_engaged=n_attach,
        sumF_head_N=sumF,
        mean_F_head_N=float(F_head.mean()) if n_attach else 0.0,
        median_F_head_N=float(np.median(F_head)) if n_attach else 0.0,
        max_F_head_N=float(F_head.max()) if n_attach else 0.0,
        p90_F_head_N=float(np.percentile(F_head, 90)) if n_attach else 0.0,
        F_band_required_N=float(F_band),
        sumF_over_band=float(sumF / F_band) if F_band > 0 else float("nan"),
    )

    # ---- (b) HILL VALIDITY (F/F_stall) -------------------------------------
    F_for_stall = F_series if F_series is not None else F_head
    if F_for_stall is not None and F_for_stall.size:
        ratio = F_for_stall / F_stall
        hill = dict(
            load_source="series" if F_series is not None else "bond",
            mean_F_over_Fstall=float(ratio.mean()),
            median_F_over_Fstall=float(np.median(ratio)),
            max_F_over_Fstall=float(ratio.max()),
            p90_F_over_Fstall=float(np.percentile(ratio, 90)),
            frac_hill_valid=float(np.mean(ratio <= 1.0)),     # ≤ 1 = sub-stall
            frac_super_stall=float(np.mean(ratio > 1.0)),     # > 1 = artefact
            F_stall_N=float(F_stall),
        )
    else:
        hill = dict(load_source="none", frac_hill_valid=float("nan"),
                    F_stall_N=float(F_stall))

    # ---- (c) REALIZED g_soft (method-of-planes, per channel) ---------------
    def _mop(mask):
        return mop_tension(pos, bg[mask], T[mask], R_cell, center=center,
                           n_planes=n_planes) if mask.any() else 0.0
    is_soft = np.ones(bg.shape[0], dtype=bool)   # every bond passed in is soft
    g_mop_total = _mop(is_soft)
    g_mop_attach = _mop(is_attach)
    g_mop_myo = _mop(is_myo)

    # IK virial decomposition (additive, per channel).
    g_ik_total = ik_hoop(L, sin2, T, R_cell)
    g_ik_attach = ik_hoop(L[is_attach], sin2[is_attach], T_attach, R_cell)
    g_ceiling_attach = ik_ceiling(L[is_attach], T_attach, R_cell)
    tension = dict(
        g_soft_gate_Npm=g_soft_gate,                # authoritative (if supplied)
        g_mop_total_Npm=g_mop_total,                # cross-check (≈ gate)
        g_mop_attach_Npm=g_mop_attach,
        g_mop_myosin_Npm=g_mop_myo,
        g_ik_total_Npm=g_ik_total,
        g_ik_attach_Npm=g_ik_attach,
        g_ceiling_attach_Npm=g_ceiling_attach,
        mean_tangential_fraction=float(sin2[is_attach].mean()) if n_attach else 0.0,
    )

    # ---- (d) AGGREGATION EFFICIENCY ----------------------------------------
    # η_agg    : generator geometry (attach realized / attach ceiling) ∈ [0,1].
    # η_medium : medium cancellation (full soft IK / attach IK).
    # η_total  : overall = realized soft tension / coherent ceiling.
    def _safe_div(a, b):
        return float(a / b) if b not in (0.0, None) and np.isfinite(b) else float("nan")
    eta_agg = _safe_div(g_ik_attach, g_ceiling_attach)
    eta_medium = _safe_div(g_ik_total, g_ik_attach)
    eta_total = _safe_div(g_ik_total, g_ceiling_attach)
    ceiling_over_band = _safe_div(g_ceiling_attach, band[0])
    efficiency = dict(
        eta_agg=eta_agg,                 # radial-projection (tangential) loss
        eta_medium=eta_medium,           # xlink/spring cancellation
        eta_total=eta_total,
        ceiling_over_band=ceiling_over_band,
        gate_over_ceiling=_safe_div(
            g_soft_gate if g_soft_gate is not None else g_mop_total,
            g_ceiling_attach),
    )

    # ---- VERDICT (transparent labelling of the continuum) ------------------
    if n_attach == 0:
        verdict = "NO-MOTORS (passive baseline; n_engaged=0)"
    elif np.isfinite(ceiling_over_band) and ceiling_over_band < 1.0:
        verdict = (f"GENERATION-LIMITED: coherent ceiling {ceiling_over_band:.2g}× "
                   "the band — per-head force × N too small even if perfectly "
                   "aggregated (lever = recruitment / per-head force)")
    elif np.isfinite(eta_agg) and eta_agg < _ETA_AGG_LOW:
        verdict = (f"AGGREGATION-LIMITED (RADIAL): η_agg={eta_agg:.2g} — myosin "
                   "pulls are radial, force funnels into local bead in/out "
                   "motion not hoop tension")
    elif np.isfinite(eta_medium) and eta_medium < _ETA_MEDIUM_LOW:
        verdict = (f"AGGREGATION-LIMITED (MEDIUM): η_medium={eta_medium:.2g} — "
                   "soft xlink gauze / minifilament springs cancel the attach "
                   "tension (compression offsetting tension)")
    else:
        verdict = ("AGGREGATION OK on these axes (ceiling ≥ band, η_agg & "
                   "η_medium ~ 1) — realized tension follows generation")

    return dict(generation=gen, hill=hill, tension=tension,
                efficiency=efficiency, verdict=verdict)


# --------------------------------------------------------------------------- #
# Live-sim adapter
# --------------------------------------------------------------------------- #
def _series_loads_from_action(action, pos_byTag: np.ndarray,
                              k_head_actin: float, k_head_spring: float):
    """Read the Hill SERIES load k_series·min(s,r) per BOUND head from the updater.

    Mirrors ``MyosinStepUpdater`` Step-3 (myosin.py:1280-1283): the head-actin
    bond and the head-backbone spring bear the load in series, so the motor
    stalls on ``k_series·min(s, r)`` not the full bond tension. Returns
    ``(F_series, bead_tags)`` aligned to the bound heads (or ``(None, None)`` if
    the action has no grip state, e.g. binned_r0 / motors off).
    """
    bound = getattr(action, "_head_bound_to_actin", None)
    grip = getattr(action, "_head_grip_s", None)
    if bound is None or grip is None:
        return None, None
    idx = np.flatnonzero(bound >= 0)
    if idx.size == 0:
        return np.empty(0), np.empty(0, dtype=np.int64)
    k_series = 1.0 / (1.0 / k_head_actin + 1.0 / k_head_spring)
    bead_tags = bound[idx].astype(np.int64)
    head_tags = np.array([action._head_global_tag(int(h)) for h in idx],
                         dtype=np.int64)
    r = np.linalg.norm(pos_byTag[head_tags] - pos_byTag[bead_tags], axis=1)
    s = grip[idx]
    F_series = k_series * np.clip(np.minimum(s, r), 0.0, None)
    return F_series, bead_tags


def aggregation_ledger(sim, *, R_cell: float, p_myo, myosin_action=None,
                       g_soft_gate: float | None = None, n_planes: int = 12) -> dict:
    """Compute the KU-3.5-active aggregation ledger for a live HOOMD snapshot.

    Extracts the harmonic-bond set (the same one the gate's method-of-planes
    sums — discovered from the integrator's bond force), builds the per-bond
    ``(k, r0, channel, T)`` arrays, optionally reads the ``MyosinStepUpdater``
    grip state for the series Hill load, and delegates to
    ``aggregation_ledger_from_bonds``.

    Args:
        sim: live HOOMD ``Simulation`` (after at least one production step).
        R_cell: cortex shell radius [m].
        p_myo: ``ResolvedCortexMyosin`` (k_head_actin, k_head_spring, F_stall).
        myosin_action: optional live ``MyosinStepUpdater`` (for the series load).
        g_soft_gate: optional authoritative gate g_soft [N/m] to report.
        n_planes: method-of-planes cut count.

    Returns:
        The ledger dict (see ``aggregation_ledger_from_bonds``).
    """
    btypes = list(sim.state.bond_types)
    # tag-ordered positions (matches _tagpos / the gate's frame).
    with sim.state.cpu_local_snapshot as s:
        pos = np.asarray(s.particles.position, dtype=np.float64)
        tg = np.asarray(s.particles.tag)
        inv = np.empty_like(tg)
        inv[tg] = np.arange(tg.size)
        pos_byTag = pos[inv].copy()
        bg = np.asarray(s.bonds.group, dtype=np.int64).copy()
        bt = np.asarray(s.bonds.typeid, dtype=np.int64).copy()

    # Per-type (k, r0) from the live integrator bond force (same discovery as
    # the gate MOP — no hard-coded params, picks up mesoscale scaling).
    bond_force = None
    for f in sim.operations.integrator.forces:
        if hasattr(f, "params") and any(t in btypes for t in getattr(f, "params", {})):
            bond_force = f
            break
    if bond_force is None:
        raise RuntimeError("no harmonic bond force found on the integrator")
    k_by_tid = np.zeros(len(btypes))
    r0_by_tid = np.zeros(len(btypes))
    for i, tname in enumerate(btypes):
        try:
            kp = bond_force.params[tname]
            k_by_tid[i] = float(kp["k"])
            r0_by_tid[i] = float(kp["r0"])
        except Exception:
            pass
    channel = np.array([classify_channel(btypes[t]) for t in bt], dtype=object)

    rA = pos_byTag[bg[:, 0]]
    rB = pos_byTag[bg[:, 1]]
    L = np.linalg.norm(rB - rA, axis=1)
    T = k_by_tid[bt] * (L - r0_by_tid[bt])     # signed bond tension

    # Series Hill load aligned to attach bonds (bond col0 = head tag).
    F_series_attach = None
    if myosin_action is not None:
        F_ser, bead_tags = _series_loads_from_action(
            myosin_action, pos_byTag, p_myo.k_head_actin, p_myo.k_head_spring)
        if F_ser is not None and F_ser.size:
            # Map series loads onto attach bonds by head tag.
            is_attach = channel == CH_ATTACH
            attach_head_tags = bg[is_attach, 0]
            bound = getattr(myosin_action, "_head_bound_to_actin")
            idx = np.flatnonzero(bound >= 0)
            head_tag_to_F = {
                int(myosin_action._head_global_tag(int(h))): float(F_ser[j])
                for j, h in enumerate(idx)
            }
            F_series_attach = np.array(
                [head_tag_to_F.get(int(ht), 0.0) for ht in attach_head_tags])

    F_stall = float(getattr(p_myo, "F_stall_per_head", F_STALL_BASE))
    return aggregation_ledger_from_bonds(
        pos=pos_byTag, bg=bg, channel=channel, T=T, R_cell=R_cell,
        F_stall=F_stall, F_series=F_series_attach, g_soft_gate=g_soft_gate,
        n_planes=n_planes)


# --------------------------------------------------------------------------- #
# Human-readable formatter
# --------------------------------------------------------------------------- #
def format_ledger(L: dict, *, indent: str = "    ") -> str:
    """Compact multi-line summary of an aggregation ledger for run logs."""
    g = L["generation"]; h = L["hill"]; t = L["tension"]; e = L["efficiency"]
    mNpm = 1e3
    pN = 1e12
    lines = []
    lines.append(f"{indent}[AGG] n_engaged={g['n_engaged']} "
                 f"Σ|F_head|={g['sumF_head_N']*pN:.3g}pN "
                 f"(mean {g['mean_F_head_N']*pN:.3g}pN, max {g['max_F_head_N']*pN:.3g}pN); "
                 f"Σ|F|/band={g['sumF_over_band']:.3g}×")
    if h.get("load_source", "none") != "none":
        lines.append(f"{indent}      Hill[{h['load_source']}] F/F_stall "
                     f"mean={h['mean_F_over_Fstall']:.3g} max={h['max_F_over_Fstall']:.3g} "
                     f"valid(≤1)={h['frac_hill_valid']*100:.0f}% "
                     f"super-stall={h['frac_super_stall']*100:.0f}%")
    gate = t["g_soft_gate_Npm"]
    gate_s = f"{gate*mNpm:.3e}" if gate is not None else "n/a"
    lines.append(f"{indent}      g_soft: gate={gate_s} mop_total={t['g_mop_total_Npm']*mNpm:.3e} "
                 f"mop_attach={t['g_mop_attach_Npm']*mNpm:.3e} "
                 f"ik_attach={t['g_ik_attach_Npm']*mNpm:.3e} "
                 f"ceiling={t['g_ceiling_attach_Npm']*mNpm:.3e} mN/m")
    lines.append(f"{indent}      η_agg={e['eta_agg']:.3g} (tangential frac "
                 f"{t['mean_tangential_fraction']:.3g}) η_medium={e['eta_medium']:.3g} "
                 f"η_total={e['eta_total']:.3g} ceiling/band={e['ceiling_over_band']:.3g}×")
    lines.append(f"{indent}      → {L['verdict']}")
    return "\n".join(lines)


# --------------------------------------------------------------------------- #
# Figure (visualize-at-closeout rule, CLAUDE.md)
# --------------------------------------------------------------------------- #
def make_aggregation_figure(ledgers: list[dict], out_path, *,
                            band: tuple[float, float] = BAND,
                            title: str | None = None) -> bool:
    """Render the KU-3.5-active aggregation ledger (one fig, 3 panels).

    Panel A: the LOSS CHAIN at the last sample — coherent ceiling → attach IK
    (after the radial/η_agg loss) → realized g_soft (after the medium/η_medium
    loss) → gate g_soft, against the KU-3.5 band. Makes the binding constraint
    visually obvious (e.g. ceiling itself far below band ⇒ generation-limited).
    Panel B: recruitment + ceiling/band over samples. Panel C: η_agg / η_medium
    / tangential fraction over samples (the aggregation-geometry health).

    Args:
        ledgers: per-sample ledger dicts (``aggregation_ledger_from_bonds``).
        out_path: PNG path.
        band: KU-3.5 active band (N/m).
        title: optional suptitle.

    Returns:
        True on success, False if matplotlib unavailable / no data.
    """
    ledgers = [L for L in ledgers if L is not None]
    if not ledgers:
        return False
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as e:  # noqa: BLE001
        print(f"    (aggregation figure skipped: {e})")
        return False
    mNpm = 1e3
    pN = 1e12
    s_idx = np.arange(1, len(ledgers) + 1)
    last = ledgers[-1]
    t = last["tension"]; e = last["efficiency"]

    fig, (axA, axB, axC) = plt.subplots(1, 3, figsize=(16, 4.8))

    # Panel A — loss chain (last sample).
    stages = ["Σ|F_head|\nceiling", "× η_agg\n(radial)", "× η_medium\n(realized)",
              "gate\ng_soft"]
    vals = [t["g_ceiling_attach_Npm"], t["g_ik_attach_Npm"], t["g_ik_total_Npm"],
            (t["g_soft_gate_Npm"] if t["g_soft_gate_Npm"] is not None
             else t["g_mop_total_Npm"])]
    vals_mNpm = [max(abs(v) * mNpm, 1e-12) for v in vals]
    axA.axhspan(band[0] * mNpm, band[1] * mNpm, color="tab:green", alpha=0.18,
                label=f"KU-3.5 band [{band[0]*mNpm:.2f}, {band[1]*mNpm:.2f}]")
    bars = axA.bar(range(len(stages)), vals_mNpm,
                   color=["tab:purple", "tab:blue", "tab:cyan", "tab:red"])
    axA.set_yscale("log")
    axA.set_xticks(range(len(stages)))
    axA.set_xticklabels(stages, fontsize=8)
    axA.set_ylabel("cortical tension  γ  (mN/m)")
    axA.set_title("Aggregation loss chain (last sample)\n"
                  f"η_agg={e['eta_agg']:.2g}  η_medium={e['eta_medium']:.2g}  "
                  f"ceiling/band={e['ceiling_over_band']:.2g}×")
    for b, v in zip(bars, vals_mNpm):
        axA.text(b.get_x() + b.get_width() / 2, v, f"{v:.1e}",
                 ha="center", va="bottom", fontsize=7)
    axA.legend(fontsize=8, loc="upper right")
    axA.grid(True, which="both", axis="y", alpha=0.25)

    # Panel B — recruitment + ceiling/band over samples.
    n_eng = [L["generation"]["n_engaged"] for L in ledgers]
    ceil_band = [L["efficiency"]["ceiling_over_band"] for L in ledgers]
    axB.plot(s_idx, n_eng, "o-", color="tab:purple", label="n_engaged heads")
    axB.set_xlabel("sample")
    axB.set_ylabel("n_engaged heads", color="tab:purple")
    axB.tick_params(axis="y", labelcolor="tab:purple")
    axB2 = axB.twinx()
    axB2.plot(s_idx, ceil_band, "s--", color="tab:orange",
              label="coherent ceiling / band")
    axB2.axhline(1.0, color="k", ls=":", lw=1.0)
    axB2.set_yscale("log")
    axB2.set_ylabel("coherent ceiling / band", color="tab:orange")
    axB2.tick_params(axis="y", labelcolor="tab:orange")
    axB.set_title("Generation: recruitment & ceiling vs band\n"
                  "(ceiling/band < 1 ⇒ generation-limited)")
    axB.grid(True, alpha=0.25)

    # Panel C — aggregation-geometry health over samples.
    eta_agg = [L["efficiency"]["eta_agg"] for L in ledgers]
    eta_med = [L["efficiency"]["eta_medium"] for L in ledgers]
    tan_frac = [L["tension"]["mean_tangential_fraction"] for L in ledgers]
    axC.plot(s_idx, eta_agg, "o-", color="tab:blue", label="η_agg (tangential)")
    axC.plot(s_idx, eta_med, "s-", color="tab:cyan", label="η_medium (no cancel)")
    axC.plot(s_idx, tan_frac, "^--", color="tab:gray",
             label="mean tangential fraction")
    axC.axhline(1.0, color="k", ls=":", lw=1.0)
    axC.set_ylim(-0.05, 1.25)
    axC.set_xlabel("sample")
    axC.set_ylabel("efficiency (dimensionless)")
    axC.set_title("Aggregation geometry health\n(≈1 ⇒ NOT the wall here)")
    axC.legend(fontsize=8, loc="lower right")
    axC.grid(True, alpha=0.25)

    fig.suptitle(title or ("KU-3.5-ACTIVE force-aggregation ledger — "
                           f"verdict: {last['verdict'].split(':')[0]}"), fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=130)
    plt.close(fig)
    print(f"    aggregation figure -> {out_path}")
    return True


# --------------------------------------------------------------------------- #
# Self-test fixtures (known force → known aggregation verdict)
# --------------------------------------------------------------------------- #
def _fib_sphere(n: int, R: float) -> np.ndarray:
    i = np.arange(n) + 0.5
    phi = np.arccos(1 - 2 * i / n)
    theta = np.pi * (1 + 5 ** 0.5) * i
    return np.stack([np.sin(phi) * np.cos(theta),
                     np.sin(phi) * np.sin(theta),
                     np.cos(phi)], axis=1) * R


def _tangential_shell(n: int, R: float, k_nn: int = 6):
    """k-NN tangential mesh on a sphere (returns pos, bg, length)."""
    pos = _fib_sphere(n, R)
    d2 = np.sum((pos[:, None, :] - pos[None, :, :]) ** 2, axis=2)
    np.fill_diagonal(d2, np.inf)
    nn = np.argsort(d2, axis=1)[:, :k_nn]
    bonds = set()
    for a in range(n):
        for b in nn[a]:
            bonds.add((min(a, int(b)), max(a, int(b))))
    return pos, np.asarray(sorted(bonds), dtype=np.int64)


def self_test(verbose: bool = True) -> bool:
    """Fixture-validate the estimator: known force → known aggregation verdict.

    Three controlled cases on a synthetic R=5µm shell prove the ledger
    localises the loss to the right branch, plus the IK-vs-MOP calibration.
    Returns True iff all checks pass.
    """
    R = 5.0e-6
    ok = True

    def check(name, cond):
        nonlocal ok
        ok = ok and bool(cond)
        if verbose:
            print(f"    [{'PASS' if cond else 'FAIL'}] {name}")

    # ---- Calibration: corrected MOP ≈ IK on an isotropic tensile shell -----
    pos, bg = _tangential_shell(1200, R)
    u, Ln, sin2 = bond_geometry(pos, bg)
    T0 = 1e-13
    T = np.full(bg.shape[0], T0)
    g_mop = mop_tension(pos, bg, T, R)
    g_ik = ik_hoop(Ln, sin2, T, R)
    if verbose:
        print(f"\n[1] CALIBRATION isotropic tensile shell ({bg.shape[0]} bonds)")
        print(f"    corrected MOP = {g_mop*1e3:.4e} mN/m   IK = {g_ik*1e3:.4e} mN/m"
              f"   MOP/IK = {g_mop/g_ik:.4f}")
    check("MOP ≈ IK (calibration, ratio∈[0.85,1.15])", 0.85 <= g_mop / g_ik <= 1.15)
    check("tangential mesh sin²ψ ≈ 1", sin2.mean() > 0.9)

    # ---- Case A: fully tangential tensile attach bonds → η≈1, gen-or-ok ----
    # T scaled so the coherent ceiling CLEARS the band — this isolates the η
    # (geometry) test from the generation check (which would otherwise fire
    # first). F_stall set equal so the rows stay Hill-valid.
    chan = np.full(bg.shape[0], CH_ATTACH, dtype=object)
    T_big = np.full(bg.shape[0], 2e-10)
    ledA = aggregation_ledger_from_bonds(
        pos=pos, bg=bg, channel=chan, T=T_big, R_cell=R, F_stall=2e-10)
    if verbose:
        print("\n[2] CASE A — tangential tensile attach (coherent):")
        print(format_ledger(ledA))
    check("A: η_agg ≈ 1 (tangential)", ledA["efficiency"]["eta_agg"] > 0.85)
    check("A: η_medium ≈ 1 (no medium)", abs(ledA["efficiency"]["eta_medium"] - 1.0) < 0.05)
    check("A: verdict not radial/medium",
          "RADIAL" not in ledA["verdict"] and "MEDIUM" not in ledA["verdict"])

    # ---- Case B: RADIAL attach bonds → η_agg≈0, radial verdict -------------
    # Build radial bonds: each from a shell point inward toward center.
    n_b = 800
    outer = _fib_sphere(n_b, R)
    inner = outer * 0.7                      # 30% inward (radial spoke)
    pos_r = np.vstack([outer, inner])
    bg_r = np.stack([np.arange(n_b), np.arange(n_b) + n_b], axis=1)
    chan_r = np.full(n_b, CH_ATTACH, dtype=object)
    Lr = np.linalg.norm(pos_r[bg_r[:, 1]] - pos_r[bg_r[:, 0]], axis=1)
    T_r = np.full(n_b, 5e-12)
    ledB = aggregation_ledger_from_bonds(
        pos=pos_r, bg=bg_r, channel=chan_r, T=T_r, R_cell=R, F_stall=5e-12)
    if verbose:
        print("\n[3] CASE B — radial attach (spokes toward centre):")
        print(format_ledger(ledB))
    check("B: η_agg ≈ 0 (radial)", ledB["efficiency"]["eta_agg"] < 0.15)
    check("B: verdict RADIAL or GENERATION (low realized)",
          "RADIAL" in ledB["verdict"] or "GENERATION" in ledB["verdict"])

    # ---- Case C: tangential attach + equal-and-opposite compressive xlinks -
    # Medium cancellation: same geometry, xlink bonds in compression cancel.
    # Generation is adequate (attach ceiling > band, as in Case A); the 90%
    # compressive xlink medium is what kills the realized tension → MEDIUM.
    bg_c = np.vstack([bg, bg])               # duplicate the tangential mesh
    chan_c = np.array([CH_ATTACH] * bg.shape[0] + [CH_XLINK] * bg.shape[0],
                      dtype=object)
    T_c = np.concatenate([np.full(bg.shape[0], 2e-10),      # attach: tension
                          np.full(bg.shape[0], -1.8e-10)])   # xlink: compression
    pos_c = pos
    ledC = aggregation_ledger_from_bonds(
        pos=pos_c, bg=bg_c, channel=chan_c, T=T_c, R_cell=R, F_stall=2e-10)
    if verbose:
        print("\n[4] CASE C — tangential attach + compressive xlink medium:")
        print(format_ledger(ledC))
    check("C: η_agg ≈ 1 (attach still tangential)", ledC["efficiency"]["eta_agg"] > 0.85)
    check("C: η_medium ≪ 1 (medium cancels)", ledC["efficiency"]["eta_medium"] < 0.5)
    check("C: verdict MEDIUM", "MEDIUM" in ledC["verdict"])

    # ---- Case D: tiny per-head force → GENERATION-LIMITED ------------------
    T_small = np.full(bg.shape[0], 1e-15)
    ledD = aggregation_ledger_from_bonds(
        pos=pos, bg=bg, channel=chan, T=T_small, R_cell=R, F_stall=5e-12)
    if verbose:
        print("\n[5] CASE D — tiny per-head force (under-generation):")
        print(format_ledger(ledD))
    check("D: ceiling < band", ledD["efficiency"]["ceiling_over_band"] < 1.0)
    check("D: verdict GENERATION-LIMITED", "GENERATION" in ledD["verdict"])

    if verbose:
        print(f"\n{'='*60}\nSELF-TEST {'PASSED' if ok else 'FAILED'}\n{'='*60}")
    return ok


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--self-test", action="store_true",
                    help="run the fixture validation (known force → verdict)")
    args = ap.parse_args()
    if args.self_test:
        return 0 if self_test() else 1
    ap.print_help()
    return 0


if __name__ == "__main__":
    sys.exit(main())
