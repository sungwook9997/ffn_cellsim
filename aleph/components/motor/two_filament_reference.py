r"""Two-filament sarcomere DECISION reference — the cleanest bug-vs-density test (host NumPy, NO Warp).

Host-side acceptance oracle (pure NumPy) for the single cleanest experiment that separates a *motor /
series-compliance bug* from an honest *density floor*: ONE production-topology NMII minifilament wedged
between TWO straight ANTIPARALLEL actin filaments (outward barbed ends, far ends clamped). This module is
the ANALYTIC ground truth the native Warp runner
(:mod:`aleph.components.motor.native_gates.ng1_two_filament_stall`) is gated against.

THE DECISION (why this topology):
  * The native gate measures the two clamp reactions and the crossbridge→arm→backbone→opposite-filament
    TRANSMISSION RATIO. In a lossless *collinear static* spring chain, every series bond carries the
    identical axial tension (force balance at each internal node), so the force generated at the active
    crossbridge arrives UNDIMINISHED at the opposite filament's clamp: **transmission ratio = 1** (proven
    numerically by :func:`collinear_transmission_proof`, not asserted). The series stiffness
    ``k_eff = (1/k_xb + 1/k_head_arm)^{-1}`` (:func:`series_stiffness`) is NOT a force loss — it only sets
    how far the chain STRETCHES to build that force (``F_side / k_eff``); it is the compliance the native
    gate must probe. Real transmission < 1 appears ONLY through GEOMETRY / misalignment (backbone bending,
    head-arm swing off the actin axis, a non-central crossbridge reaction) — modes the 3-D device topology
    + inner relaxation expose and this collinear analytic cannot.
  * VERDICT RULE (native vs this reference):
      - transmission ≪ 1  ⇒ a motor / series-compliance BUG, independent of density → fix the motor.
      - transmission ≈ 1 AND per-side reaction ≈ this analytic ``F_side`` ⇒ the motor is faithful; carry the
        low per-minifilament force to the FULL-CELL DENSITY question (:func:`gamma_ceiling`,
        :func:`sweep_gap_claims`) — the cortical tension is floored by the sparse *measured* minifilament
        areal density (Nie 2015 ~0.625 µm⁻²), NOT by a broken motor. NEVER add heads/density to lift γ (§6.2).

THE EMERGENT PER-SIDE FORCE (reuses the ensemble-stall oracle; NO reinvented stall math):
  At isometric (v_slide → 0) each BOUND head steps to its own stall ``f_stall`` (v_head = 0), so the per-side
  reaction is ``F_side = N_side · phi_b(f_stall) · f_stall`` where the engaged fraction
  ``phi_b = k_on/(k_on + k_off0 e^{f_stall/f0})`` EMERGES from the Bell kinetics — it is NOT the naive product
  ``N_side · f_stall`` (P2). :func:`per_side_stall_force` computes this by CALLING
  :func:`aleph.components.motor.ensemble_stall_analytic.ensemble_stall_meanfield` and cross-checks ``phi_b`` against
  :func:`aleph.components.motor.bell_kinetics_analytic.engaged_fraction_steady`; the per-head working-stroke strain
  at stall is taken from :func:`aleph.components.motor.powerstroke_analytic.stall_abscissa`. At claim_a
  (N_side = 10, f_stall = 0.5 pN, the AFINES values wired in ``ac.cell.assemble``) ``F_side ≈ 4.96 pN``.

THE FULL-CELL γ CEILING (reuses the sourced density + the codebase dipole→tension convention):
  A minifilament is a contractile force dipole of arm ``L_bb`` (Billington 0.301 µm) carrying per-side force
  ``F_side``. At an areal density ``rho`` its 2-D active tension is the force-dipole areal stress
  ``gamma = rho · F_side · L_bb / 2`` — the SAME convention as
  ``ff.network_contractility`` (``sigma_dipole = Σ f·L_dip / Area / 2``) and consistent with the
  method-of-planes estimator ``ff.gamma_estimator`` for a small dipole (arm ≪ R). ``rho`` is taken at the
  ONLY direct cortical NMII datum, Nie 2015 (``ff.gamma_floor.NIE2015_DENSITY_UM2 = 0.625 µm⁻²``), so the γ
  ceiling is a REPORTED consequence of the sourced sparse density, never tuned.

⚠ report-not-tune (§6.2 / CLAUDE.md): every magnitude here is a params_i0b3 GAP (N_side, f_stall, k_xb) or a
sourced datum (Nie density, L_bb, f0). :func:`sweep_gap_claims` tabulates the WHOLE GAP grid so PI sees the
finding across claims; nothing is chosen to hit a band. The host gate
(``tests/ac/motor/test_two_filament_sarcomere.py``) asserts these constants do NOT drift from ``ac.cell.assemble``
/ ``ff.gamma_floor``.

References:
  Stam-Hocky 2015 (minifilament); Billington 2013 JBC (L_bb = 301 nm, N/minifilament); Kovacs 2003 (duty);
  Bell 1978 (slip); Nie et al. 2015 Cytoskeleton 72:29 (cortical NMII areal density); Truong-Quang 2021
  (engaged overlap). Mirrors ``ac.cell.assemble`` NMII_* constants + ``ff.gamma_floor`` / ``ff.gamma_estimator``.
"""

from __future__ import annotations

import numpy as np
import numpy.typing as npt

from aleph.components.motor.bell_kinetics_analytic import bell_f0_from_x_beta, engaged_fraction_steady
from aleph.components.motor.ensemble_stall_analytic import ensemble_stall_meanfield
from aleph.components.motor.powerstroke_analytic import stall_abscissa

__all__ = [
    "N_SIDE_CLAIM_A",
    "F_STALL_CLAIM_A",
    "K_XB",
    "K_HEAD_ARM",
    "K_BACKBONE",
    "K_ON",
    "K_OFF0",
    "F0_PN",
    "L_BB_UM",
    "R_CELL_UM",
    "NIE2015_DENSITY_UM2",
    "NIE2015_DENSITY_RANGE_UM2",
    "MCF7_ACTIVE_BAND_PN_UM",
    "N_SIDE_CLAIMS",
    "F_STALL_CLAIMS",
    "per_side_stall_force",
    "series_stiffness",
    "full_transmission_path_stiffness",
    "ideal_transmission_ratio",
    "collinear_transmission_proof",
    "two_filament_isometric",
    "gamma_ceiling",
    "sweep_gap_claims",
    "format_sweep_table",
]

# ── SOURCED constants — mirror ac.cell.assemble NMII_* / params_i0b3.yaml (test asserts NO drift) ─────────
# claim_a magnitudes (params_i0b3.yaml GAP rows; the AFINES branch wired into ac.cell.assemble at t0):
N_SIDE_CLAIM_A: int = 10        # assemble.NMII_N_SIDE   (AFINES; params_i0b3 N_side claim_a)
F_STALL_CLAIM_A: float = 0.5    # assemble.NMII_F_STALL_TEST (AFINES/Kovacs; params_i0b3 F_stall_head claim_a) [pN]
# mechanical stiffnesses (assemble.py TEST values; k_xb is the params_i0b3 MASTER-knob GAP, mid 100–1000 band):
K_XB: float = 1.0e3             # assemble.NMII_K_XB_TEST      crossbridge (head↔actin) stiffness [pN/µm]
K_HEAD_ARM: float = 1.0e2       # assemble.NMII_K_HEAD_ARM_TEST head↔backbone arm stiffness [pN/µm]
K_BACKBONE: float = 1.0e3       # assemble.NMII_K_BACKBONE_TEST backbone PER-BOND stiffness [pN/µm] (per contract)
N_BB: int = 14                  # assemble.NMII_N_BB  backbone beads ⇒ (N_BB−1) bonds in series (end-to-end soft)
# Bell kinetics (assemble.py; k_on GAP-provisional, k_off0/f0 draft):
K_ON: float = 50.0              # assemble.NMII_KON_TEST   per-head attach rate [1/s]
K_OFF0: float = 0.35            # assemble.NMII_KOFF0      Bell zero-force off-rate [1/s]
F0_PN: float = 4.28e-3 / 0.6e-3  # assemble.NMII_F0 = kBT/x_beta ≈ 7.13 pN (Bell characteristic force)
# minifilament geometry (Billington 2013 EM; assemble.py):
L_BB_UM: float = 0.301          # assemble.NMII_L_BB_UM   backbone contour = force-dipole arm [µm]
R_CELL_UM: float = 7.5          # assemble.R_CELL_UM      MCF7 outer radius (I0-A ratified) [µm]

# ── SOURCED density + band — mirror ff.gamma_floor / ff.gamma_estimator (test asserts NO drift) ───────────
# The ONLY direct cortical NMII minifilament areal density (Nie et al. 2015 Cytoskeleton 72:29, HeLa medial
# cortex, intensity-calibrated; non-MCF7, LOW-MEDIUM confidence). ff.gamma_floor.NIE2015_DENSITY_UM2.
NIE2015_DENSITY_UM2: float = 0.625
NIE2015_DENSITY_RANGE_UM2: tuple[float, float] = (0.31, 0.94)
# MCF7 myosin-ACTIVE cortical-tension target band [pN/µm] = ACTIVE_FRACTION(0.70) · Chugh-2017 total band
# (350,650). ff.gamma_estimator.active_band_pn_um(SALBREUX_BAND_PN_UM). 1 pN/µm = 1e-3 mN/m.
MCF7_ACTIVE_BAND_PN_UM: tuple[float, float] = (245.0, 455.0)

# ── GAP claim grids (params_i0b3.yaml N_side / F_stall_head conflicted rows) ──────────────────────────────
N_SIDE_CLAIMS: tuple[int, int] = (10, 28)        # a: AFINES model ; b: Billington 2013 structural (28–30)
F_STALL_CLAIMS: tuple[float, float] = (0.5, 2.0)  # a: AFINES/Kovacs ; b: Billington/KB-3.18 [pN]


def per_side_stall_force(
    n_side: int = N_SIDE_CLAIM_A,
    f_stall: float = F_STALL_CLAIM_A,
    *,
    k_on: float = K_ON,
    k_off0: float = K_OFF0,
    f0: float = F0_PN,
    k_xb: float = K_XB,
) -> dict[str, float]:
    """Emergent per-side isometric stall force ``F_side`` [pN], reusing the ensemble-stall oracle.

    At isometric each bound head sits at its own stall ``f_stall`` (Hill v → 0), so the per-side reaction is
    ``F_side = N_side · phi_b(f_stall) · f_stall`` with the engaged fraction
    ``phi_b = k_on/(k_on + k_off0 e^{f_stall/f0})`` — EMERGENT from the Bell kinetics, NOT the naive product
    (P2). Computed by calling :func:`ensemble_stall_analytic.ensemble_stall_meanfield` (the oracle) and
    cross-checked against :func:`bell_kinetics_analytic.engaged_fraction_steady`; the per-head working-stroke
    strain at stall is :func:`powerstroke_analytic.stall_abscissa`.

    Args:
        n_side: explicit heads per anti-parallel half-filament ``N_side`` (params_i0b3 GAP).
        f_stall: per-head isometric stall force ``f_stall`` [pN] (params_i0b3 GAP).
        k_on: per-head attach rate [1/s].
        k_off0: Bell zero-force off-rate [1/s].
        f0: Bell characteristic force [pN].
        k_xb: crossbridge stiffness [pN/µm] (for the working-stroke strain diagnostic).

    Returns:
        dict with ``f_side`` (emergent per-side reaction [pN]), ``phi_b`` (engaged fraction),
        ``n_engaged`` (``N_side·phi_b``), ``f_naive`` (``N_side·f_stall``), ``oracle_f_ensemble`` (the
        ``ensemble_stall_meanfield`` value F_side is taken from), ``phi_b_crosscheck`` (independent
        ``engaged_fraction_steady``), and ``working_stroke_strain_nm`` (per-head ``f_stall/k_xb`` [nm]).
    """
    oracle = ensemble_stall_meanfield(n_side, f_stall, k_on, k_off0, f0)
    phi_b_crosscheck = float(engaged_fraction_steady(k_on, f_stall, k_off0, f0))
    strain_um = stall_abscissa(f_stall, k_xb, r0_xb=0.0)   # per-head working-stroke strain at stall [µm]
    return {
        "f_side": float(oracle["f_ensemble"]),
        "phi_b": float(oracle["phi_b"]),
        "n_engaged": float(oracle["n_engaged"]),
        "f_naive": float(oracle["f_naive"]),
        "oracle_f_ensemble": float(oracle["f_ensemble"]),
        "phi_b_crosscheck": phi_b_crosscheck,
        "working_stroke_strain_nm": float(strain_um * 1.0e3),
    }


def series_stiffness(k_xb: float = K_XB, k_head_arm: float = K_HEAD_ARM) -> float:
    """Per-side series compliance ``k_eff = (1/k_xb + 1/k_head_arm)^{-1}`` [pN/µm] (collinear).

    The crossbridge and head-arm springs in series along one arm's load path. With the assemble.py test
    values (k_xb = 1000, k_head_arm = 100) this is ≈ 90.9 pN/µm. It is NOT a force loss (a static series
    chain transmits ratio = 1); it is the softness that sets the chain's stretch ``F_side / k_eff`` under
    load — the compliance the native gate must probe.

    Args:
        k_xb: crossbridge stiffness [pN/µm].
        k_head_arm: head↔backbone arm stiffness [pN/µm].

    Returns:
        Series stiffness [pN/µm].
    """
    if k_xb <= 0.0 or k_head_arm <= 0.0:
        raise ValueError("stiffnesses must be positive")
    return 1.0 / (1.0 / k_xb + 1.0 / k_head_arm)


def full_transmission_path_stiffness(
    k_xb: float = K_XB, k_head_arm: float = K_HEAD_ARM, k_backbone: float = K_BACKBONE, n_bb: int = N_BB
) -> float:
    """Full clamp-to-clamp series stiffness [pN/µm] of the two-filament load path.

    The path traverses ``crossbridge · arm · backbone · arm · crossbridge`` in series. ⚠ The backbone is the
    DISCRETE ``(n_bb−1)``-bond chain the engine builds (``harmonic_bond_kernel`` applies the per-bond
    ``k_backbone`` to each of the ``n_bb−1`` bonds), so its END-TO-END stiffness is ``k_backbone/(n_bb−1)`` —
    the earlier single-``1/k_backbone`` term treated it as ONE rigid spring and was WRONG (Codex cross-audit):

        (2/k_xb + 2/k_head_arm + (n_bb−1)/k_backbone)^{-1}

    With the assemble per-bond values this is a SOFT ≈ 6.0 pN/µm (backbone-dominated), not the 43.5 the
    single-spring form gave. The transmission ratio is still 1 (series force balance is count-independent); only
    the compliance / chain stretch changes — and the discrete gate must relax over that larger stretch.

    Args:
        k_xb: crossbridge stiffness [pN/µm].
        k_head_arm: head↔backbone arm stiffness [pN/µm].
        k_backbone: backbone PER-BOND stiffness [pN/µm].
        n_bb: backbone bead count ⇒ ``n_bb−1`` bonds in series.

    Returns:
        Full-path series stiffness [pN/µm].
    """
    if k_xb <= 0.0 or k_head_arm <= 0.0 or k_backbone <= 0.0:
        raise ValueError("stiffnesses must be positive")
    if n_bb < 2:
        raise ValueError("n_bb must be >= 2")
    return 1.0 / (2.0 / k_xb + 2.0 / k_head_arm + (n_bb - 1) / k_backbone)


def ideal_transmission_ratio() -> float:
    """The lossless collinear-static transmission ratio = ``1.0`` (exact).

    A static series spring chain balances force at every internal node, so the identical axial tension
    propagates from the active crossbridge to the opposite filament's clamp — the generated force is
    transmitted UNDIMINISHED. Real loss (< 1) enters only via geometry/misalignment (backbone bending, arm
    swing, non-central reaction), which the native Warp gate exposes and this collinear reference cannot.
    :func:`collinear_transmission_proof` computes this ratio numerically rather than asserting it.

    Returns:
        ``1.0``.
    """
    return 1.0


def collinear_transmission_proof(
    k_springs: npt.ArrayLike | None = None,
    f_generated: float = 1.0,
    active_index: int = 0,
) -> dict[str, object]:
    """Numerically PROVE the ideal collinear transmission ratio = 1 via a 1-D series-spring equilibrium.

    Builds a chain of ``n`` springs (default the full clamp-to-clamp path
    ``[k_xb, k_head_arm, k_backbone, k_head_arm, k_xb]``) between two FIXED walls (the two actin clamps),
    imposes an active rest-length preload on spring ``active_index`` (the myosin power stroke) sized so that
    spring's tension equals ``f_generated``, solves the internal-node force balance ``K x = b``, and reads
    the tension in every spring. Series force balance forces a UNIFORM tension, so the tension at the far
    wall equals the tension at the near wall: ratio = 1 to machine precision.

    Args:
        k_springs: series stiffnesses [pN/µm] along the load path (default the full transmission path).
        f_generated: the active crossbridge tension to reproduce [pN] (sets the preload; ratio is
            independent of this value).
        active_index: which spring carries the active preload (the crossbridge).

    Returns:
        dict with ``ratio`` (far-wall tension / near-wall tension, = 1), ``tension_uniform_pN`` (the common
        tension), ``tension_spread_pN`` (max−min across springs, ≈ 0), ``k_series_pN_per_um`` (the full
        series stiffness), ``preload_um`` (the active rest-length offset), and ``tensions_pN`` (per-spring).
    """
    if k_springs is None:
        # the DISCRETE clamp-to-clamp chain: crossbridge · arm · (N_BB−1 backbone bonds) · arm · crossbridge
        # (the backbone is a bond CHAIN, not one spring — matches full_transmission_path_stiffness / the engine).
        k_springs = [K_XB, K_HEAD_ARM] + [K_BACKBONE] * (N_BB - 1) + [K_HEAD_ARM, K_XB]
    k = np.asarray(k_springs, dtype=np.float64)
    n = k.shape[0]
    if n < 1:
        raise ValueError("need at least one spring")
    if not (0 <= active_index < n):
        raise ValueError("active_index out of range")
    if np.any(k <= 0.0):
        raise ValueError("stiffnesses must be positive")

    k_series = 1.0 / np.sum(1.0 / k)
    preload = f_generated / k_series          # rest-length offset so the common tension == f_generated

    # rest lengths: all zero except the active spring, which "wants" to be shorter by `preload` (contractile)
    rest = np.zeros(n)
    rest[active_index] = -preload

    # internal nodes 1..n-1 (nodes 0 and n are the fixed walls at x=0 and x=0, an unstretched reference);
    # equilibrium: for node j, k_j (x_j - x_{j-1} - rest_j) = k_{j+1} (x_{j+1} - x_j - rest_{j+1}).
    m = n - 1
    if m == 0:                                # single spring between two coincident walls
        tension = k[0] * (0.0 - 0.0 - rest[0])
        tensions = np.array([-tension])
        return {
            "ratio": 1.0, "tension_uniform_pN": float(abs(tension)),
            "tension_spread_pN": 0.0, "k_series_pN_per_um": float(k_series),
            "preload_um": float(preload), "tensions_pN": tensions.tolist(),
        }
    a = np.zeros((m, m))
    b = np.zeros(m)
    x_left, x_right = 0.0, 0.0                 # both walls fixed at the same reference coordinate
    for row, j in enumerate(range(1, n)):      # node j is internal for j in 1..n-1
        a[row, row] = k[j - 1] + k[j]
        if row > 0:
            a[row, row - 1] = -k[j - 1]
        if row < m - 1:
            a[row, row + 1] = -k[j]
        # node-j force balance: (k_{j-1}+k_j) x_j − k_{j-1} x_{j-1} − k_j x_{j+1}
        #                       = k_{j-1} rest_{j-1} − k_j rest_j
        b[row] = k[j - 1] * rest[j - 1] - k[j] * rest[j]
        if j == 1:
            b[row] += k[j - 1] * x_left
        if j == n - 1:
            b[row] += k[j] * x_right
    x_internal = np.linalg.solve(a, b)
    x = np.concatenate([[x_left], x_internal, [x_right]])
    tensions = k * (np.diff(x) - rest)         # per-spring tension [pN] (signed)
    t_near = tensions[0]
    t_far = tensions[-1]
    ratio = float(t_far / t_near) if t_near != 0.0 else 1.0
    return {
        "ratio": ratio,
        "tension_uniform_pN": float(np.mean(np.abs(tensions))),
        "tension_spread_pN": float(np.max(tensions) - np.min(tensions)),
        "k_series_pN_per_um": float(k_series),
        "preload_um": float(preload),
        "tensions_pN": tensions.tolist(),
    }


def two_filament_isometric(
    n_side: int = N_SIDE_CLAIM_A,
    f_stall: float = F_STALL_CLAIM_A,
    **kwargs: float,
) -> dict[str, object]:
    """The full analytic prediction for the two-filament isometric-stall gate (the native runner's target).

    Bundles the emergent per-side reaction, the two equal-and-opposite clamp reactions, the ideal
    transmission ratio, and the series compliance the native gate probes.

    Args:
        n_side: heads per anti-parallel side (params_i0b3 GAP).
        f_stall: per-head stall force [pN] (params_i0b3 GAP).
        **kwargs: forwarded to :func:`per_side_stall_force` (``k_on``, ``k_off0``, ``f0``, ``k_xb``).

    Returns:
        dict with ``f_side`` [pN], ``clamp_reaction_A_pN`` / ``clamp_reaction_B_pN`` (equal-and-opposite,
        the two filaments' clamps hold ±F_side inward pull), ``ideal_transmission`` (= 1),
        ``k_eff_pN_per_um`` (per-side series), ``k_full_pN_per_um`` (clamp-to-clamp series),
        ``chain_stretch_at_load_nm`` (``F_side/k_eff`` — the compliance the native gate probes), plus the
        pass-through kinetics diagnostics from :func:`per_side_stall_force`.
    """
    f = per_side_stall_force(n_side, f_stall, **kwargs)
    f_side = f["f_side"]
    k_eff = series_stiffness()
    k_full = full_transmission_path_stiffness()
    return {
        "f_side": f_side,
        # the minifilament pulls both actins inward by F_side; each clamp reacts with an equal, outward-
        # directed force. Sign convention: reaction on filament A is +F_side (holding against inward pull),
        # on filament B is −F_side (the antiparallel side) — equal magnitude, opposite sense.
        "clamp_reaction_A_pN": +f_side,
        "clamp_reaction_B_pN": -f_side,
        "ideal_transmission": ideal_transmission_ratio(),
        "k_eff_pN_per_um": k_eff,
        "k_full_pN_per_um": k_full,
        "chain_stretch_at_load_nm": (f_side / k_eff) * 1.0e3,
        **{k: v for k, v in f.items() if k != "f_side"},
    }


def gamma_ceiling(
    f_side: float,
    *,
    density_um2: float = NIE2015_DENSITY_UM2,
    l_dip_um: float = L_BB_UM,
    band_pn_um: tuple[float, float] = MCF7_ACTIVE_BAND_PN_UM,
) -> dict[str, float]:
    """Full-cell cortical-tension ceiling [pN/µm] from the per-side force via the force-dipole areal stress.

    A minifilament is a contractile force dipole (arm ``l_dip``, per-side force ``f_side``); at areal density
    ``density`` the 2-D active tension is ``gamma = density · f_side · l_dip / 2`` — the codebase convention
    (``ff.network_contractility`` ``sigma_dipole = Σ f·L_dip / Area / 2``; consistent with the
    ``ff.gamma_estimator`` method-of-planes for a small dipole, arm ≪ R). At the Nie-2015 sourced density
    this is a REPORTED consequence of the sparse measured density, never tuned.

    Args:
        f_side: emergent per-side isometric force [pN].
        density_um2: minifilament areal density [µm⁻²] (default Nie 2015 0.625).
        l_dip_um: force-dipole arm = backbone contour ``L_bb`` [µm].
        band_pn_um: MCF7 myosin-active cortical-tension target band [pN/µm].

    Returns:
        dict with ``gamma_pn_um`` (ceiling [pN/µm]), ``gamma_mN_per_m`` (= gamma·1e-3),
        ``floor_factor_lo`` / ``floor_factor_hi`` (band_lo/gamma, band_hi/gamma — ×-under the band),
        ``rho_needed_band_mid_um2`` (areal density to reach the band midpoint at this ``f_side``), and
        ``density_um2`` (echoed).
    """
    if f_side < 0.0:
        raise ValueError("f_side must be non-negative")
    if density_um2 <= 0.0 or l_dip_um <= 0.0:
        raise ValueError("density and l_dip must be positive")
    gamma = density_um2 * f_side * l_dip_um / 2.0
    band_lo, band_hi = band_pn_um
    band_mid = 0.5 * (band_lo + band_hi)
    dipole_coeff = f_side * l_dip_um / 2.0                   # gamma per unit areal density [pN·µm]
    rho_needed = band_mid / dipole_coeff if dipole_coeff > 0.0 else float("inf")
    return {
        "gamma_pn_um": gamma,
        "gamma_mN_per_m": gamma * 1.0e-3,
        "floor_factor_lo": band_lo / gamma if gamma > 0.0 else float("inf"),
        "floor_factor_hi": band_hi / gamma if gamma > 0.0 else float("inf"),
        "rho_needed_band_mid_um2": rho_needed,
        "density_um2": density_um2,
    }


def sweep_gap_claims(
    n_side_claims: tuple[int, ...] = N_SIDE_CLAIMS,
    f_stall_claims: tuple[float, ...] = F_STALL_CLAIMS,
    *,
    density_um2: float = NIE2015_DENSITY_UM2,
) -> list[dict[str, float]]:
    """Tabulate ``F_side`` and the R=7.5 µm-cell γ ceiling over the (N_side, f_stall) GAP grid.

    The whole point of the DECISION: even the most generous GAP claim leaves the cortical tension floored
    at the sourced Nie-2015 density — so a faithful low per-side force is a DENSITY finding, not a motor bug.

    Args:
        n_side_claims: the ``N_side`` GAP claims to sweep (default (10, 28)).
        f_stall_claims: the ``f_stall`` GAP claims to sweep (default (0.5, 2.0)).
        density_um2: minifilament areal density [µm⁻²] (default Nie 2015 0.625).

    Returns:
        list of row dicts (one per (N_side, f_stall) combination), each with ``n_side``, ``f_stall_pN``,
        ``phi_b``, ``f_naive_pN``, ``f_side_pN``, ``gamma_pn_um``, ``gamma_mN_per_m``,
        ``floor_factor_lo`` / ``floor_factor_hi``, ``rho_needed_band_mid_um2``.
    """
    rows: list[dict[str, float]] = []
    for n_side in n_side_claims:
        for f_stall in f_stall_claims:
            fs = per_side_stall_force(n_side, f_stall)
            gc = gamma_ceiling(fs["f_side"], density_um2=density_um2)
            rows.append({
                "n_side": int(n_side),
                "f_stall_pN": float(f_stall),
                "phi_b": fs["phi_b"],
                "f_naive_pN": fs["f_naive"],
                "f_side_pN": fs["f_side"],
                "gamma_pn_um": gc["gamma_pn_um"],
                "gamma_mN_per_m": gc["gamma_mN_per_m"],
                "floor_factor_lo": gc["floor_factor_lo"],
                "floor_factor_hi": gc["floor_factor_hi"],
                "rho_needed_band_mid_um2": gc["rho_needed_band_mid_um2"],
            })
    return rows


def format_sweep_table(rows: list[dict[str, float]] | None = None) -> str:
    """Render the GAP-claim sweep as a fixed-width text table (for the host gate + the native report).

    Args:
        rows: rows from :func:`sweep_gap_claims` (default: a fresh sweep of the standard grid).

    Returns:
        A multi-line table string.
    """
    if rows is None:
        rows = sweep_gap_claims()
    band_lo, band_hi = MCF7_ACTIVE_BAND_PN_UM
    header = (
        f"  N_side  f_stall   phi_b   F_naive   F_side    gamma      gamma      floor-x    rho_need\n"
        f"          [pN]              [pN]      [pN]      [pN/um]    [mN/m]     under band [/um^2]"
    )
    lines = [
        "Two-filament sarcomere DECISION sweep — F_side (emergent) + R=7.5um-cell gamma ceiling",
        f"  density = Nie2015 {NIE2015_DENSITY_UM2}/um^2 (range {NIE2015_DENSITY_RANGE_UM2})  |  "
        f"L_bb = {L_BB_UM} um  |  MCF7 active band = {band_lo:.0f}-{band_hi:.0f} pN/um",
        header,
    ]
    for r in rows:
        lines.append(
            f"  {r['n_side']:>5d}  {r['f_stall_pN']:>6.2f}  {r['phi_b']:>7.4f}  "
            f"{r['f_naive_pN']:>7.2f}  {r['f_side_pN']:>7.3f}  {r['gamma_pn_um']:>8.4f}  "
            f"{r['gamma_mN_per_m']:>9.2e}  "
            f"{r['floor_factor_lo']:>5.0f}-{r['floor_factor_hi']:<5.0f} {r['rho_needed_band_mid_um2']:>8.1f}"
        )
    k_eff = series_stiffness()
    k_full = full_transmission_path_stiffness()
    lines.append(
        f"  k_eff (per-side crossbridge+arm series) = {k_eff:.3f} pN/um   "
        f"[full clamp-to-clamp path = {k_full:.3f} pN/um]"
    )
    lines.append(f"  ideal collinear transmission ratio = {ideal_transmission_ratio():.1f}  "
                 f"(loss < 1 is geometric only — the native Warp gate exposes it)")
    return "\n".join(lines)


def _main() -> None:
    """Print the full decision reference (runnable: ``python -m aleph.components.motor.two_filament_reference``)."""
    claim_a = two_filament_isometric()
    print("=" * 96)
    print("TWO-FILAMENT SARCOMERE DECISION REFERENCE (host NumPy analytic; NO Warp)")
    print("=" * 96)
    print(f"claim_a (N_side={N_SIDE_CLAIM_A}, f_stall={F_STALL_CLAIM_A} pN, AFINES; ac.cell.assemble t0):")
    print(f"  phi_b (Bell-emergent engaged fraction) = {claim_a['phi_b']:.6f}")
    print(f"  F_naive = N_side*f_stall               = {N_SIDE_CLAIM_A * F_STALL_CLAIM_A:.3f} pN")
    print(f"  F_side  = N_side*phi_b*f_stall          = {claim_a['f_side']:.4f} pN   <-- emergent per-side")
    print(f"  clamp reactions  A / B                  = "
          f"{claim_a['clamp_reaction_A_pN']:+.4f} / {claim_a['clamp_reaction_B_pN']:+.4f} pN "
          f"(equal & opposite)")
    print(f"  per-head working-stroke strain at stall = {claim_a['working_stroke_strain_nm']:.3f} nm "
          f"(= f_stall/k_xb)")
    print(f"  k_eff (crossbridge+arm series)          = {claim_a['k_eff_pN_per_um']:.3f} pN/um")
    print(f"  chain stretch at load  = F_side/k_eff   = {claim_a['chain_stretch_at_load_nm']:.3f} nm")
    proof = collinear_transmission_proof(f_generated=claim_a["f_side"])
    print(f"  ideal transmission ratio (proven)       = {proof['ratio']:.12f}  "
          f"(tension spread {proof['tension_spread_pN']:.2e} pN; k_series {proof['k_series_pN_per_um']:.3f})")
    print("-" * 96)
    print(format_sweep_table())
    print("=" * 96)


if __name__ == "__main__":
    _main()
