"""Resolve KU-derived parameters from the primary scales in a YAML config.

All formulas come from cited KUs; no fitting, no hand-tuning.

Geometry (KU-1.27, 2D isotropic Mikado):
    Line density:           ρ_L = N_f · L_fiber / L_box²        [1/m]
    Crossings per fiber:    ⟨n_int⟩ = 2 ρ_L L_fiber / π          (Buffon-Sylvester)
    Mikado segment length:  ℓ_c = π / (2 ρ_L)                   (Poisson spacing on a fiber)

Time integration (KU-1.26, overdamped Langevin):
    Per-bead Stokes drag:   γ_b = 6 π η_water · bead_radius
    Relaxation times:       τ_xl       = γ_b / k_xl
                            τ_stretch  = γ_b ℓ₀ / μ
                            τ_bend     = γ_b ℓ₀³ / κ
    CFL bound:              dt < α · τ_min, α = cfl_safety_factor (default 0.1).

We *match the Mikado segment length ℓ_c to a target* (typically the
biological mesh ξ from KU-1.7) — i.e. invert:
    N_f = π · L_box² / (2 · ℓ_c · L_fiber)

The emergent coordination ⟨z⟩ = backbone_z + ⟨n_int⟩ / beads_per_fiber
is then a **predicted** quantity, not a tuning knob. Derivation:
each crossing creates **one** cross-link that touches two fibers
(one bead each), so total links L = N · ⟨n_int⟩ / 2 and the per-bead
contribution to z is 2L / (N · beads) = ⟨n_int⟩ / beads.

For Phase 1 defaults (ℓ_c = ξ = 2 μm) this predicts ⟨z⟩ ≈ 2.6, which
is below the KU-1.3 biological collagen reference of 3.2. The gap is
real biology: 3D collagen has covalent cross-links (lysyl oxidase) and
bundled-fiber connectivity that a 2D geometric Mikado does not
reproduce. Matching ⟨z⟩=3.2 instead of ℓ_c=2 μm would require
ℓ_c ≈ 1.25 μm (denser than the experimental porosity). We follow the
Storm-MacKintosh convention and prioritise ξ; the gate accepts the
sub-isostatic band [2.5, 3.9] so the emergent 2.6 passes while being
honestly distinguishable from the biological 3.2.

Biological ξ (KU-1.7, c^(-1/2)) and Mikado ℓ_c (KU-1.27) are *distinct*
quantities; the config stores them separately and only equates them by
modelling convention (Storm & MacKintosh 2005, Broedersz & MacKintosh
2014).
"""

from __future__ import annotations

import math
from pathlib import Path
from typing import Any

import yaml


def _ku127_n_fibers_from_segment_length(
    L_box: float, target_segment_length: float, L_fiber: float
) -> int:
    """Invert KU-1.27 Mikado segment-length formula.

    ℓ_c = π / (2 ρ_L)  ⇒  N = π L_box² / (2 ℓ_c L_fiber)
    """
    return int(math.ceil(math.pi * L_box * L_box
                         / (2.0 * target_segment_length * L_fiber)))


def _line_density(n_fibers: int, L_fiber: float, L_box: float) -> float:
    return n_fibers * L_fiber / (L_box * L_box)


def _mikado_segment_length(line_density: float) -> float:
    return math.pi / (2.0 * line_density)


def _expected_crossings_per_fiber(line_density: float, L_fiber: float) -> float:
    return 2.0 * line_density * L_fiber / math.pi


def _backbone_z(beads_per_fiber: int) -> float:
    if beads_per_fiber < 2:
        return 0.0
    return (2.0 * (beads_per_fiber - 2) + 2.0) / beads_per_fiber


def resolve(cfg: dict[str, Any]) -> dict[str, Any]:
    """Fill the ``ecm.derived.*`` block of a Phase 1 config in place."""
    ecm = cfg["ecm"]
    L_box = float(ecm["L_box"])
    L_fiber = float(ecm["L_fiber"])
    target_ell_c = float(ecm["target_segment_length"])
    beads = int(ecm["beads_per_fiber"])
    cost_cap = int(ecm["acceptance"]["n_fibers_max"])
    demo = bool(ecm.get("demo_mode", False))

    n_fibers_kb = _ku127_n_fibers_from_segment_length(L_box, target_ell_c, L_fiber)
    if n_fibers_kb > cost_cap:
        if not demo:
            raise ValueError(
                f"KU-1.27 derivation requires n_fibers={n_fibers_kb}, exceeding "
                f"cost cap {cost_cap}. Set ecm.demo_mode=true to override."
            )
        n_fibers = cost_cap
        resolved_at = "demo_mode_cost_cap"
    else:
        n_fibers = n_fibers_kb
        resolved_at = "KU-1.27_segment_match"

    rho_L = _line_density(n_fibers, L_fiber, L_box)
    ell_c_predicted = _mikado_segment_length(rho_L)
    n_int_per_fiber = _expected_crossings_per_fiber(rho_L, L_fiber)
    backbone_z = _backbone_z(beads)
    # XL contribution per bead: each crossing → 1 link touching 1 bead
    # on each of 2 fibers ⇒ ⟨z⟩_xl = ⟨n_int⟩ / beads_per_fiber.
    expected_z = backbone_z + n_int_per_fiber / beads

    derived = ecm.setdefault("derived", {})
    derived["n_fibers"] = n_fibers
    derived["rest_length"] = L_fiber / (beads - 1)
    derived["backbone_z"] = backbone_z
    derived["line_density"] = rho_L
    derived["mikado_ell_c_predicted"] = ell_c_predicted
    derived["expected_xl_per_fiber"] = n_int_per_fiber
    derived["expected_total_z"] = expected_z
    derived["fibers_resolved_at"] = resolved_at

    # ---- KU-1.26 time-integration derivation ----
    eta = float(ecm["water_viscosity"])
    a = float(ecm["bead_radius"])
    mu = float(ecm["stretching_modulus"])
    kappa = float(ecm["bending_modulus"])
    k_xl = float(ecm["xl_stiffness"])
    ell_0 = derived["rest_length"]

    gamma_b = 6.0 * math.pi * eta * a
    tau_xl = gamma_b / k_xl
    tau_stretch = gamma_b * ell_0 / mu
    tau_bend = gamma_b * ell_0 ** 3 / kappa
    tau_min = min(tau_xl, tau_stretch, tau_bend)

    derived["gamma_b"] = gamma_b
    derived["tau_xl"] = tau_xl
    derived["tau_stretch"] = tau_stretch
    derived["tau_bend"] = tau_bend
    derived["tau_min"] = tau_min

    dyn = ecm.get("dynamics", {})
    alpha = float(dyn.get("cfl_safety_factor", 0.1))
    dt_cfl = alpha * tau_min
    if dyn.get("dt") is None:
        # Stay strictly below the CFL bound (1 % margin) so the
        # `dt < α·τ_min` gate has room.
        dyn["dt"] = 0.99 * dt_cfl
    derived["dt_cfl"] = dt_cfl
    ecm["dynamics"] = dyn

    return cfg


def load_config(path: str | Path) -> dict[str, Any]:
    with open(path, "r") as f:
        cfg = yaml.safe_load(f)
    return resolve(cfg)


# ---------------------------------------------------------------------- #
# Bridge / Unit 2.1                                                      #
# ---------------------------------------------------------------------- #


def _pereverzev_peak(
    k_off_slip: float, F_s: float, k_off_catch: float, F_c: float
) -> tuple[float, float]:
    """Closed-form (F*, τ_max) for the Pereverzev catch-slip form (KU-2.5).

    F* = F_s F_c / (F_s + F_c) · ln(k_c F_s / (k_s F_c))   [N]
    τ_max = 1 / k_off(F*)                                   [s]

    Returns (nan, nan) if the bond is slip-only (no catch peak).
    """
    arg = (k_off_catch * F_s) / (k_off_slip * F_c)
    if arg <= 1.0:
        return float("nan"), float("nan")
    F_star = (F_s * F_c) / (F_s + F_c) * math.log(arg)
    k_at_peak = (
        k_off_slip * math.exp(F_star / F_s)
        + k_off_catch * math.exp(-F_star / F_c)
    )
    return F_star, 1.0 / k_at_peak


def resolve_bridge(cfg: dict[str, Any]) -> dict[str, Any]:
    """Fill the ``bridge.derived.*`` block of a Phase 1 Unit 2.1 config in place.

    All formulas come from KU-1.21, KU-2.4, KU-2.5, KU-2.8, KU-2.18.
    No fitting, no hand-tuning.
    """
    b = cfg["bridge"]
    sub = b["substrate"]
    mc = b["motor_clutch"]
    bond = b["catch_bond"]

    # ---- Substrate stiffness (KU-1.21) ----
    E = float(sub["young_modulus"])
    nu = float(sub["poisson_ratio"])
    a = float(sub["contact_radius"])
    k_sub = math.pi * (E / (1.0 - nu * nu)) * a

    # ---- Catch-slip peak (KU-2.5 closed-form) ----
    F_star, tau_max = _pereverzev_peak(
        float(bond["k_off_slip"]), float(bond["F_s"]),
        float(bond["k_off_catch"]), float(bond["F_c"]),
    )

    # ---- Bond event rate at the slip characteristic force (sets dt budget) ----
    F_s = float(bond["F_s"])
    k_off_at_Fs = (
        float(bond["k_off_slip"]) * math.exp(1.0)
        + float(bond["k_off_catch"]) * math.exp(-F_s / float(bond["F_c"]))
    )
    bond_event_rate_max = max(float(mc["k_on"]), k_off_at_Fs)
    cfl_safe_dt = 0.1 / bond_event_rate_max

    # ---- KU-2.8 biphasic optimum ----
    # Two analytics, faithful to KU-2.8 and to the catch-bond extension
    # cited there (Alonso-Matilla 2023 Biophys J).
    #
    # (a) KU-2.8 motor-stall heuristic: k_sub* = N_m F_stall / (v_un τ_max).
    #     This is an *upper bound* derived for slip-only bonds and tends
    #     to underestimate the simulated peak when catch bonds stabilise
    #     engagement (Bangasser 2017 Nat Cell Biol, fig 2).
    # (b) Matched-stiffness heuristic (Bangasser 2013 Biophys J,
    #     Alonso-Matilla 2023): the peak occurs when the substrate
    #     stiffness matches the *effective* clutch-bundle stiffness,
    #     ~ ½ N_clutches k_int (half-engaged steady state). This is the
    #     better predictor for catch-stabilised systems and is what the
    #     biphasic regression test compares against.
    n_clutches = float(mc["n_clutches"])
    k_int = float(mc["k_int"])
    if math.isfinite(tau_max) and tau_max > 0.0:
        k_sub_star_motor = (
            float(mc["n_motors"]) * float(mc["F_stall_per_motor"])
            / (float(mc["v_unloaded"]) * tau_max)
        )
    else:
        k_sub_star_motor = float("nan")
    k_sub_star_matched = 0.5 * n_clutches * k_int   # Bangasser 2013

    def _ksub_to_E(k: float) -> float:
        return k / (math.pi * a) * (1.0 - nu * nu) if math.isfinite(k) else float("nan")

    derived = b.setdefault("derived", {})
    derived["substrate_stiffness"] = k_sub
    derived["catch_peak_force"] = F_star
    derived["catch_peak_lifetime"] = tau_max
    derived["bond_event_rate_max"] = bond_event_rate_max
    derived["cfl_safe_dt"] = cfl_safe_dt
    derived["biphasic_kSubStar_motor"] = k_sub_star_motor
    derived["biphasic_EStar_motor"] = _ksub_to_E(k_sub_star_motor)
    derived["biphasic_kSubStar_matched"] = k_sub_star_matched
    derived["biphasic_EStar_matched"] = _ksub_to_E(k_sub_star_matched)
    # Primary analytic the gate compares against (see Sanity Gate notes).
    derived["biphasic_kSubStar"] = k_sub_star_matched
    derived["biphasic_EStar"] = _ksub_to_E(k_sub_star_matched)
    return cfg


def load_bridge_config(path: str | Path) -> dict[str, Any]:
    with open(path, "r") as f:
        cfg = yaml.safe_load(f)
    return resolve_bridge(cfg)


# ---------------------------------------------------------------------- #
# Bridge / Unit 2.2 (mature FA + ECM adapter)                            #
# ---------------------------------------------------------------------- #


def resolve_bridge_22(cfg: dict[str, Any], unit1_cfg_path: str | Path) -> dict[str, Any]:
    """Resolve a Phase 1 Unit 2.2 config in place.

    Extends :func:`resolve_bridge` with the maturation derived block
    (vinculin steady state, k_int_eff, talin rate at 5 pN) and pulls
    Worker A's μ, κ from the ECM Unit 1.1 YAML so the adapter speaks the
    same ECM language as Worker A's network.
    """
    cfg = resolve_bridge(cfg)
    b = cfg["bridge"]

    # Cross-load ECM moduli from Worker A's Unit 1.1 config so the adapter
    # uses the *same* fibre mechanics. Honour explicit overrides (non-null
    # values in the Unit 2.2 yaml) — only fill blanks.
    ecm = b.setdefault("ecm_adapter", {})
    if ecm.get("stretching_modulus") is None or ecm.get("bending_modulus") is None:
        with open(unit1_cfg_path, "r") as f:
            unit1 = yaml.safe_load(f)
        if ecm.get("stretching_modulus") is None:
            ecm["stretching_modulus"] = float(unit1["ecm"]["stretching_modulus"])
        if ecm.get("bending_modulus") is None:
            ecm["bending_modulus"] = float(unit1["ecm"]["bending_modulus"])

    # Maturation derived quantities.
    vin = b["vinculin"]
    n_unfolded_target = 1                                # Phase 1: 1-state max
    N_vin_ss = (
        float(vin["k_rec"]) * n_unfolded_target * float(vin["N_free"])
        / max(float(vin["k_diss"]), 1e-30)
    )
    mc = b["motor_clutch"]
    k_int_eff_ss = float(mc["k_int"]) * (1.0 + float(vin["alpha"]) * N_vin_ss)

    talin = b["talin"]
    k_unfold_5pN = float(talin["k_u0"]) * math.exp(
        5.0e-12 * float(talin["dx_star"]) / float(talin["kT"])
    )

    derived = b.setdefault("derived", {})
    derived["vinculin_steady_state_at_n1"] = N_vin_ss
    derived["k_int_eff_at_steady_vin"] = k_int_eff_ss
    derived["talin_unfold_rate_at_5pN"] = k_unfold_5pN
    return cfg


def load_bridge_22_config(
    path: str | Path,
    unit1_cfg_path: str | Path | None = None,
) -> dict[str, Any]:
    """Load and resolve a Unit 2.2 config; cross-load Unit 1.1 by default."""
    path = Path(path)
    if unit1_cfg_path is None:
        unit1_cfg_path = path.parent / "phase1_unit1.yaml"
    with open(path, "r") as f:
        cfg = yaml.safe_load(f)
    return resolve_bridge_22(cfg, unit1_cfg_path)
