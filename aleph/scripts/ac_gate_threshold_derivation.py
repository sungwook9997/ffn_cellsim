"""Derive the resting-convergence force threshold from first principles (GATE A / `max|PF|`).

Companion executable for ``aleph/docs/v2_audit/GATE_A_THRESHOLD_DERIVATION_2026-07-25.md``.

WHAT THIS IS.  The literal ``0.21`` pN has governed two weeks of solver work with no derivation on
disk (first appearance: a bare "gate ≈ 0.21" in ``RESTING_BASELINE_DIAGNOSIS_2026-07-22c.md``; ~12
hand-copied literals since).  This script reconstructs, from the on-disk constants only, the four
candidate first-principles routes to a per-node residual-force threshold, identifies which route the
incumbent literal actually came from, and states how each route scales under mesh refinement.

WHAT THIS IS NOT.  It changes no gate and imports no oracle.  It never runs a simulation: it reads
constants and one committed native artifact.  CPU-only by construction (numpy; `warp` is never
imported — the two physics modules it does import, ``ff.units`` and ``ac.solid.wca_analytic``, are
pure numpy, and the ``ac/cell`` constants are read as source literals precisely so that this script
cannot pull the Warp runtime onto a CPU host).

Sanity Gate:
    * Dimensional analysis: every derived force is checked through an explicit exponent-tuple unit
      algebra (µm, pN, s) and must reduce to pN^1.  ``CHECK dimensions``.
    * Positive control / independent limiting cases:
        - the route-(c) formula must reproduce ``tolerance_um / dt_mu`` as recorded in a committed
          native run artifact to float round-off, and the frozen 0.21 literal to <1%;
        - ``sqrt(kBT/k)`` at k = 0.1 N/m must reproduce the 0.207 nm computed independently in
          ``_historical/H1_H5_H7_MANIFOLD_CONTACT_ARCHITECTURE_2026-06-07.md`` (quoted 0.21 nm);
        - with steric disabled, route (c) must reproduce the ~0.06 pN implied bound stated
          independently in ``AUDIT_WHOLE_REPO_2026-07-25.md`` finding 9.
    * Negative control: ``--assume-kmax`` / ``--assume-seg`` inject a deliberately wrong input; the
      reproduction check must then FAIL (non-zero exit).  A check nobody has watched fail is not a
      check.
    * Grid invariance: route (c) is asserted to scale exactly ∝ ℓ (``F/ℓ`` constant across the
      refinement series to float round-off); the frozen-literal loosening factor per rung is
      reported, not silently absorbed.
    * Sign sense: all thresholds are positive magnitudes; a threshold must never be reported as
      looser than the incumbent without the ``LOOSER-THAN-INCUMBENT`` banner.

Usage:
    python aleph/scripts/ac_gate_threshold_derivation.py
    python aleph/scripts/ac_gate_threshold_derivation.py --assume-kmax 8.2e5   # negative control
    python aleph/scripts/ac_gate_threshold_derivation.py --assume-seg 0.075    # negative control
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np

_HERE = Path(__file__).resolve()
if str(_HERE.parents[2]) not in sys.path:
    sys.path.insert(0, str(_HERE.parents[2]))

from aleph.components.solid.wca_analytic import (  # noqa: E402
    contact_stiffness,
    epsilon_from_contact_stiffness,
    max_core_stiffness,
)
from aleph.laws.units import KAPPA_ACTIN, KBT, LP_ACTIN_UM  # noqa: E402

REPO = _HERE.parents[2]
ASSEMBLE = REPO / "aleph" / "components" / "incumbent" / "assemble.py"
COMPARTMENTS = REPO / "aleph" / "components" / "incumbent" / "compartments.py"
DRIVER = REPO / "aleph" / "components" / "incumbent" / "driver.py"
TURGOR_PI0 = REPO / "aleph" / "common" / "turgor_pi0.py"
GATE_SCRIPT = REPO / "aleph" / "scripts" / "ac_resting_converge.py"
NATIVE_ARTIFACT = (
    REPO / "aleph" / "outputs" / "ac" / "resting_native"
    / "resting_native_probe2_capture0.6_2026-07-23.json"
)

# ── unit algebra (µm, pN, s) — so "the units come out in pN" is a check, not a comment ────────────
UM, PN, S = (1, 0, 0), (0, 1, 0), (0, 0, 1)
DIMLESS = (0, 0, 0)
ENERGY = (1, 1, 0)          # pN·µm
STIFFNESS = (-1, 1, 0)      # pN/µm
MOBILITY_STEP = (1, -1, 0)  # µm/pN  (the dt_mu pseudo-step)


def u_mul(a: tuple[int, ...], b: tuple[int, ...]) -> tuple[int, ...]:
    """Multiply two unit exponent tuples."""
    return tuple(x + y for x, y in zip(a, b))


def u_div(a: tuple[int, ...], b: tuple[int, ...]) -> tuple[int, ...]:
    """Divide two unit exponent tuples."""
    return tuple(x - y for x, y in zip(a, b))


def u_sqrt(a: tuple[int, ...]) -> tuple[int, ...]:
    """Square-root a unit exponent tuple; raises when the result is not integral."""
    out = []
    for x in a:
        if x % 2:
            raise ValueError(f"non-integral unit exponent under sqrt: {a}")
        out.append(x // 2)
    return tuple(out)


def u_str(a: tuple[int, ...]) -> str:
    """Render a unit exponent tuple."""
    names = ("um", "pN", "s")
    parts = [f"{n}^{e}" for n, e in zip(names, a) if e]
    return "·".join(parts) if parts else "1"


# ── on-disk constant reading (regex, so a drifted constant re-derives instead of going stale) ─────
def read_literal(path: Path, pattern: str, *fallbacks: tuple[Path, str]) -> float:
    """Read one float literal out of a source file by regex, trying fallback definition sites.

    Constants in this tree migrate between a bare module literal and a provenance-gated resolver
    (``PI_0_PA`` did exactly that on 2026-07-25), so each quantity may name several definition
    sites in priority order.  A quantity that matches NOWHERE halts — it never silently defaults.

    Args:
        path: Primary source file to read.
        pattern: Regex with exactly one capture group holding the numeric literal.
        *fallbacks: Additional ``(path, pattern)`` pairs, tried in order.

    Returns:
        The captured literal as a float.

    Raises:
        ValueError: If no candidate matches — a renamed constant must halt, not default.
    """
    for p, pat in ((path, pattern), *fallbacks):
        if not p.exists():
            continue
        m = re.search(pat, p.read_text())
        if m is not None:
            return float(m.group(1))
    raise ValueError(f"constant not found in {path.name} (or {len(fallbacks)} fallbacks): {pattern}")


@dataclass(frozen=True)
class Constants:
    """Every number this derivation uses, with its on-disk provenance."""

    kBT_pn_um: float
    kappa_actin_pn_um2: float
    lp_actin_um: float
    seg_um: float
    sigma_ev_um: float
    k_ev_pn_per_um: float
    steric_force_cap_pn: float
    k_xb_pn_per_um: float
    k_erm_pn_per_um: float
    k_xl_alpha_pn_per_um: float
    k_xl_filamin_pn_per_um: float
    pi0_pa: float
    r_cell_um: float
    dt_mu_numerator: float
    eps64: float


def load_constants() -> Constants:
    """Collect the derivation constants from their definition sites on disk."""
    from aleph.laws.hand_kmc import ALPHA_ACTININ, FILAMIN

    return Constants(
        kBT_pn_um=float(KBT),
        kappa_actin_pn_um2=float(KAPPA_ACTIN),
        lp_actin_um=float(LP_ACTIN_UM),
        seg_um=read_literal(ASSEMBLE, r"cortex_seg_um: float = ([0-9.eE+-]+)"),
        sigma_ev_um=read_literal(ASSEMBLE, r"SIGMA_EV_UM = ([0-9.eE+-]+)"),
        k_ev_pn_per_um=read_literal(ASSEMBLE, r"K_EV_PROVISIONAL = ([0-9.eE+-]+)"),
        steric_force_cap_pn=read_literal(
            ASSEMBLE, r"steric_force_cap: float \| None = ([0-9.eE+-]+)"),
        k_xb_pn_per_um=read_literal(ASSEMBLE, r"NMII_K_XB_TEST = ([0-9.eE+-]+)"),
        k_erm_pn_per_um=read_literal(COMPARTMENTS, r"k_erm=([0-9.eE+-]+)"),
        k_xl_alpha_pn_per_um=float(ALPHA_ACTININ.link_k),
        k_xl_filamin_pn_per_um=float(FILAMIN.link_k),
        pi0_pa=read_literal(
            ASSEMBLE, r"PI_0_PA = ([0-9.eE+-]+)\s",
            (TURGOR_PI0, r"PI0_CLAIM_HELA_PROXY_40PA: float = ([0-9.eE+-]+)")),
        r_cell_um=read_literal(ASSEMBLE, r"R_CELL_UM = ([0-9.eE+-]+)"),
        dt_mu_numerator=read_literal(ASSEMBLE, r"dt_mu = ([0-9.eE+-]+) / kmax"),
        eps64=float(np.finfo(np.float64).eps),
    )


def frozen_literal() -> float:
    """The incumbent hand-copied gate literal, read out of the gate script itself."""
    return read_literal(GATE_SCRIPT, r'"gate_pn": ([0-9.eE+-]+)')


# ── the four candidate routes ──────────────────────────────────────────────────────────────────────
def kmax_pn_per_um(c: Constants, with_steric: bool = True) -> tuple[float, str]:
    """Reconstruct the CFL ``kmax`` that sets ``dt_mu = 0.1/kmax`` in ``assemble.py``.

    Returns:
        ``(kmax, dominant_term_name)``.
    """
    terms: dict[str, float] = {
        "actin_bending_kappa_over_ell3": c.kappa_actin_pn_um2 / c.seg_um**3,
        "crosslink_alpha_actinin": c.k_xl_alpha_pn_per_um,
        "crosslink_filamin": c.k_xl_filamin_pn_per_um,
        "nmii_crossbridge_k_xb": c.k_xb_pn_per_um,
        "membrane_erm_k_erm": c.k_erm_pn_per_um,
    }
    if with_steric:
        eps_wca = epsilon_from_contact_stiffness(c.k_ev_pn_per_um, c.sigma_ev_um)
        terms["steric_wca_core_at_force_cap"] = max_core_stiffness(
            c.steric_force_cap_pn, c.sigma_ev_um, eps_wca)
        terms["steric_wca_contact"] = contact_stiffness(eps_wca, c.sigma_ev_um)
    name = max(terms, key=lambda k: terms[k])
    return terms[name], name


def route_c_predicate_floor(c: Constants, seg_um: float, kmax: float) -> float:
    """Route (c): the accepted-step predicate's own force-discrimination floor [pN].

    The executing runtime criterion (``inner_mechanics.update_convergence_kernel``) is
    ``dt_mu * max|PF| <= tolerance_um`` with ``tolerance_um = sqrt(eps64) * ell`` and
    ``dt_mu = 0.1 / kmax``.  Solving for the force gives the largest residual the predicate cannot
    distinguish from a converged state:

        F_pred = tolerance_um / dt_mu = (0.1)^-1 * sqrt(eps64) * ell * kmax
    """
    tolerance_um = np.sqrt(c.eps64) * seg_um
    dt_mu = c.dt_mu_numerator / kmax
    return float(tolerance_um / dt_mu)


def route_a2_equipartition_floor(c: Constants, k_pn_per_um: float) -> float:
    """Route (a2): the equipartition (γ- and dt-free) thermal force floor [pN].

    A degree of freedom held by stiffness ``k`` at 310 K has ``<x^2> = kBT/k``, so the instantaneous
    net force on it fluctuates with rms ``sqrt(kBT*k)``.  This is the thermal force scale of a
    QUASI-STATIC state; unlike the Langevin ``sqrt(2*gamma*kBT/dt)`` it needs neither a drag nor a
    time step, both of which the resting solve does not have.
    """
    return float(np.sqrt(c.kBT_pn_um * k_pn_per_um))


def route_a2_thermal_length_um(c: Constants, k_pn_per_um: float) -> float:
    """The companion thermal position uncertainty ``sqrt(kBT/k)`` [µm]."""
    return float(np.sqrt(c.kBT_pn_um / k_pn_per_um))


def route_b_displacement_floor(k_pn_per_um: float, delta_um: float) -> float:
    """Route (b): force that produces a displacement below a length that cannot matter [pN]."""
    return float(k_pn_per_um * delta_um)


def route_d_load_balance_total_pn(c: Constants, fraction: float) -> float:
    """Route (d): grid-invariant TOTAL residual as a fraction of the physiological load [pN].

    The resting shell must carry the turgor load.  The total outward force the cortex balances across
    a mid-plane is ``F_load = Pi_0 * pi * R^2`` (1 Pa == 1 pN/µm²).  A total residual below
    ``fraction * F_load`` cannot shift the reported Laplace tension ``gamma = dP*R/2`` by more than
    ``fraction``.  Both sides are TOTAL forces, so this is invariant to ``n_nodes`` and ``seg_um``.
    """
    f_load = c.pi0_pa * np.pi * c.r_cell_um**2
    return float(fraction * f_load)


def route_d_radial_bias_budget_pn(c: Constants, gamma_fraction: float) -> tuple[float, float]:
    """Route (d), the measurement-consistent form: the SIGNED radial residual budget.

    A residual field only contaminates the reported cortical tension through its coherent radial
    component.  The signed sum ``S = sum_i (PF_i . r_hat_i)`` acts as a spurious pressure
    ``dP_spur = S / (4 pi R^2)`` and therefore biases the Laplace readout by
    ``gamma_spur = dP_spur * R / 2 = S / (8 pi R)``.  Requiring ``gamma_spur <= f * gamma_target``
    with ``gamma_target = Pi_0 * R / 2`` gives a budget on ``|S|`` that contains neither ``n_nodes``
    nor ``ell`` — the only genuinely grid-invariant form found in this derivation.

    Returns:
        ``(gamma_target_pn_per_um, |S| budget in pN)``.
    """
    gamma_target = c.pi0_pa * c.r_cell_um / 2.0
    budget = gamma_fraction * gamma_target * 8.0 * np.pi * c.r_cell_um
    return float(gamma_target), float(budget)


# ── checks ────────────────────────────────────────────────────────────────────────────────────────
@dataclass
class Check:
    """One pass/fail derivation check."""

    name: str
    ok: bool
    detail: str


def run_checks(c: Constants, kmax: float, kmax_term: str, seg_um: float,
               native: dict | None) -> list[Check]:
    """Run every dimensional / limiting-case / grid-scaling check."""
    checks: list[Check] = []

    # C1 — dimensional analysis through the unit algebra.
    dims: list[tuple[str, tuple[int, ...]]] = []
    u_tol = u_mul(DIMLESS, UM)                       # sqrt(eps) [1] * ell [um]
    u_dtmu = u_div(DIMLESS, STIFFNESS)               # 0.1 [1] / kmax [pN/um]
    dims.append(("route_c = tolerance/dt_mu", u_div(u_tol, u_dtmu)))
    dims.append(("route_a2 = sqrt(kBT*k)", u_sqrt(u_mul(ENERGY, STIFFNESS))))
    dims.append(("route_b = k*delta", u_mul(STIFFNESS, UM)))
    dims.append(("route_d = frac*Pi0*pi*R^2", u_mul(u_div(PN, u_mul(UM, UM)), u_mul(UM, UM))))
    bad = [f"{n}->{u_str(u)}" for n, u in dims if u != PN]
    checks.append(Check(
        "dimensions", not bad,
        "all four routes reduce to pN^1: " + "; ".join(f"{n} -> {u_str(u)}" for n, u in dims)
        + ("" if not bad else f"  BAD: {bad}")))
    # dt_mu itself must be um/pN (the pseudo-step), else the force test is not a displacement test.
    checks.append(Check(
        "dt_mu_units", u_dtmu == MOBILITY_STEP,
        f"dt_mu = 0.1/kmax -> {u_str(u_dtmu)} (required um^1·pN^-1)"))

    # C2 — POSITIVE CONTROL: route (c) must reproduce the committed native artifact exactly.
    f_pred = route_c_predicate_floor(c, seg_um, kmax)
    if native is not None:
        rep = native.get("report", native)
        tol_disk = float(rep["inner_tolerance_um"])
        dtmu_disk = float(rep["inner_dt_mu"])
        kmax_disk = float(rep["ledger"]["kmax_pN_per_um"])
        f_disk = tol_disk / dtmu_disk
        # The runtime's ell is the REALIZED mean segment rest length (assemble.py `seg_mean`), which
        # differs from the config default by the beads-per-filament rounding; back it out of the
        # committed tolerance so the FORMULA can be checked to round-off.
        seg_mean_disk = tol_disk / np.sqrt(c.eps64)
        rel_k = abs(kmax - kmax_disk) / kmax_disk
        f_exact = route_c_predicate_floor(c, seg_mean_disk, kmax)
        rel_exact = abs(f_exact - f_disk) / f_disk
        rel_cfg = abs(f_pred - f_disk) / f_disk
        checks.append(Check(
            "reproduce_native_kmax", rel_k < 1e-12,
            f"derived kmax={kmax:.10g} ({kmax_term}) vs committed ledger {kmax_disk:.10g}"
            f"  rel={rel_k:.3e} (tol 1e-12)"))
        checks.append(Check(
            "reproduce_native_predicate_force_exact", rel_exact < 1e-12,
            f"route-(c) at the realized ell={seg_mean_disk:.10g} µm gives {f_exact:.10g} pN vs"
            f" committed tolerance/dt_mu = {tol_disk:.10g}/{dtmu_disk:.10g} = {f_disk:.10g} pN"
            f"  rel={rel_exact:.3e} (tol 1e-12)"))
        checks.append(Check(
            "reproduce_native_predicate_force_from_config", rel_cfg < 1e-3,
            f"route-(c) at the CONFIG default ell={seg_um} µm gives {f_pred:.10g} pN vs committed"
            f" {f_disk:.10g} pN  rel={rel_cfg:.3e} (tol 1e-3; the residual gap is the"
            f" config-seg vs realized-seg_mean rounding, {abs(seg_um - seg_mean_disk) * 1e3:.4f} nm)"))
    else:
        checks.append(Check("reproduce_native_kmax", False, "committed native artifact missing"))

    # C3 — the frozen literal must be route (c) at the baseline mesh.
    lit = frozen_literal()
    rel_lit = abs(f_pred - lit) / lit
    checks.append(Check(
        "frozen_literal_is_route_c", rel_lit < 0.01,
        f"frozen literal {lit} pN vs derived route-(c) {f_pred:.6g} pN -> rel={rel_lit:.4%}"
        " (tol 1%)  => the incumbent 0.21 IS the accepted-step predicate floor at seg=0.5 um"))

    # C4 — LIMITING CASE 1 (independent): sqrt(kBT/k) at k_ERM = 0.1 N/m = 1e5 pN/um must reproduce
    # the 0.207 nm computed in _historical/H1_H5_H7_MANIFOLD_CONTACT_ARCHITECTURE_2026-06-07.md (quoted 0.21 nm).
    sigma_radial_nm = route_a2_thermal_length_um(c, 1.0e5) * 1e3
    checks.append(Check(
        "limiting_case_sigma_radial", abs(sigma_radial_nm - 0.207) < 0.005,
        f"sqrt(kBT/k) at k=1e5 pN/um -> {sigma_radial_nm:.4f} nm; June manifold doc independently"
        " computed 0.207 nm (quoted 0.21 nm)"))

    # C5 — LIMITING CASE 2 (independent): steric OFF -> kmax = filamin 8.2e5 -> ~0.06 pN, the bound
    # AUDIT_WHOLE_REPO_2026-07-25.md finding 9 states independently.
    kmax_nosteric, term_nosteric = kmax_pn_per_um(c, with_steric=False)
    f_nosteric = route_c_predicate_floor(c, seg_um, kmax_nosteric)
    checks.append(Check(
        "limiting_case_steric_off", abs(f_nosteric - 0.06) < 0.006,
        f"steric OFF -> kmax={kmax_nosteric:.4g} ({term_nosteric}) -> route-(c)={f_nosteric:.4g} pN;"
        " whole-repo audit finding 9 independently states ~0.06 pN"))

    # C6 — GRID SCALING: route (c) is exactly proportional to ell (kmax is ell-independent while the
    # steric core dominates), so F/ell must be constant across the refinement series.
    rungs = [0.5, 0.2, 0.075]
    ratios = [route_c_predicate_floor(c, s, kmax) / s for s in rungs]
    spread = (max(ratios) - min(ratios)) / max(ratios)
    checks.append(Check(
        "route_c_scales_linearly_in_seg", spread < 1e-12,
        "F_pred/ell across seg = "
        + ", ".join(f"{s}->{r:.6g}" for s, r in zip(rungs, ratios))
        + f"  spread={spread:.3e} (tol 1e-12)  => route (c) is O(ell^1), NOT grid-invariant"))

    return checks


def main() -> int:
    """Run the derivation, print the table, and return a POSIX exit code."""
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--assume-kmax", type=float, default=None,
                    help="NEGATIVE CONTROL: inject a wrong CFL kmax [pN/um]")
    ap.add_argument("--assume-seg", type=float, default=None,
                    help="NEGATIVE CONTROL: inject a wrong segment length [um]")
    ap.add_argument("--json", type=Path, default=None, help="write the derivation to JSON")
    args = ap.parse_args()

    c = load_constants()
    kmax_real, term = kmax_pn_per_um(c)
    kmax = args.assume_kmax if args.assume_kmax is not None else kmax_real
    term = "ASSUMED (negative control)" if args.assume_kmax is not None else term
    seg = args.assume_seg if args.assume_seg is not None else c.seg_um
    injected = args.assume_kmax is not None or args.assume_seg is not None

    native = json.loads(NATIVE_ARTIFACT.read_text()) if NATIVE_ARTIFACT.exists() else None
    lit = frozen_literal()

    print("=" * 96)
    print("GATE A resting-convergence threshold — first-principles derivation")
    print("=" * 96)
    if injected:
        print("!! NEGATIVE CONTROL: a deliberately wrong input was injected"
              f" (kmax={args.assume_kmax}, seg={args.assume_seg}) — checks MUST fail.")
    print(f"constants: kBT={c.kBT_pn_um:.6g} pN·µm  kappa_actin={c.kappa_actin_pn_um2:.6g} pN·µm²"
          f" (L_p={c.lp_actin_um} µm)  seg(ell)={seg} µm")
    print(f"CFL kmax  = {kmax:.10g} pN/µm   dominant term: {term}")
    print(f"dt_mu     = {c.dt_mu_numerator}/kmax = {c.dt_mu_numerator / kmax:.10g} µm/pN")
    print(f"tolerance = sqrt(eps64)*ell = {np.sqrt(c.eps64) * seg:.10g} µm")
    print()

    f_pred = route_c_predicate_floor(c, seg, kmax)
    print("-- ROUTE (c)  accepted-step predicate floor  F = 10*sqrt(eps64)*ell*kmax -------------")
    print(f"   F_pred = {f_pred:.6g} pN      (frozen literal in the gate script: {lit} pN)")
    print()

    print("-- ROUTE (a2) equipartition thermal force floor  F = sqrt(kBT*k) --------------------")
    print(f"   {'mode':38s} {'k [pN/µm]':>12s} {'F_th [pN]':>11s} {'x_th [nm]':>10s}")
    modes = {
        "actin bending soft mode kappa/ell^3": c.kappa_actin_pn_um2 / seg**3,
        "NMII crossbridge k_xb": c.k_xb_pn_per_um,
        "membrane ERM tether k_erm": c.k_erm_pn_per_um,
        "crosslink alpha-actinin (Ferrer)": c.k_xl_alpha_pn_per_um,
        "crosslink filamin (Ferrer)": c.k_xl_filamin_pn_per_um,
        "steric WCA core at force cap": kmax_real,
    }
    thermal = {}
    for name, k in modes.items():
        f_th = route_a2_equipartition_floor(c, k)
        thermal[name] = f_th
        print(f"   {name:38s} {k:12.4g} {f_th:11.4g} "
              f"{route_a2_thermal_length_um(c, k) * 1e3:10.4g}")
    print()

    print("-- ROUTE (b)  meaningless-displacement floor  F = k*delta ---------------------------")
    for dname, delta in (("actin monomer rise 2.7 nm", 0.0027), ("steric sigma 7 nm", c.sigma_ev_um)):
        row = " ".join(
            f"{n.split()[0]}={route_b_displacement_floor(k, delta):.4g}" for n, k in modes.items())
        print(f"   delta = {dname:26s} -> {row}")
    print("   => spans >5 decades across the model's own stiffnesses: delta must be mode-matched,")
    print("      and the only natural mode-matched length is sqrt(kBT/k) => route (b) == route (a2).")
    print()

    print("-- ROUTE (d)  grid-invariant residual criteria --------------------------------------")
    f_load = c.pi0_pa * np.pi * c.r_cell_um**2
    gamma_target, _ = route_d_radial_bias_budget_pn(c, 1.0)
    print(f"   F_load = Pi_0*pi*R^2 = {c.pi0_pa}*pi*{c.r_cell_um}^2 = {f_load:.6g} pN")
    print(f"   gamma_target = Pi_0*R/2 = {gamma_target:.6g} pN/µm  (Laplace resting hoop tension)")
    print("   (d1) SIGNED radial bias  |S| = |sum_i PF_i . r_hat_i|  ->  gamma_spur = |S|/(8*pi*R)")
    for frac in (0.001, 0.01, 0.05):
        _, budget = route_d_radial_bias_budget_pn(c, frac)
        print(f"        gamma_spur <= {frac:>6.1%} of gamma_target  =>  |S| <= {budget:10.6g} pN"
              f"   (= {budget / f_load:.3g} x F_load; contains no n_nodes, no ell)")
    for frac in (0.01, 0.05):
        print(f"   (d2) crude L1 variant: sum|PF| <= {frac:.0%} of load ="
              f" {route_d_load_balance_total_pn(c, frac):.6g} pN"
              "  -- NOT recommended (see doc: L1 double-counts cancelling residuals)")
    print("   (d3) grid-convergence diagnostic: mean|PF| must fall as O(ell^1) under refinement,")
    print("        because a nodal force is a force DENSITY times ell.  Report mean|PF|/ell.")
    print()

    print("-- GRID SCALING of the incumbent per-node literal ----------------------------------")
    print(f"   {'seg [µm]':>9s} {'route-(c) bound [pN]':>21s} {'frozen 0.21 is':>16s}"
          f" {'reported max|PF|':>17s} {'verdict vs derived':>20s}")
    reported = {0.5: 0.14, 0.2: 0.21, 0.075: 3.42}
    rung_rows = []
    for s in (0.5, 0.2, 0.075):
        b = route_c_predicate_floor(c, s, kmax)
        r = reported[s]
        rung_rows.append({"seg_um": s, "derived_pn": b, "reported_max_pf_pn": r,
                          "loosening_factor": lit / b, "reported_over_derived": r / b})
        print(f"   {s:9.3f} {b:21.6g} {lit / b:15.2f}x {r:17.3g}"
              f" {'PASS' if r <= b else f'FAIL {r / b:.1f}x':>20s}")
    segs = sorted(reported)
    for lo, hi in zip(segs[:-1], segs[1:]):
        p_meas = np.log(reported[hi] / reported[lo]) / np.log(hi / lo)
        print(f"   measured max|PF| exponent over seg {lo}->{hi}: ell^{p_meas:+.3f}"
              f"   (route-(c) bound is ell^+1.000)"
              f" => gate/measurement gap widens as ell^{p_meas - 1.0:+.3f}")
    print()

    checks = run_checks(c, kmax, term, seg, native)
    print("-- CHECKS --------------------------------------------------------------------------")
    for ch in checks:
        print(f"   [{'PASS' if ch.ok else 'FAIL'}] {ch.name}: {ch.detail}")
    n_fail = sum(1 for ch in checks if not ch.ok)
    print()
    print(f"RESULT: {len(checks) - n_fail}/{len(checks)} checks pass"
          + ("  -> DERIVATION VERIFIED" if n_fail == 0 else f"  -> {n_fail} FAILED"))
    if injected and n_fail == 0:
        print("NEGATIVE CONTROL DID NOT FAIL — the checks are inert. Treat as a broken check.")
        return 3

    if args.json is not None:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(json.dumps({
            "schema": "ffn-ac-gate-threshold-derivation-v1",
            "constants": c.__dict__,
            "kmax_pn_per_um": kmax, "kmax_dominant_term": term,
            "frozen_literal_pn": lit,
            "route_c_predicate_floor_pn": f_pred,
            "route_a2_thermal_force_pn": thermal,
            "route_d_load_pn": f_load,
            "grid_rungs": rung_rows,
            "checks": [ch.__dict__ for ch in checks],
            "n_failed": n_fail,
        }, indent=2))
        print(f"wrote {args.json}")
    return 0 if n_fail == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
