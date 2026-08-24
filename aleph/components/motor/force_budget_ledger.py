r"""Rigorous cortical-tension force-budget ledger for the assembled Active Cell (F10, done right).

WHY THIS EXISTS
---------------
A prior consult reduced the assembled cell's myosin cortical tension to a single hand-wavy line —
"cortical tension ~210x below band" — from ONE pessimistic number: all 442 head-resolved
Stam-Hocky minifilaments of the saved R=7.5 um run treated as crossing a single diametral cut
plane, giving gamma = 4.72e-5 N/m, then divided into a ~1e-2 N/m MCF7 whole-cell cortical tension
(Moazzeni; project memory ``project-mcf7-parameter-collection``) -> ~212x.

That number is loose on BOTH ends:

  * it is a LOOSE UPPER BOUND on the method-of-planes tension (every minifilament forced onto one
    plane with full projection ``|u.n|=1``), and
  * it is normalized against the WRONG target. 1e-2 N/m (= 10 mN/m) is a total apparent whole-cell
    tension. The myosin-GENERATED cortex tension must be compared to the ACTIVE FRACTION of the
    cortical-tension band (``ff.gamma_estimator.active_band_pn_um`` ~ 0.245-0.455 mN/m; the
    FF_STAGE6P hard rule that ``gamma_active`` is compared to ``ACTIVE_FRACTION * band``), NOT the
    full total band and certainly not a 10 mN/m whole-cell value.

This module recomputes the budget THREE ways with correct normalization and sweeps the I0-B3
PI-GAP magnitudes, so the ~210x decomposes cleanly into (wrong normalization) x (real below-band
deficit), and the real deficit into (under-population) x (per-head budget). No constant is tuned
to a target; every value is sourced (see the CONSTANTS block).

THE THREE SCENARIOS
-------------------
(i)   LOOSE BOUND -- reproduces the prior number. Method-of-planes UPPER bound: every minifilament
      is assumed to cross the SAME diametral cut with full projection, each carrying its one-sided
      tension ``F_side = N_side * phi * f_stall``::

          gamma_i = n_mf * F_side / (2*pi*R_cortex)

      This overestimates: a real diametral cut through an isotropically distributed population is
      crossed by only a fraction of it, each element contributing ``T*|u.n|`` with ``<|u.n|> < 1``
      (the isotropic Kirkwood factor lowers it ~100x).

(ii)  DISTRIBUTED active (force-dipole) stress -- the areal average of the *active* stress a cut
      would release::

          sigma_a * h = (n_heads / A) * (phi * f_stall * d_step)

      Each engaged head is a contractile force-dipole of moment ``(phi*f_stall) * d_step``. This is
      far below (i) because the working stroke ``d_step`` (5-10 nm) is a minuscule lever versus the
      cell circumference; it isolates the purely *active* (power-stroke) contribution from the
      sustained-load bound of (i). NOTE: reading ``d`` as the element/overlap length (~L_bb=301 nm)
      instead of the step -- the isotropic sustained-tension Kirkwood stress -- raises (ii) by
      ~``L_bb/d_step`` (~40x) and recovers ~1% of (i); that variant is DOCUMENTED, never tuned in.

(iii) PHYSIOLOGICAL DENSITY -- re-evaluate (i) at the prior FF/KB physiological engaged-NMII areal
      density 16-21 minifilament/um^2 (gamma-floor layered-resolution, 2026-06-09) versus the saved
      run's ``442/A = 0.625/um^2``. Isolates how much of the deficit is UNDER-POPULATION.

HEADLINE (what the sweep proves — and what it does NOT)
-------------------------------------------------------
The prior ~212x = ~30x wrong-normalization (1e-2 N/m whole-cell vs ~0.33 mN/m active-cortex band)
times ~7x real below-band deficit *on the LOOSE bound*. This module reconciles the ~212x; it does
NOT establish that physiological density closes the band, because scenario (i) is an explicit UPPER
bound. ⚠ CORRECTION (do not resurrect the earlier optimism): the physical force-dipole / method-of-
planes estimator (``two_filament_reference.gamma_ceiling`` = rho*F*L_bb/2, the convention the NG-1
gate's reference uses) gives ~100x LESS than the loose bound at the same density, and stays
16-38x BELOW the band at 16-21 /um^2 (claim_a needs rho ~= 469 /um^2). So the loose-bound "2.7-3.5x
above band at physiological density" is NOT the physical answer; the two estimators BRACKET the true
gamma ~100x apart (isotropic-Kirkwood / network-propagation gap). Additionally, the 16-21 /um^2
"physiological" density is itself a back-calculated GAP target (gamma-floor layered-resolution), not
a measurement -- the only DIRECT cortical NMII density is Nie 2015 = 0.625 /um^2. Under-population is
a large lever, but "physiological density closes the band" is UNRESOLVED: the real gamma requires the
assembled-cell method-of-planes with the FIXED motor (Round-1 measured the broken motor). See
``tests/ac/motor/test_density_ledger_bridge.py`` for the honest bracket gate.

CONSTANTS -- every value sourced (HARD RULE: no empirical magic numbers, none tuned to a band)
------------------------------------------------------------------------------------------------
  R_CELL_UM / CORTEX_MEMBRANE_GAP_UM / R_CORTEX_UM  ac/cell/assemble.py:97,106,107
  CURRENT_N_MF=442, CURRENT_N_HEADS=8840            saved assembled run (442 mf * 2*N_side heads)
  N_side  claim_a=10 / claim_b=28                   ac/motor/params_i0b3.yaml:68,72
  f_stall claim_a=0.5 / claim_b=2.0 pN              ac/motor/params_i0b3.yaml:50,54
  d_step 5-10 nm                                    ac/motor/params_i0b3.yaml:305 (reference_band)
  k_on=50/s, k_off0=0.35/s, kBT, x_beta -> f0       ac/motor/params_i0b3.yaml:238,253,294,264
  PHYS_DENSITY_BAND 16-21 /um^2                     gamma-floor layered-resolution (2026-06-09)
  active band (245-455 pN/um)                       ff.gamma_estimator.active_band_pn_um (REUSED)
  PRIOR_TARGET_N_PER_M 1e-2                          Moazzeni MCF7 whole-cell tension (project memory)
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

from aleph.laws.gamma_estimator import (
    ACTIVE_FRACTION,
    MCF7_IQR_PN_UM,
    SALBREUX_BAND_PN_UM,
    active_band_pn_um,
    gamma_to_N_per_m,
)

__all__ = [
    "engaged_fraction",
    "per_side_force_pn",
    "gamma_loose_bound",
    "gamma_distributed",
    "gamma_at_density",
    "BandPosition",
    "band_position",
    "ScenarioResult",
    "SweepRow",
    "sweep",
    "ForceBudgetLedger",
    "build_ledger",
]

# ── geometry (ac/cell/assemble.py) ───────────────────────────────────────────────────────────────────
R_CELL_UM = 7.5                      # cell OUTER (membrane) radius [um] — assemble.py:97
CORTEX_MEMBRANE_GAP_UM = 0.10        # membrane->cortex-shell offset [um] — assemble.py:106
R_CORTEX_UM = R_CELL_UM - CORTEX_MEMBRANE_GAP_UM   # 7.40 um — assemble.py:107; where minifilaments sit
CELL_AREA_UM2 = 4.0 * math.pi * R_CELL_UM**2       # 706.86 um^2 — sets the count<->areal-density map
# The method-of-planes cut circumference uses the CORTEX-shell radius (7.40 um): minifilaments are seeded
# at cortex myosin-link midpoints on that shell (assemble.py:239-241), so the great circle a diametral cut
# actually traverses has radius R_CORTEX_UM, not R_CELL_UM. (Using R_CELL_UM gives 4.65e-5 vs 4.72e-5 — a
# 1.4% geometric difference, immaterial to the >=5x below-band conclusion; exposed as CELL_CIRC_UM.)
CORTEX_CIRC_UM = 2.0 * math.pi * R_CORTEX_UM       # 46.496 um — the physical diametral cut length
CELL_CIRC_UM = 2.0 * math.pi * R_CELL_UM           # 47.124 um — nominal-sphere variant (for the note)

# ── saved assembled run (the state the prior consult scored) ────────────────────────────────────────
CURRENT_N_MF = 442                                 # minifilaments on the R=7.5 um sphere (saved run)
CURRENT_N_HEADS = 8840                             # = 442 * 2 * N_side(=10) explicit heads
CURRENT_DENSITY_PER_UM2 = CURRENT_N_MF / CELL_AREA_UM2   # 0.6252 minifilament/um^2

# ── I0-B3 PI-GAP magnitudes (ac/motor/params_i0b3.yaml) — CONFLICTED claims, swept as oracle variables ─
N_SIDE_CLAIM = {"a": 10, "b": 28}                  # heads/side: AFINES(10) vs Billington(28) — yaml:68,72
F_STALL_CLAIM_PN = {"a": 0.5, "b": 2.0}            # per-head stall [pN]: AFINES(0.5) vs Billington(2.0) — yaml:50,54
D_STEP_BAND_UM = (0.005, 0.010)                    # working-stroke 5-10 nm — yaml:305 reference_band
K_ON = 50.0                                        # per-head attach rate [1/s] — yaml:238
K_OFF0 = 0.35                                      # Bell zero-force detach prefactor [1/s] — yaml:253
KBT_PN_UM = 4.28e-3                                # thermal energy at 37 C [pN*um] — yaml:294
X_BETA_UM = 0.6e-3                                 # Bell bond length (0.6 nm) — yaml:264
F0_PN = KBT_PN_UM / X_BETA_UM                      # Bell characteristic force ~7.13 pN — yaml:271

# ── physiological engaged-NMII areal density (gamma-floor layered-resolution, 2026-06-09) ───────────
PHYS_DENSITY_BAND = (16.0, 21.0)                   # minifilament/um^2 needed to lift the native gamma-floor

# ── prior mis-applied normalization target (documented, NOT the correct comparison) ─────────────────
PRIOR_TARGET_N_PER_M = 1.0e-2                      # Moazzeni MCF7 whole-cell cortical tension [N/m]


def engaged_fraction(f_stall_pn: float, *, k_on: float = K_ON, k_off0: float = K_OFF0,
                     f0_pn: float = F0_PN) -> float:
    """Bell-loaded bound-head fraction ``phi`` when a head bears the stall load ``f_stall``.

    Mechanistic (Bell-Evans, per the project hard rule): the detachment rate at load ``f`` is
    ``k_off = k_off0 * exp(f/f0)`` and the steady bound fraction is ``phi = k_on/(k_on+k_off)``.
    Evaluated at ``f = f_stall`` this gives 0.9925 at claim_a (0.5 pN) — the value the prior consult
    used — and 0.9908 at claim_b (2.0 pN); it is EMERGENT from k_on/k_off, never an imposed duty.

    Args:
        f_stall_pn: per-head isometric stall force [pN].
        k_on: per-head actin attachment rate [1/s].
        k_off0: Bell zero-force detachment prefactor [1/s].
        f0_pn: Bell characteristic force ``kBT/x_beta`` [pN].

    Returns:
        The engaged (bound) head fraction ``phi`` in ``(0, 1]``.
    """
    k_off = k_off0 * math.exp(f_stall_pn / f0_pn)
    return k_on / (k_on + k_off)


def per_side_force_pn(n_side: int, f_stall_pn: float, phi: float | None = None) -> float:
    """One-sided minifilament tension ``F_side = N_side * phi * f_stall`` [pN].

    A bipolar minifilament pulls two antiparallel actin filaments together; the tension transmitted
    across a cut between them equals the pull of ONE side's ``N_side`` heads (each at ``phi*f_stall``).

    Args:
        n_side: explicit motor heads per anti-parallel half-filament.
        f_stall_pn: per-head isometric stall force [pN].
        phi: engaged fraction; defaults to the Bell-loaded value at ``f_stall``.

    Returns:
        The one-sided contractile tension per minifilament [pN].
    """
    if phi is None:
        phi = engaged_fraction(f_stall_pn)
    return n_side * phi * f_stall_pn


def gamma_loose_bound(n_mf: float, n_side: int, f_stall_pn: float, *,
                      circ_um: float = CORTEX_CIRC_UM, phi: float | None = None) -> tuple[float, float]:
    """Scenario (i): the method-of-planes LOOSE UPPER BOUND on cortical tension [N/m, pN/um].

    Assumes every minifilament crosses the SAME diametral cut with full projection::

        gamma = n_mf * (N_side * phi * f_stall) / (2*pi*R_cortex)

    Overestimates the true (distributed) value because only a fraction of an isotropic population
    crosses any one plane and each contributes ``T*|u.n|`` with ``<|u.n|> < 1``.

    Args:
        n_mf: number of minifilaments (may be fractional = density*area).
        n_side: heads per anti-parallel side.
        f_stall_pn: per-head stall force [pN].
        circ_um: cut circumference [um] (default = the cortex-shell great circle).
        phi: engaged fraction; defaults to the Bell-loaded value at ``f_stall``.

    Returns:
        ``(gamma_N_per_m, gamma_pN_per_um)``.
    """
    f_side = per_side_force_pn(n_side, f_stall_pn, phi)
    gamma_pn_um = n_mf * f_side / circ_um
    return gamma_to_N_per_m(gamma_pn_um), gamma_pn_um


def gamma_distributed(n_heads: float, f_stall_pn: float, d_step_um: float, *,
                      area_um2: float = CELL_AREA_UM2, phi: float | None = None) -> tuple[float, float]:
    """Scenario (ii): the DISTRIBUTED active (force-dipole) stress ``sigma_a*h`` [N/m, pN/um].

    Areal density of contractile force-dipoles, each of moment ``(phi*f_stall)*d_step``::

        sigma_a * h = (n_heads / A) * (phi * f_stall * d_step)

    This is the purely *active* (power-stroke) contribution; it is far below scenario (i) because the
    working stroke ``d_step`` (5-10 nm) is a tiny lever versus the cell (see module docstring for the
    element-length Kirkwood variant).

    Args:
        n_heads: number of explicit myosin heads.
        f_stall_pn: per-head stall force [pN].
        d_step_um: myosin working-stroke displacement [um].
        area_um2: cortex surface area [um^2].
        phi: engaged fraction; defaults to the Bell-loaded value at ``f_stall``.

    Returns:
        ``(gamma_N_per_m, gamma_pN_per_um)``.
    """
    if phi is None:
        phi = engaged_fraction(f_stall_pn)
    gamma_pn_um = (n_heads / area_um2) * (phi * f_stall_pn * d_step_um)
    return gamma_to_N_per_m(gamma_pn_um), gamma_pn_um


def gamma_at_density(density_per_um2: float, n_side: int, f_stall_pn: float, *,
                     area_um2: float = CELL_AREA_UM2, circ_um: float = CORTEX_CIRC_UM) -> tuple[float, float]:
    """Scenario (iii): the loose-bound gamma [N/m, pN/um] at a chosen minifilament areal density.

    Sets ``n_mf = density * A`` and evaluates :func:`gamma_loose_bound`. Used to re-price the budget
    at the physiological engaged-NMII density (16-21/um^2) versus the saved run's 0.625/um^2.

    Args:
        density_per_um2: minifilament areal density [1/um^2].
        n_side: heads per anti-parallel side.
        f_stall_pn: per-head stall force [pN].
        area_um2: cortex surface area [um^2] (count<->density map).
        circ_um: cut circumference [um].

    Returns:
        ``(gamma_N_per_m, gamma_pN_per_um)``.
    """
    n_mf = density_per_um2 * area_um2
    return gamma_loose_bound(n_mf, n_side, f_stall_pn, circ_um=circ_um)


@dataclass(frozen=True)
class BandPosition:
    """Where a gamma sits relative to a two-edge literature band (all comparisons in N/m).

    Attributes:
        label: ``"below"``, ``"in-band"`` or ``"above"``.
        factor_lo: band-lower / gamma  (``>1`` => gamma is below the lower edge).
        factor_hi: band-upper / gamma  (``>1`` => gamma is below the upper edge).
    """

    label: str
    factor_lo: float
    factor_hi: float

    def summary(self) -> str:
        """One-line human summary, e.g. ``'5.2-9.6x below'`` or ``'2.7-3.5x above'``."""
        if self.label == "below":            # factors are band/gamma > 1; near edge = lower factor
            return f"{self.factor_lo:.1f}-{self.factor_hi:.1f}x below"
        if self.label == "above":            # gamma/band = 1/factor; near edge = upper
            return f"{1.0 / self.factor_hi:.1f}-{1.0 / self.factor_lo:.1f}x above"
        return "in-band"


def band_position(gamma_N_per_m: float, band_pn_um: tuple[float, float]) -> BandPosition:
    """Classify a gamma against a literature band given in pN/um.

    Args:
        gamma_N_per_m: the tension to classify [N/m].
        band_pn_um: ``(lower, upper)`` band edges [pN/um].

    Returns:
        A :class:`BandPosition`.
    """
    lo = gamma_to_N_per_m(band_pn_um[0])
    hi = gamma_to_N_per_m(band_pn_um[1])
    factor_lo = lo / gamma_N_per_m
    factor_hi = hi / gamma_N_per_m
    if gamma_N_per_m < lo:
        label = "below"
    elif gamma_N_per_m > hi:
        label = "above"
    else:
        label = "in-band"
    return BandPosition(label=label, factor_lo=factor_lo, factor_hi=factor_hi)


@dataclass(frozen=True)
class ScenarioResult:
    """A single scenario's gamma plus its position versus the active band.

    Attributes:
        name: scenario label.
        gamma_N_per_m: cortical tension [N/m].
        gamma_pN_per_um: same tension in engine/gamma_estimator units [pN/um].
        band: position versus the myosin-active band.
        detail: free-form provenance / assumption string.
    """

    name: str
    gamma_N_per_m: float
    gamma_pN_per_um: float
    band: BandPosition
    detail: str = ""


@dataclass(frozen=True)
class SweepRow:
    """One ``(N_side, f_stall, density)`` combination of the loose-bound sensitivity sweep."""

    n_side: int
    f_stall_pn: float
    density_per_um2: float
    n_mf: float
    phi: float
    gamma_N_per_m: float
    gamma_pN_per_um: float
    band: BandPosition


def sweep(n_sides: tuple[int, ...] = (10, 28),
          f_stalls: tuple[float, ...] = (0.5, 2.0),
          densities: tuple[float, ...] = (CURRENT_DENSITY_PER_UM2, 16.0, 21.0), *,
          band_pn_um: tuple[float, float] | None = None) -> list[SweepRow]:
    """Sweep the loose-bound (scenario i) gamma over the I0-B3 GAP magnitudes and density.

    Args:
        n_sides: heads-per-side values (claim_a=10, claim_b=28).
        f_stalls: per-head stall values [pN] (claim_a=0.5, claim_b=2.0).
        densities: minifilament areal densities [1/um^2] (current, then physiological 16 & 21).
        band_pn_um: band to score against [pN/um]; defaults to the myosin-active band.

    Returns:
        A list of :class:`SweepRow`, ordered density-major then ``(N_side, f_stall)``.
    """
    if band_pn_um is None:
        band_pn_um = active_band_pn_um()
    rows: list[SweepRow] = []
    for dens in densities:
        n_mf = dens * CELL_AREA_UM2
        for ns in n_sides:
            for fs in f_stalls:
                phi = engaged_fraction(fs)
                g_n, g_pn = gamma_at_density(dens, ns, fs)
                rows.append(SweepRow(
                    n_side=ns, f_stall_pn=fs, density_per_um2=dens, n_mf=n_mf, phi=phi,
                    gamma_N_per_m=g_n, gamma_pN_per_um=g_pn, band=band_position(g_n, band_pn_um)))
    return rows


def format_sweep_table(rows: list[SweepRow]) -> str:
    """Render the sweep as a fixed-width text table (SI + engine units + band position)."""
    head = (f"{'N_side':>6} {'f_stall':>8} {'density':>9} {'n_mf':>8} "
            f"{'gamma[N/m]':>12} {'gamma[pN/um]':>13} {'vs active band':>18}")
    sep = "-" * len(head)
    lines = [head, sep]
    for r in rows:
        lines.append(f"{r.n_side:>6d} {r.f_stall_pn:>7.1f}p {r.density_per_um2:>8.3f} {r.n_mf:>8.0f} "
                     f"{r.gamma_N_per_m:>12.3e} {r.gamma_pN_per_um:>13.4f} {r.band.summary():>18}")
    return "\n".join(lines)


@dataclass
class ForceBudgetLedger:
    """The full three-scenario budget + sensitivity sweep + the ~210x decomposition.

    Attributes:
        loose: scenario (i) at claim_a (reproduces the prior 4.72e-5 N/m).
        distributed: scenario (ii) at claim_a, d_step midpoint (with band-edge variants in ``extras``).
        physiological: scenario (iii) loose bound at the physiological density band, claim_a.
        rows: the 12-row ``(N_side, f_stall, density)`` loose-bound sweep.
        decomposition: the multiplicative split of the prior ~210x.
        extras: auxiliary numbers used by the figure / report.
    """

    loose: ScenarioResult
    distributed: ScenarioResult
    physiological: list[ScenarioResult]
    rows: list[SweepRow]
    decomposition: dict = field(default_factory=dict)
    extras: dict = field(default_factory=dict)

    def format_report(self) -> str:
        """A full text report: scenarios, the ~210x decomposition, and the sweep table."""
        act = active_band_pn_um()
        L: list[str] = []
        L.append("=" * 78)
        L.append("CORTICAL-TENSION FORCE BUDGET — assembled Active Cell (F10, done right)")
        L.append("=" * 78)
        L.append(f"saved run: {CURRENT_N_MF} minifilaments / {CURRENT_N_HEADS} heads on R={R_CELL_UM} um "
                 f"(density {CURRENT_DENSITY_PER_UM2:.3f} /um^2)")
        L.append(f"active myosin band (REUSED ff.gamma_estimator): {act[0]:.0f}-{act[1]:.0f} pN/um "
                 f"= {gamma_to_N_per_m(act[0]):.2e}-{gamma_to_N_per_m(act[1]):.2e} N/m "
                 f"(= ACTIVE_FRACTION {ACTIVE_FRACTION} * Salbreux/Chugh total band)")
        L.append(f"total Salbreux/Chugh band: {SALBREUX_BAND_PN_UM[0]:.0f}-{SALBREUX_BAND_PN_UM[1]:.0f} pN/um; "
                 f"MCF7 IQR: {MCF7_IQR_PN_UM[0]:.0f}-{MCF7_IQR_PN_UM[1]:.0f} pN/um")
        L.append("")
        for s in [self.loose, self.distributed, *self.physiological]:
            L.append(f"  {s.name:<34} {s.gamma_N_per_m:.3e} N/m  ({s.gamma_pN_per_um:.4g} pN/um)  "
                     f"{s.band.summary():>14} active band")
            if s.detail:
                L.append(f"      {s.detail}")
        L.append("")
        d = self.decomposition
        L.append("--- decomposition of the prior ~210x ---")
        L.append(f"  prior ratio (loose / 1e-2 N/m Moazzeni whole-cell) : {d['prior_ratio']:.1f}x")
        L.append(f"    = wrong-normalization {d['normalization_factor']:.1f}x "
                 f"(1e-2 N/m whole-cell vs active-band geomean {d['active_geomean_N_per_m']:.2e} N/m)")
        L.append(f"    x real below-band deficit {d['residual_below_band']:.1f}x "
                 f"(loose bound vs active-band geomean)")
        L.append(f"  under-population headroom (16-21 vs 0.625 /um^2) : {d['density_headroom']:.1f}x "
                 f"(NOTE: 16-21 is a back-calc GAP target, not measured; direct = Nie 0.625 /um^2)")
        L.append(f"  per-head-budget headroom (claim_a -> claim_b)    : {d['perhead_headroom']:.1f}x")
        L.append(f"  ⚠ all quantities above are on the LOOSE UPPER bound. The physical force-dipole /")
        L.append(f"    method-of-planes estimator (two_filament_reference.gamma_ceiling) is ~100x lower and")
        L.append(f"    stays 16-38x BELOW band at 16-21 /um^2 (claim_a needs rho~469). BAND-CROSSING at")
        L.append(f"    physiological density is UNRESOLVED (the two estimators bracket gamma ~100x apart);")
        L.append(f"    it needs the assembled-cell method-of-planes with the fixed motor, NOT this loose bound.")
        L.append("")
        L.append("--- loose-bound (UPPER) sensitivity sweep (N_side x f_stall x density) ---")
        L.append(format_sweep_table(self.rows))
        return "\n".join(L)


def build_ledger() -> ForceBudgetLedger:
    """Assemble the full force-budget ledger from the sourced constants.

    Returns:
        A :class:`ForceBudgetLedger` with the three scenarios, the 12-row loose-bound sweep, and the
        multiplicative decomposition of the prior ~210x.
    """
    act = active_band_pn_um()
    ns_a, fs_a = N_SIDE_CLAIM["a"], F_STALL_CLAIM_PN["a"]

    # (i) LOOSE BOUND at claim_a — reproduces the prior 4.72e-5 N/m.
    g_i_n, g_i_pn = gamma_loose_bound(CURRENT_N_MF, ns_a, fs_a)
    loose = ScenarioResult(
        name="(i) loose bound [claim_a]", gamma_N_per_m=g_i_n, gamma_pN_per_um=g_i_pn,
        band=band_position(g_i_n, act),
        detail=(f"all {CURRENT_N_MF} minifilaments on one cut, F_side={per_side_force_pn(ns_a, fs_a):.3f} pN; "
                f"nominal-R variant {gamma_loose_bound(CURRENT_N_MF, ns_a, fs_a, circ_um=CELL_CIRC_UM)[0]:.3e} N/m"))

    # (ii) DISTRIBUTED force-dipole stress at claim_a, d_step midpoint (+ band-edge variants).
    d_mid = 0.5 * (D_STEP_BAND_UM[0] + D_STEP_BAND_UM[1])
    g_ii_n, g_ii_pn = gamma_distributed(CURRENT_N_HEADS, fs_a, d_mid)
    g_ii_lo = gamma_distributed(CURRENT_N_HEADS, fs_a, D_STEP_BAND_UM[0])[0]
    g_ii_hi = gamma_distributed(CURRENT_N_HEADS, fs_a, D_STEP_BAND_UM[1])[0]
    distributed = ScenarioResult(
        name="(ii) distributed dipole [claim_a]", gamma_N_per_m=g_ii_n, gamma_pN_per_um=g_ii_pn,
        band=band_position(g_ii_n, act),
        detail=(f"(n_heads/A)*(phi*f*d_step), d_step={d_mid * 1e3:.1f} nm "
                f"(band {g_ii_lo:.2e}-{g_ii_hi:.2e} N/m for 5-10 nm)"))

    # (iii) PHYSIOLOGICAL DENSITY loose bound at claim_a, both band edges.
    physiological: list[ScenarioResult] = []
    for dens in PHYS_DENSITY_BAND:
        g_n, g_pn = gamma_at_density(dens, ns_a, fs_a)
        physiological.append(ScenarioResult(
            name=f"(iii) loose @ {dens:.0f}/um^2 [claim_a]", gamma_N_per_m=g_n, gamma_pN_per_um=g_pn,
            band=band_position(g_n, act),
            detail=f"physiological engaged-NMII density (gamma-floor layered-resolution 2026-06-09)"))

    rows = sweep(band_pn_um=act)

    # --- decomposition of the prior ~210x ---
    active_geomean = math.sqrt(gamma_to_N_per_m(act[0]) * gamma_to_N_per_m(act[1]))
    prior_ratio = PRIOR_TARGET_N_PER_M / g_i_n
    normalization_factor = PRIOR_TARGET_N_PER_M / active_geomean
    residual_below_band = active_geomean / g_i_n
    phys_geomean_density = math.sqrt(PHYS_DENSITY_BAND[0] * PHYS_DENSITY_BAND[1])
    density_headroom = phys_geomean_density / CURRENT_DENSITY_PER_UM2
    perhead_headroom = ((N_SIDE_CLAIM["b"] / N_SIDE_CLAIM["a"])
                        * (F_STALL_CLAIM_PN["b"] / F_STALL_CLAIM_PN["a"])
                        * (engaged_fraction(F_STALL_CLAIM_PN["b"]) / engaged_fraction(F_STALL_CLAIM_PN["a"])))
    decomposition = {
        "prior_ratio": prior_ratio,
        "normalization_factor": normalization_factor,
        "residual_below_band": residual_below_band,
        "active_geomean_N_per_m": active_geomean,
        "density_headroom": density_headroom,
        "perhead_headroom": perhead_headroom,
        "phys_geomean_density": phys_geomean_density,
    }
    extras = {
        "active_band_pn_um": act,
        "active_band_N_per_m": (gamma_to_N_per_m(act[0]), gamma_to_N_per_m(act[1])),
        "salbreux_band_N_per_m": (gamma_to_N_per_m(SALBREUX_BAND_PN_UM[0]), gamma_to_N_per_m(SALBREUX_BAND_PN_UM[1])),
        "mcf7_iqr_N_per_m": (gamma_to_N_per_m(MCF7_IQR_PN_UM[0]), gamma_to_N_per_m(MCF7_IQR_PN_UM[1])),
        "prior_target_N_per_m": PRIOR_TARGET_N_PER_M,
        "distributed_band_N_per_m": (g_ii_lo, g_ii_hi),
    }
    return ForceBudgetLedger(loose=loose, distributed=distributed, physiological=physiological,
                             rows=rows, decomposition=decomposition, extras=extras)


if __name__ == "__main__":  # pragma: no cover — convenience console dump
    print(build_ledger().format_report())
