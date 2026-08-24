"""γ-floor prototype — MD-free cortical tension from mechanical equilibrium (Stage 6d, ENGINE.md §4).

The decisive FF experiment. γ (cortical tension) is **quasi-static** (the load barely changes during
its measurement), so even the load-coupled myosin has a well-defined steady distribution that can be
sampled once (quenched). Protocol (ENGINE.md §4):

  1. MC-sample a realization of the cortex fiber network + crosslinker + myosin bound state at init.
  2. Solve the **mechanical equilibrium** (bending + crosslinker springs + myosin active prestress)
     with the inextensibility projector + reshape — NO explicit BAOAB thermostat (regime A/B).
  3. Measure γ via the method-of-planes over all load-bearing elements (actin axial tension +
     crosslinker links + myosin links).
  4. Repeat over N realizations (embarrassingly parallel quenched ensemble) → a γ distribution.
  5. **Decisive validation:** compare this MD-free γ to the BAOAB-MD γ from the archived
     ``cortical_tension.py`` (the γ-floor finding) — does an MD-free mechanical solve reproduce the
     same γ (⇒ FF is MD-equivalent for γ), or does it localize where in-time load-coupling matters?

γ is **swept as a controlled variable** over the myosin prestress f_myo (NOT tuned to a band — hard
rule). All constants are lit-anchored (ff.units / cortex_assembly / hand_kmc presets).

The myosin active prestress: a non-muscle myosin IIA minifilament between two actin points walks to
the + ends and stalls against the network at mechanical equilibrium (v=0), exerting a contractile
force ≈ n_engaged_heads · F_stall_per_head (NMIIA F_stall=0.5 pN/head, n_heads_per_side=10 → tens of
pN). That stall force is the per-link prestress f_myo; we sweep it and mark the lit-anchored value.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from aleph.laws.turgor_pi0 import (
    PI0_CLAIM_HELA_PROXY_40PA,
    TurgorPi0Unspecified,          # noqa: F401  (re-exported for callers that catch the gate)
    resolve_pi0_pa,
)
from aleph.laws import units as U
from aleph.laws.constraints import reshape, segment_axial_tension
from aleph.laws.cortex_assembly import CortexParams, build_cortex_network
from aleph.laws.fiber_network import FiberNetwork
from aleph.laws.forces_warp import make_bending_force_fn
from aleph.laws.gamma_estimator import method_of_planes_gamma
from aleph.laws.hand_kmc import ALPHA_ACTININ, FILAMIN, NMIIA_MYOSIN


#: ⚠ RETIRED-CORTEX-BUILD BANNER (audit R1, 2026-07-25) — travels with every number this module emits.
_CORTEX_BUILD_BANNER: str = (
    "measured on the pre-2026-07-23 cortex build; not re-run; see "
    "docs/v2_audit/CORTEX_STRUCTURE_AUDIT_2026-07-23.md. THREE of the five 2026-07-23 structural "
    "defects are STILL LIVE on this path: (1) build_crosslinked_cortex takes n_xl from the caller and "
    "never reads ff/architecture_spec.CORTEX, and every production call site passes n_xl = n_filaments "
    "(density_per_fil = 1.0, the under-percolating relic replaced by 20 in bc5ff3b0); "
    "(2) resolve_overlaps is never passed to build_cortex_network, so it runs at the default False; "
    "(3) length_dist stays 'mono'. bc5ff3b0's own message names under-percolation as 'the upstream "
    "root of ... no myosin transmission' — and the gamma-floor measures EXACTLY myosin transmission."
)

# NMIIA per-side contractile force — the prestress, a FORCE-MAGNITUDE lever SWEPT as a controlled
# variable (gamma_floor_sweep; never tuned to band). F_stall_per_head=0.5 pN (Kovacs 2003). The
# per-side head count has a lit RANGE: the Stam-Hocky/AFINES coarse model = 10 heads/side (→ 5 pN, the
# conservative anchor below, configs/phase1_h3.yaml); the STRUCTURAL minifilament (Billington 2013 EM,
# ~29 molecules / ~58 heads) ≈ 29 heads/side (→ ~15 pN at full engagement); the ensemble-measured
# minifilament stall ≈ 17 pN (Stachowiak 2009). The engaged fraction (duty) sits between these. ⚠️ The
# floor is robust across this whole range: at the MCF7 cortex (R=7.5), γ_active vs the MCF7-active
# target is floored 216× (5 pN) → 64× (17 pN structural) → 36× (30 pN) — PERSISTS even at the structural
# upper end (FF_STAGE6Q). The 5 pN anchor is the conservative low end; the sweep spans to 30× it.
NMIIA_F_STALL_PER_HEAD = 0.5    # pN/head (Kovacs 2003)
NMIIA_HEADS_PER_SIDE = 10       # Stam-Hocky/AFINES coarse model; STRUCTURAL (Billington 2013) ≈ 29 (swept)
NMIIA_MINIFIL_STALL_PN = NMIIA_F_STALL_PER_HEAD * NMIIA_HEADS_PER_SIDE   # 5.0 pN per side (conservative anchor)

# ENGAGED (load-bearing overlap) FRACTION — the mechanistically-correct determinant of cortical tension
# (Truong-Quang et al. 2021, PMC8586027): cortical tension is set by the myosin-actin OVERLAP, NOT the
# total minifilament count. In low-tension INTERPHASE ~35% of cortical NMII lies OUTSIDE the actin cortex
# (bound by one end → does NOT transmit contractile stress → NON-engaged); this overhang → ~0 in high-
# tension MITOSIS (full overlap). So the force-bearing (engaged) fraction ≈ 0.65 interphase → 1.0 mitotic,
# and the interphase→mitosis ~3× tension rise is driven by engagement at ~constant myosin amount. FF's
# 2D cortical shell does not resolve the RADIAL overhang geometrically, so engagement enters as this
# measured scalar on the force-bearing count; because γ_myo is exactly LINEAR in the contributing-dipole
# count (verified, FF_STAGE6Q), the scalar is exact, not a fudge. ⚠️ This is a ~1.5× lever (0.65→1.0) —
# it makes the model faithful and DEEPENS the interphase floor, it does NOT close the ~30× gap.
ENGAGED_FRACTION_INTERPHASE = 0.65   # Truong-Quang 2021 (~35% overhang, resting/adherent baseline)
ENGAGED_FRACTION_MITOTIC = 1.0       # full actin-cortex overlap (rounded/mitotic)

# Physiological turgor — STATE-DEPENDENT osmotic closure (2026-06-30 turgor workflow; the band-implied
# 133 Pa was a tuned circular value). FF units: 1 Pa = 1 pN/µm². Animal cells (no wall) hold no static
# turgor; the net excess is the sub-mM osmotic difference balanced by/coupled to cortical tension
# (Stewart 2011; Kay & Blaustein 2019 pump-leak). Resting baseline re-anchored to the MEASURED
# interphase value (Fischer-Friedrich 2014 Sci Rep 4:6213, HeLa interphase ΔP=40±30 Pa — HeLa proxy,
# no MCF7-specific datum exists, flagged to PI). The volume response uses Guo et al. 2017 (PNAS
# 114:E8618) entropic excluded-volume closure Π(V)=N·kB·T/(V−Vmin), which DERIVES the bulk modulus
# (no magic K_vol) from lit-anchored inputs: c_osm≈200 mM cytoplasmic osmolytes (cross-validated vs
# cytoplasmic salt) and Vmin≈0.30·V0 (Venkova 2022 eLife / Adar 2025 Ponder fit). Young-Laplace
# ΔP=2γ/R is now an OUTPUT/consistency-check, not the input defining ΔP0 (correct causal direction).
# ⚠ Π₀ IS GATED (PI 2026-07-25) — aleph/laws/params_turgor.yaml holds `value: null`
# (`evidence_status: "GAP — PI"`). There is NO silent default any more: `turgor_pressure`,
# `external_force`, `measure_gamma`, `gamma_floor_run` and `gamma_floor_production` all take
# `dP0=None` and RAISE `TurgorPi0Unspecified` if turgor is applied without an explicit value.
# Turgor constants moved to `ff/turgor_constants.py` (2026-07-28) so the hot path need not
# import this module; re-exported here so this file and its callers are untouched.
from aleph.laws.turgor_constants import (  # noqa: E402,F401
    OSMOLYTE_C_MM, TURGOR_DP0, TURGOR_DP0_PROVENANCE, TURGOR_PI_IN0, VMIN_FRAC,
)

# GROUNDED production operating point (configs/phase1_h3.yaml — the ×40 mesoscale cell, NOT a
# prototype). n_fil=1000 (Plan v2 §3 H.3), n_xl=1000 (KU-3.19).
#
# ⚠️ CITATION FIX (2026-06-30): the myosin count's "Salbreux 2012, 3/µm²" attribution is a CONFIRMED
# MISATTRIBUTION (project audit docs/v2_audit/H7_CORTICAL_MYOSIN_DENSITY_DATUM_2026-06-07.md — no
# Salbreux paper states a per-µm² minifilament count). The ONLY direct measurement of cortical NMII
# minifilament AREAL DENSITY is Nie et al. 2015 Cytoskeleton 72(1):29-46 (PMID 25641802): ~0.625
# minifil/µm² (HeLa medial cortex, intensity-calibrated 1 focus≈1 minifilament; range 0.31-0.94;
# LOW-MEDIUM confidence; non-MCF7 — a breast/MCF7 value does NOT exist in the literature). DECISIVE:
# the measured ~0.6/µm² is ~26-35× BELOW the ~16-21/µm² the active-γ ceiling needs, so the real datum
# CONFIRMS the γ-floor, it does not rescue it (re-anchoring lowers γ further). n_myo here is the
# ×40-mesoscale count; density is the controlled variable in gamma_floor_density_sweep (NOT tuned).
# NATIVE-SCALE production (PI 2026-07-01): the ×40 mesoscale (N=1000) was a CPU/HOOMD hardware
# constraint — RETIRED now that FF is GPU-native (the A5000 runs native scale; γ_active is N-DEPENDENT).
# LIT-FAITHFUL DENSITIES (FF_STAGE6Q, 2026-07-01): every cortical density is now set to its DIRECT
# measured areal value at the MCF7 radius R=7.5µm (area 4πR²≈707µm²), per the physiological-baseline
# HARD rule — NOT a convenient ratio:
#   • actin  ~100 /µm² (KB-3.18 cortex composition) → N_fil = round(100·4πR²) = 70686. This SUPERSEDES
#     the CLAUDE.md "~38000 native" (=30/µm² at R=10, ~3× sparser than KB-3.18). ⚠️ CLAUDE.md↔KB
#     conflict — surfaced to PI as a citation/magic-number fix (the actin COUNT does not drive the
#     γ-floor; γ is myosin-bound — so this is a faithfulness fix, not a floor lever).
#   • myosin 0.625 /µm² (Nie 2015, the only direct cortical NMII minifilament datum) → N_myo = 442.
#     This DECOUPLES myosin from the old n_fil/10 ratio (which, at the raised actin count, would put
#     myosin at ~10/µm² = 16× Nie). Using the measured Nie density is physiological-baseline-correct
#     and DEEPENS the floor (Nie is sparse) — the honest consequence, consistent with the density sweep.
#   • crosslinkers 1:1 with actin (γ-irrelevant, FF_STAGE6M) → N_xl = N_fil.
# Net: γ_active drops vs the prior ratio-based preset and the floor is DEEPER (~10³× under the MCF7-active
# band), because the measured-sparse Nie myosin density is now the operating point. The open lever is
# unchanged: the missing MCF7-ADHERENT load-engaged NMII density datum (engaged density, not this count).
PROD_N_FIL = 70686   # actin ~100/µm² (KB-3.18) at MCF7 R=7.5µm = round(100·4π·7.5²)
PROD_N_XL = 70686    # 1:1 crosslinkers (γ-irrelevant)
PROD_N_MYO = 442     # NMIIA minifilaments at Nie 2015 0.625/µm² (= round(0.625·4π·7.5²)) — physiological baseline

# Nie et al. 2015 measured cortical NMII minifilament areal density [µm⁻²] — the only direct datum.
NIE2015_DENSITY_UM2 = 0.625
NIE2015_DENSITY_RANGE = (0.31, 0.94)

# `CrosslinkedCortex` moved to `ff/fiber_network.py` (2026-07-28); re-exported so this
# module and its callers are untouched.
from aleph.laws.fiber_network import CrosslinkedCortex  # noqa: E402,F401
from aleph.laws.fiber_network import mesoscale_reach  # noqa: E402,F401  (moved 2026-07-28)


# The two helpers below moved to `ff/fiber_network.py` (2026-07-28) so callers can reach them
# without importing this module. Re-exported under their private names because this file
# uses them throughout and the rename is not the point of that move.
from aleph.laws.fiber_network import cross_fiber_pairs as _cross_fiber_pairs  # noqa: E402
from aleph.laws.fiber_network import node_fiber_map as _node_fiber_map        # noqa: E402


def build_crosslinked_cortex(params: CortexParams | None = None, *, n_filaments: int = 100,
                             n_xl: int = 300, n_myo: int = 100,
                             alpha_fraction: float = 0.30,
                             orientation: str = "isotropic", nematic_S: float = 1.0,
                             length_dist: str = "mono",
                             rng: np.random.Generator | None = None) -> CrosslinkedCortex:
    """Assemble one quenched cortex realization: fibers + crosslinker links + myosin links.

    Crosslinkers are α-actinin (``alpha_fraction``, Stricker 2010 30 %) / filamin (rest); each links
    two near nodes on different fibers (capture radius = the crosslinker Hand's ε). Myosin links
    connect cross-fiber node pairs within the NMIIA capture radius. Counts scale with the H.3 cortex
    (n_xl≈n_filaments·... ; defaults sized for a fast prototype).

    ⚠⚠ **PRE-2026-07-23 CORTEX BUILD — three declared structural defects are STILL LIVE here**
    (verified 2026-07-25 against this function and its ~25 call sites; see
    ``docs/v2_audit/CORTEX_STRUCTURE_AUDIT_2026-07-23.md`` +
    ``docs/v2_audit/AUDIT_WHOLE_REPO_2026-07-25.md`` R1):

    1. ``n_xl`` is a CALLER ARGUMENT and this function never reads ``ff/architecture_spec.CORTEX``.
       Every in-tree call site passes ``n_xl = n_filaments`` (``PROD_N_XL == PROD_N_FIL == 70686``),
       i.e. crosslink density 1.0/filament — the ×40 relic that left the native network at mean
       degree 2 (giant component 79.6 %, 11,430 fragments). The 2026-07-23 fix (``bc5ff3b0``) raised
       ``density_per_fil`` 1.0 → 20 in ``architecture_spec.py`` ONLY; that commit's claim that
       "both weave and gamma_floor read the same spec" is **false for this function**.
    2. ``resolve_overlaps`` is never passed to :func:`build_cortex_network`, so the cortex is built at
       the default ``False`` (build interpenetrations retained; ``ff.weave.weave`` is the only path
       that passes it).
    3. ``length_dist`` defaults to ``"mono"`` (single filament length), not the exponential
       distribution.

    Everything measured through this builder — including the ~530× γ-floor headline — therefore
    carries :data:`_CORTEX_BUILD_BANNER` and must be re-run on the fixed cortex before being cited.
    Fixing them here changes the build and breaks γ-floor bit-parity, so it is a PI call, not a
    drive-by edit; the banner is the interim honest statement.
    """
    if params is None:
        params = CortexParams()
    if rng is None:
        rng = np.random.default_rng(0)
    net, _ = build_cortex_network(params, rng=rng, n_filaments=n_filaments,
                                  orientation=orientation, nematic_S=nematic_S,
                                  length_dist=length_dist)

    # Mesoscale bind reach (sanctioned ×40 coarse-graining dual) — the molecular ε (60 nm / 210 nm)
    # cannot connect the sparse coarse-grained filaments; widen to √(A/n_fil) as the archived cortex.
    reach = mesoscale_reach(params.R_um, n_filaments)

    # Crosslinkers (α-actinin / filamin split). rest length = the FORMATION distance so the link is
    # force-free at construction (the archived per-r0 binning rationale: construction energy ≪ kT);
    # tension develops only as myosin contracts the network — NOT pre-loaded by an arbitrary rest.
    xl_pairs = _cross_fiber_pairs(net, reach, n_xl, rng)
    nxl = xl_pairs.shape[0]
    is_alpha = rng.random(nxl) < alpha_fraction
    xl_k = np.where(is_alpha, ALPHA_ACTININ.link_k, FILAMIN.link_k)
    if nxl:
        xl_rest = np.linalg.norm(net.pos[xl_pairs[:, 1]] - net.pos[xl_pairs[:, 0]], axis=1)
    else:
        xl_rest = np.zeros(0)

    # Myosin links: same mesoscale reach — exclude already-used crosslinker pairs.
    used = {(int(a), int(b)) for a, b in xl_pairs}
    myo_pairs = _cross_fiber_pairs(net, reach, n_myo, rng, exclude=used)

    r0_mean = float(np.linalg.norm(net.pos - net.pos.mean(axis=0), axis=1).mean())
    return CrosslinkedCortex(
        net=net, xl_i=xl_pairs[:, 0], xl_j=xl_pairs[:, 1], xl_k=xl_k, xl_rest=xl_rest,
        myo_i=myo_pairs[:, 0], myo_j=myo_pairs[:, 1], R_um=params.R_um, R0_mean=r0_mean)


def _link_spring_force(pos: np.ndarray, i: np.ndarray, j: np.ndarray, k: np.ndarray,
                       rest: np.ndarray, out: np.ndarray) -> None:
    """Accumulate Hookean link forces k·(L−rest)·û into ``out`` (in place)."""
    if i.shape[0] == 0:
        return
    d = pos[j] - pos[i]
    L = np.linalg.norm(d, axis=1)
    safe = L > 1e-12
    u = np.zeros_like(d); u[safe] = d[safe] / L[safe, None]
    f = (k * (L - rest))[:, None] * u          # on i toward j when stretched
    np.add.at(out, i, f)
    np.add.at(out, j, -f)


def _myosin_force(pos: np.ndarray, i: np.ndarray, j: np.ndarray, f_myo: float,
                  out: np.ndarray) -> None:
    """Accumulate myosin contractile forces (constant magnitude f_myo pulling i,j together)."""
    if i.shape[0] == 0 or f_myo == 0.0:
        return
    d = pos[j] - pos[i]
    L = np.linalg.norm(d, axis=1)
    safe = L > 1e-12
    u = np.zeros_like(d); u[safe] = d[safe] / L[safe, None]
    f = f_myo * u                               # on i toward j (contractile = shorten)
    np.add.at(out, i, f)
    np.add.at(out, j, -f)


def turgor_pressure(cortex: CrosslinkedCortex, pos: np.ndarray, *, dP0: float | None = None,
                    pi_in0: float = TURGOR_PI_IN0, vmin_frac: float = VMIN_FRAC) -> tuple[float, float]:
    """State-dependent osmotic turgor ΔP(V) [pN/µm²] + mean radius [µm] (Guo 2017 entropic closure).

    ⚠ ``dP0`` (Π₀) is REQUIRED — PI 2026-07-25 gated it (``aleph/laws/params_turgor.yaml``,
    ``value: null``). Omitting it raises :class:`TurgorPi0Unspecified` naming the three competing
    in-tree values; pass ``dP0=TURGOR_DP0`` to keep the historical 40 Pa CONVENIENCE claim explicitly.

    Π_in(V) = N·kB·T/(V − Vmin) with Vmin = vmin_frac·V0 and N·kB·T = Π_in0·(V0−Vmin) (so Π_in(V0)=
    Π_in0); the medium balances all but the net resting excess, Π_out = Π_in0 − dP0, giving

        ΔP(V) = Π_in0·(V0−Vmin)/(V−Vmin) − (Π_in0 − dP0),   ΔP(V0) = dP0.

    This DERIVES the bulk modulus K_vol = −V·dΔP/dV|_{V0} = Π_in0/(1−vmin_frac) (≈7e5 pN/µm² at the
    lit inputs — no magic K_vol), and ΔP rises steeply as the shell is compressed (V↓), the
    physiological osmotic counter-force. ΔP floored at 0 (a deflated shell exerts no inward suction).
    """
    dP0 = resolve_pi0_pa(dP0, caller="ff.gamma_floor.turgor_pressure").value_pa
    centre = pos.mean(axis=0)
    R_mean = float(np.linalg.norm(pos - centre, axis=1).mean())
    V = (4.0 / 3.0) * np.pi * R_mean**3
    R0 = cortex.R0_mean if cortex.R0_mean > 0.0 else cortex.R_um   # self-consistent rest reference
    V0 = (4.0 / 3.0) * np.pi * R0**3
    vmin = vmin_frac * V0
    dP = pi_in0 * (V0 - vmin) / max(V - vmin, 1e-12 * V0) - (pi_in0 - dP0)
    return max(dP, 0.0), R_mean


def _turgor_force(cortex: CrosslinkedCortex, pos: np.ndarray, dP: float, R_mean: float,
                  out: np.ndarray) -> None:
    """Accumulate the radial outward Young-Laplace turgor force ΔP·(area/node)·r̂ (in place)."""
    if dP == 0.0:
        return
    centre = pos.mean(axis=0)
    d = pos - centre
    r = np.linalg.norm(d, axis=1)
    safe = r > 1e-12
    rhat = np.zeros_like(d); rhat[safe] = d[safe] / r[safe, None]
    area_per_node = 4.0 * np.pi * R_mean**2 / pos.shape[0]
    out += dP * area_per_node * rhat


def external_force(cortex: CrosslinkedCortex, pos: np.ndarray, bending_fn, f_myo: float,
                   *, turgor: bool = True, dP0: float | None = None) -> np.ndarray:
    """Total external force (bending + crosslinker springs + myosin + turgor) on all nodes (N,3) [pN].

    Turgor is ON by default (physiological-baseline HARD rule): the resting cell is turgor-
    pressurised and myosin is the modulator on top. The radial turgor force is the passive counter-
    force that lets the actively contracted shell reach equilibrium.

    ⚠ With ``turgor=True``, ``dP0`` (Π₀) is REQUIRED (PI 2026-07-25 gate — no silent default);
    omitting it raises :class:`TurgorPi0Unspecified`. With ``turgor=False`` Π₀ is unused and not
    demanded, so the actomyosin-only load path (the γ_active channel) is unaffected by the gate.
    """
    N = cortex.net.n_nodes
    F = np.asarray(bending_fn(pos.reshape(-1)), dtype=np.float64).reshape(N, 3)
    _link_spring_force(pos, cortex.xl_i, cortex.xl_j, cortex.xl_k, cortex.xl_rest, F)
    _myosin_force(pos, cortex.myo_i, cortex.myo_j, f_myo, F)
    if turgor:
        dP, R_mean = turgor_pressure(cortex, pos, dP0=dP0)
        _turgor_force(cortex, pos, dP, R_mean, F)
    return F


def equilibrate(cortex: CrosslinkedCortex, f_myo: float = 0.0, *, n_steps: int = 600,
                reshape_every: int = 25, dt_mu: float = 0.0, turgor: bool = False,
                dP0: float | None = None,
                method: str = "explicit", device: str = "cpu",
                crosslink_turnover: bool = False) -> np.ndarray:
    """Settle the cortex by projected overdamped descent + periodic reshape (NF2007 project-then-
    reshape, explicit; sufficient + robust for the quasi-static γ prototype).

    DEFAULT (``f_myo=0, turgor=False``) settles the RESTING passive shell — bending kinks relax while
    the fibers stay on the shell (the physiological-baseline resting geometry; myosin is then added
    as a MODULATOR at measurement time, per the hard rule). Passing ``f_myo>0`` / ``turgor=True``
    relaxes the actively-loaded shell instead — but note (Stage-6d finding) the floored actomyosin
    network has NO static equilibrium against physiological turgor (it inflates/collapses); that
    instability is itself the γ-floor signature, surfaced to PI rather than tuned away.

    ``method='implicit'`` uses the unconditionally-stable NF2007 Eq 2 step (``ff.relax``) on the full
    external force — far fewer force evaluations when the configuration is stiff (the resting bending
    settle is soft, so the default 'explicit' is already cheap there; implicit is the path for
    stiff/loaded shells).

    ``crosslink_turnover=True`` is the ROBUST resting baseline for STIFF (lit-anchored, ~4.6e5 pN/µm)
    crosslinks (Stage 6M/6N): settle the resting SHELL **bending-only** (crosslink force excluded → the
    stiff crosslinks never enter the step, so no CFL throttle / overshoot), then bind the crosslinks
    FORCE-FREE at the settled geometry (``xl_rest`` ← settled lengths; Hand §10.1 attach, the
    quasi-static / instantaneous-turnover limit of real crosslink rebinding). Without this, stiff
    crosslinks amplify the bending-settle residual stretch into a large SPURIOUS passive γ_xl (and the
    explicit/implicit settle can over-stretch). The active γ (myosin) is identical with or without it.

    ⚠ ``turgor=True`` now REQUIRES an explicit ``dP0`` (Π₀ gate, PI 2026-07-25) — the default
    ``turgor=False`` resting settle is unchanged and needs no Π₀.

    Returns the settled positions (N, 3) [µm]; also writes them into ``cortex.net.pos``.
    """
    from aleph.laws.constraints import project_constraint_forces

    net = cortex.net
    if crosslink_turnover:
        from aleph.laws.network_warp import relax_on_device
        net.pos = relax_on_device(net, links=None, n_steps=n_steps, reshape_every=reshape_every,
                                  dt_mu=dt_mu, device=device)
        if cortex.xl_i.size:                                  # crosslinks rebind force-free (Hand attach)
            cortex.xl_rest = np.linalg.norm(net.pos[cortex.xl_j] - net.pos[cortex.xl_i], axis=1)
        return net.pos
    bending_fn = make_bending_force_fn(net)
    if method == "device":
        # Fully on-device (Warp) settle — bending + crosslink springs (+ myosin), reshape, all on the
        # GPU (device="cuda:0" on gbook A5000). No turgor in the device loop (the resting baseline is
        # turgor-free per the physiological-baseline rule; loaded shells use the implicit path).
        from aleph.laws.network_warp import relax_on_device
        links = (np.stack([cortex.xl_i, cortex.xl_j], axis=1).astype(np.int64)
                 if cortex.xl_i.size else None)
        myo = (np.stack([cortex.myo_i, cortex.myo_j], axis=1).astype(np.int64)
               if (cortex.myo_i.size and f_myo) else None)
        net.pos = relax_on_device(
            net, links=links, k_xl=(cortex.xl_k if links is not None else None),
            xl_rest=(cortex.xl_rest if links is not None else None),
            myo_links=myo, f_myo=f_myo, n_steps=n_steps, reshape_every=reshape_every,
            dt_mu=dt_mu, device=device)
        return net.pos
    if method == "implicit":
        from aleph.laws.relax import _default_gamma, relax_implicit
        ext_fn = lambda xf: external_force(cortex, np.asarray(xf).reshape(net.n_nodes, 3),
                                           bending_fn, f_myo, turgor=turgor, dP0=dP0).reshape(-1)
        relax_implicit(net, ext_fn, gamma=_default_gamma(net), dt=1e3,
                       n_steps=max(1, n_steps // 30), reshape_every=5, project=False)
        return net.pos
    if dt_mu <= 0.0:
        seg = float(net.seg_rest.mean())
        k_bend = float(net.kappa.max()) / seg**3
        k_max = max(k_bend, float(cortex.xl_k.max()) if cortex.xl_i.size else 0.0)
        dt_mu = 0.1 / k_max
    x = net.pos.reshape(-1, 3).copy()
    for step in range(n_steps):
        net.pos = x
        F = external_force(cortex, x, bending_fn, f_myo, turgor=turgor, dP0=dP0)
        Fp = project_constraint_forces(net, F)          # tangent to the constraint manifold
        x = x + dt_mu * Fp
        if (step + 1) % reshape_every == 0:
            net.pos = x
            x = reshape(net, n_iter=2)
    net.pos = reshape(net, n_iter=6)
    return net.pos


def gamma_passive_young_laplace(dP: float, R_mean: float) -> float:
    """Young-Laplace pressure-balance tension γ = ΔP·R/2 [pN/µm] of the turgor-pressurised shell.

    ⚠️ BOOKKEEPING CAVEAT (2026-06-30, ultracode cortical-tension-active-fraction review). This is
    NOT independent evidence that "the cortex is mostly passive / at band". Two reasons:
      1. ``TURGOR_DP0`` (133 Pa) is **band-implied / tuned, no sourced row** (PARAM_AUDIT_SIMUCELL3D
         2026-06-25), so ΔP·R/2 "reproducing the band" is CIRCULAR (the turgor was set to make it so).
      2. Young-Laplace ΔP·R/2 is the in-plane tension that BALANCES the osmotic pressure (Stewart 2011:
         ΔP = γ·2/R) — it is the pressure's *partner*, not an independent passive cortical-tension
         channel; reporting both the turgor pressure and this tension double-books.
    The genuine blebbistatin-INSENSITIVE passive cortical floor is only ~0.04 mN/m (~9-12 % of band;
    Fischer-Friedrich 2014/2016) — and the band is ~50-90 % MYOSIN-dependent (central ~70 %). So the
    active actomyosin floor is the REAL gap, not a benign "passive carries it". Kept as a diagnostic
    (never folded into the actomyosin γ), but do not cite it as a passive-cortex result.
    See docs/v2_audit/FF_STAGE6D_GAMMA_FLOOR_2026-06-29.md §3.0c.
    """
    return 0.5 * dP * R_mean


def measure_gamma(cortex: CrosslinkedCortex, f_myo: float, *, n_planes: int = 50,
                  turgor: bool = True, dP0: float | None = None,
                  engaged_fraction: float = 1.0) -> dict:
    """Measure γ at the current ``cortex.net.pos``.

    ⚠ ``turgor=True`` REQUIRES an explicit ``dP0`` (Π₀ gate, PI 2026-07-25): γ_passive = Π₀·R/2 is
    ENTIRELY set by Π₀, so a silent Π₀ would silently set the reported passive tension. Omitting it
    raises :class:`TurgorPi0Unspecified`. The returned dict carries the Π₀ provenance record so any
    artifact built from it records the value AND that it is not sourced for MCF7. ``turgor=False``
    (the γ_active-only channel) needs no Π₀ and is unaffected.

    Two channels, kept SEPARATE (archived rule — turgor is NOT actomyosin tension):
      * γ_active (= γ_total here) — method-of-planes over the actomyosin load path: actin segments
        (axial tension from constraint multipliers), crosslinker links (k·(L−rest)), myosin links
        (contractile tension +f_myo).
      * γ_passive — Young-Laplace ΔP·R/2 of the turgor-pressurised shell.
    Returns both + the actomyosin sub-channels [pN/µm].

    ``engaged_fraction`` (Truong-Quang 2021): the load-bearing myosin-actin OVERLAP fraction. Only this
    fraction of placed minifilaments transmits contractile stress (the rest overhang, one-end-bound);
    since γ_myo is exactly linear in the contributing-dipole count, the reported ``gamma_myo`` = raw
    method-of-planes γ_myo × engaged_fraction (``gamma_myo_all`` keeps the raw all-placed value).
    Default 1.0 = the raw measurement; production uses the physiological interphase 0.65.
    """
    net = cortex.net
    pos = net.pos
    R = cortex.R_um
    bending_fn = make_bending_force_fn(net)
    # actomyosin load path ONLY (turgor is NOT actomyosin tension — kept to its own γ_passive channel)
    F = external_force(cortex, pos, bending_fn, f_myo, turgor=False)

    # actin backbone axial tension (per segment)
    seg_tau = segment_axial_tension(net, F)
    seg = net.segments
    rA_a, rB_a = pos[seg[:, 0]], pos[seg[:, 1]]

    # crosslinker links
    if cortex.xl_i.size:
        d = pos[cortex.xl_j] - pos[cortex.xl_i]
        L = np.linalg.norm(d, axis=1)
        xl_tau = cortex.xl_k * (L - cortex.xl_rest)
        rA_x, rB_x = pos[cortex.xl_i], pos[cortex.xl_j]
    else:
        xl_tau = np.zeros(0); rA_x = rB_x = np.zeros((0, 3))

    # myosin links (contractile tension = +f_myo)
    if cortex.myo_i.size:
        myo_tau = np.full(cortex.myo_i.size, f_myo)
        rA_m, rB_m = pos[cortex.myo_i], pos[cortex.myo_j]
    else:
        myo_tau = np.zeros(0); rA_m = rB_m = np.zeros((0, 3))

    centre = pos.mean(axis=0)
    g_all = method_of_planes_gamma(
        np.concatenate([rA_a, rA_x, rA_m]), np.concatenate([rB_a, rB_x, rB_m]),
        np.concatenate([seg_tau, xl_tau, myo_tau]), R, n_planes=n_planes, centre=centre)
    g_actin = method_of_planes_gamma(rA_a, rB_a, seg_tau, R, n_planes=n_planes, centre=centre)
    g_xl = method_of_planes_gamma(rA_x, rB_x, xl_tau, R, n_planes=n_planes, centre=centre)
    g_myo = method_of_planes_gamma(rA_m, rB_m, myo_tau, R, n_planes=n_planes, centre=centre)
    # γ_passive = the resting physiological turgor tension dP0·R0/2 (the baseline myosin modulates
    # FROM; physiological-baseline rule). Computed at the resting dP0/R0, NOT the settle geometry —
    # the bending-only settle shrinks the shell slightly (great-circle arcs straighten toward
    # chords), which is a dynamics artifact, not the resting turgor.
    if turgor:
        pi0 = resolve_pi0_pa(dP0, caller="ff.gamma_floor.measure_gamma")
        g_passive = gamma_passive_young_laplace(pi0.value_pa, R)
        pi0_record = pi0.as_artifact_record()
    else:
        g_passive = 0.0
        pi0_record = {"turgor_PI_0_Pa": 0.0, "turgor_PI_0_provenance": "TURGOR_OFF"}
    g_myo_engaged = g_myo * engaged_fraction        # Truong-Quang overlap: only the engaged fraction transmits
    return {"gamma_total": g_all, "gamma_active": g_all, "gamma_passive": g_passive,
            "gamma_actin": g_actin, "gamma_xl": g_xl, "gamma_myo": g_myo_engaged,
            "gamma_myo_all": g_myo, "engaged_fraction": engaged_fraction, **pi0_record}


def gamma_floor_run(f_myo: float, *, n_filaments: int = 100, n_xl: int = 300, n_myo: int = 100,
                    seed: int = 0, n_steps: int = 600, method: str = "explicit",
                    device: str = "cpu", crosslink_turnover: bool = True,
                    orientation: str = "isotropic", nematic_S: float = 1.0,
                    dP0: float | None = None,
                    engaged_fraction: float = 1.0) -> dict:
    """One quenched realization → γ at prestress ``f_myo``, measured at the resting physiological
    geometry (physiological-baseline rule: settle the resting passive shell, then add myosin as the
    modulator and measure). Returns the actomyosin γ channels + the passive turgor γ + meta.

    ``method='device'`` + ``device='cuda:0'`` runs the resting-shell settle fully on the GPU (Warp,
    gbook A5000) — the production GPU-resident path. ``crosslink_turnover=True`` (default) uses the
    robust force-free-rebinding resting baseline (Stage 6N) so the lit-anchored STIFF crosslink
    stiffness gives a clean γ (no spurious passive γ_xl); the active γ is identical either way.

    ⚠ ``dP0`` (Π₀) is REQUIRED because the γ_passive channel is measured (``turgor=True``): Π₀ is
    gated (PI 2026-07-25, ``aleph/laws/params_turgor.yaml``) and has no default. Pass
    ``dP0=TURGOR_DP0`` to reproduce the historical 40 Pa CONVENIENCE (HeLa-proxy) baseline; the
    returned dict then carries ``turgor_PI_0_provenance="CONVENIENCE"``."""
    rng = np.random.default_rng(seed)
    cortex = build_crosslinked_cortex(n_filaments=n_filaments, n_xl=n_xl, n_myo=n_myo, rng=rng,
                                      orientation=orientation, nematic_S=nematic_S)
    equilibrate(cortex, 0.0, n_steps=n_steps, turgor=False,   # settle resting passive shell
                method=method, device=device, crosslink_turnover=crosslink_turnover)
    out = measure_gamma(cortex, f_myo, turgor=True, dP0=dP0, engaged_fraction=engaged_fraction)
    out.update(f_myo=f_myo, seed=seed, n_xl=int(cortex.xl_i.size), n_myo=int(cortex.myo_i.size),
               orientation=orientation, nematic_S=(nematic_S if orientation != "isotropic" else 0.0))
    return out


def gamma_floor_production(*, f_myo: float = NMIIA_MINIFIL_STALL_PN, n_real: int = 4,
                           n_steps: int = 400, base_seed: int = 0, parallel: bool = True,
                           method: str = "explicit", device: str = "cpu",
                           n_filaments: int = PROD_N_FIL, n_xl: int = PROD_N_XL,
                           n_myo: int = PROD_N_MYO, dP0: float | None = None,
                           engaged_fraction: float = ENGAGED_FRACTION_INTERPHASE) -> dict:
    """γ at the GROUNDED production operating point — NATIVE ~38000 filaments (PI 2026-07-01; the ×40
    mesoscale 1000 is RETIRED, was a CPU constraint). ensemble.

    ⚠⚠ RETIRED-BUILD BANNER (2026-07-25). The headline this function produced — γ_active ≈ 6.2–6.6e-4
    mN/m, "floored ~530×" — was measured on the **pre-2026-07-23 cortex build** and has **NOT been
    re-run** since ``docs/v2_audit/CORTEX_STRUCTURE_AUDIT_2026-07-23.md`` declared that build broken
    (5 structural defects). THREE of those defects are STILL LIVE on this exact path:
    ``build_crosslinked_cortex`` takes ``n_xl`` from the caller and never reads
    ``ff/architecture_spec.CORTEX`` (every call site here passes ``n_xl = PROD_N_XL = n_filaments``,
    i.e. density_per_fil = 1.0, the under-percolating relic the 2026-07-23 fix replaced with 20);
    it omits ``resolve_overlaps`` (so ``build_cortex_network`` runs with the default ``False``); and
    it leaves ``length_dist="mono"``. Commit ``bc5ff3b0`` — whose own message names under-percolation
    as "the upstream root of … no myosin transmission" — reached only ``architecture_spec.py``, and
    the γ-floor measures EXACTLY myosin transmission. Do not cite the ~530× number until one native
    A5000 run of this function at the FIXED cortex reproduces or revises it.

    ⚠ ``dP0`` (Π₀) is REQUIRED (gated, PI 2026-07-25) because the γ_passive channel is measured.

    The single-cell γ-floor number; at native the MD-free γ_active ≈ 6.6e-4 mN/m (floor
    ~530×, N-converged) — the ×40 (1000) under-reported it ~5× (γ is N-dependent). Returns the γ
    distribution + the floor factor vs band. ``n_filaments``/``n_xl``/``n_myo`` override the native
    defaults (e.g. a small-N smoke test). ``method='device'`` + ``device='cuda:0'`` runs the ensemble
    GPU-resident on the gbook A5000 — REQUIRED at native scale (the CPU path is impractical).
    """
    from aleph.laws.gamma_estimator import MCF7_IQR_PN_UM, SALBREUX_BAND_PN_UM, active_band_pn_um
    ens = gamma_floor_ensemble(f_myo, n_real=n_real, n_filaments=n_filaments, n_xl=n_xl,
                               n_myo=n_myo, n_steps=n_steps, base_seed=base_seed,
                               parallel=parallel, method=method, device=device, dP0=dP0,
                               engaged_fraction=engaged_fraction)
    # FLOOR METRIC = γ_myo, the CLEAN myosin-induced channel (FF_STAGE6Q, 2026-07-01). At the native
    # lit-faithful density the γ_active SUM is contaminated by a STRUCTURAL passive actin-network
    # residual (γ_actin(f_myo=0) ≈ 0.84 pN/µm, plateaus 1000–8000 steps — crosslinked-geodesic
    # frustration on the curved shell, NOT under-relaxation); myosin contraction RELAXES that residual
    # so γ_active(f5) < γ_active(f0) and γ_active is NOT a usable myosin metric at native density. The
    # γ_myo channel is artifact-free: linear in f_myo, zero-intercept (γ_myo = 0.0199·f_myo) → the
    # honest active cortical tension. Floor is reported on γ_myo; γ_active is kept (flagged) for continuity.
    g_myo = np.array([r["gamma_myo"] for r in ens["runs"]])
    g_act = np.array([r["gamma_active"] for r in ens["runs"]])
    gm = g_myo.mean()
    mcf7_active_lo = active_band_pn_um(MCF7_IQR_PN_UM)[0]      # 70% of the MCF7 suspended IQR
    salbreux_active_lo = active_band_pn_um(SALBREUX_BAND_PN_UM)[0]
    f = lambda lo: float(lo / gm) if gm > 0 else float("inf")
    return {"gamma_myo_mean": float(gm), "gamma_myo_std": float(g_myo.std()),
            "gamma_myo_mN_per_m": float(gm) * 1e-3,
            "gamma_active_mean": float(g_act.mean()),          # passive-contaminated SUM (flagged, not the floor metric)
            "gamma_active_mN_per_m": float(g_act.mean()) * 1e-3,
            "gamma_passive": ens["runs"][0]["gamma_passive"],
            "floor_vs_mcf7_active": f(mcf7_active_lo),          # PRIMARY: clean γ_myo vs MCF7-faithful active target
            "floor_vs_salbreux_active": f(salbreux_active_lo),  # clean γ_myo vs generic active fraction
            "floor_factor_under_band": f(SALBREUX_BAND_PN_UM[0]),  # legacy: clean γ_myo vs generic TOTAL band
            "f_myo": f_myo, "n_fil": n_filaments, "n_xl": n_xl, "n_myo": n_myo,
            "engaged_fraction": engaged_fraction,   # Truong-Quang overlap (0.65 interphase / 1.0 mitotic)
            "gamma_myo_samples": g_myo.tolist(),
            # Π₀ provenance travels WITH the number (PI 2026-07-25 gate).
            **{k: v for k, v in ens["runs"][0].items() if k.startswith("turgor_PI_0_")},
            # Retired-cortex-build banner travels WITH the number (audit R1, 2026-07-25).
            "cortex_build_provenance": _CORTEX_BUILD_BANNER}


def gamma_floor_ensemble(f_myo: float, *, n_real: int = 8, n_filaments: int = 100,
                         n_xl: int = 300, n_myo: int = 100, n_steps: int = 4000,
                         base_seed: int = 0, parallel: bool = True,
                         method: str = "explicit", device: str = "cpu",
                         dP0: float | None = None,
                         engaged_fraction: float = 1.0) -> dict:
    """Run ``n_real`` quenched realizations → γ distribution at prestress ``f_myo``.

    ``method='device'`` forces SERIAL execution (the realizations run GPU-resident on ``device``; a
    ProcessPool would fork CUDA, which is unsafe — and the GPU is fast enough not to need it)."""
    seeds = [base_seed + r for r in range(n_real)]
    args = dict(n_filaments=n_filaments, n_xl=n_xl, n_myo=n_myo, n_steps=n_steps,
                dP0=dP0, engaged_fraction=engaged_fraction)
    if method == "device":
        runs = [gamma_floor_run(f_myo, seed=s, method="device", device=device, **args) for s in seeds]
    elif parallel and n_real > 1:
        from concurrent.futures import ProcessPoolExecutor
        with ProcessPoolExecutor() as ex:
            runs = list(ex.map(_run_one, [(f_myo, s, args) for s in seeds]))
    else:
        runs = [gamma_floor_run(f_myo, seed=s, **args) for s in seeds]
    g = np.array([r["gamma_total"] for r in runs])
    return {"f_myo": f_myo, "n_real": n_real, "gamma_mean": float(g.mean()),
            "gamma_std": float(g.std()), "gamma_median": float(np.median(g)),
            "gamma_samples": g.tolist(), "runs": runs}


def _run_one(packed):
    """Top-level worker for ProcessPoolExecutor (picklable)."""
    f_myo, seed, args = packed
    return gamma_floor_run(f_myo, seed=seed, **args)
