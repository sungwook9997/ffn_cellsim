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

from dataclasses import dataclass

import numpy as np

from ffn_sim.ff import units as U
from ffn_sim.ff.constraints import reshape, segment_axial_tension
from ffn_sim.ff.cortex_assembly import CortexParams, build_cortex_network
from ffn_sim.ff.fiber_network import FiberNetwork
from ffn_sim.ff.forces_warp import make_bending_force_fn
from ffn_sim.ff.gamma_estimator import method_of_planes_gamma
from ffn_sim.ff.hand_kmc import ALPHA_ACTININ, FILAMIN, NMIIA_MYOSIN


# NMIIA minifilament stall force (lit anchor for the prestress sweep), FF units (pN).
# F_stall_per_head=0.5 pN, n_heads_per_side=10 (configs/phase1_h3.yaml). Two sides engage against
# the two actin tracks; a high-duty-ratio NMIIA engages a sizeable fraction under load. We take the
# per-side stall (10 heads × 0.5 pN = 5 pN) as the grounded anchor and SWEEP around it.
NMIIA_F_STALL_PER_HEAD = 0.5
NMIIA_HEADS_PER_SIDE = 10
NMIIA_MINIFIL_STALL_PN = NMIIA_F_STALL_PER_HEAD * NMIIA_HEADS_PER_SIDE   # 5.0 pN (per side)

# Physiological turgor (the passive counter-force; physiological-baseline HARD rule). FF units:
# 1 Pa = 1 pN/µm². Grounded: turgor_dP0=133 Pa (SOLID, registry LIVE — band-implied MCF7 baseline),
# K_vol=1e3 Pa osmotic bulk modulus (ResolvedDCM). Young-Laplace passive γ = ΔP·R/2.
TURGOR_DP0 = 133.0          # pN/µm²  (= 133 Pa)
TURGOR_K_VOL = 1.0e3        # pN/µm²  (ΔP per ΔV/V)

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
PROD_N_FIL = 1000
PROD_N_XL = 1000
PROD_N_MYO = 100   # ×40-mesoscale minifilament count (provenance per the citation-fix note above)

# Nie et al. 2015 measured cortical NMII minifilament areal density [µm⁻²] — the only direct datum.
NIE2015_DENSITY_UM2 = 0.625
NIE2015_DENSITY_RANGE = (0.31, 0.94)


@dataclass(slots=True)
class CrosslinkedCortex:
    """A cortex fiber network plus its crosslinker + myosin links (one quenched realization)."""

    net: FiberNetwork
    xl_i: np.ndarray            # (Nxl,) node index endpoint A
    xl_j: np.ndarray            # (Nxl,) node index endpoint B (different fiber)
    xl_k: np.ndarray            # (Nxl,) link stiffness [pN/µm]
    xl_rest: np.ndarray         # (Nxl,) link rest length [µm]
    myo_i: np.ndarray           # (Nmyo,) node index endpoint A
    myo_j: np.ndarray           # (Nmyo,) node index endpoint B (different fiber)
    R_um: float


def mesoscale_reach(R_um: float, n_filaments: int) -> float:
    """Mesoscale crosslinker/myosin reach ``√(A_shell / n_fil)`` [µm] — the geometric dual of the
    ×40 mesoscopic coarse-graining (the ONLY sanctioned coarse-graining; Plan v2 §3 H.3 v3.1). At
    ×40 the ~1000 effective filaments are far sparser than the native mesh, so the molecular
    crosslinker ε (60 nm) never connects coarse filaments; the archived connected cortex
    (``connected_mesh.mesoscale_reach``) widens the bind span to this mesoscale value so the
    coarse-grained network percolates. Faithful port of that rule."""
    return float(np.sqrt(4.0 * np.pi * R_um**2 / max(n_filaments, 1)))


def _node_fiber_map(net: FiberNetwork) -> np.ndarray:
    """(N,) fiber index of each node."""
    off = net.fiber_offsets
    return np.searchsorted(off, np.arange(net.n_nodes), side="right") - 1


def _cross_fiber_pairs(net: FiberNetwork, capture_um: float, n_links: int,
                       rng: np.random.Generator, exclude: set | None = None) -> np.ndarray:
    """Sample up to ``n_links`` node pairs on DIFFERENT fibers within ``capture_um`` of each other.

    Returns (M, 2) node-index pairs (M ≤ n_links). Uses a KD-tree pair query.
    """
    from scipy.spatial import cKDTree

    fib = _node_fiber_map(net)
    tree = cKDTree(net.pos)
    pairs = tree.query_pairs(r=capture_um, output_type="ndarray")
    if pairs.shape[0] == 0:
        return np.zeros((0, 2), np.int64)
    diff = fib[pairs[:, 0]] != fib[pairs[:, 1]]        # keep only cross-fiber pairs
    pairs = pairs[diff]
    if exclude:
        keep = np.array([(int(a), int(b)) not in exclude for a, b in pairs], bool)
        pairs = pairs[keep]
    if pairs.shape[0] == 0:
        return np.zeros((0, 2), np.int64)
    if pairs.shape[0] > n_links:
        sel = rng.choice(pairs.shape[0], size=n_links, replace=False)
        pairs = pairs[sel]
    return pairs.astype(np.int64)


def build_crosslinked_cortex(params: CortexParams | None = None, *, n_filaments: int = 100,
                             n_xl: int = 300, n_myo: int = 100,
                             alpha_fraction: float = 0.30,
                             rng: np.random.Generator | None = None) -> CrosslinkedCortex:
    """Assemble one quenched cortex realization: fibers + crosslinker links + myosin links.

    Crosslinkers are α-actinin (``alpha_fraction``, Stricker 2010 30 %) / filamin (rest); each links
    two near nodes on different fibers (capture radius = the crosslinker Hand's ε). Myosin links
    connect cross-fiber node pairs within the NMIIA capture radius. Counts scale with the H.3 cortex
    (n_xl≈n_filaments·... ; defaults sized for a fast prototype).
    """
    if params is None:
        params = CortexParams()
    if rng is None:
        rng = np.random.default_rng(0)
    net, _ = build_cortex_network(params, rng=rng, n_filaments=n_filaments)

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

    return CrosslinkedCortex(
        net=net, xl_i=xl_pairs[:, 0], xl_j=xl_pairs[:, 1], xl_k=xl_k, xl_rest=xl_rest,
        myo_i=myo_pairs[:, 0], myo_j=myo_pairs[:, 1], R_um=params.R_um)


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


def turgor_pressure(cortex: CrosslinkedCortex, pos: np.ndarray, *, dP0: float = TURGOR_DP0,
                    K_vol: float = TURGOR_K_VOL) -> tuple[float, float]:
    """Enclosed-volume turgor ΔP [pN/µm²] and mean radius [µm] at ``pos`` (sphere-proxy volume).

    ΔP = dP0 + K_vol·(V0−V)/V0 with V = 4/3·π·R_mean³ (R_mean = mean node radius), V0 from the rest
    radius. Contraction (R↓ ⇒ V↓) raises ΔP → the radial counter-force that lets the actively
    contracted shell reach a Young-Laplace equilibrium. ΔP floored at 0 (a deflated shell exerts no
    inward suction here).
    """
    centre = pos.mean(axis=0)
    R_mean = float(np.linalg.norm(pos - centre, axis=1).mean())
    V = (4.0 / 3.0) * np.pi * R_mean**3
    V0 = (4.0 / 3.0) * np.pi * cortex.R_um**3
    dP = dP0 + K_vol * (V0 - V) / V0
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
                   *, turgor: bool = True, dP0: float = TURGOR_DP0,
                   K_vol: float = TURGOR_K_VOL) -> np.ndarray:
    """Total external force (bending + crosslinker springs + myosin + turgor) on all nodes (N,3) [pN].

    Turgor is ON by default (physiological-baseline HARD rule): the resting cell is turgor-
    pressurised and myosin is the modulator on top. The radial turgor force is the passive counter-
    force that lets the actively contracted shell reach equilibrium.
    """
    N = cortex.net.n_nodes
    F = np.asarray(bending_fn(pos.reshape(-1)), dtype=np.float64).reshape(N, 3)
    _link_spring_force(pos, cortex.xl_i, cortex.xl_j, cortex.xl_k, cortex.xl_rest, F)
    _myosin_force(pos, cortex.myo_i, cortex.myo_j, f_myo, F)
    if turgor:
        dP, R_mean = turgor_pressure(cortex, pos, dP0=dP0, K_vol=K_vol)
        _turgor_force(cortex, pos, dP, R_mean, F)
    return F


def equilibrate(cortex: CrosslinkedCortex, f_myo: float = 0.0, *, n_steps: int = 600,
                reshape_every: int = 25, dt_mu: float = 0.0, turgor: bool = False,
                method: str = "explicit") -> np.ndarray:
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

    Returns the settled positions (N, 3) [µm]; also writes them into ``cortex.net.pos``.
    """
    from ffn_sim.ff.constraints import project_constraint_forces

    net = cortex.net
    bending_fn = make_bending_force_fn(net)
    if method == "implicit":
        from ffn_sim.ff.relax import _default_gamma, relax_implicit
        ext_fn = lambda xf: external_force(cortex, np.asarray(xf).reshape(net.n_nodes, 3),
                                           bending_fn, f_myo, turgor=turgor).reshape(-1)
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
        F = external_force(cortex, x, bending_fn, f_myo, turgor=turgor)
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
                  turgor: bool = True) -> dict:
    """Measure γ at the current ``cortex.net.pos``.

    Two channels, kept SEPARATE (archived rule — turgor is NOT actomyosin tension):
      * γ_active (= γ_total here) — method-of-planes over the actomyosin load path: actin segments
        (axial tension from constraint multipliers), crosslinker links (k·(L−rest)), myosin links
        (contractile tension +f_myo).
      * γ_passive — Young-Laplace ΔP·R/2 of the turgor-pressurised shell.
    Returns both + the actomyosin sub-channels [pN/µm].
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
    g_passive = gamma_passive_young_laplace(TURGOR_DP0, R) if turgor else 0.0
    return {"gamma_total": g_all, "gamma_active": g_all, "gamma_passive": g_passive,
            "gamma_actin": g_actin, "gamma_xl": g_xl, "gamma_myo": g_myo}


def gamma_floor_run(f_myo: float, *, n_filaments: int = 100, n_xl: int = 300, n_myo: int = 100,
                    seed: int = 0, n_steps: int = 600) -> dict:
    """One quenched realization → γ at prestress ``f_myo``, measured at the resting physiological
    geometry (physiological-baseline rule: settle the resting passive shell, then add myosin as the
    modulator and measure). Returns the actomyosin γ channels + the passive turgor γ + meta."""
    rng = np.random.default_rng(seed)
    cortex = build_crosslinked_cortex(n_filaments=n_filaments, n_xl=n_xl, n_myo=n_myo, rng=rng)
    equilibrate(cortex, 0.0, n_steps=n_steps, turgor=False)   # settle resting passive shell
    out = measure_gamma(cortex, f_myo, turgor=True)
    out.update(f_myo=f_myo, seed=seed, n_xl=int(cortex.xl_i.size), n_myo=int(cortex.myo_i.size))
    return out


def gamma_floor_production(*, f_myo: float = NMIIA_MINIFIL_STALL_PN, n_real: int = 4,
                           n_steps: int = 400, base_seed: int = 0, parallel: bool = True) -> dict:
    """γ at the GROUNDED production operating point (N=1000 / n_xl=1000 / n_myo=100), ensemble.

    All counts from configs/phase1_h3.yaml (the ×40 mesoscale cell) — no prototype downscaling. This
    is the defensible single-cell γ-floor number; at this point the MD-free γ_active lands on the
    archived BAOAB-MD g_soft (~1.4e-4 mN/m). Returns the γ distribution + the floor factor vs band.
    """
    from ffn_sim.ff.gamma_estimator import SALBREUX_BAND_PN_UM
    ens = gamma_floor_ensemble(f_myo, n_real=n_real, n_filaments=PROD_N_FIL, n_xl=PROD_N_XL,
                               n_myo=PROD_N_MYO, n_steps=n_steps, base_seed=base_seed,
                               parallel=parallel)
    g = np.array([r["gamma_active"] for r in ens["runs"]])
    band_lo = SALBREUX_BAND_PN_UM[0]
    return {"gamma_active_mean": float(g.mean()), "gamma_active_std": float(g.std()),
            "gamma_active_mN_per_m": float(g.mean()) * 1e-3,
            "gamma_passive": ens["runs"][0]["gamma_passive"],
            "floor_factor_under_band": float(band_lo / g.mean()) if g.mean() > 0 else float("inf"),
            "f_myo": f_myo, "n_fil": PROD_N_FIL, "n_xl": PROD_N_XL, "n_myo": PROD_N_MYO,
            "gamma_samples": g.tolist()}


def gamma_floor_ensemble(f_myo: float, *, n_real: int = 8, n_filaments: int = 100,
                         n_xl: int = 300, n_myo: int = 100, n_steps: int = 4000,
                         base_seed: int = 0, parallel: bool = True) -> dict:
    """Run ``n_real`` quenched realizations → γ distribution at prestress ``f_myo``."""
    seeds = [base_seed + r for r in range(n_real)]
    args = dict(n_filaments=n_filaments, n_xl=n_xl, n_myo=n_myo, n_steps=n_steps)
    if parallel and n_real > 1:
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
