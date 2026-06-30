"""FF single-cell actomyosin network vs Taeyoon Kim's BD actin-network model (Stage 6i, peer oracle).

The correct validation target for the FF SINGLE-CELL network is a peer fine-grained model, NOT the
multi-cell spheroid A/A0 law. Taeyoon Kim's Brownian-dynamics cross-linked actin network (Kim, MS
thesis, MIT 2007, `references/181655768-MIT.pdf`, "Simulation of Actin Cytoskeleton Structure and
Rheology") is that peer: a 3D cubic box of polymerizing actin filaments cross-linked by ACPs, from
which network MORPHOLOGY (mesh size ξ, connectivity, cross-linking angle) and rheology are measured.

FF is athermal mechanical-equilibrium (no Brownian), so the directly-comparable Kim observables are
STRUCTURAL: the **mesh size ξ vs actin concentration C_A** (the semiflexible-network scaling
ξ ∝ C_A^(−1/2)) and the **network connectivity z vs the crosslinker ratio R = C_ACP/C_A**. (Thermal
G′/G″ from MSD are not FF-native; the elastic modulus is a follow-up.)

Grounded inputs (Kim 2007 + actin literature; none invented): actin rise ≈ 2.7 nm/monomer
(370 monomers/µm), monomer diameter 7 nm; C_A in µM → length density ρ_L = C_A·N_A·rise, mesh
ξ = ρ_L^(−1/2) (Schmidt/MacKintosh; gives ξ(151 µM) ≈ 64 nm, matching Kim's ~50–100 nm cortex mesh).
R = ACP:actin ratio (Kim uses R≈0.5); ACP capture radius = the crosslinker reach.
"""

from __future__ import annotations

import numpy as np

from ffn_sim.ff.fiber_network import build_fiber_network

N_AVOGADRO = 6.02214076e23
ACTIN_RISE_UM = 2.7e-3          # µm per monomer (≈370 monomers/µm; standard F-actin)

# Actin filament AXIAL stretching modulus EA (extensional rigidity). Kojima, Ishijima & Yanagida
# 1994, PNAS 91(26):12962 — direct glass-needle stretch, 43.7 ± 4.6 pN/nm over a 1 µm segment ⇒
# EA = k·L = 4.4e-8 N = 4.4e4 pN. (⚠️ 1000× trap: EA = 4.4e4 pN, NOT 4.4e7.) Used ONLY in the
# shear_modulus measurement to set the actin segment axial spring k_axial = EA/L_seg — it is NOT a
# production runtime constant. PI/SE: Kojima 1994 is not yet a SourceEvidence row (see
# references/SE_REGISTRATION_CANDIDATES_2026-06-30.md); register before any LIVE config use.
EA_ACTIN_PN = 4.4e4            # pN  (Kojima 1994; actin axial stretching modulus)


def length_density_per_um2(C_A_uM: float) -> float:
    """Actin contour-length density ρ_L [µm/µm³ = µm⁻²] from monomer concentration C_A [µM]."""
    # C_A [µM] = C_A·1e-6 mol/L = C_A·1e-6·N_A / 1e-3 m³ = C_A·1e-6·N_A·1e3 monomers/m³
    n_per_m3 = C_A_uM * 1e-6 * N_AVOGADRO * 1e3
    n_per_um3 = n_per_m3 * 1e-18                # monomers per µm³
    return n_per_um3 * ACTIN_RISE_UM           # µm contour per µm³ = µm⁻²


def mesh_size_theory_um(C_A_uM: float) -> float:
    """Semiflexible-network mesh size ξ = ρ_L^(−1/2) [µm] (Schmidt 1989 / MacKintosh)."""
    return length_density_per_um2(C_A_uM) ** -0.5


def build_box_network(C_A_uM: float, *, box_um: float = 1.0, L_fil_um: float = 0.5,
                      seg_um: float = 0.1, R_acp: float = 0.5, acp_reach_um: float = 0.05,
                      rng: np.random.Generator | None = None):
    """Build a 3D cubic-box cross-linked actin network at concentration ``C_A_uM`` (Kim-style).

    Filaments (random position + orientation, length ``L_fil_um``) are placed so the total contour
    length matches ρ_L·V; ACP crosslinks join filament nodes on different filaments within
    ``acp_reach_um``, up to R_acp · (number of filaments-ish) — here R sets the ACP budget relative to
    actin. Returns (net, xl_pairs, meta).
    """
    if rng is None:
        rng = np.random.default_rng(0)
    V = box_um**3
    rho_L = length_density_per_um2(C_A_uM)             # µm⁻² target contour density
    L_tot_target = rho_L * V                            # total contour length [µm]
    n_fil = max(1, int(round(L_tot_target / L_fil_um)))
    n_beads = max(2, int(round(L_fil_um / seg_um)) + 1)
    fibers = []
    for _ in range(n_fil):
        start = rng.random(3) * box_um
        d = rng.standard_normal(3); d /= np.linalg.norm(d)
        s = np.arange(n_beads) * seg_um
        fibers.append(start[None, :] + s[:, None] * d[None, :])
    net = build_fiber_network(fibers, kappa=7.3e-26)   # κ unused for morphology
    # crosslinks: node pairs on DIFFERENT filaments within acp_reach (ACP budget = R·n_fil)
    from scipy.spatial import cKDTree
    off = net.fiber_offsets
    fib_of = np.searchsorted(off, np.arange(net.n_nodes), side="right") - 1
    tree = cKDTree(net.pos)
    pairs = tree.query_pairs(r=acp_reach_um, output_type="ndarray")
    if pairs.shape[0]:
        pairs = pairs[fib_of[pairs[:, 0]] != fib_of[pairs[:, 1]]]
    budget = int(round(R_acp * n_fil))
    if pairs.shape[0] > budget and budget > 0:
        pairs = pairs[rng.choice(pairs.shape[0], budget, replace=False)]
    meta = {"C_A_uM": C_A_uM, "n_fil": n_fil, "rho_L_um2": rho_L,
            "L_tot_um": n_fil * L_fil_um, "box_um": box_um, "n_xl": int(pairs.shape[0])}
    return net, (pairs if pairs.shape[0] else np.zeros((0, 2), int)), meta


def measure_mesh_size_um(net, *, box_um: float, n_probe: int = 4000,
                         rng: np.random.Generator | None = None) -> float:
    """Measured mesh ξ [µm]: mean nearest-filament-node distance from random probe points × 2.

    A random point in a network of mesh ξ sits on average ~ξ/2 from the nearest filament, so
    ξ ≈ 2·⟨nearest distance⟩ (a standard geometric estimator; the scaling with C_A is the test)."""
    if rng is None:
        rng = np.random.default_rng(1)
    from scipy.spatial import cKDTree
    # densely resample fiber SEGMENTS into sub-points (~10 nm spacing) so nearest distance is to the
    # filament LINE, not the sparse nodes (sparse nodes can't resolve a sub-segment mesh → bias).
    seg = net.segments
    a, b = net.pos[seg[:, 0]], net.pos[seg[:, 1]]
    nsub = 10
    t = (np.arange(nsub) + 0.5) / nsub
    sub = (a[:, None, :] * (1 - t)[None, :, None] + b[:, None, :] * t[None, :, None]).reshape(-1, 3)
    tree = cKDTree(sub)
    probes = rng.random((n_probe, 3)) * box_um
    d, _ = tree.query(probes, k=1)
    return float(2.0 * d.mean())


def connectivity_z(net, xl_pairs) -> float:
    """Mean crosslinks per filament (Kim connectivity): 2·n_xl / n_fil."""
    n_fil = net.n_fibers
    return float(2.0 * xl_pairs.shape[0] / n_fil) if n_fil else 0.0


def shear_modulus(C_A_uM: float = 300.0, *, R_acp: float = 1.5, box_um: float = 1.0,
                  gamma: float = 0.03, k_xl: float = 10.0, k_axial: float | None = None,
                  n_steps: int = 4000, reshape_every: int = 20, seed: int = 0):
    """Athermal elastic SHEAR modulus G [pN/µm² = Pa] of the FF cross-linked actin network (Kim rheology).

    Applies affine simple shear u_x = γ·z, pins the top/bottom boundary at the sheared position,
    relaxes the interior, and reads σ_xz = (x-reaction on the top plane)/area; G = σ_xz/γ.
    (1 pN/µm² = 1 Pa exactly.) Cross-linked networks have an elastic floppy→rigid TRANSITION with
    connectivity (Head/Levine/MacKintosh 2003 PRE; Kim 2007).

    DEFORMABLE ELEMENTS (Stage-6c fix — previously a single ``k_xl`` was — wrongly — used for BOTH the
    actin backbone and the crosslink, conflating them; and at the broken value the backbone was as soft
    as a crosslink):
      * **Actin backbone** — two modes. ``k_axial=None`` (DEFAULT) ⇒ INEXTENSIBLE (FF-native NF2007
        §5.3 reshape, like the rest of the engine): the correct stiff limit AND tractable for dense
        networks (a 1.8e5 pN/µm explicit spring needs ~1e6 steps). Valid in the crosslink-limited
        regime k_xl ≪ EA/L_seg. ``k_axial=<float>`` ⇒ a FINITE axial spring (= EA/L_seg, Kojima 1994
        EA = ``EA_ACTIN_PN`` = 4.4e4 pN) — needed in the sourced-stiffness regime where k_xl ~ EA/L_seg
        (α-actinin), so actin compliance co-determines G; tractable only for small (dilute, Kim-matched)
        networks (few nodes) because the stiff spring forces a tiny step.
      * ``k_xl`` — crosslink JUNCTION stiffness [pN/µm], the soft deformable element. Lit anchors
        (Ferrer 2008 PNAS, AFM): α-actinin 455 pN/nm = 4.6e5 pN/µm, filamin 820 pN/nm = 8.2e5 pN/µm.
        ⚠️ The default 10.0 and the ``hand_kmc`` crosslinker ``link_k=0.1`` are UNSOURCED magnitude
        knobs (a confirmed pN/µm-vs-pN/nm slip compounded with an AFINES soft surrogate) — correcting
        the LIVE production constant is PI-gated (register a Ferrer-stiffness SourceEvidence row first).

    Regimes (the robust Kim/analytic comparison — assert the TRANSITION + LINEARITY, NOT an absolute):
      * crosslink-limited (k_xl ≪ bending stiffness κ/seg³): G is LINEAR in k_xl (the junction is the
        soft series element). Analytic cross-check G ≈ k_xl·ρ_L·ℓc (primary ground-truth).
      * backbone-limited (k_xl ≫ κ/seg³): G SATURATES at the bending/inextensibility-set enthalpic
        value. FF is athermal ⇒ the enthalpic branch; the dilute in-vitro Kim/Gardel G' (0.1–1000 Pa
        low-f) is the THERMAL branch FF does not target.

    Returns (G [Pa], connectivity z, meta) with meta += rho_L_per_um2, lc_um, k_xl, G_analytic_crosslink.
    """
    from ffn_sim.ff.constraints import reshape
    from ffn_sim.ff.forces_warp import make_bending_force_fn
    net, xl, meta = build_box_network(C_A_uM, box_um=box_um, R_acp=R_acp,
                                      rng=np.random.default_rng(seed))
    N = net.n_nodes
    bfn = make_bending_force_fn(net)
    kappa = float(net.kappa.max())
    seg = float(net.seg_rest.mean())
    use_reshape = k_axial is None                              # rigid actin (reshape) vs finite axial spring
    sp = net.segments
    pos0 = net.pos.copy()
    r0 = np.linalg.norm(pos0[sp[:, 1]] - pos0[sp[:, 0]], axis=1)
    x = pos0.copy()
    x[:, 0] += gamma * x[:, 2]                                  # affine simple shear
    margin = 0.12 * box_um
    top = pos0[:, 2] > box_um - margin
    pinned = (pos0[:, 2] < margin) | top
    pinpos = x[pinned].copy()
    dt_mu = 0.2 / max(16 * kappa / seg**3, k_xl, (k_axial or 0.0))  # CFL: stiffest spring (incl finite k_axial)

    def F(xx):
        f = bfn(xx.reshape(-1)).reshape(N, 3).copy()           # bending (actin backbone resists via κ)
        if len(xl):                                            # crosslink junctions (k_xl) — soft element
            d = xx[xl[:, 1]] - xx[xl[:, 0]]; f0 = k_xl * d
            np.add.at(f, xl[:, 0], f0); np.add.at(f, xl[:, 1], -f0)
        if not use_reshape:                                    # finite actin axial spring (sourced EA/L_seg)
            d = xx[sp[:, 1]] - xx[sp[:, 0]]; L = np.linalg.norm(d, axis=1) + 1e-12
            fa = (k_axial * (L - r0))[:, None] * (d / L[:, None])
            np.add.at(f, sp[:, 0], fa); np.add.at(f, sp[:, 1], -fa)
        return f

    for step in range(n_steps):
        ff = F(x); ff[pinned] = 0.0
        if not np.isfinite(ff).all():
            return None
        x = x + dt_mu * ff; x[pinned] = pinpos
        if use_reshape and (step + 1) % reshape_every == 0:    # inextensible actin (NF2007 §5.3), then re-pin
            net.pos = x; x = reshape(net, n_iter=2); x[pinned] = pinpos
    react = -F(x)[top]
    sigma_xz = react[:, 0].sum() / (box_um * box_um)
    # analytic cross-check inputs (crosslink-limited G ≈ k_xl·ρ_L·ℓc; primary ground-truth)
    rho_L = length_density_per_um2(C_A_uM)                      # µm⁻² contour-length density
    z = connectivity_z(net, xl)
    # crosslink spacing ℓc = total contour length / (2·n_xl) (each crosslink = one point on each of 2 fibers)
    lc = float(net.seg_rest.sum()) / (2.0 * xl.shape[0]) if xl.shape[0] else float("inf")
    meta = {**meta, "rho_L_per_um2": rho_L, "lc_um": lc, "k_xl": float(k_xl),
            "G_analytic_crosslink": float(k_xl * rho_L * (lc if np.isfinite(lc) else 0.0))}
    return float(sigma_xz / gamma), z, meta
