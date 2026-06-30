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


def shear_modulus(C_A_uM: float = 300.0, *, R_acp: float = 0.5, box_um: float = 1.0,
                  gamma: float = 0.03, k_xl: float = 10.0, n_steps: int = 3000, seed: int = 0):
    """Athermal elastic SHEAR modulus G [pN/µm²] of the FF cross-linked actin network (Kim rheology).

    Applies affine simple shear u_x = γ·z, pins the top/bottom boundary at the sheared position,
    relaxes the interior (actin segment springs + crosslink springs + bending), and reads the shear
    stress σ_xz = (x-reaction on the top plane)/area; G = σ_xz/γ. Cross-linked networks have an elastic
    floppy→rigid TRANSITION with connectivity (Head/Levine/MacKintosh 2003 PRE; Kim 2007): G≈0 below
    threshold, rising steeply above. ``k_xl`` sets the element stiffness scale (the G MAGNITUDE; the
    transition/trend is the robust Kim comparison). Returns (G, z, meta).
    """
    net, xl, meta = build_box_network(C_A_uM, box_um=box_um, R_acp=R_acp,
                                      rng=np.random.default_rng(seed))
    from ffn_sim.ff.forces_warp import make_bending_force_fn
    N = net.n_nodes
    bfn = make_bending_force_fn(net)
    kappa = float(net.kappa.max())
    seg = float(net.seg_rest.mean())
    pos0 = net.pos.copy()
    x = pos0.copy()
    x[:, 0] += gamma * x[:, 2]                                  # affine simple shear
    margin = 0.12 * box_um
    top = pos0[:, 2] > box_um - margin
    pinned = (pos0[:, 2] < margin) | top
    pinpos = x[pinned].copy()
    sp = net.segments
    r0 = np.linalg.norm(pos0[sp[:, 1]] - pos0[sp[:, 0]], axis=1)
    dt_mu = 0.2 / max(16 * kappa / seg**3, k_xl)

    def F(xx):
        f = bfn(xx.reshape(-1)).reshape(N, 3).copy()
        if len(xl):
            d = xx[xl[:, 1]] - xx[xl[:, 0]]; f0 = k_xl * d
            np.add.at(f, xl[:, 0], f0); np.add.at(f, xl[:, 1], -f0)
        d = xx[sp[:, 1]] - xx[sp[:, 0]]; L = np.linalg.norm(d, axis=1) + 1e-12
        fa = (k_xl * (L - r0))[:, None] * (d / L[:, None])
        np.add.at(f, sp[:, 0], fa); np.add.at(f, sp[:, 1], -fa)
        return f

    for _ in range(n_steps):
        ff = F(x); ff[pinned] = 0.0
        if not np.isfinite(ff).all():
            return None
        x = x + dt_mu * ff; x[pinned] = pinpos
    react = -F(x)[top]
    sigma_xz = react[:, 0].sum() / (box_um * box_um)
    return float(sigma_xz / gamma), connectivity_z(net, xl), meta
