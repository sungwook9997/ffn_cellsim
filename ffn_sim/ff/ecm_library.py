"""FF ECM library — grounded, multi-material, 2D/3D, alignment-controlled matrices.

The ECM a cell lives in is not one material. This module is the going-forward FF ECM builder that
generalises ``ecm_mikado`` (collagen-only, 3D, isotropic) to the full space the PI asked for
(2026-07-10): several ECM MATERIALS (individual + mixed), 2D and 3D, several ALIGNMENT degrees,
each set from literature so the *emergent* stiffness reproduces the real Pa.

Two constitutive classes (this distinction is the physics, not a convenience):

* **Fibrillar ECM** — collagen-I, fibrin, agarose. A Mikado semiflexible-fiber network (Wilhelm &
  Frey 2003 PRL 91:108103; Head-Levine-MacKintosh 2003 PRE 68:061907; Broedersz & MacKintosh 2014
  Rev.Mod.Phys. 86:995). Straight rods placed at random with a controllable orientation
  distribution, cross-linked at near-contacts. The macroscopic modulus **emerges** from the
  microstructure (fiber length density → connectivity ⟨z⟩, bending κ=k_BT·L_p, stretch EA, crosslink
  k_xl). We set the microstructure from literature and *validate* the emergent modulus — we never
  tune the network to hit a target Pa (hard rule: no-parameter-tuning-to-outcome).

* **Continuum gel** — polyacrylamide (PAA), hyaluronic acid, Matrigel. A flexible chemical / basement-
  membrane gel whose network mesh (~nm) is unresolvable at cell scale; it behaves as a linear-elastic
  continuum with a modulus SET BY CHEMISTRY (a material input, not an emergent microstructural
  quantity). Represented as a coarse random elastic spring network whose bond stiffness is the
  analytic inverse of the target E (a constitutive discretization — Cauchy-Born affine estimate —
  not a fit). Indentation then validates that pressing returns the input E.

Units: µm · pN · s (FF convention, ``ff.units``). 1 pN/µm² = 1 Pa exactly, so every modulus this
library feeds ``ecm_mechanics`` comes out directly in Pa. κ[pN·µm²]=k_BT·L_p; EA[pN]=E_fibril[Pa]·π·r_f²[µm²].
"""

from __future__ import annotations

from dataclasses import dataclass, field, replace

import numpy as np
from scipy.spatial import cKDTree

from ffn_sim.ff import units as U
from ffn_sim.ff.fiber_network import FiberNetwork, build_fiber_network, concat_fiber_networks

# ─────────────────────────────────────────────────────────────────────────────────────────────────
# Material registry
# ─────────────────────────────────────────────────────────────────────────────────────────────────


@dataclass(frozen=True)
class ECMSpec:
    """Grounded parameter card for one ECM material (literature inputs; sources in ``sources``).

    Fibrillar materials carry the fiber + network params (``lp_um``, ``fiber_diameter_um``,
    ``E_fibril_Pa``, ``fiber_len_um``, ``k_xl_pN_um``, ``xl_contact_um``) and the concentration→mesh
    anchor (``ref_conc``/``ref_mesh_um``/``conc_exponent``). Continuum gels carry ``E_gel_Pa`` (the
    material Young's modulus input) and ``poisson``; their fiber fields are unused.
    """

    key: str
    name: str
    is_fibrillar: bool
    # validation band (Pa) at reference conditions + which modulus it refers to ("E" or "G")
    modulus_band_Pa: tuple[float, float]
    modulus_kind: str
    poisson: float
    # ── fibrillar params ──
    lp_um: float = 0.0                 # persistence length → κ = k_BT·L_p
    fiber_diameter_um: float = 0.0     # → r_f, EA = E_fibril·π·r_f²
    E_fibril_Pa: float = 0.0           # single-fiber material Young's modulus (axial stretch)
    fiber_len_um: float = 0.0          # rod contour length L_f
    seg_um: float = 0.5                # discretisation segment length
    k_xl_pN_um: float = 0.0            # crosslink Hookean stiffness
    xl_contact_um: float = 0.0         # near-contact crosslink distance (~fiber diameter scale)
    ref_conc: float = 0.0              # reference concentration
    ref_conc_unit: str = ""            # "mg/mL" | "%w/v" | "kPa-input"
    ref_mesh_um: float = 0.0           # mesh ξ at ref_conc
    conc_exponent: float = 0.5         # ξ ~ conc^(−conc_exponent) (3D collagen: 0.5)
    strain_stiffens: bool = False
    # ── continuum-gel params ──
    E_gel_Pa: float = 0.0              # target Young's modulus (material input for continuum gels)
    # ── documentation ──
    tissue_by_S: tuple[tuple[float, str], ...] = ()   # (S_order, tissue) mapping
    aliases: tuple[str, ...] = ()
    sources: tuple[str, ...] = ()
    notes: str = ""

    @property
    def kappa_pN_um2(self) -> float:
        """Fiber bending rigidity κ = k_BT·L_p [pN·µm²] (measured L_p, not the solid-rod E·I — collagen
        fibrils slide within a bundle, so bending is far softer than E·I; stretch uses the full EA)."""
        return U.KBT * self.lp_um

    @property
    def r_f_um(self) -> float:
        """Fiber radius [µm]."""
        return 0.5 * self.fiber_diameter_um

    @property
    def EA_pN(self) -> float:
        """Axial stretch modulus EA = E_fibril · π r_f² [pN] (a force; k_seg = EA/seg)."""
        return self.E_fibril_Pa * np.pi * self.r_f_um ** 2

    def k_seg_pN_um(self, seg_um: float | None = None) -> float:
        """Per-segment axial spring stiffness k = EA/ℓ₀ [pN/µm] — derived, finite-extensibility (enables
        the emergent strain-stiffening gate). Falls back to a stiff default if EA is unset."""
        seg = self.seg_um if seg_um is None else seg_um
        return self.EA_pN / seg if self.EA_pN > 0 else 5.0e4


# Grounded material cards. Fibrillar cards anchor to the KB (collagen Unit-1 ECM); the others are
# anchored to primary literature (DOIs registered via the KB-registration step). Values are literature
# inputs — none chosen to make a modulus gate pass.
REGISTRY: dict[str, ECMSpec] = {
    "collagen_I": ECMSpec(
        key="collagen_I", name="Collagen-I", is_fibrillar=True,
        modulus_band_Pa=(30.0, 100.0), modulus_kind="G", poisson=0.2,
        lp_um=17.0, fiber_diameter_um=0.10, E_fibril_Pa=1.1e6, fiber_len_um=10.0, seg_um=0.5,
        k_xl_pN_um=1000.0, xl_contact_um=0.75,
        ref_conc=1.5, ref_conc_unit="mg/mL", ref_mesh_um=2.0, conc_exponent=0.5, strain_stiffens=True,
        tissue_by_S=((0.05, "loose stroma / dermis (healthy, TACS-1 dense-isotropic)"),
                     (0.35, "tumor stroma tangential (TACS-2)"),
                     (0.6, "tumor invasion highway (TACS-3 radial-aligned)"),
                     (0.85, "tendon / ligament / aligned scar")),
        aliases=("collagen", "col1", "type-I collagen", "fibrillar collagen"),
        sources=("KB-1.1", "KB-1.2", "KB-1.3", "KB-1.7", "KB-1.30", "KB-1.V.2.1",
                 "Yang-Kaufman 2009 BiophysJ 10.1016/j.bpj.2008.10.063",
                 "Licup 2015 PNAS 10.1073/pnas.1504258112"),
        notes="Emergent G'(c): 1→5, 3→55, 7→342 Pa (n~2.0-2.1). ⟨z⟩~3.2 sub-isostatic bending-dominated."),
    "fibrin": ECMSpec(
        key="fibrin", name="Fibrin", is_fibrillar=True,
        modulus_band_Pa=(10.0, 1000.0), modulus_kind="G", poisson=0.2,
        lp_um=1.0, fiber_diameter_um=0.13, E_fibril_Pa=1.7e6, fiber_len_um=8.0, seg_um=0.4,
        k_xl_pN_um=1000.0, xl_contact_um=0.6,
        ref_conc=2.0, ref_conc_unit="mg/mL", ref_mesh_um=1.5, conc_exponent=0.5, strain_stiffens=True,
        tissue_by_S=((0.05, "blood clot / provisional wound matrix (isotropic)"),
                     (0.5, "aligned fibrin under flow / tissue-engineered constructs")),
        aliases=("fibrin gel", "fibrin clot"),
        sources=("Piechocka 2010 BiophysJ 10.1016/j.bpj.2010.01.040 (G0 0.1-2000 Pa)",
                 "Collet 2005 PNAS 10.1073/pnas.0504120102 (single-fiber E=1.7±1.3 MPa)",
                 "Storm 2005 Nature 10.1038/nature03521", "Liu 2006 Science 10.1126/science.1127317"),
        notes="Semiflexible, highly extensible fibers (strain to 180%), strong strain-stiffening. "
              "G0 0.1-2000 Pa over 0.1-8 mg/mL, G~c^2.3 (Piechocka). Fiber Lp debated (protofibril~sub-µm vs fiber~tens µm)."),
    "agarose": ECMSpec(
        key="agarose", name="Agarose", is_fibrillar=False,
        modulus_band_Pa=(1.0e3, 1.0e5), modulus_kind="E", poisson=0.5,
        E_gel_Pa=3.0e4,
        ref_conc=1.0, ref_conc_unit="kPa-input",
        tissue_by_S=((0.0, "isotropic 3D encapsulation gel; cartilage/chondrocyte mimic (E~kPa-tens of kPa)"),),
        aliases=("agarose gel",),
        sources=("Normand 2000 Biomacromolecules 10.1021/bm005583j",
                 "Pernodet 1997 Electrophoresis 10.1002/elps.1150180111",
                 "Martikainen 2020 Macromolecules 10.1021/acs.macromol.0c00601",
                 "Zignego 2014 JBiomech 10.1016/j.jbiomech.2013.10.051"),
        notes="Polysaccharide thermogel, E ~5-1000 kPa over 0.5-4 %w/v (typ 30 kPa @1%). Structurally a fine "
              "sub-isostatic bundle network (mesh ~0.24µm, Martikainen 2020) but at cell scale mechanically a "
              "near-incompressible elastic CONTINUUM — the athermal Mikado underestimates its stiffness ~60×, "
              "so its modulus is modeled as a continuum input (like PA/Matrigel), validated by indentation."),
    "pa_gel": ECMSpec(
        key="pa_gel", name="Polyacrylamide gel", is_fibrillar=False,
        modulus_band_Pa=(100.0, 5.0e4), modulus_kind="E", poisson=0.45,
        E_gel_Pa=5.0e3,
        ref_conc=5.0, ref_conc_unit="kPa-input",
        tissue_by_S=((0.0, "engineered stiffness substrate (brain 0.1-1, mammary/fat 1, muscle 10, cartilage 100 kPa)"),),
        aliases=("polyacrylamide", "PAA", "PA gel", "acrylamide gel"),
        sources=("KB-1.21", "KB-1.5", "Discher 2005 Science 10.1126/science.1116995",
                 "Tse-Engler 2010 CurrProtocCellBiol 10.1002/0471143030.cb1016s47",
                 "Palchesko 2012 PLoSONE 10.1371/journal.pone.0051499"),
        notes="Linear-elastic NON-fibrillar chemical gel. E tunable 0.1-40 kPa by acrylamide%/bis-ratio. "
              "E is a material input; indentation validates it (continuum lattice)."),
    "hyaluronic_acid": ECMSpec(
        key="hyaluronic_acid", name="Hyaluronic acid gel", is_fibrillar=False,
        modulus_band_Pa=(10.0, 3.0e3), modulus_kind="E", poisson=0.4,
        E_gel_Pa=3.0e2,
        ref_conc=0.3, ref_conc_unit="kPa-input",
        tissue_by_S=((0.0, "brain / soft neural ECM mimic (crosslinked MeHA)"),),
        aliases=("HA", "hyaluronan", "MeHA", "methacrylated HA"),
        sources=("Burdick 2011 AdvMater", "Lou 2018 Biomaterials"),
        notes="Soft flexible polyelectrolyte, tunable E ~10 Pa-3 kPa by crosslink density. Continuum lattice."),
    "matrigel": ECMSpec(
        key="matrigel", name="Matrigel (reconstituted basement membrane)", is_fibrillar=False,
        modulus_band_Pa=(30.0, 900.0), modulus_kind="E", poisson=0.45,
        E_gel_Pa=4.5e2,
        ref_conc=1.0, ref_conc_unit="kPa-input",
        tissue_by_S=((0.0, "epithelial basement membrane (laminin-111 + collagen-IV, isotropic fine mesh)"),),
        aliases=("basement membrane", "rBM", "reconstituted BM", "laminin gel"),
        sources=("Soofi 2009 JStructBiol 10.1016/j.jsb.2009.05.005 (AFM 37°C E=443± Pa)",
                 "Reed 2009 Langmuir 10.1021/la8033098",
                 "Li 2021 PNAS 10.1073/pnas.2022422118", "Fabris 2018 BiophysJ 10.1016/j.bpj.2018.09.020"),
        notes="Soft BM mimic, E ~30-900 Pa (typ 450, Soofi; batch-variable). Fine near-continuum mesh; continuum lattice."),
}


def get_spec(key: str) -> ECMSpec:
    """Look up an :class:`ECMSpec` by key or alias (case-insensitive)."""
    if key in REGISTRY:
        return REGISTRY[key]
    k = key.strip().lower()
    if k in REGISTRY:
        return REGISTRY[k]
    for spec in REGISTRY.values():
        if k == spec.name.lower() or k in (a.lower() for a in spec.aliases):
            return spec
    raise KeyError(f"unknown ECM material {key!r}; known: {list(REGISTRY)}")


# ─────────────────────────────────────────────────────────────────────────────────────────────────
# Alignment: nematic-order-controlled orientation sampling
# ─────────────────────────────────────────────────────────────────────────────────────────────────

def _watson_S_of_kappa(kappa: float, ngrid: int = 4001) -> float:
    """Nematic order S = ⟨P₂(cosθ)⟩ for the Watson axial distribution p(u) ∝ exp(κu²), u=cosθ ∈ [−1,1].
    κ>0 concentrates fibers about the director (prolate); κ=0 → isotropic (S=0)."""
    u = np.linspace(-1.0, 1.0, ngrid)
    w = np.exp(kappa * (u ** 2 - 1.0))              # subtract 1 for numerical stability (κ large)
    _trapz = np.trapezoid if hasattr(np, "trapezoid") else np.trapz   # numpy 2.x renamed trapz→trapezoid
    Z = _trapz(w, u)
    u2 = _trapz(u ** 2 * w, u) / Z
    return 1.5 * u2 - 0.5


def _watson_kappa_for_S(S: float) -> float:
    """Invert S(κ) for the Watson distribution by bisection. S∈[0,~0.99]; returns κ≥0."""
    S = float(np.clip(S, 0.0, 0.985))
    if S < 1e-4:
        return 0.0
    lo, hi = 0.0, 2.0
    while _watson_S_of_kappa(hi) < S and hi < 1e4:
        hi *= 2.0
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if _watson_S_of_kappa(mid) < S:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def _rotation_to(director: np.ndarray) -> np.ndarray:
    """Rotation matrix mapping +z to the unit ``director`` (for placing a pole-referenced sample)."""
    d = np.asarray(director, float)
    d = d / (np.linalg.norm(d) + 1e-12)
    z = np.array([0.0, 0.0, 1.0])
    v = np.cross(z, d)
    s = np.linalg.norm(v)
    c = float(np.dot(z, d))
    if s < 1e-9:
        return np.eye(3) if c > 0 else np.diag([1.0, -1.0, -1.0])
    vx = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    return np.eye(3) + vx + vx @ vx * ((1.0 - c) / (s * s))


def sample_orientations(n: int, S: float, *, director=(0.0, 0.0, 1.0), dim: int = 3,
                        rng: np.random.Generator) -> np.ndarray:
    """Sample ``n`` unit orientation vectors with nematic order parameter ≈ ``S`` about ``director``.

    3D uses the Watson axial distribution p(u)∝exp(κu²) (u=cosθ from the director) with κ solved for
    the target S; 2D uses a von-Mises on the in-plane angle (director projected into the xy-plane).
    S=0 → isotropic; S→1 → a perfectly aligned bundle. Orientations are axial (±v equivalent)."""
    d = np.asarray(director, float)
    d = d / (np.linalg.norm(d) + 1e-12)
    if dim == 2:
        # in-plane angle φ about the director's in-plane projection; p(φ)∝exp(κ cos2(φ−φ0)) (von Mises in 2φ)
        phi0 = np.arctan2(d[1], d[0])
        # S_2D = I1(κ)/I0(κ); invert numerically
        from scipy.special import i0e, i1e

        def S2d(k):
            return i1e(k) / i0e(k)
        Sc = float(np.clip(S, 0.0, 0.985))
        lo, hi = 0.0, 2.0
        while S2d(hi) < Sc and hi < 1e4:
            hi *= 2.0
        k = 0.5 * (lo + hi)
        if Sc >= 1e-4:
            for _ in range(80):
                mid = 0.5 * (lo + hi)
                if S2d(mid) < Sc:
                    lo = mid
                else:
                    hi = mid
            k = 0.5 * (lo + hi)
        else:
            k = 0.0
        # sample 2φ from von Mises(mean=2φ0, κ=k), then φ; axial → collapse to [0,π)
        two_phi = rng.vonmises(2.0 * phi0, k, size=n)
        phi = 0.5 * two_phi
        v = np.stack([np.cos(phi), np.sin(phi), np.zeros(n)], axis=1)
        return v
    # 3D Watson
    kappa = _watson_kappa_for_S(S)
    # inverse-CDF sample of u=cosθ on a grid
    ug = np.linspace(-1.0, 1.0, 4001)
    w = np.exp(kappa * (ug ** 2 - 1.0))
    cdf = np.cumsum(w)
    cdf /= cdf[-1]
    u = np.interp(rng.random(n), cdf, ug)
    theta = np.arccos(np.clip(u, -1.0, 1.0))
    phi = rng.random(n) * 2.0 * np.pi
    st = np.sin(theta)
    local = np.stack([st * np.cos(phi), st * np.sin(phi), np.cos(theta)], axis=1)
    R = _rotation_to(d)
    return local @ R.T


# ─────────────────────────────────────────────────────────────────────────────────────────────────
# ECM network container + builders
# ─────────────────────────────────────────────────────────────────────────────────────────────────


@dataclass(slots=True)
class ECMNetwork:
    """A built ECM: fiber network + crosslinks + axial segment springs + pinned BC + metadata.

    ``seg_i/seg_j/seg_k/seg_rest`` are the per-segment axial springs (finite EA — the higher-fidelity
    alternative to the inextensible reshape). ``xl_*`` are the crosslinks. Both feed
    ``link_spring_kernel`` unchanged (they are the same Hookean bond primitive)."""

    net: FiberNetwork
    xl_i: np.ndarray
    xl_j: np.ndarray
    xl_k: np.ndarray
    xl_rest: np.ndarray
    seg_i: np.ndarray
    seg_j: np.ndarray
    seg_k: np.ndarray
    seg_rest: np.ndarray
    pinned: np.ndarray
    box_lo: np.ndarray
    box_hi: np.ndarray
    dim: int
    material: str
    mesh_size_um: float
    connectivity_z: float
    S_target: float
    S_measured: float
    director: np.ndarray
    concentration: float = 0.0
    meta: dict = field(default_factory=dict)

    @property
    def links(self) -> np.ndarray:
        """All Hookean bonds (crosslinks + axial segments) stacked as (M,2) for the spring kernel."""
        xl = np.stack([self.xl_i, self.xl_j], 1) if self.xl_i.size else np.zeros((0, 2), np.int64)
        sg = np.stack([self.seg_i, self.seg_j], 1) if self.seg_i.size else np.zeros((0, 2), np.int64)
        return np.concatenate([xl, sg], 0).astype(np.int32)

    @property
    def link_k(self) -> np.ndarray:
        return np.concatenate([self.xl_k, self.seg_k]).astype(np.float64)

    @property
    def link_rest(self) -> np.ndarray:
        return np.concatenate([self.xl_rest, self.seg_rest]).astype(np.float64)


def _crosslink_near_contacts(net: FiberNetwork, xl_contact_um: float, k_xl: float):
    """Mikado near-contact crosslinks: pairs of nodes on DIFFERENT fibers within ``xl_contact_um``,
    at most one (nearest) per fiber-pair. Returns (xl_i, xl_j, xl_k, xl_rest)."""
    node_fiber = np.concatenate([np.full(net.fiber_offsets[f + 1] - net.fiber_offsets[f], f)
                                 for f in range(net.n_fibers)])
    tree = cKDTree(net.pos)
    pairs = tree.query_pairs(r=xl_contact_um, output_type="ndarray")
    if pairs.shape[0]:
        pairs = pairs[node_fiber[pairs[:, 0]] != node_fiber[pairs[:, 1]]]
    if pairs.shape[0]:
        keyfp = node_fiber[pairs[:, 0]].astype(np.int64) * net.n_fibers + node_fiber[pairs[:, 1]]
        d = np.linalg.norm(net.pos[pairs[:, 1]] - net.pos[pairs[:, 0]], axis=1)
        order = np.argsort(d)
        seen, keep = set(), []
        for idx in order:
            k = int(keyfp[idx])
            if k in seen:
                continue
            seen.add(k)
            keep.append(idx)
        pairs = pairs[np.array(keep)]
    if not pairs.shape[0]:
        return (np.zeros(0, np.int64), np.zeros(0, np.int64), np.zeros(0), np.zeros(0))
    xi, xj = pairs[:, 0], pairs[:, 1]
    xr = np.linalg.norm(net.pos[xj] - net.pos[xi], axis=1)
    return xi, xj, np.full(xi.shape[0], float(k_xl)), xr


def _segment_springs(net: FiberNetwork, k_seg: float):
    """Per-fiber consecutive-node axial springs (finite EA). Returns (seg_i, seg_j, seg_k, seg_rest)."""
    si = net.segments[:, 0].astype(np.int64)
    sj = net.segments[:, 1].astype(np.int64)
    return si, sj, np.full(si.shape[0], float(k_seg)), np.asarray(net.seg_rest, float)


def _connectivity_z(net: FiberNetwork, xl_i, xl_j) -> float:
    """Mean crosslinks per fiber ⟨z⟩ (each crosslink attaches 2 fibers)."""
    if xl_i.size == 0 or net.n_fibers == 0:
        return 0.0
    return float(2.0 * xl_i.size / net.n_fibers)


def _mesh_size(net: FiberNetwork, xl_i, xl_j, seg_um, fiber_len_um) -> float:
    """Emergent pore size ξ = median inter-crosslink spacing along fibers [µm]."""
    if xl_i.size == 0:
        return fiber_len_um
    xl_nodes = np.concatenate([xl_i, xl_j])
    spac = []
    for f in range(net.n_fibers):
        a0, b0 = net.fiber_offsets[f], net.fiber_offsets[f + 1]
        on = np.sort(xl_nodes[(xl_nodes >= a0) & (xl_nodes < b0)])
        if on.size >= 2:
            spac.extend(np.diff(on) * seg_um)
    return float(np.median(spac)) if spac else fiber_len_um


def fibers_for_concentration(spec: ECMSpec, box_lo, box_hi, conc: float, dim: int) -> int:
    """Number of rods to hit the concentration-implied mesh ξ(conc). ξ = ref_mesh·(conc/ref_conc)^(−exp).
    3D: fiber length density ρ_L = ξ^(−2) [µm⁻²], N = ρ_L·V/L_f. 2D: areal length density ρ_A = ξ^(−1)
    [µm⁻¹] over a thin slab, N = ρ_A·A/L_f."""
    box_lo = np.asarray(box_lo, float)
    box_hi = np.asarray(box_hi, float)
    xi = spec.ref_mesh_um * (conc / spec.ref_conc) ** (-spec.conc_exponent)
    if dim == 2:
        Lx, Ly = (box_hi - box_lo)[:2]
        rho_A = 1.0 / xi                      # µm of fiber per µm² of sheet
        return max(1, int(round(rho_A * Lx * Ly / spec.fiber_len_um)))
    V = float(np.prod(box_hi - box_lo))
    rho_L = 1.0 / xi ** 2                      # µm of fiber per µm³
    return max(1, int(round(rho_L * V / spec.fiber_len_um)))


def _random_rods(box_lo, box_hi, n_fibers, length, orient, seg_um, dim, rng):
    """Build ``n_fibers`` discretised straight rods with the given per-fiber ``orient`` unit axes.
    Centres uniform in the box (2D: z at the mid-plane); endpoints clipped into the box."""
    box_lo = np.asarray(box_lo, float)
    box_hi = np.asarray(box_hi, float)
    nb = max(2, int(round(length / seg_um)) + 1)
    zc = 0.5 * (box_lo[2] + box_hi[2])
    fibers = []
    for i in range(n_fibers):
        c = box_lo + rng.random(3) * (box_hi - box_lo)
        if dim == 2:
            c[2] = zc
        v = orient[i]
        if dim == 2:
            v = np.array([v[0], v[1], 0.0])
        v = v / (np.linalg.norm(v) + 1e-12)
        a = np.clip(c - 0.5 * length * v, box_lo, box_hi)
        b = np.clip(c + 0.5 * length * v, box_lo, box_hi)
        if np.linalg.norm(b - a) < seg_um:
            continue
        t = np.linspace(0.0, 1.0, nb)[:, None]
        fibers.append(a[None, :] * (1 - t) + b[None, :] * t)
    return fibers


def _pin_mask(net: FiberNetwork, box_lo, box_hi, pin_faces, margin_um) -> np.ndarray:
    """Dirichlet-pin nodes within ``margin_um`` of the named box faces (embedded in bulk matrix / clamped
    to a coverslip). ``pin_faces`` is a tuple of {"x_lo","x_hi","y_lo","y_hi","z_lo","z_hi"}."""
    pinned = np.zeros(net.n_nodes, bool)
    axmap = {"x": 0, "y": 1, "z": 2}
    for face in pin_faces:
        ax = axmap[face[0]]
        if face.endswith("lo"):
            pinned |= net.pos[:, ax] < box_lo[ax] + margin_um
        else:
            pinned |= net.pos[:, ax] > box_hi[ax] - margin_um
    return pinned


def build_fibrillar_ecm(spec: ECMSpec, box_lo, box_hi, *, n_fibers: int | None = None,
                        concentration: float | None = None, dim: int = 3, alignment_S: float = 0.0,
                        director=(1.0, 0.0, 0.0), pin_faces=("z_lo",), pin_margin_um: float = 1.0,
                        rng: np.random.Generator | None = None) -> ECMNetwork:
    """Build a fibrillar ECM (Mikado semiflexible-fiber network) for ``spec``.

    Density is set either by ``n_fibers`` directly or by ``concentration`` (→ mesh ξ → N). Orientation
    is sampled at nematic order ``alignment_S`` about ``director`` (0 = isotropic). ``dim`` is 2 (thin
    planar sheet) or 3 (embedded slab). The macroscopic modulus is NOT set here — it emerges from this
    microstructure and is measured by ``ecm_mechanics``.
    """
    if not spec.is_fibrillar:
        raise ValueError(f"{spec.key} is a continuum gel; use build_continuum_ecm")
    if rng is None:
        rng = np.random.default_rng(0)
    box_lo = np.asarray(box_lo, float)
    box_hi = np.asarray(box_hi, float)
    if n_fibers is None:
        if concentration is None:
            concentration = spec.ref_conc
        n_fibers = fibers_for_concentration(spec, box_lo, box_hi, concentration, dim)
    orient = sample_orientations(n_fibers, alignment_S, director=director, dim=dim, rng=rng)
    fibers = _random_rods(box_lo, box_hi, n_fibers, spec.fiber_len_um, orient, spec.seg_um, dim, rng)
    net = build_fiber_network(fibers, kappa=spec.kappa_pN_um2)
    xl_i, xl_j, xl_k, xl_rest = _crosslink_near_contacts(net, spec.xl_contact_um, spec.k_xl_pN_um)
    seg_i, seg_j, seg_k, seg_rest = _segment_springs(net, spec.k_seg_pN_um())
    pinned = _pin_mask(net, box_lo, box_hi, pin_faces, pin_margin_um)
    z = _connectivity_z(net, xl_i, xl_j)
    mesh = _mesh_size(net, xl_i, xl_j, spec.seg_um, spec.fiber_len_um)
    S_meas = _safe_order(net)
    return ECMNetwork(net=net, xl_i=xl_i, xl_j=xl_j, xl_k=xl_k, xl_rest=xl_rest,
                      seg_i=seg_i, seg_j=seg_j, seg_k=seg_k, seg_rest=seg_rest,
                      pinned=pinned, box_lo=box_lo, box_hi=box_hi, dim=dim, material=spec.key,
                      mesh_size_um=mesh, connectivity_z=z, S_target=alignment_S, S_measured=S_meas,
                      director=np.asarray(director, float), concentration=concentration or 0.0,
                      meta={"n_fibers": net.n_fibers, "n_nodes": net.n_nodes, "n_xl": int(xl_i.size),
                            "kappa": spec.kappa_pN_um2, "EA_pN": spec.EA_pN})


def _safe_order(net: FiberNetwork) -> float:
    from ffn_sim.ff.architecture_metrics import parallel_order_parameter
    try:
        return parallel_order_parameter(net) if net.n_fibers > 1 else 0.0
    except Exception:
        return float("nan")


def build_continuum_ecm(spec: ECMSpec, box_lo, box_hi, *, dim: int = 3, node_spacing_um: float = 2.0,
                        E_gel_Pa: float | None = None, k_bond_pN_um: float | None = None,
                        pin_faces=("z_lo",), pin_margin_um: float = 1.0,
                        rng: np.random.Generator | None = None) -> ECMNetwork:
    """Build a continuum-gel ECM (PAA / HA / Matrigel) as a coarse random elastic spring network.

    A linear-elastic isotropic solid discretised as a jittered lattice of nodes connected by
    nearest-neighbour central-force springs (a random spring network; Delaunay in the interior). The
    bond stiffness is set from the target modulus by the Cauchy-Born affine estimate for an isotropic
    random spring network — E ≈ (1/6)·k·⟨z⟩·(N/V)·ℓ² for a 3D network — but because that estimate is
    approximate, ``ecm_mechanics.calibrate_continuum_k`` refines it against the measured shear/uniaxial
    modulus once; here we accept an explicit ``k_bond_pN_um`` override (the calibrated value) when given.
    E is a material INPUT (literature), not an emergent quantity — the honest model for a chemical gel.
    """
    if spec.is_fibrillar:
        raise ValueError(f"{spec.key} is fibrillar; use build_fibrillar_ecm")
    if rng is None:
        rng = np.random.default_rng(0)
    box_lo = np.asarray(box_lo, float)
    box_hi = np.asarray(box_hi, float)
    E = spec.E_gel_Pa if E_gel_Pa is None else E_gel_Pa
    pts, _ = _lattice_points(box_lo, box_hi, node_spacing_um, dim, rng)
    # node-only carrier network (no bending); the Delaunay springs below ARE the elasticity
    net = _point_network(pts)
    bonds = _delaunay_bonds(pts, dim)
    rest = np.linalg.norm(pts[bonds[:, 1]] - pts[bonds[:, 0]], axis=1)
    # affine estimate of E per unit k, refined by override
    if k_bond_pN_um is None:
        mean_l = float(rest.mean())
        z_bar = 2.0 * len(bonds) / max(len(pts), 1)
        n_over_v = len(pts) / float(np.prod(box_hi - box_lo)) if dim == 3 else \
            len(pts) / float(np.prod((box_hi - box_lo)[:2]) * node_spacing_um)
        E_per_k = (1.0 / 6.0) * z_bar * n_over_v * mean_l ** 2      # Pa per (pN/µm)
        k_bond = E / max(E_per_k, 1e-9)
    else:
        k_bond = k_bond_pN_um
    seg_i, seg_j = bonds[:, 0], bonds[:, 1]
    seg_k = np.full(len(bonds), float(k_bond))
    pinned = _pin_mask(net, box_lo, box_hi, pin_faces, pin_margin_um)
    return ECMNetwork(net=net, xl_i=np.zeros(0, np.int64), xl_j=np.zeros(0, np.int64),
                      xl_k=np.zeros(0), xl_rest=np.zeros(0),
                      seg_i=seg_i, seg_j=seg_j, seg_k=seg_k, seg_rest=rest,
                      pinned=pinned, box_lo=box_lo, box_hi=box_hi, dim=dim, material=spec.key,
                      mesh_size_um=node_spacing_um, connectivity_z=2.0 * len(bonds) / max(len(pts), 1),
                      S_target=0.0, S_measured=0.0, director=np.array([1.0, 0.0, 0.0]),
                      concentration=0.0,
                      meta={"n_nodes": len(pts), "n_bonds": int(len(bonds)), "E_gel_Pa": E,
                            "k_bond_pN_um": float(k_bond), "continuum": True})


def _lattice_points(box_lo, box_hi, spacing, dim, rng, jitter=0.3):
    """Jittered regular grid of points filling the box (jitter breaks lattice anisotropy → isotropic)."""
    box_lo = np.asarray(box_lo, float)
    box_hi = np.asarray(box_hi, float)
    axes = []
    for a in range(3):
        n = max(2, int(round((box_hi[a] - box_lo[a]) / spacing)) + 1)
        if dim == 2 and a == 2:
            axes.append(np.array([0.5 * (box_lo[2] + box_hi[2])]))
        else:
            axes.append(np.linspace(box_lo[a], box_hi[a], n))
    gx, gy, gz = np.meshgrid(*axes, indexing="ij")
    pts = np.stack([gx.ravel(), gy.ravel(), gz.ravel()], 1)
    interior = np.ones(len(pts), bool)
    for a in range(dim):
        interior &= (pts[:, a] > box_lo[a] + 1e-6) & (pts[:, a] < box_hi[a] - 1e-6)
    pts[interior] += (rng.random((interior.sum(), 3)) - 0.5) * jitter * spacing * (
        np.array([1, 1, 1]) if dim == 3 else np.array([1, 1, 0]))
    return pts, None


def _point_network(pts) -> FiberNetwork:
    """A FiberNetwork of isolated point-nodes (no segments/bending) — the continuum-lattice carrier."""
    return FiberNetwork(pos=np.ascontiguousarray(pts, np.float64),
                        fiber_offsets=np.arange(len(pts) + 1, dtype=np.int64),
                        segments=np.zeros((0, 2), np.int64), bend_triples=np.zeros((0, 3), np.int64),
                        seg_rest=np.zeros(0), kappa=np.zeros(len(pts)))


def _delaunay_bonds(pts, dim) -> np.ndarray:
    """Unique edges of the Delaunay triangulation (the random-spring-network connectivity)."""
    from scipy.spatial import Delaunay
    P = pts[:, :2] if dim == 2 else pts
    tri = Delaunay(P)
    edges = set()
    simp = tri.simplices
    for s in simp:
        for a in range(len(s)):
            for b in range(a + 1, len(s)):
                i, j = int(s[a]), int(s[b])
                edges.add((i, j) if i < j else (j, i))
    return np.array(sorted(edges), np.int64)


def build_ecm(material: str, box_lo, box_hi, **kw) -> ECMNetwork:
    """Dispatch: build a fibrillar or continuum ECM for ``material`` (registry key/alias)."""
    spec = get_spec(material)
    if spec.is_fibrillar:
        return build_fibrillar_ecm(spec, box_lo, box_hi, **kw)
    # continuum gels ignore fibrillar-only kwargs
    for k in ("n_fibers", "concentration", "alignment_S", "director"):
        kw.pop(k, None)
    return build_continuum_ecm(spec, box_lo, box_hi, **kw)


def build_composite(components, box_lo, box_hi, *, dim: int = 3, interlink_um: float = 0.0,
                    interlink_k: float = 100.0, pin_faces=("z_lo",), pin_margin_um: float = 1.0,
                    rng: np.random.Generator | None = None) -> ECMNetwork:
    """Build a MIXED / interpenetrating ECM from several components.

    ``components`` is a list of dicts, each ``{"material": key, ...builder kwargs}`` (e.g.
    ``{"material":"collagen_I","concentration":1.5,"alignment_S":0.4}``). Networks are built
    independently in the same box, stacked with ``concat_fiber_networks``, and (optionally) inter-
    crosslinked at near-contacts within ``interlink_um`` (interpenetrating-network coupling, e.g.
    collagen fibers entangled with a Matrigel lattice). Returns one merged ECMNetwork.
    """
    if rng is None:
        rng = np.random.default_rng(0)
    subs = []
    for comp in components:
        comp = dict(comp)
        mat = comp.pop("material")
        comp.setdefault("dim", dim)
        comp.setdefault("pin_faces", ())            # pin globally after merge
        sub = build_ecm(mat, box_lo, box_hi, rng=rng, **comp)
        subs.append(sub)
    merged_net, node_off = concat_fiber_networks([s.net for s in subs])
    # remap and concatenate each sub's crosslink + segment bonds into the merged index space
    xl_i, xl_j, xl_k, xl_rest = [], [], [], []
    sg_i, sg_j, sg_k, sg_rest = [], [], [], []
    for s, base in zip(subs, node_off[:-1]):
        if s.xl_i.size:
            xl_i.append(s.xl_i + base); xl_j.append(s.xl_j + base); xl_k.append(s.xl_k); xl_rest.append(s.xl_rest)
        if s.seg_i.size:
            sg_i.append(s.seg_i + base); sg_j.append(s.seg_j + base); sg_k.append(s.seg_k); sg_rest.append(s.seg_rest)
    xi = np.concatenate(xl_i) if xl_i else np.zeros(0, np.int64)
    xj = np.concatenate(xl_j) if xl_j else np.zeros(0, np.int64)
    xk = np.concatenate(xl_k) if xl_k else np.zeros(0)
    xr = np.concatenate(xl_rest) if xl_rest else np.zeros(0)
    si = np.concatenate(sg_i) if sg_i else np.zeros(0, np.int64)
    sj = np.concatenate(sg_j) if sg_j else np.zeros(0, np.int64)
    sk = np.concatenate(sg_k) if sg_k else np.zeros(0)
    sr = np.concatenate(sg_rest) if sg_rest else np.zeros(0)
    # inter-network crosslinks (couple different components at near-contact)
    if interlink_um > 0.0 and len(subs) > 1:
        comp_of = np.concatenate([np.full(node_off[k + 1] - node_off[k], k) for k in range(len(subs))])
        tree = cKDTree(merged_net.pos)
        pairs = tree.query_pairs(r=interlink_um, output_type="ndarray")
        if pairs.shape[0]:
            pairs = pairs[comp_of[pairs[:, 0]] != comp_of[pairs[:, 1]]]
        if pairs.shape[0]:
            ii, jj = pairs[:, 0], pairs[:, 1]
            rr = np.linalg.norm(merged_net.pos[jj] - merged_net.pos[ii], axis=1)
            xi = np.concatenate([xi, ii]); xj = np.concatenate([xj, jj])
            xk = np.concatenate([xk, np.full(ii.size, interlink_k)]); xr = np.concatenate([xr, rr])
    pinned = _pin_mask(merged_net, np.asarray(box_lo, float), np.asarray(box_hi, float), pin_faces, pin_margin_um)
    z = _connectivity_z(merged_net, xi, xj)
    mesh = float(np.median([s.mesh_size_um for s in subs]))
    S_meas = _safe_order(merged_net)
    return ECMNetwork(net=merged_net, xl_i=xi, xl_j=xj, xl_k=xk, xl_rest=xr,
                      seg_i=si, seg_j=sj, seg_k=sk, seg_rest=sr, pinned=pinned,
                      box_lo=np.asarray(box_lo, float), box_hi=np.asarray(box_hi, float), dim=dim,
                      material="+".join(c["material"] for c in components),
                      mesh_size_um=mesh, connectivity_z=z, S_target=float("nan"), S_measured=S_meas,
                      director=np.array([1.0, 0.0, 0.0]),
                      meta={"components": [c["material"] for c in components], "n_nodes": merged_net.n_nodes})


__all__ = ["ECMSpec", "REGISTRY", "get_spec", "ECMNetwork", "sample_orientations",
           "build_fibrillar_ecm", "build_continuum_ecm", "build_ecm", "build_composite",
           "fibers_for_concentration"]
