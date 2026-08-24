r"""Region builders on the membrane-adjacent manifold — one dispatcher for the whole region set — I4 / P3.

``build_region(spec, rng)`` turns one :class:`RegionSpec` into a :class:`~ac.weave.woven_cell.RegionLeaf`.
Dispatch by manifold:

  * ``"sphere"`` (CORTEX) -> DELEGATE to ``ff.weave.weave(CORTEX)`` unchanged. The sphere seed is already the
    authoritative isotropic shell, and delegating (with untouched RNG order) is what makes Gate-1 bit-identical.
  * ``"patch"`` (LAMELLIPODIUM) -> the dendritic Arp2/3 rebuild (``ac.weave.lamellipodium``).
  * ``"bundle"`` (ventral/dorsal SF, transverse arc, perinuclear cap, filopodium) -> an EMERGENT
    isotropic-conditional seed on the region manifold (filaments placed on the declared geometry with
    orientations sampled isotropically-conditional, NOT a pre-made sarcomeric bundle). ``emergent=False`` on a
    bundle region reproduces the DEMOTED scripted ``ff.weave._build_bundle`` as a labeled control only.

HONEST SCOPE (I4). Cortex = FULL FIDELITY (bit-identical parity). Lamellipodium = FULL FIDELITY (dendritic
Arp2/3, angle-harmonic branch, N-fixed activation). Ventral SF / dorsal SF / transverse arc / perinuclear cap /
filopodium = honest REGION BUILDERS / SEEDS: they place the correct manifold, seed filaments isotropic-
conditional, place alpha-actinin + NMII + FA/LINC anchor SITES, and concat correctly into the one network — but
the bundle CONDENSATION (SF/arc/cap/filopodium emerging from the isotropic seed) is the I5 NATIVE proof, NOT
claimed here. The cap->LINC->nucleus tether stubs (``linc_sites``) are the I7 load-path seed; the nucleus-side
endpoints are supplied by ``ac.nucleus`` at I7 (consumed read-only, not built here).

engine units: length um, stiffness pN/um.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import numpy.typing as npt

from aleph.laws.architecture_spec import ArchitectureSpec, CORTEX, LAMELLIPODIUM
# Both crosslinker Hands ORIGINATE in `ff.hand_kmc`; `gamma_floor` only re-exported them, and reaching
# them through it dragged a 684-line module — whose own full-native results are retired — plus, through
# its function-local `ff.relax` import, the nine-file parked `dcm/` engine into this run's build
# verification. Same objects, same values; only the path changed (2026-07-28).
from aleph.laws.fiber_network import cross_fiber_pairs as _cross_fiber_pairs
from aleph.laws.hand_kmc import ALPHA_ACTININ, FILAMIN
from aleph.components.weave.branch_angle import ARP23_K_THETA, ARP23_THETA0_RAD, sample_branch_angles
from aleph.components.weave.lamellipodium import build_lamellipodium
from aleph.components.weave.woven_cell import RegionLeaf

__all__ = [
    "RegionSpec",
    "build_region",
    "CORTEX_REGION",
    "LAMELLIPODIUM_REGION",
    "VENTRAL_SF_REGION",
    "DORSAL_SF_REGION",
    "TRANSVERSE_ARC_REGION",
    "PERINUCLEAR_CAP_REGION",
    "FILOPODIUM_REGION",
    "FULL_REGION_SET",
    "derive_cortex_arp23_split",
    "cortex_arp23_region",
]


@dataclass(frozen=True, slots=True)
class RegionSpec:
    """One region to weave: an ff ``ArchitectureSpec`` + placement + emergence flags + anchor role.

    Attributes:
        arch: the ff architecture row (manifold / filament / crosslinker / motor).
        manifold_kind: the placement manifold for a bundle region
            (``"basal_footprint"`` | ``"dome"`` | ``"arc"`` | ``"finger"``); ignored for sphere/patch.
        centre_um: manifold centre [um].
        scale_um: manifold characteristic size [um] (footprint radius / dome radius / arc radius / finger len).
        z_basal_um: basal plane z for FA sites [um].
        emergent: bundle regions seed isotropic-conditional (True) or the demoted scripted bundle (False).
        seed_concentration: von-Mises-Fisher concentration of the isotropic-conditional orientation
            (0 = fully isotropic; the emergence-honoring default is small).
        anchor_role: ``"both_ends"`` (ventral SF FA-FA) | ``"one_end"`` (dorsal SF) | ``"none"`` (arc) |
            ``"linc"`` (perinuclear cap) | ``"tip"`` (filopodium).
        n_daughter_pool: lamellipodium dormant daughter capacity (patch only).
        npf_kappa: lamellipodium NPF orientation bias (patch only; 0 = isotropic).
    """

    arch: ArchitectureSpec
    manifold_kind: str = "basal_footprint"
    centre_um: tuple[float, float, float] = (0.0, 0.0, 0.0)
    scale_um: float = 7.5
    z_basal_um: float = 0.0
    emergent: bool = True
    seed_concentration: float = 0.0
    anchor_role: str = "none"
    n_daughter_pool: int = 0
    npf_kappa: float = 0.0

    @property
    def name(self) -> str:
        return self.arch.name


# ── the region table (I4 full region set) ─────────────────────────────────────────────────────────
CORTEX_REGION = RegionSpec(arch=CORTEX)
LAMELLIPODIUM_REGION = RegionSpec(
    arch=LAMELLIPODIUM, centre_um=(0.0, 7.5, 0.0), scale_um=8.0, n_daughter_pool=100, npf_kappa=0.0)


def _sf_arch(name: str, n_fil: int, length: float) -> ArchitectureSpec:
    """A ventral/dorsal-SF-like architecture row: formin bundle, alpha-actinin, NMIIA, mixed polarity."""
    from aleph.laws.architecture_spec import CrosslinkerSpec, FilamentSpec, MotorSpec
    from aleph.laws.hand_kmc import NMIIA_MYOSIN
    return ArchitectureSpec(
        name=name, manifold="bundle", R_um=length,
        filament=FilamentSpec(nucleator="formin", length_um=length, seg_um=0.5, n_filaments=n_fil,
                              orientation="parallel", polarity="mixed"),
        crosslinker=CrosslinkerSpec(hand=ALPHA_ACTININ, bind_mode="parallel", density_per_fil=2.0),
        motor=MotorSpec(hand=NMIIA_MYOSIN, mode="sarcomeric", density_per_fil=1.0))


VENTRAL_SF_REGION = RegionSpec(
    arch=_sf_arch("ventral_sf", 20, 10.0), manifold_kind="basal_footprint", scale_um=6.0,
    z_basal_um=0.0, anchor_role="both_ends")
DORSAL_SF_REGION = RegionSpec(
    arch=_sf_arch("dorsal_sf", 15, 8.0), manifold_kind="basal_footprint", scale_um=6.0,
    centre_um=(0.0, 4.0, 0.0), z_basal_um=0.0, anchor_role="one_end")
TRANSVERSE_ARC_REGION = RegionSpec(
    arch=_sf_arch("transverse_arc", 15, 8.0), manifold_kind="arc", scale_um=6.0,
    centre_um=(0.0, 5.0, 1.0), anchor_role="none")
PERINUCLEAR_CAP_REGION = RegionSpec(
    arch=_sf_arch("perinuclear_cap", 12, 9.0), manifold_kind="dome", scale_um=3.5,
    centre_um=(0.0, 0.0, 0.0), anchor_role="linc")


def _filopodium_arch() -> ArchitectureSpec:
    from aleph.laws.architecture_spec import CrosslinkerSpec, FilamentSpec, MotorSpec
    return ArchitectureSpec(
        name="filopodium", manifold="bundle", R_um=3.0,
        filament=FilamentSpec(nucleator="formin", length_um=3.0, seg_um=0.5, n_filaments=20,
                              orientation="parallel", polarity="uniform"),
        crosslinker=CrosslinkerSpec(hand=FILAMIN, bind_mode="parallel", density_per_fil=3.0),
        motor=MotorSpec(hand=None, mode="none"))


FILOPODIUM_REGION = RegionSpec(
    arch=_filopodium_arch(), manifold_kind="finger", scale_um=3.0, centre_um=(0.0, 8.0, 0.0),
    anchor_role="tip")

FULL_REGION_SET = [
    CORTEX_REGION, LAMELLIPODIUM_REGION, VENTRAL_SF_REGION, DORSAL_SF_REGION,
    TRANSVERSE_ARC_REGION, PERINUCLEAR_CAP_REGION, FILOPODIUM_REGION,
]


# ── manifold samplers (base point + local mean axis) ─────────────────────────────────────────────
def _sample_manifold(kind: str, centre, scale: float, z0: float, n: int, rng: np.random.Generator):
    """Sample ``n`` (base_point, local_axis) pairs on a region manifold. The local axis is the DECLARED
    geometric orientation the isotropic-conditional seed is blended toward (SEEDED, ablatable)."""
    c = np.asarray(centre, float)
    if kind == "basal_footprint":                              # disk at z0 (ventral/dorsal SF footprint)
        r = scale * np.sqrt(rng.random(n))
        a = rng.uniform(0, 2 * np.pi, n)
        base = np.column_stack([c[0] + r * np.cos(a), c[1] + r * np.sin(a), np.full(n, z0)])
        axis = np.tile(np.array([1.0, 0.0, 0.0]), (n, 1))      # in-plane declared axis
    elif kind == "dome":                                       # hemispherical cap over the nucleus
        u = rng.random(n)
        phi = np.arccos(1.0 - 0.6 * u)                          # apical band (theta small)
        az = rng.uniform(0, 2 * np.pi, n)
        base = np.column_stack([c[0] + scale * np.sin(phi) * np.cos(az),
                                c[1] + scale * np.sin(phi) * np.sin(az),
                                c[2] + scale * np.cos(phi)])
        # tangent-to-dome meridional axis (cap bundles arc OVER the nucleus)
        axis = np.column_stack([np.cos(phi) * np.cos(az), np.cos(phi) * np.sin(az), -np.sin(phi)])
    elif kind == "arc":                                        # curved band parallel to the leading edge
        t = rng.uniform(-1.0, 1.0, n)
        ang = t * (np.pi / 3.0)
        base = np.column_stack([c[0] + scale * np.sin(ang), c[1] + scale * (1.0 - np.cos(ang)), np.full(n, c[2])])
        axis = np.column_stack([np.cos(ang), np.sin(ang), np.zeros(n)])   # tangent to the arc
    elif kind == "finger":                                     # elongated protrusion (filopodium)
        base = np.column_stack([c[0] + rng.uniform(-0.05, 0.05, n), np.full(n, c[1]),
                                c[2] + rng.uniform(-0.05, 0.05, n)])
        axis = np.tile(np.array([0.0, 1.0, 0.0]), (n, 1))      # protrusion axis
    else:
        raise ValueError(f"unknown manifold_kind {kind!r}")
    return base, axis


def _isotropic_conditional_dirs(axis: npt.NDArray[np.float64], kappa: float, rng: np.random.Generator):
    """Directions sampled isotropic-conditional: a random unit vector blended toward ``axis`` by ``kappa``.

    ``kappa = 0`` => fully isotropic (the bundle must CONDENSE, not be scripted); larger ``kappa`` => a mild
    declared bias toward the manifold axis (the SEEDED, separately-ablatable orientation prior of §5).
    """
    n = axis.shape[0]
    v = rng.standard_normal((n, 3))
    v /= np.linalg.norm(v, axis=1, keepdims=True) + 1e-12
    if kappa <= 0.0:
        return v
    d = axis + kappa * v
    return d / (np.linalg.norm(d, axis=1, keepdims=True) + 1e-12)


def _build_emergent_bundle(spec: RegionSpec, rng: np.random.Generator) -> RegionLeaf:
    """Isotropic-conditional seed of a bundle region on its manifold (SF / arc / cap / filopodium)."""
    from aleph.laws.fiber_network import build_fiber_network
    from aleph.laws import units as U

    fs = spec.arch.filament
    n_fil = fs.n_filaments
    nb = max(2, int(round(fs.length_um / fs.seg_um)) + 1)
    base, axis = _sample_manifold(spec.manifold_kind, spec.centre_um, spec.scale_um, spec.z_basal_um, n_fil, rng)
    dirs = _isotropic_conditional_dirs(axis, spec.seed_concentration, rng)
    fibers = [base[f][None, :] + np.arange(nb)[:, None] * fs.seg_um * dirs[f][None, :] for f in range(n_fil)]
    net = build_fiber_network(fibers, kappa=U.KAPPA_ACTIN)

    # seed-scale mesoscale reach ~ 1.5 * inter-filament spacing (the coarse seed is far sparser than the
    # native full-density bundle; the physiological fascin/alpha-actinin reach applies at the native run).
    reach = max(1.5 * spec.scale_um / max(np.sqrt(n_fil), 1.0), 0.6)
    n_xl = int(round(n_fil * spec.arch.crosslinker.density_per_fil))
    xl_pairs = _cross_fiber_pairs(net, reach, n_xl, rng)
    nxl = xl_pairs.shape[0]
    xl_i, xl_j = (xl_pairs[:, 0], xl_pairs[:, 1]) if nxl else (np.zeros(0, np.int64), np.zeros(0, np.int64))
    xl_k = np.full(nxl, spec.arch.crosslinker.hand.link_k)
    xl_rest = np.linalg.norm(net.pos[xl_j] - net.pos[xl_i], axis=1) if nxl else np.zeros(0)

    if spec.arch.motor.hand is not None:
        n_myo = int(round(n_fil * spec.arch.motor.density_per_fil))
        used = {(int(a), int(b)) for a, b in xl_pairs}
        myo_pairs = _cross_fiber_pairs(net, reach, n_myo, rng, exclude=used)
    else:
        myo_pairs = np.zeros((0, 2), np.int64)
    myo_i = myo_pairs[:, 0] if myo_pairs.size else np.zeros(0, np.int64)
    myo_j = myo_pairs[:, 1] if myo_pairs.size else np.zeros(0, np.int64)

    off = net.fiber_offsets
    ends0 = off[:-1].astype(np.int64)
    ends1 = (off[1:] - 1).astype(np.int64)
    fa_sites = np.zeros(0, np.int64)
    linc_sites = np.zeros(0, np.int64)
    if spec.anchor_role == "both_ends":
        fa_sites = np.concatenate([ends0, ends1])
    elif spec.anchor_role in ("one_end", "tip"):
        fa_sites = ends0 if spec.anchor_role == "one_end" else ends1
    elif spec.anchor_role == "linc":
        # LINC tether stubs = the fiber node nearest the nucleus centre (dome apex side)
        c = np.asarray(spec.centre_um, float)
        d = np.linalg.norm(net.pos - c, axis=1)
        linc_sites = np.array([int(off[f] + np.argmin(d[off[f]:off[f + 1]])) for f in range(n_fil)], np.int64)

    polarity = np.where(np.arange(n_fil) % 2 == 0, 1, -1).astype(np.int64) if fs.polarity == "mixed" \
        else np.ones(n_fil, np.int64)

    return RegionLeaf(
        name=spec.name, manifold=spec.manifold_kind, nucleator=fs.nucleator, pos=net.pos.copy(),
        fiber_offsets=net.fiber_offsets.astype(np.int64), xl_i=xl_i.astype(np.int64), xl_j=xl_j.astype(np.int64),
        xl_k=xl_k, xl_rest=xl_rest, myo_i=myo_i.astype(np.int64), myo_j=myo_j.astype(np.int64),
        polarity=polarity, provenance="isotropic_conditional_seed", fa_sites=fa_sites, linc_sites=linc_sites)


def _build_cortex_delegated(spec: RegionSpec, rng: np.random.Generator,
                            overlap_free: bool = False,
                            overlap_mode: str = "transverse", overlap_span: int = 2) -> RegionLeaf:
    """Delegate to ``ff.weave.weave(CORTEX)`` (Gate-1 parity when overlap_free=False)."""
    from aleph.laws.weave import weave
    cc = weave(spec.arch, rng=rng, overlap_free=overlap_free,
               overlap_mode=overlap_mode, overlap_span=overlap_span)
    n_fib = cc.net.fiber_offsets.shape[0] - 1
    return RegionLeaf(
        name=spec.name, manifold="sphere", nucleator=spec.arch.filament.nucleator, pos=cc.net.pos.copy(),
        fiber_offsets=cc.net.fiber_offsets.astype(np.int64),
        xl_i=cc.xl_i.astype(np.int64), xl_j=cc.xl_j.astype(np.int64), xl_k=cc.xl_k, xl_rest=cc.xl_rest,
        myo_i=cc.myo_i.astype(np.int64), myo_j=cc.myo_j.astype(np.int64),
        polarity=np.ones(n_fib, np.int64), provenance="delegated_ff_weave",
        branch_triples=cc.branch_triples.astype(np.int64))


def _build_lamellipodium_region(spec: RegionSpec, rng: np.random.Generator) -> RegionLeaf:
    """Dendritic Arp2/3 rebuild -> RegionLeaf (crosslinks = the branch anchors; NMII off in the core)."""
    fs = spec.arch.filament
    seed = build_lamellipodium(
        n_mothers=fs.n_filaments, n_daughter_pool=spec.n_daughter_pool, length_um=fs.length_um,
        seg_um=fs.seg_um, patch_half_um=spec.scale_um, rng=rng,
        theta0=fs.branch_angle_rad or ARP23_THETA0_RAD, k_theta=ARP23_K_THETA, npf_kappa=spec.npf_kappa)
    # crosslinks = the stiff daughter-base <-> branch anchors (force-free at the branch geometry)
    a = seed.branch_anchors
    if a.size:
        xl_i, xl_j = a[:, 0].astype(np.int64), a[:, 1].astype(np.int64)
        xl_k = np.full(a.shape[0], ALPHA_ACTININ.link_k)
        xl_rest = np.linalg.norm(seed.pos[xl_j] - seed.pos[xl_i], axis=1)
    else:
        xl_i = xl_j = np.zeros(0, np.int64); xl_k = xl_rest = np.zeros(0)
    return RegionLeaf(
        name=spec.name, manifold="patch", nucleator="arp23", pos=seed.pos, fiber_offsets=seed.fiber_offsets,
        xl_i=xl_i, xl_j=xl_j, xl_k=xl_k, xl_rest=xl_rest, myo_i=np.zeros(0, np.int64),
        myo_j=np.zeros(0, np.int64), polarity=seed.polarity, provenance="isotropic_conditional_seed",
        branch_triples=seed.branch_triples, branch_anchors=seed.branch_anchors,
        branch_active=seed.branch_active, active_mask=seed.active_mask)


# ── mixed formin + Arp2/3 cortex sub-population (FLAGGED prototype, 2026-07-24) ────────────────────
# Native cortex is NOT formin-only: it is ~2/3 formin (long) + ~1/3 Arp2/3-nucleated (short, branched) BY MASS
# (Bovellan 2014, KB-3.18). The default cortex (``architecture_spec.CORTEX``) is formin-only; this adds the
# missing Arp2/3 branched sub-population as an OPT-IN split of the fixed cortex filament BUDGET, so total areal
# density (~100/µm²) is preserved. The Arp2/3 sub-population REUSES the Arp2/3 angle-harmonic branch path
# (``branch_angle.sample_branch_angles`` + ``ARP23_THETA0_RAD``/``ARP23_K_THETA``; the ``branch_triples`` it
# emits feed the SAME device ``branch_angle_kernel`` the lamellipodium uses) — only the sphere-native placement
# is new. See docs/v2_audit/CORTEX_ARP23_POPULATION_2026-07-24.md.


def derive_cortex_arp23_split(
    n_total: int,
    arp23_mass_fraction: float,
    formin_length_um: float,
    arp23_length_um: float,
) -> tuple[int, int]:
    r"""Split a fixed cortex filament budget into (n_formin, n_arp23) from an Arp2/3 MASS fraction.

    The physiological datum is a MASS fraction (Bovellan 2014: ~1/3 of cortical F-actin is Arp2/3-nucleated),
    NOT a count fraction — and mass ∝ n·length, so because Arp2/3 filaments are ~10–20× shorter than the
    representative formin filaments, the Arp2/3 COUNT fraction is far larger than the mass fraction. Solving

        f = N_a·L_a / (N_a·L_a + N_f·L_f)   (mass fraction)   and   N_a + N_f = N_total   (fixed budget),

    gives ``M = N_total / (f/L_a + (1−f)/L_f)`` (total contour "mass" in length units), then
    ``N_a = round(f·M/L_a)`` and ``N_f = N_total − N_a`` (so the total is EXACT ⇒ areal density preserved).

    Worked example (the task's illustrative numbers): L_a=0.1 µm, L_f=1.0 µm, f=1/3 ⇒ N_a/N_f = 5 (≈5× more
    Arp2/3 filaments BY COUNT than formin, even though Arp2/3 is only 1/3 by mass).

    Args:
        n_total: fixed cortex filament budget (e.g. 70,686) — preserved exactly across the split.
        arp23_mass_fraction: Arp2/3 fraction of cortical actin BY MASS, in (0, 1) (~0.33; Bovellan 2014).
        formin_length_um: representative formin (long) filament contour length [µm].
        arp23_length_um: representative Arp2/3 (short) filament contour length [µm].

    Returns:
        ``(n_formin, n_arp23)`` with ``n_formin + n_arp23 == n_total`` (both ≥ 1).
    """
    if not (0.0 < arp23_mass_fraction < 1.0):
        raise ValueError("arp23_mass_fraction must be in (0, 1)")
    if formin_length_um <= 0.0 or arp23_length_um <= 0.0:
        raise ValueError("filament lengths must be > 0")
    if n_total < 2:
        raise ValueError("n_total must be >= 2 to split")
    f = float(arp23_mass_fraction)
    total_mass = n_total / (f / arp23_length_um + (1.0 - f) / formin_length_um)
    n_arp23 = int(round(f * total_mass / arp23_length_um))
    n_arp23 = max(1, min(n_total - 1, n_arp23))
    n_formin = n_total - n_arp23
    return n_formin, n_arp23


def _cortex_arp23_arch(
    n_mothers: int, *, length_um: float, seg_um: float, R_um: float, name: str = "cortex_arp23",
) -> ArchitectureSpec:
    """An Arp2/3 cortical sub-population architecture row (short branched filaments on the cortex sphere).

    Distinct manifold ``"sphere_dendritic"`` so ``build_region`` dispatches to the sphere-native dendritic
    builder (NOT the formin sphere delegation). Branch geometry mirrors the lamellipodium row (θ₀=70°,
    σ_θ=9°; Fäßler 2020) — the branch path is REUSED, only the placement manifold differs."""
    from aleph.laws.architecture_spec import (
        ARP23_BRANCH_ANGLE_RAD, ARP23_BRANCH_SIGMA_DEG, CrosslinkerSpec, FilamentSpec, MotorSpec,
    )
    return ArchitectureSpec(
        name=name, manifold="sphere_dendritic", R_um=R_um,
        filament=FilamentSpec(nucleator="arp23", length_um=length_um, seg_um=seg_um, n_filaments=n_mothers,
                              orientation="branched_twomode", polarity="uniform",
                              branch_angle_rad=ARP23_BRANCH_ANGLE_RAD, branch_sigma_deg=ARP23_BRANCH_SIGMA_DEG),
        crosslinker=CrosslinkerSpec(hand=FILAMIN, bind_mode="any", density_per_fil=0.0),
        motor=MotorSpec(hand=None, mode="none"))


def cortex_arp23_region(
    n_arp23: int, *, length_um: float, seg_um: float, mother_fraction: float, R_um: float,
    name: str = "cortex_arp23",
) -> RegionSpec:
    """A cortex Arp2/3 sub-population :class:`RegionSpec` of ``n_arp23`` short branched filaments.

    ``mother_fraction`` of the population are branch ROOTS (mothers); the rest are DAUGHTERS, each attached at
    an Arp2/3 branch junction (so each daughter carries one branch). ``mother_fraction`` and the Arp2/3
    ``length_um``/``seg_um`` are FLAGGED modeling choices (PI A/B/C) — see the memo. Encoded on the standard
    RegionSpec: mothers → ``arch.filament.n_filaments``; daughters → ``n_daughter_pool`` (mirrors the
    lamellipodium ``n_mothers`` + ``n_daughter_pool`` convention)."""
    if n_arp23 < 2:
        raise ValueError("n_arp23 must be >= 2 (need at least one mother + one branched daughter)")
    if not (0.0 < mother_fraction < 1.0):
        raise ValueError("mother_fraction must be in (0, 1)")
    n_mothers = max(1, min(n_arp23 - 1, int(round(mother_fraction * n_arp23))))
    n_daughters = n_arp23 - n_mothers
    arch = _cortex_arp23_arch(n_mothers, length_um=length_um, seg_um=seg_um, R_um=R_um, name=name)
    return RegionSpec(arch=arch, manifold_kind="sphere_dendritic", scale_um=R_um, n_daughter_pool=n_daughters)


def _sphere_tangent_basis(normal: npt.NDArray[np.float64]) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]:
    """Two orthonormal tangent-plane vectors at a unit sphere ``normal`` (for placing filaments on the shell)."""
    ref = np.array([0.0, 0.0, 1.0]) if abs(normal[2]) < 0.9 else np.array([1.0, 0.0, 0.0])
    e1 = ref - np.dot(ref, normal) * normal
    e1 = e1 / (np.linalg.norm(e1) + 1e-12)
    e2 = np.cross(normal, e1)
    return e1, e2


def _relax_arp23_overlaps(
    pos: npt.NDArray[np.float64],
    fiber_of: npt.NDArray[np.int64],
    anchor_pairs: npt.NDArray[np.int64],
    sigma_ev_um: float = 0.007,
    max_sweeps: int = 60,
    relax: float = 1.0,
    overlap_mode: str = "transverse",
    overlap_span: int = 2,
) -> tuple[npt.NDArray[np.float64], dict]:
    r"""Host WCA overlap-relax of the Arp2/3 leaf, holding each daughter-base↔branch-vertex anchor welded.

    Mirrors :func:`ff.cortex_assembly._resolve_cross_fiber_overlaps` (same soft-sphere push to the WCA contact
    distance ``r_c = 2^(1/6)·σ_EV``, a param-free geometric target NOT tuned to a force gate), with two Arp2/3-
    specific constraints so the branch topology survives the relax:

      1. the intentional daughter-base↔branch-vertex coincidences (``anchor_pairs``: co-located nodes on DIFFERENT
         fibers — the rest-0 α-actinin anchor) are EXCLUDED from the push; they are welds, not EV contacts, so
         pushing them apart would fight the weld and never converge; and
      2. after every sweep each daughter base ``d0`` is re-welded onto its (sub-σ-nudged) mother branch vertex
         ``bn``, so the branch stays force-free and the θ₀=70° geometry is preserved (all node moves are ≪ the
         0.5 µm mesh, so the branch angle is unchanged to <1°).

    Each sweep moves every remaining different-fiber sub-``r_c`` node pair symmetrically apart by half their
    overlap. Converges in a handful of sweeps. Returns ``(relaxed_pos, stats)``.
    """
    from scipy.spatial import cKDTree

    if overlap_mode == "radial_span":
        # SMOOTHNESS-PRESERVING out-of-plane resolution (shared core; centre = origin for the cortex sphere).
        # Excludes + re-welds the branch anchors exactly like the transverse path so the θ₀=70° geometry survives.
        from aleph.laws.cortex_assembly import _resolve_cross_fiber_overlaps_radial_span
        F = int(fiber_of.max()) + 1 if fiber_of.size else 0
        lens = np.bincount(fiber_of, minlength=F) if F else np.zeros(0, np.int64)
        offsets = np.concatenate([[0], np.cumsum(lens)]).astype(np.int64)
        p, sweeps, move = _resolve_cross_fiber_overlaps_radial_span(
            pos, offsets, np.zeros(3), sigma_ev_um, max_sweeps, span=overlap_span,
            anchor_pairs=(anchor_pairs if anchor_pairs.size else None))
        return p, {"arp23_overlap_relax_sweeps": int(sweeps),
                   "arp23_overlap_relax_max_move_um": float(move), "overlap_mode": overlap_mode}
    if overlap_mode != "transverse":
        raise ValueError(f"overlap_mode must be 'transverse' or 'radial_span', got {overlap_mode!r}")

    r_c = 2.0 ** (1.0 / 6.0) * float(sigma_ev_um)
    r_target = r_c * 1.01                                       # settle just past the cutoff, not on it
    p = np.ascontiguousarray(pos, np.float64).copy()
    if anchor_pairs.size:
        d0 = anchor_pairs[:, 0].astype(np.int64)
        bn = anchor_pairs[:, 1].astype(np.int64)
        anchor_key = {(min(int(a), int(b)), max(int(a), int(b))) for a, b in anchor_pairs}
    else:
        d0 = bn = np.zeros(0, np.int64)
        anchor_key = set()
    sweeps = 0
    for sweeps in range(1, int(max_sweeps) + 1):
        pairs = cKDTree(p).query_pairs(r_c, output_type="ndarray")           # (i<j) node pairs within r_c
        if pairs.shape[0]:
            pairs = pairs[fiber_of[pairs[:, 0]] != fiber_of[pairs[:, 1]]]     # cross-filament only
        if pairs.shape[0] and anchor_key:                                    # drop the welded anchor pairs
            keep = np.fromiter(((int(a), int(b)) not in anchor_key for a, b in pairs), bool, pairs.shape[0])
            pairs = pairs[keep]
        if pairs.shape[0] == 0:
            break
        i, j = pairs[:, 0], pairs[:, 1]
        dvec = p[i] - p[j]
        r = np.maximum(np.linalg.norm(dvec, axis=1), 1e-9)
        push = (float(relax) * 0.5 * (r_target - r) / r)[:, None] * dvec
        disp = np.zeros_like(p)
        np.add.at(disp, i, push)
        np.add.at(disp, j, -push)
        p = p + disp
        if d0.size:                                                          # re-weld daughter base ↔ vertex
            p[d0] = p[bn]
    return p, {"arp23_overlap_relax_sweeps": int(sweeps),
               "arp23_overlap_relax_max_move_um": float(np.linalg.norm(p - pos, axis=1).max())}


def _build_cortex_arp23_region(spec: RegionSpec, rng: np.random.Generator,
                               overlap_free: bool = False,
                               overlap_mode: str = "transverse", overlap_span: int = 2) -> RegionLeaf:
    """Short branched Arp2/3 cortex sub-population woven onto the cortex SPHERE (REUSES the Arp2/3 branch path).

    Placement (the only new geometry): mother BASES are isotropic on the shell at ``spec.scale_um``; each mother
    lays ``nb`` beads along a random tangent direction. Branching (REUSED): daughters branch off a random mother
    at ``θ₀ ± σ_θ`` drawn by :func:`branch_angle.sample_branch_angles` (θ₀=70°, k_θ=``ARP23_K_THETA``), rotated
    IN THE LOCAL TANGENT PLANE. The ``branch_triples``/``branch_anchors`` are built exactly as
    :func:`ac.weave.lamellipodium.build_lamellipodium` and feed the SAME device ``branch_angle_kernel``. Crosslinks
    = the stiff daughter-base↔branch anchors (force-free at the branch geometry), like the lamellipodium region.

    ``overlap_free`` (opt-in, mirrors the formin cortex path in :func:`_build_cortex_delegated` →
    ``ff.weave.weave(..., overlap_free=True)``): when True, the dense short Arp2/3 sub-population — placed on the
    shell without a WCA pass — is overlap-relaxed in-build by :func:`_relax_arp23_overlaps` (iterative soft-sphere
    push of different-fiber sub-``r_c`` pairs + a per-sweep re-weld of each daughter base onto its branch vertex),
    so no sub-σ inter-filament contact survives while the θ₀=70° branch geometry and rest-0 anchors are preserved.
    The crosslink rest lengths below are then computed from the RESOLVED positions, so they stay force-free at
    rest. Default OFF leaves the raw seed (the historical behaviour)."""
    fs = spec.arch.filament
    R = spec.scale_um
    seg = fs.seg_um
    nb = int(round(fs.length_um / seg)) + 1
    if nb < 3:
        raise ValueError(
            "Arp2/3 cortex filament needs >=3 nodes to carry a branch; refine cortex_arp23_seg_um "
            f"(got length_um={fs.length_um}, seg_um={seg} → {nb} nodes)")
    n_mothers = fs.n_filaments
    n_daughters = spec.n_daughter_pool
    n_fib = n_mothers + n_daughters
    theta0 = fs.branch_angle_rad or ARP23_THETA0_RAD

    # mothers: isotropic base on the sphere, random tangent-plane direction
    normals = rng.standard_normal((n_mothers, 3))
    normals /= np.linalg.norm(normals, axis=1, keepdims=True) + 1e-12
    fibers: list[npt.NDArray[np.float64]] = []
    mdirs: list[npt.NDArray[np.float64]] = []
    mnorms: list[npt.NDArray[np.float64]] = []
    for m in range(n_mothers):
        nrm = normals[m]
        e1, e2 = _sphere_tangent_basis(nrm)
        a = rng.uniform(0.0, 2.0 * np.pi)
        d = np.cos(a) * e1 + np.sin(a) * e2
        base = R * nrm
        fibers.append(base[None, :] + np.arange(nb)[:, None] * seg * d[None, :])
        mdirs.append(d); mnorms.append(nrm)

    # daughters: branch off a random mother at θ₀ ± σ_θ (REUSED arp23 sampling), rotated in the tangent plane
    thetas = sample_branch_angles(theta0, ARP23_K_THETA, n_daughters, rng) if n_daughters else np.zeros(0)
    parent_of: list[tuple[int, int, int]] = []
    for d_local in range(n_daughters):
        mi = int(rng.integers(n_mothers))
        b = int(rng.integers(1, nb - 1))                        # interior branch bead on the mother
        mdir = mdirs[mi]; nrm = mnorms[mi]
        perp = np.cross(nrm, mdir)                              # in-tangent-plane perpendicular
        perp = perp / (np.linalg.norm(perp) + 1e-12)
        side = 1.0 if rng.random() < 0.5 else -1.0
        th = side * float(thetas[d_local])
        ddir = np.cos(th) * mdir + np.sin(th) * perp
        branch_pos = fibers[mi][b]
        fibers.append(branch_pos[None, :] + np.arange(nb)[:, None] * seg * ddir[None, :])
        parent_of.append((n_mothers + d_local, mi, b))

    offsets = np.zeros(n_fib + 1, np.int64)
    for f in range(n_fib):
        offsets[f + 1] = offsets[f] + fibers[f].shape[0]
    pos = np.concatenate(fibers, axis=0)

    triples, anchors = [], []
    for (di, mi, b) in parent_of:
        m_after = int(offsets[mi]) + min(b + 1, nb - 1)         # mother node just past the branch
        bn = int(offsets[mi]) + b                               # branch vertex (on the mother)
        d0 = int(offsets[di]); d1 = int(offsets[di]) + 1        # daughter base (co-located) + first arm node
        triples.append([m_after, bn, d1]); anchors.append([d0, bn])
    branch_triples = np.array(triples, np.int64) if triples else np.zeros((0, 3), np.int64)
    branch_anchors = np.array(anchors, np.int64) if anchors else np.zeros((0, 2), np.int64)
    branch_active = np.ones(branch_triples.shape[0], bool)

    # overlap-free build (opt-in; mirrors the formin cortex overlap_free path): WCA-relax the short Arp2/3 leaf so
    # the dense sub-population starts WITHOUT the ~sub-σ inter-filament interpenetrations that otherwise saturate
    # the build-time steric force, holding every daughter-base↔branch-vertex anchor welded (θ₀=70° preserved).
    if overlap_free and pos.shape[0] > 1:
        fiber_of = np.repeat(np.arange(n_fib, dtype=np.int64), np.diff(offsets))
        pos, _ = _relax_arp23_overlaps(pos, fiber_of, branch_anchors,
                                       overlap_mode=overlap_mode, overlap_span=overlap_span)

    if branch_anchors.size:
        xl_i = branch_anchors[:, 0].astype(np.int64); xl_j = branch_anchors[:, 1].astype(np.int64)
        xl_k = np.full(branch_anchors.shape[0], ALPHA_ACTININ.link_k)
        xl_rest = np.linalg.norm(pos[xl_j] - pos[xl_i], axis=1)
    else:
        xl_i = xl_j = np.zeros(0, np.int64); xl_k = xl_rest = np.zeros(0)

    polarity = np.ones(n_fib, np.int64)                        # barbed distal (last node) for all
    return RegionLeaf(
        name=spec.name, manifold="sphere_dendritic", nucleator="arp23", pos=pos, fiber_offsets=offsets,
        xl_i=xl_i, xl_j=xl_j, xl_k=xl_k, xl_rest=xl_rest, myo_i=np.zeros(0, np.int64), myo_j=np.zeros(0, np.int64),
        polarity=polarity, provenance="sphere_dendritic_arp23_seed", branch_triples=branch_triples,
        branch_anchors=branch_anchors, branch_active=branch_active)


def _build_demoted_bundle(spec: RegionSpec, rng: np.random.Generator) -> RegionLeaf:
    """The DEMOTED scripted ``ff.weave._build_bundle`` — a NON-AUTHORITATIVE control only (emergent=False)."""
    from aleph.laws.weave import _build_bundle
    net = _build_bundle(spec.arch, rng)
    n_fib = net.fiber_offsets.shape[0] - 1
    return RegionLeaf(
        name=spec.name, manifold="bundle", nucleator=spec.arch.filament.nucleator, pos=net.pos.copy(),
        fiber_offsets=net.fiber_offsets.astype(np.int64), xl_i=np.zeros(0, np.int64), xl_j=np.zeros(0, np.int64),
        xl_k=np.zeros(0), xl_rest=np.zeros(0), myo_i=np.zeros(0, np.int64), myo_j=np.zeros(0, np.int64),
        polarity=np.ones(n_fib, np.int64), provenance="demoted_scripted_bundle")


def build_region(spec: RegionSpec, rng: np.random.Generator, overlap_free: bool = False,
                 overlap_mode: str = "transverse", overlap_span: int = 2) -> RegionLeaf:
    """Build one region into a :class:`RegionLeaf` (dispatch by manifold; see module docstring)."""
    if spec.arch.manifold == "sphere":
        return _build_cortex_delegated(spec, rng, overlap_free=overlap_free,
                                       overlap_mode=overlap_mode, overlap_span=overlap_span)
    if spec.arch.manifold == "sphere_dendritic":              # mixed-cortex Arp2/3 sub-population (FLAGGED)
        return _build_cortex_arp23_region(spec, rng, overlap_free=overlap_free,
                                          overlap_mode=overlap_mode, overlap_span=overlap_span)
    if spec.arch.manifold == "patch":
        return _build_lamellipodium_region(spec, rng)
    if spec.arch.manifold == "bundle":
        return _build_emergent_bundle(spec, rng) if spec.emergent else _build_demoted_bundle(spec, rng)
    raise ValueError(f"unsupported manifold {spec.arch.manifold!r} (sphere|patch|bundle)")
