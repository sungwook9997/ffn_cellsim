r"""Lamellipodium REBUILD — dendritic Arp2/3 network (host NumPy builder) — I4 / §3A.b headline net-new.

The old "lamellipodium" was ``ff.motility_warp.directed_front_growth_kernel`` — directed barbed-end growth on
CORTEX fibers, a ratchet BODY-FORCE proxy with NO dendritic topology. This is the PI-mandated rebuild
("갈아엎기"): a genuine branched network —

  * TOPOLOGY: mother filament + Arp2/3 side-branch -> daughter filament, a mother/daughter TREE seeded at the
    leading edge (a flat membrane-adjacent patch). Every junction is an angle-harmonic 3-body branch
    (``ac.weave.branch_angle``), NOT a rigid 72 deg.
  * SEED (emergent, not scripted): mother orientations are ISOTROPIC conditional on the declared leading-edge
    manifold + NPF polarity field (a von-Mises bias toward the protrusion axis whose concentration is the NPF
    field — ablatable; NPF uniform => isotropic). The classic +-35 deg two-mode is NOT hard-set — it is the
    EMERGENT consequence of theta0 = 70 deg branching once mothers align natively (I5 proof). Daughter branch
    angles are THERMAL draws theta0 +- sigma_theta (``branch_angle.sample_branch_angles``), so the measured
    branch-angle distribution matches theta0 +- sigma_theta by construction (Gate).
  * N-FIXED (§5-confirm-1 / growing-N reconciliation): a pre-allocated DORMANT daughter pool is ACTIVATED by
    the flux-limited nucleation KMC (``ac.weave.nucleation.activate_dormant_daughters``) — the dendritic
    topology emerges WITHOUT allocating new nodes. Capping terminates growth (``nucleation.capping_prob``).
  * FORCE: barbed ends push the membrane (Mogilner-Oster ratchet, reused from ``ff.motility_warp``);
    retrograde flow is the emergent reaction. (Force wiring is the lead's; here we build the network + gate.)

Requires the I1c monomer field LIVE (the nucleation flux limit consumes it) — consumed READ-ONLY as a scalar/
array concentration; this module never runs the Warp field.

engine units: length um, angle rad.

Sanity Gate (self-tested in tests/ac/weave/test_lamellipodium_oracle.py):
  * dendritic topology: every daughter has exactly one parent branch (a tree, not a ring); branch_triples and
    branch_anchors reference valid mother/daughter nodes; N stays fixed under activation (no new nodes).
  * branch-angle distribution: the MEASURED junction angles (label-blind ``branch_angle``) match theta0 +-
    sigma_theta (mean within the sin-measure shift; SD within tolerance of sqrt(kT/k_theta)).
  * flux limit: activating with a depleted monomer field leaves daughters dormant (fewer active);
    a rich field activates more -> monotone in monomer concentration.
  * seed is not pre-aligned: with a UNIFORM (isotropic) NPF field the mother nematic order S is low (the
    two-mode is not scripted); the branch angles are still theta0 (topology is dendritic regardless).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import numpy.typing as npt

from aleph.components.weave.branch_angle import ARP23_K_THETA, ARP23_THETA0_RAD, sample_branch_angles, thermal_sigma
from aleph.components.weave.membrane_ratchet import ensemble_protrusion_velocity
from aleph.components.weave.nucleation import activate_dormant_daughters, branch_nucleation_rate, capping_prob

__all__ = [
    "LamellipodiumSeed",
    "build_lamellipodium",
    "LamellipodiumTick",
    "AdaptiveLamellipodium",
]


@dataclass(slots=True)
class LamellipodiumSeed:
    """The dendritic Arp2/3 network build product (local node indices)."""

    pos: npt.NDArray[np.float64]                # (n, 3)
    fiber_offsets: npt.NDArray[np.int64]        # (F+1,)
    branch_triples: npt.NDArray[np.int64]       # (Nbr, 3) [mother_after, branch_node, daughter_first]
    branch_anchors: npt.NDArray[np.int64]       # (Nbr, 2) [daughter_base, branch_node]
    branch_active: npt.NDArray[np.bool_]        # (Nbr,) active branch junctions
    active_mask: npt.NDArray[np.bool_]          # (F,) mothers active; dormant daughters False until nucleated
    polarity: npt.NDArray[np.int64]             # (F,) barbed end at last node (+1) for all (distal = barbed)
    mother_axes: npt.NDArray[np.float64]        # (n_mothers, 3) unit mother directions (for the seed-order gate)


def _dir_in_plane(phi: float) -> npt.NDArray[np.float64]:
    """In-plane unit direction at angle ``phi`` from +y (protrusion axis), patch in x-y."""
    return np.array([np.sin(phi), np.cos(phi), 0.0])


def _rotate_in_plane(v: npt.NDArray[np.float64], angle: float) -> npt.NDArray[np.float64]:
    """Rotate an in-plane (x-y) unit vector by ``angle`` about +z."""
    c, s = np.cos(angle), np.sin(angle)
    return np.array([c * v[0] - s * v[1], s * v[0] + c * v[1], 0.0])


def build_lamellipodium(
    *,
    n_mothers: int,
    n_daughter_pool: int,
    length_um: float,
    seg_um: float,
    patch_half_um: float,
    rng: np.random.Generator,
    theta0: float = ARP23_THETA0_RAD,
    k_theta: float = ARP23_K_THETA,
    npf_kappa: float = 0.0,
    monomer_c: float = 10.0,
    npf_activity: float = 1.0,
    k_arp0: float = 5.0,
    k_m: float = 5.0,
    tau: float = 1.0,
) -> LamellipodiumSeed:
    """Build the dendritic Arp2/3 lamellipodium seed with a pre-allocated dormant daughter pool.

    Args:
        n_mothers: seeded mother filaments (active).
        n_daughter_pool: pre-allocated daughter filaments (dormant until nucleated; N-FIXED capacity).
        length_um: filament contour length [um].
        seg_um: segment rest length [um].
        patch_half_um: leading-edge patch half-extent [um] (mother base placement).
        rng: NumPy generator.
        theta0: rest branch angle [rad] (default Arp2/3 70 deg).
        k_theta: angle stiffness [pN*um/rad^2] (default the equipartition anchor for sigma = 9 deg).
        npf_kappa: von-Mises concentration of the mother-orientation bias toward +y (0 => isotropic;
            larger => more aligned to the protrusion axis). This is the DECLARED NPF field (SEEDED, ablatable).
        monomer_c: local free-monomer concentration from the I1c field [uM] (flux limit; read-only).
        npf_activity: NPF activity in [0,1] gating nucleation.
        k_arp0, k_m: nucleation rate constants [1/s], [uM] (I0-B4 GAP).
        tau: nucleation KMC tick [s].

    Returns:
        A :class:`LamellipodiumSeed`.
    """
    if n_mothers < 1 or n_daughter_pool < 0:
        raise ValueError("n_mothers >= 1 and n_daughter_pool >= 0 required")
    nb = max(2, int(round(length_um / seg_um)) + 1)
    n_fib = n_mothers + n_daughter_pool
    sigma = thermal_sigma(k_theta)

    fibers: list[npt.NDArray[np.float64]] = []
    mother_axes: list[npt.NDArray[np.float64]] = []

    # ---- mothers: isotropic-conditional on the leading-edge manifold + NPF field --------------------
    for _ in range(n_mothers):
        # von-Mises about +y (phi=0); kappa=0 => uniform (isotropic seed, two-mode NOT scripted)
        phi = rng.vonmises(0.0, npf_kappa)
        d = _dir_in_plane(phi)
        base = np.array([rng.uniform(-patch_half_um, patch_half_um),
                         rng.uniform(-patch_half_um, patch_half_um), 0.0])
        fibers.append(base[None, :] + np.arange(nb)[:, None] * seg_um * d[None, :])
        mother_axes.append(d)

    # ---- daughter pool: pre-allocated, branch off a random mother at theta0 +- sigma_theta ----------
    branch_triples, branch_anchors, parent_of = [], [], []
    thetas = sample_branch_angles(theta0, k_theta, n_daughter_pool, rng) if n_daughter_pool else np.zeros(0)
    for d_local in range(n_daughter_pool):
        mi = int(rng.integers(n_mothers))
        b = int(rng.integers(1, nb - 1))                       # interior branch bead on the mother
        mdir = mother_axes[mi]
        side = 1.0 if rng.random() < 0.5 else -1.0
        ddir = _rotate_in_plane(mdir, side * float(thetas[d_local]))
        branch_pos = fibers[mi][b]
        fibers.append(branch_pos[None, :] + np.arange(nb)[:, None] * seg_um * ddir[None, :])
        parent_of.append((n_mothers + d_local, mi, b))

    fibers = np.array(fibers) if False else fibers  # keep list (variable-safe); build offsets below
    offsets = np.zeros(n_fib + 1, np.int64)
    for f in range(n_fib):
        offsets[f + 1] = offsets[f] + fibers[f].shape[0]
    pos = np.concatenate(fibers, axis=0)

    for (di, mi, b) in parent_of:
        m_after = int(offsets[mi]) + min(b + 1, nb - 1)        # mother node just past the branch (mother arm)
        bn = int(offsets[mi]) + b                              # branch vertex (on the mother)
        d0 = int(offsets[di])                                  # daughter base (co-located with the branch)
        d1 = int(offsets[di]) + 1                              # daughter first segment node (daughter arm)
        branch_triples.append([m_after, bn, d1])
        branch_anchors.append([d0, bn])

    branch_triples = np.array(branch_triples, np.int64) if branch_triples else np.zeros((0, 3), np.int64)
    branch_anchors = np.array(branch_anchors, np.int64) if branch_anchors else np.zeros((0, 2), np.int64)

    # ---- N-FIXED activation: flux-limited nucleation activates dormant daughters (no new nodes) -----
    active_mask = np.zeros(n_fib, bool)
    active_mask[:n_mothers] = True                             # mothers active
    rate = np.zeros(n_fib)
    rate[n_mothers:] = branch_nucleation_rate(monomer_c, npf_activity, k_arp0, k_m)
    active_mask = activate_dormant_daughters(active_mask, rate, tau, rng)
    # a branch junction is live iff its daughter is active
    daughter_of_branch = branch_triples[:, 2] if branch_triples.size else np.zeros(0, np.int64)
    branch_active = np.array(
        [bool(active_mask[np.searchsorted(offsets, dn, side="right") - 1]) for dn in daughter_of_branch],
        bool) if branch_triples.size else np.zeros(0, bool)

    polarity = np.ones(n_fib, np.int64)                        # barbed end = distal (last node) for all
    return LamellipodiumSeed(
        pos=pos, fiber_offsets=offsets, branch_triples=branch_triples, branch_anchors=branch_anchors,
        branch_active=branch_active, active_mask=active_mask, polarity=polarity,
        mother_axes=np.array(mother_axes))


# =====================================================================================================
# Adaptive Arp2/3 dynamics — nucleation / capping / Brownian-ratchet protrusion over time (host oracle)
# =====================================================================================================
#
# The :func:`build_lamellipodium` product is a static SEED. The dendritic net is genuinely ADAPTIVE: the
# active/growing filament population evolves under three competing mechanistic channels the roadmap §3.4
# names ("nucleation -> branch angle-harmonic + thermal, capping, retrograde"):
#
#   * NUCLEATION activates dormant daughters at the flux-limited Arp2/3 rate (``branch_nucleation_rate``,
#     Michaelis in the I1c monomer field, gated by NPF) — N-FIXED (pops from the dormant pool, no new nodes).
#   * CAPPING terminates a barbed end's growth at the first-order rate ``k_cap`` (``capping_prob``); a capped
#     filament stays in the network but no longer elongates or pushes the membrane.
#   * DISSOLUTION recycles a capped filament back to the dormant pool at the first-order rate ``k_dissolve``
#     (pointed-end depolymerisation / debranching), CLOSING the treadmilling cycle dormant -> growing ->
#     capped -> dormant. Without this recycle a finite dormant pool exhausts and the array would cap out; the
#     closed cycle is what gives a genuine population STEADY STATE (nucleation flux = capping flux = dissolution
#     flux), independent of the pool size once it exceeds the steady occupancy (§3 free-list discipline).
#   * PROTRUSION advances the growing barbed ends against the membrane load by the elastic Brownian ratchet
#     (``ac.weave.membrane_ratchet``): the load is SHARED across the growing ends, so recruiting branches drops
#     the per-filament load and holds velocity up (the Mogilner-Oster adaptive force-velocity). RETROGRADE flow
#     is the complementary reaction: barbed-end polymer that is NOT converted to protrusion slides rearward.
#     The protrusion/retrograde partition is the motor-clutch engagement ``clutch_engagement`` in [0, 1] — a
#     SEEDED, ablatable input (it EMERGES from the FA clutch at native, the lead's I6 wiring); not tuned here.
#
# The steady state (growing count ~ nucleation-flux / k_cap, mean growing-tip length ~ V/k_cap) is the
# distribution the ensemble harness (``ac.engine.ensemble``) validates over N realizations — never a single
# run. All magnitudes (k_arp0, k_m, k_cap, k_on_c, k_off) are I0-B4 GAP inputs, swept as oracle variables;
# the SHAPE gates (N-fixed bound, monotone load response, nucleation/capping fixed point) are
# magnitude-independent (report-not-tune, §6.2).


@dataclass(slots=True)
class LamellipodiumTick:
    """One physical-tick snapshot of the adaptive lamellipodium (all counts label-blind)."""

    t: float                                      # elapsed time [s]
    n_active: int                                 # active filaments (mothers + nucleated daughters)
    n_growing: int                                # active AND uncapped barbed ends (load-sharing)
    mean_growing_length_um: float                 # mean contour length of growing filaments [um]
    protrusion_velocity_um_s: float               # membrane advance speed this tick [um/s]
    retrograde_velocity_um_s: float               # rearward network flow this tick [um/s]
    membrane_front_um: float                      # accumulated leading-edge position [um]


@dataclass(slots=True)
class AdaptiveLamellipodium:
    """Mutable adaptive state over a fixed :class:`LamellipodiumSeed` topology (N-FIXED).

    Args:
        seed: the built dendritic topology (immutable; supplies mothers + dormant daughter pool).
        seg_um: segment rest length [um] (a nucleated filament starts one segment long).
        n_mothers: number of mother filaments (always active, never capped-off from the count of live net).
        active: (F,) current active flags (evolves as daughters nucleate).
        capped: (F,) barbed-end capped flags (a capped end no longer grows or pushes).
        length_um: (F,) current contour length [um]; dormant filaments are 0 until nucleated.
        membrane_front_um: accumulated leading-edge position [um].
        time_s: elapsed physical time [s].
    """

    seed: LamellipodiumSeed
    seg_um: float
    n_mothers: int
    active: npt.NDArray[np.bool_]
    capped: npt.NDArray[np.bool_]
    length_um: npt.NDArray[np.float64]
    membrane_front_um: float = 0.0
    time_s: float = 0.0

    @classmethod
    def from_seed(cls, seed: LamellipodiumSeed, *, seg_um: float, length_um: float, n_mothers: int) -> "AdaptiveLamellipodium":
        """Initialise the adaptive state: active filaments start at ``length_um``, dormant at 0.

        ``n_mothers`` is retained as topological metadata (which filaments carry the seed's branch triples);
        the dynamics treat every filament uniformly — a lamellipodial array has no permanently privileged
        barbed end, all cycle through dormant -> growing -> capped -> dormant by turnover.
        """
        n_fib = seed.fiber_offsets.shape[0] - 1
        active = np.asarray(seed.active_mask, bool).copy()
        length = np.where(active, float(length_um), 0.0)
        return cls(
            seed=seed, seg_um=float(seg_um), n_mothers=int(n_mothers),
            active=active, capped=np.zeros(n_fib, bool), length_um=length)

    @property
    def growing(self) -> npt.NDArray[np.bool_]:
        """(F,) filaments that are active AND uncapped (the load-sharing, elongating barbed ends)."""
        return self.active & ~self.capped

    def step(
        self,
        *,
        dt: float,
        membrane_load_pN: float,
        k_arp0: float,
        k_m: float,
        k_cap: float,
        k_dissolve: float,
        k_on_c: float,
        k_off: float,
        monomer_c: float,
        npf_activity: float,
        clutch_engagement: float,
        rng: np.random.Generator,
    ) -> LamellipodiumTick:
        """Advance one physical tick: nucleate, dissolve capped, cap growing, protrude by the ratchet.

        Args:
            dt: physical tick [s].
            membrane_load_pN: total opposing membrane load on the growing front [pN] (>= 0).
            k_arp0, k_m: Arp2/3 nucleation kinetics [1/s], [uM] (I0-B4 GAP).
            k_cap: capping-protein termination rate [1/s] (I0-B4 GAP).
            k_dissolve: capped-filament dissolution/recycle rate [1/s] (I0-B4 GAP; closes the turnover cycle).
            k_on_c, k_off: barbed-end elongation/dissociation kinetics [1/s] (I0-B4 GAP).
            monomer_c: local free-monomer concentration from the I1c field [uM] (read-only flux limit).
            npf_activity: NPF activity in [0, 1] gating nucleation (SEEDED membrane field).
            clutch_engagement: protrusion fraction in [0, 1] of barbed-end polymer (SEEDED; emerges from the
                FA clutch at native). ``1`` => all growth becomes protrusion; ``0`` => all becomes retrograde.
            rng: NumPy generator (device RNG mirror).

        Returns:
            A :class:`LamellipodiumTick` snapshot of the post-step state.
        """
        if dt <= 0.0:
            raise ValueError("dt must be > 0")
        if not 0.0 <= clutch_engagement <= 1.0:
            raise ValueError("clutch_engagement must be in [0, 1]")
        n_fib = self.active.shape[0]

        # 1. NUCLEATION — flux-limited activation of dormant filaments (N-FIXED; pops the dormant pool).
        rate = np.zeros(n_fib)
        dormant = ~self.active
        rate[dormant] = branch_nucleation_rate(monomer_c, npf_activity, k_arp0, k_m)
        new_active = activate_dormant_daughters(self.active, rate, dt, rng)
        just_born = new_active & ~self.active
        self.active = new_active
        self.length_um[just_born] = self.seg_um            # a nucleated barbed end starts one segment long

        # 2. DISSOLUTION — capped filaments recycle to the dormant pool (closes the turnover cycle).
        p_dis = capping_prob(dt, k_dissolve)               # 1 - exp(-dt*k_dissolve), same Poisson form
        dissolve_fire = self.capped & (rng.random(n_fib) < p_dis)
        self.active[dissolve_fire] = False
        self.capped[dissolve_fire] = False
        self.length_um[dissolve_fire] = 0.0

        # 3. CAPPING — a first-order KMC channel terminates growth on the currently-growing ends.
        p_cap = capping_prob(dt, k_cap)
        growing = self.growing
        cap_fire = growing & (rng.random(n_fib) < p_cap)
        self.capped[cap_fire] = True

        # 4. PROTRUSION — the growing ends share the membrane load; the ratchet sets the common speed.
        growing = self.growing
        n_growing = int(growing.sum())
        v_poly = 0.0
        if n_growing > 0:
            v_poly = ensemble_protrusion_velocity(membrane_load_pN, n_growing, k_on_c, k_off)
            self.length_um[growing] = np.maximum(self.length_um[growing] + v_poly * dt, 0.0)
        v_protr = clutch_engagement * v_poly
        v_retro = (1.0 - clutch_engagement) * v_poly
        self.membrane_front_um += max(v_protr, 0.0) * dt   # the membrane does not recede on super-stall
        self.time_s += dt

        mean_len = float(self.length_um[growing].mean()) if n_growing > 0 else 0.0
        return LamellipodiumTick(
            t=self.time_s, n_active=int(self.active.sum()), n_growing=n_growing,
            mean_growing_length_um=mean_len, protrusion_velocity_um_s=float(v_protr),
            retrograde_velocity_um_s=float(v_retro), membrane_front_um=float(self.membrane_front_um))

    def run(self, n_steps: int, **step_kwargs: object) -> list[LamellipodiumTick]:
        """Run ``n_steps`` ticks with fixed per-step kwargs; return the per-tick trajectory."""
        if n_steps < 1:
            raise ValueError("n_steps must be >= 1")
        return [self.step(**step_kwargs) for _ in range(n_steps)]  # type: ignore[arg-type]
