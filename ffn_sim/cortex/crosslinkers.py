"""D2 dynamic xlink_head↔actin_cortex bond updater (KU-3.19).

Off-rate per species (mechanism audit 2026-05-30, PI-ratified):
α-actinin = Bell SLIP; filamin = Pereverzev two-pathway CATCH-SLIP
(A1 correction — filamin is a documented catch bond, not a slip bond).

Phase 1 H.3 brief §Crosslinkers (D2). The updater is a
``hoomd.custom.Action`` wrapped in a periodic ``CustomUpdater`` that runs
every ``xlink_batch_steps`` integration steps (D2 batched-updater
convention; gated by the same batch CFL as H.4 IntegrinBondUpdater).

Architecture (per brief §Crosslinkers)
--------------------------------------
Each crosslinker is a TWO-PARTICLE head pair (`xlink_head`) joined by a
permanent intra-xlink harmonic bond (`xlink_intra`):

* α-actinin (KU-3.19 30 %): rest length 35 nm
* filamin (KU-3.19 70 %): rest length 150 nm
* Same intra-xlink stiffness k_xl = 0.1 pN/μm (Furuike 2001 / Ferrer 2008)

Each head can be **bound** to one nearby actin_cortex bead via a dynamic
harmonic bond (`xlink_attach`), OR **unbound** (free to diffuse with the
other head, tethered only by the intra-xlink bond).

Per batch tick the action:

1. Reads the current particle positions via ``sim.state.get_snapshot()``.
2. For every **currently engaged** ``xlink_head ↔ actin_cortex`` bond:
   computes the bond force ``F = k_attach · max(0, |r_head − r_actin| −
   r0)``, samples the species-specific Bell-Evans off-rate
   ``k_off(F) = k_off⁰ · exp(F · x_β / kT)`` (slip bond — monotonically
   increasing in F), draws ``p_break = 1 − exp(−k_off · Δt_batch)``,
   removes the bond on fire.
3. For every **currently unbound** head: queries a freshly-built KDTree
   over all actin_cortex bead positions for any acceptor within
   ``max_bind_dist`` (60 nm per brief). If at least one acceptor exists,
   samples ``p_bind = 1 − exp(−k_on · Δt_batch)`` (same single-head
   acceptance probability regardless of acceptor count, by D2
   contract); on fire, registers the bond to the nearest acceptor.
4. Rebuilds ``bonds.group`` / ``bonds.typeid`` and writes a fresh
   Snapshot.

Bell-Evans batch CFL contract
-----------------------------
``xlink_batch_steps · dt · k_off_max ≤ 10⁻³``, where ``k_off_max`` is
the Bell-Evans rate evaluated at the maximum credible force per the
brief (~1 pN at WCA scale). ``ResolvedCrosslinkers.resolve_xlinks``
enforces this at config-resolve time by shrinking
``xlink_batch_steps`` if violated. The Updater's ``__init__``
re-asserts the bound.

Sanity Gate
-----------
*Per CLAUDE.md Hard Rule. STATIC checks in
``ffn_sim/tests/test_crosslinkers.py``.*

1. **Dimensional analysis**
   - ``Δt_batch = batch_steps · dt`` [s].
   - ``k_off(F) = k_off⁰ · exp(F · x_β / kT)``: ``x_β`` [m], ``F`` [N],
     ``kT`` [J] → ``F · x_β / kT`` dimensionless → ``k_off`` [1/s]. ✓
   - ``k_off · Δt_batch`` dimensionless → break probability. ✓
   - ``k_on · Δt_batch`` dimensionless → binding probability. ✓
   - All positions in metres, forces in Newtons.

2. **Boundary cases**
   - No xlinks at startup (n_xl = 0): generator returns an empty layout;
     the updater idles. RUNTIME pass.
   - Cortex too sparse (no actin within max_bind_dist of any head):
     binding loop finds zero candidates; the updater idles. RUNTIME pass.
   - ``k_off⁰ · Δt_batch ≥ 10⁻³`` (D2 violation): caught by
     ``resolve_xlinks`` → ``xlink_batch_steps`` is shrunk; ``__init__``
     re-asserts.
   - ``x_β ≤ 0`` or ``k_off⁰ ≤ 0``: RUNTIME ValueError.
   - ``max_bind_dist ≤ 0``: RUNTIME ValueError.

3. **Conservation invariants**
   - Particle count after extend is ``n_cortex_beads + 2·n_xl``
     (cortex beads unchanged; xlink_head particles added).
   - Intra-xlink bond count is invariant ``= n_xl`` (each xlink has one
     permanent head-to-head bond).
   - Dynamic head-to-actin bond count is bounded by ``2·n_xl`` (each of
     the 2N heads can attach at most one actin bond at any time).
   - Newton 3rd law: dynamic ``xlink_attach`` bonds use
     ``md.bond.Harmonic`` (symmetric force pair); break/bind events
     mutate topology but never inject net momentum (no thermostat
     coupling on the kinetic step).

4. **Numerical sanity**
   - dt of the host simulation is provided at Updater construction; the
     CFL gate is checked once at ``__init__``.
   - Per-batch act() is O(n_engaged + n_unbound · log(n_actin)) via
     scipy.spatial.cKDTree for the binding-acceptor query. For a single
     cortex (1000 xlinks, 7000 actin beads), per-tick wall-time ≈ tens
     of ms — well under the integration cost of 100 BAOAB steps between
     ticks (~100s of ms).
   - RNG isolated per Action (``np.random.default_rng`` with
     ``simulation.seed + seed_offset``).

5. **Sign / sense** (A1 correction 2026-05-30: filamin is a CATCH bond)
   - α-actinin: Bell SLIP bond — ``k_off(F)`` monotonically INCREASES
     with F (force accelerates unbinding). STATIC test asserts
     ``d k_off / d F > 0`` at all F for α-actinin params.
   - filamin: Pereverzev two-pathway CATCH-SLIP bond — ``k_off(F)``
     FALLS with F (force stabilises) up to a peak force F*, then RISES
     (slip branch). Filamin is a documented catch bond (Ehrlicher 2011
     Nature; Rognoni 2012 PNAS; Gieseke/Rief 2013), same catch-slip
     family as KU-2.5 Pereverzev for H.4 integrin↔ligand. STATIC test
     asserts the catch-slip signature (off-rate falls then rises;
     ``k_off(F*) < k_off(0)``).
   - Bond force magnitude ``F = k_attach · max(0, |Δr| − r0)``; with
     binned ``r0 = bin_center`` and ``|Δr| ≈ bin_center`` at construction,
     F ≈ 0 → ``k_off ≈ k_off⁰`` → binding equilibrium dominated by
     k_on/k_off⁰ ratio.

6. **Measurement-protocol consistency**
   - The equilibrium fraction of engaged xlink-heads, in the absence of
     external force (F ≈ 0), should approach ``k_on/(k_on + k_off⁰)``
     per single-state two-state kinetics. The demo gate runs a short
     equilibration and asserts the empirical fraction is within 50 %
     of this prediction (loose tolerance for the smoke run; production
     gate ``KU-3.19_PRODUCTION=1`` will tighten).
   - The xlink-attach-bond population should equilibrate within a few
     1/k_off⁰ characteristic times.

References
----------
- Brief: ``ffn_sim/docs/briefs/H3_cortex.md`` §Crosslinkers (D2).
- KU-3.19 (Stricker 2010 α-actinin / filamin ratio, Furuike 2001 filamin
  rates, Wachsstock 1994 / Goldmann 2000 α-actinin rates).
- Bell-Evans: Bell 1978 (Science 200:618), Evans-Ritchie 1997.
- D2 batched-updater contract: shared with H.4
  ``ffn_sim/bridge/integrin_bonds.py``.
- ``ffn_sim/integrator/baoab.py`` (D3 BAOAB, frozen 2026-05-20).
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Any

import numpy as np

import hoomd
import hoomd.custom
import hoomd.md as md


# ---------------------------------------------------------------------------
# Resolved parameters
# ---------------------------------------------------------------------------
@dataclass(slots=True)
class ResolvedCrosslinkers:
    """KU-3.19 crosslinker parameters + Bell-Evans species rates."""

    n_xl: int                    # total xlinks per cortex
    alpha_fraction: float        # 0.30 per KU-3.19
    alpha_length: float          # m   α-actinin rest length (35 nm)
    filamin_length: float        # m   filamin rest length (150 nm)
    k_intra: float               # N/m intra-xlink head-to-head k
    k_attach: float              # N/m head-to-actin attach bond k
    max_bind_dist: float         # m   acceptor search radius (60 nm)

    # α-actinin: Bell SLIP (k_off(F) = k_off⁰ · exp(+F·x_β/kT)).
    # C4 correction 2026-05-30: re-anchored to Ferrer 2008 PNAS single-
    # molecule value 0.066 /s (was 1.0 /s, bulk Wachsstock 1994).
    alpha_k_off0: float          # 1/s α-actinin unloaded off-rate (Ferrer 2008 PNAS)
    alpha_x_beta: float          # m   α-actinin Bell strength length (0.4 nm)
    # filamin: Pereverzev two-pathway CATCH-SLIP (A1 correction 2026-05-30).
    # k_off(F) = k_catch0·exp(−F·x_catch/kT) + k_slip0·exp(+F·x_slip/kT).
    # Filamin is a documented catch bond (Ehrlicher 2011; Rognoni 2012;
    # Gieseke/Rief 2013) — force stabilises the bond up to peak F* then
    # destabilises it. Previously mis-modelled as a pure slip bond.
    filamin_k_catch0: float      # 1/s filamin catch-pathway zero-force rate (Furuike 2001)
    filamin_x_catch: float       # m   filamin catch-pathway distance (Ehrlicher 2011 / Rognoni 2012)
    filamin_k_slip0: float       # 1/s filamin slip-pathway zero-force prefactor (Pereverzev 2005)
    filamin_x_slip: float        # m   filamin slip-pathway distance (Furuike 2001 Bell length)

    k_on: float                  # 1/s single-head binding rate when acceptor present

    # D2 batch
    batch_steps: int             # BAOAB steps between Updater ticks
    dt: float                    # host sim dt [s]

    # Per-bin r0 binning for force-free initial attachment (matches H.1 cortex.py)
    n_bins: int
    seed: int

    # Derived
    n_alpha: int = 0
    n_filamin: int = 0
    batch_dt: float = 0.0        # batch_steps · dt
    k_off_max: float = 0.0       # max(k_off(F_max_credible)) — for CFL
    extras: dict[str, Any] = field(default_factory=dict)


def _require_finite_positive(name: str, x: float) -> None:
    if not (math.isfinite(x) and x > 0.0):
        raise ValueError(f"{name} must be finite and > 0, got {x!r}")


def resolve_crosslinkers(cfg: dict, *, dt: float) -> ResolvedCrosslinkers:
    """Resolve crosslinker config block with D2 batch CFL enforcement.

    ``cfg`` is the ``cortex.dynamic_crosslinkers`` sub-dict OR the raw
    cortex config (will look up ``dynamic_crosslinkers`` key). ``dt`` is
    the host simulation timestep (typically ``p_cortex.dt_cfl``).
    """
    if "cortex" in cfg:
        cfg = cfg["cortex"]
    if "dynamic_crosslinkers" in cfg:
        cfg = cfg["dynamic_crosslinkers"]

    p = ResolvedCrosslinkers(
        n_xl=int(cfg["n_xl"]),
        alpha_fraction=float(cfg["alpha_fraction"]),
        alpha_length=float(cfg["alpha_length"]),
        filamin_length=float(cfg["filamin_length"]),
        k_intra=float(cfg["k_intra"]),
        k_attach=float(cfg["k_attach"]),
        max_bind_dist=float(cfg["max_bind_dist"]),
        alpha_k_off0=float(cfg["alpha_k_off0"]),
        alpha_x_beta=float(cfg["alpha_x_beta"]),
        filamin_k_catch0=float(cfg["filamin_k_catch0"]),
        filamin_x_catch=float(cfg["filamin_x_catch"]),
        filamin_k_slip0=float(cfg["filamin_k_slip0"]),
        filamin_x_slip=float(cfg["filamin_x_slip"]),
        k_on=float(cfg["k_on"]),
        batch_steps=int(cfg["batch_steps"]),
        dt=float(dt),
        n_bins=int(cfg.get("n_bins", 10)),
        seed=int(cfg.get("seed", 42)),
    )

    # ---- §2 boundary checks ----
    for name, x in [
        ("alpha_length", p.alpha_length),
        ("filamin_length", p.filamin_length),
        ("k_intra", p.k_intra), ("k_attach", p.k_attach),
        ("max_bind_dist", p.max_bind_dist),
        ("alpha_k_off0", p.alpha_k_off0),
        ("alpha_x_beta", p.alpha_x_beta),
        ("filamin_k_catch0", p.filamin_k_catch0),
        ("filamin_x_catch", p.filamin_x_catch),
        ("filamin_k_slip0", p.filamin_k_slip0),
        ("filamin_x_slip", p.filamin_x_slip),
        ("k_on", p.k_on), ("dt", p.dt),
    ]:
        _require_finite_positive(name, x)
    if not (0.0 <= p.alpha_fraction <= 1.0):
        raise ValueError(f"alpha_fraction must be in [0, 1]; got {p.alpha_fraction}")
    if p.n_xl < 0:
        raise ValueError(f"n_xl must be ≥ 0; got {p.n_xl}")
    if p.batch_steps < 1:
        raise ValueError(f"batch_steps must be ≥ 1; got {p.batch_steps}")
    if p.n_bins < 1:
        raise ValueError(f"n_bins must be ≥ 1; got {p.n_bins}")

    # ---- §1 derived: counts + batch CFL ----
    p.n_alpha = int(round(p.n_xl * p.alpha_fraction))
    p.n_filamin = p.n_xl - p.n_alpha
    p.batch_dt = p.batch_steps * p.dt

    # ---- §2 batch CFL: shrink batch_steps until k_off_max · batch_dt ≤ 1e-3 ----
    # Max credible bond force: roughly k_attach · max_bind_dist ~ 6 pN
    # for the brief defaults (1e-7 · 60e-9 = 6e-15 N — far smaller than 1 pN
    # = 1e-12 N because k_attach is tiny). Use a more meaningful bound:
    # at thermal-scale displacement |Δr| ~ √(kT/k_attach) ≈ 200 nm, force ≈
    # k · 200 nm = 2e-14 N ≈ 0.02 pN. The force-dependent amplification at
    # this F is ≈ 1 for both species (α-actinin slip exp(+F·x_β/kT) and
    # filamin catch-slip). So at smoke scale k_off_max ≈ max(α_k_off0,
    # filamin k_catch0 + k_slip0). With α_k_off0 = 0.066 /s (Ferrer 2008)
    # and filamin ≈ 0.12 /s (catch0 0.1 + slip0 0.02), k_off_max ≈ 0.12 /s.
    # CFL: batch_dt · 0.12 ≤ 1e-3 → batch_dt ≤ 8.3 ms. With dt ≈ 13 ns,
    # batch_steps ≤ 6.4e5 — any reasonable batch_steps (100) is well within.
    # F at fully-stretched dynamic bond at the bind-radius edge (envelope).
    _F_env = p.k_attach * p.max_bind_dist
    _kT = 4.28e-21
    # α-actinin slip envelope: k_off0 · exp(+F·x_β/kT).
    _k_off_alpha = p.alpha_k_off0 * math.exp(p.alpha_x_beta * _F_env / _kT)
    # filamin catch-slip envelope: catch branch decays, slip branch grows;
    # at large F the slip branch dominates → take catch0 (F=0 ceiling of the
    # catch branch) + slip0·exp(+F·x_slip/kT) as a conservative upper bound.
    _k_off_filamin = (
        p.filamin_k_catch0
        + p.filamin_k_slip0 * math.exp(p.filamin_x_slip * _F_env / _kT)
    )
    p.k_off_max = max(_k_off_alpha, _k_off_filamin)
    cfl_product = p.batch_dt * p.k_off_max
    if cfl_product > 1.0e-3:
        # Shrink batch_steps to satisfy the CFL.
        target_batch_dt = 1.0e-3 / p.k_off_max
        new_batch_steps = max(1, int(math.floor(target_batch_dt / p.dt)))
        if new_batch_steps < 1:
            raise RuntimeError(
                f"D2 batch CFL impossible to satisfy: dt={p.dt:.3e} > "
                f"1e-3 / k_off_max ({1e-3 / p.k_off_max:.3e})."
            )
        p.batch_steps = new_batch_steps
        p.batch_dt = p.batch_steps * p.dt
        p.extras["batch_steps_shrunk_from"] = int(cfg["batch_steps"])

    # ---- §4 finite checks on derived ----
    for name in ("batch_dt", "k_off_max"):
        v = getattr(p, name)
        if not (math.isfinite(v) and v > 0.0):
            raise RuntimeError(
                f"Derived {name}={v!r} is not finite-positive."
            )

    return p


# ---------------------------------------------------------------------------
# Per-bin r0 helpers (force-free attachment construction; mirrors H.1 pattern)
# ---------------------------------------------------------------------------
def xlink_attach_bin_names(n_bins: int) -> list[str]:
    """HOOMD bond type names for per-r0 binning of dynamic attach bonds."""
    return [f"xlink_attach_b{i}" for i in range(n_bins)]


def xlink_attach_bin_rest_lengths(n_bins: int, max_bind_dist: float) -> np.ndarray:
    """Bin-center rest lengths uniformly across (0, max_bind_dist].

    First bin center > 0 so that any well-formed dynamic attach bond
    placed within max_bind_dist has a positive ``r0`` to anchor to;
    a length-0 r0 would cause a zero-force WCA-like degeneracy.
    """
    edges = np.linspace(0.0, max_bind_dist, n_bins + 1, dtype=np.float64)
    return 0.5 * (edges[:-1] + edges[1:])


def _bin_index_for_distance(dist: float, n_bins: int, max_bind_dist: float) -> int:
    bin_width = max_bind_dist / n_bins
    idx = int(min(n_bins - 1, max(0, int(dist / bin_width))))
    return idx


# ---------------------------------------------------------------------------
# Topology extension (add xlink_head particles to a cortex frame)
# ---------------------------------------------------------------------------
@dataclass(slots=True)
class XlinkLayout:
    """Per-xlink bookkeeping: tag ranges + species assignment.

    Attributes
    ----------
    head_tag_pairs : ndarray, shape (n_xl, 2)
        For xlink i, the (head_a_tag, head_b_tag) particle tags in the
        post-extension frame.
    species : ndarray of "alpha"/"filamin", shape (n_xl,)
    intra_r0 : ndarray, shape (n_xl,)
        Intra-xlink head-to-head rest length per xlink (set by species).
    head_positions : ndarray, shape (2·n_xl, 3)
        Initial head positions in the post-extension frame.
    actin_tag_seed : int
        Tag offset where the first xlink_head particle was inserted
        (= n_cortex_beads); subsequent head_tags are this + 1, 2, ...
    """

    head_tag_pairs: np.ndarray
    species: np.ndarray
    intra_r0: np.ndarray
    head_positions: np.ndarray
    actin_tag_seed: int


def generate_xlink_layout(
    cortex_positions: np.ndarray,
    cortex_filament_idx: np.ndarray,
    p_xl: ResolvedCrosslinkers,
    *,
    n_cortex_beads: int,
    rng: np.random.Generator | None = None,
) -> XlinkLayout:
    """Place ``n_xl`` xlinks: each as a head-pair near a pair of actin beads
    from different filaments.

    Algorithm (per xlink):

    1. Pick a random actin bead A.
    2. Find all actin beads B from a DIFFERENT filament within
       ``max_bind_dist``.
    3. If at least one candidate exists, sample B uniformly. Place
       head_a at ``A + 1/3 · (B − A)`` and head_b at ``A + 2/3 · (B − A)``
       — both inside the bind-radius shells of A and B respectively, so
       initial attach-bonds (handled by the dynamic Updater on its first
       tick) will find candidates immediately.
    4. If no candidate exists, the xlink is "homeless": place its
       head-pair at A + offset along a random direction at the
       intra-xlink rest length. Both heads will be unbound at start;
       the Updater will bind them when actin diffuses into range.

    Species (α-actinin / filamin) assigned by random with fractions per
    KU-3.19.
    """
    if rng is None:
        rng = np.random.default_rng(p_xl.seed)

    if p_xl.n_xl == 0:
        return XlinkLayout(
            head_tag_pairs=np.empty((0, 2), dtype=np.int64),
            species=np.empty((0,), dtype="<U7"),
            intra_r0=np.empty((0,), dtype=np.float64),
            head_positions=np.empty((0, 3), dtype=np.float64),
            actin_tag_seed=n_cortex_beads,
        )

    n_actin = cortex_positions.shape[0]
    head_positions = np.empty((2 * p_xl.n_xl, 3), dtype=np.float64)
    head_tag_pairs = np.empty((p_xl.n_xl, 2), dtype=np.int64)
    species = np.empty((p_xl.n_xl,), dtype="<U7")
    intra_r0 = np.empty((p_xl.n_xl,), dtype=np.float64)

    # Jitter σ = 50 nm (~3× WCA σ/2) so two coincidentally-placed heads
    # are typically >70 nm apart — outside the WCA cutoff (r_cut ≈ 67 nm
    # for σ_LJ = 60 nm). 5 nm was empirically too small (test_cell_full
    # 단계 6 showed 12 nm separation between collided heads → LJ ~1e-12 J
    # → BAOAB runaway). 50 nm is still well below the bind radius (60 nm
    # for production) so subsequent dynamic attach binding remains
    # geometrically plausible.
    jitter_sigma = 50.0e-9

    n_homeless = 0
    for i in range(p_xl.n_xl):
        # Species assignment
        sp = "alpha" if rng.random() < p_xl.alpha_fraction else "filamin"
        species[i] = sp
        intra_r0[i] = (
            p_xl.alpha_length if sp == "alpha" else p_xl.filamin_length
        )

        # Pick anchor bead A
        a_idx = int(rng.integers(0, n_actin))
        a_fil = cortex_filament_idx[a_idx]
        r_a = cortex_positions[a_idx]
        # Find different-filament candidates within max_bind_dist
        diff = cortex_positions - r_a
        d2 = np.einsum("ij,ij->i", diff, diff)
        mask = (
            (cortex_filament_idx != a_fil)
            & (d2 <= p_xl.max_bind_dist ** 2)
        )
        candidates = np.nonzero(mask)[0]
        if candidates.size > 0:
            b_idx = int(rng.choice(candidates))
            r_b = cortex_positions[b_idx]
            # Place heads at 1/3 and 2/3 between A and B, with jitter
            # (see jitter_sigma rationale above).
            jitter = rng.standard_normal((2, 3)) * jitter_sigma
            head_positions[2 * i] = r_a + (r_b - r_a) / 3.0 + jitter[0]
            head_positions[2 * i + 1] = (
                r_a + 2.0 * (r_b - r_a) / 3.0 + jitter[1]
            )
        else:
            # Homeless: place head_a near r_a (with jitter to avoid
            # stacking when multiple homeless xlinks pick same anchor);
            # head_b at intra-xlink rest length from head_a along a
            # random direction.
            n_homeless += 1
            dir_random = rng.standard_normal(3)
            dir_random /= np.linalg.norm(dir_random)
            jitter_a = rng.standard_normal(3) * jitter_sigma
            head_positions[2 * i] = r_a + jitter_a
            head_positions[2 * i + 1] = (
                r_a + jitter_a + intra_r0[i] * dir_random
            )

        head_tag_pairs[i, 0] = n_cortex_beads + 2 * i
        head_tag_pairs[i, 1] = n_cortex_beads + 2 * i + 1

    p_xl.extras["n_homeless"] = n_homeless
    p_xl.extras["n_alpha_realised"] = int((species == "alpha").sum())
    p_xl.extras["n_filamin_realised"] = int((species == "filamin").sum())

    return XlinkLayout(
        head_tag_pairs=head_tag_pairs,
        species=species,
        intra_r0=intra_r0,
        head_positions=head_positions,
        actin_tag_seed=n_cortex_beads,
    )


# ---------------------------------------------------------------------------
# Dynamic xlink_head ↔ actin_cortex Bell-Evans Updater
# ---------------------------------------------------------------------------
def _bell_evans_k_off(F_mag: np.ndarray, k_off0: float, x_beta: float, kT: float) -> np.ndarray:
    """Bell-Evans SLIP off-rate ``k_off(F) = k_off⁰ · exp(+F · x_β / kT)``.

    Used for α-actinin (the fluid, force-released crosslinker). Off-rate
    rises monotonically with force.
    """
    return k_off0 * np.exp(F_mag * x_beta / kT)


def _pereverzev_catch_slip_k_off(
    F_mag: np.ndarray,
    k_catch0: float,
    x_catch: float,
    k_slip0: float,
    x_slip: float,
    kT: float,
) -> np.ndarray:
    """Pereverzev (2005) two-pathway CATCH-SLIP off-rate.

    ``k_off(F) = k_catch0 · exp(−F · x_catch / kT)
               + k_slip0  · exp(+F · x_slip  / kT)``

    Used for filamin, a documented catch bond (Ehrlicher 2011 Nature;
    Rognoni 2012 PNAS; Gieseke/Rief 2013). The catch branch DECREASES the
    off-rate with force (force stabilises the bond); the slip branch
    INCREASES it. With ``k_catch0 > k_slip0`` and ``x_catch > x_slip`` the
    catch branch dominates at low force, giving the catch-slip signature:
    ``k_off`` falls to a minimum at the peak force

        ``F* = (kT/(x_catch+x_slip))·ln((k_catch0·x_catch)/(k_slip0·x_slip))``

    then rises. This is the same two-pathway functional structure the H.4
    integrin updater uses via ``validation/pereverzev.pereverzev_k_off``;
    re-implemented locally here because runtime modules may not import the
    validation oracles (CLAUDE.md hard rule). The oracle remains the
    acceptance authority in the tests.
    """
    F_mag = np.asarray(F_mag, dtype=np.float64)
    return (
        k_catch0 * np.exp(-F_mag * x_catch / kT)
        + k_slip0 * np.exp(+F_mag * x_slip / kT)
    )


class XlinkBondUpdater(hoomd.custom.Action):
    """D2 Bell-Evans slip dynamic xlink_head↔actin_cortex bond updater.

    Parameters
    ----------
    p : ResolvedCrosslinkers
        Resolved crosslinker config (provides batch_dt, k_attach,
        Bell-Evans rates).
    layout : XlinkLayout
        Per-xlink head tag pairs + species assignment (built by
        ``generate_xlink_layout``).
    kT : float
        Thermal energy for Bell-Evans exponent.
    n_cortex_actin : int
        Number of actin_cortex particles (tags [0, n_cortex_actin) are
        actin; tags [n_cortex_actin, n_cortex_actin + 2·n_xl) are heads).
    seed_offset : int, default 2
        Offset added to ``simulation.seed`` for this Action's RNG.
        Default 2 to avoid collision with BAOAB Updater (offset 0) and
        any other dynamic-bond Updater (e.g. H.4 IntegrinBondUpdater
        uses offset 1).

    Notes
    -----
    The action does NOT mutate particle positions; only the
    ``xlink_attach_*`` bond entries in ``bonds.group`` /
    ``bonds.typeid``. All other bonds (cortex-bond, xlink_intra) are
    preserved across calls.
    """

    def __init__(
        self,
        *,
        p: ResolvedCrosslinkers,
        layout: XlinkLayout,
        kT: float,
        n_cortex_actin: int,
        seed_offset: int = 2,
    ) -> None:
        super().__init__()
        self.p = p
        self.layout = layout
        self.kT = float(kT)
        self.n_cortex_actin = int(n_cortex_actin)
        self._rng = np.random.default_rng(p.seed + seed_offset)

        # Re-assert D2 batch CFL.
        cfl_product = p.batch_dt * p.k_off_max
        if cfl_product > 1.0e-3 + 1.0e-12:
            raise RuntimeError(
                f"XlinkBondUpdater violates D2 batch CFL: "
                f"batch_dt·k_off_max = {p.batch_dt:.3e}·{p.k_off_max:.3e}"
                f" = {cfl_product:.3e} > 1e-3 ceiling. "
                "Shrink batch_steps or pre-resolve via resolve_crosslinkers."
            )

        # Per-head state: bound or unbound, and (if bound) the actin tag.
        n_heads = 2 * p.n_xl
        self._head_bound_to_actin = np.full(n_heads, -1, dtype=np.int64)
        # Per-head species (lookup-table for Bell-Evans rate selection).
        self._head_species = np.empty(n_heads, dtype="<U7")
        for i in range(p.n_xl):
            self._head_species[2 * i] = layout.species[i]
            self._head_species[2 * i + 1] = layout.species[i]

        self._sim_ref: hoomd.Simulation | None = None
        self._steps_run = 0
        self._n_break_total = 0
        self._n_bind_total = 0

    def attach(self, simulation: hoomd.Simulation) -> None:  # noqa: D401
        super().attach(simulation)
        self._sim_ref = simulation

    def act(self, timestep: int) -> None:  # noqa: D401
        sim = self._sim_ref
        assert sim is not None

        read_snap = sim.state.get_snapshot()
        if read_snap.communicator.rank != 0:
            return

        pos = np.asarray(read_snap.particles.position, dtype=np.float64).copy()
        bg = np.asarray(read_snap.bonds.group, dtype=np.int64).copy()
        bt = np.asarray(read_snap.bonds.typeid, dtype=np.uint32).copy()
        bond_type_names = list(read_snap.bonds.types)

        # Indices of xlink_attach_* bin types in the bond type list.
        attach_bin_typeids = [
            i for i, name in enumerate(bond_type_names)
            if name.startswith("xlink_attach_b")
        ]
        if not attach_bin_typeids:
            raise RuntimeError(
                "snap.bonds.types must include xlink_attach_b* types; "
                f"got {bond_type_names}."
            )

        is_attach = np.isin(bt, np.asarray(attach_bin_typeids, dtype=bt.dtype))
        attach_bonds = bg[is_attach]
        other_bg = bg[~is_attach]
        other_bt = bt[~is_attach]

        # Convention: every attach bond has columns [head_tag, actin_tag].
        # head_tag ≥ self.n_cortex_actin (heads come AFTER actin in the tag space).
        if attach_bonds.shape[0] > 0:
            head_tags = attach_bonds[:, 0]
            actin_tags = attach_bonds[:, 1]
            # Validate ordering (defensive — catches misuse early).
            if (head_tags < self.n_cortex_actin).any():
                raise RuntimeError(
                    "xlink_attach bond column 0 must be a head tag "
                    f"(≥ {self.n_cortex_actin}); got tag(s) below this."
                )
            r_head = pos[head_tags]
            r_actin = pos[actin_tags]
            dr = r_head - r_actin
            r = np.linalg.norm(dr, axis=1)
            # Per-bond r0: the bin center this bond was placed in.
            bond_bins = bt[is_attach] - attach_bin_typeids[0]
            bin_r0 = xlink_attach_bin_rest_lengths(
                self.p.n_bins, self.p.max_bind_dist
            )
            r0_per_bond = bin_r0[bond_bins]
            F_mag = self.p.k_attach * np.clip(r - r0_per_bond, 0.0, None)

            # Per-bond species (look up head tag → species → Bell-Evans rates).
            head_local_idx = head_tags - self.n_cortex_actin
            sp_per_bond = self._head_species[head_local_idx]
            alpha_mask = sp_per_bond == "alpha"
            k_off = np.empty(attach_bonds.shape[0], dtype=np.float64)
            # α-actinin: Bell SLIP (k_off rises with F).
            k_off[alpha_mask] = _bell_evans_k_off(
                F_mag[alpha_mask], self.p.alpha_k_off0,
                self.p.alpha_x_beta, self.kT,
            )
            # filamin: Pereverzev two-pathway CATCH-SLIP (A1 correction;
            # k_off falls then rises with F, peak at F*).
            k_off[~alpha_mask] = _pereverzev_catch_slip_k_off(
                F_mag[~alpha_mask],
                self.p.filamin_k_catch0, self.p.filamin_x_catch,
                self.p.filamin_k_slip0, self.p.filamin_x_slip,
                self.kT,
            )
            p_break = 1.0 - np.exp(-k_off * self.p.batch_dt)
            u = self._rng.uniform(0.0, 1.0, size=attach_bonds.shape[0])
            broke = u < p_break
            n_broke = int(broke.sum())
            self._n_break_total += n_broke
            if n_broke > 0:
                broken_heads = head_local_idx[broke]
                self._head_bound_to_actin[broken_heads] = -1
                attach_bonds = attach_bonds[~broke]
                bond_bins = bond_bins[~broke]

        else:
            bond_bins = np.empty((0,), dtype=np.int64)

        # ---- Step 2: bind unbound heads via KDTree query over actin ----
        unbound_local = np.flatnonzero(self._head_bound_to_actin < 0)
        if unbound_local.size > 0:
            unbound_tags = unbound_local + self.n_cortex_actin
            r_heads = pos[unbound_tags]

            # KDTree over actin_cortex positions; query within max_bind_dist.
            try:
                from scipy.spatial import cKDTree
            except ImportError:
                raise RuntimeError(
                    "scipy.spatial.cKDTree required for XlinkBondUpdater "
                    "binding-acceptor query."
                )
            r_actin_all = pos[:self.n_cortex_actin]
            tree = cKDTree(r_actin_all)
            # Returns list[ndarray] of indices per head.
            nbr_lists = tree.query_ball_point(
                r_heads, r=self.p.max_bind_dist
            )

            # Per-head: sample binding probability if at least one candidate.
            p_bind = 1.0 - np.exp(-self.p.k_on * self.p.batch_dt)
            u2 = self._rng.uniform(0.0, 1.0, size=unbound_local.size)
            new_bonds_list = []
            new_bins_list = []
            for k, nbrs in enumerate(nbr_lists):
                if len(nbrs) == 0:
                    continue
                if u2[k] >= p_bind:
                    continue
                # Bind to NEAREST candidate.
                nbrs_arr = np.asarray(nbrs, dtype=np.int64)
                d_nbrs = np.linalg.norm(
                    r_actin_all[nbrs_arr] - r_heads[k], axis=1
                )
                nearest_local = int(np.argmin(d_nbrs))
                actin_tag = int(nbrs_arr[nearest_local])
                d_use = float(d_nbrs[nearest_local])
                bin_idx = _bin_index_for_distance(
                    d_use, self.p.n_bins, self.p.max_bind_dist
                )
                head_tag = int(unbound_tags[k])
                # Validate "different particle" — head and actin tags are
                # by construction distinct (heads are appended AFTER actin),
                # but a defensive check is cheap.
                if head_tag == actin_tag:
                    continue
                new_bonds_list.append((head_tag, actin_tag))
                new_bins_list.append(bin_idx)
                self._head_bound_to_actin[unbound_local[k]] = actin_tag

            if new_bonds_list:
                new_bonds = np.array(new_bonds_list, dtype=np.int64)
                new_bins = np.array(new_bins_list, dtype=np.int64)
                attach_bonds = np.concatenate(
                    [attach_bonds, new_bonds], axis=0
                )
                bond_bins = np.concatenate([bond_bins, new_bins], axis=0)
                self._n_bind_total += int(new_bonds.shape[0])

        # ---- Step 3: rebuild Snapshot with updated attach topology ----
        attach_typeids = (
            np.asarray(bond_bins, dtype=np.uint32) + attach_bin_typeids[0]
        )
        if other_bg.shape[0] > 0:
            new_bg = np.concatenate(
                [attach_bonds.astype(np.uint32), other_bg.astype(np.uint32)],
                axis=0,
            )
            new_bt = np.concatenate(
                [attach_typeids, other_bt.astype(np.uint32)]
            )
        else:
            new_bg = attach_bonds.astype(np.uint32)
            new_bt = attach_typeids

        write_snap = hoomd.Snapshot()
        N_part = int(read_snap.particles.N)
        write_snap.particles.N = N_part
        write_snap.particles.types = list(read_snap.particles.types)
        write_snap.particles.typeid[:] = np.asarray(read_snap.particles.typeid)
        write_snap.particles.position[:] = pos
        write_snap.particles.velocity[:] = np.asarray(read_snap.particles.velocity)
        write_snap.particles.mass[:] = np.asarray(read_snap.particles.mass)
        write_snap.particles.image[:] = np.asarray(read_snap.particles.image)
        box = read_snap.configuration.box
        write_snap.configuration.box = list(box)

        write_snap.bonds.N = int(new_bg.shape[0])
        write_snap.bonds.types = list(bond_type_names)
        if new_bg.shape[0] > 0:
            write_snap.bonds.group[:] = new_bg
            write_snap.bonds.typeid[:] = new_bt

        for grp_name in ("angles", "dihedrals", "impropers"):
            src = getattr(read_snap, grp_name)
            dst = getattr(write_snap, grp_name)
            if int(src.N) > 0:
                dst.N = int(src.N)
                dst.types = list(src.types)
                dst.group[:] = np.asarray(src.group)
                dst.typeid[:] = np.asarray(src.typeid)

        sim.state.set_snapshot(write_snap)
        self._steps_run += 1

    # ------------------------------------------------------------------
    # Introspection
    # ------------------------------------------------------------------
    @property
    def n_engaged(self) -> int:
        return int((self._head_bound_to_actin >= 0).sum())

    @property
    def n_engaged_alpha(self) -> int:
        m = (self._head_bound_to_actin >= 0) & (self._head_species == "alpha")
        return int(m.sum())

    @property
    def n_engaged_filamin(self) -> int:
        m = (self._head_bound_to_actin >= 0) & (self._head_species == "filamin")
        return int(m.sum())

    @property
    def n_break_total(self) -> int:
        return self._n_break_total

    @property
    def n_bind_total(self) -> int:
        return self._n_bind_total

    @property
    def steps_run(self) -> int:
        return self._steps_run


def make_xlink_updater(
    *,
    p: ResolvedCrosslinkers,
    layout: XlinkLayout,
    kT: float,
    n_cortex_actin: int,
    seed_offset: int = 2,
) -> tuple[XlinkBondUpdater, hoomd.update.CustomUpdater]:
    """Build the XlinkBondUpdater Action wrapped in a periodic CustomUpdater."""
    action = XlinkBondUpdater(
        p=p, layout=layout, kT=kT,
        n_cortex_actin=n_cortex_actin, seed_offset=seed_offset,
    )
    updater = hoomd.update.CustomUpdater(
        action=action, trigger=hoomd.trigger.Periodic(p.batch_steps)
    )
    return action, updater


# ---------------------------------------------------------------------------
# Frame builder: extend a cortex frame with xlink_head particles + intra bonds
# ---------------------------------------------------------------------------
def extend_cortex_state_with_xlinks(
    cortex_snap,
    layout: XlinkLayout,
    p_xl: ResolvedCrosslinkers,
):
    """Append xlink_head particles + xlink_intra bonds + attach-bin bond types.

    No initial attach bonds are placed; the dynamic XlinkBondUpdater will
    bind on its first tick. Initial state has all heads unbound.

    Parameters
    ----------
    cortex_snap : gsd.hoomd.Frame
        The cortex frame from ``build_cortex_state``.
    layout : XlinkLayout
        Per-xlink tag pairs + species (from ``generate_xlink_layout``).
    p_xl : ResolvedCrosslinkers
        For n_bins / max_bind_dist used to register attach-bin types.

    Returns
    -------
    gsd.hoomd.Frame
        New frame with xlink_head particles + intra bonds.
    """
    import gsd.hoomd

    snap_old = cortex_snap
    n_actin = int(snap_old.particles.N)
    n_xl = int(p_xl.n_xl)
    n_heads = 2 * n_xl
    n_part_new = n_actin + n_heads

    snap = gsd.hoomd.Frame()
    # ---- Particles: actin_cortex + xlink_head ----
    snap.particles.N = n_part_new
    snap.particles.types = list(snap_old.particles.types) + ["xlink_head"]
    head_typeid = len(snap_old.particles.types)
    typeid_new = np.empty(n_part_new, dtype=np.uint32)
    typeid_new[:n_actin] = np.asarray(snap_old.particles.typeid)
    typeid_new[n_actin:] = head_typeid
    snap.particles.typeid = typeid_new

    pos_new = np.empty((n_part_new, 3), dtype=np.float64)
    pos_new[:n_actin] = np.asarray(snap_old.particles.position)
    if n_heads > 0:
        pos_new[n_actin:] = layout.head_positions
    snap.particles.position = pos_new

    mass_new = np.empty(n_part_new, dtype=np.float64)
    mass_new[:n_actin] = np.asarray(snap_old.particles.mass)
    mass_new[n_actin:] = 1.0
    snap.particles.mass = mass_new

    # ---- Bonds: cortex backbone + (legacy static xl bins from cortex.py) +
    #             xlink_intra + dynamic attach bins ----
    old_bg = np.asarray(snap_old.bonds.group, dtype=np.uint32)
    old_bt = np.asarray(snap_old.bonds.typeid, dtype=np.uint32)
    old_types = list(snap_old.bonds.types)

    # Append new bond types: xlink_intra + xlink_attach_b{i}
    bond_types_new = list(old_types) + ["xlink_intra"]
    attach_typeid_start = len(bond_types_new)
    bond_types_new += xlink_attach_bin_names(p_xl.n_bins)
    xlink_intra_typeid = len(old_types)

    # Build n_xl intra bonds (head_a ↔ head_b)
    if n_xl > 0:
        intra_groups = layout.head_tag_pairs.astype(np.uint32)
        intra_typeids = np.full(n_xl, xlink_intra_typeid, dtype=np.uint32)
    else:
        intra_groups = np.empty((0, 2), dtype=np.uint32)
        intra_typeids = np.empty((0,), dtype=np.uint32)

    new_bg = np.concatenate([old_bg, intra_groups], axis=0)
    new_bt = np.concatenate([old_bt, intra_typeids])
    snap.bonds.N = int(new_bg.shape[0])
    snap.bonds.types = bond_types_new
    snap.bonds.group = new_bg
    snap.bonds.typeid = new_bt

    # ---- Angles: cortex backbone angles only (no xlink angles) ----
    if int(snap_old.angles.N) > 0:
        snap.angles.N = int(snap_old.angles.N)
        snap.angles.types = list(snap_old.angles.types)
        snap.angles.typeid = np.asarray(snap_old.angles.typeid)
        snap.angles.group = np.asarray(snap_old.angles.group)

    snap.configuration.box = list(snap_old.configuration.box)
    return snap


def build_cortex_xlink_simulation(
    p_cortex,
    p_xl: ResolvedCrosslinkers,
    *,
    device=None,
    with_baoab: bool = True,
    rng: np.random.Generator | None = None,
):
    """End-to-end builder: cortex frame + xlinks + HOOMD wiring + Updater.

    Mirrors ``cortex.build_cortex_simulation`` but extends with the
    dynamic xlink layer.

    Returns
    -------
    sim, baoab_updater, baoab_action, xlink_updater, xlink_action, topology, layout
    """
    # Lazy import to avoid circular dependency at module load.
    from ffn_sim.cortex.cortex import (
        build_cortex_state,
        generate_cortex_topology,
    )
    from ffn_sim.integrator.baoab import make_baoab_updater

    topology = generate_cortex_topology(p_cortex, rng=rng)
    cortex_snap, _, _ = build_cortex_state(
        p_cortex, with_crosslinkers=False, rng=rng
    )

    # Build xlink layout from cortex topology.
    n_cortex_beads = p_cortex.n_filaments * p_cortex.beads_per_filament
    cortex_positions = topology.positions.reshape(n_cortex_beads, 3)
    cortex_filament_idx = np.repeat(
        np.arange(p_cortex.n_filaments, dtype=np.int64),
        p_cortex.beads_per_filament,
    )
    layout = generate_xlink_layout(
        cortex_positions, cortex_filament_idx, p_xl,
        n_cortex_beads=n_cortex_beads, rng=rng,
    )

    snap = extend_cortex_state_with_xlinks(cortex_snap, layout, p_xl)

    sim = hoomd.Simulation(
        device=device or hoomd.device.CPU(), seed=p_cortex.seed
    )
    sim.create_state_from_snapshot(snap)

    # ---- Forces ----
    bond = md.bond.Harmonic()
    bond.params["cortex-bond"] = dict(k=p_cortex.bond_k, r0=p_cortex.rest_length)
    # Legacy static xl bins from cortex.py (in case xl_enabled there) are
    # not active when xl_enabled=False in cortex.py, so cortex_snap has
    # only cortex-bond. extend_cortex_state preserves any existing types.
    # The dynamic xlink_intra + xlink_attach_b{i} need params:
    bond.params["xlink_intra"] = dict(k=p_xl.k_intra, r0=0.0)
    # Per-xlink intra rest length varies (alpha 35 nm vs filamin 150 nm)
    # but md.bond.Harmonic.params is per-TYPE, not per-bond. We use a
    # single average rest length, accepting that this is a refinement
    # gap (next iteration: split xlink_intra into xlink_intra_alpha +
    # xlink_intra_filamin types).
    avg_intra_r0 = float(
        p_xl.alpha_fraction * p_xl.alpha_length
        + (1.0 - p_xl.alpha_fraction) * p_xl.filamin_length
    )
    bond.params["xlink_intra"] = dict(k=p_xl.k_intra, r0=avg_intra_r0)
    bin_r0 = xlink_attach_bin_rest_lengths(p_xl.n_bins, p_xl.max_bind_dist)
    for i, name in enumerate(xlink_attach_bin_names(p_xl.n_bins)):
        bond.params[name] = dict(k=p_xl.k_attach, r0=float(bin_r0[i]))

    angle = md.angle.Harmonic()
    angle.params["cortex-angle"] = dict(k=p_cortex.angle_k, t0=p_cortex.angle_t0)

    nlist = md.nlist.Tree(buffer=0.5 * p_cortex.lj_sigma)
    lj = md.pair.LJ(nlist=nlist, default_r_cut=0.0)
    # actin_cortex × actin_cortex (D7 WCA, same as cortex.py)
    lj.params[("actin_cortex", "actin_cortex")] = dict(
        epsilon=p_cortex.lj_epsilon, sigma=p_cortex.lj_sigma
    )
    lj.r_cut[("actin_cortex", "actin_cortex")] = (
        p_cortex.lj_r_cut if p_cortex.lj_enabled else 0.0
    )
    # xlink_head × xlink_head: WCA repulsive (same params as actin × actin)
    lj.params[("xlink_head", "xlink_head")] = dict(
        epsilon=p_cortex.lj_epsilon, sigma=p_cortex.lj_sigma
    )
    lj.r_cut[("xlink_head", "xlink_head")] = (
        p_cortex.lj_r_cut if p_cortex.lj_enabled else 0.0
    )
    # xlink_head × actin_cortex: NO LJ (heads bind to actin via dynamic
    # attach bonds; LJ would prevent close approach for binding).
    lj.params[("xlink_head", "actin_cortex")] = dict(
        epsilon=p_cortex.lj_epsilon, sigma=p_cortex.lj_sigma
    )
    lj.r_cut[("xlink_head", "actin_cortex")] = 0.0
    lj.mode = "shift"

    ig = md.Integrator(dt=p_cortex.dt_cfl)
    ig.forces.append(bond)
    ig.forces.append(angle)
    ig.forces.append(lj)
    sim.operations.integrator = ig

    baoab_updater = None
    baoab_action = None
    if with_baoab:
        baoab_action, baoab_updater = make_baoab_updater(
            kT=p_cortex.kT,
            gamma={"actin_cortex": p_cortex.gamma_b, "xlink_head": p_cortex.gamma_b},
            dt=p_cortex.dt_cfl,
            seed=p_cortex.seed,
        )
        sim.operations.updaters.append(baoab_updater)

    xlink_action, xlink_updater = make_xlink_updater(
        p=p_xl, layout=layout, kT=p_cortex.kT,
        n_cortex_actin=n_cortex_beads,
    )
    sim.operations.updaters.append(xlink_updater)

    return sim, baoab_updater, baoab_action, xlink_updater, xlink_action, topology, layout
