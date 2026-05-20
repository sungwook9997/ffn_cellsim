"""H.1 ECM Mikado cross-link seeding — HOOMD wiring (KU-1.27, KU-1.28).

REUSES the v1 oracle's geometric segment-intersection finder
``ffn_sim/validation/oracles/ecm/cross_links._segment_intersections``
and ``generate_cross_links`` (Phase-1, primary-copy only — no periodic
image sweep, per oracle's documented constraint). The oracle's
**force kernel** ``compute_xl_energy_and_forces`` is NOT imported here
(CLAUDE.md hard rule).

Each cross-link is realised as a HOOMD bond of type ``"xl"`` between
the nearest bead on fiber A and the nearest bead on fiber B at the
intersection. Force form (KU-1.28):

    U_xl = ½ k_xl (|r_b − r_a| − r0)²,   k_xl = 1·10⁻³ N/m  (1 pN/nm)

H.1 brief §Topology specifies ``r0 = distance at intersection (≈ 0)``.
We follow the brief's literal simplification (``r0 = 0`` for all xl
bonds, single bond type) rather than the oracle's per-link rest_length
= current MI bead-bead distance:

- HOOMD's ``md.bond.Harmonic`` exposes only *per-type* (k, r0); per-link
  r0 would require either ~20 binned bond types or a custom Python
  force compute (slow per step).
- With Δt = α·τ_min ≈ 3.8 ns and ``k_xl/γ_b ≈ 1.5 ns⁻¹``, an
  initially-stretched xl at the maximum bead-bead offset ~ℓ₀/√2 ≈
  350 nm relaxes to within thermal noise in O(τ_xl/Δt) = O(170) steps —
  i.e. the construction-time stress is wiped out within the KU-1.30
  equilibration prelude. So the simplification does NOT change the
  long-time mechanics.
- Energy-oracle agreement at construction (the H.1 brief §Validation
  KU-1.24 gate) is therefore tested on the *ecm-bond + angle* subset
  only; the xl subset has its own topology smoke gate (count + bead
  indices) but not an energy-vs-oracle gate. This mismatch with the
  oracle convention is surfaced to PI in the M2 closeout.

Sanity Gate
-----------
*Per CLAUDE.md hard rule "Sanity Gate Protocol mandatory before first
execution of any physics/numerics module."*

1. **Dimensional analysis**

   - ``k_xl`` [N/m] · (length [m])² → [N·m] = [J]. ✓
   - ``F_xl = k_xl · |r_b − r_a|`` at r0=0 has units [N]. ✓
   - Per-bead Stokes velocity ``F/γ_b`` is [m/s]. ``F/γ_b · Δt`` is [m]
     so the per-step xl-driven displacement is well-typed. ✓
   - τ_xl = γ_b / k_xl, already verified in mikado.py's resolve_derived.

2. **Boundary cases**

   - **No intersections** (e.g. n_fibers=1 or all parallel): returns an
     empty bond list; the snapshot is unchanged. RUNTIME pass-through.
   - **Self-intersection** (fiber crossing itself across a periodic
     boundary): excluded by ``i < j`` in the oracle finder. ✓
   - **Degenerate intersection** on a bead (s_a or s_b lands exactly at
     a discrete bead): the bead index is well-defined by ``np.round``.
     Both fibers may share the same bead → forbidden by ``i ≠ j``
     fiber-pair constraint at the oracle level; never seen for
     non-self-touching Mikado.
   - **Two intersections placing the same (fiber, bead) twice**: would
     create a multi-bond on the same vertex. RUNTIME §6 measurement
     detects duplicates and raises (we deduplicate by sorted
     (global_idx_a, global_idx_b) tuple; if dedup count differs from
     raw count, we surface to PI).

3. **Conservation invariants**

   - **Newton 3rd law**: HOOMD's ``md.bond.Harmonic`` is symmetric per
     bond. Σ F_xl over the two bonded beads of a single xl is exactly
     zero by construction. STATIC test ``test_h1_cross_links.py
     ::TestSignSense::test_xl_pair_newton_third`` (single-bond
     two-bead simulation).
   - **Topology mutation count**: ``N_bonds_total = N_ecm_bonds +
     N_xl_bonds``. ``N_xl_bonds`` must equal the oracle's
     ``len(generate_cross_links(...))``. STATIC.
   - **Particle count unchanged**: cross_links add bonds only; no new
     particles. STATIC.

4. **Numerical sanity**

   - Bond-group indices stay in ``[0, N_particles)``: enforced when
     converting (fiber, bead) → global = fiber·beads_per_fiber + bead;
     we assert ``max(group) < N`` before assigning.
   - Rest length stored as float64.
   - Stiffness ``k_xl`` finite-positive (oracle config + resolve_derived).
   - No xl bond between a bead and itself (i.e. same global index in
     both columns of a bond.group row). Trivially excluded by the
     ``fiber_a ≠ fiber_b`` constraint at the oracle level.

5. **Sign / sense**

   - At construction, two xl-bonded beads are at MI distance ``d``;
     with r0=0, the harmonic force pulls them together along the MI
     displacement (attractive). STATIC test:
     ``test_xl_pair_attractive_force`` — single xl bond between two
     widely-separated beads, confirm force on each points along the
     bond axis toward the other.
   - For ``d → 0`` (overlap), the harmonic force → 0 (linear regime).
     No singular repulsion at contact — the LJ WCA already handles
     hard-core exclusion separately (D7).

6. **Measurement protocol**

   - **Topology smoke**: ``len(xl_bonds)`` matches
     ``len(generate_cross_links(net, k_xl))``. STATIC.
   - **Coordination check**: ``measure_coordination(xl_bonds,
     n_beads_total, beads_per_fiber)`` reports the achieved ⟨z⟩.
     If outside the brief's acceptance range
     ``[2.5, 3.9]`` (KU-1.3 sub-isostatic), raise unless demo_mode.
     RUNTIME.
   - The xl bond count is bounded by the oracle's
     ``expected_xl_per_fiber = 2 ρ_L L_fiber / π`` × n_fibers / 2,
     within a finite-size correction. STATIC tolerance ±15% (Phase-1
     primary-copy-only construction misses ~5% of intersections that
     cross periodic boundaries; see oracle docstring).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import gsd.hoomd

from ffn_sim.ecm.mikado import ResolvedH1
from ffn_sim.validation.oracles.ecm.cross_links import (
    CrossLink as OracleCrossLink,
    generate_cross_links as oracle_generate_cross_links,
    measure_coordination as oracle_measure_coordination,
)
from ffn_sim.validation.oracles.ecm.fiber_network import (
    FiberNetwork,
    generate_2d_fiber_network,
)


XL_BOND_TYPE_PREFIX: str = "xl_b"
# Number of r0 bins for cross-link rest-length quantization (PI 2026-05-20,
# option B). Bin centers at {0, w, 2w, ..., n·w} with w = XL_BIN_WIDTH_M.
# With L_fiber/2 = 5 μm as the geometric ceiling on bead-to-intersection
# offset, every realised per-link MI distance falls in [0, ~ℓ₀] = [0, 500 nm].
# 11 bin centers at 50 nm spacing give quantization error ≤ 25 nm per link
# and ~0.07 kT residual per-link energy at construction (≈100× reduction
# vs the brief's r0=0 simplification).
XL_BIN_WIDTH_M: float = 50.0e-9          # 50 nm
XL_N_BINS: int = 11                       # bin centers 0, 50, ..., 500 nm

# Public-facing single name kept for backward-compat with M2-pre tests that
# wired one xl type; modern code should iterate xl_bin_type_names(...).
XL_BOND_TYPE_NAME: str = XL_BOND_TYPE_PREFIX + "0"


def xl_bin_type_names() -> list[str]:
    """Return the canonical ordered bin type names ``[xl_b0, xl_b1, ...]``."""
    return [f"{XL_BOND_TYPE_PREFIX}{i}" for i in range(XL_N_BINS)]


def xl_bin_rest_lengths() -> np.ndarray:
    """Quantized r0 (m) for each bin index 0..XL_N_BINS-1."""
    return XL_BIN_WIDTH_M * np.arange(XL_N_BINS, dtype=np.float64)


def quantize_rest_lengths(rest: np.ndarray) -> np.ndarray:
    """Map each per-link rest length to a bin index in [0, XL_N_BINS).

    Standard "round to nearest bin center" quantization; clamps the upper
    tail to ``XL_N_BINS - 1`` so an over-range link (rare, > XL_N_BINS·w)
    still has a valid bin assignment. Max quantization error ≤ w/2.
    """
    idx = np.round(np.asarray(rest, dtype=np.float64) / XL_BIN_WIDTH_M)
    idx = np.clip(idx, 0, XL_N_BINS - 1).astype(np.int64)
    return idx


@dataclass(slots=True)
class XLBonds:
    """Cross-link bond data ready to merge into a HOOMD GSD frame.

    Attributes
    ----------
    group : np.ndarray, shape (N_xl, 2), dtype uint32
        Per-link (global_idx_a, global_idx_b) pairs.
    rest_lengths : np.ndarray, shape (N_xl,), dtype float64
        Per-link MI bead-bead distance at construction (oracle
        convention). Diagnostic / quantization input; the actual r0
        written to HOOMD is the *bin center* `XL_BIN_WIDTH_M · bin_idx`.
    bin_idx : np.ndarray, shape (N_xl,), dtype int64
        Per-link bin index in ``[0, XL_N_BINS)``. Maps directly to bond
        ``typeid`` when the frame is assembled.
    k_xl : float
        Per-bin-type harmonic stiffness (KU-1.28, default 1e-3 N/m). The
        same ``k`` is used across all xl bins; only r0 varies.
    """

    group: np.ndarray
    rest_lengths: np.ndarray
    bin_idx: np.ndarray
    k_xl: float


def _global_index(fiber: np.ndarray, bead: np.ndarray, beads_per_fiber: int) -> np.ndarray:
    return (fiber.astype(np.int64) * beads_per_fiber
            + bead.astype(np.int64)).astype(np.uint32)


def generate_xl_bonds(p: ResolvedH1) -> tuple[XLBonds, FiberNetwork]:
    """Build the cross-link bond list for the H.1 Mikado defined by `p`.

    Internally re-runs the oracle's ``generate_2d_fiber_network`` with
    the exact seed used by ``build_mikado_state`` (so the topology is
    identical bit-for-bit) and then calls the oracle's
    ``generate_cross_links`` for the segment-intersection geometry.

    Returns
    -------
    XLBonds
        Aggregated bond group + per-link rest lengths + k_xl.
    FiberNetwork
        The same network object the cross-link finder ran on (returned
        so that callers — typically build_mikado_simulation — don't have
        to re-generate it).
    """
    net = generate_2d_fiber_network(
        L_box=p.L_box,
        n_fibers=p.n_fibers,
        L_fiber=p.L_fiber,
        beads_per_fiber=p.beads_per_fiber,
        S_order=0.0,
        theta0=0.0,
        seed=p.seed,
        params={"resolved_from": "ffn_sim/configs/phase1_h1.yaml"},
    )
    oracle_links: list[OracleCrossLink] = oracle_generate_cross_links(
        net, stiffness=p.xl_stiffness
    )
    if not oracle_links:
        empty_group = np.empty((0, 2), dtype=np.uint32)
        empty_rest = np.empty((0,), dtype=np.float64)
        empty_idx = np.empty((0,), dtype=np.int64)
        return XLBonds(
            group=empty_group,
            rest_lengths=empty_rest,
            bin_idx=empty_idx,
            k_xl=p.xl_stiffness,
        ), net

    fa = np.fromiter((xl.fiber_a for xl in oracle_links), dtype=np.int64,
                     count=len(oracle_links))
    ba = np.fromiter((xl.bead_a for xl in oracle_links), dtype=np.int64,
                     count=len(oracle_links))
    fb = np.fromiter((xl.fiber_b for xl in oracle_links), dtype=np.int64,
                     count=len(oracle_links))
    bb = np.fromiter((xl.bead_b for xl in oracle_links), dtype=np.int64,
                     count=len(oracle_links))
    rest = np.fromiter((xl.rest_length for xl in oracle_links),
                       dtype=np.float64, count=len(oracle_links))

    ia = _global_index(fa, ba, p.beads_per_fiber)
    ib = _global_index(fb, bb, p.beads_per_fiber)

    # Sanity Gate §4: indices in range; no self-bonds; dedup duplicates.
    N = p.n_fibers * p.beads_per_fiber
    if (ia >= N).any() or (ib >= N).any():
        raise RuntimeError(
            "Cross-link global index out of range [0, N_particles). "
            f"max(ia)={int(ia.max())}, max(ib)={int(ib.max())}, N={N}."
        )
    if (ia == ib).any():
        raise RuntimeError(
            "Cross-link with identical endpoints (self-bond); oracle "
            "should have excluded these via fiber_a ≠ fiber_b."
        )

    # Sanity Gate §2/§6: deduplicate symmetric (a,b) ≡ (b,a) bonds.
    pair_lo = np.minimum(ia, ib)
    pair_hi = np.maximum(ia, ib)
    unique_view = pair_lo.astype(np.uint64) * np.uint64(N) + pair_hi.astype(np.uint64)
    _, unique_idx = np.unique(unique_view, return_index=True)
    if unique_idx.size != ia.size:
        # Surface to PI per §2 — duplicate xls on the same vertex pair.
        raise RuntimeError(
            f"Duplicate cross-link pairs detected ({ia.size} raw vs "
            f"{unique_idx.size} unique); surface to PI."
        )

    group = np.stack([ia, ib], axis=-1).astype(np.uint32)
    bin_idx = quantize_rest_lengths(rest)
    return (
        XLBonds(
            group=group,
            rest_lengths=rest,
            bin_idx=bin_idx,
            k_xl=p.xl_stiffness,
        ),
        net,
    )


def add_xl_to_frame(snap: gsd.hoomd.Frame, xl: XLBonds) -> gsd.hoomd.Frame:
    """Append cross-link bonds (one bond type per r0 bin) to ``snap``.

    The frame's existing bonds (ecm-bond) are preserved at type-id 0.
    xl bin types follow as 'xl_b0' (id 1), 'xl_b1' (id 2), ... up to
    'xl_b{XL_N_BINS-1}'. All xl bin types are registered even if some
    bins receive zero links, so that downstream code can always look up
    ``bond.params["xl_b<i>"]`` without conditional checks.
    """
    n_existing = int(snap.bonds.N)
    existing_group = np.asarray(snap.bonds.group, dtype=np.uint32)
    existing_typeid = np.asarray(snap.bonds.typeid, dtype=np.uint32)
    existing_types = list(snap.bonds.types)

    bin_names = xl_bin_type_names()
    # Type IDs for the xl bins: assigned right after existing types.
    base_id = len(existing_types)
    new_types = existing_types + bin_names

    if xl.group.shape[0] == 0:
        snap.bonds.types = new_types
        return snap

    # Per-link typeid = base_id + bin_idx.
    xl_typeids = (base_id + xl.bin_idx).astype(np.uint32)

    snap.bonds.N = n_existing + xl.group.shape[0]
    snap.bonds.types = new_types
    snap.bonds.typeid = np.concatenate([existing_typeid, xl_typeids])
    snap.bonds.group = np.concatenate([existing_group, xl.group], axis=0)
    return snap


def n_xl_bonds(xl: XLBonds) -> int:
    """Cross-link bond count (independent of bin distribution)."""
    return int(xl.group.shape[0])


def measure_coordination(xl: XLBonds, p: ResolvedH1) -> float:
    """Achieved ⟨z⟩ including ecm-bond backbone + xl contributions."""
    fake_links: list[OracleCrossLink] = [
        OracleCrossLink(0, 0, 0, 0, rest_length=0.0, stiffness=xl.k_xl)
        for _ in range(xl.group.shape[0])
    ]
    return oracle_measure_coordination(
        fake_links,
        n_beads_total=p.n_fibers * p.beads_per_fiber,
        beads_per_fiber=p.beads_per_fiber,
    )
