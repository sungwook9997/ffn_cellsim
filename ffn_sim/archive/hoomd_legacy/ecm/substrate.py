"""ECM substrate compliance — tunable ``k_sub`` anchor spring (KU-1.V.1).

ADDITIVE, DEFAULT-OFF module (Track A; see
``docs/ECM_PLATFORM_EXTENSIBILITY.md`` §"The key unlock"). Today the FA
substrate-ligand layer is held immobile by
:class:`ffn_sim.archive.hoomd_legacy.cell.cell.SubstrateLigandPin`, which **resets every ligand
to its construction position every step** — i.e. an *infinitely* rigid
backing (GPa-equivalent), the closest representation to the PI experimental
condition (b) collagen-I coated on a rigid dish. That makes substrate
stiffness a *binary* axis (rigid pin OR free), with no way to represent the
in-vitro→in-vivo compliance continuum (PA/PEG gels 0.1–40 kPa, tissue,
tumor stroma).

This module supplies the missing **finite-stiffness** endpoint of that axis:
a tunable anchor harmonic spring of stiffness ``k_sub`` [N/m] from each
ligand to its captured construction (z=0) site::

    U_i = ½ k_sub |r_i − a_i|²
    F_i = −k_sub (r_i − a_i)

where ``a_i`` is the construction position of ligand ``i`` (the z=0 anchor
the rigid pin would otherwise reset to). Sweeping ``k_sub`` makes substrate
compliance a **continuous axis** with the rigid coating and the soft gel as
the two endpoints.

Relationship to the rigid pin (CRITICAL — not bit-for-bit identical)
--------------------------------------------------------------------
The rigid default stays the EXISTING ``SubstrateLigandPin`` path. This
module is ONLY the finite-``k_sub`` spring and is a SEPARATE mechanism:

* :class:`SubstrateLigandPin` is a *post-BAOAB position reset* — it
  overwrites the ligand row's position each step, adding NO force to
  ``net_force``. A pinned ligand's per-step net motion is exactly zero.
* :class:`SubstrateAnchorSpring` is a *force* contribution to
  ``net_force`` (a ``md.force.Custom``). Even a very large ``k_sub`` only
  makes the ligand *stiff*, not *immobile*: the BAOAB heat bath still
  integrates it, so a finite-large ``k_sub`` is **NOT** bit-for-bit
  identical to the position-reset pin (it leaves a residual thermal
  excursion ``σ = √(kT/k_sub)`` and injects a force the pin never did).

The ``k_sub → ∞`` limit recovers the rigid-pin behaviour *in the physics
sense* (σ → 0, the ligand is asymptotically clamped at its anchor), and the
two endpoints bracket the compliance axis — but the BIT-IDENTICAL rigid
backing is ``pin: true`` / the ``SubstrateLigandPin`` path, NOT a large
``k_sub`` here. Attaching this spring is therefore a deliberate
"compliant-substrate" selection: a caller wanting the bit-identical legacy
run must NOT attach it and must keep the pin. (Mechanism-recon note,
2026-05-30 / 31.)

The two paths are also mutually exclusive on the same ligand tags: a pin
that resets position AND a spring that pulls toward the same anchor would
double-handle the ground. The preset layer selects exactly one
(rigid-pin XOR finite-spring) per substrate.

Dimensional bridge ``k_sub → E_sub`` (documented, NOT hard-coded)
-----------------------------------------------------------------
A per-ligand spring stiffness ``k_sub`` [N/m] maps to an effective
substrate Young's modulus via the ligand areal density (KU-1.V.1)::

    E_sub  ≈  k_sub · n_ligand / A_substrate          [N/m · 1/m² ·? ]

i.e. ``E_sub`` scales like the spring stiffness per unit ligand-projected
area (a discrete-spring → continuum-modulus bridge; the exact prefactor
depends on the indentation / contact model used as the calibration oracle,
e.g. flat-punch vs Hertzian, against the KU-1.V.1 bands). This module
deliberately does **NOT** hard-code an ``E_sub`` or the bridge prefactor:
``k_sub`` is the runtime knob, ``E_sub`` is a *derived diagnostic* the
preset/calibration layer computes once it fixes ``n_ligand``,
``A_substrate`` and the contact-model oracle. A first-order bridge helper
:func:`effective_E_sub` is provided for that diagnostic, flagged as a TODO
pending the calibration-oracle ratification (no gate depends on it).

Magic-Number Block — k_sub (no empirical magic number)
------------------------------------------------------
``k_sub`` is a **REQUIRED** parameter with **NO default**. It is the ECM
condition's substrate stiffness and is supplied by the preset layer per
condition (rigid dish → ∞ / pin; PA-gel 0.1–40 kPa → finite ``k_sub``
calibrated against KU-1.V.1). There is deliberately no in-module default:
inventing one would be an empirical magic number, and the rigid endpoint is
already owned by the pin path. The resolver therefore *requires* ``k_sub``
and raises if it is absent or non-finite-positive.

Sanity Gate
-----------
*Per CLAUDE.md Hard Rule. STATIC + analytical checks in
``ffn_sim/tests/test_substrate_spring.py``.*

1. **Dimensional analysis**
   - ``k_sub`` [N/m]; ``(r_i − a_i)`` [m]; ``F_i = k_sub · Δr`` [N]. ✓
   - ``U_i = ½ k_sub |Δr|²`` [N/m · m² = J]. ✓
   - At construction (``r_i = a_i``): ``Δr = 0`` → ``F = 0``, ``U = 0``.
     Force-free at construction (this is the zero-force configuration).
   - Bridge: ``E_sub ≈ k_sub · n_ligand / A_substrate``
     [N/m · 1/m²] → modulus-like; documented, not hard-coded.

2. **Boundary cases**
   - ``k_sub ≤ 0`` or non-finite: caught by :func:`resolve_substrate` /
     the force constructor ValueError.
   - ``k_sub`` absent: REQUIRED — :func:`resolve_substrate` raises
     (no magic default).
   - Empty ligand set (n = 0): the custom force iterates zero ligand
     particles; trivially zero force.
   - A ligand whose tag is NOT in the captured anchor map: the force
     constructor raises (every masked tag must have an anchor).
   - ``r_i = a_i`` (the zero-force config): ``F = 0`` exactly (no
     direction singularity — unlike the radial ERM/EV fields, the anchor
     spring uses the full displacement vector ``r − a``, which is simply
     the zero vector here).

3. **Conservation invariants**
   - Net force from this spring on the system: NOT zero (it is an
     EXTERNAL field anchored at fixed lab-frame sites ``a_i`` — exactly
     like ERM, and exactly the mechanical "ground" the pin represented).
     Intentional: the substrate IS the inertial reference.
   - Total potential ``Σ ½ k_sub |r_i − a_i|²`` — additive across ligand
     tags, monotonically non-negative.
   - Force applies ONLY to particles whose tag is in ``ligand_tags``;
     all other particles (integrin, actin, …) are untouched (mask).

4. **Numerical sanity / CFL**
   - All positions / anchors read in float64.
   - CFL: the spring's per-bead stiffness is exactly ``k_sub``. The relax
     time ``τ = γ_ligand / k_sub`` must satisfy ``dt ≤ α · τ``; the attach
     helper raises (parity with ``erm.py`` / ``enclosed_volume.py``) if
     violated. A stiff ``k_sub`` therefore tightens CFL — the caller must
     pick ``k_sub`` (and/or ``dt``) so ``dt ≤ cfl_safety_factor·γ/k_sub``.
   - Force compute matches the SI analytic ``k_sub·|Δr|`` to float64
     precision (STATIC test: displace one ligand by δ, read its force).

5. **Sign / sense**
   - Ligand displaced AWAY from anchor (``r_i ≠ a_i``): ``F_i = −k_sub
     (r_i − a_i)`` points back TOWARD ``a_i`` — restoring. STATIC test
     (projection of F on the outward displacement is negative).
   - Ligand AT anchor: ``F = 0``. STATIC.
   - ``k_sub → ∞`` recovers the rigid-pin limit (σ = √(kT/k_sub) → 0):
     documented; note finite-large ``k_sub`` ≠ pin bit-for-bit (§above).

6. **Measurement protocol**
   - Per-ligand displacement at equilibrium is centered at ``a_i`` with
     std-dev ``σ = √(kT/k_sub)`` (the spring is a 3D harmonic trap). A
     stiffer ``k_sub`` → smaller σ → more rigid substrate; the ``k_sub →
     ∞`` (σ → 0) endpoint is the rigid-dish (b) limit. STATIC test checks
     the σ formula + that the BAOAB-on spring keeps ligands bounded near
     their anchors over a short run.

References
----------
- Plan: ``ffn_sim/docs/ECM_PLATFORM_EXTENSIBILITY.md`` §"The key unlock",
  §"Condition matrix", §"Platform prep" (``ecm/substrate.py`` row).
- KU-1.V.1 substrate-stiffness ECM-V layer (rigid coated dish ~GPa cap;
  PA/PEG gel 0.1–40 kPa sweep; tissue 0.1–50 kPa; YAP transition ~5 kPa) —
  the ``k_sub → E_sub`` calibration bands (Notion KU v2 ECM-V layer).
- Rigid endpoint / legacy bit-identical path:
  ``ffn_sim/cell/cell.py`` :class:`SubstrateLigandPin` (position reset).
- Structural template: ``ffn_sim/cortex/erm.py``
  (:class:`ERMHarmonic` — external harmonic field on a tag range +
  Sanity Gate §1–6 + CFL attach gate) and
  ``ffn_sim/cortex/enclosed_volume.py`` (CFL gate prose).
- HOOMD 7 ``md.force.Custom`` API.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

import hoomd
import hoomd.md as md


# ---------------------------------------------------------------------------
# Resolved parameters
# ---------------------------------------------------------------------------
@dataclass(slots=True)
class ResolvedSubstrate:
    """Tunable-compliance substrate-anchor parameters (all SI).

    Attributes:
        k_sub: Per-ligand anchor-spring stiffness [N/m]. REQUIRED — the ECM
            condition's substrate stiffness (no magic default; the rigid
            endpoint is the ``SubstrateLigandPin`` path, not a large value
            here). Must be finite and > 0.
        sigma_thermal: Derived per-ligand thermal trap std-dev
            ``√(kT/k_sub)`` [m]. A stiffer ``k_sub`` → smaller σ → more
            rigid substrate; σ → 0 is the rigid-dish (b) limit. Populated by
            :func:`resolve_substrate`.
    """

    k_sub: float                       # N/m  REQUIRED anchor-spring stiffness
    sigma_thermal: float = 0.0         # m    √(kT/k_sub)  (derived diagnostic)


def _require_finite_positive(name: str, x: float) -> None:
    """Raise ``ValueError`` unless ``x`` is finite and strictly positive."""
    if not (math.isfinite(x) and x > 0.0):
        raise ValueError(f"{name} must be finite and > 0; got {x!r}")


def resolve_substrate(cfg: dict, *, kT: float) -> ResolvedSubstrate:
    """Resolve the ``substrate`` (tunable-compliance) config block.

    ``k_sub`` is **REQUIRED** (no default — see the Magic-Number Block in the
    module docstring): it is the ECM condition's substrate stiffness, supplied
    by the preset layer. The rigid (b)-dish endpoint is the
    ``SubstrateLigandPin`` position-reset path, NOT a large ``k_sub`` here, so
    there is deliberately no in-module default to invent.

    Args:
        cfg: The YAML root, an ``ecm``/``ecm_condition`` sub-dict, or a
            ``substrate`` sub-dict. Must contain ``k_sub`` [N/m].
        kT: Thermal energy [J] (from the host run's resolved params), used to
            derive the diagnostic trap width ``σ = √(kT/k_sub)``.

    Returns:
        The resolved substrate parameters with ``sigma_thermal`` populated.

    Raises:
        ValueError: If ``k_sub`` is absent, non-finite, or ≤ 0, or if ``kT``
            is non-finite-positive.
    """
    if "ecm_condition" in cfg:
        cfg = cfg["ecm_condition"]
    if "ecm" in cfg:
        cfg = cfg["ecm"]
    if "substrate" in cfg:
        cfg = cfg["substrate"]

    if "k_sub" not in cfg or cfg["k_sub"] is None:
        raise ValueError(
            "k_sub is REQUIRED (no magic default): supply the ECM "
            "condition's substrate stiffness [N/m]. The rigid (b)-dish "
            "endpoint is the SubstrateLigandPin position-reset path, not a "
            "large k_sub here."
        )
    k_sub = float(cfg["k_sub"])
    _require_finite_positive("k_sub", k_sub)
    _require_finite_positive("kT", kT)

    p = ResolvedSubstrate(k_sub=k_sub)
    p.sigma_thermal = math.sqrt(kT / k_sub)
    return p


# ---------------------------------------------------------------------------
# Dimensional-bridge diagnostic (k_sub → E_sub) — NOT used by the runtime
# ---------------------------------------------------------------------------
def effective_E_sub(
    k_sub: float, n_ligand: int, A_substrate: float
) -> float:
    """First-order ``k_sub`` → substrate stiffness DENSITY [N/m³] (DIAGNOSTIC).

    Maps the per-ligand spring stiffness to an effective substrate modulus
    via the ligand areal density (KU-1.V.1)::

        E_sub ≈ k_sub · n_ligand / A_substrate.

    .. warning::
        DIAGNOSTIC ONLY — never called by the runtime force. The exact
        prefactor depends on the contact / indentation model chosen as the
        calibration oracle (flat-punch vs Hertzian) against the KU-1.V.1
        bands; this returns the leading-order spring-stiffness-per-unit-area
        estimate. **TODO (deeper physics):** ratify the calibration oracle
        and prefactor with the PI before using ``E_sub`` in any gate. No gate
        currently depends on this value.

    Args:
        k_sub: Per-ligand anchor-spring stiffness [N/m].
        n_ligand: Number of substrate ligands sharing the contact area.
        A_substrate: Substrate contact (projected) area [m²].

    Returns:
        Leading-order substrate stiffness DENSITY ``k_sub·n_ligand/A`` [N/m³]
        (= N/m · 1/m²). NOTE: this is NOT yet a Pa Young's modulus — it needs a
        contact length-scale factor; ratify the calibration oracle (flat-punch
        vs Hertzian, KU-1.V.1) with the PI before reading it as E_sub [Pa].
        (PI-flag 2026-05-31: the old "[Pa]" label was dimensionally wrong.)

    Raises:
        ValueError: If ``k_sub`` or ``A_substrate`` is non-finite-positive,
            or ``n_ligand`` ≤ 0.
    """
    _require_finite_positive("k_sub", k_sub)
    _require_finite_positive("A_substrate", A_substrate)
    if n_ligand <= 0:
        raise ValueError(f"n_ligand must be > 0; got {n_ligand}")
    return k_sub * n_ligand / A_substrate


# ---------------------------------------------------------------------------
# Custom force compute
# ---------------------------------------------------------------------------
class SubstrateAnchorSpring(md.force.Custom):
    """Per-ligand anchor harmonic spring to captured construction sites.

    Applies a restoring spring force ``F_i = −k_sub (r_i − a_i)`` to every
    particle whose tag is in ``ligand_tags``, where ``a_i`` is that ligand's
    captured construction (z=0) position. All other particles are untouched.
    This is the **finite-stiffness** substrate endpoint; the rigid endpoint
    is the :class:`ffn_sim.archive.hoomd_legacy.cell.cell.SubstrateLigandPin` position-reset path
    (see module docstring — finite-large ``k_sub`` is NOT bit-for-bit the
    pin).

    Mirrors :class:`ffn_sim.archive.hoomd_legacy.cortex.erm.ERMHarmonic` exactly (tag mask,
    per-bead restoring force toward a fixed lab-frame anchor,
    ``cpu_local_force_arrays`` write, ``force[~mask]=0``) — the only
    difference is the anchor is a per-ligand 3D site ``a_i`` (captured at
    construction) rather than a common radius ``R_cell``.

    The anchor is looked up BY TAG (a tag→anchor map gathered onto rows),
    NOT by row index, because HOOMD does not guarantee the local-snapshot row
    order matches the order of the supplied ``anchor_positions`` array. This
    mirrors :class:`SubstrateLigandPin`'s ``_anchor_by_tag`` lookup.

    Args:
        p: Resolved substrate params (provides ``k_sub``).
        ligand_tags: Global HOOMD tags of the substrate-ligand particles.
        anchor_positions: Construction-time anchor sites, shape
            ``(len(ligand_tags), 3)`` [m], row-aligned with ``ligand_tags``.
        aniso: HOOMD anisotropic-force flag (default ``False``; this is a
            translational point force).

    Raises:
        ValueError: If ``k_sub`` is non-finite-positive, if
            ``anchor_positions`` is not ``(n_ligand, 3)``, or if a ligand tag
            is duplicated.

    Notes:
        Implemented as ``md.force.Custom`` so the contribution participates in
        HOOMD's ``net_force`` accumulator, which the frozen L-M BAOAB Updater
        reads each step. Computed in float64 from positions via
        ``cpu_local_snapshot`` → ``cpu_local_force_arrays`` (row-indexed; same
        single-rank guarantee ``erm.py`` documents).
    """

    def __init__(
        self,
        p: ResolvedSubstrate,
        ligand_tags: np.ndarray,
        anchor_positions: np.ndarray,
        *,
        aniso: bool = False,
    ) -> None:
        super().__init__(aniso=aniso)
        _require_finite_positive("k_sub", p.k_sub)
        self.p = p
        self._ligand_tags = np.asarray(ligand_tags, dtype=np.int64)
        anchor = np.asarray(anchor_positions, dtype=np.float64).copy()
        if anchor.shape != (self._ligand_tags.shape[0], 3):
            raise ValueError(
                "anchor_positions must have shape (n_ligand, 3); got "
                f"{anchor.shape} for {self._ligand_tags.shape[0]} tags."
            )
        if not np.isfinite(anchor).all():
            raise ValueError("anchor_positions contains non-finite values.")
        # Tag→anchor map (robust to HOOMD row reordering — mirrors
        # SubstrateLigandPin). Duplicate tags would silently drop anchors.
        if np.unique(self._ligand_tags).size != self._ligand_tags.size:
            raise ValueError("ligand_tags contains duplicate tags.")
        self._anchor = anchor
        self._anchor_by_tag = dict(zip(self._ligand_tags.tolist(), anchor))
        self._tag_set = set(self._ligand_tags.tolist())

    @staticmethod
    def predicted_sigma(k_sub: float, kT: float) -> float:
        """Analytic per-ligand displacement std-dev at thermal equilibrium.

        The anchor spring is a 3D harmonic trap, so each Cartesian component
        is Gaussian with std-dev ``σ = √(kT/k_sub)`` [m]. σ → 0 as
        ``k_sub → ∞`` (the rigid-dish limit).

        Args:
            k_sub: Anchor-spring stiffness [N/m].
            kT: Thermal energy [J].

        Returns:
            Per-component thermal trap width ``√(kT/k_sub)`` [m].
        """
        return math.sqrt(kT / k_sub)

    def set_forces(self, timestep: int) -> None:  # noqa: D401
        # Read positions + tags from the read-only state snapshot, compute the
        # restoring force analytically, then write into the per-row force
        # arrays. Both contexts use ROW indexing; HOOMD guarantees row i in
        # cpu_local_snapshot is the same particle as row i in
        # cpu_local_force_arrays on a single rank (same contract erm.py uses).
        with self._state.cpu_local_snapshot as snap:
            tag = np.asarray(snap.particles.tag).copy()
            pos = np.asarray(snap.particles.position).copy()

        mask = np.isin(tag, self._ligand_tags)
        F_vec = np.zeros_like(pos)
        U_per = np.zeros(pos.shape[0], dtype=np.float64)

        rows = np.flatnonzero(mask)
        if rows.size > 0:
            # Gather each masked row's anchor BY TAG (row order ≠ tag order).
            anchors = np.array(
                [self._anchor_by_tag[int(tag[r])] for r in rows],
                dtype=np.float64,
            )
            dx = pos[rows] - anchors            # r_i − a_i   [m]
            F_vec[rows] = -self.p.k_sub * dx     # restoring toward a_i
            U_per[rows] = 0.5 * self.p.k_sub * np.einsum("ij,ij->i", dx, dx)

        with self.cpu_local_force_arrays as arrays:
            # Replace the force / energy arrays entirely (HOOMD accumulates
            # these into net_force / net_energy for THIS force compute).
            arrays.force[:] = F_vec
            arrays.potential_energy[:] = U_per


# ---------------------------------------------------------------------------
# Public helper: attach the substrate-anchor spring to an existing simulation
# ---------------------------------------------------------------------------
def attach_substrate_spring(
    sim: hoomd.Simulation,
    k_sub: float,
    ligand_tags: np.ndarray,
    anchor_positions: np.ndarray,
    gamma_ligand: float,
    *,
    cfl_safety_factor: float = 0.1,
    cfl_strict: bool = True,
) -> SubstrateAnchorSpring:
    """Append a :class:`SubstrateAnchorSpring` to an existing MD Integrator.

    This is the **finite-compliance** substrate attach. Attaching it is a
    deliberate "compliant-substrate" selection: a caller wanting the
    bit-for-bit legacy rigid backing must NOT call this and must keep
    :class:`ffn_sim.archive.hoomd_legacy.cell.cell.SubstrateLigandPin` instead (the two are
    mutually exclusive on the same ligand tags — see module docstring). The
    Lead wires exactly one (rigid-pin XOR this spring) per substrate; when
    this is NOT attached, existing runs are bit-for-bit identical.

    CFL gate (parity with ``erm.py`` / ``enclosed_volume.py``): the spring's
    per-ligand stiffness is exactly ``k_sub``, so the relax time is
    ``τ = γ_ligand / k_sub`` and we require
    ``dt ≤ cfl_safety_factor · τ``; otherwise raise (a stiff ``k_sub``
    tightens CFL).

    Args:
        sim: Already-built simulation (must have an ``md.Integrator`` set).
        k_sub: Per-ligand anchor-spring stiffness [N/m] (REQUIRED; the ECM
            condition's substrate stiffness — no magic default).
        ligand_tags: Global HOOMD tags of the substrate-ligand particles.
        anchor_positions: Construction-time anchor sites, shape
            ``(len(ligand_tags), 3)`` [m], row-aligned with ``ligand_tags``.
        gamma_ligand: Per-ligand Stokes drag [N·s/m], used for the CFL gate
            ``dt ≤ cfl_safety_factor · γ_ligand / k_sub``.
        cfl_safety_factor: Same convention as the bond/angle/ERM/EV CFL gates
            (D3 BAOAB). Default ``0.1``.
        cfl_strict: If ``True`` (default) raise on CFL violation; set ``False``
            for diagnostic-only runs.

    Returns:
        The attached :class:`SubstrateAnchorSpring` (kept for introspection).

    Raises:
        RuntimeError: If ``sim`` has no integrator, or the CFL gate fails with
            ``cfl_strict=True``.
        ValueError: If ``k_sub`` or ``gamma_ligand`` is non-finite-positive
            (the force constructor / gate validates ``k_sub``; the gate
            validates ``gamma_ligand``).
    """
    ig = sim.operations.integrator
    if ig is None:
        raise RuntimeError(
            "sim.operations.integrator must be set before attaching the "
            "substrate-anchor spring."
        )
    _require_finite_positive("k_sub", k_sub)
    _require_finite_positive("gamma_ligand", gamma_ligand)

    # CFL gate — copied from enclosed_volume.py / erm.py. The spring's
    # per-ligand stiffness IS k_sub (not a derived softer k_eff), so a stiff
    # substrate genuinely tightens CFL.
    dt = float(ig.dt)
    tau = gamma_ligand / k_sub
    dt_cfl = cfl_safety_factor * tau
    if dt > dt_cfl and cfl_strict:
        raise RuntimeError(
            f"Substrate-spring CFL violated: dt = {dt:.3e} s > "
            f"{cfl_safety_factor:.2f} · τ = {dt_cfl:.3e} s "
            f"(τ = γ_ligand / k_sub = {tau:.3e} s, k_sub = {k_sub:.3e} N/m, "
            f"γ_ligand = {gamma_ligand:.3e} N·s/m). Either soften k_sub "
            "(less rigid substrate — changes the ECM condition; the rigid "
            "endpoint is the SubstrateLigandPin path, not a huge k_sub) OR "
            f"reduce dt (need dt ≤ {dt_cfl:.3e} s, "
            f"{(dt / dt_cfl):.1f}× smaller). Pass cfl_strict=False to skip "
            "this gate for diagnostic runs."
        )

    p = ResolvedSubstrate(k_sub=float(k_sub))
    spring = SubstrateAnchorSpring(p, ligand_tags, anchor_positions)
    ig.forces.append(spring)
    return spring
