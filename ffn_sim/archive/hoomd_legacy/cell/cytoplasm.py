"""H.10 cytoplasm — Tier-1 effective-viscosity drag contrast (KU-3.B3.1).

ADDITIVE + DEFAULT-OFF + FDT-SAFE cytoplasm medium model. The cytoplasm is the
**background medium** the cell moves through, not a force; Tier-1 captures the
single most robust metastatic discriminator in the dataset — the cytoplasm
**viscosity contrast** (MCF7 ≈ 5× MDA-MB-231) — by replacing the water viscosity
``η_water`` with a per-cell-type effective cytoplasm viscosity ``η_eff`` for the
cytoplasm-immersed bead types, via the **existing** BAOAB ``gamma_map`` contract.

This is a *value* change to the per-type Stokes drag ``γ_b = 6π η R``, fed to the
Leimkuhler-Matthews BAOAB integrator (``ffn_sim/integrator/baoab.py``). It does
**NOT** touch the integrator code, does **NOT** add particles, does **NOT** enter
the method-of-planes / virial (cytoplasm is dynamics-only — drift & diffusion —
orthogonal to the static cortex-γ tension sum). Tier-2 (Markovian-GLE storage
modulus) and Tier-3 (poroelastic filler particles) reach into the frozen
``integrator/`` and are **PI-GATED** — designed in ``docs/H10_CYTOPLASM_DESIGN.md``
but **NOT** implemented here.

Design reference
----------------
``ffn_sim/docs/H10_CYTOPLASM_DESIGN.md`` §"Tier 1" + ``docs/briefs/H10_cytoplasm.md``
+ ``docs/CANCER_CELLTYPE_PARAM_MAP.md`` §Cytoplasm-viscosity.

How it attaches (Lead integration — ONE call site)
---------------------------------------------------
``cell.py`` builds ``gamma_map`` with each cytoplasm-immersed type at the water
drag ``p_cortex.gamma_b``. After that dict is assembled and **before** it is
handed to ``make_baoab_updater`` / ``make_constrained_baoab_updater``, the Lead
applies::

    from ffn_sim.archive.hoomd_legacy.cell.cytoplasm import apply_cytoplasm_drag
    gamma_map = apply_cytoplasm_drag(gamma_map, p_cytoplasm)  # p_cytoplasm=None → no-op

When ``p_cytoplasm is None`` (the default) the dict is returned **unchanged**, so
pre-H.10 runs are bit-for-bit identical. When provided, the immersed types' drag
is scaled by ``η_eff / η_water`` (the override is expressed as a *ratio* of the
already-resolved ``γ_b`` so that ``η_eff == η_water`` is exactly the input float —
no radius round-trip, no rounding drift). Membrane / cortex-shell (water-exposed)
types are left at water drag unless explicitly listed as immersed.

FDT-safety (load-bearing)
-------------------------
BAOAB sets each type's thermal-noise amplitude
``bd_prefactor = √(kT / (2 γ Δt))`` from the **same** per-type ``γ`` that drives
the deterministic drift ``F/γ·Δt`` (``baoab.py`` lines ~315 & ~462-469). Therefore
the equilibrium diffusion ``D = kT/γ`` and the velocity variance ``⟨v²⟩ = kT/m``
stay **exact** for any per-type ``γ`` — raising ``η_eff`` lowers ``D`` by exactly
``η_water/η_eff`` (slower diffusion, the correct sign) with **no** fluctuation-
dissipation violation. This is why Tier-1 is the cheapest faithful option: FDT is
preserved by the integrator's own per-type construction, not by anything here.

Magic-Number Block (η bands — Dessard et al. 2024, provenance)
-------------------------------------------------------------
The cytoplasm viscosities are **measured acceptance oracles**, never fitting
targets, never grid-tuned. Microrheology viscosity is probe-length-dependent;
for MCF-7 the source reports η = 65.9 ± 11.4 Pa·s as the ALL-WIRE ensemble mean
(n=60) and η = 56.4 ± 16.6 Pa·s for the L = 3 ± 1 µm wire-length subset.

  - MCF10A (normal)            η = 41.6 Pa·s   (MRS;  KU-3.B3.1 / H.10 brief)
  - MCF-7   (low-invasion)     η = 65.9 ± 11.4 Pa·s  (all-wire mean; 56.4 L=3µm)
  - MDA-MB-231 (high-invasion) η = 12.0 Pa·s (dict)  ⚠ Dessard source = 10.7 ± 5.4;
                               the "12.0 refined" provenance is unclear → flag PI
  - contrast: MCF-7 ≈ 5× MDA — the validation target.

Source: Dessard, Manneville & Berret 2024, *Nanoscale Adv* 6(6):1727-1738
(DOI 10.1039/d4na00003j, PMC10929591), magnetic rotational spectroscopy on
adherent cells at 37 °C; cross-checked against H.10 brief KU-3.B3.1. (Citation
fix 2026-06-12: the value+PMC are correct; the prior "Hu 2024" attribution was
wrong — first author Marie Dessard, corresp. J.-F. Berret. See
docs/v2_audit/MCF7_PARAMETER_COLLECTION_2026-06-12.md §6.) ``η_water = 6.913e-4
Pa·s`` (NIST, 310 K; KU-1.26) is the default → bit-identity. Defaults live in
``ffn_sim/cell/cytoplasm.py`` constants below and the additive ``cytoplasm:``
block in ``ffn_sim/configs/phase1_h10.yaml``.

PI-gate (why this is DEFAULT-OFF, not auto-on)
----------------------------------------------
``η_eff`` is ~5×10⁴× the water drag. In the overdamped BAOAB limit a larger drag
*intentionally* slows the cell (it is deeply overdamped in cytoplasm) and
*relaxes* every τ = γ_b/k_eff CFL gate (see :func:`cfl_relaxation_at_drag_jump`),
but it is a large, physics-altering, frozen-medium-policy change that interacts
with gates tuned at η_water. Per ``docs/H10_CYTOPLASM_DESIGN.md`` §"The PI-gate":
**the η_eff magnitude + a re-run of the existing CFL/timestep gates must be
surfaced to PI before enabling.** Not a magic number (η_eff is the Dessard-2024
measurement); a medium-policy gate. Tier-2/3 are separate frozen-``integrator/``
PI-gates and are not built here.

Sanity Gate
-----------
*Per CLAUDE.md Hard Rule "Sanity Gate Protocol mandatory before first execution
of any physics/numerics module." STATIC/analytical checks live in
``ffn_sim/tests/test_h10_cytoplasm.py``; the FDT/diffusion check is the Brownian
smoke there.*

1. **Dimensional.** ``γ_b = 6π η R`` → [Pa·s · m] = [N·s/m] ✓.
   ``D = kT/γ_b`` → [J] / [N·s/m] = [N·m·m/(N·s)] = [m²/s] ✓. The override is a
   dimensionless ratio ``η_eff/η_water`` × an existing [N·s/m] → [N·s/m] ✓.
2. **Boundary.** ``p_cytoplasm is None`` **or** ``η_eff == η_water`` returns the
   gamma_map **unchanged** (exact float equality) → the run is bit-for-bit
   identical (the default-off contract). Tested with ``==``, not ``isclose``.
3. **Conservation / FDT (critical).** BAOAB's noise amplitude reads γ **per-type**
   from the same gamma_map (verified against ``baoab.py``), so ``D = kT/γ_b`` and
   ``⟨v²⟩ = kT/m`` stay exact under per-type ``η_eff``. The Brownian smoke asserts
   the measured ``D`` matches ``kT/γ_b`` per type, hence the MCF7/MDA D-ratio.
4. **Numerical / CFL.** ``τ_v = m/γ_b`` shrinks ∝ 1/η_eff (BAOAB overdamped limit
   deepens — fine); the *stiffness* CFL gates are τ = γ_b/k_eff and **grow** ∝
   η_eff, so a larger drag **relaxes** them. :func:`cfl_relaxation_at_drag_jump`
   re-asserts (does not assume) that the dt-ceiling at the raised drag is ≥ the
   ceiling at water drag, at the worst-case ~5×10⁴× jump.
5. **Sign-sense.** Higher ``η_eff`` ⇒ lower ``D`` ⇒ slower diffusion. MCF7 (more
   viscous) diffuses **slower** than MDA. Wrong sign = bug.
6. **Measurement-protocol consistency.** η reported at the stated probe length
   L ≈ 3 µm (η ∝ L²); the contrast (ratio), not the absolute Pa·s, is the gate.

References
----------
- Dessard, Manneville & Berret 2024, *Nanoscale Adv* 6(6):1727-1738 (10.1039/d4na00003j, PMC10929591) — cytoplasm viscosity by MRS.
- ``ffn_sim/docs/H10_CYTOPLASM_DESIGN.md`` (Tier-1/2/3 staging, PI-gates).
- ``ffn_sim/docs/briefs/H10_cytoplasm.md`` (KU-3.B3 anchors).
- ``ffn_sim/integrator/baoab.py`` (the gamma_map contract + per-type FDT noise).
"""

from __future__ import annotations

import logging
import math
from dataclasses import dataclass, field
from typing import Mapping

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Literature-anchored constants (Magic-Number Block — see module docstring)
# ---------------------------------------------------------------------------
#: Water viscosity at 310 K [Pa·s] (NIST; KU-1.26). The DEFAULT η_eff →
#: bit-for-bit identity with the pre-H.10 gamma_map.
ETA_WATER: float = 6.913e-4

#: Per-cell-type cytoplasm effective viscosities [Pa·s], MRS probe L ≈ 3 µm.
#: Dessard 2024 (PMC10929591) for MCF7/MDA; MRS for MCF10A. Acceptance
#: oracles, NEVER fitting targets. ``"water"`` is the explicit no-op default.
ETA_CYTO_BY_CELLTYPE: dict[str, float] = {
    "water": ETA_WATER,   # default / no-op (bit-identity)
    "MCF10A": 41.6,       # normal epithelial  (MRS; KU-3.B3.1)
    "MCF7": 65.9,         # low-invasion       (Dessard 2024 all-wire; 56.4 L=3µm)
    "MDA-MB-231": 12.0,   # high-invasion      (dict; Dessard 2024 = 10.7 — flag PI)
}

#: Cytoplasm-immersed interior particle types (water-exposed membrane / cortex
#: shell types are NOT scaled by default). These are the interior actin /
#: motor / crosslinker / lamellipodium beads bathed in the cytosol. The cortex
#: actin shell itself sits at the cytoplasm/membrane interface; whether it is
#: "immersed" is an open PI item (H10 design §Open items), so it is included by
#: default as interior (conservative for the viscosity-contrast discriminator)
#: and can be overridden per-cell via ``CytoplasmTier1.immersed_types``.
DEFAULT_IMMERSED_TYPES: tuple[str, ...] = (
    "actin_cortex",
    "xlink_head",
    "cortex_myosin_backbone",
    "cortex_myosin_head",
    "actin_lamel",
    "wave_particle",
    "nucleus_bead",
)


@dataclass(slots=True, frozen=True)
class CytoplasmTier1:
    """Resolved Tier-1 cytoplasm: a per-type ``η_eff`` drag-override spec.

    DEFAULT-OFF: ``eta_eff == ETA_WATER`` (the default) is a no-op — the
    gamma_map override is the identity (exact float equality), so the run is
    bit-for-bit identical to pre-H.10. Construct via :func:`resolve_cytoplasm`
    (config / cell-type) and apply via :func:`apply_cytoplasm_drag`.

    Attributes
    ----------
    eta_eff
        Cytoplasm effective viscosity [Pa·s] for the immersed types. Default
        ``ETA_WATER`` → no-op. A literature value (Dessard 2024) when enabled.
    eta_water
        Reference water viscosity [Pa·s] that the resolved ``gamma_b`` already
        encodes. The override scales each immersed type's drag by
        ``eta_eff / eta_water``, so ``eta_eff == eta_water`` is the identity.
    immersed_types
        Particle-type names whose Stokes drag is scaled to cytoplasm viscosity.
        Types absent from this set keep their water drag. Defaults to
        :data:`DEFAULT_IMMERSED_TYPES`.
    cell_type
        Optional provenance label (``"MCF7"`` …) for logging / reporting.
    """

    eta_eff: float = ETA_WATER
    eta_water: float = ETA_WATER
    immersed_types: frozenset[str] = field(
        default_factory=lambda: frozenset(DEFAULT_IMMERSED_TYPES)
    )
    cell_type: str | None = None

    def __post_init__(self) -> None:
        for name, v in (("eta_eff", self.eta_eff), ("eta_water", self.eta_water)):
            if not (math.isfinite(v) and v > 0.0):
                raise ValueError(f"{name} must be finite and > 0, got {v!r}")

    @property
    def viscosity_ratio(self) -> float:
        """``η_eff / η_water`` — the dimensionless drag-scale factor (= 1 ⇒ no-op)."""
        return self.eta_eff / self.eta_water

    @property
    def is_noop(self) -> bool:
        """True iff the override is the exact identity (``η_eff == η_water``)."""
        return self.eta_eff == self.eta_water

    def diffusion_ratio_vs_water(self) -> float:
        """``D_eff / D_water = η_water / η_eff`` for an immersed bead (sign check).

        Higher viscosity ⇒ lower diffusion (the metastatic discriminator's
        correct sign). For MCF7 vs MDA this returns the ~1/5.5 slowdown.
        """
        return self.eta_water / self.eta_eff


def resolve_cytoplasm(
    *,
    eta_eff: float | None = None,
    cell_type: str | None = None,
    eta_water: float = ETA_WATER,
    immersed_types: Mapping[str, object] | tuple[str, ...] | None = None,
) -> CytoplasmTier1:
    """Resolve a :class:`CytoplasmTier1` from an explicit ``η_eff`` or a cell-type.

    Exactly one viscosity source is used, in priority order:
      1. explicit ``eta_eff`` (a measured Pa·s value), else
      2. ``cell_type`` looked up in :data:`ETA_CYTO_BY_CELLTYPE`
         (``MCF10A`` / ``MCF7`` / ``MDA-MB-231`` / ``water``), else
      3. ``ETA_WATER`` (the no-op default).

    Args:
        eta_eff: Cytoplasm viscosity [Pa·s]; overrides ``cell_type`` if given.
        cell_type: Cell-type key into the Dessard-2024 viscosity table.
        eta_water: Reference water viscosity the resolved ``γ_b`` encodes
            [Pa·s] (default NIST 310 K). The override is a ratio to this.
        immersed_types: Iterable of immersed type names; default
            :data:`DEFAULT_IMMERSED_TYPES`.

    Returns:
        A :class:`CytoplasmTier1`. With no viscosity source it is a no-op
        (``η_eff == η_water``), preserving bit-identity.

    Raises:
        KeyError: ``cell_type`` not in the viscosity table.
        ValueError: non-finite / non-positive viscosity.
    """
    if eta_eff is not None:
        eta = float(eta_eff)
    elif cell_type is not None:
        if cell_type not in ETA_CYTO_BY_CELLTYPE:
            raise KeyError(
                f"cell_type {cell_type!r} not in ETA_CYTO_BY_CELLTYPE "
                f"{sorted(ETA_CYTO_BY_CELLTYPE)}; supply an explicit eta_eff "
                "or add a literature-anchored entry (Dessard 2024)."
            )
        eta = float(ETA_CYTO_BY_CELLTYPE[cell_type])
    else:
        eta = float(eta_water)

    if immersed_types is None:
        types = frozenset(DEFAULT_IMMERSED_TYPES)
    else:
        types = frozenset(immersed_types)

    return CytoplasmTier1(
        eta_eff=eta,
        eta_water=float(eta_water),
        immersed_types=types,
        cell_type=cell_type,
    )


def cytoplasm_gamma_override(
    gamma_map: Mapping[str, float],
    cyto: CytoplasmTier1,
) -> dict[str, float]:
    """Return a NEW gamma_map with immersed types scaled to cytoplasm viscosity.

    Pure & non-mutating: builds a fresh dict. Each present type in
    ``cyto.immersed_types`` has its drag multiplied by ``η_eff / η_water``;
    every other entry is copied verbatim. Because the scale is a ratio of the
    **already-resolved** water ``γ_b``, ``η_eff == η_water`` reproduces the
    input float **exactly** (no ``6πηR`` round-trip), so the default is a
    bit-for-bit identity (Sanity Gate §2). Immersed types absent from
    ``gamma_map`` are ignored (a subset of cell features may be enabled).

    Args:
        gamma_map: The water-drag per-type Stokes drag dict (BAOAB contract).
        cyto: Resolved Tier-1 spec (:func:`resolve_cytoplasm`).

    Returns:
        A new ``dict[str, float]`` — the override to hand to BAOAB.
    """
    ratio = cyto.viscosity_ratio
    out: dict[str, float] = dict(gamma_map)
    if cyto.is_noop:
        # Exact identity: return the unchanged copy (no float arithmetic at all,
        # so values are bit-for-bit the input — the default-off contract).
        return out
    scaled = []
    for typ in gamma_map:
        if typ in cyto.immersed_types:
            out[typ] = float(gamma_map[typ]) * ratio
            scaled.append(typ)
    logger.info(
        "Cytoplasm Tier-1 ENABLED (cell_type=%s): η_eff=%.4g Pa·s, "
        "ratio η_eff/η_water=%.4g× applied to %d immersed type(s) %s; "
        "D_eff/D_water=%.4g (slower diffusion). PI-gate: re-run CFL/timestep "
        "gates at this drag.",
        cyto.cell_type, cyto.eta_eff, ratio, len(scaled), sorted(scaled),
        cyto.diffusion_ratio_vs_water(),
    )
    return out


def apply_cytoplasm_drag(
    gamma_map: Mapping[str, float],
    p_cytoplasm: CytoplasmTier1 | None,
) -> dict[str, float]:
    """Lead-facing single call site: apply the Tier-1 override, or no-op.

    DEFAULT-OFF: ``p_cytoplasm is None`` returns ``dict(gamma_map)`` unchanged
    (bit-identity). Otherwise delegates to :func:`cytoplasm_gamma_override`.
    This is the **one** line ``cell.py`` adds, immediately before building the
    BAOAB updater::

        gamma_map = apply_cytoplasm_drag(gamma_map, p_cytoplasm)

    Args:
        gamma_map: Per-type water Stokes-drag dict (BAOAB contract).
        p_cytoplasm: Resolved Tier-1 spec, or ``None`` to disable.

    Returns:
        A new ``dict[str, float]`` (always a copy, never the input object).
    """
    if p_cytoplasm is None:
        return dict(gamma_map)
    return cytoplasm_gamma_override(gamma_map, p_cytoplasm)


# ---------------------------------------------------------------------------
# CFL re-check at the drag jump (Sanity Gate §4)
# ---------------------------------------------------------------------------
@dataclass(slots=True, frozen=True)
class CFLReCheck:
    """Result of the τ = γ_b/k_eff CFL re-assertion at the drag jump.

    Attributes
    ----------
    ratio
        Drag scale ``η_eff/η_water`` applied.
    dt_ceiling_water
        ``cfl_safety · γ_water / k_eff`` — the stiffness-CFL dt ceiling at water
        drag [s] (the gate the run is already tuned to).
    dt_ceiling_cyto
        ``cfl_safety · γ_cyto / k_eff`` — the dt ceiling at the cytoplasm drag.
    relaxes
        True iff ``dt_ceiling_cyto >= dt_ceiling_water`` (the larger drag is
        SAFE — it relaxes, not tightens, the stiffness CFL).
    tau_v_water, tau_v_cyto
        Velocity-relaxation times ``m/γ`` [s] — these SHRINK with drag (the
        overdamped BAOAB limit only deepens; informational, not a gate).
    """

    ratio: float
    dt_ceiling_water: float
    dt_ceiling_cyto: float
    relaxes: bool
    tau_v_water: float
    tau_v_cyto: float


def cfl_relaxation_at_drag_jump(
    *,
    gamma_water: float,
    k_eff: float,
    cfl_safety: float,
    eta_water: float = ETA_WATER,
    eta_eff: float,
    bead_mass: float = 1.0,
) -> CFLReCheck:
    """Re-assert (do not assume) the stiffness CFL at the drag jump.

    The BAOAB stiffness CFL is a relaxation-time gate ``dt ≤ cfl_safety · τ``
    with ``τ = γ_b / k_eff``. Scaling drag by ``η_eff/η_water`` scales ``τ`` by
    the same factor, so the dt **ceiling grows** with viscosity — a larger drag
    cannot violate a gate tuned at water drag. This function computes both
    ceilings and the boolean ``relaxes`` at the actual (worst-case ~5×10⁴×)
    jump so it can be asserted in the gate test rather than assumed.

    The velocity-relaxation time ``τ_v = m/γ`` SHRINKS ∝ 1/η_eff (the system
    sits even deeper in the overdamped limit BAOAB targets) — reported but not
    a CFL constraint for the overdamped integrator.

    Args:
        gamma_water: Water-drag Stokes γ_b [N·s/m] (e.g. ``p_cortex.gamma_b``).
        k_eff: Stiffest effective spring the CFL gate guards [N/m].
        cfl_safety: CFL safety factor (e.g. ``p_cortex.cfl_safety_factor``).
        eta_water: Reference water viscosity [Pa·s].
        eta_eff: Cytoplasm viscosity [Pa·s].
        bead_mass: Bead mass [kg] for the informational τ_v (default 1.0).

    Returns:
        A :class:`CFLReCheck`.

    Raises:
        ValueError: non-finite / non-positive inputs.
    """
    for name, v in (
        ("gamma_water", gamma_water), ("k_eff", k_eff),
        ("cfl_safety", cfl_safety), ("eta_water", eta_water),
        ("eta_eff", eta_eff), ("bead_mass", bead_mass),
    ):
        if not (math.isfinite(v) and v > 0.0):
            raise ValueError(f"{name} must be finite and > 0, got {v!r}")

    ratio = eta_eff / eta_water
    gamma_cyto = gamma_water * ratio
    dt_water = cfl_safety * gamma_water / k_eff
    dt_cyto = cfl_safety * gamma_cyto / k_eff
    return CFLReCheck(
        ratio=ratio,
        dt_ceiling_water=dt_water,
        dt_ceiling_cyto=dt_cyto,
        relaxes=(dt_cyto >= dt_water),
        tau_v_water=bead_mass / gamma_water,
        tau_v_cyto=bead_mass / gamma_cyto,
    )


def stokes_drag(eta: float, radius: float) -> float:
    """``γ_b = 6π η R`` [N·s/m] — Stokes drag (diagnostic / reporting helper).

    Mirrors ``ffn_sim/cell/dt_reconcile.py::_stokes_drag`` and the cortex
    ``γ_b`` derivation. The runtime override uses the ratio form (so it is
    radius-agnostic and bit-identical at η=water); this helper exists only for
    reporting the absolute γ at a stated radius.
    """
    if not (math.isfinite(eta) and eta > 0.0):
        raise ValueError(f"eta must be finite and > 0, got {eta!r}")
    if not (math.isfinite(radius) and radius > 0.0):
        raise ValueError(f"radius must be finite and > 0, got {radius!r}")
    return 6.0 * math.pi * eta * radius
