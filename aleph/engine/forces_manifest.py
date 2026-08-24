"""Force-channel OBSERVER — the on/off state + numeric value of every force channel and flag-gated path.

WHY THIS EXISTS (designed 2026-06-21, never built; ranked #1 on the 2026-07-25 audit's
convert-to-structure list).

The most expensive mistake in project history: a lumped ``SettlingForce`` proxy running at
19,000-47,000x real buoyant weight was ON in every committed spreading run while the artifact
recorded it as ``proxy_off``.  It performed the flattening that was credited to the lamellipodium,
and the project's own regression gates locked the confound in — so a CORRECT fix would have
registered as a regression.  It invalidated every pre-06-21 DCM spreading run including the
capstone A/A0 law, and it was caught only by PI pushback.  The same class recurred four weeks
later in a different engine (Biot FSI ``default-off`` in every landed native run).  The 06-20
closeout credits a permanent ``forces_manifest`` recurrence gate; the 2026-07-25 audit checked and
found that gate exists ONLY in that document and its Obsidian mirror — never in code.  This module
is that gate's missing observer half.

WHAT IT IS AND IS NOT
  * It is an OBSERVER.  It changes no physics, no default, and no value.  It imports the physics
    modules read-only and reports what they hold.
  * It is NOT a source of truth for values.  Every magnitude is READ from the live module
    constants / the live config / the live composed object.  Nothing is hardcoded here — a
    hardcoded value table would let the observer and the physics drift apart while continuing to
    agree with each other, which is the "gate a broken build cannot fail" antipattern.
  * It does NOT infer.  A channel that cannot be observed is reported ``enabled=None`` with
    ``observation="unobserved"``.  Absence is NEVER rendered as "off" — reading a missing record
    as "off" is precisely the mislabel that produced ``proxy_off``.

THE TWO HISTORICAL SHAPES IT MUST EXPOSE
  (i)  a channel declared off that is actually on  -> ``config_intent`` vs ``enabled`` disagree
       (recorded per channel as ``intent_mismatch``), plus env-var overrides that flip a
       default-off numeric path ON without appearing in any config record.
  (ii) a channel declared on at a physiological value that is actually running at a
       convenience/default value  -> ``value`` + ``provenance`` are always reported together, so a
       CONVENIENCE or PI_GAP magnitude on an enabled production channel is visible.

OBSERVATION STRENGTH.  Each channel carries ``observation``:
  ``runtime``     read off the live composed :class:`~aleph.components.incumbent.assemble.AssembledCell` /
                  solver object — the strongest evidence (what the run actually did).
  ``config``      read off the run's ``CellConfig`` record (intent; what was asked for).
  ``code-default`` read off a live module constant because the config does not carry it.
  ``env``         read off the process environment (the driver-route override channels).
  ``unobserved``  not determinable from the inputs supplied.  Never coerced to a boolean.

Verification against the checked-in production declaration is a SEPARATE tool
(``aleph/outputs/tag_kb/verify_forces.py``); this module only observes.

Sanity Gate
-----------
* Dimensional: every numeric channel carries an explicit ``unit`` string in engine units
  (um-pN-s) or the literal ``"bool"``/``"count"``/``"enum"``; no bare numbers.
* Boundary: a channel absent from a config record yields ``enabled=None`` /
  ``observation="unobserved"``, never ``False`` (the ``proxy_off`` failure mode).
* Conservation/invariant: the FORCE-kind channel set is asserted against the actual
  ``ac.cell.driver._accumulate_all`` assembly contract by
  :func:`assert_force_channels_cover_accumulate`, so a force channel added to the assembly
  without being declared here is a hard error rather than a silent omission.
* Sign-sense: not applicable — no force is computed here (read-only observer).
* Measurement-protocol: ``observation`` records HOW each state was determined, so a
  config-only dump can never be reported as runtime evidence.
"""

from __future__ import annotations

import ast
import json
import os
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Callable

__all__ = [
    "SCHEMA",
    "PROVENANCE",
    "KINDS",
    "OBSERVATIONS",
    "ForceChannel",
    "dump_force_channels",
    "dump_from_artifact",
    "write_dump",
    "assert_force_channels_cover_accumulate",
    "ENV_OVERRIDE_VARS",
]

SCHEMA = "ffn-ac-forces-manifest-v1"

#: Provenance vocabulary.  SOURCED = traceable to a literature/KB anchor; DERIVED = computed from
#: a sourced anchor or from geometry (Magic-Number-Block clean); CONVENIENCE = a numerical guard or
#: backward-compatibility default that is explicitly NOT a physiological magnitude; PI_GAP = the
#: physiological value is unknown and the code says so (``None``/TEST/provisional).
# Re-exported. The definition moved to `aleph/units/provenance.py` on 2026-08-20 so that
# `world/` — canonical since that date — can read it without importing the tree it replaces.
from aleph.units.provenance import PROVENANCE  # noqa: E402,F401

#: FORCE      — a term summed into the node force array by ``_accumulate_all``.
#: SUBCHANNEL — a force term inside a composed compartment's own ``accumulate``.
#: KINETIC    — an event/rate path that commits only on an accepted step.
#: BUILD      — a build-time choice that changes the initial state the physics runs from.
#: NUMERIC    — a solver/preconditioner/accelerator path that changes a convergence claim.
KINDS = ("FORCE", "SUBCHANNEL", "KINETIC", "BUILD", "NUMERIC")

OBSERVATIONS = ("runtime", "config", "code-default", "env", "unobserved")

#: Environment variables that flip an otherwise default-off numeric path ON.  These are read by
#: ``ac.cell.implicit_mechanics.ProjectedAnalyticCG.__init__`` and are NOT recorded in the driver's
#: report JSON — so a run with ``AC_MG=1`` is indistinguishable from one without, after the fact.
#: That is a live instance of shape (i); see the module docstring.
ENV_OVERRIDE_VARS = (
    "AC_FQ_COARSE", "AC_FQ_MODE", "AC_FQ_ITERS", "AC_MG", "AC_MG_FQ", "AC_MG_ASM",
)

_MISSING = object()


@dataclass
class ForceChannel:
    """One observed force channel or flag-gated path.

    Attributes:
        name: Stable channel id (the key the declaration is written against).
        kind: One of :data:`KINDS`.
        enabled: ``True``/``False`` when observable, ``None`` when it is not.  ``None`` is a
            first-class value and must never be collapsed to ``False``.
        value: The numeric/enum magnitude actually in force, or ``None`` for a pure on/off path
            or an unsourced (PI-GAP) magnitude the code leaves ``None``.
        unit: Explicit unit in engine units (um-pN-s), or ``bool``/``count``/``enum``/``str``.
        provenance: One of :data:`PROVENANCE`.
        code_site: ``path::symbol`` of the site that owns the value (symbols, not line numbers,
            so the reference does not rot on edits).
        gate: The expression that decides ``enabled``, quoted from the code.
        observation: One of :data:`OBSERVATIONS` — HOW ``enabled``/``value`` were determined.
        config_intent: What the config record asked for, when that is separately knowable.
        intent_mismatch: ``True`` when a runtime observation contradicts ``config_intent``.
        notes: Free text carried from the code comment that documents the channel.
    """

    name: str
    kind: str
    enabled: bool | None
    value: Any
    unit: str
    provenance: str
    code_site: str
    gate: str
    observation: str
    config_intent: Any = None
    intent_mismatch: bool = False
    notes: str = ""

    def __post_init__(self) -> None:
        if self.kind not in KINDS:
            raise ValueError(f"{self.name}: kind {self.kind!r} not in {KINDS}")
        if self.provenance not in PROVENANCE:
            raise ValueError(f"{self.name}: provenance {self.provenance!r} not in {PROVENANCE}")
        if self.observation not in OBSERVATIONS:
            raise ValueError(f"{self.name}: observation {self.observation!r} not in {OBSERVATIONS}")
        if self.observation == "unobserved" and self.enabled is not None:
            raise ValueError(f"{self.name}: unobserved channel must not claim enabled={self.enabled!r}")


# ─────────────────────────────────────────────────────────────────────────────────────────────────
# Config access that distinguishes ABSENT from None.  An artifact whose config record predates a
# channel does not record "off" for it — it records nothing, and that difference is load-bearing.
# ─────────────────────────────────────────────────────────────────────────────────────────────────
def _cfg_get(cfg: Any, key: str) -> Any:
    """Return ``cfg[key]`` / ``cfg.key``, or :data:`_MISSING` when the record has no such field."""
    if cfg is None:
        return _MISSING
    if isinstance(cfg, dict):
        return cfg.get(key, _MISSING)
    return getattr(cfg, key, _MISSING)


def _bool_channel(
    cfg: Any, key: str, *, name: str, kind: str, provenance: str, code_site: str,
    notes: str = "", runtime: bool | None = None, runtime_gate: str = "",
) -> ForceChannel:
    """Observe a boolean flag from the config, optionally cross-checked against runtime truth."""
    raw = _cfg_get(cfg, key)
    intent = None if raw is _MISSING else bool(raw)
    if runtime is not None:
        return ForceChannel(
            name=name, kind=kind, enabled=bool(runtime), value=bool(runtime), unit="bool",
            provenance=provenance, code_site=code_site, gate=runtime_gate or f"cfg.{key}",
            observation="runtime", config_intent=intent,
            intent_mismatch=(intent is not None and bool(runtime) != intent), notes=notes,
        )
    if raw is _MISSING:
        return ForceChannel(
            name=name, kind=kind, enabled=None, value=None, unit="bool", provenance=provenance,
            code_site=code_site, gate=f"cfg.{key}", observation="unobserved",
            notes=(notes + " | NOT RECORDED in this config: state cannot be certified.").strip(" |"),
        )
    return ForceChannel(
        name=name, kind=kind, enabled=bool(raw), value=bool(raw), unit="bool",
        provenance=provenance, code_site=code_site, gate=f"cfg.{key}", observation="config",
        config_intent=intent, notes=notes,
    )


def _value_channel(
    *, name: str, kind: str, enabled: bool | None, value: Any, unit: str, provenance: str,
    code_site: str, gate: str, observation: str, config_intent: Any = None, notes: str = "",
) -> ForceChannel:
    return ForceChannel(
        name=name, kind=kind, enabled=enabled, value=value, unit=unit, provenance=provenance,
        code_site=code_site, gate=gate, observation=observation, config_intent=config_intent,
        notes=notes,
    )


# ─────────────────────────────────────────────────────────────────────────────────────────────────
# The FORCE channels — one per term in ``ac.cell.driver._accumulate_all`` (the inner force-assembly
# contract).  ``gate`` quotes the actual guard; ``runtime_attr`` is the composed object whose
# presence IS the truth of whether the term is summed.
# ─────────────────────────────────────────────────────────────────────────────────────────────────
_ACCUM_SITE = "aleph/components/incumbent/driver.py::_accumulate_all"

#: name -> (config key or None, runtime predicate on AssembledCell, gate text, provenance, notes)
_FORCE_TERMS: tuple[tuple[str, str | None, Callable[[Any], bool], str, str, str], ...] = (
    (
        "actin_bending_cytosim", None, lambda c: int(getattr(c, "n_tri", 0)) > 0,
        "if cell.n_tri: cytosim_bending_kernel", "DERIVED",
        "Cytosim bending, alpha = kappa/seg^3 per triple; always on when the cortex has triples.",
    ),
    (
        "actin_crosslink_link_spring", None, lambda c: int(getattr(c, "n_xl", 0)) > 0,
        "if cell.n_xl: link_spring_kernel", "SOURCED",
        "Hookean crosslink spring; k_xl from Ferrer 2008 alpha-actinin (per-link array kxl_d).",
    ),
    (
        "arp23_branch_angle", None,
        lambda c: getattr(c, "branch_triples_d", None) is not None and int(getattr(c, "n_branch", 0)) > 0,
        "if cell.branch_triples_d is not None and cell.n_branch: branch_angle_kernel", "SOURCED",
        "Arp2/3 70deg angle-harmonic (Faessler 2020 theta0/k_theta). Fires only for the mixed "
        "formin+Arp2/3 cortex (cortex_arp23_fraction > 0); the formin-only default launches nothing.",
    ),
    (
        "nmii_minifilament_force", "with_myosin", lambda c: getattr(c, "myosin", None) is not None,
        "if cell.myosin is not None: cell.myosin.accumulate", "PI_GAP",
        "Head-resolved Stam-Hocky NMII (I3). Its per-head magnitudes (F_stall_head, k_xb, v0, "
        "kappa_hill, N_side) are all params_i0b3.yaml PI GAPs running at TEST values.",
    ),
    (
        "nucleus_compartment_force", "with_nucleus", lambda c: getattr(c, "nucleus", None) is not None,
        "if cell.nucleus is not None: cell.nucleus.accumulate", "PI_GAP",
        "Deformable-mesh nucleus: envelope bending + lamina (chromatin/lamin-B/lamin-A-C) + "
        "nucleoplasm volume + LINC. EVERY lamina modulus is an I0-B2 PI GAP at a TEST value.",
    ),
    (
        "membrane_compartment_force", "with_membrane", lambda c: getattr(c, "membrane", None) is not None,
        "if cell.membrane is not None: cell.membrane.accumulate", "SOURCED",
        "Plasma-membrane Helfrich sheet: bending + area tension + ERM tethers.",
    ),
    (
        "membrane_pressure_traction", None,
        lambda c: getattr(c, "membrane_pressure", None) is not None,
        "if cell.membrane_pressure is not None: accumulate", "DERIVED",
        "Live spatial pore pressure -> membrane normal traction (NG-3). Composed with the "
        "membrane + pressure pair; not independently flagged.",
    ),
    (
        "steric_wca", "with_steric", lambda c: getattr(c, "steric", None) is not None,
        "if cell.steric is not None: cell.steric.accumulate", "SOURCED",
        "All-fiber WCA excluded volume. CLAUDE.md hard rule: LJ repulsive ON from Phase 1.",
    ),
    (
        "biot_pressure_coupling", "with_pressure", lambda c: getattr(c, "pressure", None) is not None,
        "if cell.pressure is not None: cell.pressure.accumulate", "DERIVED",
        "Poroelastic fluid->solid coupling -alpha*p*I (I1a). Turning this off is the Biot-FSI "
        "default-off recurrence: the cell then runs with no pore-pressure load at all.",
    ),
)


def assert_force_channels_cover_accumulate() -> list[str]:
    """Assert the FORCE channel set matches the real ``_accumulate_all`` assembly contract.

    Parses ``ac/cell/driver.py`` statically (no import, no CUDA) and collects every ``wp.launch``
    kernel name and every ``<attr>.accumulate(...)`` call inside ``_accumulate_all``.  Raises
    ``AssertionError`` when the assembly contains a force term this module does not declare.

    This is the invariant that stops the observer from rotting: adding a new force term to the
    inner assembly without declaring it here fails loudly instead of going unobserved.

    Returns:
        The list of assembly terms found, for reporting.

    Raises:
        AssertionError: when a term in the assembly has no declared channel.
    """
    # `ac/cell/driver.py` moved to `aleph/components/incumbent/driver.py` in the 1b895ec0 rename and
    # this path did not follow it, so this function has raised FileNotFoundError on every call since —
    # i.e. the invariant this docstring calls "what stops the observer from rotting" has been DEAD, and
    # silently, because nothing that called it could report why. It also blocked every `observe`
    # artifact write. Restored, and it PASSES: 9 assembly terms, all declared.
    src = Path(__file__).resolve().parents[1] / "components" / "incumbent" / "driver.py"
    tree = ast.parse(src.read_text())
    fn = next(
        (n for n in ast.walk(tree)
         if isinstance(n, ast.FunctionDef) and n.name == "_accumulate_all"), None)
    if fn is None:  # pragma: no cover - the assembly function is the contract; its absence is fatal
        raise AssertionError("ac/cell/driver.py::_accumulate_all not found — assembly contract moved")
    terms: list[str] = []
    for node in ast.walk(fn):
        if not isinstance(node, ast.Call):
            continue
        f = node.func
        if isinstance(f, ast.Attribute) and f.attr == "launch":
            if node.args and isinstance(node.args[0], ast.Name):
                terms.append(node.args[0].id)
        elif isinstance(f, ast.Attribute) and f.attr == "accumulate":
            owner = f.value
            if isinstance(owner, ast.Attribute):
                terms.append(f"{owner.attr}.accumulate")
    # ``_zero`` clears the array; it contributes no force.
    terms = [t for t in terms if t != "_zero"]
    declared = {
        "cytosim_bending_kernel", "link_spring_kernel", "branch_angle_kernel",
        "myosin.accumulate", "nucleus.accumulate", "membrane.accumulate",
        "membrane_pressure.accumulate", "steric.accumulate", "pressure.accumulate",
    }
    undeclared = sorted(set(terms) - declared)
    if undeclared:
        raise AssertionError(
            "force terms in _accumulate_all with no declared channel in forces_manifest.py: "
            f"{undeclared} — declare them (and their production state) before running production.")
    return terms


# ─────────────────────────────────────────────────────────────────────────────────────────────────
# Observers
# ─────────────────────────────────────────────────────────────────────────────────────────────────
def _force_channels(cfg: Any, cell: Any) -> list[ForceChannel]:
    out: list[ForceChannel] = []
    for name, cfg_key, pred, gate, prov, notes in _FORCE_TERMS:
        raw = _cfg_get(cfg, cfg_key) if cfg_key else _MISSING
        intent = None if raw is _MISSING else bool(raw)
        if cell is not None:
            actual = bool(pred(cell))
            out.append(ForceChannel(
                name=name, kind="FORCE", enabled=actual, value=None, unit="bool",
                provenance=prov, code_site=_ACCUM_SITE, gate=gate, observation="runtime",
                config_intent=intent,
                intent_mismatch=(intent is not None and actual != intent), notes=notes))
        elif intent is not None:
            out.append(ForceChannel(
                name=name, kind="FORCE", enabled=intent, value=None, unit="bool",
                provenance=prov, code_site=_ACCUM_SITE, gate=gate, observation="config",
                config_intent=intent, notes=notes))
        else:
            # No runtime object and no config record.  The term is gated on built topology
            # (n_tri / n_xl / n_branch) or on a field this config predates: NOT observable.
            out.append(ForceChannel(
                name=name, kind="FORCE", enabled=None, value=None, unit="bool",
                provenance=prov, code_site=_ACCUM_SITE, gate=gate, observation="unobserved",
                notes=(notes + " | NOT RECORDED / no runtime object: state cannot be certified.")))
    return out


def _constant_channels() -> list[ForceChannel]:
    """Magnitudes read LIVE from the physics modules (they are module constants, not config)."""
    from aleph.components.incumbent import assemble as A
    from aleph.components.incumbent import compartments as C

    site_a = "aleph/components/incumbent/assemble.py"
    site_c = "aleph/components/incumbent/compartments.py"
    mem = C.MEMBRANE_DEFAULTS
    nuc = C.NUCLEUS_DEFAULTS
    v = _value_channel

    def const(name, value, unit, prov, site, notes, kind="SUBCHANNEL", enabled=True):
        return v(name=name, kind=kind, enabled=enabled, value=value, unit=unit, provenance=prov,
                 code_site=site, gate="module constant (always in force when its channel is on)",
                 observation="code-default", notes=notes)

    return [
        # ── geometry / turgor ────────────────────────────────────────────────────────────────
        const("cell_radius", A.R_CELL_UM, "um", "SOURCED", f"{site_a}::R_CELL_UM",
              "MCF7 suspended outer (membrane) radius; Wagner 2011 Coulter measurement."),
        const("cortex_membrane_gap", A.CORTEX_MEMBRANE_GAP_UM, "um", "DERIVED",
              f"{site_a}::CORTEX_MEMBRANE_GAP_UM",
              "= 1/2 * h_cortex (h_cortex ~200 nm, KB-3.1/3.5)."),
        const("osmotic_turgor", A.PI_0_PA, "Pa", "CONVENIENCE", f"{site_a}::PI_0_PA",
              "Resting osmotic turgor Pi_0. Constant Delta-pi (PI 2026-07-23 Option A: sigma=1, "
              "not coupled to volume/RVD). Pi_0=0 is the 2026-06-04 force-free-shell failure, so "
              "a NONZERO value is required — but the MAGNITUDE is NOT sourced for MCF7. "
              "AC_DECISION_CARDS_2026-07-22 §77, verbatim: '40 Pa stays a HeLa diagnostic proxy "
              "only'. Classified CONVENIENCE, not SOURCED, so no artifact can carry this number as "
              "though it were an MCF7 measurement. Pi_0 sets the ENTIRE resting cortical tension "
              "through gamma = Delta-P*R/2, so its provenance is load-bearing."),
        # ── Biot / fluid ─────────────────────────────────────────────────────────────────────
        const("biot_storage_S", A.BIOT_STORAGE_S, "um^2/pN", "DERIVED", f"{site_a}::BIOT_STORAGE_S",
              "S = 1/M; derived so mobility/S = c_v = 50 um^2/s (params_i0b1 anchor)."),
        const("biot_mobility", A.BIOT_MOBILITY, "um^4/(pN*s)", "DERIVED", f"{site_a}::BIOT_MOBILITY",
              "k/mu; derived jointly with S from the c_v anchor."),
        const("biot_alpha", A.BIOT_ALPHA, "dimensionless", "SOURCED", f"{site_a}::BIOT_ALPHA",
              "Biot-Willis coupling, PI-ratified 2026-07-16."),
        const("membrane_hydraulic_Lp", A.L_P, "um/(s*Pa)", "PI_GAP", f"{site_a}::L_P",
              "Membrane hydraulic conductivity; code labels this draft/MCF7, not an audited row."),
        # ── steric ───────────────────────────────────────────────────────────────────────────
        const("steric_sigma_ev", A.SIGMA_EV_UM, "um", "SOURCED", f"{site_a}::SIGMA_EV_UM",
              "F-actin steric diameter 7 nm. The physical candidate; the node-vs-segment "
              "discretization of EV remains a PI GAP (params_i0b2b)."),
        const("steric_k_ev", A.K_EV_PROVISIONAL, "pN/um", "CONVENIENCE",
              f"{site_a}::K_EV_PROVISIONAL",
              "WCA contact stiffness, self-labelled PROVISIONAL Magic-Number-Block scale; to be "
              "closed at the native gate from the measured passive F_op."),
        # ── membrane sub-channels ────────────────────────────────────────────────────────────
        const("membrane_bending_kappa", mem["kappa_m"], "pN*um", "SOURCED",
              f"{site_c}::MEMBRANE_DEFAULTS[kappa_m]", "Helfrich kappa_m = 20 kBT (KB-3.B1.2)."),
        const("membrane_area_tension_gamma", mem["gamma_mem"], "pN/um", "SOURCED",
              f"{site_c}::MEMBRANE_DEFAULTS[gamma_mem]",
              "In-plane bilayer tension plateau (KB-3.B1.1). The K_A reservoir upturn stays OFF "
              "(hard-truth #8) — that is a declared model limit, not a value."),
        const("membrane_erm_stiffness", mem["k_erm"], "pN/um", "SOURCED",
              f"{site_c}::MEMBRANE_DEFAULTS[k_erm]",
              "Single ezrin-F-actin linker 4.6 pN/nm, PI-ratified 2026-07-21 (Braunger 2014 "
              "single-bond analysis). KB-3.B1.6 registration PENDING, not yet SoT. Replaces the "
              "retired 1e2 bleb-PDE continuum value (~46x too soft)."),
        const("membrane_erm_reach", mem["erm_reach_um"], "um", "DERIVED",
              f"{site_c}::MEMBRANE_DEFAULTS[erm_reach_um]", "Membrane<->cortex ERM pairing reach."),
        # ── nucleus sub-channels: every lamina modulus is a declared TEST value ──────────────
        const("nucleus_chromatin_modulus", nuc["k_chrom"], "pN/um", "PI_GAP",
              f"{site_c}::NUCLEUS_DEFAULTS[k_chrom]", "I0-B2 PI GAP, TEST value."),
        const("nucleus_lamin_b_modulus", nuc["k_lamin_b"], "pN/um", "PI_GAP",
              f"{site_c}::NUCLEUS_DEFAULTS[k_lamin_b]", "I0-B2 PI GAP, TEST value."),
        const("nucleus_lamin_ac_modulus", nuc["k_lamin_ac"], "pN/um", "PI_GAP",
              f"{site_c}::NUCLEUS_DEFAULTS[k_lamin_ac]", "I0-B2 PI GAP, TEST value (5x ratio)."),
        const("nucleus_lamin_ac_knee", nuc["knee_strain"], "dimensionless", "PI_GAP",
              f"{site_c}::NUCLEUS_DEFAULTS[knee_strain]",
              "Lamin-A/C engagement strain; KB-3.B2.2 ~0.10 but code marks it TEST."),
        const("nucleus_envelope_bending", nuc["kappa_ne"], "pN*um", "PI_GAP",
              f"{site_c}::NUCLEUS_DEFAULTS[kappa_ne]", "~20 kBT envelope bending; TEST."),
        const("nucleus_volume_penalty", nuc["k_vol"], "pN/um^2", "PI_GAP",
              f"{site_c}::NUCLEUS_DEFAULTS[k_vol]",
              "Nucleoplasm bulk penalty (nu->1/2); TEST. p_vol = 0 at V = V0."),
        const("nucleus_linc_stiffness", nuc["k_linc"], "pN/um", "PI_GAP",
              f"{site_c}::NUCLEUS_DEFAULTS[k_linc]",
              "LINC tether STIFFNESS. PI GAP: the sourced 8 pN datum is a TENSION, not a stiffness."),
        const("nucleus_rupture_strain", nuc["eps_rupture"], "dimensionless", "PI_GAP",
              f"{site_c}::NUCLEUS_DEFAULTS[eps_rupture]",
              "Emergent envelope rupture strain; PI GAP, TEST, explicitly never tuned."),
        # ── NMII reference topology + TEST mechanics ─────────────────────────────────────────
        const("nmii_k_xb", A.NMII_K_XB_TEST, "pN/um", "PI_GAP", f"{site_a}::NMII_K_XB_TEST",
              "MASTER force knob, params_i0b3 PI GAP. Running at the mid physical-band TEST value; "
              "does not enter the resting result only because heads are UNBOUND at t0."),
        const("nmii_f_stall_head", A.NMII_F_STALL_TEST, "pN", "PI_GAP",
              f"{site_a}::NMII_F_STALL_TEST", "Per-head stall; params_i0b3 conflicted 0.5 vs 2.0."),
        const("nmii_v0", A.NMII_V0_TEST, "um/s", "PI_GAP", f"{site_a}::NMII_V0_TEST",
              "Unloaded head velocity; params_i0b3 conflicted 0.12 vs 0.2."),
        const("nmii_kappa_hill", A.NMII_KAPPA_TEST, "dimensionless", "PI_GAP",
              f"{site_a}::NMII_KAPPA_TEST", "Hill FV curvature; params_i0b3 THREE-way conflict."),
        const("nmii_k_on", A.NMII_KON_TEST, "1/s", "PI_GAP", f"{site_a}::NMII_KON_TEST",
              "Per-head attach rate; config-chosen, not an audited head-level row."),
        const("nmii_k_off0", A.NMII_KOFF0, "1/s", "PI_GAP", f"{site_a}::NMII_KOFF0",
              "Bell prefactor; params_i0b3 draft/provisional."),
        const("nmii_bell_force_f0", A.NMII_F0, "pN", "DERIVED", f"{site_a}::NMII_F0",
              "f0 = kBT/x_beta ~ 7.13 pN (physical, from the Veigel x_beta)."),
        const("nmii_capture_reach", A.NMII_CAPTURE_UM, "um", "DERIVED",
              f"{site_a}::NMII_CAPTURE_UM",
              "Head->actin perpendicular capture = r0_head + ~10 nm slack; reconciled 2026-07-24 "
              "with ff/hand_kmc capture_radius_um. NOT the loose 0.6 diagnostic probe reach."),
        const("nmii_n_side", A.NMII_N_SIDE, "count", "PI_GAP", f"{site_a}::NMII_N_SIDE",
              "Heads per anti-parallel side; params_i0b3 conflicted 10 (AFINES) vs 28-30 "
              "(Billington structural)."),
        const("nmii_backbone_length", A.NMII_L_BB_UM, "um", "PI_GAP", f"{site_a}::NMII_L_BB_UM",
              "Minifilament backbone contour (Billington 2013 EM candidate); params_i0b3 GAP."),
    ]


def _config_channels(cfg: Any) -> list[ForceChannel]:
    """BUILD/KINETIC/SUBCHANNEL channels carried by the run's ``CellConfig`` record."""
    site = "aleph/components/incumbent/assemble.py::CellConfig"
    out: list[ForceChannel] = []
    b = _bool_channel

    out.append(b(cfg, "overlap_free_cortex", name="overlap_free_cortex", kind="BUILD",
                 provenance="DERIVED", code_site=site,
                 notes="Build the cortex overlap-free (radial ~0.2 um thickness + WCA relaxation) "
                       "so the shell starts at ZERO steric force. Code: 'Default ON = the "
                       "PHYSIOLOGICAL BASELINE'. OFF reproduces the legacy zero-thickness "
                       "gamma-validation build and reintroduces the ~64k build interpenetrations "
                       "(steric-only max 2634 pN)."))
    out.append(b(cfg, "nmii_straddle_placement", name="nmii_straddle_placement", kind="BUILD",
                 provenance="DERIVED", code_site=site,
                 notes="Orient each bipolar minifilament to STRADDLE two anti-parallel filaments so "
                       "its heads land ON actin. Code: 'ON = the PHYSIOLOGICAL BASELINE'. OFF is "
                       "the legacy fixed perpendicular projection that left heads 0.13-0.6 um off "
                       "actin (only ~17% seedable at native)."))
    out.append(b(cfg, "erm_radial_pairing", name="membrane_erm_radial_pairing", kind="BUILD",
                 provenance="DERIVED", code_site=site,
                 notes="Pair each membrane node to its most-RADIAL cortex node so the tether "
                       "transmits the radial turgor load. Code default is False for "
                       "backward-compat and the comment states 'the resting physiological "
                       "baseline turns it on'."))
    out.append(b(cfg, "erm_density_mcf7_production", name="membrane_erm_density_mcf7_latch",
                 kind="BUILD", provenance="CONVENIENCE", code_site=site,
                 notes="Explicit provenance latch: asserts a supplied ERM density is an MCF7 "
                       "production contract, not a cross-cell proxy."))

    def num(key, name, unit, prov, notes, kind="BUILD"):
        raw = _cfg_get(cfg, key)
        if raw is _MISSING:
            out.append(ForceChannel(
                name=name, kind=kind, enabled=None, value=None, unit=unit, provenance=prov,
                code_site=site, gate=f"cfg.{key}", observation="unobserved",
                notes=notes + " | NOT RECORDED in this config: value cannot be certified."))
        else:
            out.append(ForceChannel(
                name=name, kind=kind, enabled=(raw is not None), value=raw, unit=unit,
                provenance=prov, code_site=site, gate=f"cfg.{key}", observation="config",
                config_intent=raw, notes=notes))

    num("n_filaments", "cortex_filament_count", "count", "SOURCED",
        "Full native cortical F-actin population = 100 um^-2 * 4*pi*(7.5 um)^2. PI HARD rule: "
        "never validate or debug at a reduced population.")
    num("cortex_seg_um", "cortex_segment_length", "um", "CONVENIENCE",
        "Segment rest length; caps the mesh/pore size. 0.5 is the COARSE baseline. Physiological "
        "cortical mesh is 50-100 nm (Morone 2006; Bovellan 2014; Chugh & Paluch 2018) => ~0.075.")
    num("cortex_length_um", "cortex_filament_length", "um", "CONVENIENCE",
        "Representative cortical filament contour. Native cortical filaments are ~hundreds of nm, "
        "so 3.0 um are COARSE representative (Cytosim-style) filaments, not a length census.")
    num("cortex_density_per_fil", "cortex_crosslinks_per_filament", "count", "DERIVED",
        "Crosslinks per filament => contour spacing = length/density. 20 keeps the network "
        "single-spanning (percolation fix 2026-07-23; degree-2 at 1.0 fragmented it).")
    num("cortex_arp23_fraction", "cortex_arp23_mass_fraction", "dimensionless", "SOURCED",
        "Arp2/3 fraction of cortical actin BY MASS. Native is ~1/3 (Bovellan 2014, KB-3.18); the "
        "code default 0.0 is a pure-formin cortex, i.e. a declared model simplification.")
    num("membrane_subdivisions", "membrane_mesh_subdivisions", "count", "CONVENIENCE",
        "Plasma-membrane icosphere level. Code default 3 = 642 verts is the backward-compat / "
        "gamma-validation smooth resting sphere. PI-ratified 2026-07-24: 8 = 655,362 verts (~33 nm) "
        "is the production RESTING config (drops the Helfrich residual 0.7766 -> 0.0768 by grid "
        "convergence). PI 2026-07-22: 6 = dynamic baseline, 7 = validation.")
    num("nucleus_subdivisions", "nucleus_mesh_subdivisions", "count", "CONVENIENCE",
        "Nuclear-envelope icosphere level. PI 2026-07-22: 6 dynamic baseline / 7 validation.")
    num("cortex_overlap_mode", "cortex_overlap_resolution_mode", "enum", "DERIVED",
        "How overlap_free_cortex resolves build interpenetrations. 'transverse' (default, "
        "byte-identical to history) shoves crossing nodes IN-plane, which knocks a node off its "
        "great-circle arc = a bending KINK whose force blows up x288 at the fine 75 nm mesh. "
        "'radial_span' separates OUT-of-plane over a smooth cosine bump.")
    num("cortex_overlap_span", "cortex_overlap_span", "count", "CONVENIENCE",
        "radial_span bump half-window in arc-neighbours.")
    num("steric_force_cap", "steric_force_cap", "pN", "CONVENIENCE",
        "WCA near-core CFL guard. A DECLARED numerical guard (>> the ~0.02 pN physical operating "
        "load), never a physics magnitude tuned to an outcome. It is NOT the steric verdict.")
    num("dx_um", "field_grid_resolution", "um", "CONVENIENCE",
        "Conservative field-grid resolution (numerical sizing; not literature).")
    num("grid_pad_um", "field_grid_pad", "um", "CONVENIENCE", "Grid half-extent beyond R_cell.")
    num("steric_grid_dim", "steric_hashgrid_dim", "count", "CONVENIENCE",
        "Hash-grid dimension; accelerator only, grid-invariant.")
    num("motor_segment_grid_dim", "motor_segment_hashgrid_dim", "count", "CONVENIENCE",
        "Segment-midpoint HashGrid dimension; accelerator only.")
    num("erm_density_per_um2", "membrane_erm_areal_density", "um^-2", "PI_GAP",
        "No MCF7 ERM areal-density datum is registered. None keeps a diagnostic "
        "one-per-membrane-node topology and NG-3 must remain open.")
    num("nmii_backbone_lp_um", "nmii_backbone_persistence_length", "um", "PI_GAP",
        "Backbone persistence length is a physical magnitude, not a rigid-rod convenience. No "
        "ratified NMII-minifilament value exists, so production-form bending is opt-in and "
        "source-gated; None leaves the gap explicit.")
    num("seed", "rng_seed", "count", "CONVENIENCE", "Build RNG seed.")

    # ── all-or-none KINETIC contracts: enabled only when EVERY required field is supplied ──
    def all_or_none(name, keys, prov, notes):
        vals = {k: _cfg_get(cfg, k) for k in keys}
        if all(x is _MISSING for x in vals.values()):
            out.append(ForceChannel(
                name=name, kind="KINETIC", enabled=None, value=None, unit="bool", provenance=prov,
                code_site=site, gate=f"all of {list(keys)} not None", observation="unobserved",
                notes=notes + " | NOT RECORDED in this config: state cannot be certified."))
            return
        present = {k: (None if x is _MISSING else x) for k, x in vals.items()}
        on = all(present[k] is not None for k in keys)
        out.append(ForceChannel(
            name=name, kind="KINETIC", enabled=on, value=present, unit="mixed", provenance=prov,
            code_site=site, gate=f"all of {list(keys)} not None", observation="config",
            config_intent=present, notes=notes))

    all_or_none(
        "membrane_erm_bell_kinetics",
        ("erm_k_on_s", "erm_k_off0_s", "erm_bell_force_pn", "erm_capture_radius_um"),
        "PI_GAP",
        "ERM on/off kinetics are an all-or-none, source-gated contract with NO biological "
        "defaults: MCF7 k_on/k_off0/F0/capture are Contract-Graph GAPs (2026-07-21 audit). With "
        "them unset the ERM tethers are permanent springs, not a kinetic linkage.")
    all_or_none(
        "resting_bound_myosin_setpoint",
        ("resting_bound_myosin_fraction", "resting_bound_myosin_force_pn"),
        "PI_GAP",
        "The resting cortical-tension source (fraction x per-head force sets gamma_cortex = "
        "Delta-P*R/2 - gamma_mem). BOTH magnitudes are PI GAPs absent from the Contract-Graph and "
        "are never defaulted. When ON the resting path also engages radial ERM pairing. Any "
        "non-null pair MUST carry resting_bound_myosin_source.")

    for key, name in (("erm_density_source", "membrane_erm_density_source"),
                      ("erm_kinetics_source", "membrane_erm_kinetics_source"),
                      ("nmii_backbone_lp_source", "nmii_backbone_lp_source"),
                      ("resting_bound_myosin_source", "resting_bound_myosin_source")):
        raw = _cfg_get(cfg, key)
        if raw is _MISSING:
            out.append(ForceChannel(
                name=name, kind="BUILD", enabled=None, value=None, unit="str",
                provenance="CONVENIENCE", code_site=site, gate=f"cfg.{key}",
                observation="unobserved",
                notes="Provenance label. NOT RECORDED in this config."))
        else:
            out.append(ForceChannel(
                name=name, kind="BUILD", enabled=bool(raw), value=raw, unit="str",
                provenance="CONVENIENCE", code_site=site, gate=f"cfg.{key}",
                observation="config", config_intent=raw,
                notes="Provenance label required alongside its source-gated magnitude. An empty "
                      "label on a non-null magnitude means the magnitude is unattributed."))
    return out


def _solver_channels(solver_kwargs: dict[str, Any] | None, env: dict[str, str] | None,
                     env_observed: bool) -> list[ForceChannel]:
    """NUMERIC solver paths.  These live in ``ProjectedAnalyticCG`` kwargs + env vars, NOT in
    ``CellConfig`` — so a run's report JSON does not record them at all."""
    site = "aleph/components/incumbent/implicit_mechanics.py::ProjectedAnalyticCG.__init__"
    out: list[ForceChannel] = []
    kw = solver_kwargs or {}
    env = env if env is not None else {}

    specs = (
        ("fiber_quotient_coarse", "fiber_quotient_coarse", "bool", False, "AC_FQ_COARSE",
         "Fiber-quotient inter-fiber coarse deflation (Path A/B). Preconditioning only (SPD, "
         "cannot move the fixed point). DEFAULT OFF. GATE A 2026-07-24 finding: it does NOT close "
         "the native resting plateau — the residual is discrete-myosin local non-equilibrium, not "
         "a collective mode. Kept flagged OFF."),
        ("fq_coarse_mode", "fq_coarse_mode", "enum", "pathA", "AC_FQ_MODE",
         "pathA = matrix-free P^T A P block-Jacobi inner CG; pathB = explicit crosslink-weighted "
         "BSR coarse operator with deep block-Jacobi CG."),
        ("fq_coarse_iterations", "fq_coarse_iterations", "count", 40, "AC_FQ_ITERS",
         "Inner iterations of the fiber-quotient coarse solve."),
        ("multigrid", "multigrid", "bool", False, "AC_MG",
         "Fiber-arclength geometric-multigrid V-cycle preconditioner. DEFAULT OFF. When ON it "
         "REPLACES the incumbent smoother+coarse stack; SPD so the fixed point and residual gate "
         "are identical. This is the path that broke the fine 75 nm (2.9M node) plateau."),
        ("mg_fiber_quotient", "mg_fiber_quotient", "bool", True, "AC_MG_FQ",
         "Use the fiber-quotient level as the multigrid coarsest solver. Note the env override is "
         "INVERTED (AC_MG_FQ=0 disables); default True."),
        ("mg_assembled_coarse", "mg_assembled_coarse", "bool", False, "AC_MG_ASM",
         "SPEED FIX2: assembled per-fiber coarse operator (~40x with FIX1). DEFAULT OFF pending a "
         "native A/B comparison."),
        ("cg_check_every", "cg_check_every", "count", 0, "",
         "SPEED FIX1: host early-exit cadence for the outer PCG. 0 = incumbent fixed-budget loop. "
         "The returned dx is BIT-IDENTICAL to running to max_iterations. DEFAULT OFF."),
        ("multigrid_smooth", "mg_smooth", "count", 2, "", "Multigrid pre/post smoothing sweeps."),
        ("multigrid_omega", "mg_omega", "dimensionless", 0.8, "", "Multigrid smoother relaxation."),
        ("implicit_coarse_modes", "coarse_modes", "count", 0, "",
         "Global rigid-body + constant-strain (l<=2) coarse deflation modes. 0 = off, 12 = full. "
         "Preconditioner-only; never moves the fixed point or relaxes a residual gate."),
    )
    for name, kwarg, unit, default, env_var, notes in specs:
        supplied = kwarg in kw
        env_set = bool(env_var) and env_var in env
        if not supplied and not env_observed:
            out.append(ForceChannel(
                name=name, kind="NUMERIC", enabled=None, value=None, unit=unit,
                provenance="CONVENIENCE", code_site=site,
                gate=(f"kwarg {kwarg} or env {env_var}" if env_var else f"kwarg {kwarg}"),
                observation="unobserved",
                notes=notes + " | NOT OBSERVABLE: this path is a solver kwarg + env var, and the "
                              "driver report JSON records neither. State cannot be certified."))
            continue
        if supplied:
            value, obs = kw[kwarg], "runtime"
        elif env_set:
            raw = env[env_var]
            if unit == "bool":
                value = (raw != "0") if env_var == "AC_MG_FQ" else (raw == "1")
            elif unit == "count":
                value = int(raw)
            else:
                value = raw
            obs = "env"
        else:
            value, obs = default, "code-default"
        enabled = bool(value) if unit == "bool" else (value != default if unit == "count" else True)
        if unit == "enum":
            enabled = True
        out.append(ForceChannel(
            name=name, kind="NUMERIC", enabled=enabled, value=value, unit=unit,
            provenance="CONVENIENCE", code_site=site,
            gate=(f"kwarg {kwarg} or env {env_var}" if env_var else f"kwarg {kwarg}"),
            observation=obs, config_intent=default,
            intent_mismatch=(obs == "env" and value != default),
            notes=notes + (f" | ENV OVERRIDE {env_var}={env[env_var]!r} ACTIVE" if obs == "env" else "")))
    return out


def dump_force_channels(
    cfg: Any = None, *, cell: Any = None, solver_kwargs: dict[str, Any] | None = None,
    env: dict[str, str] | None = None, capture_env: bool = False,
    profile_claim: str = "", run_label: str = "",
) -> dict[str, Any]:
    """Observe every force channel and flag-gated path and return a JSON-able dump.

    Args:
        cfg: A ``CellConfig`` instance or the ``config`` dict from a run's report JSON.  ``None``
            reports every config-carried channel as ``unobserved``.
        cell: The live composed ``AssembledCell``.  When supplied, FORCE-channel state is read off
            the ACTUALLY COMPOSED objects rather than the config's intent, and any disagreement is
            recorded as ``intent_mismatch``.  This is the only observation strong enough to have
            caught ``SettlingForce``.
        solver_kwargs: The kwargs actually passed to ``ProjectedAnalyticCG`` for this run.
        env: Environment mapping to read override vars from.  Defaults to ``os.environ`` when
            ``capture_env`` is True.
        capture_env: Set True only when observing a LIVE process whose environment is the run's
            environment.  Reconstructing a dump from a committed artifact must leave this False —
            this process's env is not the historical run's env.
        profile_claim: Which declared production profile this run claims to satisfy (e.g.
            ``"resting_production"``).  An empty claim is itself reportable.
        run_label: Free-text identifier for the run/artifact being observed.

    Returns:
        A dict with ``schema``, ``profile_claim``, ``run_label``, ``coverage`` and ``channels``.
    """
    if env is None and capture_env:
        env = {k: v for k, v in os.environ.items() if k in ENV_OVERRIDE_VARS}
    env_observed = capture_env or env is not None

    channels: list[ForceChannel] = []
    channels += _force_channels(cfg, cell)
    channels += _constant_channels()
    channels += _config_channels(cfg)
    channels += _solver_channels(solver_kwargs, env, env_observed)

    names = [c.name for c in channels]
    dupes = sorted({n for n in names if names.count(n) > 1})
    if dupes:
        raise ValueError(f"duplicate channel names in dump: {dupes}")

    coverage = {
        "n_channels": len(channels),
        "by_observation": {o: sum(1 for c in channels if c.observation == o) for o in OBSERVATIONS},
        "by_kind": {k: sum(1 for c in channels if c.kind == k) for k in KINDS},
        "by_provenance": {p: sum(1 for c in channels if c.provenance == p) for p in PROVENANCE},
        "n_intent_mismatch": sum(1 for c in channels if c.intent_mismatch),
        "env_observed": bool(env_observed),
        "runtime_observed": cell is not None,
        "assembly_terms": assert_force_channels_cover_accumulate(),
    }
    return {
        "schema": SCHEMA,
        "profile_claim": profile_claim,
        "run_label": run_label,
        "coverage": coverage,
        "channels": [asdict(c) for c in channels],
    }


def dump_from_artifact(path: str | Path, *, profile_claim: str = "") -> dict[str, Any]:
    """Reconstruct a dump from a committed run report JSON (``--json`` output of the driver).

    The artifact carries ``config`` only, so solver + env channels come back ``unobserved`` and
    FORCE channels are ``config``-strength at best.  That weakness is REPORTED, not papered over:
    a config record is intent, and the ``SettlingForce`` failure was precisely intent diverging
    from runtime.
    """
    p = Path(path)
    doc = json.loads(p.read_text())
    dump = dump_force_channels(
        doc.get("config"), profile_claim=profile_claim, run_label=str(p))
    dump["artifact"] = {"path": str(p), "schema": doc.get("schema", "")}
    return dump


def write_dump(path: str | Path, dump: dict[str, Any]) -> Path:
    """Write ``dump`` as indented JSON and return the path."""
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(dump, indent=2, default=str) + "\n")
    return p


def _main() -> None:
    import argparse

    ap = argparse.ArgumentParser(
        description="Observe every force channel / flag-gated path (read-only; changes no physics).")
    ap.add_argument("--artifact", type=str, default="",
                    help="run report JSON to reconstruct a dump from (config-strength observation)")
    ap.add_argument("--defaults", action="store_true",
                    help="dump the CURRENT CellConfig code defaults (what an unargumented run does)")
    ap.add_argument("--profile-claim", type=str, default="",
                    help="the declared production profile this run claims to satisfy")
    ap.add_argument("--capture-env", action="store_true",
                    help="record this process's override env vars (LIVE runs only)")
    ap.add_argument("--out", type=str, default="", help="write the dump JSON here")
    args = ap.parse_args()

    if args.artifact:
        dump = dump_from_artifact(args.artifact, profile_claim=args.profile_claim)
    elif args.defaults:
        from aleph.components.incumbent.assemble import CellConfig
        dump = dump_force_channels(
            CellConfig(), profile_claim=args.profile_claim, capture_env=args.capture_env,
            run_label="CellConfig() code defaults")
    else:
        dump = dump_force_channels(
            None, profile_claim=args.profile_claim, capture_env=args.capture_env,
            run_label="no config supplied")

    if args.out:
        print(f"wrote {write_dump(args.out, dump)}")
    else:
        print(json.dumps(dump, indent=2, default=str))


if __name__ == "__main__":
    _main()
