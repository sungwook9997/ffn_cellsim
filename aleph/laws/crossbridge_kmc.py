r"""NMII head->actin crossbridge kinetics: the catch-slip off-rate, and the constants that are absent.

THE LAW, AND WHY IT IS NOT BELL.  ``components/motor/hand.py`` records the fidelity correction of
2026-07-23 in as many words: *"the NMII head-actin bond is CATCH-SLIP, not pure Bell slip. Kovacs,
Thirumurugan, Knight & Sellers 2007, PNAS 104(24):9994"* — resistive load SLOWS ADP release and
therefore SLOWS detachment, up to a peak-lifetime load ``F*``, past which forced unbinding takes over.
That is the OPPOSITE sign of a slip bond at low load. ``hand.step_detach_kernel`` stays a pure Bell
slip and a test holds it that way, because the ERM / alpha-actinin / crosslink consumers depend on it;
it is not the NMII head's law and this module does not use it.

WHY THIS MODULE EXISTS AT ALL, given three Pereverzev implementations already in the tree.  The device
form lives in ``components/motor/hand.py`` and ``world/`` **may not import** ``components/``
(``tests/architecture/test_layer_directions.py``); ``laws/fa_clutch_warp._pereverzev`` is the integrin
clutch's, private, and takes a SIGNED load where a crossbridge slips under its magnitude; and
``laws/hand_kmc.NMIIA_MYOSIN`` — the one preset actually named for this molecule — is
``catch_slip=False``, i.e. **the exact law the 2026-07-23 correction says is wrong for this head.**
⚠ That preset is reported here rather than edited: it is another lane's file and other consumers read
it. The host arithmetic below is not re-derived either — :meth:`CrossbridgeKinetics.off_rate` calls
``hand_kmc.pereverzev_off_rate``, so there is one host source of truth and this module adds the NMII
scope, the device mirror and the refusals.

WHAT IS ANSWERED HERE: the FORM and the SIGN.  ``d k_off / d|F| < 0`` at low load and ``> 0`` past
``F*`` is a claim about mechanism, it is sourced, and it is testable without any magnitude. It is
recorded as a law rather than as a number.

WHAT IS REFUSED: EVERY MAGNITUDE.  ``k_catch0``, ``x_catch``, ``k_slip0``, ``x_slip``, ``k_on`` and
``capture_radius_um`` are I0-B3 PI-GAPs — ``hand.NMIICatchSlipParams``' own docstring says *"Every
value is an I0-B3 GAP - PI (magnitudes constrained by Kovacs 5x/12x slowdown + NM2B duty 0.2-0.3, but
no direct single-molecule fit); NEVER chosen to lift a force/gamma band."* So
:class:`CrossbridgeKinetics` has **no defaults**: every field is required and the refusal names the
gap. ``scripts/ac_gate_b_cortex_motor_native.CATCH_SLIP_PROVISIONAL`` is what a run has been passing;
it is labelled provisional there and this module will not adopt it as a default.

⚠ AND ONE STRUCTURAL REFUSAL THAT IS NOT A MISSING NUMBER.  A parameter set whose catch branch does
not dominate at low load is a SLIP bond wearing a catch bond's name — its off-rate rises monotonically
and there is no interior ``F*`` at all. ``k_catch0*x_catch > k_slip0*x_slip`` is therefore checked at
construction, not at use, so a set that cannot exhibit the Kovacs sign is rejected by name rather than
silently producing a Bell curve under a catch-slip label. This is the same inconsistency
``hand_kmc.INTEGRIN_A5B1`` carries as a KB-CONSISTENCY FLAG, caught one layer earlier.

THE KERNEL PROPOSES; IT DOES NOT COMMIT.  ``CLAUDE.md`` §Architectural principle: *"a kinetic connector
commits ONLY on an accepted physical step under one device-resident transaction."* PI decision 18
(2026-08-21) leaves the acceptance predicate UNDEFINED. So
:func:`propose_crossbridge_transitions_kernel` reads the committed ``bound`` array and writes a
SEPARATE ``proposal`` array — it is structurally incapable of mutating committed state, which is the
cheapest possible form of that transaction and needs no snapshot twin to roll back. Whether the
proposal is ever applied is the caller's, and the family's, business.

Sanity Gate (recorded per CLAUDE.md, before first execution):
  * dimensional — ``k_catch0``/``k_slip0``/``k_on`` [1/s]; ``x_catch``/``x_slip``/``capture_radius_um``
    [µm]; ``kT`` [pN·µm]; loads [pN]; ``tau`` [s]. ``|F|*x/kT`` is dimensionless, which is what makes
    the exponent legal, and a Bell length given in nm instead of µm is a 1000x error in the exponent.
  * boundary — a zero load gives ``k_catch0 + k_slip0``, the unloaded off-rate, finite and positive; a
    load large enough to overflow ``exp`` is not guarded because it is not reachable: the slip exponent
    at 100 pN with a 1 nm bond length is 23. ``tau = 0`` gives probability 0, not a division.
  * conservation/invariant — a proposal changes at most one bond's state per tick and never touches
    another bond's slot; the kernel writes ``proposal[t]`` and ``tension_pn[t]`` for its own ``t``
    alone, with no atomics and no cross-thread read.
  * CFL/precision — no integration; float64 throughout. The KMC tick ``tau`` is the OUTER physical
    step, not the inner mechanical iteration, so nothing here enters a CFL bound.
  * sign sense — THE claim of this module. ``off_rate`` FALLS with load below ``F*`` and RISES above
    it. The pure-slip degeneracy is refused at construction so the sign cannot be silently lost.
  * measurement protocol — host arithmetic is numpy and is for oracles and tests; the runtime path is
    one kernel launch of ``dim = n_bonds`` writing two device arrays and reading none back.

engine units: force pN, length µm, rate 1/s, energy pN·µm.  Runtime: NVIDIA Warp on CUDA; the host
helpers are CPU-importable so a card is not needed to check the law's shape.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import warp as wp

from aleph.laws import units as U
from aleph.laws.hand_kmc import pereverzev_off_rate

__all__ = [
    "NMII_CROSSBRIDGE_GAPS", "CrossbridgeKinetics", "crossbridge_off_rate",
    "propose_crossbridge_transitions_kernel",
]

#: What is unresolved behind each required field, quoted in the refusal so a caller cannot pick a
#: value without reading what is and is not known about it. Keyed by field name.
NMII_CROSSBRIDGE_GAPS: dict[str, str] = {
    "k_catch0": "zero-force catch-pathway rate [1/s]. I0-B3 PI-GAP — hand.NMIICatchSlipParams: "
                "'Every value is an I0-B3 GAP - PI ... NEVER chosen to lift a force/gamma band'. "
                "Constrained by the Kovacs 2007 5x/12x load slowdown, not fitted to it. "
                "ac_gate_b_cortex_motor_native.CATCH_SLIP_PROVISIONAL passes 0.35 and says provisional.",
    "x_catch_um": "catch-pathway bond length [um]. I0-B3 PI-GAP; the provisional run value is 1.0e-3 "
                  "(1 nm). Larger => stronger load-strengthening.",
    "k_slip0": "zero-force slip-pathway rate [1/s]. I0-B3 PI-GAP; provisional run value 0.35.",
    "x_slip_um": "slip-pathway (Bell) bond length [um]. I0-B3 PI-GAP; provisional run value 0.6e-3, "
                 "which is the Veigel 2002 x_beta that hand_kmc's Bell NMIIA preset uses — a Bell "
                 "length borrowed into the slip branch of a catch-slip law is a transfer, not a source.",
    "k_on": "per-head actin attachment rate [1/s]. params_i0b3.yaml lists k_on in "
            "'still_gap_under_nm2b'; its recorded source is the ff/hand_kmc NMIIA preset 50/s from "
            "configs/phase1_h3.yaml, i.e. a configuration file, not a measurement.",
    "capture_radius_um": "head->actin capture radius epsilon [um]. hand_kmc's NMIIA preset carries "
                         "210 nm as head_actin_capture_perp; PI-GAP, and it is a KINETIC radius, not "
                         "the geometric reach a builder pairs with.",
}


def _require(name: str, value: float | None) -> float:
    """Return a declared PI-GAP magnitude, or refuse with what is unresolved about it named."""
    if value is None:
        raise ValueError(
            f"{name} has no default and must be declared. PI-GAP — {NMII_CROSSBRIDGE_GAPS[name]}"
        )
    v = float(value)
    if not np.isfinite(v) or v <= 0.0:
        raise ValueError(f"{name} must be finite and positive; got {value!r}")
    return v


@dataclass(frozen=True, slots=True)
class CrossbridgeKinetics:
    """The NMII head->actin catch-slip law, with its scope and every magnitude declared.

    No field has a default. ``scope`` and ``provenance`` are required for the same reason
    :class:`~aleph.world.bond.BondCount` requires them: a magnitude with no stated origin is how a
    provisional value acquires a physiological label.

    Attributes:
        k_catch0: zero-force catch-pathway rate [1/s].
        x_catch_um: catch-pathway bond length [µm].
        k_slip0: zero-force slip-pathway rate [1/s].
        x_slip_um: slip-pathway (Bell) bond length [µm].
        k_on: per-head attachment rate while an acceptor is inside ``capture_radius_um`` [1/s].
        capture_radius_um: the KINETIC capture radius ε [µm]. Distinct from a builder's geometric
            reach: the reach decides what could ever be paired, this decides what binds this tick.
        scope: the cell line, state, assay and temperature these magnitudes are true FOR.
        provenance: the citation or derivation, per field if they differ.
        kT: thermal energy [pN·µm]. The one value here that is NOT a gap — it is the physiological
            temperature, so it defaults to :data:`aleph.laws.units.KBT`.
    """

    k_catch0: float
    x_catch_um: float
    k_slip0: float
    x_slip_um: float
    k_on: float
    capture_radius_um: float
    scope: str
    provenance: str
    kT: float = U.KBT

    def __post_init__(self) -> None:
        for name in ("k_catch0", "x_catch_um", "k_slip0", "x_slip_um", "k_on", "capture_radius_um"):
            _require(name, getattr(self, name))
        if not np.isfinite(self.kT) or self.kT <= 0.0:
            raise ValueError(f"kT must be finite and positive; got {self.kT!r}")
        if not str(self.scope).strip():
            raise ValueError(
                "scope is required: a rate measured on skeletal myosin or on NM2B is not thereby a "
                "rate for the head this bond carries, and without saying which, a transfer cannot be "
                "distinguished from a source."
            )
        if not str(self.provenance).strip():
            raise ValueError("provenance is required and must be non-empty")
        if self.k_catch0 * self.x_catch_um <= self.k_slip0 * self.x_slip_um:
            raise ValueError(
                "k_catch0*x_catch <= k_slip0*x_slip: this parameter set is SLIP-DOMINATED at every "
                "load, so its off-rate rises monotonically and there is no interior peak F*. That is "
                "a Bell bond wearing a catch bond's name — the Kovacs 2007 correction of 2026-07-23 "
                "is precisely that the NMII head is NOT that. Refused at construction rather than at "
                "use, because at use it produces a plausible curve with the wrong sign."
            )

    def off_rate(self, f_pn):
        """Off-rate ``p_off(|F|)`` [1/s]. Host oracle; delegates to the ``laws`` source of truth.

        Args:
            f_pn: crossbridge tension magnitude [pN], scalar or array. The MAGNITUDE: a crossbridge
                slips under its total force, and a signed load would make the catch branch a slip
                branch for one sign of it.
        """
        return pereverzev_off_rate(np.abs(np.asarray(f_pn, np.float64)), self.k_catch0,
                                   self.x_catch_um, self.k_slip0, self.x_slip_um, self.kT)

    @property
    def unloaded_off_rate(self) -> float:
        """``p_off(0) = k_catch0 + k_slip0`` [1/s]."""
        return float(self.k_catch0 + self.k_slip0)

    @property
    def peak_force_pn(self) -> float:
        r"""The peak-lifetime load ``F* = kT/(x_c+x_s) * ln[(k_c0 x_c)/(k_s0 x_s)]`` [pN].

        Where the catch branch hands over to the slip branch, i.e. where the bond is longest-lived.
        Guaranteed to exist and be positive because ``__post_init__`` refuses a slip-dominated set.

        ⚠ **This is an ORACLE, not a gate that has been run.** When the PI closes the four constants,
        ``F*`` becomes comparable against the ratified ``F_stall_head`` = 2.0 pN, and a set whose peak
        lifetime sits far from the stall force is saying something about the motor that its authors
        may not have intended. The comparison is written down here BEFORE the constants exist, which
        is the only order in which it is a contract.
        """
        return float(self.kT / (self.x_catch_um + self.x_slip_um)
                     * np.log((self.k_catch0 * self.x_catch_um) / (self.k_slip0 * self.x_slip_um)))

    @property
    def unloaded_bound_fraction(self) -> float:
        """``k_on / (k_on + p_off(0))`` — the EMERGENT unloaded engaged fraction, dimensionless.

        Reported as a cross-check against the NM2B duty band 0.2-0.3 (Nagy 2013), never imposed:
        ``params_i0b3.yaml`` is explicit that "the engaged fraction is EMERGENT from k_on/k_off ...
        NOT imposed as this duty".

        Provenance confirmed against the Contract-Graph 2026-08-21, per ``CLAUDE.md`` ("before citing
        a KB source in a deliverable, confirm its ``source_audit.verdict`` is OK"): **SE549**
        ``Nagy2013_JBC``, source_type *Direct measurement*, ``source_audit.verdict = OK``, DOI
        ``10.1074/jbc.M112.424671`` — the same DOI ``params_i0b3.yaml`` carries, and it is LINKED: 8 edges
        to 4 KnowledgeClaims.

        ⚠ **A retraction is recorded here rather than quietly overwritten.** This docstring said for
        one commit that SE549 had ZERO edges and was "audited and orphan". That was wrong, and the
        cause was a join, not the graph: ``edges`` keys on the Notion page UUID in
        ``source_evidence.id``, while ``SE549`` is the ``uid`` column. Joining on ``uid`` returns
        zero for EVERY row — 521 of the 565 sources are in fact linked. A query that returns zero
        for every input is not a finding about any input, and the shape of the answer should have
        said so before it was reported.

        What remains true and is a different statement: there is no **code** edge from this module to
        a contract. That needs an ``Implements: KU-x.y`` line naming a KnowledgeClaim, and choosing
        which of SE549's four to name is a gate-contract decision, so it is PI-authored
        (``CLAUDE.md`` §Knowledge base) and is reported rather than created here.
        """
        return float(self.k_on / (self.k_on + self.unloaded_off_rate))

    def record(self) -> dict[str, object]:
        """The artifact row — every magnitude beside its scope and provenance, and the two oracles."""
        return {
            "law": "pereverzev_catch_slip", "chemistry": "nmii_head_actin_kovacs2007",
            "k_catch0": self.k_catch0, "x_catch_um": self.x_catch_um,
            "k_slip0": self.k_slip0, "x_slip_um": self.x_slip_um,
            "k_on": self.k_on, "capture_radius_um": self.capture_radius_um, "kT": self.kT,
            "scope": self.scope, "provenance": self.provenance,
            "unloaded_off_rate_per_s": self.unloaded_off_rate,
            "peak_force_pn": self.peak_force_pn,
            "unloaded_bound_fraction": self.unloaded_bound_fraction,
        }


@wp.func
def crossbridge_off_rate(f: wp.float64, k_catch0: wp.float64, x_catch: wp.float64,
                         k_slip0: wp.float64, x_slip: wp.float64, kT: wp.float64) -> wp.float64:
    """Device mirror of :meth:`CrossbridgeKinetics.off_rate` — same form, same magnitude convention.

    ``p_off = k_catch0*exp(-|f|*x_catch/kT) + k_slip0*exp(+|f|*x_slip/kT)`` [1/s].
    """
    fa = wp.abs(f)
    return k_catch0 * wp.exp(-fa * x_catch / kT) + k_slip0 * wp.exp(fa * x_slip / kT)


@wp.kernel
def propose_crossbridge_transitions_kernel(
    pos: wp.array(dtype=wp.vec3d),
    node_head: wp.array(dtype=wp.int32),
    node_actin: wp.array(dtype=wp.int32),
    bound: wp.array(dtype=wp.int32),
    k_xb: wp.float64,
    r0_xb: wp.float64,
    tau: wp.float64,
    seed: wp.int32,
    k_catch0: wp.float64,
    x_catch: wp.float64,
    k_slip0: wp.float64,
    x_slip: wp.float64,
    kT: wp.float64,
    k_on: wp.float64,
    capture_radius: wp.float64,
    proposal: wp.array(dtype=wp.int32),
    tension_pn: wp.array(dtype=wp.float64),
):
    """Propose one KMC transition per crossbridge. Reads ``bound``; writes ``proposal``, never ``bound``.

    A bound head detaches with probability ``1 - exp(-tau * p_off(|F|))`` where ``|F| = k_xb*|L - r0_xb|``
    is the crossbridge tension magnitude; a free head attaches with ``1 - exp(-tau * k_on)`` when its
    actin site is inside ``capture_radius``. Both are Poisson per-tick probabilities (NF2007 §10.1).

    ⚠ **The separate output array IS the transaction.** ``CLAUDE.md`` requires a kinetic connector to
    commit only on an accepted physical step, and PI decision 18 (2026-08-21) leaves the acceptance
    predicate undefined — so this kernel cannot write committed state at all, and a rejected step is
    handled by not copying rather than by rolling back. No snapshot twin, because there is nothing to
    restore.

    ⚠ **The tension is AXIAL only.** ``STATE.md`` (f), 2026-08-15: *"The crossbridge BOND has no
    transverse channel; the head is not thereby free."* Recorded on the kernel that computes it rather
    than in a note beside it — ``tension_pn`` is the magnitude along the head-actin line and is the
    load the off-rate sees, so anything the transverse channel would add is missing from BOTH.
    """
    t = wp.tid()
    h = node_head[t]
    a = node_actin[t]
    d = pos[a] - pos[h]
    L = wp.length(d)

    if bound[t] == 0:
        tension_pn[t] = wp.float64(0.0)
        p_on = wp.float64(1.0) - wp.exp(-tau * k_on)
        state = wp.rand_init(seed, t)
        if L <= capture_radius and wp.float64(wp.randf(state)) < p_on:
            proposal[t] = wp.int32(1)
        else:
            proposal[t] = wp.int32(0)
        return

    f = k_xb * wp.abs(L - r0_xb)
    tension_pn[t] = f
    p_off = crossbridge_off_rate(f, k_catch0, x_catch, k_slip0, x_slip, kT)
    p_detach = wp.float64(1.0) - wp.exp(-tau * p_off)
    state = wp.rand_init(seed, t)
    if wp.float64(wp.randf(state)) < p_detach:
        proposal[t] = wp.int32(0)
    else:
        proposal[t] = wp.int32(1)


def _codegen_check() -> int:
    """Type-check this module's kernel by emitting its CUDA source, with no device and no execution.

    The dev machine has no card, so the alternative is shipping kernel syntax nobody has put in front
    of a compiler. This generates the CUDA C++ Warp would compile, which type-checks every expression,
    and produces NO number — so it is not a CPU result under any reading.
    """
    from warp._src.context import ModuleBuilder  # private, deliberately: no public codegen API

    module = wp.get_module(__name__)
    return len(ModuleBuilder(module, module.options).codegen("cuda"))


def _demo() -> None:
    """Self-check: the SIGN is the claim, every magnitude refuses, and a slip set is rejected by name."""
    provisional = dict(
        k_catch0=0.35, x_catch_um=1.0e-3, k_slip0=0.35, x_slip_um=0.6e-3,
        k_on=50.0, capture_radius_um=0.210,
        scope="PROVISIONAL — no scope. Self-check fixture only, NOT a run configuration.",
        provenance="ac_gate_b_cortex_motor_native.CATCH_SLIP_PROVISIONAL (labelled provisional there) "
                   "+ hand_kmc NMIIA k_on/epsilon. Every one is an I0-B3 PI-GAP.")

    # Every magnitude is required, and the refusal names what is unresolved about that one.
    for missing, expect in (("k_catch0", "NEVER chosen to lift"), ("x_catch_um", "load-strengthening"),
                            ("k_slip0", "provisional run value 0.35"), ("x_slip_um", "a transfer"),
                            ("k_on", "not a measurement"), ("capture_radius_um", "KINETIC radius")):
        try:
            CrossbridgeKinetics(**{**provisional, missing: None})
        except ValueError as exc:
            assert "PI-GAP" in str(exc) and expect in str(exc), (missing, exc)
        else:  # pragma: no cover
            raise AssertionError(f"{missing} must be refused, not defaulted")

    for bad, expect in ((dict(scope="  "), "scope is required"),
                        (dict(provenance=""), "provenance is required"),
                        (dict(kT=0.0), "kT must be finite")):
        try:
            CrossbridgeKinetics(**{**provisional, **bad})
        except ValueError as exc:
            assert expect in str(exc), (bad, exc)
        else:  # pragma: no cover
            raise AssertionError(f"CrossbridgeKinetics({bad}) must refuse")

    # A slip-dominated set is a Bell bond under a catch bond's name. Refused at CONSTRUCTION.
    try:
        CrossbridgeKinetics(**{**provisional, "k_catch0": 0.01, "k_slip0": 1.0})
    except ValueError as exc:
        assert "SLIP-DOMINATED" in str(exc) and "wrong sign" in str(exc)
    else:  # pragma: no cover
        raise AssertionError("a slip-dominated parameter set must not construct")

    cs = CrossbridgeKinetics(**provisional)

    # THE CLAIM, and it is a sign rather than a magnitude: the off-rate FALLS with load below F* and
    # RISES above it. A Bell slip is monotone increasing everywhere; this is not.
    f_star = cs.peak_force_pn
    assert f_star > 0.0
    f = np.linspace(0.0, 4.0 * f_star, 4001)
    p = np.asarray(cs.off_rate(f), np.float64)
    below, above = f < 0.98 * f_star, f > 1.02 * f_star
    assert np.all(np.diff(p[below]) < 0.0), "catch branch: detachment must SLOW under load"
    assert np.all(np.diff(p[above]) > 0.0), "slip branch: detachment must ACCELERATE past F*"
    assert abs(f[int(np.argmin(p))] - f_star) < 2.0 * (f[1] - f[0]), "the numeric minimum IS F*"
    assert np.isclose(float(cs.off_rate(0.0)), cs.unloaded_off_rate)
    # The magnitude convention: a crossbridge slips under |F|, so the law is even in the load.
    assert np.allclose(cs.off_rate(-f), p), "the off-rate must read the tension MAGNITUDE"

    # The engaged fraction is EMERGENT and is reported as a cross-check, never imposed.
    phi = cs.unloaded_bound_fraction
    assert 0.0 < phi < 1.0 and np.isclose(phi, 50.0 / (50.0 + 0.7))

    # The record carries every magnitude beside its scope, so a provisional cannot read as a datum.
    row = cs.record()
    assert "PROVISIONAL" in str(row["scope"]) and "PI-GAP" in str(row["provenance"])
    assert row["peak_force_pn"] == f_star

    print(f"crossbridge_kmc self-check OK — catch->slip crossover F* = {f_star:.4g} pN on a "
          f"PROVISIONAL set (no magnitude claimed); codegen {_codegen_check()} chars")


if __name__ == "__main__":
    _demo()
