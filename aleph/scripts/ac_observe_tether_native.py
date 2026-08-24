#!/usr/bin/env python3
"""Wire :mod:`aleph.observe.tether` to the membrane the engine actually builds.

Why this driver exists
----------------------
``aleph/observe/tether.py`` is 682 lines with **zero call sites**
(``PARALLEL_AUDIT_PARAMS_DETECTORS_2026-08-20.md`` §1). It is the forward operator for a
membrane-tether pull, and the audit flags it as the most urgent of the five unwired detectors
because the throughput budget is sized from an AFM protocol.

What the wiring found, and it changes the framing
-------------------------------------------------
1. **The forward map is not missing from the runtime — it is present without its protocol.**
   ``aleph.laws.membrane_surface.erm_rupture_force`` is the SAME closed form the operator
   implements, ``f = 2*pi*sqrt(2*kappa*(gamma_mem + gamma_MCA))``, and
   ``components/incumbent/compartments.py:644`` evaluates it once as a scalar and hands it to
   ``erm_tether_kernel`` as a **rupture threshold**. So the number exists on every native run. It
   exists with no manifest, no error bar, no identifiability statement and no record — which is
   precisely the gap ``aleph/observe`` was written to close.

2. **The AFM harness that IS wired is an indenter, not a tether.**
   ``outputs/ac/afm_cortical_tension_readiness/`` implements a *spherical indentation* protocol.
   ``tether.py`` reads a *bead-pull tether plateau*. They are different experiments and this driver
   does not pretend otherwise: no indentation forward operator exists anywhere in
   ``aleph/observe`` (checked by grep, 2026-08-20). The audit's phrase "the operator that would
   read that instrument" holds for the membrane-cortex composite, not for the indenter.

3. **On today's runtime the tether force is an echo of three constants.** ``gamma_mem`` is the
   constant Raucher-Sheetz plateau (``compartments.py:579``, the ``K_A`` upturn is default-OFF),
   ``kappa`` is a resolved constant and ``GAMMA_MCA_PN_UM`` is a literal. A per-frame series built
   from them has zero scatter, so the operator reports ``standard_error = None`` and this driver
   stamps ``INPUT_ECHO_NOT_A_MEASUREMENT``. That verdict is the point of wiring it: the detector
   was always going to say this, and until it was called nothing did.

   The one input that can make the force move is the membrane's own apparent AREA, through
   ``aleph.engine.membrane_area.area_tension`` — itself an uncalled module today. This driver calls
   it, so a run whose area varies produces a varying tether force and a real error bar, with no
   edit here.

What it reads
-------------
The **membrane** population, and only that: its per-frame total apparent area [um^2], which is what
``membrane_area_reduce_kernel`` already sums on device over ``faces_d``. In arena terms that is the
FACE block of the membrane strand-set. Everything else on the record is a resolved constant.

Host-side by construction. Reducing an accepted trajectory to an observable belongs BETWEEN accepted
steps (``CLAUDE.md`` §GPU-only inside the loop); nothing here evaluates a law or takes a step, so
nothing here is a physics measurement and no number it prints may be quoted as one.

Usage:
    ac_observe_tether_native.py --declaration decl.json --out record.json

The declaration carries every protocol field and every state series. **Nothing is defaulted.** A
missing field is a refusal with the field named, because a defaulted protocol field is exactly the
rumour ``aleph/observe/manifest.py`` opens by describing.

Implements: KU-3.B1.4
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any

import numpy as np

from aleph.engine.membrane_area import MembraneAreaCard, area_tension
from aleph.laws.membrane_surface import GAMMA_MCA_PN_UM
from aleph.observe.manifest import ProtocolManifest
from aleph.observe.operator import ArrayState, ObservationContext, Refusal, observe
from aleph.observe.tether import AfmTetherForceOperator, afm_tether_manifest

SCHEMA = "ac.observe.tether/forward-operator@1"

#: Protocol fields the declaration must carry. These are exactly the arguments
#: :func:`~aleph.observe.tether.afm_tether_manifest` refuses to default, plus the four evidence
#: fields. The list is written out rather than derived from the signature so that adding a field
#: upstream fails here loudly instead of being silently dropped.
REQUIRED_PROTOCOL_FIELDS: tuple[str, ...] = (
    "temperature",
    "medium",
    "probe_geometry",
    "control_mode",
    "rate",
    "sampling_cadence",
    "fit_window",
    "resolution_limit",
    "specimen_identity",
    "batch",
    "citation_audit",
    "applicability_statement",
    "evidence_role",
    "operator_version",
)

#: State the declaration must carry alongside the protocol.
REQUIRED_STATE_FIELDS: tuple[str, ...] = (
    "membrane_area_um2",            # (T,) per-frame apparent area, from the accepted steps
    "rest_area_um2",                # A0
    "reservoir_capacity",           # MembraneAreaCard.capacity -- required, deliberately
    "reservoir_capacity_provenance",
    "gamma_mem_pn_per_um",
    "k_area_pn_per_um",
    "tau_lysis_pn_per_um",
    "bending_rigidity_pn_um",
    "sample_interval_s",
    "temperature_k",
    "retraction_speed_um_per_s",
)


def record_membrane_area_um2(membrane_position_d: Any, faces_d: Any, device: str) -> float:
    """Sum the membrane's total apparent area on device and read the one scalar back [um^2].

    For the live loop: call this **between accepted physical steps**, never inside the inner solve.
    It is one 8-byte D2H per accepted frame, which is what the residency gate permits there and
    forbids inside it.

    ``warp`` is imported here rather than at module scope only to keep the kernel symbol out of the
    reduction path. It is NOT an isolation claim: this module imports
    ``aleph.laws.membrane_surface`` for ``GAMMA_MCA_PN_UM``, which pulls ``warp`` in anyway. The
    constant is imported rather than repeated because one literal with two homes is how the two
    copies drift.

    Args:
        membrane_position_d: Device ``vec3d`` array of membrane vertex positions.
        faces_d: Device ``int32`` ``(n_faces, 3)`` array of outward-wound triangles.
        device: The CUDA device the arrays live on.

    Returns:
        Total membrane area [um^2].
    """
    import warp as wp

    from aleph.laws.membrane_surface import membrane_area_reduce_kernel

    area_d = wp.zeros(1, dtype=wp.float64, device=device)
    wp.launch(
        membrane_area_reduce_kernel,
        dim=int(faces_d.shape[0]),
        inputs=[membrane_position_d, faces_d],
        outputs=[area_d],
        device=device,
    )
    return float(area_d.numpy()[0])


def _missing(declaration: dict[str, Any], keys: tuple[str, ...]) -> list[str]:
    """Names of ``keys`` absent from ``declaration`` or present as null."""
    return [k for k in keys if declaration.get(k) is None]


def tether_from_declaration(declaration: dict[str, Any]) -> dict[str, Any]:
    """Run the tether operator over a declared membrane trajectory and return the record body.

    Args:
        declaration: The protocol and state fields, all required. See the module docstring.

    Returns:
        A JSON-able record. Its ``verdict`` is one of ``REFUSED`` (the operator declined and said
        why), ``INPUT_ECHO_NOT_A_MEASUREMENT`` (the force had no scatter, so the number restates
        its inputs) or ``MEASURED`` (the force varied and carries a correlated error bar).

    Raises:
        ValueError: If a required field is missing, or if the assembled manifest is incomplete.
    """
    gaps = _missing(declaration, REQUIRED_PROTOCOL_FIELDS + REQUIRED_STATE_FIELDS)
    if gaps:
        raise ValueError(
            "declaration is missing required field(s): "
            + ", ".join(gaps)
            + ". None of these is defaulted: a protocol field filled in on the caller's behalf is "
            "an assumption that becomes invisible one step downstream, and a state field filled in "
            "is a physics constant chosen by this script"
        )

    manifest: ProtocolManifest = afm_tether_manifest(
        **{k: declaration[k] for k in REQUIRED_PROTOCOL_FIELDS},
        **{
            k: declaration[k]
            for k in (
                "probe_location", "probe_direction", "adhesion_state", "dwell", "loading_history",
                "preconditioning", "calibration", "osmotic_state", "substrate", "confinement",
                "pharmacology",
            )
            if declaration.get(k) is not None
        },
    )
    # `manifest` wired, not merely built: nothing else in the tree calls `is_complete()` or
    # `not_recorded_fields()`, so until now the thirty-four fields were a type nobody was forced to
    # satisfy. It is NOT a hard gate, deliberately -- `is_complete()` is false while ANY field is
    # NOT_RECORDED, and refusing to write on that would push a caller to assert NOT_APPLICABLE
    # instead, which is the one substitution the two-member enum exists to prevent. So the
    # confession rides on the record and downgrades the verdict below.
    manifest_complete = manifest.is_complete()

    card = MembraneAreaCard(
        gamma_mem_pn_per_um=float(declaration["gamma_mem_pn_per_um"]),
        k_area_pn_per_um=float(declaration["k_area_pn_per_um"]),
        tau_lysis_pn_per_um=float(declaration["tau_lysis_pn_per_um"]),
        capacity=float(declaration["reservoir_capacity"]),
        capacity_provenance=str(declaration["reservoir_capacity_provenance"]),
    )
    areas = np.asarray(declaration["membrane_area_um2"], dtype=np.float64)
    rest_area = float(declaration["rest_area_um2"])

    # sigma_bilayer(t) from the membrane's own area. This is the only per-frame input; kappa and W
    # are constants and are stamped as such below.
    sigma_bilayer = np.array(
        [area_tension(float(a), rest_area, card) for a in areas], dtype=np.float64
    )
    kappa = np.full(areas.shape, float(declaration["bending_rigidity_pn_um"]), dtype=np.float64)
    # W: the membrane-cortex attachment energy. A literal in `laws/membrane_surface.py`, NOT a
    # quantity the ERM population produces -- `erm_density_per_um2` is ABSENT_FROM_CONTRACT_GRAPH
    # (audit §3). Overridable so a sourced value or a per-frame series can replace it without an
    # edit; the record says which was used.
    attachment_declared = declaration.get("membrane_attachment_energy_pn_per_um")
    attachment = np.full(areas.shape, GAMMA_MCA_PN_UM, dtype=np.float64)
    if attachment_declared is not None:
        attachment = np.broadcast_to(
            np.asarray(attachment_declared, dtype=np.float64), areas.shape
        ).astype(np.float64)

    state = ArrayState(
        arrays={
            "membrane": {
                "tension_pn_per_um": sigma_bilayer,
                "bending_rigidity_pn_um": kappa,
            },
            "cortex": {"membrane_attachment_energy_pn_per_um": attachment},
        }
    )
    context = ObservationContext(
        sample_interval_s=float(declaration["sample_interval_s"]),
        temperature_k=float(declaration["temperature_k"]),
        extras={"retraction_speed_um_per_s": float(declaration["retraction_speed_um_per_s"])},
    )
    operator = AfmTetherForceOperator(manifest=manifest)
    result = observe(operator, state, context)

    body: dict[str, Any] = {
        "schema": SCHEMA,
        "kind": "diagnostic",
        "quantitative_claim": "BLOCKED",
        "note": (
            "host-side reduction of an accepted trajectory; no law was evaluated and no step was "
            "taken here, so no number below is a physics measurement on its own"
        ),
        "manifest_hash": manifest.manifest_hash(),
        "manifest_complete": manifest_complete,
        "manifest_not_recorded_fields": list(manifest.not_recorded_fields()),
        "manifest_not_applicable_fields": list(manifest.not_applicable_fields()),
        "inputs": {
            "n_frames": int(areas.size),
            "membrane_area_um2_first": float(areas[0]) if areas.size else None,
            "membrane_area_um2_last": float(areas[-1]) if areas.size else None,
            "rest_area_um2": rest_area,
            "reservoir_capacity": card.capacity,
            "reservoir_capacity_provenance": card.capacity_provenance,
            "bilayer_tension_is_constant": bool(
                sigma_bilayer.size and np.all(sigma_bilayer == sigma_bilayer[0])
            ),
            "attachment_energy_source": (
                "declared" if attachment_declared is not None
                else "laws.membrane_surface.GAMMA_MCA_PN_UM literal (KU-3.B1.4)"
            ),
        },
    }

    if isinstance(result, Refusal):
        body["verdict"] = "REFUSED"
        body["refusal"] = result.as_dict()
        return body

    body["apparent_quantity"] = result.as_dict()
    if result.detail.get("constant_plateau"):
        body["verdict"] = "INPUT_ECHO_NOT_A_MEASUREMENT"
    elif manifest_complete:
        body["verdict"] = "MEASURED"
    else:
        # A number under a protocol with confessed holes is still a number; it is not comparable to
        # an external record, because `ExternalDatasetManifest.protocol_manifest_digest` joins on a
        # manifest hash and these two manifests do not describe the same protocol.
        body["verdict"] = "MEASURED_UNDER_INCOMPLETE_PROTOCOL"
    return body


def main(argv: list[str] | None = None) -> int:
    """Read a declaration, run the operator, write the record. Returns the process exit code."""
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--declaration", required=True, type=Path, help="JSON protocol + state")
    parser.add_argument("--out", required=True, type=Path, help="where to write the record")
    args = parser.parse_args(argv)

    declaration = json.loads(args.declaration.read_text())
    body = tether_from_declaration(declaration)
    body["declaration_source"] = str(args.declaration)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(body, indent=2) + "\n")
    print(f"{body['verdict']}  ->  {args.out}")
    # A refusal is a scientific result, not a crash: exit 0 and let the record carry the finding.
    # An echo is not a failure either. Only an unreadable declaration raises, above.
    return 0


if __name__ == "__main__":
    sys.exit(main())
