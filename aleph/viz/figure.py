"""One command from a run artefact to a figure.

`ALEPH-PORT-3602`.

    python -m aleph.viz.figure runs/indentation/reference_level3.npz -o figure.svg

No hand-editing, no notebook, no manual step between an artefact on disk and a picture that can be
reviewed. That is the whole requirement, and the reason it is a requirement is that every manual step
between a run and its figure is a step where the figure can stop describing the run.

What counts as an artefact
--------------------------
Two shapes, both of which already exist in this repository:

``.npz``
    A bundle of named arrays. An ``(N, 3)`` float array is read as an owner's node positions; an
    ``(M, 2)`` integer array named ``<owner>_segments`` or ``<owner>_crosslinks`` is read as that
    owner's topology. **Anything else is listed as undrawn with the reason**, which is not a
    workaround -- it is the output. ``runs/indentation/reference_level3.npz`` carries a ``meta`` array
    of shape ``(5,)``, and the honest thing for a renderer to do with it is say so.

a capture directory
    Anything holding ``capture_manifest.json`` as written by this package's evidence capture. The
    manifest's census names the owners; geometry is read from an optional ``geometry/`` sidecar of
    ``<owner>_positions.npy`` and ``<owner>_segments.npy``. An owner named by the census with no
    geometry on disk is listed undrawn -- again, that is the answer and not a failure.

Refusal
-------
A refusal exits **2**, never 0. ``scripts/gpu_preflight.py`` set that convention here and the reason
carries over exactly: a tool that refuses and exits 0 lets ``make-figure && publish`` treat "nothing
was drawn" as success. Exit 1 is reserved for an unexpected fault so a caller can tell a refusal from
a crash.

Units: um throughout.
"""

from __future__ import annotations

import argparse
import json
import sys
from collections.abc import Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np

from aleph.viz.emergence import observe_emergence
from aleph.viz.render import DrawnScene, OwnerElements, draw_scene, render_scene_svg
from aleph.viz.scene import ComponentNode, IsolationState, SceneTree

__all__ = [
    "EXIT_OK",
    "EXIT_REFUSED",
    "EXIT_FAULT",
    "ArtefactError",
    "ReadArtefact",
    "figure_from_artefact",
    "main",
    "read_artefact",
]

EXIT_OK = 0
EXIT_FAULT = 1
EXIT_REFUSED = 2

#: Palette by owner name, so two filament populations never read as one colour. An owner not named
#: here gets a neutral stroke rather than a colour reused from another owner.
_OWNER_COLOUR: dict[str, str] = {
    "membrane": "#4fc3f7",
    "cortex": "#7e9cc4",
    "cortex_filaments": "#7e9cc4",
    "ecm": "#ffb74d",
    "lamellipodium": "#81c784",
    "filopodium": "#aed581",
    "sf_arc": "#f06292",
    "nmii": "#ba68c8",
    "microtubule": "#4db6ac",
    "intermediate_filament": "#9575cd",
    "nucleus": "#e57373",
    "focal_adhesion": "#fff176",
    "field": "#4fc3f7",
    "anchor": "#ffb74d",
}
_DEFAULT_COLOUR = "#90a4ae"


class ArtefactError(RuntimeError):
    """The artefact could not be read as one. Carries the reason the caller should see."""


@dataclass(slots=True)
class ReadArtefact:
    """An artefact resolved into owners, plus everything about it that could not be read.

    Attributes:
        label: Human identity of the artefact.
        owners: Per-owner arrays.
        unreadable: ``(name, reason)`` for every named thing that is not an owner's geometry.
        accepted_step: Accepted step, when the artefact records one.
        topology_epoch: Topology epoch, when the artefact records one.
    """

    label: str
    owners: dict[str, OwnerElements] = field(default_factory=dict)
    unreadable: list[tuple[str, str]] = field(default_factory=list)
    accepted_step: int = 0
    topology_epoch: int = 0


def _colour(name: str) -> str:
    return _OWNER_COLOUR.get(name, _DEFAULT_COLOUR)


def _read_npz(path: Path) -> ReadArtefact:
    """Read a bundle of named arrays into owners, naming everything that is not one."""
    artefact = ReadArtefact(label=path.name)
    with np.load(path, allow_pickle=False) as bundle:
        names = list(bundle.files)
        arrays = {name: bundle[name] for name in names}

    topology: dict[str, dict[str, np.ndarray]] = {}
    for name, array in arrays.items():
        for suffix, kind in (("_segments", "segments"), ("_crosslinks", "crosslinks")):
            if name.endswith(suffix):
                owner = name[: -len(suffix)]
                if array.ndim == 2 and array.shape[1] >= 2:
                    topology.setdefault(owner, {})[kind] = array[:, :2].astype(np.int64)
                else:
                    artefact.unreadable.append(
                        (name, f"shape {array.shape} is not (M, 2) index pairs")
                    )
                break

    consumed = {
        f"{owner}{suffix}"
        for owner in topology
        for suffix in ("_segments", "_crosslinks")
        if f"{owner}{suffix}" in arrays
    }
    for name, array in arrays.items():
        if name in consumed:
            continue
        if array.ndim == 2 and array.shape[1] == 3 and array.shape[0] > 0:
            artefact.owners[name] = OwnerElements(
                name=name,
                positions=array.astype(np.float64),
                segments=topology.get(name, {}).get("segments"),
                crosslinks=topology.get(name, {}).get("crosslinks"),
                colour=_colour(name),
            )
        else:
            artefact.unreadable.append(
                (
                    name,
                    f"array of shape {array.shape} is not (N, 3) positions in um, so it is not an "
                    "owner's geometry",
                )
            )
    if not artefact.owners:
        raise ArtefactError(
            f"{path} holds no (N, 3) position array, so there is no owner geometry to draw. "
            f"Arrays present: {', '.join(f'{n}{arrays[n].shape}' for n in names) or 'none'}"
        )
    return artefact


def _read_capture(directory: Path) -> ReadArtefact:
    """Read a capture directory: the manifest names the owners, a sidecar may carry their geometry."""
    manifest_path = directory / "capture_manifest.json"
    try:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise ArtefactError(f"{manifest_path} could not be read as JSON: {error}") from error

    schema = str(manifest.get("schema", ""))
    if not schema.startswith("aleph-viz-capture@"):
        raise ArtefactError(
            f"{manifest_path} declares schema {schema!r}, which this reader does not know"
        )

    artefact = ReadArtefact(label=str(manifest.get("run_label", directory.name)))
    frames = manifest.get("frames") or []
    if frames:
        artefact.accepted_step = int(frames[-1].get("accepted_step", 0))
        artefact.topology_epoch = int(frames[-1].get("topology_epoch", 0))

    census = manifest.get("census") or {}
    owners = [str(name) for name in (census.get("owners") or [])]
    if not owners:
        raise ArtefactError(f"{manifest_path} names no owners in its census, so there is nothing to draw")

    sidecar = directory / "geometry"
    for owner in owners:
        positions_path = sidecar / f"{owner}_positions.npy"
        if not positions_path.is_file():
            artefact.unreadable.append(
                (
                    owner,
                    "named by the capture census but no geometry/"
                    f"{owner}_positions.npy is on disk, so its arrays could not be read",
                )
            )
            continue
        positions = np.load(positions_path, allow_pickle=False)
        segments_path = sidecar / f"{owner}_segments.npy"
        segments = (
            np.load(segments_path, allow_pickle=False) if segments_path.is_file() else None
        )
        artefact.owners[owner] = OwnerElements(
            name=owner,
            positions=positions,
            segments=segments,
            colour=_colour(owner),
        )
    return artefact


def read_artefact(path: Path | str) -> ReadArtefact:
    """Resolve an artefact path into owners, or refuse with a reason.

    Args:
        path: A ``.npz`` bundle, a capture directory, or a ``capture_manifest.json``.

    Returns:
        The :class:`ReadArtefact`.

    Raises:
        ArtefactError: If the path is not an artefact this reader knows. The message names what was
            found, because "unsupported" without the finding is not actionable.
    """
    target = Path(path)
    if not target.exists():
        raise ArtefactError(f"{target} does not exist")
    if target.is_dir():
        if (target / "capture_manifest.json").is_file():
            return _read_capture(target)
        raise ArtefactError(
            f"{target} is a directory with no capture_manifest.json in it, so it is not a capture"
        )
    if target.name == "capture_manifest.json":
        return _read_capture(target.parent)
    if target.suffix == ".npz":
        return _read_npz(target)
    raise ArtefactError(
        f"{target} has suffix {target.suffix!r}; this reader knows .npz bundles and capture "
        "directories holding capture_manifest.json"
    )


def _scene_for(artefact: ReadArtefact) -> SceneTree:
    """A scene tree over the artefact's owners.

    Every owner named by the artefact becomes a component, **including the ones whose geometry could
    not be read**. Leaving those out would make the renderer's undrawn list shorter and the picture
    look more complete than the artefact is, which is the failure this whole entry is against.
    """
    components: list[ComponentNode] = []
    for name, owner in sorted(artefact.owners.items()):
        positions = np.asarray(owner.positions)
        entities: list[str] = ["node"]
        counts: list[int] = [int(positions.shape[0])]
        if owner.segments is not None:
            entities.append("segment")
            counts.append(int(np.asarray(owner.segments).shape[0]))
        if owner.crosslinks is not None:
            entities.append("crosslink")
            counts.append(int(np.asarray(owner.crosslinks).shape[0]))
        components.append(
            ComponentNode(name=name, entities=tuple(entities), element_counts=tuple(counts))
        )
    for name, reason in artefact.unreadable:
        if any(node.name == name for node in components):
            continue
        components.append(
            ComponentNode(
                name=name,
                entities=(),
                element_counts=(),
                isolation=IsolationState.EXCLUDED,
                isolation_reason=reason,
            )
        )
    return SceneTree(
        accepted_step=artefact.accepted_step,
        topology_epoch=artefact.topology_epoch,
        components=tuple(components),
        connectors=(),
    )


def figure_from_artefact(
    path: Path | str,
    *,
    axes: str = "xy",
    measure_emergence: bool = True,
    draws: int = 2000,
) -> tuple[str, DrawnScene]:
    """Read an artefact and return ``(svg_text, scene)``.

    Args:
        path: The artefact.
        axes: Which world plane to project onto.
        measure_emergence: Whether to run the observers on every owner that has segments.
        draws: Monte-Carlo draws for the null band.

    Returns:
        The SVG text and the drawn scene, so a caller can assert on the scene rather than on pixels.
    """
    artefact = read_artefact(path)
    tree = _scene_for(artefact)

    emergence: dict[str, Any] = {}
    if measure_emergence:
        for name, owner in sorted(artefact.owners.items()):
            if owner.segments is None or owner.positions is None:
                continue
            try:
                emergence[name] = observe_emergence(
                    owner.positions, owner.segments, draws=draws, local_draws=500
                )
            except (ValueError, RuntimeError) as error:
                # An observer that cannot measure says so; it never guesses, and it never stops the
                # figure. The reason travels into the picture with the owner's name attached.
                artefact.unreadable.append((name, f"emergence not measured: {error}"))

    scene = draw_scene(tree, artefact.owners, emergence=emergence)
    for name, reason in artefact.unreadable:
        if name not in scene.undrawn_owners and name not in scene.drawn_owners:
            scene.undrawn.append((name, reason))
    svg = render_scene_svg(
        scene, axes=axes, title=f"Aleph — {artefact.label}"
    )
    return svg, scene


def main(argv: Sequence[str] | None = None) -> int:
    """The one command. Returns an exit code; never raises for an unreadable artefact."""
    parser = argparse.ArgumentParser(
        prog="python -m aleph.viz.figure",
        description="Render an Aleph run artefact to an SVG figure in one command.",
    )
    parser.add_argument("artefact", help=".npz bundle, or a directory holding capture_manifest.json")
    parser.add_argument("-o", "--out", default=None, help="output .svg path (default: alongside)")
    parser.add_argument(
        "--axes", default="xy", choices=("xy", "xz", "yz"), help="world plane to project onto"
    )
    parser.add_argument(
        "--no-emergence", action="store_true", help="skip the nematic / condensation observers"
    )
    parser.add_argument("--draws", type=int, default=2000, help="null-band Monte-Carlo draws")
    parser.add_argument("--json", action="store_true", help="print a JSON summary on stdout")
    args = parser.parse_args(argv)

    source = Path(args.artefact)
    try:
        svg, scene = figure_from_artefact(
            source,
            axes=args.axes,
            measure_emergence=not args.no_emergence,
            draws=args.draws,
        )
    except ArtefactError as error:
        print(f"refused: {error}", file=sys.stderr)
        return EXIT_REFUSED

    destination = Path(args.out) if args.out else source.with_suffix(".svg")
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(svg, encoding="utf-8")

    problems = scene.coverage_problems()
    summary = {
        "artefact": str(source),
        "figure": str(destination),
        "elements": len(scene.elements),
        "drawn_owners": sorted(scene.drawn_owners),
        "undrawn": [{"owner": name, "reason": reason} for name, reason in scene.undrawn],
        "coverage_problems": problems,
        "faults": list(scene.tree.faults),
        "emergence": {
            name: observation.as_json_obj()
            for name, observation in sorted(scene.emergence.items())
        },
        "quantitative_status": "BLOCKED",
        "quantitative_status_reason": (
            "a figure is a projection of accepted state selected for legibility; an order parameter "
            "computed from it is a statement about this artefact's geometry and not a measurement "
            "of a cell"
        ),
    }
    if args.json:
        print(json.dumps(summary, indent=2))
    else:
        print(f"wrote {destination} — {len(scene.elements)} elements, {len(scene.undrawn)} undrawn")
        for name, reason in scene.undrawn:
            print(f"  undrawn  {name}: {reason}")
        for name, observation in sorted(scene.emergence.items()):
            verdict = "ORDERED" if observation.is_ordered else "not ordered"
            print(
                f"  emergence {name}: S={observation.order:.4f} z={observation.z:.2f} {verdict}, "
                f"{observation.n_bundles} bundle(s)"
            )
        for problem in problems:
            print(f"  COVERAGE DEFECT  {problem}", file=sys.stderr)

    # A coverage defect means an owner was silently omitted, which is the one thing this tool exists
    # to prevent. It is a refusal, not a warning.
    return EXIT_REFUSED if problems else EXIT_OK


if __name__ == "__main__":  # pragma: no cover - exercised through `main` in the tests
    raise SystemExit(main())
