#!/usr/bin/env python
"""Render journal-quality frames of a GSD trajectory with fresnel.

`fresnel <https://fresnel.readthedocs.io>`_ is the Glotzer-Lab path tracer that
ships alongside HOOMD-blue; it renders beads as spheres and bonds as cylinders
with ambient-occlusion lighting, using the Embree CPU backend on macOS and the
OptiX GPU backend on CUDA hosts (the Linux A5000 later in the project). This CLI
reads a GSD produced by :mod:`ffn_sim.scripts.npz_to_gsd` (or by a live
``hoomd.write.GSD`` writer) and renders each frame, coloring beads by a chosen
per-particle scalar (``log/particles/<scalar>``) through a perceptually-uniform
matplotlib colormap, with a rendered colorbar legend, a um scale bar, and
ambient-occlusion path tracing. Output is a PNG sequence; multi-frame runs are
assembled into an mp4 via ffmpeg.

Visualization-integrity (CLAUDE.md):

* SI units: the scale bar is derived from the sigma->um anchor stored in the GSD
  log (``sigma_um``, written by ``npz_to_gsd.py`` when the source ``.npz`` has it;
  0.1 um/sigma project default, CLI-overridable via ``--sigma-um``). The colorbar
  spans the full data range (no truncation).
* signed scalars use a diverging colormap (``coolwarm``) centred at zero;
  non-negative scalars use a sequential one (``viridis``; any matplotlib map via
  ``--cmap``).
* the colorbar legend and scale bar are composited onto every frame so the
  rendered PNG is self-describing.

fresnel is installed into ``ffn_sim`` via ``conda install -c conda-forge fresnel``
(0.13.8 np2py313 + embree 4.4.1 CPU backend on this Mac). It is imported lazily;
if it is absent the script fails with an actionable message. ``gsd`` and
``matplotlib`` are hard deps already present in ``ffn_sim``.

Run::

    python ffn_sim/scripts/render_fresnel.py --in traj.gsd --scalar curvature \\
        --out frames_dir [--samples 64] [--single-frame K]
"""

from __future__ import annotations

import argparse
import logging
import shutil
import subprocess
from pathlib import Path

import gsd.hoomd
import matplotlib

matplotlib.use("Agg")  # headless: no DISPLAY on gbook / CI

import matplotlib.cm as mcm
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np

logger = logging.getLogger(__name__)

# Project sigma->um anchor (CLAUDE.md: 1 sigma = 0.1 um). Used only when the GSD
# carries no 'sigma_um' log constant and no --sigma-um override is given.
_DEFAULT_SIGMA_UM = 0.1
_SEQUENTIAL_DEFAULT = "viridis"
_DIVERGING_DEFAULT = "coolwarm"

# Per-scalar colorbar labels with units (kT bending energy; reduced tension).
_SCALAR_LABELS = {
    "curvature": "bending energy  (kT)",
    "tension": "tension  (kT/sigma)",
}


def _import_fresnel():
    """Import fresnel, raising an actionable error if unavailable.

    Returns:
        The imported ``fresnel`` module.

    Raises:
        ImportError: With install guidance if fresnel is not importable.
    """
    try:
        import fresnel  # noqa: F401
    except Exception as exc:  # pragma: no cover - environment dependent
        raise ImportError(
            "fresnel is not installed/importable in this environment. Install with "
            "`conda install -c conda-forge fresnel` (a py313 build, 0.13.8, exists "
            "and pulls the embree CPU backend; OptiX is used automatically on a "
            "CUDA host). The GSD itself (positions + bonds + log/particles/*) is "
            "renderer-agnostic and already complete, so OVITO/Blender remain "
            "downstream options if fresnel is unavailable."
        ) from exc
    return fresnel


def _read_gsd_frames(path: str) -> list:
    """Load every frame of a GSD trajectory into memory.

    Args:
        path: Path to the ``.gsd`` file.

    Returns:
        List of frames (H.3 renders are typically strided, so this is small).
    """
    with gsd.hoomd.open(name=path, mode="r") as traj:
        return [traj[i] for i in range(len(traj))]


def _frame_scalar(frame, scalar: str) -> np.ndarray | None:
    """Fetch a per-particle scalar from a frame's log, trying key spellings.

    Args:
        frame: A GSD frame.
        scalar: Field name, e.g. ``'curvature'`` or ``'tension'``.

    Returns:
        The ``(N,)`` array, or ``None`` if not present under any known key.
    """
    log = frame.log or {}
    for key in (f"particles/{scalar}", f"log/particles/{scalar}", scalar):
        if key in log:
            return np.asarray(log[key]).ravel()
    return None


def _gsd_sigma_um(frame, override: float | None) -> float:
    """Resolve um-per-sigma from a CLI override, then the GSD log, then default.

    Args:
        frame: A GSD frame (its ``log`` may carry a ``sigma_um`` constant).
        override: Explicit CLI value, or ``None``.

    Returns:
        The sigma->um conversion factor (um per sigma).
    """
    if override is not None:
        return float(override)
    log = frame.log or {}
    for key in ("sigma_um", "log/sigma_um"):
        if key in log:
            return float(np.asarray(log[key]).ravel()[0])
    logger.warning(
        "no 'sigma_um' in GSD log and no --sigma-um; using project default "
        "%.3g um/sigma",
        _DEFAULT_SIGMA_UM,
    )
    return _DEFAULT_SIGMA_UM


def _choose_cmap(values: np.ndarray, cmap_name: str | None) -> tuple[str, bool]:
    """Pick a colormap and whether the field is signed (diverging).

    Args:
        values: All scalar values across the frames being rendered.
        cmap_name: Explicit colormap override, or ``None`` to auto-select.

    Returns:
        ``(cmap_name, is_diverging)``.
    """
    finite = values[np.isfinite(values)]
    signed = bool(finite.size and finite.min() < 0.0 < finite.max())
    if cmap_name is not None:
        return cmap_name, signed
    return (_DIVERGING_DEFAULT if signed else _SEQUENTIAL_DEFAULT), signed


def _make_norm(values: np.ndarray, is_diverging: bool) -> mcolors.Normalize:
    """Build a color normalization over the full (untruncated) data range.

    Diverging fields are centred at 0 and made symmetric so the midpoint colour
    means "zero"; sequential fields span ``[min, max]``.

    Args:
        values: All scalar values across the frames.
        is_diverging: Whether to centre at zero symmetrically.

    Returns:
        A matplotlib ``Normalize``.
    """
    finite = values[np.isfinite(values)]
    if finite.size == 0:
        return mcolors.Normalize(0.0, 1.0)
    vmin, vmax = float(finite.min()), float(finite.max())
    if vmin == vmax:
        vmax = vmin + 1.0
    if is_diverging:
        m = max(abs(vmin), abs(vmax))
        return mcolors.Normalize(-m, m)
    return mcolors.Normalize(vmin, vmax)


def _nice_round(x: float) -> float:
    """Round ``x`` down to a 1/2/5 x 10^k 'nice' number for a scale bar."""
    if x <= 0:
        return 1.0
    exp = np.floor(np.log10(x))
    base = 10.0**exp
    for m in (5.0, 2.0, 1.0):
        if m * base <= x:
            return m * base
    return base


def _fmt_um(um: float) -> str:
    """Format a length in um with a sensible unit (nm below 1 um)."""
    if um < 1.0:
        return f"{um * 1000.0:.0f} nm"
    if um == int(um):
        return f"{int(um)} um"
    return f"{um:.1f} um"


def _composite_overlays(
    raw_png: Path,
    out_png: Path,
    *,
    norm: mcolors.Normalize,
    cmap_name: str,
    scalar_label: str,
    box_edge_um: float,
) -> None:
    """Composite a colorbar legend + um scale bar onto a rendered PNG.

    Loads the raw fresnel PNG, draws it on a headless matplotlib canvas with no
    axes, then overlays (a) a horizontal scale bar annotated in um sized to the
    box edge and (b) a colorbar keyed to the same normalization/colormap used to
    color the beads. SI throughout (viz-integrity).

    Args:
        raw_png: Raw fresnel render (RGBA PNG).
        out_png: Destination PNG with overlays.
        norm: The color normalization used for the beads.
        cmap_name: The colormap used for the beads.
        scalar_label: Colorbar label (includes units where known).
        box_edge_um: Simulation box edge length in um, used to size the scale bar.
    """
    img = plt.imread(str(raw_png))
    h, w = img.shape[0], img.shape[1]
    cbar_px = 150
    fig = plt.figure(figsize=((w + cbar_px) / 100.0, h / 100.0), dpi=100)
    img_frac = w / (w + cbar_px)
    ax = fig.add_axes([0.0, 0.0, img_frac, 1.0])
    ax.imshow(img)
    ax.set_xlim(0, w)
    ax.set_ylim(h, 0)
    ax.set_axis_off()

    # --- scale bar (um): ~20% of the box edge, rounded to a nice number ------
    target_um = box_edge_um * 0.2 if box_edge_um > 0 else 1.0
    bar_um = _nice_round(target_um)
    frac = float(np.clip((bar_um / box_edge_um) if box_edge_um > 0 else 0.2, 0.05, 0.6))
    x0, y0 = 0.06, 0.92
    ax.plot(
        [x0 * w, (x0 + frac) * w],
        [y0 * h, y0 * h],
        color="white",
        lw=4,
        solid_capstyle="butt",
    )
    ax.text(
        (x0 + frac / 2.0) * w,
        (y0 - 0.025) * h,
        _fmt_um(bar_um),
        color="white",
        ha="center",
        va="bottom",
        fontsize=12,
        fontweight="bold",
    )

    # --- colorbar legend (full data range) ----------------------------------
    cax = fig.add_axes([img_frac + 0.04, 0.12, 0.035, 0.76])
    sm = mcm.ScalarMappable(norm=norm, cmap=cmap_name)
    sm.set_array([])
    cb = fig.colorbar(sm, cax=cax)
    cb.set_label(scalar_label, fontsize=11)
    cb.ax.tick_params(labelsize=9)

    fig.savefig(str(out_png), dpi=100, facecolor="black")
    plt.close(fig)


def _save_buffer_png(buf, path: Path) -> None:
    """Persist a fresnel image buffer to a PNG (PIL if present, else matplotlib).

    fresnel returns a :class:`fresnel.util.ImageArray`, a thin wrapper whose
    ``np.asarray`` yields a 0-d object array; slicing it with ``[:]`` returns the
    real ``(H, W, 4)`` uint8 RGBA ndarray, which is what the PNG writers need.

    Args:
        buf: A fresnel image buffer (``fresnel.util.ImageArray`` or any object
            sliceable to an ``(H, W, 4)`` uint8 RGBA ndarray).
        path: Output PNG path.
    """
    arr = np.asarray(buf[:], dtype=np.uint8)  # ImageArray[:] -> RGBA ndarray
    try:
        import PIL.Image  # type: ignore

        PIL.Image.fromarray(arr, mode="RGBA").save(str(path))
    except Exception:
        plt.imsave(str(path), arr)


def _assemble_movie(frames_dir: Path, out_mp4: Path, fps: int = 24) -> bool:
    """Assemble ``frame_*.png`` in ``frames_dir`` into an mp4 via ffmpeg.

    Args:
        frames_dir: Directory containing ``frame_%05d.png``.
        out_mp4: Output mp4 path.
        fps: Frames per second.

    Returns:
        ``True`` if ffmpeg ran successfully, else ``False`` (with a warning).
    """
    ffmpeg = shutil.which("ffmpeg")
    if ffmpeg is None:
        logger.warning("ffmpeg not found on PATH; skipping mp4 assembly")
        return False
    cmd = [
        ffmpeg,
        "-y",
        "-framerate",
        str(fps),
        "-i",
        str(frames_dir / "frame_%05d.png"),
        "-c:v",
        "libx264",
        "-pix_fmt",
        "yuv420p",
        "-crf",
        "18",
        str(out_mp4),
    ]
    logger.info("assembling mp4: %s", " ".join(cmd))
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0:
        logger.warning("ffmpeg failed (%d): %s", res.returncode, res.stderr[-500:])
        return False
    logger.info("wrote %s", out_mp4)
    return True


def render_gsd(
    in_path: str,
    out_dir: str,
    scalar: str = "curvature",
    *,
    samples: int = 64,
    single_frame: int | None = None,
    cmap_name: str | None = None,
    bead_scale: float = 0.5,
    bond_radius_frac: float = 0.35,
    width: int = 900,
    sigma_um: float | None = None,
    use_path_tracer: bool = True,
    make_movie: bool = True,
) -> list[str]:
    """Render frames of a GSD trajectory, colored by a per-particle scalar.

    Args:
        in_path: Input ``.gsd``.
        out_dir: Output directory (created); receives ``frame_*.png`` and, for
            multi-frame renders, ``movie.mp4``.
        scalar: Per-particle field under ``log/particles/*`` to color by (default
            ``'curvature'``). If absent, beads are colored flat grey and a warning
            is logged.
        samples: Path-tracer light samples per pixel (AA / AO quality). Ignored by
            the preview tracer.
        single_frame: Render only this frame index instead of the whole
            trajectory. Negative indexes from the end.
        cmap_name: Override the colormap; default auto (``viridis`` sequential /
            ``coolwarm`` diverging).
        bead_scale: Sphere radius as a fraction of the GSD particle diameter (falls
            back to ``bead_scale`` sigma if no diameter is stored).
        bond_radius_frac: Cylinder radius as a fraction of the bead radius.
        width: Render width in pixels (square output).
        sigma_um: Override um-per-sigma for the scale bar; default reads the GSD
            log (or the project default 0.1).
        use_path_tracer: Use the path tracer (AO, soft shadows); ``False`` forces
            the fast preview tracer. The path tracer is attempted first and falls
            back to preview on any failure.
        make_movie: Assemble an mp4 from the PNG sequence (multi-frame only) if
            ffmpeg is available.

    Returns:
        List of written PNG paths (plus the mp4 path appended if produced).

    Raises:
        ImportError: If fresnel is unavailable.
        IndexError: If ``single_frame`` is out of range.
    """
    fresnel = _import_fresnel()

    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)

    frames = _read_gsd_frames(in_path)
    n_total = len(frames)
    if n_total == 0:
        logger.warning("GSD %s has no frames; nothing to render", in_path)
        return []

    if single_frame is not None:
        idx = single_frame if single_frame >= 0 else n_total + single_frame
        if not (0 <= idx < n_total):
            raise IndexError(f"--single-frame {single_frame} out of range [0,{n_total})")
        indices = [idx]
    else:
        indices = list(range(n_total))

    # --- global color normalization over the frames being rendered ----------
    all_vals: list[np.ndarray] = []
    have_scalar = True
    for i in indices:
        s = _frame_scalar(frames[i], scalar)
        if s is None:
            have_scalar = False
            break
        all_vals.append(s)

    if have_scalar and all_vals:
        stacked = np.concatenate(all_vals)
        cmap_resolved, is_div = _choose_cmap(stacked, cmap_name)
        norm = _make_norm(stacked, is_div)
        scalar_label = _SCALAR_LABELS.get(scalar, scalar)
    else:
        logger.warning(
            "scalar %r not found in GSD log; rendering flat grey beads", scalar
        )
        cmap_resolved = cmap_name or _SEQUENTIAL_DEFAULT
        norm = mcolors.Normalize(0.0, 1.0)
        scalar_label = scalar
    cmap = matplotlib.colormaps[cmap_resolved]

    # --- tracer -------------------------------------------------------------
    device = fresnel.Device()
    if use_path_tracer:
        tracer = fresnel.tracer.Path(device=device, w=width, h=width)
    else:
        tracer = fresnel.tracer.Preview(device=device, w=width, h=width)

    box_edge_sigma = float(np.asarray(frames[0].configuration.box)[0])
    sig_um = _gsd_sigma_um(frames[0], sigma_um)
    box_edge_um = box_edge_sigma * sig_um
    logger.info(
        "box edge = %.3g sigma = %.3g um (sigma_um=%.3g)",
        box_edge_sigma,
        box_edge_um,
        sig_um,
    )

    written: list[str] = []

    for n, i in enumerate(indices):
        frame = frames[i]
        pos = np.asarray(frame.particles.position, dtype=np.float32)
        n_part = pos.shape[0]

        diam = getattr(frame.particles, "diameter", None)
        if diam is not None and np.size(diam) == n_part and np.any(np.asarray(diam) > 0):
            radius = float(np.mean(np.asarray(diam))) * bead_scale
        else:
            radius = float(bead_scale)

        scene = fresnel.Scene(device=device)

        # beads as spheres, colored by the scalar
        sph = fresnel.geometry.Sphere(scene, N=n_part, radius=radius)
        sph.position[:] = pos
        if have_scalar:
            vals = _frame_scalar(frame, scalar)
            sph.color[:] = fresnel.color.linear(cmap(norm(vals))[:, :3])
        else:
            sph.color[:] = fresnel.color.linear([[0.6, 0.6, 0.6]] * n_part)
        sph.material = fresnel.material.Material(roughness=0.5, primitive_color_mix=1.0)

        # bonds as cylinders between bonded bead pairs
        bgroup = (
            np.asarray(frame.bonds.group) if frame.bonds.N else np.empty((0, 2), int)
        )
        if bgroup.shape[0] > 0:
            cyl = fresnel.geometry.Cylinder(scene, N=bgroup.shape[0])
            cyl.points[:] = np.stack([pos[bgroup[:, 0]], pos[bgroup[:, 1]]], axis=1)
            cyl.radius[:] = radius * bond_radius_frac
            if have_scalar:
                vals = _frame_scalar(frame, scalar)
                c0 = fresnel.color.linear(cmap(norm(vals[bgroup[:, 0]]))[:, :3])
                c1 = fresnel.color.linear(cmap(norm(vals[bgroup[:, 1]]))[:, :3])
                cyl.color[:] = np.stack([c0, c1], axis=1)
            else:
                cyl.color[:] = fresnel.color.linear([0.5, 0.5, 0.5])
            cyl.material = fresnel.material.Material(
                roughness=0.5, primitive_color_mix=1.0
            )

        # ambient-occlusion lighting + fitted orthographic camera
        scene.lights = fresnel.light.cloudy()
        scene.background_color = fresnel.color.linear([0.0, 0.0, 0.0])
        scene.camera = fresnel.camera.Orthographic.fit(scene, view="front", margin=0.12)

        raw = out / f"_raw_{n:05d}.png"
        try:
            buf = (
                tracer.sample(scene, samples=samples)
                if use_path_tracer
                else tracer.render(scene)
            )
        except Exception as exc:  # pragma: no cover - backend dependent
            logger.warning("path tracer failed (%s); falling back to preview", exc)
            buf = fresnel.tracer.Preview(device=device, w=width, h=width).render(scene)
        _save_buffer_png(buf, raw)

        final = out / f"frame_{n:05d}.png"
        _composite_overlays(
            raw,
            final,
            norm=norm,
            cmap_name=cmap_resolved,
            scalar_label=scalar_label,
            box_edge_um=box_edge_um,
        )
        raw.unlink(missing_ok=True)
        written.append(str(final))
        logger.info("rendered frame %d/%d -> %s", n + 1, len(indices), final.name)

    if make_movie and len(written) > 1:
        mp4 = out / "movie.mp4"
        if _assemble_movie(out, mp4):
            written.append(str(mp4))

    return written


def main(argv: list[str] | None = None) -> int:
    """CLI entry point."""
    p = argparse.ArgumentParser(
        description="Render a GSD trajectory with fresnel (spheres+cylinders, "
        "scalar-colored, colorbar + um scale bar)."
    )
    p.add_argument("--in", dest="in_path", required=True, help="input .gsd path")
    p.add_argument(
        "--scalar",
        default="curvature",
        help="per-particle field to color by (default: curvature)",
    )
    p.add_argument("--out", dest="out_dir", required=True, help="output directory")
    p.add_argument(
        "--samples", type=int, default=64, help="path-tracer samples (default 64)"
    )
    p.add_argument(
        "--single-frame",
        type=int,
        default=None,
        help="render only this frame index (negative = from end)",
    )
    p.add_argument("--cmap", default=None, help="matplotlib colormap override")
    p.add_argument(
        "--sigma-um",
        type=float,
        default=None,
        help="override um-per-sigma for the scale bar (default: GSD log or 0.1)",
    )
    p.add_argument(
        "--preview",
        action="store_true",
        help="use the fast preview tracer instead of the path tracer",
    )
    p.add_argument("--width", type=int, default=900, help="render width px (square)")
    p.add_argument("-v", "--verbose", action="store_true", help="DEBUG logging")
    args = p.parse_args(argv)

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    )
    written = render_gsd(
        args.in_path,
        args.out_dir,
        scalar=args.scalar,
        samples=args.samples,
        single_frame=args.single_frame,
        cmap_name=args.cmap,
        width=args.width,
        sigma_um=args.sigma_um,
        use_path_tracer=not args.preview,
    )
    logger.info("done: %d output file(s) in %s", len(written), args.out_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
