"""V2 visualization package: dev-grade stubs that consume HDF5 frame dumps.

The stub renderers in this package are explicitly *not* the
confocal-quality production path. The Blender confocal pipeline is a
separate design unit and reads the same HDF5 frames produced by
:mod:`acs.v2.output.frame_dump`. This package only provides fast PNG
and lightweight HTML output suitable for review/debug on the dev host.
"""

from acs.v2.viz.stub3d import (
    render_frame_html,
    render_frame_png,
)

__all__ = [
    "render_frame_html",
    "render_frame_png",
]
