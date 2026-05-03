"""V2 output package: HDF5 frame dumps and JSON metadata writers."""

from acs.v2.output.frame_dump import (
    FrameDumpReadResult,
    read_frame,
    write_frame,
)

__all__ = [
    "FrameDumpReadResult",
    "read_frame",
    "write_frame",
]
