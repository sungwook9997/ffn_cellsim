r"""Per-head ``walk_dir`` fill from actin barbed-end polarity — the I3 power-stroke hand-off (host oracle) — I4.

The I3 motor (``ac.motor.hand`` / ``minifilament_warp``) advances a bound head's crossbridge attachment point
along a per-head unit vector ``walk_dir`` by the walked ``abscissa`` (the myosin power stroke). ``hand.py``
allocates ``walk_dir`` ZERO (= passive head, the OFF/regression default) and its INTEGRATION.md §1 makes it a
HARD downstream obligation: **I4-weave MUST overwrite ``state["walk_dir"][h]`` from the true actin polarity
whenever a head (re)binds** — otherwise a head walks along the minifilament axis rather than along the actin it
grabbed, and the directed contraction is wrong.

Myosin II walks toward the BARBED (+) end of actin. So for a head bound to actin node ``a`` on a filament whose
barbed end is node ``b``, ``walk_dir = unit(pos[b] - pos[a])``. For a bipolar minifilament straddling two
ANTI-PARALLEL actin filaments (barbed ends pointing outward, as in a stress-fiber sarcomere or a contractile
cortex patch), both heads walk toward their respective outward barbed ends -> the reaction pulls the two actin
filaments' pointed ends together -> NET CONTRACTION emerges (the sign the powerstroke oracle certifies). This
module fills ``walk_dir``; the contraction sign / self-limiting force-velocity is proven in
``ac.motor.powerstroke_analytic`` (referenced, not duplicated).

This is a pure GEOMETRY function of the actin the head is on (its barbed-end node) — a SEEDED physiological
polarity, not a tuned magnitude. A zero ``walk_dir`` stays passive (production must fill it; physiological-
baseline rule).

Sanity Gate (self-tested in tests/ac/weave/test_walk_dir_oracle.py):
  * direction: for a head bound to actin node a with barbed end b, walk_dir == unit(pos[b]-pos[a]); |walk_dir|=1.
  * passive default: an unbound head (anchor = -1) gets walk_dir = 0 (stays passive; no directed force).
  * reattach update: rebinding a head to a filament of OPPOSITE polarity flips walk_dir (the hand-off tracks
    the actual actin, not a frozen geometry).
  * bipolar contraction sign: two anti-parallel actin filaments with barbed ends outward -> the two heads'
    walk_dirs point outward, and the minifilament's net directed action is contractile (inward), matching
    ac.motor.powerstroke_analytic's engaged-head sign.
"""

from __future__ import annotations

import numpy as np
import numpy.typing as npt

__all__ = ["barbed_end_node", "walk_dir_from_polarity", "fill_walk_dir"]

_AXIS_EPS = 1.0e-12


def barbed_end_node(
    fiber_offsets: npt.NDArray[np.int64], polarity: npt.NDArray[np.int64]
) -> npt.NDArray[np.int64]:
    """Global node index of every fiber's BARBED (+) end, from its per-fiber polarity flag.

    Args:
        fiber_offsets: (F+1,) contiguous node-ownership offsets (fiber f = nodes [off[f], off[f+1])).
        polarity: (F,) barbed-end flag: +1 = barbed end at the LAST node, -1 = barbed end at the FIRST node.

    Returns:
        (F,) global node index of each fiber's barbed end.
    """
    fiber_offsets = np.asarray(fiber_offsets, np.int64)
    polarity = np.asarray(polarity, np.int64)
    n_fib = fiber_offsets.shape[0] - 1
    if polarity.shape[0] != n_fib:
        raise ValueError("polarity must have one entry per fiber")
    first = fiber_offsets[:-1]
    last = fiber_offsets[1:] - 1
    return np.where(polarity >= 0, last, first).astype(np.int64)


def walk_dir_from_polarity(
    anchor_node: int, barbed_node: int, pos: npt.NDArray[np.float64]
) -> npt.NDArray[np.float64]:
    """Unit ``walk_dir`` for a head bound at ``anchor_node`` walking toward the filament barbed end ``barbed_node``.

    Returns a zero vector (passive) if the head is at (or coincident with) the barbed end so no direction is
    defined, or if ``anchor_node < 0`` (unbound).
    """
    if anchor_node < 0:
        return np.zeros(3)
    v = np.asarray(pos, float)[barbed_node] - np.asarray(pos, float)[anchor_node]
    n = float(np.linalg.norm(v))
    if n < _AXIS_EPS:
        return np.zeros(3)
    return v / n


def fill_walk_dir(
    bound: npt.NDArray[np.bool_],
    anchor: npt.NDArray[np.int64],
    node_fiber: npt.NDArray[np.int64],
    barbed_node_of_fiber: npt.NDArray[np.int64],
    pos: npt.NDArray[np.float64],
) -> npt.NDArray[np.float64]:
    """Fill ``walk_dir`` for every head from the barbed-end polarity of the actin it is bound to.

    This is the host mirror of the on-device attach-preset overwrite (INTEGRATION.md §1): a bound head's
    ``walk_dir`` points from its actin anchor toward that filament's barbed end; a free head (bound = False)
    gets zero (passive).

    Args:
        bound: (H,) per-head bound flags.
        anchor: (H,) actin node each bound head grabbed (-1 if free).
        node_fiber: (N,) fiber id of every actin node (maps an anchor to its filament).
        barbed_node_of_fiber: (F,) global barbed-end node of each fiber (from :func:`barbed_end_node`).
        pos: (N, 3) node positions [um].

    Returns:
        (H, 3) per-head walk_dir; unit vectors for bound heads, zero for free heads.
    """
    bound = np.asarray(bound, bool)
    anchor = np.asarray(anchor, np.int64)
    node_fiber = np.asarray(node_fiber, np.int64)
    barbed_node_of_fiber = np.asarray(barbed_node_of_fiber, np.int64)
    pos = np.asarray(pos, float)
    out = np.zeros((bound.shape[0], 3))
    for h in range(bound.shape[0]):
        if not bound[h] or anchor[h] < 0:
            continue
        fib = int(node_fiber[anchor[h]])
        out[h] = walk_dir_from_polarity(int(anchor[h]), int(barbed_node_of_fiber[fib]), pos)
    return out
