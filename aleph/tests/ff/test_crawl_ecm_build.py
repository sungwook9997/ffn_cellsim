r"""S4 integration step-1 (build-side) coverage — `build(..., ecm=True)` constructs a collagen-I Mikado slab
under the cell and attaches the basal clutches to fiber nodes, additively (ecm=False ⇒ no ECM in S).
"""

import numpy as np

from aleph.laws.ecm_mikado import MikadoNetwork
from aleph.scripts.ff_crawl_on_substrate import build


def test_ecm_off_is_additive():
    """Default (ecm=False): no ECM is built — S carries None (the validated crawl is untouched)."""
    S = build(n_cortex_fil=60, ecm=False)
    assert S["ecm"] is None and S["ecm_node"] is None


def test_ecm_builds_slab_and_attaches_clutches():
    """ecm=True: a Mikado collagen slab is built and the basal clutches bind fiber nodes within capture. (A
    dense-enough matrix — a too-sparse ECM leaves the FA out of reach, which is physically correct: no matrix,
    no adhesion.)"""
    S = build(n_cortex_fil=120, ecm=True, ecm_fibers=500, ecm_capture=1.5)
    mk = S["ecm"]
    assert isinstance(mk, MikadoNetwork)
    assert (len(mk.net.fiber_offsets) - 1) == 500                # the requested collagen fiber count
    assert mk.net.pos.shape[0] > 500                             # discretised into many nodes
    en = S["ecm_node"]
    assert en.shape[0] == S["basal"].size                        # one entry per basal clutch
    bound = en[en >= 0]
    assert bound.size > 0                                        # clutches reached the (dense) collagen
    assert bound.min() >= 0 and bound.max() < mk.net.pos.shape[0]   # every bound index is a valid ECM node


def test_ecm_slab_sits_under_the_cell():
    """The collagen slab spans below the basal cap (z from a depth under z_sub up to just above it) so the
    ventral FA clutches can reach it; the far (bottom) face is Dirichlet-pinned as bulk matrix."""
    S = build(n_cortex_fil=60, ecm=True, ecm_fibers=120, ecm_depth=4.0)
    mk = S["ecm"]; z = mk.net.pos[:, 2]
    assert z.min() < S["z_sub"]                                  # slab extends below the cell's basal plane
    assert z.max() <= S["z_sub"] + 1.0                           # and not above the ventral cap
    assert mk.pinned.any()                                       # a pinned bulk-collagen boundary exists
