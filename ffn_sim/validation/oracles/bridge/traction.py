"""Traction reducer over the engaged clutches of an FA (KU-2.12).

For Phase 1 Unit 2.1 the FA's loading axis is 1-D (the +x direction
matches the actin retrograde flow direction). The traction vector is
the sum of engaged-clutch forces projected along that axis, exposed as a
``(2,)`` ndarray in Newtons so the 2-D extension at Unit 2.2 can fill in
the second component without breaking the API.

Sanity Gate
-----------
1. Dimensional analysis: ``clutch_forces`` in N; output in N.
2. Boundary cases: zero engaged clutches → ``(0, 0)``.
3. Conservation: traction equals the negative of the force the FA
   exerts back on the actin bundle (Newton 3rd law; not enforced here,
   just stated).
4. Numerical sanity: float64; pure NumPy.
5. Sign / sense: positive +x traction means the FA pulls the substrate
   along the actin flow direction.
6. Measurement protocol: the validation notebook plots per-clutch force
   histograms separately to confirm the 5–20 pN range (KU-2.12), so
   this reducer is just the integrator.
"""

from __future__ import annotations

import numpy as np

from ffn_sim.validation.oracles.bridge.types import FocalAdhesion


def compute_traction(fa: FocalAdhesion) -> np.ndarray:
    """Sum of engaged-clutch forces along the FA's loading axis (1-D, +x).

    Returns
    -------
    traction : np.ndarray, shape (2,)
        ``(F_x, 0.0)`` in Newtons. The y-component is reserved for the
        Unit 2.2 2-D extension.
    """
    if fa.clutches_engaged.any():
        Fx = float(fa.clutch_forces[fa.clutches_engaged].sum())
    else:
        Fx = 0.0
    return np.array([Fx, 0.0], dtype=np.float64)
