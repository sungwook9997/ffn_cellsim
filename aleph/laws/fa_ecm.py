r"""FF focal-adhesion ↔ ECM-fiber coupling (Warp, µm·pN·s) — the clutch binds a LIVE collagen-I (Mikado) fiber
node instead of a FIXED substrate point, so traction REMODELS the matrix (Newton reaction on the fiber).

The crawl/adhesion driver anchors each FA clutch to a fixed dish point (`fa_clutch_warp.clutch_spring_kernel`,
one-sided — the reaction goes to the immovable substrate). The PI's actual matrix is collagen-I: a disordered
Mikado fiber network (`ecm_mikado`) whose nodes are live DOFs. This module is the S4 primitive: a TWO-sided
clutch spring between an actin end node (in the cell's position array) and a collagen fiber node (in the ECM's
position array). The same Hookean law `f = k_int·(L−rest)/L·(ECM−actin)` now applies **+f to the actin (traction
on the cell) and −f to the fiber node (the Newton reaction)** — so the cell pulls its collagen substrate toward
it, recruiting/aligning fibers (the emergent 1/r remodel validated for DCM⊗Kim). No lumped traction law: the
force is the fine-grained clutch spring, and matrix remodelling is the reaction, molecule by molecule.

The cell and the ECM are SEPARATE position/force arrays (distinct networks, each with its own bending +
crosslink + drag dynamics); the clutch is the only mechanical bridge. `attach_clutches_to_ecm` forms the
nascent bonds (each clutch grabs the nearest fiber node within a capture radius); `clutch_ecm_spring_kernel`
applies the coupled force each step. Both are additive — a cell with no ECM (`ecm_node<0` everywhere) is
untouched, and the fixed-substrate path is unchanged.
"""
from __future__ import annotations

import numpy as np
import warp as wp

from scipy.spatial import cKDTree


@wp.kernel
def clutch_ecm_spring_kernel(
    cell_pos: wp.array(dtype=wp.vec3d),          # cell node positions [µm]
    actin: wp.array(dtype=wp.int32),             # (M,) actin end-node index per clutch (into cell_pos)
    ecm_pos: wp.array(dtype=wp.vec3d),           # ECM (collagen) node positions [µm]
    ecm_node: wp.array(dtype=wp.int32),          # (M,) bound ECM node index per clutch (into ecm_pos); <0 = unbound
    k_int: wp.float64,                           # clutch stiffness [pN/µm]
    rest: wp.float64,                            # clutch rest length [µm]
    cell_force: wp.array(dtype=wp.vec3d),        # += traction on the cell
    ecm_force: wp.array(dtype=wp.vec3d),         # += Newton reaction on the collagen fiber (remodels the matrix)
):
    """Two-sided Hookean clutch spring actin↔collagen-node: f = k_int·(L−rest)/L·(ECM−actin). +f on the actin
    (traction), −f on the fiber node (reaction). Unbound clutch (ecm_node<0) exerts nothing on either side."""
    t = wp.tid()
    j = ecm_node[t]
    if j < 0:
        return
    a = actin[t]
    d = ecm_pos[j] - cell_pos[a]
    L = wp.length(d)
    if L > wp.float64(1e-12):
        f = (k_int * (L - rest) / L) * d
        wp.atomic_add(cell_force, a, f)          # traction on the cell (toward the fiber)
        wp.atomic_add(ecm_force, j, -f)          # Newton pair: pull the collagen fiber toward the cell


@wp.kernel
def clutch_ecm_slip_traction_kernel(
    cell_force: wp.array(dtype=wp.vec3d),        # += forward retrograde-flow traction on the cell actin node
    ecm_force: wp.array(dtype=wp.vec3d),         # += equal-opposite reaction on the LIVE collagen node (two-way)
    actin: wp.array(dtype=wp.int32),             # (M,) cell actin-end node index per clutch
    ecm_node: wp.array(dtype=wp.int32),          # (M,) bound collagen node index per clutch (<0 = unbound)
    slip: wp.array(dtype=wp.float64),            # (M,) accumulated retrograde slip s [µm]
    phat: wp.vec3d,                              # crawl polarity axis (unit)
    k_int: wp.float64,                           # clutch stiffness [pN/µm]
    f_star: wp.float64,                          # catch-slip release force [pN]: sustained per-clutch traction is capped here
):
    """Molecular-clutch RETROGRADE-FLOW traction on the LIVE collagen — the two-way twin of
    ``motility_warp.clutch_slip_traction_kernel`` (which assumes a FIXED dish anchor). A bound clutch grips actin
    flowing REARWARD; the receded material (slip ``s``, accumulated once per step by
    ``clutch_slip_accumulate_kernel``) makes the clutch pull the cell node FORWARD by ``k_int·s·phat`` — and
    because the collagen is LIVE (not an immovable dish), the equal-and-opposite reaction pushes the engaged
    collagen node REARWARD (``−k_int·s·phat``). Momentum-conserving: the pinned collagen (z_lo BC) resists, so the
    cell crawls forward against a matrix it can also displace ⇒ traction/propulsion vs matrix stiffness emerges
    (molecular-clutch stiffness-sensing) instead of a body force. Same staggered launch as the ``clutch_ecm_spring``
    pair: cell side in the force-eval (``ecm_force``=dummy), collagen side in the substep (``cell_force``=dummy).

    The per-clutch force is CAPPED at ``f_star`` (Bell-Evans catch-slip release, KB): a real clutch releases at ~F*,
    so its sustained traction cannot exceed it. Without the cap the discrete kmc release (every kmc_every steps)
    lets ``k·s`` overshoot F* between ticks; the cap makes the per-clutch load physiological (≤F*) at any cadence."""
    t = wp.tid()
    j = ecm_node[t]
    if j < 0:
        return
    a = actin[t]
    f = wp.min(k_int * slip[t], f_star)          # sustained per-clutch traction ≤ catch-slip release force F*
    wp.atomic_add(cell_force, a, wp.vec3d(f * phat[0], f * phat[1], f * phat[2]))      # cell pulled FORWARD (crawl)
    wp.atomic_add(ecm_force, j, wp.vec3d(-f * phat[0], -f * phat[1], -f * phat[2]))    # collagen reaction (two-way)


@wp.kernel
def mask_ecm_by_bound_kernel(base: wp.array(dtype=wp.int32), bound: wp.array(dtype=wp.int32),
                             out: wp.array(dtype=wp.int32)):
    """``out[t] = base[t]`` if the clutch is BOUND (``bound[t]==1``) else ``-1`` — so only ENGAGED clutches grip
    the collagen. Makes the ECM traction respect the catch-slip bound state (a released clutch stops pulling the
    matrix), and makes the traction-OFF control (all clutches unbound) truly detach → the remodel is traction-driven."""
    t = wp.tid()
    if bound[t] == wp.int32(1):
        out[t] = base[t]
    else:
        out[t] = wp.int32(-1)


def attach_clutches_to_ecm(actin_pos: np.ndarray, ecm_pos: np.ndarray, capture_um: float) -> np.ndarray:
    """Form nascent FA↔collagen bonds: each clutch (at ``actin_pos[i]``) binds the NEAREST ECM node within
    ``capture_um``; returns the per-clutch ECM node index (``-1`` = no fiber in reach). Host, at attach time
    (occasional), via a KD-tree — cheap relative to the per-step force kernel."""
    actin_pos = np.ascontiguousarray(actin_pos, np.float64)
    ecm_pos = np.ascontiguousarray(ecm_pos, np.float64)
    if ecm_pos.shape[0] == 0 or actin_pos.shape[0] == 0:
        return np.full(actin_pos.shape[0], -1, np.int64)
    dist, idx = cKDTree(ecm_pos).query(actin_pos, k=1)
    return np.where(dist <= float(capture_um), idx, -1).astype(np.int64)


def _inplane_unit(director) -> np.ndarray:
    """The in-plane (xy) unit direction of ``director`` (contact guidance is an in-plane fiber-axis effect)."""
    d = np.asarray(director, float).copy()
    d[2] = 0.0
    n = np.linalg.norm(d)
    return d / n if n > 1e-12 else np.array([1.0, 0.0, 0.0])


def nematic_order_2d(vectors: np.ndarray, director) -> float:
    """2D nematic order S = ⟨2cos²φ − 1⟩ of the in-plane parts of ``vectors`` about the in-plane ``director``.

    φ is the angle between each vector's in-plane projection and the director. S=1 all-parallel, S=0 isotropic,
    S=−1 all-perpendicular (a headless/axial order parameter, sign of the vector ignored)."""
    v = np.asarray(vectors, float)
    if v.ndim != 2 or v.shape[0] == 0:
        return 0.0
    dhat = _inplane_unit(director)
    ehat = np.array([-dhat[1], dhat[0], 0.0])                 # in-plane ⟂ (= ẑ × d̂)
    vp = v.copy(); vp[:, 2] = 0.0
    n = np.linalg.norm(vp, axis=1)
    keep = n > 1e-12
    if not keep.any():
        return 0.0
    cos = (vp[keep] @ dhat) / n[keep]                          # cosφ (in-plane)
    return float(np.mean(2.0 * cos * cos - 1.0))


def traction_director_anisotropy(actin_pos: np.ndarray, ecm_pos_bound: np.ndarray, k_int: float,
                                 rest: float, director) -> dict:
    """Contact-guidance readout (ROADMAP NEAR #6): decompose the engaged-clutch traction on the fiber director.

    For each ENGAGED clutch (actin end node ``actin_pos[i]`` bound to fiber node ``ecm_pos_bound[i]``) the
    traction VECTOR is the SAME isotropic Hookean clutch law used by :func:`clutch_ecm_spring_kernel`,
    ``f = k_int·max(L−rest,0)/L·(ecm − actin)`` — no directional prefactor is introduced; any anisotropy must
    EMERGE from the aligned microstructure resisting that isotropic pull. The vectors are projected onto the
    in-plane director d̂ (‖) and the in-plane perpendicular ê=ẑ×d̂ (⟂):

    - ``F_par``  = Σ |f·d̂|,   ``F_perp`` = Σ |f·ê|  (sum-abs)
    - ``A_F``    = F_par / F_perp  — the traction anisotropy ratio (the Ray-2017 "force anisotropy"). It is a
      RATIO of two projections of the SAME clutch population, so the seed-to-seed few-clutch noise that made the
      absolute traction curve noisy cancels in numerator/denominator.
    - ``A_F_rms`` uses RMS instead of sum-abs (robust to a few large clutches).
    - ``S_bound`` = 2D nematic order of the engaged-clutch traction directions about the director.

    Returns zeros / A_F=nan for an empty engaged set. Pure host numpy; call once at the end of a quasi-static
    run (not per step)."""
    a = np.ascontiguousarray(actin_pos, np.float64)
    e = np.ascontiguousarray(ecm_pos_bound, np.float64)
    out = dict(F_par_pN=0.0, F_perp_pN=0.0, A_F=float("nan"), A_F_rms=float("nan"),
               n_engaged=int(a.shape[0]), S_bound=0.0)
    if a.shape[0] == 0:
        return out
    d = e - a
    L = np.linalg.norm(d, axis=1)
    ext = np.maximum(L - float(rest), 0.0)
    ok = L > 1e-12
    f = np.zeros_like(d)
    f[ok] = (float(k_int) * ext[ok] / L[ok])[:, None] * d[ok]  # per-clutch traction vector [pN]
    dhat = _inplane_unit(director)
    ehat = np.array([-dhat[1], dhat[0], 0.0])                  # in-plane ⟂
    fpar = f @ dhat
    fperp = f @ ehat
    F_par = float(np.abs(fpar).sum())
    F_perp = float(np.abs(fperp).sum())
    out["F_par_pN"] = F_par
    out["F_perp_pN"] = F_perp
    out["A_F"] = float(F_par / F_perp) if F_perp > 1e-12 else float("nan")
    rms_par = float(np.sqrt(np.mean(fpar ** 2)))
    rms_perp = float(np.sqrt(np.mean(fperp ** 2)))
    out["A_F_rms"] = float(rms_par / rms_perp) if rms_perp > 1e-12 else float("nan")
    out["S_bound"] = nematic_order_2d(f, director)
    return out


__all__ = ["clutch_ecm_spring_kernel", "mask_ecm_by_bound_kernel", "attach_clutches_to_ecm",
           "traction_director_anisotropy", "nematic_order_2d"]
