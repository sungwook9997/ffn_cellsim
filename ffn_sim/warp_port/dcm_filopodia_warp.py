"""Warp @wp.kernel device force kernels for explicit FILOPODIA (Phase C · B3).

A filopodium is a thin actin-bundle finger that protrudes from a cell-surface
node, extends outward at a polymerization velocity, and whose TIP probes and
adheres on contact (the host :mod:`dcm_filopodia_host` owns extension + adhesion
bookkeeping; these kernels apply the per-step forces of the formed tip
adhesions). The mechanism is faithful to the unified actin-architecture view
(KB-3.8 filopodium structure / KB-3.6 Brownian-ratchet protrusion): the cortex,
the lamellipodium, and the filopodium are ONE category — filament bundles + a
polymerizing tip — differing only in their weaving and where the tip grips.

Two adhesion classes, each reusing an existing committed geometry pattern so the
node-FACE / node-to-plane philosophy is shared across the engine:

  * **tip-vs-cell-FACE** (:func:`filopodia_tip_face_force_kernel`): the tip adheres
    on a face of ANOTHER cell at the closest point (Ericson ``closest_bary``, the
    SAME geometry as :func:`dcm_neighbor_warp.contact_grid_kernel`). The spring
    pulls the ORIGINATING cell's base node toward that face point and, by Newton-3,
    scatters the reaction onto the gripped face's three vertices by barycentric
    weights — node-FACE, never node-node. This is the protrusion-mediated
    cell-cell traction that draws two cells together.

  * **tip-vs-SUBSTRATE** (:func:`filopodia_tip_plane_force_kernel`): the tip grips
    the dish at a fixed ligand site (z = z0) and pulls the base node toward it,
    ``F = k_tip·(anchor − r_base)`` — the SAME node-to-plane clutch spring as
    :func:`dcm_neighbor_warp.ecm_clutch_force_kernel`. The anchor is the rigid
    dish (no reaction), so the filopodium transmits its pull as traction.

The base node is the cell-surface node the filopodium grew FROM (the membrane
attachment point); the tip itself is a host-managed advancing anchor point (like
the lamellipodium actin bead — NOT a BAOAB DOF), so the only integrated body is
the membrane node it pulls. ``force_cap`` bounds each spring so a single adhesion
cannot inject a non-physical impulse (matches the lamellipodium tether cap).

Sanity Gate
-----------
Dimensional analysis:
  * ``k_tip`` [N/m] · displacement [m] = force [N] — consistent. ``k_tip`` is a
    filopodial-bundle tip-complex stiffness; PROVISIONAL anchor in the host.
  * ``force_cap`` [N] bounds |F|; the cap branch rescales by ``cap/|F|`` (unitless),
    preserving direction → still [N]. Setting ``force_cap`` huge ⇒ no cap.

Boundary cases:
  * ``t >= n`` thread returns (over-launch guard, like cadherin/ecm kernels).
  * face kernel: a degenerate target face (zero area / coincident base & closest
    point, ``L ≤ 1e-30``) contributes no force (guarded), never a divide-by-zero.
  * plane kernel: anchor == base ⇒ zero displacement ⇒ zero force (no NaN).
  * a tip with ``face_id < 0`` / detached (host sets sentinel) is skipped, so a
    broken adhesion injects nothing until re-formed.

Sign-sense:
  * face: ``r_vec = cpa − p_base`` points FROM the base TOWARD the face point, so
    ``+k_tip·r_vec`` on the base is ATTRACTIVE (pulls the originating cell toward
    the neighbour's face). The Newton-3 reaction ``−bary·F`` on the face vertices
    pulls the neighbour toward the base — momentum-clean, the two cells approach.
  * plane: ``anchor − r_base`` points toward the gripped dish site → attractive
    traction, identical sign to the validated ECM clutch.
  * forces are ``wp.atomic_add``-accumulated into the shared ``force`` array (a node
    is shared by many faces' threads), composing with turgor/contact/cohesion.

SI units throughout (positions ~1e-6 m, forces ~1e-9 N).
"""

from __future__ import annotations

import warp as wp

from ffn_sim.warp_port.dcm_contact_warp import closest_bary

wp.init()


@wp.kernel
def filopodia_tip_face_force_kernel(
    base_idx: wp.array(dtype=wp.int32),      # (M,) originating cell's base (membrane) node
    face_id: wp.array(dtype=wp.int32),       # (M,) adhered target face (−1 = detached/skip)
    bary: wp.array(dtype=wp.vec3d),          # (M,) frozen barycentric of the tip's grip point
    n_tip: wp.int32,
    faces: wp.array(dtype=wp.int32, ndim=2), # (n_faces,3) node ids of each face
    pos: wp.array(dtype=wp.vec3d),
    k_tip: wp.float64, force_cap: wp.float64,
    force: wp.array(dtype=wp.vec3d),
):
    """Tip-vs-cell-FACE adhesion spring (node-FACE, Newton-3).

    For each formed tip adhesion ``t``: the grip point on the target face is the
    barycentric combination ``cpa = Σ bary[k]·faces_vertex[k]`` of the face's three
    CURRENT vertex positions (so the grip tracks the deforming neighbour membrane,
    exactly like :func:`dcm_neighbor_warp.cadherin_bond_force_kernel` tracks its
    endpoints). The originating base node is pulled toward ``cpa`` with
    ``F = k_tip·(cpa − r_base)`` (capped); the equal-and-opposite ``−F`` is scattered
    onto the face's three vertices by ``bary`` (Newton-3), so total momentum is
    conserved and the two cells are drawn together. ``bary`` was set at adhesion time
    by the host from ``closest_bary`` — re-derivable here but frozen to keep the grip
    a fixed material point on the neighbour membrane (a real tip-complex sticks to
    the membrane patch it bound, not to the moving closest point)."""
    t = wp.tid()
    if t >= n_tip:
        return
    fj = face_id[t]
    if fj < wp.int32(0):
        return
    bvec = bary[t]
    ia = faces[fj, 0]
    ib = faces[fj, 1]
    ic = faces[fj, 2]
    cpa = pos[ia] * bvec[0] + pos[ib] * bvec[1] + pos[ic] * bvec[2]
    pb = pos[base_idx[t]]
    r_vec = cpa - pb                            # base → grip point (attractive dir)
    L = wp.length(r_vec)
    if L <= wp.float64(1.0e-30):
        return
    fmag = k_tip * L
    if fmag > force_cap:
        fmag = force_cap
    fvec = r_vec * (fmag / L)                    # |fvec| = capped, dir = toward face
    wp.atomic_add(force, base_idx[t], fvec)      # pull originating cell toward face
    wp.atomic_add(force, ia, fvec * (-bvec[0]))  # Newton-3 reaction on the gripped face
    wp.atomic_add(force, ib, fvec * (-bvec[1]))
    wp.atomic_add(force, ic, fvec * (-bvec[2]))


@wp.kernel
def filopodia_tip_plane_force_kernel(
    base_idx: wp.array(dtype=wp.int32),      # (K,) base node of each substrate-gripped filopodium
    anchor: wp.array(dtype=wp.vec3d),        # (K,) fixed dish ligand site the tip gripped (z=z0)
    n_tip: wp.int32,
    k_tip: wp.float64, force_cap: wp.float64,
    pos: wp.array(dtype=wp.vec3d),
    force: wp.array(dtype=wp.vec3d),
):
    """Tip-vs-SUBSTRATE clutch spring (node-to-plane).

    A filopodium whose tip reached the dish grips a fixed ligand site ``anchor``
    (z = z0). The originating base node is pulled toward that site with
    ``F = k_tip·(anchor − r_base)`` (capped). The anchor is the RIGID dish (no
    reaction), so the filopodium transmits its pull to the substrate as traction —
    identical in form and sign to :func:`dcm_neighbor_warp.ecm_clutch_force_kernel`
    (the host breaks the grip at a documented distance/catch-slip rate). Over-launch
    guard ``t >= n_tip``; zero displacement ⇒ zero force (no NaN)."""
    t = wp.tid()
    if t >= n_tip:
        return
    i = base_idx[t]
    d = anchor[t] - pos[i]
    L = wp.length(d)
    if L <= wp.float64(1.0e-30):
        return
    fmag = k_tip * L
    if fmag > force_cap:
        fmag = force_cap
    wp.atomic_add(force, i, d * (fmag / L))
