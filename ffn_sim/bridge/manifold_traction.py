"""Per-contact-patch FA traction FIELD on the surface manifold (GEOMETRY-only bin).

This is the H.7 surface-manifold design's **first-class output** for an adherent
spread cell: it turns the explicit focal-adhesion integrin↔ligand bond forces
into a spatial **traction vector field**, one entry per cell-surface contact
patch, decomposed into the patch's local normal/tangential frame.

Division of labour (HARD — CLAUDE.md / AGENTS.md / the H.7 manifold design)
---------------------------------------------------------------------------
* The :class:`~ffn_sim.cortex.surface_manifold.SurfaceManifold` supplies
  **GEOMETRY ONLY**: the contact patches (triangles), their centroids, areas,
  outward normals ``n̂`` and tangent frames ``(e1, e2)``, and the bead→patch map
  (:meth:`SurfaceManifold.nearest_patch`). It owns **no mechanics** — no force,
  no tension, no γ, no spring. This module never asks it for any.
* The **FORCES** come entirely from the explicit FA integrin↔ligand bonds
  (bond type ``integrin_ligand``, :data:`BOND_TYPE_INTEGRIN`). Per the runtime
  ``md.bond.Harmonic`` + the Pereverzev updater (``bridge/integrin_bonds.py``),
  each engaged bond carries

      ``F_mag = k_int_bare · clip(|Δr| − integrin_r0, 0)``   along  û = (lig − int)/|·|

  i.e. the integrin→ligand bond unit vector. This module reads those bonds and
  **re-bins** their force vectors onto patches — it creates no force of its own.

What "traction" means here
--------------------------
The integrin end of each bond sits on the cell's basal cortex/membrane surface;
the ligand end is pinned on the rigid substrate plane (z=0). The bond force the
substrate exerts back on the cell (the traction the cell feels) is along
``û_int→lig``; equivalently the force the cell exerts on the substrate is its
negative. We report the **force the substrate ligand exerts on the cell**
(``+F_mag · û_int→lig``, integrin pulled toward its ligand) summed per patch,
then split into:

* **normal traction**   ``t_n = n̂ · ΣF``      (signed; + along outward n̂, i.e.
  pulling the surface OUT toward the substrate / detaching, − pressing IN),
* **tangential traction** ``t_t = |ΣF − t_n·n̂|`` (the in-plane shear magnitude —
  the load that resists cell spreading / migration).

Both are reported as a **force [N]** (the patch's summed bond force component)
AND a **traction [Pa]** (that force divided by the patch's geometric area).

Contact-conservation gate (HARD)
---------------------------------
The manifold only RE-BINS the explicit bond forces; it must neither create nor
lose force. :func:`measure_fa_traction_field` therefore asserts that the sum of
the per-patch accumulated force vectors equals the sum of all engaged integrin
bond force vectors to floating-point tolerance (``np.allclose``). A failure means
the binning lost/created force (a bug) — the caller must treat it as fatal, not
silently accept it. The returned dict carries the gate result so a driver can
exit nonzero.

Sanity Gate (per CLAUDE.md Sanity Gate Protocol)
------------------------------------------------
* **Dimensional.** ``k_int_bare`` [N/m] · ``clip(|Δr|−r₀,0)`` [m] → ``F_mag`` [N].
  ``ΣF`` [N]; patch area [m²]; traction = [N]/[m²] = [Pa]. ✓
* **Boundary.** Zero engaged bonds → all-zero fields, empty basal mask, gate
  trivially passes (0 == 0). A patch with no bonds carries exactly zero force.
* **Conservation.** Σ(per-patch force) ≡ Σ(per-bond force) — the contact gate.
  Newton's 3rd law is inherited from ``md.bond.Harmonic`` (the substrate side
  carries −F); we report the cell-side traction only, consistently.
* **Sign / sense.** ``F_mag = k·clip(|Δr|−r₀, 0) ≥ 0`` (a stretched bond pulls;
  a compressed one is clipped to zero, matching ``integrin_bonds.py:271``). The
  sign of the *normal* traction is physical (outward-pull + / inward-press −),
  carried by ``n̂·ΣF``; the tangential magnitude is unsigned by construction.
* **Numerical.** float64 throughout; positions read tag-ordered (a HOOMD
  snapshot is tag-ordered) exactly as
  ``cortex/cortical_tension.py:_read_tag_ordered_positions_and_bonds``.
* **Measurement-protocol consistency.** ``F_mag`` is recomputed from the same
  closed form the runtime updater applies (``bridge/integrin_bonds.py``), so the
  field is the engaged-clutch force with no baseline subtraction (construction
  r₀ is the realised separation ⇒ force-free baseline).

This module imports **no HOOMD force/integrator construction** and creates **no
spring/tension** — grep-clean of force creation beyond reading the explicit bond
geometry. ``numpy`` only (plus the manifold's geometry).
"""

from __future__ import annotations

from typing import Any

import numpy as np

from ffn_sim.cortex.surface_manifold import SurfaceManifold

__all__ = [
    "BOND_TYPE_INTEGRIN",
    "TYPE_INTEGRIN",
    "TYPE_LIGAND",
    "measure_fa_traction_field",
]

# Bond / particle type strings — the SINGLE source of truth is
# ``ffn_sim.bridge.fa`` (BOND_TYPE_INTEGRIN, TYPE_INTEGRIN, TYPE_LIGAND) and the
# cell-integration mirror in ``ffn_sim.cell.cell`` (FA_BOND_INTEGRIN_LIGAND,
# FA_TYPE_INTEGRIN, FA_TYPE_SUBSTRATE_LIGAND). They agree on these literals; we
# re-declare them here (rather than import the heavy fa/cell modules) so this
# pure-measurement module stays import-light. A startup assert in the driver can
# cross-check against fa.py if desired.
BOND_TYPE_INTEGRIN: str = "integrin_ligand"
TYPE_INTEGRIN: str = "integrin"
TYPE_LIGAND: str = "ligand"

# Number of azimuthal (φ) bins for the FA-polarity distribution. A geometry-free
# reporting resolution (histogram bin count), NOT a physics tuning knob — it does
# not enter any force or gate. 36 bins → 10° resolution.
_N_PHI_BINS_DEFAULT: int = 36


def _read_tag_ordered_positions_and_bonds(
    sim: Any,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, list[str]]:
    """Read tag-ordered positions + bond topology from a HOOMD simulation.

    Mirrors ``cortex/cortical_tension.py:_read_tag_ordered_positions_and_bonds``:
    a HOOMD snapshot's particle rows are tag-ordered, so a bond's endpoint tag
    indexes directly into the returned ``pos_by_tag`` array.

    Args:
        sim: A built HOOMD ``Simulation`` (its state holds particles + bonds).

    Returns:
        ``(pos_by_tag, bond_group, bond_typeid, bond_types)``:

        * ``pos_by_tag``: ``(N, 3)`` positions re-indexed into tag order [m].
        * ``bond_group``: ``(B, 2)`` per-bond endpoint tag indices.
        * ``bond_typeid``: ``(B,)`` per-bond type id.
        * ``bond_types``: bond-type name strings (type-id order).
    """
    bond_types = list(sim.state.bond_types)
    with sim.state.cpu_local_snapshot as s:
        pos = np.asarray(s.particles.position)
        tg = np.asarray(s.particles.tag)
        inv = np.empty_like(tg)
        inv[tg] = np.arange(tg.size)
        pos_by_tag = pos[inv].copy()
        # Particle type id per tag (to identify which bond end is the integrin).
        ptypeid = np.asarray(s.particles.typeid)
        ptypeid_by_tag = ptypeid[inv].copy()
        bond_group = np.asarray(s.bonds.group).copy()
        bond_typeid = np.asarray(s.bonds.typeid).copy()
    particle_types = list(sim.state.particle_types)
    return pos_by_tag, ptypeid_by_tag, particle_types, bond_group, bond_typeid, bond_types


def measure_fa_traction_field(
    sim: Any,
    manifold: SurfaceManifold,
    *,
    k_int_bare: float,
    integrin_r0: float,
    n_phi_bins: int = _N_PHI_BINS_DEFAULT,
    conservation_atol: float = 1.0e-18,
    conservation_rtol: float = 1.0e-6,
) -> dict[str, Any]:
    """Bin the explicit FA integrin-bond forces into a per-patch traction field.

    The manifold supplies geometry (patches + local frames + areas); the FA
    ``integrin_ligand`` bonds supply the forces
    (``F_mag = k_int_bare·clip(|Δr|−integrin_r0, 0)`` along the integrin→ligand
    unit vector — the same closed form the runtime
    :class:`~ffn_sim.bridge.integrin_bonds.IntegrinBondUpdater` applies). Each
    bond's INTEGRIN end is mapped to its home patch
    (:meth:`SurfaceManifold.nearest_patch`); per patch the bond force vectors are
    summed and decomposed into the patch's outward-normal and tangential
    components, reported as both force [N] and traction [Pa].

    Args:
        sim: A built, FA-adhered HOOMD ``Simulation`` (e.g. ``cell.simulation``).
        manifold: A :class:`SurfaceManifold` fitted to the cell surface (GEOMETRY
            ONLY — supplies patches/frames/areas, no physics).
        k_int_bare: Integrin bond stiffness [N/m] (``cell.p_fa.k_int_bare``).
        integrin_r0: Integrin bond rest length [m] (``cell.p_fa.integrin_r0``).
        n_phi_bins: Number of azimuthal (φ) histogram bins for the FA-polarity
            distribution (reporting resolution only).
        conservation_atol: Absolute tolerance [N] for the contact-conservation
            gate (Σ per-patch force == Σ per-bond force).
        conservation_rtol: Relative tolerance for the same gate.

    Returns:
        A dict with the per-patch fields and gate result::

            {
              # --- per-patch arrays (length n_tri; non-contact patches are 0) ---
              "patch_centroids": (n_tri, 3) [m],
              "patch_normal_force_N":    (n_tri,) signed [N]  (n̂·ΣF),
              "patch_tangential_force_N":(n_tri,) [N]         (|ΣF − (n̂·ΣF)n̂|),
              "patch_force_vec_N":       (n_tri, 3) [N]       (ΣF per patch),
              "patch_normal_traction_Pa":    (n_tri,) signed [Pa],
              "patch_tangential_traction_Pa":(n_tri,) [Pa],
              "patch_traction_mag_Pa":       (n_tri,) [Pa]    (|ΣF|/area),
              "patch_area_m2":               (n_tri,) [m²],
              # --- basal contact mask (patches carrying ≥1 integrin bond) ---
              "basal_patch_mask": (n_tri,) bool,
              "basal_patch_ids":  (n_basal,) int,
              # --- totals ---
              "total_normal_force_N": float,      # Σ n̂·ΣF over contact patches
              "total_tangential_force_N": float,  # Σ tangential force
              "total_force_vec_N": (3,) [N],      # vector sum of all bond forces
              "total_normal_traction_Pa": float,  # contact-area-weighted
              "total_tangential_traction_Pa": float,
              "contact_area_m2": float,           # Σ area of basal patches
              # --- azimuthal (φ) distribution of |traction| (FA polarity) ---
              "phi_bin_edges_rad": (n_phi_bins+1,),
              "phi_traction_mag_Pa": (n_phi_bins,) [Pa],  # |force|/area per φ wedge
              "phi_force_mag_N": (n_phi_bins,) [N],
              "polarity_vector_xy": (2,),   # Σ (|F| · n̂_xy) — net in-plane heading
              "polarity_magnitude": float,  # |polarity_vector_xy| / Σ|F|  in [0,1]
              "polarity_angle_rad": float,  # atan2 of the polarity vector
              # --- bookkeeping / gate ---
              "n_integrin_bonds": int,        # engaged integrin_ligand bonds
              "n_bonds_loaded": int,          # of those, with F_mag > 0
              "n_basal_patches": int,
              "subdivisions_hint": str,       # patch scale vs FA footprint note
              "contact_conservation": {
                  "passed": bool,
                  "sum_patch_force_N": (3,) [N],
                  "sum_bond_force_N": (3,) [N],
                  "abs_residual_N": (3,) [N],
                  "max_abs_residual_N": float,
                  "atol": float, "rtol": float,
              },
            }
    """
    if not (np.isfinite(k_int_bare) and k_int_bare > 0.0):
        raise ValueError(f"k_int_bare must be finite > 0; got {k_int_bare!r}")
    if not (np.isfinite(integrin_r0) and integrin_r0 >= 0.0):
        raise ValueError(f"integrin_r0 must be finite ≥ 0; got {integrin_r0!r}")
    if n_phi_bins < 1:
        raise ValueError(f"n_phi_bins must be ≥ 1; got {n_phi_bins}")

    n_tri = manifold.n_tri
    centroids = np.asarray(manifold.tri_centroids, dtype=np.float64)
    normals = np.asarray(manifold.tri_normals, dtype=np.float64)
    # Patch areas: prefer the manifold's own tri_areas (flat-triangle area), else
    # compute from the vertices. GEOMETRY ONLY — never a physical observable.
    areas = getattr(manifold, "tri_areas", None)
    if areas is None:
        v0 = manifold.verts[manifold.tris[:, 0]]
        v1 = manifold.verts[manifold.tris[:, 1]]
        v2 = manifold.verts[manifold.tris[:, 2]]
        areas = 0.5 * np.linalg.norm(np.cross(v1 - v0, v2 - v0), axis=1)
    areas = np.asarray(areas, dtype=np.float64)

    (
        pos_by_tag,
        ptypeid_by_tag,
        particle_types,
        bg,
        bt,
        bond_types,
    ) = _read_tag_ordered_positions_and_bonds(sim)

    # --- select integrin_ligand bonds ------------------------------------
    try:
        int_bond_typeid = bond_types.index(BOND_TYPE_INTEGRIN)
    except ValueError:
        raise RuntimeError(
            f"sim has no {BOND_TYPE_INTEGRIN!r} bond type; this is not an "
            f"FA-adhered cell. bond types present: {bond_types}."
        )
    is_int = bt == int_bond_typeid
    int_bonds = bg[is_int]
    n_int_bonds = int(int_bonds.shape[0])

    # Identify, per bond, which endpoint is the integrin and which the ligand.
    # The builder always inserts integrin_ligand bonds as (integrin_tag,
    # ligand_tag) (bridge/integrin_bonds.py act() col0=integrin, col1=ligand),
    # but we resolve it ROBUSTLY by particle type so a swapped column cannot
    # silently flip the bond direction.
    try:
        integrin_ptype = particle_types.index(TYPE_INTEGRIN)
        ligand_ptype = particle_types.index(TYPE_LIGAND)
    except ValueError:
        raise RuntimeError(
            f"sim must define both {TYPE_INTEGRIN!r} and {TYPE_LIGAND!r} "
            f"particle types; present: {particle_types}."
        )

    # Pre-allocate per-patch accumulators (length n_tri).
    patch_force_vec = np.zeros((n_tri, 3), dtype=np.float64)
    basal_mask = np.zeros(n_tri, dtype=bool)

    sum_bond_force = np.zeros(3, dtype=np.float64)
    n_loaded = 0
    # Per-bond arrays kept for the φ distribution + polarity (integrin-end frame).
    bond_force_mag_list: list[float] = []
    bond_integrin_pos_list: list[np.ndarray] = []

    if n_int_bonds > 0:
        col0 = int_bonds[:, 0]
        col1 = int_bonds[:, 1]
        t0 = ptypeid_by_tag[col0]
        # col0 is the integrin when its particle type is "integrin"; else col1.
        col0_is_int = t0 == integrin_ptype
        integrin_tag = np.where(col0_is_int, col0, col1)
        ligand_tag = np.where(col0_is_int, col1, col0)
        # Validate: every bond must be exactly one integrin + one ligand.
        end_ok = (
            (ptypeid_by_tag[integrin_tag] == integrin_ptype)
            & (ptypeid_by_tag[ligand_tag] == ligand_ptype)
        )
        if not np.all(end_ok):
            n_bad = int(np.count_nonzero(~end_ok))
            raise RuntimeError(
                f"{n_bad}/{n_int_bonds} {BOND_TYPE_INTEGRIN!r} bonds are not a "
                "clean (integrin, ligand) pair — the FA bond topology is "
                "inconsistent; refusing to bin an ambiguous force direction."
            )

        r_int = pos_by_tag[integrin_tag]   # (B, 3)
        r_lig = pos_by_tag[ligand_tag]     # (B, 3)
        dr = r_lig - r_int                 # integrin → ligand
        r = np.linalg.norm(dr, axis=1)
        r_safe = np.where(r > 0.0, r, 1.0)
        u_hat = dr / r_safe[:, None]       # integrin→ligand unit vector
        # SAME closed form as bridge/integrin_bonds.py:271 — F_mag is the bond
        # tension; the force the ligand exerts on the cell points integrin→ligand.
        F_mag = k_int_bare * np.clip(r - integrin_r0, 0.0, None)  # (B,) ≥ 0 [N]
        F_vec = F_mag[:, None] * u_hat                            # (B, 3) [N]

        sum_bond_force = F_vec.sum(axis=0)
        n_loaded = int(np.count_nonzero(F_mag > 0.0))

        # Map each bond's INTEGRIN-end position → home patch (geometry only).
        home_patch = manifold.nearest_patch(r_int)  # (B,) int

        # Accumulate ΣF per patch (np.add.at handles repeated patch ids).
        np.add.at(patch_force_vec, home_patch, F_vec)
        basal_mask[home_patch] = True

        bond_force_mag_list = list(np.linalg.norm(F_vec, axis=1))
        bond_integrin_pos_list = list(r_int)

    # --- per-patch normal / tangential decomposition ----------------------
    # t_n = n̂·ΣF (signed); tangential vector = ΣF − t_n·n̂; t_t = |tangential|.
    patch_normal_force = np.einsum("ij,ij->i", patch_force_vec, normals)  # (n_tri,)
    tang_vec = patch_force_vec - patch_normal_force[:, None] * normals
    patch_tangential_force = np.linalg.norm(tang_vec, axis=1)
    patch_force_mag = np.linalg.norm(patch_force_vec, axis=1)

    area_safe = np.where(areas > 0.0, areas, np.inf)  # 0-area → 0 traction
    patch_normal_traction = patch_normal_force / area_safe
    patch_tangential_traction = patch_tangential_force / area_safe
    patch_traction_mag = patch_force_mag / area_safe

    basal_ids = np.flatnonzero(basal_mask)
    n_basal = int(basal_ids.size)
    contact_area = float(np.sum(areas[basal_mask]))

    total_normal_force = float(np.sum(patch_normal_force[basal_mask]))
    total_tangential_force = float(np.sum(patch_tangential_force[basal_mask]))
    contact_area_safe = contact_area if contact_area > 0.0 else np.inf
    total_normal_traction = total_normal_force / contact_area_safe
    total_tangential_traction = total_tangential_force / contact_area_safe

    # --- CONTACT-CONSERVATION GATE ---------------------------------------
    # The manifold only RE-BINS the explicit bond forces. Σ over patches of the
    # accumulated per-patch force MUST equal Σ over all integrin bonds of F_vec.
    sum_patch_force = patch_force_vec.sum(axis=0)
    residual = sum_patch_force - sum_bond_force
    max_abs_residual = float(np.max(np.abs(residual))) if residual.size else 0.0
    passed = bool(
        np.allclose(
            sum_patch_force,
            sum_bond_force,
            atol=conservation_atol,
            rtol=conservation_rtol,
        )
    )

    # --- azimuthal (φ) distribution of |traction| (FA polarity) -----------
    # φ taken about the cell's z-axis (the substrate-normal axis) at each bond's
    # integrin end. This is the basal-plane heading of the adhesion load — the
    # design's first-class FA-polarity readout.
    phi_edges = np.linspace(-np.pi, np.pi, n_phi_bins + 1)
    phi_force_mag = np.zeros(n_phi_bins, dtype=np.float64)
    phi_area = np.zeros(n_phi_bins, dtype=np.float64)
    polarity_vec_xy = np.zeros(2, dtype=np.float64)
    sum_force_mag_total = 0.0
    if bond_force_mag_list:
        r_int_arr = np.asarray(bond_integrin_pos_list, dtype=np.float64)
        fmag_arr = np.asarray(bond_force_mag_list, dtype=np.float64)
        phi = np.arctan2(r_int_arr[:, 1], r_int_arr[:, 0])  # (B,)
        bin_idx = np.clip(
            np.digitize(phi, phi_edges) - 1, 0, n_phi_bins - 1
        )
        np.add.at(phi_force_mag, bin_idx, fmag_arr)
        # Per-bin contact area: area of the home patch each bond fell on, summed
        # (so a wedge's |traction| = Σ|F| / Σ(home-patch area) is dimensionally Pa).
        home_patch = manifold.nearest_patch(r_int_arr)
        np.add.at(phi_area, bin_idx, areas[home_patch])
        # In-plane polarity heading: |F|-weighted unit-xy of each integrin.
        rxy = r_int_arr[:, :2]
        rxy_norm = np.linalg.norm(rxy, axis=1)
        rxy_safe = np.where(rxy_norm > 0.0, rxy_norm, 1.0)
        rxy_hat = rxy / rxy_safe[:, None]
        polarity_vec_xy = np.sum(fmag_arr[:, None] * rxy_hat, axis=0)
        sum_force_mag_total = float(np.sum(fmag_arr))
    phi_area_safe = np.where(phi_area > 0.0, phi_area, np.inf)
    phi_traction_mag = phi_force_mag / phi_area_safe
    polarity_magnitude = (
        float(np.linalg.norm(polarity_vec_xy) / sum_force_mag_total)
        if sum_force_mag_total > 0.0
        else 0.0
    )
    polarity_angle = float(np.arctan2(polarity_vec_xy[1], polarity_vec_xy[0]))

    # --- resolution-vs-FA-footprint note (no magic number) ----------------
    # The mean patch edge length is the manifold's geometric scale; the FA
    # footprint is ~√(A_mature) (the design's FA patch scale). We REPORT how they
    # compare so a reader can judge whether the icosphere resolves the contact
    # (resolution-sensitivity hook), but nothing here is tuned to it.
    mean_edge = float(getattr(manifold, "mean_edge_length", 0.0))
    subdivisions_hint = (
        f"mean_patch_edge={mean_edge:.3e} m; n_tri={n_tri}; "
        f"basal_patches={n_basal}; bonds/basal_patch="
        f"{(n_loaded / n_basal):.2f}" if n_basal else
        f"mean_patch_edge={mean_edge:.3e} m; n_tri={n_tri}; no basal patches"
    )

    return {
        "patch_centroids": centroids,
        "patch_normal_force_N": patch_normal_force,
        "patch_tangential_force_N": patch_tangential_force,
        "patch_force_vec_N": patch_force_vec,
        "patch_normal_traction_Pa": patch_normal_traction,
        "patch_tangential_traction_Pa": patch_tangential_traction,
        "patch_traction_mag_Pa": patch_traction_mag,
        "patch_area_m2": areas,
        "basal_patch_mask": basal_mask,
        "basal_patch_ids": basal_ids,
        "total_normal_force_N": total_normal_force,
        "total_tangential_force_N": total_tangential_force,
        "total_force_vec_N": sum_bond_force,
        "total_normal_traction_Pa": float(total_normal_traction),
        "total_tangential_traction_Pa": float(total_tangential_traction),
        "contact_area_m2": contact_area,
        "phi_bin_edges_rad": phi_edges,
        "phi_traction_mag_Pa": phi_traction_mag,
        "phi_force_mag_N": phi_force_mag,
        "polarity_vector_xy": polarity_vec_xy,
        "polarity_magnitude": polarity_magnitude,
        "polarity_angle_rad": polarity_angle,
        "n_integrin_bonds": n_int_bonds,
        "n_bonds_loaded": n_loaded,
        "n_basal_patches": n_basal,
        "subdivisions_hint": subdivisions_hint,
        "contact_conservation": {
            "passed": passed,
            "sum_patch_force_N": sum_patch_force,
            "sum_bond_force_N": sum_bond_force,
            "abs_residual_N": np.abs(residual),
            "max_abs_residual_N": max_abs_residual,
            "atol": float(conservation_atol),
            "rtol": float(conservation_rtol),
        },
    }
