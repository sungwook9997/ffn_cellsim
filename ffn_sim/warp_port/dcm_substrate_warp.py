"""Warp port of the DCM substrate forces (Phase C · G2).

Ports the two substrate-adhesion forces from ``cell.dcm_gpu_forces`` — the route-a
spreading drivers that replace the deleted ``SettlingForce`` body-force proxy:

1. ``DcmSubstrateForceGPU`` — the z-WELL: a capped-harmonic adhesive well at z=z0
   (``plane_well_forces``) plus a stiff non-capped rigid-dish floor. z-only::

       |dz| <= rng : F_z = -k_well·dz          (harmonic well, toward z0)
       dz  < -rng  : F_z = +k_well·rng          (capped steric push-up)
       dz  > +rng  : F_z = 0                     (out of adhesion reach)
       always      : F_z += k_floor·max(0, z0-z) (rigid dish, z<z0 only)

   with ``k_well = 2·W_cs/rng²`` (series-combined with ``k_sub`` if a soft substrate).
   One thread per node, own-row write — no atomics, machine-eps.

2. ``DcmSubstrateWettingGPU`` — the in-plane WETTING drive: the xy-gradient of the
   substrate adhesion energy ``U_adh = -W_cs_Jm2·A_contact`` over basal contact
   triangles (the conservative, cortex-balanced analogue of soap-film wetting — NOT
   the rejected body-force proxy). Per face::

       zc   = mean vertex z;  w = clip(1 - (zc-z0)/rng, 0, 1)   (contact engagement)
       S    = ½·signed xy-area;  coef = ½·W·w·sign(S)
       ∂A_xy/∂v0 = (y1-y2, x2-x1)  (cyclic for v1,v2)  → xy force, scatter to nodes
   then a per-node cap on the in-plane magnitude (BAOAB guard). z is left to the well.

Parity contract (mirrors the DCM turgor piece):
* well: single own-row kernel → machine-eps vs the committed HOOMD reference.
* wetting ``reduce='host'`` (GATED): the per-face force law runs in Warp, the
  node scatter is done in numpy (``np.add.at``) exactly as the reference → isolates
  the force law (machine-eps).  ``reduce='warp'`` (diagnostic): scatter via
  ``wp.atomic_add`` (order differs → reduction-order tolerance < 1e-8).

Reference = the committed HOOMD fixture ``fixtures/dcm_substrate_ref.npz`` (the REAL
forces run on CPU; see ``fixtures/generate_dcm_substrate_fixture.py``) — guard-rail 2.
Fixed mesh (faces static for the step); the remesh is the separate Phase-C risk piece.
"""

from __future__ import annotations

import numpy as np

import warp as wp

wp.init()


# ──────────────────────────────────────────────────────────────────────────
# 1. z-well (DcmSubstrateForceGPU = plane_well_forces + rigid-dish floor)
# ──────────────────────────────────────────────────────────────────────────
@wp.kernel
def dcm_substrate_well_kernel(
    pos: wp.array(dtype=wp.vec3d),
    z0: wp.float64,
    k_well: wp.float64,
    rng: wp.float64,
    k_floor: wp.float64,
    force: wp.array(dtype=wp.vec3d),     # (N,) out, own-row write (no atomics)
):
    i = wp.tid()
    z = pos[i][2]
    dz = z - z0
    fz = wp.float64(0.0)
    if wp.abs(dz) <= rng:
        fz = -k_well * dz                 # harmonic well toward z0
    elif dz < -rng:
        fz = k_well * rng                 # capped steric push-up
    # rigid-dish floor (z < z0 only), non-capped — added on top of the well term
    below = z0 - z
    if below > wp.float64(0.0):
        fz = fz + k_floor * below
    force[i] = wp.vec3d(wp.float64(0.0), wp.float64(0.0), fz)


def run_dcm_substrate_well_warp(
    *,
    pos: np.ndarray,            # (N, 3)
    z0: float,
    W_cs: float,
    adh_range: float,
    k_sub: float | None = None,
    k_floor: float = 1.0,
    device: str = "cpu",
) -> dict:
    """DCM substrate z-well force (N,3) in Warp, matching DcmSubstrateForceGPU."""
    pos = np.ascontiguousarray(pos, dtype=np.float64)
    N = pos.shape[0]
    k_well = 2.0 * W_cs / (adh_range ** 2)
    if k_sub is not None:
        k_well = (k_well * k_sub) / (k_well + k_sub)

    pos_d = wp.array(pos, dtype=wp.vec3d, device=device)
    force_d = wp.zeros(N, dtype=wp.vec3d, device=device)
    wp.launch(dcm_substrate_well_kernel, dim=N,
              inputs=[pos_d, wp.float64(z0), wp.float64(k_well), wp.float64(adh_range),
                      wp.float64(k_floor), force_d], device=device)
    wp.synchronize_device(device)
    return {"force": force_d.numpy().astype(np.float64), "k_well": float(k_well)}


# ──────────────────────────────────────────────────────────────────────────
# 2. in-plane wetting (DcmSubstrateWettingGPU)
# ──────────────────────────────────────────────────────────────────────────
@wp.func
def _wetting_coef(v0: wp.vec3d, v1: wp.vec3d, v2: wp.vec3d,
                  z0: wp.float64, W: wp.float64, rng: wp.float64) -> wp.float64:
    """Per-face scalar coef = ½·W·w·sign(S) (engagement w, signed xy-area S)."""
    zc = (v0[2] + v1[2] + v2[2]) / wp.float64(3.0)
    w = wp.float64(1.0) - (zc - z0) / rng
    if w < wp.float64(0.0):
        w = wp.float64(0.0)
    if w > wp.float64(1.0):
        w = wp.float64(1.0)
    S = wp.float64(0.5) * ((v1[0] - v0[0]) * (v2[1] - v0[1])
                           - (v1[1] - v0[1]) * (v2[0] - v0[0]))
    sgn = wp.float64(0.0)
    if S > wp.float64(0.0):
        sgn = wp.float64(1.0)
    if S < wp.float64(0.0):
        sgn = wp.float64(-1.0)
    return (wp.float64(0.5) * W) * w * sgn


@wp.kernel
def dcm_wetting_perface_kernel(
    pos: wp.array(dtype=wp.vec3d),
    faces: wp.array(dtype=wp.int32, ndim=2),
    z0: wp.float64, W: wp.float64, rng: wp.float64,
    fout: wp.array(dtype=wp.float64, ndim=2),    # (M, 6): f0x f0y f1x f1y f2x f2y
):
    """host-reduce path: emit per-face vertex xy-forces; numpy does the scatter."""
    fi = wp.tid()
    v0 = pos[faces[fi, 0]]
    v1 = pos[faces[fi, 1]]
    v2 = pos[faces[fi, 2]]
    c = _wetting_coef(v0, v1, v2, z0, W, rng)
    fout[fi, 0] = c * (v1[1] - v2[1]);  fout[fi, 1] = c * (v2[0] - v1[0])
    fout[fi, 2] = c * (v2[1] - v0[1]);  fout[fi, 3] = c * (v0[0] - v2[0])
    fout[fi, 4] = c * (v0[1] - v1[1]);  fout[fi, 5] = c * (v1[0] - v0[0])


@wp.kernel
def dcm_wetting_scatter_kernel(
    pos: wp.array(dtype=wp.vec3d),
    faces: wp.array(dtype=wp.int32, ndim=2),
    z0: wp.float64, W: wp.float64, rng: wp.float64,
    force: wp.array(dtype=wp.vec3d),             # (N,) out, atomic accumulate
):
    """warp-reduce path: compute per-face forces and atomic-scatter to the 3 nodes."""
    fi = wp.tid()
    i0 = faces[fi, 0]
    i1 = faces[fi, 1]
    i2 = faces[fi, 2]
    v0 = pos[i0]
    v1 = pos[i1]
    v2 = pos[i2]
    c = _wetting_coef(v0, v1, v2, z0, W, rng)
    z = wp.float64(0.0)
    wp.atomic_add(force, i0, wp.vec3d(c * (v1[1] - v2[1]), c * (v2[0] - v1[0]), z))
    wp.atomic_add(force, i1, wp.vec3d(c * (v2[1] - v0[1]), c * (v0[0] - v2[0]), z))
    wp.atomic_add(force, i2, wp.vec3d(c * (v0[1] - v1[1]), c * (v1[0] - v0[0]), z))


@wp.kernel
def dcm_wetting_scatter_integrin_kernel(
    pos: wp.array(dtype=wp.vec3d),
    faces: wp.array(dtype=wp.int32, ndim=2),
    fcell: wp.array(dtype=wp.int32),
    integrin: wp.array(dtype=wp.float64),        # (n_cells,) per-cell substrate gain
    z0: wp.float64, W: wp.float64, rng: wp.float64,
    force: wp.array(dtype=wp.vec3d),             # (N,) out, atomic accumulate
):
    """M3 junction-switch variant of :func:`dcm_wetting_scatter_kernel`: the per-face
    in-plane wetting (cell→substrate traction) is scaled by the owning cell's
    ``integrin`` gain (raised to ``integrin_strong_factor`` when the cell crowd-switches),
    the cadherin→integrin clutch's substrate-strengthening half."""
    fi = wp.tid()
    i0 = faces[fi, 0]
    i1 = faces[fi, 1]
    i2 = faces[fi, 2]
    v0 = pos[i0]
    v1 = pos[i1]
    v2 = pos[i2]
    g = integrin[fcell[fi]]
    c = _wetting_coef(v0, v1, v2, z0, W, rng) * g
    z = wp.float64(0.0)
    wp.atomic_add(force, i0, wp.vec3d(c * (v1[1] - v2[1]), c * (v2[0] - v1[0]), z))
    wp.atomic_add(force, i1, wp.vec3d(c * (v2[1] - v0[1]), c * (v0[0] - v2[0]), z))
    wp.atomic_add(force, i2, wp.vec3d(c * (v0[1] - v1[1]), c * (v1[0] - v0[0]), z))


# ──────────────────────────────────────────────────────────────────────────
# Device-loop variants (ACCUMULATE into a shared force buffer; for the
# device-resident hybrid loop, not the standalone parity above which overwrites).
# Force LAWS are identical to the parity-verified kernels — only the write changes
# (own-row add for the well; separate wetting buffer + per-node cap-then-add so the
# cap applies to the wetting contribution ONLY, never the accumulated total).
# ──────────────────────────────────────────────────────────────────────────
@wp.kernel
def dcm_substrate_well_accum_kernel(
    pos: wp.array(dtype=wp.vec3d),
    z0: wp.float64, k_well: wp.float64, rng: wp.float64, k_floor: wp.float64,
    force: wp.array(dtype=wp.vec3d),     # (N,) ACCUMULATE (own-row read-add-write)
):
    i = wp.tid()
    z = pos[i][2]
    dz = z - z0
    fz = wp.float64(0.0)
    if wp.abs(dz) <= rng:
        fz = -k_well * dz
    elif dz < -rng:
        fz = k_well * rng
    below = z0 - z
    if below > wp.float64(0.0):
        fz = fz + k_floor * below
    fo = force[i]
    force[i] = wp.vec3d(fo[0], fo[1], fo[2] + fz)


@wp.kernel
def dcm_wetting_cap_add_kernel(
    wbuf: wp.array(dtype=wp.vec3d),      # (N,) accumulated wetting force (xy)
    cap: wp.float64,
    force: wp.array(dtype=wp.vec3d),     # (N,) total force to add the capped wetting into
):
    i = wp.tid()
    w = wbuf[i]
    mag = wp.sqrt(w[0] * w[0] + w[1] * w[1])
    sx = w[0]
    sy = w[1]
    if mag > cap:
        s = cap / mag
        sx = w[0] * s
        sy = w[1] * s
    fo = force[i]
    force[i] = wp.vec3d(fo[0] + sx, fo[1] + sy, fo[2])


def _cap_xy(F: np.ndarray, cap: float) -> np.ndarray:
    """Per-node cap on the in-plane magnitude (exactly the reference's clip)."""
    mag = np.sqrt((F[:, :2] ** 2).sum(axis=1))
    scale = np.where(mag > cap, cap / np.where(mag > 0, mag, 1.0), 1.0)
    F[:, 0] *= scale
    F[:, 1] *= scale
    return F


def run_dcm_substrate_wetting_warp(
    *,
    pos: np.ndarray,            # (N, 3), node order (= tag order in the fixture)
    faces: np.ndarray,          # (M, 3) node-id triplets
    z0: float,
    W_cs_Jm2: float,
    adh_range: float,
    force_cap: float = 5.0e-8,
    reduce: str = "host",
    device: str = "cpu",
) -> dict:
    """DCM in-plane substrate wetting force (N,3) in Warp, matching DcmSubstrateWettingGPU.

    ``reduce='host'`` isolates the per-face force law (machine-eps); ``reduce='warp'``
    uses atomic scatter (reduction-order tol). The per-node cap is applied after.
    """
    pos = np.ascontiguousarray(pos, dtype=np.float64)
    faces = np.ascontiguousarray(faces, dtype=np.int32)
    N = pos.shape[0]
    M = faces.shape[0]
    pos_d = wp.array(pos, dtype=wp.vec3d, device=device)
    faces_d = wp.array(faces, dtype=wp.int32, device=device)

    if reduce == "host":
        fout_d = wp.zeros((M, 6), dtype=wp.float64, device=device)
        wp.launch(dcm_wetting_perface_kernel, dim=M,
                  inputs=[pos_d, faces_d, wp.float64(z0), wp.float64(W_cs_Jm2),
                          wp.float64(adh_range), fout_d], device=device)
        wp.synchronize_device(device)
        ff = fout_d.numpy()
        F = np.zeros((N, 3), dtype=np.float64)
        np.add.at(F[:, :2], faces[:, 0], ff[:, 0:2])
        np.add.at(F[:, :2], faces[:, 1], ff[:, 2:4])
        np.add.at(F[:, :2], faces[:, 2], ff[:, 4:6])
    elif reduce == "warp":
        force_d = wp.zeros(N, dtype=wp.vec3d, device=device)
        wp.launch(dcm_wetting_scatter_kernel, dim=M,
                  inputs=[pos_d, faces_d, wp.float64(z0), wp.float64(W_cs_Jm2),
                          wp.float64(adh_range), force_d], device=device)
        wp.synchronize_device(device)
        F = force_d.numpy().astype(np.float64)
    else:
        raise ValueError(f"reduce must be 'host' or 'warp', got {reduce!r}")

    F = _cap_xy(F, force_cap)
    return {"force": F}


# ──────────────────────────────────────────────────────────────────────────
# C5. ECM MECHANO-FEEDBACK  (ADDITIVE · DEFAULT-OFF)
# ──────────────────────────────────────────────────────────────────────────
# Two mechanisms layered on the rigid-dish substrate above, both off by default so
# that an existing rigid-dish run is byte-identical (E_sub=None and strain-stiffening
# disabled both skip the new code paths and fall through to the original force law):
#
#   (1) SUBSTRATE-STIFFNESS MECHANOSENSING — the well/clutch effective stiffness the
#       substrate presents scales physically with a substrate Young's modulus E_sub,
#       and the transmitted traction follows the Chan-Odde rigidity-sensing curve
#       (rises with stiffness up to a plateau). The E_sub→∞ limit recovers the rigid
#       dish exactly.
#   (2) STRAIN-STIFFENING HOOK — a phenomenological local substrate stiffening under
#       accumulated normal traction (collagen/ECM strain-stiffens under load).
#
# ── KB anchors (queried 2026-06-22; do NOT invent magic numbers) ────────────
# KB-2.4  Chan-Odde motor-clutch (ChanOdde2008_Science; Bangasser2013_BiophysJ):
#         v_actin = v_unloaded·(1 − F_total/(N_m·F_stall)); clutch F_i = k_int·Δx;
#         defaults N_total=50, k_on=1 s⁻¹, v_unloaded=100 nm/s, F_stall=2 pN, k_int=1 pN/nm.
# KB-2.8  Optimal stiffness for max traction (biphasic): k_sub* ≈ N_m·F_stall/(v_unloaded·τ_bond),
#         k_sub* ~ 1–10 kPa. ⚠ This project's own Phase-1 finding (in KB-2.8): the curve is
#         SATURATING not sharply peaked — peak prominence 0.34%, ⟨F⟩ asymptotes to
#         ~100 pN = N_m·F_stall for E ≳ 500 Pa. So we model a saturating (monotone) rise
#         to a plateau, the physically-faithful form for this parameter regime.
# KB-2.12 Traction magnitudes: per-clutch 5–20 pN; per-FA 1–10 nN (~5 nN); per-cell 10–100 nN.
# KB-1.5  Substrate E range: brain 0.1–1, fat/mammary ~1, muscle 10, cartilage 100 kPa, bone 10 GPa;
#         Phase-1 default E=5 kPa, ν=0.45.
# KB-1.4 / KB-1.V.2.5  Storm-MacKintosh collagen strain-stiffening (Storm2005_Nature;
#         Licup2015_PNAS; Jansen2018_BiophysJ): differential modulus K(γ) ~ G0/(1−γ/γc)²
#         diverges above a critical strain γc≈0.16; in the nonlinear regime K ~ σ (stiffness
#         rises ~linearly with stress, 10–100× above γc). KB-1.V.2.5 marks γc as an EMERGENT,
#         geometry-set quantity, NOT a free knob — here it parametrises the phenomenological
#         hook only (the emergent Mikado-network result lives in the ecm/ module, not here).

# KB-2.8 optimal-stiffness crossover, expressed as a spring constant (N/m). The Chan-Odde
# crossover k_sub* ~ 1–10 kPa is a *modulus*; converted to a clutch spring over the contact
# patch by the same E→k map used below (e_sub_to_k_sub) so the crossover and the substrate
# spring live in one consistent unit. We pin the modulus crossover and convert on demand.
_CHAN_ODDE_E_OPT_PA: float = 5.0e3       # KB-2.8 k_sub* midpoint ~1–10 kPa → 5 kPa
# Chan-Odde clutch defaults (KB-2.4) used to set the traction plateau scale.
_KINT_N_PER_M: float = 1.0e-3            # k_int = 1 pN/nm = 1e-3 N/m (KB-2.4)


def e_sub_to_k_sub(E_sub: float, contact_radius: float, poisson: float = 0.45) -> float:
    """Map a substrate Young's modulus to the spring constant it presents to a clutch/well.

    Hertzian-contact form for an elastic half-space indented over a circular patch of
    radius ``a``: the substrate behaves as a linear spring of stiffness

        k_sub = 2·E*·a ,   E* = E_sub / (1 − ν²)

    (the standard flat-punch / Hertz half-space compliance; E* is the plane-strain reduced
    modulus). This is the physical bridge from the KB-1.5 modulus axis to the N/m spring the
    existing ``run_dcm_substrate_well_warp`` already series-combines into ``k_well``. As
    ``E_sub → ∞`` the spring stiffens without bound, so the series combination recovers the
    rigid dish (k_well unchanged) exactly — the rigid-dish limit.

    Args:
        E_sub: substrate Young's modulus (Pa). KB-1.5 biological range 0.1 kPa–100 kPa.
        contact_radius: clutch/well contact patch radius ``a`` (m); the ×40 mesh contact scale.
        poisson: substrate Poisson ratio ν (KB-1.5 default 0.45).

    Returns:
        k_sub (N/m), the substrate spring constant.

    Sanity Gate:
        * Dimensional: [Pa]·[m] = [N/m²]·[m] = [N/m] ✓.
        * E_sub→∞ ⇒ k_sub→∞ ⇒ series ``(k_well·k_sub)/(k_well+k_sub)`` → k_well (rigid limit).
        * Monotone increasing in E_sub (a stiffer substrate is a stiffer spring).
    """
    e_star = E_sub / (1.0 - poisson * poisson)
    return 2.0 * e_star * contact_radius


def chan_odde_traction_factor(E_sub: float, E_opt: float = _CHAN_ODDE_E_OPT_PA) -> float:
    """Chan-Odde rigidity-sensing traction modifier T(E_sub)/T(∞) ∈ (0, 1].

    The transmitted traction rises with substrate stiffness up to a plateau (KB-2.8: in this
    project's regime the biphasic curve is SATURATING, peak prominence 0.34%, not sharply
    peaked). We model that saturating rise with a Michaelis-Menten / Hill-1 form anchored to
    the KB-2.8 crossover modulus ``E_opt`` (~5 kPa, the 1–10 kPa midpoint)::

        T(E)/T(∞) = E / (E + E_opt)

    so that:
        * E → 0   ⇒ factor → 0   (soft substrate, frictional-slippage, clutches never load),
        * E = E_opt ⇒ factor = ½ (half-maximal traction at the KB-2.8 crossover stiffness),
        * E → ∞   ⇒ factor → 1   (stiff/rigid dish, full traction plateau = the original
                                   force law) — recovers the current rigid-dish behaviour.

    This is the SATURATING reading of KB-2.8 (the project's own Phase-1 result), chosen over a
    sharply-peaked biphasic form deliberately; see the module docstring. The descending (very-
    stiff "load-and-fail") arm of the canonical Chan-Odde curve has 0.34% prominence in this
    regime and is intentionally omitted — flagged as a PROVISIONAL modelling choice.

    Args:
        E_sub: substrate Young's modulus (Pa).
        E_opt: half-max crossover modulus (Pa); KB-2.8 k_sub* ~1–10 kPa midpoint = 5 kPa.

    Returns:
        Dimensionless traction factor in (0, 1].

    Sanity Gate:
        * Dimensionless (ratio of Pa/Pa) ✓.
        * Monotone increasing in E_sub; →1 as E_sub→∞ (rigid limit); =½ at E_sub=E_opt.
    """
    if E_sub == float("inf"):
        return 1.0
    return E_sub / (E_sub + E_opt)


def run_dcm_substrate_well_mechano_warp(
    *,
    pos: np.ndarray,            # (N, 3)
    z0: float,
    W_cs: float,
    adh_range: float,
    k_floor: float = 1.0,
    k_sub: float | None = None,
    E_sub: float | None = None,
    contact_radius: float | None = None,
    poisson: float = 0.45,
    device: str = "cpu",
) -> dict:
    """z-well with C5 substrate-stiffness mechanosensing (ADDITIVE, DEFAULT-OFF).

    With ``E_sub=None`` (default) this is byte-identical to
    :func:`run_dcm_substrate_well_warp` — same k_well, same kernel, same output (the explicit
    ``k_sub`` legacy path is preserved untouched). With ``E_sub`` set, the Chan-Odde rigidity-
    sensing physics is engaged:

        1. the substrate spring ``k_sub`` is DERIVED from E_sub (``e_sub_to_k_sub``) and
           series-combined into k_well exactly as the legacy soft-substrate path does, so a
           stiffer substrate transmits a stiffer well;
        2. the well force is additionally scaled by the Chan-Odde saturating traction factor
           ``chan_odde_traction_factor(E_sub)`` (rises with stiffness toward a plateau).

    The E_sub→∞ limit gives k_sub→∞ (series → k_well) and factor→1, recovering the rigid dish.

    Args:
        pos, z0, W_cs, adh_range, k_floor, device: as :func:`run_dcm_substrate_well_warp`.
        k_sub: legacy explicit substrate spring (N/m); used only when E_sub is None.
        E_sub: substrate Young's modulus (Pa). None ⇒ mechanosensing OFF (rigid dish).
        contact_radius: clutch contact patch radius ``a`` (m); REQUIRED when E_sub is set.
            The ×40 mesh contact scale — default-derive from adh_range if not given.
        poisson: substrate Poisson ratio (KB-1.5 default 0.45).

    Returns:
        dict with ``force`` (N,3), ``k_well`` (effective), ``k_sub`` (derived or None),
        ``traction_factor`` (1.0 when off).

    Sanity Gate:
        * E_sub=None ⇒ identical to the legacy well (force unchanged) — off-path invariance.
        * E_sub→∞ ⇒ k_sub→∞, factor→1 ⇒ rigid-dish force recovered.
        * Monotone: well force magnitude increases with E_sub (via both k_sub and the factor).
    """
    pos = np.ascontiguousarray(pos, dtype=np.float64)
    N = pos.shape[0]
    k_well = 2.0 * W_cs / (adh_range ** 2)

    k_sub_used: float | None = k_sub
    factor = 1.0
    if E_sub is not None:
        a = contact_radius if contact_radius is not None else adh_range
        k_sub_used = e_sub_to_k_sub(E_sub, a, poisson=poisson)
        factor = chan_odde_traction_factor(E_sub)
    if k_sub_used is not None and np.isfinite(k_sub_used):
        # series spring; the k_sub→∞ (rigid) limit is k_well, handled by skipping the
        # combination when k_sub is non-finite (avoids inf/inf = nan).
        k_well = (k_well * k_sub_used) / (k_well + k_sub_used)
    k_well_eff = k_well * factor

    pos_d = wp.array(pos, dtype=wp.vec3d, device=device)
    force_d = wp.zeros(N, dtype=wp.vec3d, device=device)
    # The traction factor scales the *adhesive* well stiffness only; the rigid steric
    # floor (k_floor) is geometric, not a transmitted traction, so it is left unfactored.
    wp.launch(dcm_substrate_well_kernel, dim=N,
              inputs=[pos_d, wp.float64(z0), wp.float64(k_well_eff), wp.float64(adh_range),
                      wp.float64(k_floor), force_d], device=device)
    wp.synchronize_device(device)
    return {"force": force_d.numpy().astype(np.float64),
            "k_well": float(k_well_eff), "k_sub": (None if k_sub_used is None else float(k_sub_used)),
            "traction_factor": float(factor)}


# ── (2) Strain-stiffening hook ──────────────────────────────────────────────
@wp.kernel
def dcm_substrate_well_strainstiff_kernel(
    pos: wp.array(dtype=wp.vec3d),
    z0: wp.float64,
    k_well_base: wp.float64,
    k_well_local: wp.array(dtype=wp.float64),   # (N,) per-node stiffened well stiffness
    rng: wp.float64,
    k_floor: wp.float64,
    force: wp.array(dtype=wp.vec3d),            # (N,) out, own-row write (no atomics)
):
    """z-well with a PER-NODE adhesive stiffness ``k_well_local`` (strain-stiffening hook).

    Identical to :func:`dcm_substrate_well_kernel` except the harmonic-well spring constant
    is read per-node from ``k_well_local`` (≥ ``k_well_base``) instead of a single scalar, so
    a node whose local substrate has stiffened under accumulated load pulls harder. With
    ``k_well_local[i] == k_well_base`` everywhere the kernel is byte-identical to the scalar one.
    """
    i = wp.tid()
    z = pos[i][2]
    dz = z - z0
    kw = k_well_local[i]
    fz = wp.float64(0.0)
    if wp.abs(dz) <= rng:
        fz = -kw * dz
    elif dz < -rng:
        fz = kw * rng
    below = z0 - z
    if below > wp.float64(0.0):
        fz = fz + k_floor * below
    force[i] = wp.vec3d(wp.float64(0.0), wp.float64(0.0), fz)


class SubstrateStrainStiffening:
    """Host-side per-node substrate strain-stiffening accumulator (ADDITIVE, DEFAULT-OFF).

    Phenomenological collagen/ECM strain-stiffening (KB-1.4 / KB-1.V.2.5, Storm2005_Nature):
    the local substrate differential stiffness rises with the stress the cell imposes on it.
    The canonical Storm-MacKintosh result is ``K ~ σ`` in the nonlinear regime (stiffness ~
    linear in stress, 10–100× above the critical strain). We implement that as a per-node
    multiplier on the adhesive well stiffness driven by the time-integrated downward (−z)
    force the node has exerted on the dish::

        I_i      ← I_i + max(0, −F_z,i)·dt            (accumulated downward impulse, N·s)
        K_i/K0   = 1 + β · (I_i / I_ref)              (linear K∼σ stiffening, capped at K_max/K0)

    where ``I_ref`` non-dimensionalises the impulse (a load scale) and ``β`` is the stiffening
    slope. This is a PHENOMENOLOGICAL hook (KB-1.V.2.5 marks the *emergent* fibre-network γc as
    not-a-knob — the genuine Mikado result lives in the ``ecm/`` module; this is the substrate-
    plane surrogate for the Warp DCM loop), so ``β``, ``I_ref`` and ``stiffen_max`` are flagged
    PROVISIONAL modelling parameters, default-OFF (``enabled=False``).

    Sanity Gate:
        * Dimensional: I [N·s]; I/I_ref dimensionless; K/K0 dimensionless ⇒ K [N/m] ✓.
        * enabled=False ⇒ k_well_local ≡ k_well_base (off-path invariance, byte-identical well).
        * Monotone: I_i is non-decreasing (only downward force accumulates) ⇒ K_i non-decreasing
          up to the cap ⇒ k rises with load, never falls.
        * Bounded: K_i/K0 ∈ [1, stiffen_max] (no runaway).
    """

    def __init__(
        self,
        *,
        n_nodes: int,
        k_well_base: float,
        dt: float,
        enabled: bool = False,
        beta: float = 1.0,
        I_ref: float = 1.0e-12,        # PROVISIONAL: N·s load scale (≈ per-clutch 5 pN × ~0.2 s)
        stiffen_max: float = 100.0,    # KB-1.4: K rises up to ~10–100× above γc → cap at 100×
    ):
        """Build the accumulator.

        Args:
            n_nodes: number of substrate (basal) nodes tracked.
            k_well_base: the scalar adhesive well stiffness (N/m) when unstiffened.
            dt: integrator timestep (s), for the impulse integral.
            enabled: master switch; False ⇒ no stiffening (default, byte-identical).
            beta: stiffening slope (dimensionless K∼σ gain). PROVISIONAL.
            I_ref: impulse non-dimensionalisation scale (N·s). PROVISIONAL — anchored to a
                per-clutch traction (KB-2.12 5–20 pN) over a bond lifetime (~0.1–1 s).
            stiffen_max: cap on K_i/K0 (KB-1.4 10–100× stiffening ceiling).
        """
        self.enabled = bool(enabled)
        self.k_well_base = float(k_well_base)
        self.dt = float(dt)
        self.beta = float(beta)
        self.I_ref = float(I_ref)
        self.stiffen_max = float(stiffen_max)
        self.impulse = np.zeros(int(n_nodes), dtype=np.float64)   # accumulated downward impulse
        self.k_well_local = np.full(int(n_nodes), self.k_well_base, dtype=np.float64)

    def accumulate(self, force_z: np.ndarray) -> None:
        """Add this step's downward (−z) traction to the per-node impulse and restiffen.

        Args:
            force_z: (N,) the z-component of the well force this node exerted (N). Only the
                downward part (the cell pushing INTO the dish, F_z < 0 ⇒ −F_z > 0) loads the
                substrate; upward/adhesive pull does not strain-stiffen it.
        """
        if not self.enabled:
            return
        fz = np.asarray(force_z, dtype=np.float64)
        self.impulse += np.maximum(0.0, -fz) * self.dt
        ratio = 1.0 + self.beta * (self.impulse / self.I_ref)
        np.clip(ratio, 1.0, self.stiffen_max, out=ratio)
        self.k_well_local = self.k_well_base * ratio

    def upload(self, device):
        """Push the per-node stiffness to a Warp array (for the strain-stiff kernel)."""
        return wp.array(np.ascontiguousarray(self.k_well_local), dtype=wp.float64, device=device)


def run_dcm_substrate_well_strainstiff_warp(
    *,
    pos: np.ndarray,
    z0: float,
    k_well_local: np.ndarray,   # (N,) per-node well stiffness (≥ base)
    adh_range: float,
    k_floor: float = 1.0,
    device: str = "cpu",
) -> dict:
    """Run the per-node strain-stiffened z-well kernel (helper for the loop / self-test).

    With ``k_well_local`` uniform (= base) this is byte-identical to the scalar well.
    """
    pos = np.ascontiguousarray(pos, dtype=np.float64)
    N = pos.shape[0]
    pos_d = wp.array(pos, dtype=wp.vec3d, device=device)
    kwl_d = wp.array(np.ascontiguousarray(k_well_local, dtype=np.float64),
                     dtype=wp.float64, device=device)
    force_d = wp.zeros(N, dtype=wp.vec3d, device=device)
    wp.launch(dcm_substrate_well_strainstiff_kernel, dim=N,
              inputs=[pos_d, wp.float64(z0), wp.float64(float(k_well_local.min())),
                      kwl_d, wp.float64(adh_range), wp.float64(k_floor), force_d],
              device=device)
    wp.synchronize_device(device)
    return {"force": force_d.numpy().astype(np.float64)}


# ──────────────────────────────────────────────────────────────────────────
# Self-test (C5 ECM mechano-feedback)
# ──────────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    import sys

    dev = "cpu"
    rng = np.random.default_rng(0)
    N = 64
    z0, W_cs, adh_range, k_floor = 0.0, 2.85e-3, 0.5e-6, 1.0
    # nodes scattered around the well (some above, some below z0)
    pos = rng.uniform(-1.0e-6, 1.0e-6, size=(N, 3))
    pos[:, 2] = rng.uniform(-0.4e-6, 0.4e-6, size=N)  # within adh_range band

    print("=" * 72)
    print("C5 ECM mechano-feedback self-test")
    print("=" * 72)

    # --- (A) OFF ⇒ unchanged forces (byte-identical to the legacy well) ----------
    base = run_dcm_substrate_well_warp(
        pos=pos, z0=z0, W_cs=W_cs, adh_range=adh_range, k_floor=k_floor, device=dev)
    mech_off = run_dcm_substrate_well_mechano_warp(
        pos=pos, z0=z0, W_cs=W_cs, adh_range=adh_range, k_floor=k_floor,
        E_sub=None, device=dev)
    off_match = np.array_equal(base["force"], mech_off["force"])
    print(f"\n[A] OFF byte-identical to legacy well : {off_match}  "
          f"(max|Δ|={np.abs(base['force']-mech_off['force']).max():.3e})")

    # --- (B) E_sub sweep ⇒ monotone traction toward a plateau --------------------
    Es = [1.0e2, 5.0e2, 1.0e3, 5.0e3, 1.0e4, 5.0e4, 1.0e5, 1.0e7, 1.0e9]
    print("\n[B] E_sub sweep (contact_radius=adh_range): traction factor + |F_z| sum")
    fz_sums, factors = [], []
    for E in Es:
        r = run_dcm_substrate_well_mechano_warp(
            pos=pos, z0=z0, W_cs=W_cs, adh_range=adh_range, k_floor=k_floor,
            E_sub=E, contact_radius=adh_range, device=dev)
        fzs = float(np.abs(r["force"][:, 2]).sum())
        fz_sums.append(fzs); factors.append(r["traction_factor"])
        print(f"    E={E:9.1e} Pa  factor={r['traction_factor']:.4f}  "
              f"k_sub={r['k_sub']:.3e}  Σ|F_z|={fzs:.4e}")
    monotone = all(b >= a - 1e-18 for a, b in zip(factors, factors[1:]))
    # rigid limit (E→∞): factor→1, force → legacy (since k_sub→∞ ⇒ series→k_well)
    rigid = run_dcm_substrate_well_mechano_warp(
        pos=pos, z0=z0, W_cs=W_cs, adh_range=adh_range, k_floor=k_floor,
        E_sub=float("inf"), contact_radius=adh_range, device=dev)
    rigid_match = np.allclose(rigid["force"], base["force"], rtol=0, atol=1e-18)
    print(f"    monotone traction factor          : {monotone}")
    print(f"    E_sub→∞ recovers rigid dish        : {rigid_match}  "
          f"(factor={rigid['traction_factor']:.6f}, max|Δ|="
          f"{np.abs(rigid['force']-base['force']).max():.3e})")

    # --- (C) strain-stiffening ⇒ k rises with load -------------------------------
    k_base = base["k_well"]
    ss = SubstrateStrainStiffening(n_nodes=N, k_well_base=k_base, dt=1.0e-4,
                                   enabled=True, beta=1.0)
    # OFF-invariance: a disabled accumulator never moves k_well_local
    ss_off = SubstrateStrainStiffening(n_nodes=N, k_well_base=k_base, dt=1.0e-4,
                                       enabled=False)
    ss_off.accumulate(np.full(N, -5.0e-9))      # downward load, but disabled
    off_inv = np.array_equal(ss_off.k_well_local, np.full(N, k_base))

    # apply a steady downward (−z) load for many steps and watch k climb
    k_hist = [ss.k_well_local.mean()]
    for _ in range(500):
        ss.accumulate(np.full(N, -5.0e-9))      # 5 nN downward per node
        k_hist.append(ss.k_well_local.mean())
    k_rises = k_hist[-1] > k_hist[0] and all(
        b >= a - 1e-30 for a, b in zip(k_hist, k_hist[1:]))
    capped = ss.k_well_local.max() <= k_base * ss.stiffen_max + 1e-30
    print(f"\n[C] strain-stiffening (enabled):")
    print(f"    OFF (disabled) leaves k unchanged  : {off_inv}")
    print(f"    k_well rises with load             : {k_rises}  "
          f"(k0={k_hist[0]:.3e} → k={k_hist[-1]:.3e} N/m, ×{k_hist[-1]/k_hist[0]:.2f})")
    print(f"    bounded by stiffen_max cap         : {capped}  "
          f"(max ×{ss.k_well_local.max()/k_base:.2f} ≤ {ss.stiffen_max})")
    # the stiffened kernel must produce LARGER |F_z| than the base for loaded nodes
    ss_force = run_dcm_substrate_well_strainstiff_warp(
        pos=pos, z0=z0, k_well_local=ss.k_well_local, adh_range=adh_range,
        k_floor=k_floor, device=dev)
    # compare on nodes inside the harmonic band with dz≠0 (well term active)
    dz = pos[:, 2] - z0
    band = np.abs(dz) <= adh_range
    stiff_more = np.all(
        np.abs(ss_force["force"][band, 2]) >= np.abs(base["force"][band, 2]) - 1e-30)
    print(f"    stiffened |F_z| ≥ base |F_z|        : {stiff_more}")

    ok = off_match and monotone and rigid_match and off_inv and k_rises and capped and stiff_more
    print("\n" + "=" * 72)
    print(f"SELF-TEST {'PASS' if ok else 'FAIL'}")
    print("=" * 72)
    sys.exit(0 if ok else 1)
