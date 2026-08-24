"""D1 — nucleus + MT/IF/LINC internal force-path rig (Lane D, NON-AUTHORITATIVE DIAGNOSTIC).

Scene (built on top of the D0 membrane shell)::

    outer membrane  (own node block, tension + Helfrich)
        │  diagnostic capture patch  (anchor nodes on R_cortex, tied 1:1 to nearest membrane node)
        ▼
    explicit MT aster ── cortical/membrane-side capture (plus-ends → anchor)
        │  nesprin LINC (nonlinear tension-only, equal-and-opposite)   MT MTOC-side ↔ nucleus surface
        ▼
    deformable nucleus  (Helfrich + lamina areal tension + volume)
        ▲
    explicit IF cables  (radial nucleus↔anchor spokes, Hookean slice)  ── plectin/anchor

Read-only primitives reused
    * ``aleph.laws.membrane_surface``         — membrane mesh + tension/bending
    * ``aleph.components.nucleus.envelope``         — deformable nucleus (build_nucleus, lamina/volume/Helfrich, linc_tether)
    * ``aleph.components.nucleus.lamina_analytic``  — LaminaParams, linc_tether_force/energy host oracles
    * ``aleph.components.solid.microtubule``        — explicit MT aster (Cytosim bending, Euler-buckling oracle)
    * ``aleph.components.solid.intermediate_filament`` — explicit IF cage (Hookean spokes)

Lane-local glue (this file): a global-array assembler, a Hookean capture-spring kernel (diagnostic
membrane↔anchor / MT-tip↔anchor tether with a per-bond bound/detached flag), the two diagnostic
perturbations, and per-channel reaction readbacks for the force ledger.

⚠ DIAGNOSTIC PARAMETERS. Several nucleus/LINC constants are PI GAPs (lamin-A/C split, eps_rupture,
k_linc). The values below are **declared diagnostic placeholders** so the load-path can be exercised and
gated; they are NOT production values and are surfaced to PI. The capture-spring stiffness is a numerical
transmission coupling, not a biological constant. The 'capture patch' is a diagnostic patch, NOT a cortex
production proxy. Units µm-pN-s, float64, Warp-CUDA only.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import warp as wp

from aleph.laws.membrane_surface import (
    build_membrane_mesh, build_membrane_hinges, membrane_area_kernel,
    helfrich_bending_kernel, calibrate_kappa_tilde,
)
from aleph.components.nucleus.envelope import (
    build_nucleus, lamina_areal_tension_kernel, nucleoplasm_volume_reduce_kernel,
    nucleoplasm_volume_force_kernel, linc_tether_kernel,
)
from aleph.components.nucleus.lamina_analytic import LaminaParams
from aleph.components.solid.microtubule import build_microtubule_compartment
from aleph.components.solid.intermediate_filament import build_if_compartment

# --- diagnostic parameters (engine µm-pN-s) ---
R_MEM_UM = 7.5
R_CORTEX_UM = 7.0        # capture-patch (anchor) radius, just inside the membrane
R_NUC_UM = 5.0
GAMMA_MEM_EFF = 150.0    # membrane pre-tension (Young-Laplace, as D0) [pN/µm]
KAPPA_M_PN_UM = 0.0828
KAPPA_NE_PN_UM = 0.0828  # nuclear-envelope bending [pN·µm] (~20 kBT)
E_NUC_PA = 399.0         # MCF7 in-situ nuclear modulus [Pa]
H_LAMINA_UM = 0.015
K_VOL_NUC = 2.0e3        # numerical nucleoplasm incompressibility penalty [pN/µm²] (grid-invariant, diagnostic)
# LINC (nesprin) — k_linc is a PI GAP (HALT→PI). Diagnostic placeholder in the sourced 0.1–1 pN/nm band.
K_LINC_PN_UM = 200.0     # ⚠ DIAGNOSTIC placeholder (PI GAP) [pN/µm]
LINC_STIFFENING = 0.5    # nonlinear strain-stiffening [1/µm²]
# capture spring — numerical transmission coupling (NOT biological)
K_CAPTURE_PN_UM = 400.0


# ============================ lane-local glue kernels ============================
@wp.kernel
def _capture_spring_kernel(
    pos: wp.array(dtype=wp.vec3d),
    ia: wp.array(dtype=wp.int32),
    ib: wp.array(dtype=wp.int32),
    rest: wp.array(dtype=wp.float64),
    k: wp.float64,
    bound: wp.array(dtype=wp.int32),
    f: wp.array(dtype=wp.vec3d),
):
    """Hookean capture tether between global nodes ia[e]↔ib[e]. ``bound[e]==0`` (detached) removes only
    that edge's contribution (equal-and-opposite otherwise). Diagnostic transmission, not a lumped mechanism."""
    e = wp.tid()
    if bound[e] == 0:
        return
    a = ia[e]; b = ib[e]
    d = pos[b] - pos[a]
    L = wp.length(d)
    if L < wp.float64(1.0e-9):
        return
    u = d / L
    fv = (k * (L - rest[e])) * u
    wp.atomic_add(f, a, fv)
    wp.atomic_add(f, b, -fv)


@wp.kernel
def _overdamped_relax_kernel(
    pos: wp.array(dtype=wp.vec3d), force: wp.array(dtype=wp.vec3d),
    inv_gamma: wp.float64, dt: wp.float64, fixed: wp.array(dtype=wp.int32),
    pos_new: wp.array(dtype=wp.vec3d),
):
    """Overdamped relaxation; ``fixed[n]==1`` pins a node (prescribed-displacement boundary)."""
    n = wp.tid()
    if fixed[n] == 1:
        pos_new[n] = pos[n]
        return
    pos_new[n] = pos[n] + force[n] * (inv_gamma * dt)


def _fibonacci_sphere(n, r, centre=(0, 0, 0)):
    i = np.arange(n) + 0.5
    phi = np.arccos(1.0 - 2.0 * i / n)
    gold = np.pi * (1.0 + 5 ** 0.5)
    theta = gold * i
    x = np.cos(theta) * np.sin(phi); y = np.sin(theta) * np.sin(phi); z = np.cos(phi)
    return (np.stack([x, y, z], 1) * r + np.asarray(centre)).astype(np.float64)


@dataclass
class ForcePathRig:
    """Assembles membrane + nucleus + MT + IF + LINC into one global node array on CUDA."""

    device: str
    mem_subdiv: int = 2
    nuc_subdiv: int = 2
    n_mt: int = 12
    n_if: int = 40      # = n_anchor so every capture anchor has an IF spoke to the nucleus
    n_anchor: int = 40

    def __post_init__(self):
        dev = self.device
        # ---- membrane ----
        self.mem = build_membrane_mesh(R_MEM_UM, subdivisions=self.mem_subdiv)
        mem_v = np.ascontiguousarray(self.mem.verts, np.float64)
        self.mem_faces = np.ascontiguousarray(self.mem.faces, np.int32)
        self.mem_hinges = np.ascontiguousarray(build_membrane_hinges(self.mem.faces), np.int32)
        self.mem_kappa_tilde = float(calibrate_kappa_tilde(KAPPA_M_PN_UM, mem_v, self.mem_hinges))
        Nm = mem_v.shape[0]
        # ---- nucleus ----
        k_areal = E_NUC_PA * H_LAMINA_UM  # ≈5.985 pN/µm total areal modulus
        self.lamina = LaminaParams(k_chrom=0.3 * k_areal, k_lamin_b=0.7 * k_areal,
                                   k_lamin_ac=10.0 * k_areal, knee_strain=0.10,
                                   eps_rupture=0.5, kappa_ne=KAPPA_NE_PN_UM)  # ⚠ diagnostic split (PI GAP)
        self.nuc = build_nucleus(R_NUC_UM, self.lamina, aspect=1.0, subdivisions=self.nuc_subdiv)
        nuc_v = np.ascontiguousarray(self.nuc.verts, np.float64)
        Nn = nuc_v.shape[0]
        self.nuc_faces = np.ascontiguousarray(self.nuc.faces, np.int32)
        self.nuc_hinges = np.ascontiguousarray(self.nuc.hinges, np.int32)
        self.nuc_a0 = np.ascontiguousarray(self.nuc.A0_face, np.float64)
        self.nuc_V0 = float(self.nuc.V0)
        self.nuc_kappa_tilde = float(self.nuc.kappa_tilde)
        # ---- anchor (capture patch) ----
        anc = _fibonacci_sphere(self.n_anchor, R_CORTEX_UM)
        Na = anc.shape[0]
        # ---- block offsets ----
        self.off_mem, self.off_nuc, self.off_anc = 0, Nm, Nm + Nn
        off_mt = Nm + Nn + Na
        # ---- MT aster (node_off into global array) ----
        self.mt = build_microtubule_compartment(centre=(0, 0, 0), n_mt=self.n_mt,
                                                 reach_R_um=R_CORTEX_UM, node_off=off_mt, device=dev)
        mt_p = np.ascontiguousarray(self.mt.pos0, np.float64)
        Nmt = mt_p.shape[0]
        self.off_mt = off_mt
        off_if = off_mt + Nmt
        # ---- IF cage (nucleus↔anchor radial spokes) ----
        self.iff = build_if_compartment(centre=(0, 0, 0), R_nuc_um=R_NUC_UM, R_cortex_um=R_CORTEX_UM,
                                        nuc_pos=nuc_v, cortex_pos=anc, n_fil=self.n_if, cell_type="keratin",
                                        nuc_offset=self.off_nuc, cortex_offset=self.off_anc,
                                        if_offset=off_if, crosslink=True, device=dev)
        if_p = np.ascontiguousarray(self.iff.pos0, np.float64)
        Nif = if_p.shape[0]
        self.off_if = off_if
        # ---- global node array ----
        pos_all = np.concatenate([mem_v, nuc_v, anc, mt_p, if_p], axis=0)
        self.pos0 = np.ascontiguousarray(pos_all, np.float64)
        self.N = self.pos0.shape[0]
        self.Nm, self.Nn, self.Na, self.Nmt, self.Nif = Nm, Nn, Na, Nmt, Nif
        self.mem_slice = slice(0, Nm)
        self.nuc_slice = slice(self.off_nuc, self.off_nuc + Nn)
        self.anc_slice = slice(self.off_anc, self.off_anc + Na)

        # ---- capture edges: anchor↔nearest membrane node, and MT tip↔nearest anchor ----
        cap_a, cap_b = [], []
        for j in range(Na):
            m = int(np.argmin(np.linalg.norm(mem_v - anc[j], axis=1)))
            cap_a.append(self.off_anc + j); cap_b.append(self.off_mem + m)
        # MT plus-end (tip) node per fiber ↔ nearest anchor
        foff = np.asarray(self.mt.fiber_off, np.int64)
        self._mt_tip_global, self._mt_mtoc_global = [], []
        for i in range(self.n_mt):
            tip_local = int(foff[i + 1] - 1); mtoc_local = int(foff[i])
            tip_g = off_mt + tip_local; mtoc_g = off_mt + mtoc_local
            self._mt_tip_global.append(tip_g); self._mt_mtoc_global.append(mtoc_g)
            j = int(np.argmin(np.linalg.norm(anc - mt_p[tip_local], axis=1)))
            cap_a.append(tip_g); cap_b.append(self.off_anc + j)
        cap_a = np.asarray(cap_a, np.int32); cap_b = np.asarray(cap_b, np.int32)
        cap_rest = np.linalg.norm(self.pos0[cap_b] - self.pos0[cap_a], axis=1).astype(np.float64)
        self.n_capture = cap_a.shape[0]

        # ---- LINC (nesprin) edges: each anchor ↔ nearest nucleus surface node (IF/actin-cap LINC).
        # The MT primitive here is bending-only (no axial stretch), so axial load travels the IF+LINC
        # cabling; MT contributes lateral bending/buckling resistance and is gated on its own strut.
        linc_a, linc_n = [], []
        for j in range(Na):
            k = int(np.argmin(np.linalg.norm(nuc_v - anc[j], axis=1)))
            linc_a.append(self.off_anc + j)             # anchor side
            linc_n.append(self.off_nuc + k)             # nucleus side
        self._linc_n = np.asarray(linc_n, np.int32)
        self._linc_a = np.asarray(linc_a, np.int32)
        linc_rest = np.linalg.norm(self.pos0[self._linc_a] - self.pos0[self._linc_n], axis=1).astype(np.float64)
        self.n_linc = self._linc_n.shape[0]

        # ---- device arrays ----
        self.pos_d = wp.array(self.pos0, dtype=wp.vec3d, device=dev)
        self.pos_new_d = wp.zeros(self.N, dtype=wp.vec3d, device=dev)
        self.f_d = wp.zeros(self.N, dtype=wp.vec3d, device=dev)
        self.fixed_d = wp.zeros(self.N, dtype=wp.int32, device=dev)
        self.mem_faces_d = wp.array(self.mem_faces, dtype=wp.int32, device=dev)
        self.mem_hinges_d = wp.array(self.mem_hinges, dtype=wp.int32, device=dev)
        self.nuc_faces_d = wp.array(self.nuc_faces, dtype=wp.int32, device=dev)
        self.nuc_hinges_d = wp.array(self.nuc_hinges, dtype=wp.int32, device=dev)
        self.nuc_a0_d = wp.array(self.nuc_a0, dtype=wp.float64, device=dev)
        self.nuc_ruptured_d = wp.zeros(self.nuc_faces.shape[0], dtype=wp.int32, device=dev)
        self._nuc_vol_d = wp.zeros(1, dtype=wp.float64, device=dev)
        self.cap_ia_d = wp.array(cap_a, dtype=wp.int32, device=dev)
        self.cap_ib_d = wp.array(cap_b, dtype=wp.int32, device=dev)
        self.cap_rest_d = wp.array(cap_rest, dtype=wp.float64, device=dev)
        self.cap_bound_d = wp.ones(self.n_capture, dtype=wp.int32, device=dev)
        self.linc_n_d = wp.array(self._linc_n, dtype=wp.int32, device=dev)
        self.linc_a_d = wp.array(self._linc_a, dtype=wp.int32, device=dev)
        self.linc_rest_d = wp.array(linc_rest, dtype=wp.float64, device=dev)
        self.linc_bound = np.ones(self.n_linc, dtype=bool)   # host-side detach mask (per-joint)
        self.topology_epoch = 0                              # generation counter (MT remap)
        self.f_rest_d = wp.zeros(self.N, dtype=wp.vec3d, device=dev)
        # The membrane is a PINNED diagnostic frame — with no cytosol pressure (D1 isolates the internal
        # load path) its tension would collapse it. Pinning makes it a boundary the capture patch loads
        # against; the actuated patch nodes stay pinned and are relocated by the perturbations.
        fx = self.fixed_d.numpy(); fx[self.mem_slice] = 1
        self.fixed_d.assign(np.ascontiguousarray(fx))
        # Reference V0 to the DISCRETE mesh volume (build_nucleus uses the analytic oblate volume, which
        # differs from the subdivided-icosphere volume → a phantom volume-penalty force at build). Measuring
        # it here makes the nucleus genuinely force-free at its built shape.
        self._nuc_vol_d.zero_()
        wp.launch(nucleoplasm_volume_reduce_kernel, dim=self.nuc_faces.shape[0],
                  inputs=[self.pos_d, self.nuc_faces_d, self._nuc_vol_d], device=dev)
        self.nuc_V0 = float(self._nuc_vol_d.numpy()[0])

    # ---------------- force channels ----------------
    def _acc_nucleus(self, pos, f):
        wp.launch(helfrich_bending_kernel, dim=self.nuc_hinges.shape[0],
                  inputs=[pos, self.nuc_hinges_d, wp.float64(self.nuc_kappa_tilde), f], device=self.device)
        wp.launch(lamina_areal_tension_kernel, dim=self.nuc_faces.shape[0],
                  inputs=[pos, self.nuc_faces_d, self.nuc_a0_d, self.nuc_ruptured_d,
                          wp.float64(self.lamina.k_soft), wp.float64(self.lamina.k_lamin_ac),
                          wp.float64(self.lamina.knee_strain), wp.float64(self.lamina.eps_rupture), f],
                  device=self.device)
        # volume penalty
        self._nuc_vol_d.zero_()
        wp.launch(nucleoplasm_volume_reduce_kernel, dim=self.nuc_faces.shape[0],
                  inputs=[pos, self.nuc_faces_d, self._nuc_vol_d], device=self.device)
        # NOTE: p_vol needs a host read of V; acceptable post-loop / per-outer-step (diagnostic)
        V = float(self._nuc_vol_d.numpy()[0])
        p_vol = -K_VOL_NUC * (V - self.nuc_V0) / self.nuc_V0
        wp.launch(nucleoplasm_volume_force_kernel, dim=self.nuc_faces.shape[0],
                  inputs=[pos, self.nuc_faces_d, wp.float64(p_vol), f], device=self.device)

    def _acc_membrane(self, pos, f):
        wp.launch(membrane_area_kernel, dim=self.mem_faces.shape[0],
                  inputs=[pos, self.mem_faces_d, wp.float64(GAMMA_MEM_EFF), f], device=self.device)
        wp.launch(helfrich_bending_kernel, dim=self.mem_hinges.shape[0],
                  inputs=[pos, self.mem_hinges_d, wp.float64(self.mem_kappa_tilde), f], device=self.device)

    def _acc_capture(self, pos, f):
        wp.launch(_capture_spring_kernel, dim=self.n_capture,
                  inputs=[pos, self.cap_ia_d, self.cap_ib_d, self.cap_rest_d, wp.float64(K_CAPTURE_PN_UM),
                          self.cap_bound_d, f], device=self.device)

    def _acc_linc(self, pos, f):
        wp.launch(linc_tether_kernel, dim=self.n_linc,
                  inputs=[pos, pos, self.linc_n_d, self.linc_a_d, self.linc_rest_d,
                          wp.float64(K_LINC_PN_UM), wp.float64(LINC_STIFFENING), f, f], device=self.device)

    def accumulate_all(self, pos=None, f=None, subtract_rest=True):
        """Accumulate every channel into ``f`` (defaults to self.f_d, zeroed first)."""
        pos = self.pos_d if pos is None else pos
        f = self.f_d if f is None else f
        f.zero_()
        # membrane is a pinned frame (see __post_init__) — its tension is not applied in D1 (no cytosol
        # to balance it); the capture springs transmit membrane-node motion into the internal network.
        self._acc_nucleus(pos, f)
        self.mt.accumulate(pos, f)
        self.iff.accumulate(pos, f)
        self._acc_capture(pos, f)
        self._acc_linc(pos, f)
        if subtract_rest:
            wp.launch(_sub_ref_kernel, dim=self.N, inputs=[f, self.f_rest_d], device=self.device)
        return f

    # ---------------- rest reference ----------------
    def relax_to_rest(self, n_iter: int = 400, step: float = 2.0e-4) -> float:
        self.f_rest_d.zero_()
        for _ in range(n_iter):
            self.accumulate_all(subtract_rest=False)
            wp.launch(_overdamped_relax_kernel, dim=self.N,
                      inputs=[self.pos_d, self.f_d, wp.float64(step), wp.float64(1.0),
                              self.fixed_d, self.pos_new_d], device=self.device)
            self.pos_d.assign(self.pos_new_d)
        self.pos0 = self.pos_d.numpy().copy()
        self.accumulate_all(subtract_rest=False)
        self.f_rest_d.assign(self.f_d)
        self.accumulate_all(subtract_rest=True)
        return float(np.max(np.linalg.norm(self.f_d.numpy(), axis=1)))

    def relax_and_return(self, n_iter: int = 300):
        """relax_to_rest then return self (convenience for one-line gate construction)."""
        self.relax_to_rest(n_iter=n_iter)
        return self

    def settle(self, n_iter: int = 600, step: float = 2.0e-4) -> float:
        """Overdamped quasi-static relaxation to equilibrium (respecting pinned nodes). Returns max residual."""
        for _ in range(n_iter):
            self.accumulate_all(subtract_rest=True)
            wp.launch(_overdamped_relax_kernel, dim=self.N,
                      inputs=[self.pos_d, self.f_d, wp.float64(step), wp.float64(1.0),
                              self.fixed_d, self.pos_new_d], device=self.device)
            self.pos_d.assign(self.pos_new_d)
        self.accumulate_all(subtract_rest=True)
        free = self.fixed_d.numpy() == 0
        return float(np.max(np.linalg.norm(self.f_d.numpy()[free], axis=1))) if free.any() else 0.0

    # ---------------- perturbations ----------------
    def membrane_patch_nodes(self, axis_frac=0.6):
        c = self.pos0[self.mem_slice].mean(axis=0)
        v = self.pos0[self.mem_slice]
        return np.where((v[:, 0] - c[0]) > axis_frac * R_MEM_UM)[0]  # local membrane indices

    def nucleus_nodes(self):
        return np.arange(self.Nn) + self.off_nuc

    def set_fixed(self, global_idx, fixed=True):
        f = self.fixed_d.numpy(); f[global_idx] = 1 if fixed else 0
        self.fixed_d.assign(np.ascontiguousarray(f))

    def prescribe_displacement(self, global_idx, delta):
        """Move nodes by ``delta`` (µm) and pin them — an EXTERNAL input logged in the ledger (pulse B)."""
        p = self.pos_d.numpy(); p[global_idx] = p[global_idx] + np.asarray(delta)[None, :]
        self.pos_d.assign(np.ascontiguousarray(p))
        self.set_fixed(global_idx, True)

    # ---------------- reaction readbacks (ledger) ----------------
    def channel_force(self, channel: str):
        """Return the per-node force from a single channel [pN] (for the force ledger / sign gates)."""
        f = wp.zeros(self.N, dtype=wp.vec3d, device=self.device)
        if channel == "linc":
            self._acc_linc(self.pos_d, f)
        elif channel == "capture":
            self._acc_capture(self.pos_d, f)
        elif channel == "mt":
            self.mt.accumulate(self.pos_d, f)
        elif channel == "if":
            self.iff.accumulate(self.pos_d, f)
        elif channel == "nucleus":
            self._acc_nucleus(self.pos_d, f)
        else:
            raise ValueError(channel)
        return f.numpy()

    def linc_detach(self, joint: int):
        """Detach one LINC joint — remove ONLY that graph path (its rest set huge ⇒ tension-only ⇒ zero)."""
        self.linc_bound[joint] = False
        n = self.linc_n_d.numpy(); a = self.linc_a_d.numpy()
        # route the detached joint to a self-pair (a==n) so its contribution is exactly zero, others intact
        n[joint] = a[joint]
        self.linc_n_d.assign(np.ascontiguousarray(n))

    def bump_topology_epoch(self):
        self.topology_epoch += 1
        return self.topology_epoch

    # ---------------- snapshot / restore (transactional reject) ----------------
    def snapshot(self):
        return {
            "pos": self.pos_d.numpy().copy(),
            "linc_n": self.linc_n_d.numpy().copy(),
            "cap_bound": self.cap_bound_d.numpy().copy(),
            "linc_bound": self.linc_bound.copy(),
            "epoch": self.topology_epoch,
            "fixed": self.fixed_d.numpy().copy(),
        }

    def restore(self, snap):
        self.pos_d.assign(np.ascontiguousarray(snap["pos"]))
        self.linc_n_d.assign(np.ascontiguousarray(snap["linc_n"]))
        self.cap_bound_d.assign(np.ascontiguousarray(snap["cap_bound"]))
        self.linc_bound = snap["linc_bound"].copy()
        self.topology_epoch = snap["epoch"]
        self.fixed_d.assign(np.ascontiguousarray(snap["fixed"]))

    def capture_detach(self, edge: int):
        b = self.cap_bound_d.numpy(); b[edge] = 0
        self.cap_bound_d.assign(np.ascontiguousarray(b))


@wp.kernel
def _sub_ref_kernel(force: wp.array(dtype=wp.vec3d), f_ref: wp.array(dtype=wp.vec3d)):
    n = wp.tid()
    force[n] = force[n] - f_ref[n]
