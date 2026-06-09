"""H.7 ADHERENT pivot Stage 2b/2c — ventral stress-fiber TRACTION simulation.

Builds the fine-grained ventral stress fiber (cortex/ventral_stress_fiber layout: aligned
mixed-polarity bead-spring bundle, FA-anchored ends) into a runnable HOOMD system and measures
the TRACTION = the net axial reaction force the contracting fiber delivers to its FA anchors
(pinned = rigid-substrate limit). This is the direct test of the adherent-pivot claim
(H7_ADHERENT_VENTRAL_PIVOT §2): an aligned + END-ANCHORED fiber transmits myosin tension to the
substrate, where the suspended isotropic cortex did not.

Stage 2b-1 (here): assemble the bundle + backbone bonds/angles + WCA + BAOAB, pin the FA anchor
beads (separate ``fa_anchor`` type, excluded from BAOAB → fixed = rigid substrate), equilibrate,
and validate the traction-measurement plumbing (net force on the anchor set). A passive baseline
(no myosin) must give ~0 net anchor reaction at rest. Stage 2b-2 adds the Stam-Hocky myosin
(actin-aware placement) + α-actinin; Stage 2c reads the active traction and (separate track) the
stiffness-corrected traction.

NO lumped mechanism: explicit bead-spring actin, real harmonic backbone/bend, WCA EV. Rigid-
substrate pin = the sanctioned mcf7_baseline FA limit (consumed, not re-derived).
"""
from __future__ import annotations

import argparse
import json
from copy import deepcopy
from pathlib import Path

import numpy as np

from ffn_sim.cell.manifest import load_manifest, _deep_merge
from ffn_sim.cortex.cortex import resolve_h3_derived
from ffn_sim.cortex.ventral_stress_fiber import generate_ventral_sf_layout
from ffn_sim.integrator.baoab import make_baoab_updater

_PN = 1.0e12
_UM = 1.0e6


def _angle_triplets(layout):
    """Backbone bending triplets (i-1, i, i+1) per filament (straight rest angle)."""
    tri = []
    nb = layout.n_beads_per_fil
    for f in range(layout.n_fil):
        base = f * nb
        for j in range(1, nb - 1):
            tri.append((base + j - 1, base + j, base + j + 1))
    return np.array(tri, dtype=np.int64).reshape(-1, 3)


def build_sf_sim(*, n_fil, fiber_length, bundle_radius, device, seed,
                 anchor_drag_factor=1.0e4, with_myosin=False, n_motors=None):
    """Assemble the ventral SF HOOMD sim; optionally wire continuous_stroke myosin.

    Stage 2b-2: with_myosin=True places Stam-Hocky bipolar minifilaments (actin-aware,
    along the SF filament tangents ±x̂) on the bundle and drives them with the §9
    continuous per-head force (MyosinHeadForce; the cross-bridge stiffness fix). The
    SF actin beads are the binding substrate (n_cortex_actin = n_sf). α-actinin is the
    cortex dynamic crosslinker (separate, optional; not added here — the SF backbone
    bonds already provide the axial load-path to the FA anchors, which is the 2c lever)."""
    import hoomd
    import hoomd.md as md
    import gsd.hoomd

    manifest = deepcopy(load_manifest("mcf7_baseline.yaml"))
    base = load_manifest(manifest["base_cortex_config"])
    cfg = _deep_merge(base, manifest.get("cortex_overrides"))
    cfg.setdefault("cortex", {})["R_cell"] = float(manifest["R_cell"])
    p = resolve_h3_derived(cfg)
    ell0 = float(p.rest_length)
    R = float(p.R_cell)
    z_basal = -0.5 * R

    rng = np.random.default_rng(seed)
    lay = generate_ventral_sf_layout(
        n_fil=n_fil, fiber_length=fiber_length, ell0=ell0, z_basal=z_basal,
        bundle_radius=bundle_radius, rng=rng,
    )
    n = lay.positions.shape[0]
    anchors = set(int(b) for b in lay.anchor_beads)

    # particle types: actin bead vs pinned FA anchor (excluded from BAOAB → fixed)
    typeid = np.zeros(n, dtype=np.uint32)
    for b in anchors:
        typeid[b] = 1
    angles = _angle_triplets(lay)

    L_box = max(4.0 * fiber_length, 4.0 * R)
    snap = gsd.hoomd.Frame()
    snap.particles.N = n
    snap.particles.types = ["sf_actin", "fa_anchor"]
    snap.particles.typeid = typeid
    snap.particles.position = lay.positions
    snap.particles.mass = np.ones(n)
    snap.bonds.N = int(lay.backbone_bonds.shape[0])
    snap.bonds.types = ["sf-bond"]
    snap.bonds.typeid = np.zeros(lay.backbone_bonds.shape[0], dtype=np.uint32)
    snap.bonds.group = lay.backbone_bonds.astype(np.uint32)
    snap.angles.N = int(angles.shape[0])
    snap.angles.types = ["sf-angle"]
    snap.angles.typeid = np.zeros(angles.shape[0], dtype=np.uint32)
    snap.angles.group = angles.astype(np.uint32)
    snap.configuration.box = [L_box, L_box, L_box, 0.0, 0.0, 0.0]

    # ---- Stage 2b-2: continuous_stroke myosin on the SF bundle (actin-aware) ----
    p_myo = myo_layout = None
    n_sf = n
    if with_myosin:
        from ffn_sim.cortex.myosin import (
            resolve_cortex_myosin, generate_cortex_myosin_layout,
            extend_state_with_cortex_myosin,
        )
        myo_cfg = deepcopy(cfg)
        m = myo_cfg.setdefault("cortex", {}).setdefault("myosin", {})
        m["stepping_mode"] = "continuous_stroke"   # §9 capped per-head custom force
        m["mesoscale_force_scaling"] = False        # SF bundle: native motors, no sphere-density scaling
        if n_motors is not None:
            m["n_motors_per_cell"] = int(n_motors)
        p_myo = resolve_cortex_myosin(myo_cfg, dt=p.dt_cfl, R_cell=R)
        # Per-filament tangent = polarity·x̂ (the SF filaments run along ±x̂).
        fil_tangents = (lay.polarity[:, None].astype(np.float64)
                        * np.array([1.0, 0.0, 0.0])[None, :])
        myo_layout = generate_cortex_myosin_layout(
            p_myo, R, motor_tag_start=n_sf,
            rng=np.random.default_rng(p_myo.seed),
            cortex_positions=lay.positions,
            cortex_tangents=fil_tangents,
            beads_per_filament=lay.n_beads_per_fil,
            cortex_filament_idx=lay.filament_idx,
        )
        snap = extend_state_with_cortex_myosin(snap, myo_layout, p_myo)

    sim = hoomd.Simulation(device=device, seed=seed)
    sim.create_state_from_snapshot(snap)

    bond = md.bond.Harmonic()
    bond.params["sf-bond"] = dict(k=p.bond_k, r0=ell0)
    angle = md.angle.Harmonic()
    angle.params["sf-angle"] = dict(k=p.angle_k, t0=np.pi)  # straight aligned fiber
    nlist = md.nlist.Tree(buffer=0.5 * p.lj_sigma)
    lj = md.pair.LJ(nlist=nlist, default_r_cut=0.0)
    ptypes = list(sim.state.particle_types)   # includes myosin types when added
    for ta in ptypes:
        for tb in ptypes:
            lj.params[(ta, tb)] = dict(epsilon=p.lj_epsilon, sigma=p.lj_sigma)
            lj.r_cut[(ta, tb)] = p.lj_r_cut if p.lj_enabled else 0.0
    lj.mode = "shift"

    # CFL: the §9 stiff cross-bridge gives myosin k_backbone≈1e-2 → τ_backbone≈9 ns
    # < the bare-cortex dt_cfl 13 ns. Re-derive the SF integrator dt locally (same
    # cfl_safety_factor; NO edit to the frozen integrator/) so the stiff backbone is
    # CFL-safe. For the passive (no-myosin) build dt_used = p.dt_cfl unchanged.
    dt_used = float(p.dt_cfl)
    if with_myosin and p_myo is not None and getattr(p_myo, "k_backbone", 0.0) > 0.0:
        tau_backbone = p.gamma_b / p_myo.k_backbone
        dt_used = min(dt_used, p.cfl_safety_factor * tau_backbone)

    ig = md.Integrator(dt=dt_used)
    ig.forces += [bond, angle, lj]
    sim.operations.integrator = ig
    # BAOAB integrates all beads. The FA anchors get a HIGH drag (overdamped rigid-substrate
    # limit: γ_anchor = anchor_drag_factor · γ_b) so they are ~fixed on the fiber timescale
    # WITHOUT a re-pin step (a hard re-pin after BAOAB is energetically inconsistent with the
    # L-M scheme — it pumps energy — and the integrator is PI-frozen, so no filter). The bond
    # tension delivered to these near-fixed anchors is the traction. integrator/ untouched.
    gamma_map = {"sf_actin": p.gamma_b, "fa_anchor": anchor_drag_factor * p.gamma_b}
    if with_myosin:
        gamma_map["cortex_myosin_backbone"] = p.gamma_b
        gamma_map["cortex_myosin_head"] = p.gamma_b
    baoab_action, baoab_updater = make_baoab_updater(
        kT=p.kT, gamma=gamma_map, dt=dt_used, seed=seed,
    )
    sim.operations.updaters.append(baoab_updater)

    myosin_action = None
    if with_myosin:
        from ffn_sim.cortex.myosin import (
            register_cortex_myosin_bond_params, make_cortex_myosin_updater,
            MyosinHeadForce,
        )
        register_cortex_myosin_bond_params(bond, p_myo)
        myosin_action, myosin_updater = make_cortex_myosin_updater(
            p_myo=p_myo, layout=myo_layout, kT=p.kT, n_cortex_actin=n_sf,
            cortex_bond_groups=lay.backbone_bonds,
            ell0_cortex=ell0, cortex_beads_per_filament=lay.n_beads_per_fil,
        )
        sim.operations.updaters.append(myosin_updater)
        ig.forces.append(MyosinHeadForce(action=myosin_action, p_myo=p_myo))

    anchor_tags = np.array(sorted(anchors), dtype=np.int64)
    return dict(sim=sim, bond=bond, layout=lay, anchors=anchor_tags,
                p=p, ell0=ell0, R=R, z_basal=z_basal,
                p_myo=p_myo, myosin_action=myosin_action, dt_used=dt_used)


def anchor_traction(handles):
    """Traction = the net axial backbone-bond tension the fiber delivers to its FA anchors.

    Computed from bond GEOMETRY (T = k·(r − r0) along each anchor's backbone bond), robust and
    independent of HOOMD per-particle force-access timing. Returns (Σ|axial tension| [N],
    max |bond tension| [N])."""
    snap = handles["sim"].state.get_snapshot()
    if snap.communicator.rank != 0:
        return 0.0, 0.0
    pos = np.asarray(snap.particles.position, dtype=np.float64)
    axis = handles["layout"].axis
    k = float(handles["p"].bond_k)
    r0 = float(handles["ell0"])
    anchor_set = set(int(a) for a in handles["anchors"])
    bonds = handles["layout"].backbone_bonds
    axial_sum = 0.0
    max_T = 0.0
    for a, b in bonds:
        a, b = int(a), int(b)
        anc = a if a in anchor_set else (b if b in anchor_set else None)
        if anc is None:
            continue
        other = b if anc == a else a
        d = pos[other] - pos[anc]          # from anchor toward the fiber interior
        r = float(np.linalg.norm(d))
        if r <= 0:
            continue
        T = k * (r - r0)                   # backbone tension [N] (+ = stretched/pulling in)
        axial_sum += abs(T * float((d / r) @ axis))
        max_T = max(max_T, abs(T))
    return float(axial_sum), float(max_T)


def _time_avg_traction(h, *, equilibrate, contract, n_samples):
    """Equilibrate, then TIME-AVERAGE the axial anchor traction over n_samples chunks.

    Returns (mean_pN, std_pN, n_engaged_mean). Time-averaging is the 2b-1-mandated fix
    for the noisy single-snapshot |T|; the differential (active − passive) cancels the
    intrinsic taut-WLC thermal tension so what remains is the myosin-generated traction."""
    sim = h["sim"]
    sim.run(equilibrate)
    samples, engaged = [], []
    chunk = max(1, contract // n_samples)
    for _ in range(n_samples):
        sim.run(chunk)
        t, _ = anchor_traction(h)
        samples.append(t * _PN)
        act = h.get("myosin_action")
        engaged.append(int(act.n_engaged) if act is not None else 0)
    arr = np.array(samples)
    return float(arr.mean()), float(arr.std()), float(np.mean(engaged))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--n-filaments", type=int, default=12)
    ap.add_argument("--fiber-length-um", type=float, default=5.0)
    ap.add_argument("--bundle-radius-nm", type=float, default=200.0)
    ap.add_argument("--equilibrate", type=int, default=2000)
    ap.add_argument("--contract", type=int, default=20000,
                    help="post-equilibration steps over which traction is time-averaged (2c)")
    ap.add_argument("--n-samples", type=int, default=10)
    ap.add_argument("--n-motors", type=int, default=20)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--anchor-drag-factor", type=float, default=1.0e4,
                    help="FA-anchor drag multiple of gamma_b (overdamped rigid-substrate limit)")
    ap.add_argument("--device", choices=["cpu", "gpu"], default="cpu")
    ap.add_argument("--out-json", type=str,
                    default="ffn_sim/outputs/h7/production/h7_ventral_sf_traction.json")
    args = ap.parse_args()

    import hoomd

    def _dev():
        return (hoomd.device.GPU(notice_level=0) if args.device == "gpu"
                else hoomd.device.CPU(notice_level=0))

    common = dict(
        n_fil=args.n_filaments, fiber_length=args.fiber_length_um * 1e-6,
        bundle_radius=args.bundle_radius_nm * 1e-9, seed=args.seed,
        anchor_drag_factor=args.anchor_drag_factor,
    )
    # Stage 2c: DIFFERENTIAL time-averaged traction (active − passive). The passive
    # baseline carries the intrinsic taut-WLC thermal tension (2b-1 finding); the
    # differential isolates the myosin-generated traction with the §9 corrected motor.
    h_p = build_sf_sim(device=_dev(), with_myosin=False, **common)
    n = h_p["layout"].positions.shape[0]
    print(f"  SF: {h_p['layout'].n_fil} fil × {h_p['layout'].n_beads_per_fil} beads = {n} "
          f"({h_p['anchors'].size} FA anchors), L={h_p['layout'].fiber_length*_UM:.2f}µm",
          flush=True)
    t_pass, sd_pass, _ = _time_avg_traction(
        h_p, equilibrate=args.equilibrate, contract=args.contract, n_samples=args.n_samples)
    print(f"  PASSIVE  traction = {t_pass:.3f} ± {sd_pass:.3f} pN (time-avg)", flush=True)

    h_a = build_sf_sim(device=_dev(), with_myosin=True, n_motors=args.n_motors, **common)
    print(f"  myosin: {h_a['p_myo'].n_motors_per_cell} minifilaments, "
          f"continuous_stroke, dt={h_a['dt_used']:.3e}s (CFL-rederived)", flush=True)
    t_act, sd_act, eng = _time_avg_traction(
        h_a, equilibrate=args.equilibrate, contract=args.contract, n_samples=args.n_samples)
    print(f"  ACTIVE   traction = {t_act:.3f} ± {sd_act:.3f} pN (time-avg, "
          f"{eng:.0f} engaged heads)", flush=True)
    diff = t_act - t_pass
    print(f"  ⇒ DIFFERENTIAL (active − passive) = {diff:+.3f} pN  "
          f"[{'CONTRACTILE +traction' if diff > 0 else 'no net active traction'}]", flush=True)

    out = {
        "stage": "2b-2 + 2c: continuous_stroke myosin on ventral SF + differential time-avg traction",
        "n_fil": h_a["layout"].n_fil, "n_beads": int(n), "n_anchors": int(h_a["anchors"].size),
        "fiber_length_um": h_a["layout"].fiber_length * _UM,
        "n_motors": h_a["p_myo"].n_motors_per_cell,
        "stepping_mode": "continuous_stroke",
        "dt_used_s": h_a["dt_used"],
        "equilibrate": args.equilibrate, "contract": args.contract, "n_samples": args.n_samples,
        "passive_traction_pN": t_pass, "passive_traction_std_pN": sd_pass,
        "active_traction_pN": t_act, "active_traction_std_pN": sd_act,
        "engaged_heads_mean": eng,
        "differential_traction_pN": diff,
        "note": "differential cancels the intrinsic taut-WLC thermal tension → myosin-generated "
                "traction with the §9 corrected (capped, stiff cross-bridge) motor.",
    }
    Path(args.out_json).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out_json).write_text(json.dumps(out, indent=2))
    print(f"  json → {args.out_json}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
