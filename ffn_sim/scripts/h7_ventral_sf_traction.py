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
                 anchor_drag_factor=1.0e4, with_myosin=False, n_motors=None,
                 myosin_force_scale=1.0, myosin_backbone_bending=True,
                 sarcomeric=False, n_sarcomeres=3, overlap_frac=0.3):
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
    if sarcomeric:
        # n_fil = filaments per cross-section (bundle thickness); total =
        # n_sarcomeres·2·n_fil half-filaments.
        from ffn_sim.cortex.ventral_stress_fiber import generate_sarcomeric_sf_layout
        lay = generate_sarcomeric_sf_layout(
            n_cross=n_fil, fiber_length=fiber_length, ell0=ell0, z_basal=z_basal,
            bundle_radius=bundle_radius, rng=rng,
            n_sarcomeres=n_sarcomeres, overlap_frac=overlap_frac,
        )
        # α-actinin Z-disc stiffness = the re-anchored crosslink k (KB-1.28).
        k_zdisc = float(cfg["cortex"]["dynamic_crosslinkers"]["k_attach"])
    else:
        lay = generate_ventral_sf_layout(
            n_fil=n_fil, fiber_length=fiber_length, ell0=ell0, z_basal=z_basal,
            bundle_radius=bundle_radius, rng=rng,
        )
        k_zdisc = 0.0
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
    n_bb = int(lay.backbone_bonds.shape[0])
    n_xl = int(lay.crosslink_bonds.shape[0])   # α-actinin Z-disc (sarcomeric only)
    if n_xl > 0:
        snap.bonds.N = n_bb + n_xl
        snap.bonds.types = ["sf-bond", "sf_zdisc"]
        snap.bonds.typeid = np.concatenate([
            np.zeros(n_bb, dtype=np.uint32), np.ones(n_xl, dtype=np.uint32)])
        snap.bonds.group = np.concatenate([
            lay.backbone_bonds.astype(np.uint32),
            lay.crosslink_bonds.astype(np.uint32)], axis=0)
    else:
        snap.bonds.N = n_bb
        snap.bonds.types = ["sf-bond"]
        snap.bonds.typeid = np.zeros(n_bb, dtype=np.uint32)
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
        m["backbone_bending"] = bool(myosin_backbone_bending)  # rigid-rod fidelity: keep the dipole arm straight
        if n_motors is not None:
            m["n_motors_per_cell"] = int(n_motors)
        p_myo = resolve_cortex_myosin(myo_cfg, dt=p.dt_cfl, R_cell=R, kT=p.kT)
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
    if n_xl > 0:
        # α-actinin Z-disc crosslink: barbed ends meeting at a Z-band are held
        # together (r0=0 idealised Z-disc; the bond's nlist exclusion stops the WCA
        # between the co-located pair). k = KB-1.28 re-anchored crosslink stiffness.
        bond.params["sf_zdisc"] = dict(k=k_zdisc, r0=0.0)
    angle = md.angle.Harmonic()
    angle.params["sf-angle"] = dict(k=p.angle_k, t0=np.pi)  # straight aligned fiber
    if with_myosin and p_myo is not None and p_myo.backbone_bending:
        from ffn_sim.cortex.myosin import register_cortex_myosin_angle_params
        register_cortex_myosin_angle_params(angle, p_myo)   # minifilament rigid-rod
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
        if getattr(p_myo, "backbone_bending", False) and p_myo.k_backbone_angle > 0.0:
            # τ_bend = γ_b·ℓ0³/κ_B (κ_B = k_angle·ℓ0); usually slacker than stretch.
            kappa_B = p_myo.k_backbone_angle * p_myo.backbone_segment_length
            tau_bend = p.gamma_b * p_myo.backbone_segment_length ** 3 / kappa_B
            dt_used = min(dt_used, p.cfl_safety_factor * tau_bend)

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
        ig.forces.append(MyosinHeadForce(
            action=myosin_action, p_myo=p_myo, force_scale=myosin_force_scale))

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
    # COHERENT contractile traction (the physical observable) vs Σ|T| (fluctuation
    # magnitude). For each anchor bond, the SIGNED axial force the bond exerts ON the
    # anchor: f_in = k(r−r0)·(d̂·axis) with d toward the fiber interior, so f_in>0 =
    # pulled INWARD (contractile). Split anchors into the two ends about the fiber
    # centre; the contractile tension pulls the two ends TOWARD each other. Summing the
    # signed inward force per end makes thermal fluctuations cancel (random sign), so
    # the differential isolates the myosin's COHERENT contraction — unlike Σ|axial T|,
    # which sums absolute values and is dominated by the thermal fluctuation magnitude.
    xc = float(np.mean(pos[np.asarray(sorted(anchor_set), dtype=np.int64), 0]))
    axial_abs = 0.0
    max_T = 0.0
    f_left = 0.0   # net inward (+axis) force on the LEFT-end anchors
    f_right = 0.0  # net inward (−axis) force on the RIGHT-end anchors
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
        f_axial = T * float((d / r) @ axis)   # signed force on the anchor along +axis
        axial_abs += abs(f_axial)
        max_T = max(max_T, abs(T))
        if pos[anc, 0] < xc:               # left end: inward = +axis
            f_left += f_axial
        else:                              # right end: inward = −axis
            f_right += -f_axial
    # Coherent contractile traction = mean inward pull of the two ends (>0 = contractile).
    coherent = 0.5 * (f_left + f_right)
    return float(axial_abs), float(max_T), float(coherent)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--n-filaments", type=int, default=24)
    ap.add_argument("--fiber-length-um", type=float, default=6.0)
    ap.add_argument("--bundle-radius-nm", type=float, default=400.0)
    ap.add_argument("--equilibrate", type=int, default=2000)
    ap.add_argument("--contract", type=int, default=20000,
                    help="post-equilibration steps over which traction is time-averaged (2c)")
    ap.add_argument("--n-samples", type=int, default=10)
    ap.add_argument("--n-motors", type=int, default=20)
    ap.add_argument("--sarcomeric", action="store_true",
                    help="graded-polarity SARCOMERIC SF (Z/M bands + α-actinin) instead of "
                         "the random mixed-polarity bundle — rectifies myosin sliding into traction")
    ap.add_argument("--n-sarcomeres", type=int, default=3)
    ap.add_argument("--overlap-frac", type=float, default=0.3)
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
        with_myosin=True, n_motors=args.n_motors,
        sarcomeric=args.sarcomeric, n_sarcomeres=args.n_sarcomeres,
        overlap_frac=args.overlap_frac,
    )
    # Stage 2c REDESIGN — SAME-SEED PAIRED differential. The first design differenced
    # two DIFFERENT random realizations (myosin-present vs myosin-absent), so the huge
    # taut-WLC thermal tension (~10 nN) did NOT cancel and buried the ~tens-of-pN myosin
    # signal. Here BOTH builds are bit-identical (same seed → same bundle + same myosin
    # placement + same binding + same thermostat counter), differing ONLY in the myosin
    # FORCE: force_scale=0 (myosin present, force OFF) vs 1 (force ON). HOOMD's
    # counter-based RNG makes the thermal kicks identical at every (timestep, tag)
    # regardless of position, so the differential cancels the thermal noise and the taut
    # baseline EXACTLY — what remains is the pure myosin-generated traction. Sampled at
    # MATCHED timesteps and differenced per-sample.
    h_off = build_sf_sim(device=_dev(), myosin_force_scale=0.0, **common)
    h_on = build_sf_sim(device=_dev(), myosin_force_scale=1.0, **common)
    n = h_on["layout"].positions.shape[0]
    print(f"  SF: {h_on['layout'].n_fil} fil × {h_on['layout'].n_beads_per_fil} beads = {n} "
          f"({h_on['anchors'].size} FA anchors), L={h_on['layout'].fiber_length*_UM:.2f}µm; "
          f"{h_on['p_myo'].n_motors_per_cell} minifilaments, dt={h_on['dt_used']:.3e}s",
          flush=True)
    # Equilibrate both identically (same seed → identical trajectories up to here).
    h_off["sim"].run(args.equilibrate)
    h_on["sim"].run(args.equilibrate)
    chunk = max(1, args.contract // args.n_samples)
    # COHERENT contractile traction is the physical observable (signed; thermal cancels).
    cdiffs, c_offs, c_ons, engs = [], [], [], []
    for _ in range(args.n_samples):
        h_off["sim"].run(chunk)
        h_on["sim"].run(chunk)
        _, _, c_off = anchor_traction(h_off)
        _, _, c_on = anchor_traction(h_on)
        cdiffs.append((c_on - c_off) * _PN)
        c_offs.append(c_off * _PN); c_ons.append(c_on * _PN)
        engs.append(int(h_on["myosin_action"].n_engaged))
    cdiffs = np.array(cdiffs)
    diff_mean = float(cdiffs.mean())
    diff_sem = float(cdiffs.std() / max(1, len(cdiffs) ** 0.5))
    eng = float(np.mean(engs))
    sig = abs(diff_mean) > 2 * diff_sem
    verdict = ("CONTRACTILE (+traction)" if (sig and diff_mean > 0)
               else "EXPANSILE/slackening (−)" if (sig and diff_mean < 0)
               else "within noise")
    print(f"  COHERENT force-OFF = {np.mean(c_offs):+.2f} ± {np.std(c_offs):.2f} pN", flush=True)
    print(f"  COHERENT force-ON  = {np.mean(c_ons):+.2f} ± {np.std(c_ons):.2f} pN "
          f"({eng:.0f} engaged heads)", flush=True)
    print(f"  ⇒ SAME-SEED COHERENT DIFFERENTIAL (ON − OFF) = {diff_mean:+.3f} ± {diff_sem:.3f} pN "
          f"[{verdict}]", flush=True)

    out = {
        "stage": "2c REDESIGN: same-seed paired COHERENT-traction differential (force ON vs OFF)",
        "n_fil": h_on["layout"].n_fil, "n_beads": int(n), "n_anchors": int(h_on["anchors"].size),
        "fiber_length_um": h_on["layout"].fiber_length * _UM,
        "n_motors": h_on["p_myo"].n_motors_per_cell,
        "stepping_mode": "continuous_stroke", "dt_used_s": h_on["dt_used"],
        "equilibrate": args.equilibrate, "contract": args.contract, "n_samples": args.n_samples,
        "coherent_off_traction_pN": float(np.mean(c_offs)),
        "coherent_on_traction_pN": float(np.mean(c_ons)),
        "engaged_heads_mean": eng,
        "coherent_differential_pN": diff_mean,
        "coherent_differential_sem_pN": diff_sem,
        "per_sample_coherent_differential_pN": cdiffs.tolist(),
        "significant": bool(sig), "verdict": verdict,
        "note": "same-seed paired (force_scale 1 vs 0): identical bundle/binding/thermostat → "
                "thermal noise cancels. COHERENT signed traction (two ends pulled together) is "
                "the contractile observable (vs Σ|T| fluctuation magnitude). +ve = the §9 "
                "corrected myosin generates net contractile traction at the FA anchors.",
    }
    Path(args.out_json).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out_json).write_text(json.dumps(out, indent=2))
    print(f"  json → {args.out_json}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
