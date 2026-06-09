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
                 anchor_drag_factor=1.0e4):
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

    sim = hoomd.Simulation(device=device, seed=seed)
    sim.create_state_from_snapshot(snap)

    bond = md.bond.Harmonic()
    bond.params["sf-bond"] = dict(k=p.bond_k, r0=ell0)
    angle = md.angle.Harmonic()
    angle.params["sf-angle"] = dict(k=p.angle_k, t0=np.pi)  # straight aligned fiber
    nlist = md.nlist.Tree(buffer=0.5 * p.lj_sigma)
    lj = md.pair.LJ(nlist=nlist, default_r_cut=0.0)
    for ta in ("sf_actin", "fa_anchor"):
        for tb in ("sf_actin", "fa_anchor"):
            lj.params[(ta, tb)] = dict(epsilon=p.lj_epsilon, sigma=p.lj_sigma)
            lj.r_cut[(ta, tb)] = p.lj_r_cut if p.lj_enabled else 0.0
    lj.mode = "shift"

    ig = md.Integrator(dt=p.dt_cfl)
    ig.forces += [bond, angle, lj]
    sim.operations.integrator = ig
    # BAOAB integrates all beads. The FA anchors get a HIGH drag (overdamped rigid-substrate
    # limit: γ_anchor = anchor_drag_factor · γ_b) so they are ~fixed on the fiber timescale
    # WITHOUT a re-pin step (a hard re-pin after BAOAB is energetically inconsistent with the
    # L-M scheme — it pumps energy — and the integrator is PI-frozen, so no filter). The bond
    # tension delivered to these near-fixed anchors is the traction. integrator/ untouched.
    baoab_action, baoab_updater = make_baoab_updater(
        kT=p.kT,
        gamma={"sf_actin": p.gamma_b, "fa_anchor": anchor_drag_factor * p.gamma_b},
        dt=p.dt_cfl, seed=seed,
    )
    sim.operations.updaters.append(baoab_updater)
    anchor_tags = np.array(sorted(anchors), dtype=np.int64)
    return dict(sim=sim, bond=bond, layout=lay, anchors=anchor_tags,
                p=p, ell0=ell0, R=R, z_basal=z_basal)


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


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--n-filaments", type=int, default=12)
    ap.add_argument("--fiber-length-um", type=float, default=5.0)
    ap.add_argument("--bundle-radius-nm", type=float, default=200.0)
    ap.add_argument("--equilibrate", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--anchor-drag-factor", type=float, default=1.0e4,
                    help="FA-anchor drag multiple of gamma_b (overdamped rigid-substrate limit)")
    ap.add_argument("--device", choices=["cpu", "gpu"], default="cpu")
    ap.add_argument("--out-json", type=str,
                    default="ffn_sim/outputs/h7/production/h7_ventral_sf_traction.json")
    args = ap.parse_args()

    import hoomd
    dev = hoomd.device.GPU(notice_level=0) if args.device == "gpu" else hoomd.device.CPU(notice_level=0)
    h = build_sf_sim(
        n_fil=args.n_filaments, fiber_length=args.fiber_length_um * 1e-6,
        bundle_radius=args.bundle_radius_nm * 1e-9, device=dev, seed=args.seed,
        anchor_drag_factor=args.anchor_drag_factor,
    )
    n = h["layout"].positions.shape[0]
    print(f"  SF: {h['layout'].n_fil} fil × {h['layout'].n_beads_per_fil} beads = {n} "
          f"({h['anchors'].size} FA anchors), L={h['layout'].fiber_length*_UM:.2f}µm, "
          f"z_basal={h['z_basal']*_UM:.2f}µm", flush=True)
    h["sim"].run(1)
    t0, fmax0 = anchor_traction(h)
    print(f"  [t=1] passive anchor traction = {t0*_PN:.3f} pN  (max |F|={fmax0*_PN:.3f} pN)",
          flush=True)
    h["sim"].run(args.equilibrate)
    t1, fmax1 = anchor_traction(h)
    print(f"  [equilibrated {args.equilibrate}] passive anchor traction = {t1*_PN:.3f} pN  "
          f"(max |F|={fmax1*_PN:.3f} pN)", flush=True)
    # NOTE (2b-1 finding): a fiber placed at FULL CONTOUR extension is TAUT, so it carries an
    # intrinsic thermal/entropic tension at rest (a taut WLC has tension); the single-snapshot
    # |T| also fluctuates widely. So "passive ≈ 0" is NOT expected — the SCIENCE measurement
    # (Stage 2c) must be DIFFERENTIAL (myosin_ON − myosin_OFF) and TIME-AVERAGED, and/or place
    # the FA separation with slack. This run validates only the BUILD + RUN + traction READOUT.
    out = {
        "stage": "2b-1 SCAFFOLD: bundle build + FA anchor + traction readout (WIP)",
        "n_fil": h["layout"].n_fil, "n_beads": int(n), "n_anchors": int(h["anchors"].size),
        "fiber_length_um": h["layout"].fiber_length * _UM,
        "passive_traction_pN_t0": t0 * _PN,
        "passive_traction_pN_equilibrated": t1 * _PN,
        "known_issues_for_2b2": [
            "taut full-contour placement → intrinsic thermal tension (add slack or differential)",
            "single-snapshot |T| is noisy → time-average over samples",
            "overdamped high-drag anchor is a stand-in for the FA SubstrateLigandPin (consume it)",
            "myosin + alpha-actinin not yet added (Stage 2b-2)",
        ],
        "note": "BUILD+RUN+READOUT validated; quantitative traction = Stage 2c (differential, time-avg)",
    }
    Path(args.out_json).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out_json).write_text(json.dumps(out, indent=2))
    print(f"  json → {args.out_json}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
