"""S6 remodel field cross-check vs the newly-arrived KB-1 ECM claims (KB-1.10 1/r, KB-1.9 nematic S)."""
import sys
import numpy as np

npz = sys.argv[1] if len(sys.argv) > 1 else "ffn_sim/outputs/ff/figs/ff_ecm_s6_long_on.npz"
d = np.load(npz)
p0 = np.asarray(d["ecm_pos0"], float)      # (Ne,3) frame-0 collagen
pf = np.asarray(d["ecm_posf"], float)      # (Ne,3) remodeled collagen
en = np.asarray(d["ecm_node"])             # (M,) clutch->collagen node
foff = np.asarray(d["ecm_foff"], int)      # fiber offsets
bnd = en >= 0
gr = en[bnd]
fc = p0[gr, :2].mean(0)                     # footprint centroid (xy)
disp = np.linalg.norm(pf - p0, axis=1)      # per-node displacement magnitude [µm]

# ---- KB-1.10: long-range stress/displacement propagation σ(r) ~ r^-n; n~1 fibrous, n=3 continuum ----
r = np.linalg.norm(p0[:, :2] - fc[None, :], axis=1)                 # in-plane distance from footprint [µm]
m = (disp > 1e-4) & (r > 0.5) & (r < 12.0)                          # moved nodes, in the near-to-mid field
if m.sum() > 30:
    lr, ld = np.log(r[m]), np.log(disp[m])
    n_fit = -np.polyfit(lr, ld, 1)[0]                               # slope of log-disp vs log-r → exponent n
    # binned profile for a robust read
    bins = np.linspace(0.5, 12, 13); idx = np.digitize(r[m], bins)
    prof = [(0.5 * (bins[i] + bins[i + 1]), disp[m][idx == i + 1].mean() * 1e3)
            for i in range(len(bins) - 1) if (idx == i + 1).sum() > 3]
    print(f"[KB-1.10 σ(r)~r^-n]  fitted n = {n_fit:.2f}  (KB: n~1 fibrous force-chains, n=3 continuum)")
    print("   radial displacement profile (r µm → disp nm):")
    for rr, dd in prof:
        print(f"     r={rr:4.1f}  {dd:6.1f} nm  " + "#" * int(min(dd / 5, 50)))
else:
    print("[KB-1.10] too few moved nodes for a decay fit")

# ---- KB-1.9: nematic order S = <cos 2(θ-θ0)>; radial (θ0 = direction to footprint) alignment near the cell ----
def fiber_nematic_radial(pos, rmax):
    """S for near-footprint fibers, director = radial-to-centroid. S=1 fully radial, 0 isotropic."""
    cos2 = []
    for f in range(len(foff) - 1):
        a, b = int(foff[f]), int(foff[f + 1])
        if b - a < 2:
            continue
        t = pos[b - 1] - pos[a]
        n = np.linalg.norm(t[:2])
        mid = pos[(a + b) // 2, :2]
        rv = mid - fc
        rn = np.linalg.norm(rv)
        if n < 1e-9 or rn < 1e-9 or rn > rmax:
            continue
        ct = (t[:2] @ rv) / (n * rn)                               # cos θ between fibre tangent and radial
        cos2.append(2.0 * ct * ct - 1.0)                          # cos 2θ = 2cos²θ − 1
    return float(np.mean(cos2)) if cos2 else float("nan"), len(cos2)

for rmax in (8.0, 12.0):
    s0, n0 = fiber_nematic_radial(p0, rmax)
    sf, nf = fiber_nematic_radial(pf, rmax)
    print(f"[KB-1.9 nematic S, radial, r<{rmax:.0f}µm]  frame0 S={s0:+.3f} → final S={sf:+.3f}  (ΔS={sf - s0:+.3f}, {nf} fibres)")
print("   (KB-1.9: tumor stroma S=0.3-0.7 aligned; healthy S<0.1; our ΔS = alignment BUILT by the cell over the run)")
