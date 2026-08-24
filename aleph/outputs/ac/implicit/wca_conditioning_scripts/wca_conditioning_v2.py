"""WCA steric ill-conditioning diagnostic v2 — pure NumPy, CPU-only (analysis, NOT a Warp run).

Corrects v1: (1) probe the UNCAPPED overlap band where the WCA tangent is actually assembled;
(2) demonstrate the force-cap residual pathology (deep overlap -> standing f_cap force with ZERO
tangent); (3) randomize node phase in the coverage test. Uses the repo wca_analytic oracle so all
magnitudes match production. No CUDA, no simulation runtime.
"""
from __future__ import annotations
import numpy as np
from ffn_sim.ac.solid.wca_analytic import (
    epsilon_from_contact_stiffness, contact_stiffness, wca_cutoff,
    wca_force_magnitude, wca_local_stiffness, max_core_stiffness, compressed_pair_separation,
)
EPS64 = np.finfo(np.float64).eps
SQRT_EPS = np.sqrt(EPS64)
FLOOR = SQRT_EPS
SIGMA, K_EV, F_CAP, ELL = 0.007, 1.0e3, 1.0e3, 0.5
EI = 4.114e-3 * 17.7
ALPHA_BEND = EI / ELL**3
eps = epsilon_from_contact_stiffness(K_EV, SIGMA)
rc = wca_cutoff(SIGMA)
r_eq_cap = compressed_pair_separation(F_CAP, SIGMA, eps)   # boundary of the capped region
kmax_steric = max_core_stiffness(F_CAP, SIGMA, eps)

# ---------------------------------------------------------------------------- model machinery
def make_fiber(o, ax, M, ell):
    ax = np.asarray(ax, float); ax = ax/np.linalg.norm(ax)
    return np.array([np.asarray(o, float)+i*ell*ax for i in range(M)])
def bend_triples(base, M): return [(base+i, base+i+1, base+i+2) for i in range(M-2)]
def cjac(pos, segs):
    N = pos.shape[0]; J = np.zeros((len(segs), 3*N))
    for r,(i,j) in enumerate(segs):
        u = pos[j]-pos[i]; u/=np.linalg.norm(u); J[r,3*i:3*i+3]=-u; J[r,3*j:3*j+3]=u
    return J
def proj(J): return np.eye(J.shape[1]) - J.T@np.linalg.solve(J@J.T, J)
def Kbend(pos, tri, a):
    N=pos.shape[0]; K=np.zeros((3*N,3*N)); I3=np.eye(3)
    for (p,q,r) in tri:
        w={p:1.,q:-2.,r:1.}
        for u,wu in w.items():
            for v,wv in w.items(): K[3*u:3*u+3,3*v:3*v+3]+=a*wu*wv*I3
    return K
def Kwca(pos, pairs, sigma, eps, f_cap):
    N=pos.shape[0]; K=np.zeros((3*N,3*N)); added=0
    for (i,j) in pairs:
        d=pos[i]-pos[j]; L=np.linalg.norm(d)
        if L<=1e-12 or L>=wca_cutoff(sigma): continue
        if f_cap>0 and float(wca_force_magnitude(L,sigma,eps))>=f_cap: continue  # capped -> zero tangent
        k=float(wca_local_stiffness(L,sigma,eps)); u=(d/L).reshape(3,1); blk=k*(u@u.T)
        K[3*i:3*i+3,3*i:3*i+3]+=blk; K[3*j:3*j+3,3*j:3*j+3]+=blk
        K[3*i:3*i+3,3*j:3*j+3]-=blk; K[3*j:3*j+3,3*i:3*i+3]-=blk; added+=1
    return K, added
def range_spec(A, P):
    w,V=np.linalg.eigh(P); Q=V[:,w>1e-9]; return np.linalg.eigvalsh(Q.T@A@Q)
def pcg(A,b,Minv,tol,maxit):
    x=np.zeros_like(b); r=b-A@x; z=Minv(r); p=z.copy(); rz=r@z; b0=np.linalg.norm(b)
    for k in range(1,maxit+1):
        Ap=A@p; al=rz/(p@Ap); x+=al*p; r-=al*Ap
        if np.linalg.norm(r)/b0<=tol: return k,np.linalg.norm(r)/b0
        z=Minv(r); rz2=r@z; p=z+(rz2/rz)*p; rz=rz2
    return maxit,np.linalg.norm(r)/b0
def jac(A):
    d=np.diag(A).copy(); d[d<=0]=1.; return lambda r:r/d
def coarse(A, fibers, base):
    N=A.shape[0]//3; Z=np.zeros((3*N,3*len(fibers)))
    for f,ns in enumerate(fibers):
        for dd in range(3):
            for n in ns: Z[3*n+dd,3*f+dd]=1.
    Ac_inv=np.linalg.pinv(Z.T@A@Z)
    return lambda r: base(r)+Z@(Ac_inv@(Z.T@r))

def build_cross(M, gap, phaseA=0.0, phaseB=0.0):
    """Two perpendicular fibers crossing with centerline gap 'gap'; node phase shifts crossing off-node."""
    fa = make_fiber([-ELL*(M//2)+phaseA,0,0],[1,0,0],M,ELL)
    fb = make_fiber([0,-ELL*(M//2)+phaseB,gap],[0,1,0],M,ELL)
    pos=np.vstack([fa,fb])
    tri=bend_triples(0,M)+bend_triples(M,M)
    seg=[(i,i+1) for i in range(M-1)]+[(M+i,M+i+1) for i in range(M-1)]
    fib=[list(range(M)),list(range(M,2*M))]
    return pos,tri,seg,fib

print("="*98)
print("PART A — closed-form: k_max, the force-cap boundary, and the kappa = 1/sqrt(eps64) floor")
print("="*98)
print(f"sigma=7nm  r_c={rc*1e3:.4f}nm  eps={eps:.4e}pN.um  k_ev=U''(r_min)={contact_stiffness(eps,SIGMA):.3e}")
print(f"force-cap boundary r_eq(f_cap)= {r_eq_cap*1e3:.4f}nm  (overlap {100*(rc-r_eq_cap)/rc:.1f}% of r_c)")
print(f"  -> overlaps DEEPER than this are CAPPED: force pinned at {F_CAP:.0f}pN, tangent = 0")
print(f"  -> overlaps SHALLOWER carry a live tangent up to k_max=U''(r_eq)= {kmax_steric:.4e} pN/um")
print(f"bending alpha=EI/l^3={ALPHA_BEND:.3f}  ->  WCA/bending stiffness contrast = {kmax_steric/ALPHA_BEND:.2e}")
print(f"regularizer a = sqrt(eps64)*kmax; if lambda_min sits at a and lambda_max=kmax then")
print(f"   kappa = kmax/a = 1/sqrt(eps64) = {1/SQRT_EPS:.3e}  (a HARD floor, cancels in k_ev)")

print("\n"+"="*98)
print("PART B1 — spectrum in the UNCAPPED band: lambda_max tracks the live WCA tangent (Hyp A)")
print("="*98)
M=7; CP=[(M//2,M+M//2)]
print(f"{'overlap%':>8} {'r[nm]':>8} {'Fpair[pN]':>10} {'capped?':>7} {'k_ax':>10} "
      f"{'lam_min':>10} {'lam_max':>11} {'kappa':>10} {'CGjac':>6} {'CGcrs':>6}")
for frac in [0.02,0.10,0.20,0.30,0.39,0.60,0.95]:
    gap=rc*(1-frac); pos,tri,seg,fib=build_cross(M,gap)
    r=np.linalg.norm(pos[CP[0][0]]-pos[CP[0][1]])
    Fp=float(wca_force_magnitude(r,SIGMA,eps)); capped=(Fp>=F_CAP)
    Kc,_=Kwca(pos,CP,SIGMA,eps,F_CAP)
    P=proj(cjac(pos,seg)); K=Kbend(pos,tri,ALPHA_BEND)+Kc
    kmax=max(np.max(np.abs(np.diag(K))),ALPHA_BEND); a=FLOOR*kmax
    A=a*np.eye(P.shape[0])+P@K@P
    lam=range_spec(A,P); rng=np.random.default_rng(0); b=P@rng.standard_normal(A.shape[0])
    ij,_=pcg(A,b,jac(A),SQRT_EPS,5000); ic,_=pcg(A,b,coarse(A,fib,jac(A)),SQRT_EPS,5000)
    kax=float(wca_local_stiffness(r,SIGMA,eps))
    print(f"{frac*100:>7.0f}% {r*1e3:>8.3f} {Fp:>10.2e} {str(capped):>7} {kax:>10.2e} "
          f"{lam[0]:>10.2e} {lam[-1]:>11.3e} {lam[-1]/lam[0]:>10.2e} {ij:>6} {ic:>6}")

print("\n"+"="*98)
print("PART B2 — force-cap RESIDUAL pathology: deep overlap -> standing force, ZERO tangent (the 2204 pN)")
print("="*98)
# a node with a few deep (capped) cross-fiber overlaps: net force is O(f_cap), tangent is exactly zero
for ncontact,spread in [(1,0.0),(2,0.3),(3,0.5)]:
    r_deep=r_eq_cap*0.5  # well inside the capped region
    # forces from ncontact capped neighbours, directions spread by 'spread' rad from +z
    net=np.zeros(3);
    for c in range(ncontact):
        ang=spread*(c-(ncontact-1)/2)
        nvec=np.array([np.sin(ang),0,np.cos(ang)]); net+=F_CAP*nvec
    tangent=0.0  # capped branch contributes exactly zero stiffness (kernel + oracle)
    print(f"  {ncontact} capped contacts, spread {spread:.1f}rad: |net force|={np.linalg.norm(net):7.1f} pN  "
          f"assembled tangent={tangent:.1f} pN/um  -> Newton correction in this dir = force/tangent = INF/undefined")
print("  => the solver has a standing residual it CANNOT reduce: no curvature where the force lives.")
print(f"  => matches REPORT: WCA-on resting max projected residual 2204 pN (~2x f_cap) vs WCA-off 21 pN.")

print("\n"+"="*98)
print("PART B3 — is lambda_min (Hyp B) a pre-existing rigid-mode floor, independent of WCA?")
print("="*98)
gapU=rc*0.39  # uncapped, near-peak tangent
pos,tri,seg,fib=build_cross(M,gapU)
for lbl,withc,fixed_a in [("bending only",False,None),("bending+WCA(uncapped)",True,None),
                          ("bending+WCA, a fixed=1.0 (decoupled reg)",True,1.0)]:
    Kc=Kwca(pos,CP,SIGMA,eps,F_CAP)[0] if withc else np.zeros((6*M,6*M))
    P=proj(cjac(pos,seg)); K=Kbend(pos,tri,ALPHA_BEND)+Kc
    kmax=max(np.max(np.abs(np.diag(K))),ALPHA_BEND); a=fixed_a if fixed_a else FLOOR*kmax
    A=a*np.eye(P.shape[0])+P@K@P; lam=range_spec(A,P)
    print(f"  {lbl:42} a={a:.2e} lam_min={lam[0]:.2e} lam_max={lam[-1]:.2e} kappa={lam[-1]/lam[0]:.2e}")
print("  => lambda_min is the near-null fiber mode at the sqrt(eps)*kmax floor: kappa~1/sqrt(eps) with or")
print("     without WCA. Decoupling the regularizer (fixed a) collapses kappa -> the floor is the lever.")

print("\n"+"="*98)
print("PART C1 — near-null mode COUNT grows with fiber count (why native CG stalls; coarse-space need)")
print("="*98)
for nf in [2,4,8,16]:
    poss=[]; tris=[]; segs=[]; fibs=[]; off=0
    rng=np.random.default_rng(nf)
    for k in range(nf):
        o=rng.uniform(-2,2,3); ax=rng.standard_normal(3)
        f=make_fiber(o,ax,M,ELL); poss.append(f)
        tris+=bend_triples(off,M); segs+=[(off+i,off+i+1) for i in range(M-1)]
        fibs.append(list(range(off,off+M))); off+=M
    pos=np.vstack(poss); P=proj(cjac(pos,segs)); K=Kbend(pos,tris,ALPHA_BEND)
    kmax=max(np.max(np.abs(np.diag(K))),ALPHA_BEND); a=FLOOR*kmax
    A=a*np.eye(P.shape[0])+P@K@P; lam=range_spec(A,P)
    nnull=int(np.sum(lam < 1e3*a))   # modes near the regularizer floor
    b=P@rng.standard_normal(A.shape[0])
    ij,_=pcg(A,b,jac(A),SQRT_EPS,20000); ic,_=pcg(A,b,coarse(A,fibs,jac(A)),SQRT_EPS,20000)
    print(f"  {nf:2d} fibers ({nf*M} nodes): near-floor modes={nnull:3d}  kappa={lam[-1]/lam[0]:.2e}  "
          f"CG(jac)={ij:5d}  CG(+coarse)={ic:4d}")
print("  => #near-floor modes ~ O(#fibers); Jacobi-CG iters scale with them; the per-fiber coarse space")
print("     deflates the translation part. Native = 70,686 fibers -> ~1e5 near-null modes -> CG768 stalls.")

print("\n"+"="*98)
print("PART C2 — node-node vs seg-seg: coverage of real crossings (randomized node phase)")
print("="*98)
def segseg_d(p0,p1,q0,q1):
    d1=p1-p0; d2=q1-q0; r=p0-q0; a=d1@d1; e=d2@d2; f=d2@r; c=d1@r; b=d1@d2
    den=a*e-b*b; s=0. if den<1e-18 else np.clip((b*f-c*e)/den,0,1)
    t=np.clip((b*s+f)/e,0,1) if e>1e-18 else 0.; s=np.clip((b*t-c)/a,0,1) if a>1e-18 else 0.
    return np.linalg.norm((p0+s*d1)-(q0+t*d2))
rng=np.random.default_rng(11)
for sg,lbl in [(0.007,"physical 7nm"),(0.025,"MT 25nm"),(0.05,"coarse 50nm")]:
    rct=wca_cutoff(sg); nn=ss=0
    for _ in range(6000):
        c=rng.uniform(-3,3,3); a1=rng.standard_normal(3); a1/=np.linalg.norm(a1)
        a2=rng.standard_normal(3); a2/=np.linalg.norm(a2)
        perp=rng.standard_normal(3); perp-=(perp@a1)*a1; perp/=np.linalg.norm(perp)
        gap=rng.uniform(0,rct)                       # centerline gap < r_c => a REAL crossing
        pa=rng.uniform(-ELL/2,ELL/2); pb=rng.uniform(-ELL/2,ELL/2)   # random node phase
        fa=make_fiber(c-3*ELL*a1+pa*a1,a1,M,ELL)
        fb=make_fiber(c+gap*perp-3*ELL*a2+pb*a2,a2,M,ELL)
        sa=[(fa[i],fa[i+1]) for i in range(M-1)]; sb=[(fb[i],fb[i+1]) for i in range(M-1)]
        dmin=min(segseg_d(*x,*y) for x in sa for y in sb)
        if dmin>=rct: continue
        ss+=1
        ndm=min(np.linalg.norm(x-y) for x in fa for y in fb)
        if ndm<rct: nn+=1
    print(f"  sigma={lbl:13} r_c={rct*1e3:7.2f}nm : real crossings={ss:4d}  node-node HITS={nn:4d}  "
          f"coverage={100*nn/max(ss,1):5.1f}%")
print("  => at physical sigma, node-node MISSES essentially all real crossings (nodes 0.5um apart, r_c~8nm);")
print("     the contacts it DOES fire are accidental node coincidences (discretization artifacts).")

print("\n"+"="*98)
print("PART C3 — peak nodal stiffness & k_max: seg-seg spreads load; larger effective sigma lowers k_max")
print("="*98)
print("  node-node: rank-1 on ONE pair, peak nodal = k_max.  seg-seg: split over 4 endpoints (w<=1).")
for sg,mode in [(0.007,"nn"),(0.05,"nn"),(0.007,"ss"),(0.05,"ss")]:
    ept=epsilon_from_contact_stiffness(K_EV,sg); kmx=max_core_stiffness(F_CAP,sg,ept)
    if mode=="ss":
        peak=kmx*0.5**2; lbl=f"seg-seg sigma={sg*1e3:.0f}nm"
    else:
        peak=kmx; lbl=f"node-node sigma={sg*1e3:.0f}nm"
    print(f"  {lbl:22} k_max={kmx:.3e}  peak_nodal_stiffness={peak:.3e}  explicit dt_mu=0.1/kmax={0.1/kmx:.2e}s")
