import numpy as np, sys
def cell_vol_area(P,faces,fcell,c):
    fm=fcell==c; f=faces[fm]
    v=P[f]; 
    # signed volume via divergence (sum tri (a x b)·c /6), area = sum tri area
    a,b,cc=v[:,0],v[:,1],v[:,2]
    vol=np.abs(np.einsum('ij,ij->i',np.cross(a,b),cc).sum()/6.0)
    ar=0.5*np.linalg.norm(np.cross(b-a,cc-a),axis=1).sum()
    return vol,ar
def analyze(npz,label,frame=0):
    d=np.load(npz,allow_pickle=True)
    P=d["frames"][frame].astype(np.float64); faces=d["faces"]; cof=d["cof"]
    fcell=cof[faces[:,0]]
    nc=int(cof[cof>=0].max())+1; R=7.5e-6; UM=1e6
    cen=np.array([P[cof==c].mean(0) for c in range(nc)])
    # NN center spacing
    from scipy.spatial import cKDTree
    tr=cKDTree(cen); dd,_=tr.query(cen,k=2); nn=dd[:,1]/R
    # coordination: neighbors within 2.2R
    z=np.array([len(tr.query_ball_point(cen[i],2.2*R))-1 for i in range(nc)])
    # per-cell sphericity Psi = pi^(1/3)(6V)^(2/3)/A  (=1 sphere, <1 flattened/polyhedral)
    psi=[]; 
    for c in range(nc):
        V,A=cell_vol_area(P,faces,fcell,c)
        psi.append(np.pi**(1/3)*(6*V)**(2/3)/A)
    psi=np.array(psi)
    # contact fraction: per cell, frac of its nodes within c_adh(1.8um) of a DIFFERENT cell's node
    c_adh=1.8e-6
    allt=cKDTree(P[cof>=0]); idx=np.flatnonzero(cof>=0); cofl=cof[idx]
    cfrac=[]
    for c in range(nc):
        nm=np.flatnonzero(cof==c); 
        nbrs=allt.query_ball_point(P[nm],c_adh)
        incontact=[any(cofl[j]!=c for j in nb) for nb in nbrs]
        cfrac.append(np.mean(incontact))
    cfrac=np.array(cfrac)
    print(f"\n===== {label} (frame {frame}) =====")
    print(f"  cells={nc}  NN center spacing: mean={nn.mean():.2f}R  min={nn.min():.2f}R  (gap built=2.05R)")
    print(f"  coordination z (nbrs<2.2R): mean={z.mean():.1f}  (FCC bulk=12; SimuCell3D dense z~12-14)")
    print(f"  cell sphericity Psi: mean={psi.mean():.3f}  (1.0=sphere/no-flatten; <1=flattened tissue)")
    print(f"  contact fraction phi: mean={cfrac.mean():.2f}  (frac of surface near another cell; dense tissue→high)")
    return dict(nn=nn,z=z,psi=psi,cfrac=cfrac)
analyze("aleph/outputs/warp_decohesion/n100_implicit_decisive.npz","DECISIVE (OLD weak cohesion 0.18nN/junction)",0)
