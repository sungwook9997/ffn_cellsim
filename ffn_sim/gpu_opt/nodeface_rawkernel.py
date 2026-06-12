"""SimuCell3D node-vs-face contact as a thread-per-node CUDA RawKernel.

The cupy-array port (``cell/dcm_face_contact.node_face_contact_forces_vec``) is
correct but materialises the full (node,face) candidate array → memory-bound, 36×
slower than node-node at N=100. SimuCell3D is fast because each OpenMP thread owns
one node and loops only the faces co-resident in its voxel, with the closest-point
test in registers and NO materialised candidate array. This module reproduces that:
one CUDA thread per active node, a host-built uniform grid over face centroids
(CSR per-voxel ranges), the Ericson closest-point + bilinear-tent/repulsion law in
registers, and atomicAdd scatter to the face's 3 nodes (Newton-3).

Force law + closest point are BIT-IDENTICAL to ``dcm_face_contact`` (the numpy
ground truth); ``tests``/gbook parity pins it. GPU-only (guarded cupy import).
SI units.
"""

from __future__ import annotations

import numpy as np

_SRC = r"""
extern "C" {

__device__ void closest_pt(
    double px,double py,double pz,
    double ax,double ay,double az, double bx,double by,double bz,
    double cx,double cy,double cz, double* b0,double* b1,double* b2) {
  double abx=bx-ax, aby=by-ay, abz=bz-az;
  double acx=cx-ax, acy=cy-ay, acz=cz-az;
  double apx=px-ax, apy=py-ay, apz=pz-az;
  double d1=abx*apx+aby*apy+abz*apz;
  double d2=acx*apx+acy*apy+acz*apz;
  if (d1<=0.0 && d2<=0.0){*b0=1.0;*b1=0.0;*b2=0.0;return;}
  double bpx=px-bx, bpy=py-by, bpz=pz-bz;
  double d3=abx*bpx+aby*bpy+abz*bpz;
  double d4=acx*bpx+acy*bpy+acz*bpz;
  if (d3>=0.0 && d4<=d3){*b0=0.0;*b1=1.0;*b2=0.0;return;}
  double vc=d1*d4-d3*d2;
  if (vc<=0.0 && d1>=0.0 && d3<=0.0){double v=d1/(d1-d3);*b0=1.0-v;*b1=v;*b2=0.0;return;}
  double cpx=px-cx, cpy=py-cy, cpz=pz-cz;
  double d5=abx*cpx+aby*cpy+abz*cpz;
  double d6=acx*cpx+acy*cpy+acz*cpz;
  if (d6>=0.0 && d5<=d6){*b0=0.0;*b1=0.0;*b2=1.0;return;}
  double vb=d5*d2-d1*d6;
  if (vb<=0.0 && d2>=0.0 && d6<=0.0){double w=d2/(d2-d6);*b0=1.0-w;*b1=0.0;*b2=w;return;}
  double va=d3*d6-d5*d4;
  if (va<=0.0 && (d4-d3)>=0.0 && (d5-d6)>=0.0){
      double w=(d4-d3)/((d4-d3)+(d5-d6));*b0=0.0;*b1=1.0-w;*b2=w;return;}
  double denom=1.0/(va+vb+vc); double v=vb*denom, w=vc*denom;
  *b0=1.0-v-w;*b1=v;*b2=w;
}

__global__ void node_face_contact(
    const double* pos, const long* cof, const long* faces, const long* fcell,
    const long* sorted_fidx, const long* vox_start, const double* origin,
    double vsize, long nx, long ny, long nz,
    double rep, double adh, double c_rep, double c_adh,
    const double* cad, int use_cad, int N, double* F) {
  long ni = (long)blockIdx.x * blockDim.x + threadIdx.x;
  if (ni >= N) return;
  long c1 = cof[ni];
  if (c1 < 0) return;
  double px=pos[3*ni], py=pos[3*ni+1], pz=pos[3*ni+2];
  long vx=(long)floor((px-origin[0])/vsize);
  long vy=(long)floor((py-origin[1])/vsize);
  long vz=(long)floor((pz-origin[2])/vsize);
  double half=0.5*c_adh;
  double Fnx=0.0, Fny=0.0, Fnz=0.0;
  for (int dx=-1; dx<=1; dx++) for (int dy=-1; dy<=1; dy++) for (int dz=-1; dz<=1; dz++) {
    long wx=vx+dx, wy=vy+dy, wz=vz+dz;
    if (wx<0||wx>=nx||wy<0||wy>=ny||wz<0||wz>=nz) continue;
    long vid=(wx*ny+wy)*nz+wz;
    for (long k=vox_start[vid]; k<vox_start[vid+1]; k++) {
      long fi=sorted_fidx[k];
      if (fcell[fi]==c1) continue;
      long a=faces[3*fi], b=faces[3*fi+1], c=faces[3*fi+2];
      double ax=pos[3*a],ay=pos[3*a+1],az=pos[3*a+2];
      double bx=pos[3*b],by=pos[3*b+1],bz=pos[3*b+2];
      double cx=pos[3*c],cy=pos[3*c+1],cz=pos[3*c+2];
      double b0,b1,b2;
      closest_pt(px,py,pz, ax,ay,az, bx,by,bz, cx,cy,cz, &b0,&b1,&b2);
      double cpx=ax*b0+bx*b1+cx*b2, cpy=ay*b0+by*b1+cy*b2, cpz=az*b0+bz*b1+cz*b2;
      double rvx=px-cpx, rvy=py-cpy, rvz=pz-cpz;
      double min_d=sqrt(rvx*rvx+rvy*rvy+rvz*rvz);
      double e1x=bx-ax,e1y=by-ay,e1z=bz-az, e2x=cx-ax,e2y=cy-ay,e2z=cz-az;
      double nfx=e1y*e2z-e1z*e2y, nfy=e1z*e2x-e1x*e2z, nfz=e1x*e2y-e1y*e2x;
      double nnorm=sqrt(nfx*nfx+nfy*nfy+nfz*nfz);
      double area=0.5*nnorm;
      double sign=(nnorm>0.0)?(rvx*nfx+rvy*nfy+rvz*nfz)/nnorm:0.0;
      double amp=0.0;
      if (sign<0.0 && min_d<c_rep) {
        amp=rep*area;
      } else if (adh>0.0 && sign>0.0 && min_d<c_adh) {
        double mult=1.0;
        if (use_cad) mult=sqrt(cad[c1]*cad[fcell[fi]]);
        double md=(min_d>0.0)?min_d:1.0;
        amp=(min_d>=half)?(adh*mult*(c_adh/md-1.0)*area):(adh*mult*area);
      }
      if (amp!=0.0) {
        double Fx=rvx*amp, Fy=rvy*amp, Fz=rvz*amp;
        Fnx-=Fx; Fny-=Fy; Fnz-=Fz;
        atomicAdd(&F[3*a], b0*Fx); atomicAdd(&F[3*a+1], b0*Fy); atomicAdd(&F[3*a+2], b0*Fz);
        atomicAdd(&F[3*b], b1*Fx); atomicAdd(&F[3*b+1], b1*Fy); atomicAdd(&F[3*b+2], b1*Fz);
        atomicAdd(&F[3*c], b2*Fx); atomicAdd(&F[3*c+1], b2*Fy); atomicAdd(&F[3*c+2], b2*Fz);
      }
    }
  }
  atomicAdd(&F[3*ni], Fnx); atomicAdd(&F[3*ni+1], Fny); atomicAdd(&F[3*ni+2], Fnz);
}
}
"""

_MODULE = None
_KERNEL = None


def _kernel():
    global _MODULE, _KERNEL
    if _KERNEL is None:
        import cupy as cp
        _MODULE = cp.RawModule(code=_SRC, options=("--std=c++11",))
        _KERNEL = _MODULE.get_function("node_face_contact")
    return _KERNEL


def node_face_contact_raw(pos, cell_of_node, faces, face_cell, *,
                          rep_strength, adh_strength, c_rep, c_adh, max_edge,
                          cad_mult=None):
    """Thread-per-node node-face contact (cupy). Returns per-node force (N,3)."""
    import cupy as cp
    pos = cp.ascontiguousarray(cp.asarray(pos, cp.float64))
    cof = cp.ascontiguousarray(cp.asarray(cell_of_node, cp.int64))
    faces = cp.ascontiguousarray(cp.asarray(faces, cp.int64))
    fcell = cp.ascontiguousarray(cp.asarray(face_cell, cp.int64))
    N = int(pos.shape[0]); M = int(faces.shape[0])
    F = cp.zeros((N, 3), cp.float64)
    if M < 1 or N < 1:
        return F
    v0 = pos[faces[:, 0]]; v1 = pos[faces[:, 1]]; v2 = pos[faces[:, 2]]
    fc = (v0 + v1 + v2) / 3.0
    r_search = float(max(c_rep, c_adh) + max_edge)
    act = cof >= 0
    origin = cp.minimum(pos[act].min(axis=0), fc.min(axis=0))
    fcl = cp.floor((fc - origin) / r_search).astype(cp.int64)
    ndl = cp.floor((pos[act] - origin) / r_search).astype(cp.int64)
    dims = cp.maximum(fcl.max(axis=0), ndl.max(axis=0)) + 1
    nx, ny, nz = int(dims[0]), int(dims[1]), int(dims[2])
    vid = (fcl[:, 0] * ny + fcl[:, 1]) * nz + fcl[:, 2]
    order = cp.argsort(vid).astype(cp.int64)
    n_vox = nx * ny * nz
    counts = cp.bincount(vid, minlength=n_vox)
    vox_start = cp.zeros(n_vox + 1, cp.int64)
    vox_start[1:] = cp.cumsum(counts)
    use_cad = 1 if cad_mult is not None else 0
    cad = (cp.asarray(cad_mult, cp.float64) if cad_mult is not None
           else cp.zeros(1, cp.float64))
    origin = cp.ascontiguousarray(origin.astype(cp.float64))
    threads = 128
    blocks = (N + threads - 1) // threads
    _kernel()((blocks,), (threads,), (
        pos, cof, faces, fcell, order, vox_start, origin,
        cp.float64(r_search), cp.int64(nx), cp.int64(ny), cp.int64(nz),
        cp.float64(rep_strength), cp.float64(adh_strength),
        cp.float64(c_rep), cp.float64(c_adh), cad, cp.int32(use_cad),
        cp.int32(N), F))
    return F
