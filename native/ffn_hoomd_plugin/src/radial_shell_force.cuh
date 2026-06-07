#pragma once

#include <hoomd/HOOMDMath.h>

#include "hip/hip_runtime.h"

namespace ffn_native
    {
// Native radial-shell compartment force (turgor / membrane / nucleus). All 3
// production compartment md.force.Custom share one shape: a global shell
// centroid + mean radius (reductions over the tagged shell beads), then a
// per-bead radial force from a law. Computed fully on-device via ArrayHandle
// (NO gpu_local_snapshot — removes the ~233 us/force HOOMD-API floor the cupy
// port could not). Law selector:
//   0 = nucleus  : d=r-R0 bilinear (a=k_chrom, b=k_lamin, c=d_knee, d=F_knee)
//   1 = membrane : S=4piR^2, gamma=a+b*(S-c)/c, dP=2*gamma/R, F=-dP*(S/n)*nhat  (a=gamma_mem,b=K_A,c=A0)
//   2 = turgor   : V=4/3 pi R^3, dP=a-b*(V-c)/c, F=+dP*(S/n)*nhat              (a=dP0,b=K_vol,c=V0)
//
// d_acc layout (5 Scalars): [sum_x, sum_y, sum_z, count, sum_r].
hipError_t gpu_radial_shell_force(const hoomd::Scalar4* d_pos,
                                  const unsigned int* d_tag,
                                  hoomd::Scalar4* d_force,
                                  hoomd::Scalar* d_acc,
                                  unsigned int N,
                                  unsigned int tag_start,
                                  unsigned int tag_end,
                                  int law,
                                  double R0,
                                  double pa,
                                  double pb,
                                  double pc,
                                  double pd,
                                  unsigned int block_size);
    } // namespace ffn_native
