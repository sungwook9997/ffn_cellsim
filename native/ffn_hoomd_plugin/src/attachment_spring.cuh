#pragma once

#include <hoomd/HOOMDMath.h>

#include "hip/hip_runtime.h"

namespace ffn_native
    {
// Native attachment-spring force (fixed-pool binder enabler). The dynamic
// myosin/xlink head<->actin springs are a FIXED pool of K slots, each a tuple
// (head_tag, actin_tag, k, r0) in device arrays the binder updater toggles
// (k=0 / head_tag<0 = inactive). This ForceCompute reads those arrays + the
// positions via ArrayHandle and sums the harmonic springs — so binding is a
// per-slot array write (gpu_local), NEVER a HOOMD bond mutation / global
// set_snapshot (the 51.9 ms/firing wall measured in h7_native_gate_a_profile).
//
// d_row_of_tag (N int, scratch) maps tag->row each step (ParticleSorter-safe).
hipError_t gpu_attachment_spring(const hoomd::Scalar4* d_pos,
                                 const unsigned int* d_tag,
                                 hoomd::Scalar4* d_force,
                                 int* d_row_of_tag,
                                 const int* d_head_tag,
                                 const int* d_actin_tag,
                                 const double* d_k,
                                 const double* d_r0,
                                 unsigned int N,
                                 unsigned int K,
                                 double Lx,
                                 double Ly,
                                 double Lz,
                                 unsigned int block_size);
    } // namespace ffn_native
