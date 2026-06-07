#pragma once

#include <hoomd/HOOMDMath.h>

#include "hip/hip_runtime.h"

namespace ffn_native
    {
// Native overdamped Leimkuhler-Matthews BAOAB-limit step (Stage 1a).
//
// Bit-for-bit port of the validated cupy RawKernel
// (``ffn_sim/integrator/baoab_device.py::baoab_step``): same splitmix64 +
// Box-Muller counter-based RNG keyed by (seed, timestep, tag), same two-Gaussian
// noise ``pref*(W_n + W_{n-1})*dt``, same triclinic fractional-coordinate wrap.
// Indexing:
//   * d_pos / d_force are HOOMD Scalar4 in CURRENT index order (.w of pos = type,
//     preserved on write); d_image is int3 in index order.
//   * d_tag[i] is the stable tag of the particle at index i.
//   * d_gamma_by_tag / d_bd_prefactor_by_tag / d_prv are indexed by TAG (the
//     dense [0,N) tag space OverdampedBAOABDevice requires), so per-particle
//     state follows the particle across HOOMD's ParticleSorter reordering.
hipError_t gpu_baoab_step(hoomd::Scalar4* d_pos,
                          const hoomd::Scalar4* d_force,
                          int3* d_image,
                          const unsigned int* d_tag,
                          const hoomd::Scalar* d_gamma_by_tag,
                          const hoomd::Scalar* d_bd_prefactor_by_tag,
                          hoomd::Scalar* d_prv,
                          unsigned int n_particles,
                          unsigned long long seed,
                          unsigned long long timestep,
                          double dt,
                          double Lx,
                          double Ly,
                          double Lz,
                          double xy,
                          double xz,
                          double yz,
                          unsigned int block_size);
    } // namespace ffn_native
