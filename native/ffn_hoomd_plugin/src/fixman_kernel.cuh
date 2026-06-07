#pragma once

#include "hip/hip_runtime.h"

namespace ffn_native
    {
// Native Fixman metric pseudo-force — Stage 1b port of
// ``ffn_sim/integrator/constrained_baoab.py::fixman_logdet_and_force`` (uniform
// fast path). ONE THREAD PER CHAIN: builds the m×m tridiagonal mobility-metric
// Gram matrix G, computes ln det G + the needed entries of G^{-1} via the
// tridiagonal continuant (θ/φ) recurrences, and the analytic ∇ln det G force.
// Per-chain independent → no cross-thread reduction (numerically matches the
// Python einsum/linalg path to tight tolerance; LAPACK vs continuant differ in
// the last ULPs, so parity is tolerance-based, not bitwise — acceptable per the
// gates-based acceptance, and the force is what feeds the dynamics).
//
//   pos        : (N,3) row-major double
//   inv_gamma  : (N,) double (per-bead mobility, = inv_mass convention)
//   chains     : (F, m+1) int32 bead row indices (uniform length)
//   force      : (N,3) out (Fixman force, += per disjoint chain)
//   logdet     : (F,) out (per-chain ln|det G|; U_F = 0.5 kT Σ logdet on host)
//   bad_sign   : (F,) out (1 if det G <= 0 → metric non-positive)
constexpr unsigned int FFN_FIXMAN_MMAX = 64;

// Device-pointer launcher (no host copy). force is OVERWRITTEN for chain beads
// (caller must zero non-chain beads if needed). logdet/bad_sign device arrays.
hipError_t gpu_fixman_launch(const double* d_pos,
                             const double* d_inv_gamma,
                             const int* d_chains,
                             double* d_force,
                             double* d_logdet,
                             int* d_bad_sign,
                             unsigned int F,
                             unsigned int m,
                             double Lx,
                             double Ly,
                             double Lz);

hipError_t fixman_force_host(const double* pos,
                             const double* inv_gamma,
                             const int* chains,
                             double* force,
                             double* logdet,
                             int* bad_sign,
                             unsigned int N,
                             unsigned int F,
                             unsigned int m,
                             double Lx,
                             double Ly,
                             double Lz);
    } // namespace ffn_native
