#pragma once

#include "hip/hip_runtime.h"

namespace ffn_native
    {
// Native Matrix-SHAKE (tridiagonal Newton) for linear-chain bond constraints —
// Stage 1b constraint-projection port of
// ``ffn_sim/integrator/constrained_baoab.py::shake_project_chains`` (uniform
// fast path). ONE THREAD PER CHAIN (chains are disjoint → no inter-thread race,
// per-chain serial Newton + Thomas in registers/local memory), so the float op
// order matches the Python einsum/Thomas exactly → bit-identical.
//
// Host-pointer wrapper: allocates device buffers, copies H2D, launches, copies
// D2H. ``pos`` is updated in place (projected positions); ``lambda`` receives
// the accumulated per-bond Lagrange multipliers (F*m). Orthorhombic box only
// (the cortex constraint path; matches ``_min_image_orthorhombic``).
//
//   pos / ref     : (N, 3) row-major double
//   inv_mass      : (N,) double
//   chains        : (F, m+1) int32 bead row indices (uniform length)
//   lambda        : (F, m) double  out
//   returns hipSuccess, or the launch error.
constexpr unsigned int FFN_SHAKE_MMAX = 64; // max bonds per chain

// Device-pointer launcher (no host copy) — for the device-resident constrained
// updater. pos/ref/inv_mass/chains/lambda/nonconv are ALL device pointers.
hipError_t gpu_shake_launch(double* d_pos,
                            const double* d_ref,
                            const double* d_inv_mass,
                            const int* d_chains,
                            double* d_lambda,
                            int* d_nonconv,
                            unsigned int F,
                            unsigned int m,
                            double rest_length,
                            double Lx,
                            double Ly,
                            double Lz,
                            double tol,
                            unsigned int max_iter);

// Relaxed (unilateral) M-SHAKE with the τ_bend-EMA Euler buckling gate — the H.7
// Gate-B port (see shake_kernel.cu). ``d_load_ema`` is a persistent (F*m) device
// array (the per-bond sustained-load EMA, zero-initialised once); ``alpha`` =
// dt/τ_bend; ``Fcrit`` is the Euler buckling load π²κ/ℓ₀². All pointers device.
hipError_t gpu_shake_relaxed_launch(double* d_pos,
                                    const double* d_ref,
                                    const double* d_inv_mass,
                                    const int* d_chains,
                                    double* d_lambda,
                                    int* d_nonconv,
                                    double* d_load_ema,
                                    unsigned int F,
                                    unsigned int m,
                                    double rest_length,
                                    double Lx,
                                    double Ly,
                                    double Lz,
                                    double tol,
                                    unsigned int max_iter,
                                    double Fcrit,
                                    double dt,
                                    double alpha);

hipError_t shake_project_host(double* pos,
                              const double* ref,
                              const double* inv_mass,
                              const int* chains,
                              double* lambda,
                              unsigned int N,
                              unsigned int F,
                              unsigned int m,
                              double rest_length,
                              double Lx,
                              double Ly,
                              double Lz,
                              double tol,
                              unsigned int max_iter,
                              int* out_nonconverged);
    } // namespace ffn_native
