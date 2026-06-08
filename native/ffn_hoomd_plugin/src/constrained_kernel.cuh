#pragma once

#include <hoomd/HOOMDMath.h>

#include "hip/hip_runtime.h"

namespace ffn_native
    {
// Device-resident constrained Leimkuhler-Matthews BAOAB step (Stage 1b assembly).
// Orchestrates, all on the GPU with NO host round-trip, the validated kernels in
// constrained_baoab.act() order: extract → map chains tag→row → Fixman force →
// L-M predictor (no wrap) → M-SHAKE projection → wrap+image → store prv.
//
// HOOMD arrays (Scalar4 pos/force, int3 image, uint tag) + persistent device
// buffers (inv_gamma_by_tag, bd_pref_by_tag, prv[3N], chains_tag[F*(m+1)]) +
// scratch (pos_d/ref_d/ffix_d/pred_d[3N], inv_gamma_row[N], row_of_tag[N],
// chains_row[F*(m+1)], lambda[F*m], nonconv[F], bad_sign[F], logdet[F]).
hipError_t gpu_constrained_step(hoomd::Scalar4* d_pos,
                                const hoomd::Scalar4* d_force,
                                int3* d_image,
                                const unsigned int* d_tag,
                                const double* d_inv_gamma_by_tag,
                                const double* d_bd_pref_by_tag,
                                double* d_prv,
                                const int* d_chains_tag,
                                double* d_pos_d,
                                double* d_ref_d,
                                double* d_ffix_d,
                                double* d_pred_d,
                                double* d_inv_gamma_row,
                                int* d_row_of_tag,
                                int* d_chains_row,
                                double* d_lambda,
                                int* d_nonconv,
                                int* d_bad_sign,
                                double* d_logdet,
                                unsigned int N,
                                unsigned int F,
                                unsigned int m,
                                unsigned long long seed,
                                unsigned long long timestep,
                                double dt,
                                double half_kT,
                                double rest_length,
                                double Lx,
                                double Ly,
                                double Lz,
                                double xy,
                                double xz,
                                double yz,
                                double tol,
                                unsigned int max_iter,
                                unsigned int block_size,
                                int compression_release,
                                double Fcrit,
                                double load_alpha,
                                double* d_load_ema);
    } // namespace ffn_native
