#include "constrained_kernel.cuh"

#ifdef ENABLE_HIP

#include "fixman_kernel.cuh"
#include "shake_kernel.cuh"

namespace ffn_native
    {

// --- counter RNG (identical to baoab_kernel / baoab_device.baoab_step) ------ //
__device__ inline unsigned long long ck_splitmix64(unsigned long long x)
    {
    x += 0x9E3779B97F4A7C15ULL;
    x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ULL;
    x = (x ^ (x >> 27)) * 0x94D049BB133111EBULL;
    return x ^ (x >> 31);
    }
__device__ inline double ck_uniform01(unsigned long long x)
    {
    const unsigned long long r = ck_splitmix64(x);
    return ((double)(r >> 11) + 0.5) * 1.1102230246251565e-16;
    }
__device__ inline void ck_normal_pair(unsigned long long base, double* z0, double* z1)
    {
    double u1 = ck_uniform01(base);
    const double u2 = ck_uniform01(base ^ 0xD1B54A32D192ED03ULL);
    if (u1 < 1.0e-300)
        u1 = 1.0e-300;
    const double r = sqrt(-2.0 * log(u1));
    const double theta = 6.283185307179586476925286766559 * u2;
    *z0 = r * cos(theta);
    *z1 = r * sin(theta);
    }

// 1. Scalar4 pos → row-indexed double pos_d (=ref_d); row_of_tag; inv_gamma_row.
__global__ void ck_extract(const hoomd::Scalar4* __restrict__ pos,
                           const unsigned int* __restrict__ tag,
                           const double* __restrict__ inv_gamma_by_tag,
                           double* __restrict__ pos_d,
                           double* __restrict__ ref_d,
                           double* __restrict__ inv_gamma_row,
                           int* __restrict__ row_of_tag,
                           unsigned int N)
    {
    const unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= N)
        return;
    const hoomd::Scalar4 p = pos[i];
    pos_d[3 * i + 0] = (double)p.x;
    pos_d[3 * i + 1] = (double)p.y;
    pos_d[3 * i + 2] = (double)p.z;
    ref_d[3 * i + 0] = (double)p.x;
    ref_d[3 * i + 1] = (double)p.y;
    ref_d[3 * i + 2] = (double)p.z;
    const unsigned int t = tag[i];
    inv_gamma_row[i] = inv_gamma_by_tag[t];
    row_of_tag[t] = (int)i;
    }

// 2. chains_tag → chains_row via row_of_tag.
__global__ void ck_map_chains(const int* __restrict__ chains_tag,
                              const int* __restrict__ row_of_tag,
                              int* __restrict__ chains_row,
                              unsigned int total)
    {
    const unsigned int k = blockIdx.x * blockDim.x + threadIdx.x;
    if (k >= total)
        return;
    chains_row[k] = row_of_tag[chains_tag[k]];
    }

// 4. L-M predictor with Fixman force (no wrap). pred_d = pos_d + (F+Ffix)/γ dt
//    + bd_pref (W_n + prv) dt ; store prv[tag] = W_n.
__global__ void ck_predict(const double* __restrict__ pos_d,
                           const hoomd::Scalar4* __restrict__ force,
                           const double* __restrict__ ffix,
                           const unsigned int* __restrict__ tag,
                           const double* __restrict__ inv_gamma_row,
                           const double* __restrict__ bd_pref_by_tag,
                           double* __restrict__ prv,
                           double* __restrict__ pred_d,
                           unsigned int N,
                           unsigned long long seed,
                           unsigned long long timestep,
                           double dt,
                           double half_kT)
    {
    const unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= N)
        return;
    const unsigned int t = tag[i];
    const unsigned long long tp = 3ULL * (unsigned long long)t;
    const double ig = inv_gamma_row[i];
    const double pref = bd_pref_by_tag[t];
    const unsigned long long rng_base = seed ^ (timestep * 0xD2B74407B1CE6E93ULL)
                                        ^ (((unsigned long long)t + 1ULL) * 0xCA5A826395121157ULL);
    double W0, W1, W2, Wu;
    ck_normal_pair(rng_base, &W0, &W1);
    ck_normal_pair(rng_base ^ 0x9E3779B97F4A7C15ULL, &W2, &Wu);

    // Fixman force = half_kT * (d ln det / d r); the kernel wrote the raw grad.
    const hoomd::Scalar4 f = force[i];
    pred_d[3 * i + 0] = pos_d[3 * i + 0] + ((double)f.x + half_kT * ffix[3 * i + 0]) * ig * dt
                        + pref * (W0 + prv[tp + 0]) * dt;
    pred_d[3 * i + 1] = pos_d[3 * i + 1] + ((double)f.y + half_kT * ffix[3 * i + 1]) * ig * dt
                        + pref * (W1 + prv[tp + 1]) * dt;
    pred_d[3 * i + 2] = pos_d[3 * i + 2] + ((double)f.z + half_kT * ffix[3 * i + 2]) * ig * dt
                        + pref * (W2 + prv[tp + 2]) * dt;
    prv[tp + 0] = W0;
    prv[tp + 1] = W1;
    prv[tp + 2] = W2;
    }

// 6. wrap projected pred_d into the (triclinic) box → Scalar4 pos + image delta.
__global__ void ck_wrap(const double* __restrict__ pred_d,
                        hoomd::Scalar4* __restrict__ pos,
                        int3* __restrict__ image,
                        unsigned int N,
                        double Lx,
                        double Ly,
                        double Lz,
                        double xy,
                        double xz,
                        double yz)
    {
    const unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= N)
        return;
    const double x = pred_d[3 * i + 0];
    const double y = pred_d[3 * i + 1];
    const double z = pred_d[3 * i + 2];
    double fz = z / Lz;
    double fy = (y - yz * Lz * fz) / Ly;
    double fx = (x - xy * Ly * fy - xz * Lz * fz) / Lx;
    const double nx_d = floor(fx + 0.5);
    const double ny_d = floor(fy + 0.5);
    const double nz_d = floor(fz + 0.5);
    fx -= nx_d;
    fy -= ny_d;
    fz -= nz_d;
    hoomd::Scalar4 p = pos[i]; // preserve .w (type)
    p.x = (hoomd::Scalar)(Lx * fx + xy * Ly * fy + xz * Lz * fz);
    p.y = (hoomd::Scalar)(Ly * fy + yz * Lz * fz);
    p.z = (hoomd::Scalar)(Lz * fz);
    pos[i] = p;
    int3 im = image[i];
    im.x += (int)nx_d;
    im.y += (int)ny_d;
    im.z += (int)nz_d;
    image[i] = im;
    }

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
                                unsigned int block_size)
    {
    if (N == 0)
        return hipSuccess;
    const unsigned int gN = (N + block_size - 1) / block_size;

    hipLaunchKernelGGL(ck_extract, dim3(gN), dim3(block_size), 0, 0, d_pos, d_tag,
                       d_inv_gamma_by_tag, d_pos_d, d_ref_d, d_inv_gamma_row, d_row_of_tag, N);

    if (F > 0)
        {
        const unsigned int total = F * (m + 1);
        const unsigned int gT = (total + block_size - 1) / block_size;
        hipLaunchKernelGGL(ck_map_chains, dim3(gT), dim3(block_size), 0, 0, d_chains_tag,
                           d_row_of_tag, d_chains_row, total);
        hipMemsetAsync(d_ffix_d, 0, (size_t)3 * N * sizeof(double), 0);
        gpu_fixman_launch(d_pos_d, d_inv_gamma_row, d_chains_row, d_ffix_d, d_logdet, d_bad_sign,
                          F, m, Lx, Ly, Lz);
        }
    else
        {
        hipMemsetAsync(d_ffix_d, 0, (size_t)3 * N * sizeof(double), 0);
        }

    hipLaunchKernelGGL(ck_predict, dim3(gN), dim3(block_size), 0, 0, d_pos_d, d_force, d_ffix_d,
                       d_tag, d_inv_gamma_row, d_bd_pref_by_tag, d_prv, d_pred_d, N, seed, timestep,
                       dt, half_kT);

    if (F > 0)
        gpu_shake_launch(d_pred_d, d_ref_d, d_inv_gamma_row, d_chains_row, d_lambda, d_nonconv, F,
                         m, rest_length, Lx, Ly, Lz, tol, max_iter);

    hipLaunchKernelGGL(ck_wrap, dim3(gN), dim3(block_size), 0, 0, d_pred_d, d_pos, d_image, N, Lx,
                       Ly, Lz, xy, xz, yz);

    return hipPeekAtLastError();
    }

    } // namespace ffn_native

#endif
