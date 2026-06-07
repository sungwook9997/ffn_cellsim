#include "fixman_kernel.cuh"

#ifdef ENABLE_HIP

#include <stdexcept>
#include <string>

namespace ffn_native
    {

__device__ inline double ffn_fx_minimage(double d, double L)
    {
    return d - L * rint(d / L);
    }

// One thread per chain. m bonds, m+1 beads with row indices P[0..m].
// Mirrors fixman_logdet_and_force's uniform fast path (FIXMAN_SIGN = +1).
__global__ void fixman_kernel(const double* __restrict__ pos,
                              const double* __restrict__ inv_gamma,
                              const int* __restrict__ chains,
                              double* __restrict__ force,
                              double* __restrict__ logdet,
                              int* __restrict__ bad_sign,
                              const unsigned int F,
                              const unsigned int m,
                              const double Lx,
                              const double Ly,
                              const double Lz)
    {
    const unsigned int f = blockIdx.x * blockDim.x + threadIdx.x;
    if (f >= F)
        return;

    const int* P = chains + (size_t)f * (m + 1);

    double bx[FFN_FIXMAN_MMAX], by[FFN_FIXMAN_MMAX], bz[FFN_FIXMAN_MMAX];
    double M[FFN_FIXMAN_MMAX + 1];
    double dg[FFN_FIXMAN_MMAX];       // G diagonal
    double e[FFN_FIXMAN_MMAX];        // G[a,a+1], a=0..m-2
    double th[FFN_FIXMAN_MMAX + 1];   // theta_0..theta_m
    double ph[FFN_FIXMAN_MMAX + 2];   // phi_1..phi_{m+1}
    double gx[FFN_FIXMAN_MMAX], gy[FFN_FIXMAN_MMAX], gz[FFN_FIXMAN_MMAX]; // grad per bond

    for (unsigned int c = 0; c <= m; ++c)
        M[c] = inv_gamma[P[c]];
    for (unsigned int a = 0; a < m; ++a)
        {
        const int p0 = P[a], p1 = P[a + 1];
        bx[a] = ffn_fx_minimage(pos[3 * p1 + 0] - pos[3 * p0 + 0], Lx);
        by[a] = ffn_fx_minimage(pos[3 * p1 + 1] - pos[3 * p0 + 1], Ly);
        bz[a] = ffn_fx_minimage(pos[3 * p1 + 2] - pos[3 * p0 + 2], Lz);
        const double b2 = bx[a] * bx[a] + by[a] * by[a] + bz[a] * bz[a];
        dg[a] = 4.0 * b2 * (M[a] + M[a + 1]);
        }
    for (unsigned int a = 0; a + 1 < m; ++a)
        {
        const double bdot = bx[a] * bx[a + 1] + by[a] * by[a + 1] + bz[a] * bz[a + 1];
        e[a] = -4.0 * M[a + 1] * bdot;
        }

    // tridiagonal continuants: theta (leading minors), phi (trailing minors)
    th[0] = 1.0;
    th[1] = dg[0];
    for (unsigned int k = 2; k <= m; ++k)
        th[k] = dg[k - 1] * th[k - 1] - e[k - 2] * e[k - 2] * th[k - 2];
    ph[m + 1] = 1.0;
    ph[m] = dg[m - 1];
    for (int k = (int)m - 1; k >= 1; --k)
        ph[k] = dg[k - 1] * ph[k + 1] - e[k - 1] * e[k - 1] * ph[k + 2];

    const double det = th[m];
    bad_sign[f] = (det > 0.0) ? 0 : 1;
    logdet[f] = log(fabs(det));
    const double inv_det = 1.0 / det;

    // grad_a = 0.5kT factor applied on host; here accumulate d ln det / d b_a
    for (unsigned int a = 0; a < m; ++a)
        {
        const double gii = th[a] * ph[a + 2] * inv_det;           // Ginv[a,a]
        double dx = gii * 8.0 * (M[a] + M[a + 1]) * bx[a];
        double dy = gii * 8.0 * (M[a] + M[a + 1]) * by[a];
        double dz = gii * 8.0 * (M[a] + M[a + 1]) * bz[a];
        if (a >= 1)
            {
            const double glo = -e[a - 1] * th[a - 1] * ph[a + 2] * inv_det; // Ginv[a-1,a]
            const double c = 2.0 * glo * (-4.0 * M[a]);
            dx += c * bx[a - 1];
            dy += c * by[a - 1];
            dz += c * bz[a - 1];
            }
        if (a + 1 < m)
            {
            const double gup = -e[a] * th[a] * ph[a + 3] * inv_det;          // Ginv[a,a+1]
            const double c = 2.0 * gup * (-4.0 * M[a + 1]);
            dx += c * bx[a + 1];
            dy += c * by[a + 1];
            dz += c * bz[a + 1];
            }
        gx[a] = dx;  // = d ln det / d b_a  (× 0.5kT on host → grad)
        gy[a] = dy;
        gz[a] = dz;
        }

    // F_{p_a} += grad_a ; F_{p_{a+1}} -= grad_a  (grad scaled by 0.5kT on host)
    for (unsigned int b = 0; b <= m; ++b)
        {
        double fx = 0.0, fy = 0.0, fz = 0.0;
        if (b < m) { fx += gx[b]; fy += gy[b]; fz += gz[b]; }
        if (b > 0) { fx -= gx[b - 1]; fy -= gy[b - 1]; fz -= gz[b - 1]; }
        const int pb = P[b];
        force[3 * pb + 0] = fx; // disjoint chains → plain write (caller zeroed)
        force[3 * pb + 1] = fy;
        force[3 * pb + 2] = fz;
        }
    }

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
                             double Lz)
    {
    if (F == 0)
        return hipSuccess;
    const unsigned int block = 128;
    const unsigned int grid = (F + block - 1) / block;
    hipLaunchKernelGGL(fixman_kernel, dim3(grid), dim3(block), 0, 0, d_pos, d_inv_gamma, d_chains,
                       d_force, d_logdet, d_bad_sign, F, m, Lx, Ly, Lz);
    return hipPeekAtLastError();
    }

static void ffn_fx_check(hipError_t e, const char* what)
    {
    if (e != hipSuccess)
        throw std::runtime_error(std::string("fixman_force_host ") + what + ": "
                                 + hipGetErrorString(e));
    }

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
                             double Lz)
    {
    if (m >= FFN_FIXMAN_MMAX)
        throw std::runtime_error("fixman_force_host: chain too long (m >= FFN_FIXMAN_MMAX)");
    if (F == 0)
        return hipSuccess;

    double *d_pos = nullptr, *d_invg = nullptr, *d_force = nullptr, *d_logdet = nullptr;
    int *d_chains = nullptr, *d_bad = nullptr;
    const size_t pos_bytes = (size_t)3 * N * sizeof(double);
    const size_t invg_bytes = (size_t)N * sizeof(double);
    const size_t chain_bytes = (size_t)F * (m + 1) * sizeof(int);
    const size_t ld_bytes = (size_t)F * sizeof(double);
    const size_t bad_bytes = (size_t)F * sizeof(int);

    ffn_fx_check(hipMalloc(&d_pos, pos_bytes), "malloc pos");
    ffn_fx_check(hipMalloc(&d_invg, invg_bytes), "malloc invg");
    ffn_fx_check(hipMalloc(&d_chains, chain_bytes), "malloc chains");
    ffn_fx_check(hipMalloc(&d_force, pos_bytes), "malloc force");
    ffn_fx_check(hipMalloc(&d_logdet, ld_bytes), "malloc logdet");
    ffn_fx_check(hipMalloc(&d_bad, bad_bytes), "malloc bad");

    ffn_fx_check(hipMemcpy(d_pos, pos, pos_bytes, hipMemcpyHostToDevice), "H2D pos");
    ffn_fx_check(hipMemcpy(d_invg, inv_gamma, invg_bytes, hipMemcpyHostToDevice), "H2D invg");
    ffn_fx_check(hipMemcpy(d_chains, chains, chain_bytes, hipMemcpyHostToDevice), "H2D chains");
    ffn_fx_check(hipMemset(d_force, 0, pos_bytes), "memset force");

    const unsigned int block = 128;
    const unsigned int grid = (F + block - 1) / block;
    hipLaunchKernelGGL(fixman_kernel, dim3(grid), dim3(block), 0, 0, d_pos, d_invg, d_chains,
                       d_force, d_logdet, d_bad, F, m, Lx, Ly, Lz);
    const hipError_t kerr = hipDeviceSynchronize();

    ffn_fx_check(hipMemcpy(force, d_force, pos_bytes, hipMemcpyDeviceToHost), "D2H force");
    ffn_fx_check(hipMemcpy(logdet, d_logdet, ld_bytes, hipMemcpyDeviceToHost), "D2H logdet");
    ffn_fx_check(hipMemcpy(bad_sign, d_bad, bad_bytes, hipMemcpyDeviceToHost), "D2H bad");

    hipFree(d_pos);
    hipFree(d_invg);
    hipFree(d_chains);
    hipFree(d_force);
    hipFree(d_logdet);
    hipFree(d_bad);
    return kerr;
    }

    } // namespace ffn_native

#endif
