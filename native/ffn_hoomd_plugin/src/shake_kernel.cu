#include "shake_kernel.cuh"

#ifdef ENABLE_HIP

#include <stdexcept>
#include <string>

namespace ffn_native
    {

__device__ inline double ffn_minimage(double d, double L)
    {
    return d - L * rint(d / L); // orthorhombic min-image (rint = round-half-even = numpy round)
    }

// One thread per chain. m bonds, m+1 beads with row indices P[0..m].
// Mirrors shake_project_chains' uniform fast path op-for-op (bit-identical).
__global__ void shake_kernel(double* __restrict__ pos,
                             const double* __restrict__ ref,
                             const double* __restrict__ inv_mass,
                             const int* __restrict__ chains,
                             double* __restrict__ lambda,
                             int* __restrict__ nonconv,
                             const unsigned int F,
                             const unsigned int m,
                             const double L2,
                             const double Lx,
                             const double Ly,
                             const double Lz,
                             const double tol,
                             const unsigned int max_iter)
    {
    const unsigned int f = blockIdx.x * blockDim.x + threadIdx.x;
    if (f >= F)
        return;

    const int* P = chains + (size_t)f * (m + 1);

    double d0x[FFN_SHAKE_MMAX], d0y[FFN_SHAKE_MMAX], d0z[FFN_SHAKE_MMAX];
    double Mass[FFN_SHAKE_MMAX + 1];
    double lamt[FFN_SHAKE_MMAX];
    for (unsigned int b = 0; b <= m; ++b)
        Mass[b] = inv_mass[P[b]];
    for (unsigned int a = 0; a < m; ++a)
        {
        const int p0 = P[a], p1 = P[a + 1];
        d0x[a] = ffn_minimage(ref[3 * p0 + 0] - ref[3 * p1 + 0], Lx);
        d0y[a] = ffn_minimage(ref[3 * p0 + 1] - ref[3 * p1 + 1], Ly);
        d0z[a] = ffn_minimage(ref[3 * p0 + 2] - ref[3 * p1 + 2], Lz);
        lamt[a] = 0.0;
        }

    double sx[FFN_SHAKE_MMAX], sy[FFN_SHAKE_MMAX], sz[FFN_SHAKE_MMAX], g[FFN_SHAKE_MMAX];
    double diag[FFN_SHAKE_MMAX], sub[FFN_SHAKE_MMAX], sup[FFN_SHAKE_MMAX];
    double cc[FFN_SHAKE_MMAX], dp[FFN_SHAKE_MMAX], lam[FFN_SHAKE_MMAX];

    bool converged = false;
    for (unsigned int iter = 0; iter < max_iter; ++iter)
        {
        double maxg = 0.0;
        for (unsigned int a = 0; a < m; ++a)
            {
            const int p0 = P[a], p1 = P[a + 1];
            sx[a] = ffn_minimage(pos[3 * p0 + 0] - pos[3 * p1 + 0], Lx);
            sy[a] = ffn_minimage(pos[3 * p0 + 1] - pos[3 * p1 + 1], Ly);
            sz[a] = ffn_minimage(pos[3 * p0 + 2] - pos[3 * p1 + 2], Lz);
            g[a] = sx[a] * sx[a] + sy[a] * sy[a] + sz[a] * sz[a] - L2;
            const double ag = fabs(g[a]);
            if (ag > maxg)
                maxg = ag;
            }
        if (maxg / L2 <= tol)
            {
            converged = true;
            break;
            }

        for (unsigned int a = 0; a < m; ++a)
            {
            const double sd = sx[a] * d0x[a] + sy[a] * d0y[a] + sz[a] * d0z[a];
            diag[a] = -2.0 * (Mass[a] + Mass[a + 1]) * sd;
            sub[a] = 0.0;
            sup[a] = 0.0;
            }
        for (unsigned int a = 1; a < m; ++a)
            sub[a] = 2.0 * Mass[a]
                     * (sx[a] * d0x[a - 1] + sy[a] * d0y[a - 1] + sz[a] * d0z[a - 1]);
        for (unsigned int a = 0; a + 1 < m; ++a)
            sup[a] = 2.0 * Mass[a + 1]
                     * (sx[a] * d0x[a + 1] + sy[a] * d0y[a + 1] + sz[a] * d0z[a + 1]);

        // Thomas solve: diag/sub/sup, rhs = -g  ->  lam
        cc[0] = sup[0] / diag[0];
        dp[0] = (-g[0]) / diag[0];
        for (unsigned int k = 1; k < m; ++k)
            {
            const double den = diag[k] - sub[k] * cc[k - 1];
            cc[k] = sup[k] / den;
            dp[k] = ((-g[k]) - sub[k] * dp[k - 1]) / den;
            }
        lam[m - 1] = dp[m - 1];
        for (int k = (int)m - 2; k >= 0; --k)
            lam[k] = dp[k] - cc[k] * lam[k + 1];

        for (unsigned int a = 0; a < m; ++a)
            lamt[a] += lam[a];

        // apply Δr[p_b] = Mass[b]*(λ_{b-1} d0_{b-1} − λ_b d0_b)
        for (unsigned int b = 0; b <= m; ++b)
            {
            double dxx = 0.0, dyy = 0.0, dzz = 0.0;
            if (b > 0)
                {
                dxx += Mass[b] * lam[b - 1] * d0x[b - 1];
                dyy += Mass[b] * lam[b - 1] * d0y[b - 1];
                dzz += Mass[b] * lam[b - 1] * d0z[b - 1];
                }
            if (b < m)
                {
                dxx -= Mass[b] * lam[b] * d0x[b];
                dyy -= Mass[b] * lam[b] * d0y[b];
                dzz -= Mass[b] * lam[b] * d0z[b];
                }
            const int pb = P[b];
            pos[3 * pb + 0] += dxx;
            pos[3 * pb + 1] += dyy;
            pos[3 * pb + 2] += dzz;
            }
        }

    for (unsigned int a = 0; a < m; ++a)
        lambda[(size_t)f * m + a] = lamt[a];
    nonconv[f] = converged ? 0 : 1;
    }

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
                            unsigned int max_iter)
    {
    if (F == 0)
        return hipSuccess;
    const unsigned int block = 128;
    const unsigned int grid = (F + block - 1) / block;
    hipLaunchKernelGGL(shake_kernel, dim3(grid), dim3(block), 0, 0, d_pos, d_ref, d_inv_mass,
                       d_chains, d_lambda, d_nonconv, F, m, rest_length * rest_length, Lx, Ly, Lz,
                       tol, max_iter);
    return hipPeekAtLastError();
    }

static void ffn_check(hipError_t e, const char* what)
    {
    if (e != hipSuccess)
        throw std::runtime_error(std::string("shake_project_host ") + what + ": "
                                 + hipGetErrorString(e));
    }

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
                              int* out_nonconverged)
    {
    if (m >= FFN_SHAKE_MMAX)
        throw std::runtime_error("shake_project_host: chain too long (m >= FFN_SHAKE_MMAX)");
    if (F == 0)
        return hipSuccess;

    double *d_pos = nullptr, *d_ref = nullptr, *d_invm = nullptr, *d_lam = nullptr;
    int *d_chains = nullptr, *d_nonconv = nullptr;
    const size_t pos_bytes = (size_t)3 * N * sizeof(double);
    const size_t invm_bytes = (size_t)N * sizeof(double);
    const size_t chain_bytes = (size_t)F * (m + 1) * sizeof(int);
    const size_t lam_bytes = (size_t)F * m * sizeof(double);
    const size_t nc_bytes = (size_t)F * sizeof(int);

    ffn_check(hipMalloc(&d_pos, pos_bytes), "malloc pos");
    ffn_check(hipMalloc(&d_ref, pos_bytes), "malloc ref");
    ffn_check(hipMalloc(&d_invm, invm_bytes), "malloc invm");
    ffn_check(hipMalloc(&d_chains, chain_bytes), "malloc chains");
    ffn_check(hipMalloc(&d_lam, lam_bytes), "malloc lam");
    ffn_check(hipMalloc(&d_nonconv, nc_bytes), "malloc nonconv");

    ffn_check(hipMemcpy(d_pos, pos, pos_bytes, hipMemcpyHostToDevice), "H2D pos");
    ffn_check(hipMemcpy(d_ref, ref, pos_bytes, hipMemcpyHostToDevice), "H2D ref");
    ffn_check(hipMemcpy(d_invm, inv_mass, invm_bytes, hipMemcpyHostToDevice), "H2D invm");
    ffn_check(hipMemcpy(d_chains, chains, chain_bytes, hipMemcpyHostToDevice), "H2D chains");

    gpu_shake_launch(d_pos, d_ref, d_invm, d_chains, d_lam, d_nonconv, F, m, rest_length, Lx, Ly,
                     Lz, tol, max_iter);
    const hipError_t kerr = hipDeviceSynchronize();

    ffn_check(hipMemcpy(pos, d_pos, pos_bytes, hipMemcpyDeviceToHost), "D2H pos");
    ffn_check(hipMemcpy(lambda, d_lam, lam_bytes, hipMemcpyDeviceToHost), "D2H lam");
    ffn_check(hipMemcpy(out_nonconverged, d_nonconv, nc_bytes, hipMemcpyDeviceToHost), "D2H nc");

    hipFree(d_pos);
    hipFree(d_ref);
    hipFree(d_invm);
    hipFree(d_chains);
    hipFree(d_lam);
    hipFree(d_nonconv);
    return kerr;
    }

    } // namespace ffn_native

#endif
