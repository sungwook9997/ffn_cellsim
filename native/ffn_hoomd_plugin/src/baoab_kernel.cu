#include "baoab_kernel.cuh"

#ifdef ENABLE_HIP

namespace ffn_native
    {

// --- counter-based RNG (identical to baoab_device.py::baoab_step) ---------- //
__device__ inline unsigned long long ffn_splitmix64(unsigned long long x)
    {
    x += 0x9E3779B97F4A7C15ULL;
    x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ULL;
    x = (x ^ (x >> 27)) * 0x94D049BB133111EBULL;
    return x ^ (x >> 31);
    }

__device__ inline double ffn_uniform01(unsigned long long x)
    {
    const unsigned long long r = ffn_splitmix64(x);
    return ((double)(r >> 11) + 0.5) * 1.1102230246251565e-16;
    }

__device__ inline void ffn_normal_pair(unsigned long long base, double* z0, double* z1)
    {
    double u1 = ffn_uniform01(base);
    const double u2 = ffn_uniform01(base ^ 0xD1B54A32D192ED03ULL);
    if (u1 < 1.0e-300)
        {
        u1 = 1.0e-300;
        }
    const double r = sqrt(-2.0 * log(u1));
    const double theta = 6.283185307179586476925286766559 * u2;
    *z0 = r * cos(theta);
    *z1 = r * sin(theta);
    }

__global__ void baoab_step_kernel(hoomd::Scalar4* __restrict__ d_pos,
                                  const hoomd::Scalar4* __restrict__ d_force,
                                  int3* __restrict__ d_image,
                                  const unsigned int* __restrict__ d_tag,
                                  const hoomd::Scalar* __restrict__ d_gamma_by_tag,
                                  const hoomd::Scalar* __restrict__ d_bd_prefactor_by_tag,
                                  hoomd::Scalar* __restrict__ d_prv,
                                  const unsigned int N,
                                  const unsigned long long seed,
                                  const unsigned long long timestep,
                                  const double dt,
                                  const double Lx,
                                  const double Ly,
                                  const double Lz,
                                  const double xy,
                                  const double xz,
                                  const double yz)
    {
    const unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= N)
        {
        return;
        }

    const unsigned int t = d_tag[i];
    const unsigned long long tp = 3ULL * (unsigned long long)t;
    const double inv_gamma = 1.0 / (double)d_gamma_by_tag[t];
    const double pref = (double)d_bd_prefactor_by_tag[t];

    const unsigned long long rng_base = seed ^ (timestep * 0xD2B74407B1CE6E93ULL)
                                        ^ (((unsigned long long)t + 1ULL) * 0xCA5A826395121157ULL);
    double W0, W1, W2, W_unused;
    ffn_normal_pair(rng_base, &W0, &W1);
    ffn_normal_pair(rng_base ^ 0x9E3779B97F4A7C15ULL, &W2, &W_unused);

    hoomd::Scalar4 p = d_pos[i];      // .w = type (preserved)
    const hoomd::Scalar4 f = d_force[i];

    double x = (double)p.x + (double)f.x * inv_gamma * dt + pref * (W0 + (double)d_prv[tp + 0]) * dt;
    double y = (double)p.y + (double)f.y * inv_gamma * dt + pref * (W1 + (double)d_prv[tp + 1]) * dt;
    double z = (double)p.z + (double)f.z * inv_gamma * dt + pref * (W2 + (double)d_prv[tp + 2]) * dt;

    // triclinic fractional-coordinate wrap (identical to the cupy kernel)
    double fz = z / Lz;
    double fy = (y - yz * Lz * fz) / Ly;
    double fx = (x - xy * Ly * fy - xz * Lz * fz) / Lx;

    const double nx_d = floor(fx + 0.5);
    const double ny_d = floor(fy + 0.5);
    const double nz_d = floor(fz + 0.5);

    fx -= nx_d;
    fy -= ny_d;
    fz -= nz_d;

    p.x = (hoomd::Scalar)(Lx * fx + xy * Ly * fy + xz * Lz * fz);
    p.y = (hoomd::Scalar)(Ly * fy + yz * Lz * fz);
    p.z = (hoomd::Scalar)(Lz * fz);
    d_pos[i] = p;

    int3 img = d_image[i];
    img.x += (int)nx_d;
    img.y += (int)ny_d;
    img.z += (int)nz_d;
    d_image[i] = img;

    d_prv[tp + 0] = (hoomd::Scalar)W0;
    d_prv[tp + 1] = (hoomd::Scalar)W1;
    d_prv[tp + 2] = (hoomd::Scalar)W2;
    }

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
                          unsigned int block_size)
    {
    if (n_particles == 0)
        {
        return hipSuccess;
        }
    const unsigned int grid = (n_particles + block_size - 1) / block_size;
    hipLaunchKernelGGL(baoab_step_kernel,
                       dim3(grid),
                       dim3(block_size),
                       0,
                       0,
                       d_pos,
                       d_force,
                       d_image,
                       d_tag,
                       d_gamma_by_tag,
                       d_bd_prefactor_by_tag,
                       d_prv,
                       n_particles,
                       seed,
                       timestep,
                       dt,
                       Lx,
                       Ly,
                       Lz,
                       xy,
                       xz,
                       yz);
    return hipPeekAtLastError();
    }

    } // namespace ffn_native

#endif
