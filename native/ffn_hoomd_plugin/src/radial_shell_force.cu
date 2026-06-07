#include "radial_shell_force.cuh"

#ifdef ENABLE_HIP

namespace ffn_native
    {

__global__ void rsf_zero(hoomd::Scalar* acc)
    {
    if (threadIdx.x < 5)
        acc[threadIdx.x] = 0.0;
    }

// Pass 1: shell centroid sum + count.
__global__ void rsf_reduce_centroid(const hoomd::Scalar4* __restrict__ pos,
                                    const unsigned int* __restrict__ tag,
                                    hoomd::Scalar* __restrict__ acc,
                                    unsigned int N,
                                    unsigned int t0,
                                    unsigned int t1)
    {
    const unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= N)
        return;
    const unsigned int t = tag[i];
    if (t < t0 || t >= t1)
        return;
    const hoomd::Scalar4 p = pos[i];
    atomicAdd(&acc[0], (hoomd::Scalar)p.x);
    atomicAdd(&acc[1], (hoomd::Scalar)p.y);
    atomicAdd(&acc[2], (hoomd::Scalar)p.z);
    atomicAdd(&acc[3], (hoomd::Scalar)1.0);
    }

// Pass 2: shell mean-radius sum (needs the centroid from pass 1).
__global__ void rsf_reduce_radius(const hoomd::Scalar4* __restrict__ pos,
                                  const unsigned int* __restrict__ tag,
                                  hoomd::Scalar* __restrict__ acc,
                                  unsigned int N,
                                  unsigned int t0,
                                  unsigned int t1)
    {
    const unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= N)
        return;
    const unsigned int t = tag[i];
    if (t < t0 || t >= t1)
        return;
    const double cnt = acc[3] > 0.0 ? acc[3] : 1.0;
    const double cx = acc[0] / cnt, cy = acc[1] / cnt, cz = acc[2] / cnt;
    const hoomd::Scalar4 p = pos[i];
    const double dx = (double)p.x - cx, dy = (double)p.y - cy, dz = (double)p.z - cz;
    atomicAdd(&acc[4], (hoomd::Scalar)sqrt(dx * dx + dy * dy + dz * dz));
    }

// Pass 3: per-bead radial force + per-bead potential (stored in force.w).
__global__ void rsf_apply(const hoomd::Scalar4* __restrict__ pos,
                          const unsigned int* __restrict__ tag,
                          hoomd::Scalar4* __restrict__ force,
                          const hoomd::Scalar* __restrict__ acc,
                          unsigned int N,
                          unsigned int t0,
                          unsigned int t1,
                          int law,
                          double R0,
                          double pa,
                          double pb,
                          double pc,
                          double pd)
    {
    const unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= N)
        return;
    const unsigned int t = tag[i];
    if (t < t0 || t >= t1)
        {
        hoomd::Scalar4 fz;
        fz.x = 0;
        fz.y = 0;
        fz.z = 0;
        fz.w = 0;
        force[i] = fz;
        return;
        }
    const double cnt = acc[3] > 0.0 ? acc[3] : 1.0;
    const double cx = acc[0] / cnt, cy = acc[1] / cnt, cz = acc[2] / cnt;
    const double R_mean = acc[4] / cnt;
    const hoomd::Scalar4 p = pos[i];
    const double dx = (double)p.x - cx, dy = (double)p.y - cy, dz = (double)p.z - cz;
    const double r = sqrt(dx * dx + dy * dy + dz * dz);
    const double rs = r > 0.0 ? r : 1.0;
    double nx = dx / rs, ny = dy / rs, nz = dz / rs;
    if (r <= 0.0)
        {
        nx = 0.0;
        ny = 0.0;
        nz = 0.0;
        }
    const double PI = 3.14159265358979323846;
    double Fmag = 0.0, U = 0.0;
    if (law == 0)
        { // nucleus: bilinear about R0 (pa=k_chrom, pb=k_lamin, pc=d_knee, pd=F_knee)
        const double d = r - R0, ad = fabs(d), sgn = (double)((d > 0) - (d < 0));
        const double e = ad - pc, uknee = 0.5 * pa * pc * pc;
        if (ad <= pc)
            {
            Fmag = -pa * d;
            U = 0.5 * pa * d * d;
            }
        else
            {
            Fmag = -sgn * (pd + (pa + pb) * e);
            U = uknee + pd * e + 0.5 * (pa + pb) * e * e;
            }
        }
    else if (law == 1)
        { // membrane: S=4piR^2 Laplace (pa=gamma_mem, pb=K_A, pc=A0); inward
        const double S = 4.0 * PI * R_mean * R_mean;
        const double gamma_tot = pa + pb * (S - pc) / pc;
        const double Rm = R_mean > 0.0 ? R_mean : 1e-300;
        const double dP = 2.0 * gamma_tot / Rm;
        const double A_i = S / cnt;
        Fmag = -(dP * A_i);
        const double Uc = fmax(pa * (S - pc), 0.0);
        U = (Uc + 0.5 * pb * (S - pc) * (S - pc) / pc) / cnt;
        }
    else
        { // turgor: V=4/3 pi R^3 (pa=dP0, pb=K_vol, pc=V0); outward
        const double V = (4.0 / 3.0) * PI * R_mean * R_mean * R_mean;
        const double S = 4.0 * PI * R_mean * R_mean;
        const double dP = pa - pb * (V - pc) / pc;
        const double A_i = S / cnt;
        Fmag = (dP * A_i);
        U = (0.5 * (pb / pc) * (V - pc) * (V - pc)) / cnt;
        }
    hoomd::Scalar4 fo;
    fo.x = (hoomd::Scalar)(Fmag * nx);
    fo.y = (hoomd::Scalar)(Fmag * ny);
    fo.z = (hoomd::Scalar)(Fmag * nz);
    fo.w = (hoomd::Scalar)U;
    force[i] = fo;
    }

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
                                  unsigned int block_size)
    {
    if (N == 0)
        return hipSuccess;
    const unsigned int grid = (N + block_size - 1) / block_size;
    hipLaunchKernelGGL(rsf_zero, dim3(1), dim3(32), 0, 0, d_acc);
    hipLaunchKernelGGL(rsf_reduce_centroid, dim3(grid), dim3(block_size), 0, 0, d_pos, d_tag,
                       d_acc, N, tag_start, tag_end);
    hipLaunchKernelGGL(rsf_reduce_radius, dim3(grid), dim3(block_size), 0, 0, d_pos, d_tag, d_acc,
                       N, tag_start, tag_end);
    hipLaunchKernelGGL(rsf_apply, dim3(grid), dim3(block_size), 0, 0, d_pos, d_tag, d_force, d_acc,
                       N, tag_start, tag_end, law, R0, pa, pb, pc, pd);
    return hipPeekAtLastError();
    }

    } // namespace ffn_native

#endif
