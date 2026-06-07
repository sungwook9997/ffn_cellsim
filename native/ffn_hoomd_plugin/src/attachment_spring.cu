#include "attachment_spring.cuh"

#ifdef ENABLE_HIP

namespace ffn_native
    {

__global__ void as_row_of_tag(const unsigned int* __restrict__ tag, int* __restrict__ row_of_tag,
                              unsigned int N)
    {
    const unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= N)
        return;
    row_of_tag[tag[i]] = (int)i;
    }

__global__ void as_zero_force(hoomd::Scalar4* __restrict__ force, unsigned int N)
    {
    const unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= N)
        return;
    hoomd::Scalar4 z;
    z.x = 0;
    z.y = 0;
    z.z = 0;
    z.w = 0;
    force[i] = z;
    }

__device__ inline double as_minimage(double d, double L)
    {
    return d - L * rint(d / L);
    }

// One thread per attachment slot. Inactive slots (k<=0 or head_tag<0) skip.
// Spring on head = +k(r-r0) toward actin; on actin = -that. Per-bead force is
// accumulated with atomicAdd (a bead carries many springs).
__global__ void as_spring(const hoomd::Scalar4* __restrict__ pos,
                          const int* __restrict__ row_of_tag,
                          const int* __restrict__ head_tag,
                          const int* __restrict__ actin_tag,
                          const double* __restrict__ kk,
                          const double* __restrict__ r0,
                          hoomd::Scalar4* __restrict__ force,
                          unsigned int K,
                          double Lx,
                          double Ly,
                          double Lz)
    {
    const unsigned int s = blockIdx.x * blockDim.x + threadIdx.x;
    if (s >= K)
        return;
    const int ht = head_tag[s];
    const int at = actin_tag[s];
    const double k = kk[s];
    if (ht < 0 || at < 0 || k <= 0.0)
        return;
    const int hr = row_of_tag[ht];
    const int ar = row_of_tag[at];
    const hoomd::Scalar4 ph = pos[hr];
    const hoomd::Scalar4 pa = pos[ar];
    double dx = as_minimage((double)pa.x - (double)ph.x, Lx);
    double dy = as_minimage((double)pa.y - (double)ph.y, Ly);
    double dz = as_minimage((double)pa.z - (double)ph.z, Lz);
    const double r = sqrt(dx * dx + dy * dy + dz * dz);
    if (r <= 0.0)
        return;
    const double fmag = k * (r - r0[s]) / r; // force/length along d (toward actin on head)
    const double fx = fmag * dx, fy = fmag * dy, fz = fmag * dz;
    // U = 0.5 k (r-r0)^2, split half to each bead's .w (potential).
    const double half_u = 0.25 * k * (r - r0[s]) * (r - r0[s]);
    atomicAdd(&force[hr].x, (hoomd::Scalar)fx);
    atomicAdd(&force[hr].y, (hoomd::Scalar)fy);
    atomicAdd(&force[hr].z, (hoomd::Scalar)fz);
    atomicAdd(&force[hr].w, (hoomd::Scalar)half_u);
    atomicAdd(&force[ar].x, (hoomd::Scalar)(-fx));
    atomicAdd(&force[ar].y, (hoomd::Scalar)(-fy));
    atomicAdd(&force[ar].z, (hoomd::Scalar)(-fz));
    atomicAdd(&force[ar].w, (hoomd::Scalar)half_u);
    }

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
                                 unsigned int block_size)
    {
    if (N == 0)
        return hipSuccess;
    const unsigned int gN = (N + block_size - 1) / block_size;
    hipLaunchKernelGGL(as_row_of_tag, dim3(gN), dim3(block_size), 0, 0, d_tag, d_row_of_tag, N);
    hipLaunchKernelGGL(as_zero_force, dim3(gN), dim3(block_size), 0, 0, d_force, N);
    if (K > 0)
        {
        const unsigned int gK = (K + block_size - 1) / block_size;
        hipLaunchKernelGGL(as_spring, dim3(gK), dim3(block_size), 0, 0, d_pos, d_row_of_tag,
                           d_head_tag, d_actin_tag, d_k, d_r0, d_force, K, Lx, Ly, Lz);
        }
    return hipPeekAtLastError();
    }

    } // namespace ffn_native

#endif
