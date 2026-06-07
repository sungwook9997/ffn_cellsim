#include "position_kick_kernel.cuh"

#ifdef ENABLE_HIP

namespace ffn_native
    {

__global__ void position_kick_kernel(hoomd::Scalar4* d_pos,
                                     unsigned int n_particles,
                                     hoomd::Scalar dx,
                                     hoomd::Scalar dy,
                                     hoomd::Scalar dz)
    {
    const unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n_particles)
        {
        return;
        }

    hoomd::Scalar4 pos = d_pos[idx];
    pos.x += dx;
    pos.y += dy;
    pos.z += dz;
    d_pos[idx] = pos;
    }

hipError_t gpu_position_kick(hoomd::Scalar4* d_pos,
                             unsigned int n_particles,
                             hoomd::Scalar dx,
                             hoomd::Scalar dy,
                             hoomd::Scalar dz,
                             unsigned int block_size)
    {
    const unsigned int grid = (n_particles + block_size - 1) / block_size;
    hipLaunchKernelGGL(position_kick_kernel,
                       dim3(grid),
                       dim3(block_size),
                       0,
                       0,
                       d_pos,
                       n_particles,
                       dx,
                       dy,
                       dz);
    return hipPeekAtLastError();
    }

    } // namespace ffn_native

#endif
