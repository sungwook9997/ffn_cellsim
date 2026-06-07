#pragma once

#include <hoomd/HOOMDMath.h>

#include "hip/hip_runtime.h"

namespace ffn_native
    {
hipError_t gpu_position_kick(hoomd::Scalar4* d_pos,
                             unsigned int n_particles,
                             hoomd::Scalar dx,
                             hoomd::Scalar dy,
                             hoomd::Scalar dz,
                             unsigned int block_size);
    }
