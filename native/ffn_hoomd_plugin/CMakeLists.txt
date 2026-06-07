cmake_minimum_required(VERSION 3.23)
project(ffn_hoomd_plugin LANGUAGES CXX)

find_package(HOOMD REQUIRED CONFIG)

hoomd_add_module(_ffn_native src/module.cc src/position_kick_kernel.cu src/baoab_kernel.cu src/shake_kernel.cu src/fixman_kernel.cu src/constrained_kernel.cu)
target_link_libraries(_ffn_native PRIVATE HOOMD::_hoomd HOOMD::_md)
target_compile_features(_ffn_native PRIVATE cxx_std_17)
