#//////////////////////////////////////////////////////////////////////////////////////
#// This file is distributed under the University of Illinois/NCSA Open Source License.
#// See LICENSE file in top directory for details.
#//
#// Copyright (c) 2020 QMCPACK developers.
#//
#// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
#//
#// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
#//////////////////////////////////////////////////////////////////////////////////////

set(CUDA_SRCS cuBLAS_missing_functions.cu CUDAfill.cpp CUDAallocator.cpp CUDAruntime.cpp)

cuda_add_library(platform_cuda ${CUDA_SRCS})
cuda_add_cublas_to_target(platform_cuda)
target_link_libraries(platform_cuda PRIVATE ${CUDA_LIBRARIES} ${CUDA_cusolver_LIBRARY})
