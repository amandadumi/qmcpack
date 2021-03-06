#//////////////////////////////////////////////////////////////////////////////////////
#// This file is distributed under the University of Illinois/NCSA Open Source License.
#// See LICENSE file in top directory for details.
#//
#// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
#//
#// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
#//                    Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
#//
#// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
#//////////////////////////////////////////////////////////////////////////////////////

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${QMCPACK_UNIT_TEST_DIR})

set(SRC_DIR type_traits)
set(UTEST_EXE test_${SRC_DIR})
set(UTEST_NAME deterministic-unit_test_${SRC_DIR})
set(UTEST_DIR ${qmcpack_BINARY_DIR}/tests/type_traits)

execute_process(COMMAND ${CMAKE_COMMAND} -E make_directory "${UTEST_DIR}")

set(TEST_SRCS test_qmctypes.cpp test_template_types.cpp)
if(QMC_CUDA)
  set(TEST_SRCS test_cuda_types.cpp)
endif(QMC_CUDA)
add_executable(${UTEST_EXE} ${TEST_SRCS})
target_link_libraries(${UTEST_EXE} catch_main containers)

add_unit_test(${UTEST_NAME} 1 1 "${QMCPACK_UNIT_TEST_DIR}/${UTEST_EXE}")
