#//////////////////////////////////////////////////////////////////////////////////////
#// This file is distributed under the University of Illinois/NCSA Open Source License.
#// See LICENSE file in top directory for details.
#//
#// Copyright (c) 2021 QMCPACK developers.
#//
#// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
#//
#// File created by: Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
#//////////////////////////////////////////////////////////////////////////////////////

set(SRC_DIR lammps)
set(UTEST_EXE test_${SRC_DIR})
set(UTEST_NAME deterministic-unit_test_${SRC_DIR})

add_executable(
  ${UTEST_EXE}
  test_lammps.cpp
  )

foreach(fname ${FILES_TO_COPY})
execute_process(COMMAND ${CMAKE_COMMAND} -E copy "${CMAKE_CURRENT_SOURCE_DIR}/${fname}" ${UTEST_DIR}/lammps)
endforeach()


target_link_libraries(${UTEST_EXE} catch_main qmcparticle qmcwfs Lammps::lammps)
target_link_options(${UTEST_EXE}  PUBLIC "${LAMMPS_LDFLAGS}")
target_compile_options(${UTEST_EXE} PUBLIC "${LAMMPS_CFLAGS}")

if(USE_OBJECT_TARGET)
  target_link_libraries(${UTEST_EXE} qmcutil qmcparticle_omptarget qmcwfs_omptarget  )
endif()

add_unit_test(${UTEST_NAME} 1 1 $<TARGET_FILE:${UTEST_EXE}>)
