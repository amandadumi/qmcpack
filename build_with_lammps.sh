#!/bin/bash
#export GCC_PATH="/usr/local/Cellar/gcc/13.2.0"
export LDFLAGS="-L/usr/local/opt/libomp/lib"
export CPPFLAGS="-I/usr/local/opt/libomp/include"
#export CC="$GCC_PATH/bin/gcc-13"
#export CXX="$GCC_PATH/bin/g++-13"
export CXX="/usr/local/Cellar/llvm/17.0.6_1/bin/clang++"
export CC="/usr/local/Cellar/llvm/17.0.6_1/bin/clang"
export OMPI_CC=$CC
export OMPI_CXX=$CXX

#export CLANG_CXX_LANGUAGE_STANDARD="C++17"
#export CMAKE_CXX_FLAGS="-Wl,-ld_classic,-g"
export MPICC="/usr/local/Cellar/open-mpi/5.0.2_1/bin/mpicc"
export MPICXX="/usr/local/Cellar/open-mpi/5.0.2_1/bin/mpic++"
#export MPIEXEC_EXECUTABLE="/usr/local/Cellar/open-mpi/4.1.5/bin/mpiexec"
export FFTW_HOME="/usr/local/opt/fftw"
export LAMMPS_ROOT="/Users/aedumi/software/lammps/build"
export PKG_CONFIG_PATH="/usr/local/opt/readline/lib/pkgconfig"
export DYLD_LIBRARY_PATH="/Users/aedumi/software/lammps_serial/build"
export CMAKE_INCLUDE_PATH="/usr/local/Cellar/libomp/18.1.1/include"
echo "configured the environment"
#-DCMAKE_CXX_FLAGS="-fopenmp=libomp" \
mpic++ -v
cd build
cmake  -DCMAKE_C_COMPILER=$MPICC \
      -DCMAKE_CXX_COMPILER=$MPICXX \
      -DFFTW_INCLUDE_DIRS=$FFTW_HOME/include \
      -DFFTW_LIBRARY_DIRS=$FFTW_HOME/lib \
      -DCMAKE_MACOSX_RPATH=ON \
      -DLAMMPS_ROOT=$LAMMPS_ROOT \
      -DENABLE_OFFLOAD=OFF \
      -DQMC_COMPLEX=0 ..
make -j 3  # Adjust for available core count


