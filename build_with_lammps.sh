#!/bin/bash
export GCC_PATH="/usr/local/Cellar/gcc/13.2.0"
export LDFLAGS="-L/usr/local/opt/libomp/lib"
export CPPFLAGS="-I/usr/local/opt/libomp/include"
export CC="$GCC_PATH/bin/gcc-13"
export CXX="$GCC_PATH/bin/g++-13"

#export CLANG_CXX_LANGUAGE_STANDARD="C++17"
export CMAKE_CXX_FLAGS="-Wl,-ld_classic,-g"
#export CXX="/usr/bin/clang++"
#export MPICC="/usr/local/Cellar/open-mpi/4.1.5/bin/mpicc"
#export MPICXX="/usr/local/Cellar/open-mpi/4.1.5/bin/mpic++"
#export MPIEXEC_EXECUTABLE="/usr/local/Cellar/open-mpi/4.1.5/bin/mpiexec"
export HDF5_ROOT="/usr/local/opt/hdf5"
export FFTW_HOME="/usr/local/opt/fftw"
export LAMMPS_ROOT="/Users/aedumi/software/lammps/build"
export PKG_CONFIG_PATH="/usr/local/opt/readline/lib/pkgconfig"
export DYLD_LIBRARY_PATH="$GCC_PATH/gcc/13:/Users/aedumi/software/lammps_serial/build"
echo "configured the environment"
#-DCMAKE_CXX_FLAGS="-fopenmp=libomp" \
cd build
cmake  -DCMAKE_C_COMPILER=$CC \
      -DCMAKE_CXX_COMPILER=$CXX \
      -DCMAKE_BUILD_TYPE=debug \
      -DCMAKE_CXX_FLAGS="-Wl,-ld_classic -g" \
      -DCMAKE_CXX_FLAGS_DEBUG="-Wl,-ld_classic -g" \
      -DFFTW_INCLUDE_DIRS=$FFTW_HOME/include \
      -DHDF5_ROOT=$HDF5_ROOT \
      -DFFTW_LIBRARY_DIRS=$FFTW_HOME/lib \
      -DCMAKE_MACOSX_RPATH=ON \
      -DLAMMPS_ROOT=$LAMMPS_ROOT \
      -DQMC_OMP=OFF \
      -DDEBUG_PSIBUFFER_ON=OFF \
      -DQMC_COMPLEX=0 ..
make -j 3  # Adjust for available core count


