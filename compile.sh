. ~/spack/share/spack/setup-env.sh
#spack load gcc@13.2.0
spack load llvm/2ttrz2x
spack load cmake
spack load openblas
spack load openmpi
spack load hdf5+mpi/s352ljk
spack load boost
spack load fftw
#export CXX=$(spack location -i gcc)/bin/g++
export CXX=$(spack location -i llvm/2ttrz2x)/bin/clang++
export CC=$(spack location -i llvm/2ttrz2x)/bin/clang  
export MKL_ROOT=$(spack location -i openblas)/lib
mkdir -p build_clang
cd build_clang
cmake -DCMAKE_C_COMPILER=$CC \
      -DCMAKE_CXX_COMPILER=$CXX \
      -DCMAKE_CXX_FLAGS="-g -O0" \
      -DQMC_COMPLEX=0 ..
make -j6 # Adjust for available core count

