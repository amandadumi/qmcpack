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
mkdir -p build_clang_complex
cd build_clang_complex
cmake -DCMAKE_C_COMPILER=$CC \
      -DCMAKE_CXX_COMPILER=$CXX \
      -DQMC_COMPLEX=1  ..
make -j15 # Adjust for available core count

