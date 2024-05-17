## instructions for linking lammps to qmcpack (WIP)
1. compile lammps
 a shared library lammps can be compiled with:
```
#!/bin/bash

export MPICC=/usr/local/bin/mpicc
export MPICXX=/usr/local/bin/mpicxx
rm -r build ; mkdir build; cd build
cmake -DCMAKE_CXX_COMPILER=$MPICXX \
      -DCMAKE_C_COMPILER=$MPICC \
      -DCMAKE_INSTALL_PREFIX=$HOME/software/lammps/ \ // ensuring lammps pkg-config points to right place for qmcpack cmake usage.
      -DBUILD_MPI=yes \
      -DPKG_MOLECULE=yes \
      -DBUILD_SHARED_LIBS=yes ../cmake
cmake --build .          # compilation (or type "make")
```
2. point qmcpack to lammps
```
export CC=/usr/local/Cellar/gcc/11.3.0/bin/gcc-11
export CXX=/usr/local/Cellar/gcc/11.3.0/bin/g++-11
export MPICC=/usr/local/bin/mpicc
export MPICXX=/usr/local/bin/mpic++
export HDF5_ROOT=/usr/local/opt/hdf5
export FFTW_HOME=/usr/local/opt/fftw
export LAMMPS_ROOT=/Users/aedumi/software/lammps/build // important switch
export PKG_CONFIG_PATH="/usr/local/opt/readline/lib/pkgconfig"
export DYLD_LIBRARY_PATH=/usr/local/Cellar/gcc/11.3.0/lib/gcc/11:/Users/aedumi/software/lammps/build // this is mac specific. 
echo "configured the environment"
cd build
cmake -DCMAKE_C_COMPILER=$MPICC \
      -DCMAKE_CXX_COMPILER=$MPICXX \
      -DFFTW_INCLUDE_DIRS=$FFTW_HOME/include \
      -DHDF5_ROOT=$HDF5_ROOT \
      -DFFTW_LIBRARY_DIRS=$FFTW_HOME/lib \
      -DCMAKE_MACOSX_RPATH=ON \
      -DQMC_EXTRA_LIBS=/usr/local/Cellar/gcc/11.3.0/lib/gcc/11/libgomp.dylib  \
      -DLAMMPS_ROOT=$LAMMPS_ROOT \
      -DQMC_COMPLEX=1 ..
make -j 1 VERBOSE=1  # Adjust for available core count
```


4. checking include paths are right in cmake. Helpful information is printed during the cmake step.
5. running ctest for the unit tests in the build directory will run the test_lammps function if lammps was found. Currently, test is breaking compilation, but this is ideal end goal and current design.

