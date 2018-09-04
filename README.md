# SpLLT - sparse Cholesky solver

|OpenMP|
|------|
|[![Build Status](https://travis-ci.com/NLAFET/SpLLT.svg?token=UwhpFb953M8N7PyHRDWG&branch=master)](https://travis-ci.com/NLAFET/SpLLT)|

SpLLT is a sparse direct solver for computing the solution of
[symmetric positive
definite](https://en.wikipedia.org/wiki/Positive-definite_matrix)
linear systems. The factorization phase, which is the most
computationally intensive part, is based on a task-based Cholesky
algorithm and the parallel code is implemented using a runtime
system. The code supports three different runtime systems:
[StarPU](http://starpu.gforge.inria.fr/) developed at INRIA Bordeaux
Sud-Ouest, [PaRSEC](https://bitbucket.org/icldistcomp/parsec) from
ICL, University of Tennessee and the [OpenMP
standard](http://openmp.org/) Version 4.0 or above.

The compilation is handled by [CMake](https://cmake.org/) tools and
the solver can be built as following instructions:

```bash
mkdir build # create build directory
cd build
cmake <path-to-source>
```

This will build the sequential version of the code by default. In
order to build the parallel code you need to select a runtime system
using the option `-DRUNTIME` when running the cmake `cmake
<path-to-source>` command.

# Prerequisites

In order to install properly SpLLT, some libraries are required.
The analysis of the matrix is performed through
[SPRAL](http://www.numerical.rl.ac.uk/spral/).
This software requires itself a graph partitioner as
[Metis](http://glaros.dtc.umn.edu/gkhome/) library, 
and [HWLOC](https://www.open-mpi.org/projects/hwloc/).

## Metis 5.1

During the analyse, a graph partitioner as Metis needs to be linked to SPRAL.
The latest current version supported by SPRAL is Metis-5.1 available 
[here](http://glaros.dtc.umn.edu/gkhome/metis/metis/download).
Download the source code, open Makefile.in, and set the variable CC to gcc.
Then,

```bash
make config prefix=<path_to_install> shared=1 cc=<compiler_requested>
make
make install
```

The same installation can be done with Intel compiler (changing CC variable in 
Makefile.in and adapting the name of the created folder).

## HWLOC 2.0

SPRAL requires Hardware Locality library where the sources are [here](https://www.open-mpi.org/software/hwloc/v2.0/).

```bash
CC=<C_compiler> CXX=<CXX_compiler> ./configure --prefix=<path_to_install>
make
make install
```

## SPRAL

To install SPRAL, download the sources [here](https://github.com/ralna/spral).
Then

```bash
./autogen.sh
CC=icc CXX=icpc FC=ifort ./configure --prefix=<path_to_install> --disable-openmp --disable-gpu --with-blas="-L$MKL_LIB -lmkl_gf_lp64 -lmkl_core -lmkl_intel_thread -liomp5 -lm" --with-lapack="-L$MKL_LIB -lmkl_gf_lp64 -lmkl_core -lmkl_intel_thread -liomp5 -lm" --with-metis="-L$METIS_LIB -lmetis"
make
make install
cp *.mod <path_to_install>/include
```

# RUNTIME

SpLLT is designed to be used with different runtimes.

## OMP

A parallel version of the code using the [OpenMP](https://openmp.org/)
standard with tasking capabilities can be obtained as following:

```bash
cmake -DRUNTIME=OMP <path-to-source>

```

## StarPU

A parallel version of the code using the
[StarPU](http://starpu.gforge.inria.fr/) runtime system can be
obtained as following:

```bash
cmake -DRUNTIME=StarPU <path-to-source>

```

## PaRSEC

A parallel version of the code using the
[PaRSEC](https://bitbucket.org/icldistcomp/parsec) runtime system can
be obtained as following:

```bash
cmake -DRUNTIME=Parsec <path-to-source>

```

# Examples of installation of SpLLT

To install SpLLT, for each prerequisite you can either use the environment or 
use defined variables, and provide either the path to the directory or the paths
to the library and include folders.
For example, considering GNU compilers and MKL, we can explicitly define to 
cmake the path of the library and include folders for each prerequisite, 
assuming that these paths are set in the environment.

```bash
export SPRAL_LIB=<path_to_spral_library_folder>
export SPRAL_INC=<path_to_spral_include_folder>
export METIS_LIB=<path_to_metis_library_folder>
export METIS_INC=<path_to_metis_include_folder>
export HWLOC_LIB=<path_to_hwloc_library_folder>
export HWLOC_INC=<path_to_hwloc_include_folder>

mkdir build
mkdir build/build_omp
cd build/build_omp

CC=gcc FC=gfortran CXX=g++ cmake -DRUNTIME=OMP -DSPRAL_LIB=${SPRAL_LIB} -DSPRAL_INC=${SPRAL_INC} -DMETIS_LIB=${METIS_LIB} -DMETIS_INC=${METIS_INC} -DHWLOC_LIB=${HWLOC_LIB} -DHWLOC_INC=${HWLOC_INC} -DLBLAS="${MKL_LIB}/libmkl_gf_lp64.a;${MKL_LIB}/libmkl_sequential.a;${MKL_LIB}/libmkl_core.a" -DLLAPACK="${MKL_LIB}/libmkl_gf_lp64.a;${MKL_LIB}/libmkl_sequential.a;${MKL_LIB}/libmkl_core.a" ../..
make
```
or
```bash
export SPLLT_MKL_BLAS_LAPACK_LIBS="${MKL_LIB}/libmkl_gf_lp64.a;${MKL_LIB}/libmkl_sequential.a;${MKL_LIB}/libmkl_core.a"

mkdir build
mkdir build/build_omp
cd build/build_omp

CC=gcc FC=gfortran CXX=g++ cmake -DRUNTIME=OMP -DLBLAS=${SPLLT_MKL_BLAS_LAPACK_LIBS} -DLLAPACK=${SPLLT_MKL_BLAS_LAPACK_LIBS} ../..
make
```
You can also consider only the directory of each prerequisite, as follow
```bash
export SPRAL_DIR=<path_to_spral_folder>
export METIS_DIR=<path_to_metis_folder>
export HWLOC_DIR=<path_to_hwloc_folder>
export SPLLT_MKL_BLAS_LAPACK_LIBS="${MKL_LIB}/libmkl_gf_lp64.a;${MKL_LIB}/libmkl_sequential.a;${MKL_LIB}/libmkl_core.a"

mkdir build
mkdir build/build_omp
cd build/build_omp

CC=gcc FC=gfortran CXX=g++ cmake -DRUNTIME=OMP -DLBLAS=${SPLLT_MKL_BLAS_LAPACK_LIBS} -DLLAPACK=${SPLLT_MKL_BLAS_LAPACK_LIBS} ../..
make
```
