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
[SPRAL](http://www.numerical.rl.ac.uk/spral/), where the sources
are [here](https://github.com/ralna/spral).
This software requires itself a graph partitioner as
[Metis](http://glaros.dtc.umn.edu/gkhome/) library, and [HWLOC](https://
www.open-mpi.org/projects/hwloc/).

## Metis 4.0

During the analyse, a graph partitioner as Metis needs to be linked to SPRAL.
The latest current version supported by SPRAL is [metis-4.0.3](http://glaros.
dtc.umn.edu/gkhome/fetch/sw/metis/OLD/metis-4.0.3.tar.gz).
Download the source code, open Makefile.in, and set the variable CC to gcc.
Then,

```bash
make
mkdir lib_gcc-<compiler_version>
cp libmetis.a lib_gcc-<compiler_version>
```

The same installation can be done with Intel compiler (changing CC variable in 
Makefile.in and adapting the name of the created folder).

## HWLOC 2.0

SPRAL requires Hardware Locality library where the sources are [here](https://
www.open-mpi.org/software/hwloc/v2.0/).

```bash
./configure
make
mkdir lib_gcc-<compiler_version>
cp hwloc/.libs/* lib_gcc-6.2.0/
```

## SPRAL

To install SPRAL, download the sources [here](https://github.com/ralna/spral).
Then

```bash
./configure --disable-openmp --disable-gpu --disable-openmp --with-blas="-L$MKL_LIB -lmkl_core -lmkl_intel_lp64" --with-lapack="-L$MKL_LIB -lmkl_core -lmkl_intel_lp64" --with-metis="-L$METIS_LIB -lmetis"
make
mkdir lib_gcc-<compiler_version>
cp libspral.a lib_gcc-<compiler_version>
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
