# SpLLT - sparse Cholesky sovler

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

