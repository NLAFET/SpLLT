************
Installation
************

Quick Start
===========

Under Linux, or Mac OS X:

.. code-block:: bash
   
   # Get latest development version from github
   git clone https://github.com/NLAFET/SpLLT
   # Move to source directory
   cd SpLLT 
   # Create build directory
   mkdir build 
   cd build
   # Create Makefile with cmake command. The -DRUNTIME option can be 
   # used to select a runtime system.
   cmake <path-to-source> -DRUNTIME=<StarPU|OMP|Parsec>
   # Build SpLLT software
   make

Third-party libraries
=====================

BLAS and LAPACK
---------------

The solver depends on high performance BLAS and LAPACK libraries to
perform the dense linear algebra operations that are abundant in our
kernels. For best performance, please use the library recommended by
your computer manufacturer (normally the Intel MKL). If this is not
available, use an optimized alternative, such as `OpenBLAS
<http://www.openblas.net/>`_.  The `reference BLAS
<http://www.netlib.org/blas/>`_ and `reference LAPACK
<http://www.netlib.org/lapack/>`_ libraries from netlib are at least
an order of magnitude slower than modern optimized BLAS, and should be
avoided. If bit-compatible results are desired, a bit-compatible BLAS
library must be used.

The BLAS library can be passed to `cmake` using the :code:`-DLBLAS`
option and the LAPACK library can be passed using the
:code:`-DLLAPACK` option as following

.. code-block:: bash

   cmake <path-to-source> -DLBLAS=/path/to/blas -DLLAPACK=/path/to/lapack

SPRAL
-----

SPRAL is an open source (BSD) library for sparse linear algebra and
associated algorithms. It can be downloaded directly from the SPRAL
GitHub repository: `<https://github.com/ralna/spral>`_.


.. code-block:: bash


Runtime system
--------------

In this package we used a runtime system for implementing the parallel
version of our code. SpLLT currently supports three runtime systems
among `OpenMP <http://www.openmp.org/>`_, `StarPU
<http://starpu.gforge.inria.fr/>`_ and `Parsec
<https://bitbucket.org/icldistcomp/parsec>`_ that can be set using the
:code:`-DRUNTIME` option when running the cmake command.

Compilers and compiler options
==============================
If no compiler is specified, `cmake` will pick a default
compiler to use. If `cmake` cannot find an appropriate compiler, or
you wish to specify a different compiler you can do so by setting the following
variables:

CC
  specifies the C compiler to use.
FC
  specifies the Fortran 90/95/2003/2008 compiler to use.
NVCC
  specifies the CUDA compiler to use.

Additionally, compiler flags can be specified using the following variables:

CFLAGS
   specifies options passed to the C compiler.
FCFLAGS
   specifies options passed to the Fortran compiler

NVCCFLAGS
   specifies options passed to the CUDA compiler.

For example, to compile with ``ifort -g -O3 -ip`` we could use:

.. code-block:: bash

   FC=ifort FCFLAGS="-g -O3 -ip" cmake <path-to-source>

      
Support
=======
Feeback may be sent to `florent.lopez@stfc.ac.uk <florent@stfc.ac.uk>`_ or by filing
an issue on our github: `<https://github.com/NLAFET/SpLLT/issues>`_.
