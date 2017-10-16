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

BLAS
----

The solver depends on high performance BLAS library to perform the
dense linear algebra operations that are abundant in our kernels. For
best performance, please use the library recommended by your computer
manufacturer (normally the Intel MKL). If this is not available, use
an optimized alternative, such as OpenBLAS.  The reference BLAS from
netlib are at least an order of magnitude slower than modern optimized
BLAS, and should be avoided. If bit-compatible results are desired, a
bit-compatible BLAS library must be used.

SPRAL
-----

SPRAL is an open source (BSD) library for sparse linear algebra and
associated algorithms. It can be downloaded directly from the SPRAL
GitHub repository: `<https://github.com/ralna/spral>`_.

Runtime system
--------------

In this package we used a runtime system for implementing the parallel
version of our code. SpLLT currently supports three runtime systems
among `OpenMP <http://www.openmp.org/>`_, `StarPU
<http://starpu.gforge.inria.fr/>`_ and `Parsec
<https://bitbucket.org/icldistcomp/parsec>`_ that can be set using the
:code:`-DRUNTIME` option when running the cmake command.

Support
=======
Feeback may be sent to `florent.lopez@stfc.ac.uk <florent@stfc.ac.uk>`_ or by filing
an issue on our github: `<https://github.com/NLAFET/SpLLT/issues>`_.
