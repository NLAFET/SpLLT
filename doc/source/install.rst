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
available, use an optimized alternative, such as OpenBLAS.  The
reference BLAS from netlib are at least an order of magnitude slower
than modern optimized BLAS, and should be avoided. If bit-compatible
results are desired, a bit-compatible BLAS library must be used.

Support
=======
Feeback may be sent to `florent.lopez@stfc.ac.uk <florent@stfc.ac.uk>`_ or by filing
an issue on our github: `<https://github.com/NLAFET/SpLLT/issues>`_.
