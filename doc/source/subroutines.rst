***********
Subroutines
***********

Basic subroutines
=================

   
.. f:subroutine:: spllt_analyse(akeep, fkeep, options, n, ptr, row, info, order)

   Perform the analyse phase of the factorization (referred to as
   symbolic factorization) for a matrix supplied in `Compressed Sparse
   Column (CSC) format
   <http://www.numerical.rl.ac.uk/spral/doc/latest/Fortran/csc_format.html>`_. The
   resulting symbolic factors stored in `akeep` should be passed
   unaltered in the subsequent calls to :f:subr:`spllt_factor()`. This
   routines also initializes the numeric factorization data stored in
   `fkeep` that should also be passed in the subsequent calls to
   :f:subr:`spllt_factor()`.

   :p spllt_akeep akeep [inout]: symbolic factorization data.
   :p spllt_fkeep fkeep [inout]: numeric factorization data. 
   :p spllt_options options [in]: user-supplied options to be used.
   :p integer n [in]: order of the system.
   :p integer ptr(n+1) [in]: column pointers for lower triangular part.
   :p integer row(ptr(n+1)-1) [in]: row indices of lower triangular part.
   :p spllt_inform info [out]: exit status.
   :p integer order(n) [out]: permutation array used for matrix
      ordering.

.. f:subroutine:: spllt_factor(akeep, fkeep, options, val, info)

   Perform the numerical factorization of the matrix whose structure
   is kept in `akeep` that is determined with the :f:subr:`spllt_analyse()`
   routine. The numerical factors are stored in `fkeep` and should be
   passed unaltered in the subsequent calls to :f:subr:`spllt_solve()` for
   computing the solution to the system.

   :p spllt_akeep akeep [inout]: symbolic factorization data.
   :p spllt_fkeep fkeep [inout]: numeric factorization data. 
   :p spllt_options options [in]: user-supplied options to be used.
   :p real val(*) [in]: non-zero values for :math:`A` in same format as for
      the call to :f:subr:`ssids_analyse()`
   :p spllt_inform info [out]: exit status.

.. note::

   This routine call is asynchronous, the routine
   :f:subr:`spllt_wait()` should be call afterwards to make sure that
   the factorization has been completed.

.. f:subroutine:: spllt_solve(fkeep, options, order, x, info[, job])

   Solves for multiple right-hand sides on the following problems:

   +---------------+--------------------------+
   | `job`         | Equation solved          |
   +===============+==========================+
   | 0 (or absent) | :math:`AX=B`             |
   +---------------+--------------------------+
   | 1             | :math:`PLX=B`            |
   +---------------+--------------------------+
   | 2             | :math:`(PL)^TX=B`        |
   +---------------+--------------------------+

   :p spllt_fkeep fkeep [in]: numeric factorization data. 
   :p spllt_options options [in]: user-supplied options to be used.
   :p integer order(n) [out]: permutation array used for matrix
      ordering.
   :p real x(n,nrhs) [inout]: right-hand sides :math:`B` on entry,
      solutions :math:`X` on exit. `n` represents the order of the
      system to solve.
   :p spllt_inform info [out]: exit status.
   :o integer job [in]: specifies equation to solve, as per above table.

.. f:subroutine:: spllt_wait()

   Wait for all the tasks submitted in a previous function call to be
   completed.
