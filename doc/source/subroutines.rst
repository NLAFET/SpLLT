***********
Subroutines
***********

Basic subroutines
=================

   
.. f:subroutine:: spllt_analyse(akeep, fkeep, options, n, ptr, row, info, order)

   Perform the analysis of the input matrix.

   :p spllt_akeep akeep [in]: symbolic factorization data.
   :p spllt_fkeep fkeep [in]: factorization data.
   :p spllt_options options [in]: user-supplied options to be used.
   :p integer n [in]:
   :p integer ptr(n+1) [in]:
   :p integer row(ptr(n+1)-1) [in]:
   :p spllt_inform info [out]: 
   :p integer order(n) [out]: The permutation array used after ordering.

.. f:subroutine:: spllt_factor(akeep, fkeep, options, val, info)

   Factorizes the input matrix.

.. f:subroutine:: spllt_solve(fkeep, options, order, x, info, job)

   Solves for multiple right-hand sides on the following problems:

   +---------------+--------------------------+
   | `job`         | Equation solved          |
   +===============+==========================+
   | 0 (or absent) | :math:`Ax=b`             |
   +---------------+--------------------------+
