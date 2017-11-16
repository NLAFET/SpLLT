***********
Subroutines
***********

Basic subroutines
=================

   
.. f:subroutine:: spllt_analyse(akeep, fkeep, options, n, ptr, row, info, order)

   Perform the analysis of the input matrix.

.. f:subroutine:: spllt_factor(akeep, fkeep, options, val, info)

   Factorizes the input matrix.

.. f:subroutine:: spllt_solve()

   Solves for multiple right-hand sides on the following problems:

   +---------------+--------------------------+
   | `job`         | Equation solved          |
   +===============+==========================+
   | 0 (or absent) | :math:`Ax=b`             |
   +---------------+--------------------------+
