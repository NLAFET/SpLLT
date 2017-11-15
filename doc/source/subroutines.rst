***********
Subroutines
***********

Basic subroutines
=================

   
.. f:subroutine:: spllt_analyse(adata, fdata, options, n, ptr, row, info, order)

   Perform the analysis of the input matrix.

.. f:subroutine:: spllt_factor(adata, fdata, options, val, info)

   Factorizes the input matrix.

.. f:subroutine:: spllt_solve()

   Solves for multiple right-hand sides on the following problems:

   +---------------+--------------------------+
   | `job`         | Equation solved          |
   +===============+==========================+
   | 0 (or absent) | :math:`Ax=b`             |
   +---------------+--------------------------+
