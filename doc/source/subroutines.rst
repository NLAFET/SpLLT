***********
Subroutines
***********

Basic subroutines
=================

   
.. f:subroutine:: spllt_analyse()

   Perform the analysis of the input matrix.

.. f:subroutine:: spllt_factor()

   Factorizes the input matrix.

.. f:subroutine:: spllt_solve()

   Solves for multiple right-hand sides on the following problems:

   +---------------+--------------------------+
   | `job`         | Equation solved          |
   +===============+==========================+
   | 0 (or absent) | :math:`Ax=b`             |
   +---------------+--------------------------+
