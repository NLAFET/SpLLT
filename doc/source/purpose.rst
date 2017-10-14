*******
Purpose
*******

SpLLT is a direct solver for solving **large sparse positive-definite
symmetric linear systems** of equations:

.. math::

   AX=B
          
this is done by computing the Cholesky decomposition of the input
matrix:

.. math::
   PAP^T=LL^T

where :math:`P` is a permutation matrix and the factor :math:`L` is a
lower triangular matrix. Following the matrix factorization the
solution can be retrieved by successively solving the system
:math:`LY=PB` (forward substitution) and :math:`L^{T}PX=Y` (backward
substitutions).

