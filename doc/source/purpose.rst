*******
Purpose
*******

SpLLT is a direct solver for solving **large sparse symmetric
positive-definite linear systems** of equations:

.. math::

   AX=B
          
This is done by computing the `Cholesky decomposition
<https://en.wikipedia.org/wiki/Cholesky_decomposition>`_ of the input
matrix:

.. math::
   PAP^T=LL^T

where the factor :math:`L` is a lower triangular matrix and the matrix
:math:`P` is a permutation matrix used to reduce the `fill-in
<https://en.wikipedia.org/wiki/Sparse_matrix#Reducing_fill-in>`_
generated during the factorization. Following the matrix factorization
the solution can be retrieved by successively solving the system
:math:`LY=PB` (forward substitution) and :math:`L^{T}PX=Y` (backward
substitutions).

