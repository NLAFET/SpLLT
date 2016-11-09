# SpLLT

This is a DAG-based sparse Cholesky solver which uses a runtime-based
approach. It currently supports the
[StarPU](http://starpu.gforge.inria.fr/) and
[PaRSEC](https://bitbucket.org/icldistcomp/parsec) runtime systems as
well as the [OpenMP standard](http://openmp.org/) Version 4.0 or
above.

# Runtime systems

By default the code is compiled in sequential mode but the choice of
the runtime system can be specified by setting the option
`-DRUNTIME`. For example, the sequential code can be configured as
following:

```bash
cmake -DRUNTIME=STF <path-to-source>

```

## StarPU

A parallel version of the code using the
[StarPU](http://starpu.gforge.inria.fr/) runtime system can be
obtained as following:

```bash
cmake -DRUNTIME=StarPU <path-to-source>

```

## Parsec

A parallel version of the code using the
[PaRSEC](https://bitbucket.org/icldistcomp/parsec) runtime system can
be obtained as following:

```bash
cmake -DRUNTIME=Parsec <path-to-source>

```

