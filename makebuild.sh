#! /bin/sh

# build OpenMP code
mkdir builds
mkdir builds/omp
cd builds/omp
cmake ../.. -DRUNTIME=OMP -DBLAS_LAPACK=Ref
make

# make
RESULT=$?
[ $RESULT -ne 0 ] && exit 1

exit 0
