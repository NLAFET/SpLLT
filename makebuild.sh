#! /bin/sh

# build OpenMP code
mkdir builds
mkdir builds/omp
cd builds/omp
cmake ../.. -DRUNTIME=OMP
make spllt

# make
RESULT=$?
[ $RESULT -ne 0 ] && exit 1

exit 0
