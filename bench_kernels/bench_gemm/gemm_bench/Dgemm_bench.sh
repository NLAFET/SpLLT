#!/bin/bash
# This shows how to link with SMP ESSL on Panther.
# Tested with GNU Fortran (GCC) 6.1.1 20160728
# Notice -fno-underscoring option!
# Cf. -qnoextname in XL Fortran documentation.
gfortran -DGEMM_D -fno-underscoring gemm_bench.F90 -o Dgemm_bench.exe -lesslsmp -Wl,-rpath=/gpfs/panther/local/apps/ibm/lib -L/gpfs/panther/local/apps/ibm/xlsmp/4.1.5/lib -lxlsmp -L/gpfs/panther/local/apps/ibm/xlf/15.1.5/lib -lxlf90_r -lxlfmath -lxl
