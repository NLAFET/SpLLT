#!/bin/bash
for X in S D C Z
do
	xlf2008_r -WF,-DGEMM_$X -I/gpfs/cds/local/SCD/jpf02/vxn61-jpf02/ppc64le/include -qfloat=subnormals -qmaxmem=-1 -qnosave -qsclk=micro -qarch=pwr8 -qtune=pwr8:smt8 -qsmp=omp gemm_bench.F90 -o "$X"gemm_bench-essl.exe -lessl /gpfs/cds/local/SCD/jpf02/vxn61-jpf02/ppc64le/lib/libpapi.a
	xlf2008_r -WF,-DGEMM_$X -I/gpfs/cds/local/SCD/jpf02/vxn61-jpf02/ppc64le/include -qfloat=subnormals -qmaxmem=-1 -qnosave -qsclk=micro -qarch=pwr8 -qtune=pwr8:smt8 -qsmp=omp gemm_bench.F90 -o "$X"gemm_bench-essl_smp.exe -lesslsmp /gpfs/cds/local/SCD/jpf02/vxn61-jpf02/ppc64le/lib/libpapi.a
done
