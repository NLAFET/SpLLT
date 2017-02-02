#!/bin/bash
for X in S D C Z
do
	xlf2008_r -WF,-DGEMM_$X -qfloat=subnormals -qmaxmem=-1 -qnosave -qsclk=micro -qarch=pwr8 -qtune=pwr8:smt8 -qsmp=omp gemm_bench.F90 -o "$X"gemm_bench-essl.exe -lessl
	xlf2008_r -WF,-DGEMM_$X -qfloat=subnormals -qmaxmem=-1 -qnosave -qsclk=micro -qarch=pwr8 -qtune=pwr8:smt8 -qsmp=omp gemm_bench.F90 -o "$X"gemm_bench-essl_smp.exe -lesslsmp
done
