#!/bin/bash
for X in S D C Z
do
	ifort -O3 -xHost -fpp -qopenmp -mkl=sequential -DGEMM_$X gemm_bench.F90 -o "$X"gemm_bench-mkl_seq.exe
	ifort -O3 -xHost -fpp -qopenmp -mkl=parallel   -DGEMM_$X gemm_bench.F90 -o "$X"gemm_bench-mkl_par.exe
done
