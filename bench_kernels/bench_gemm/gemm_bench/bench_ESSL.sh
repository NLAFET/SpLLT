#!/bin/bash
for X in S D C Z
do
	for ((N=1024; N<=16384; N+=1024))
	do
#		OMP_NUM_THREADS=1 OMP_PROC_BIND=close OMP_PLACES=cores /usr/bin/numactl -m 0 -N 0 ./"$X"gemm_bench-essl.exe 10 N N $N $N $N 1.0 1.0 >> "$X"gemm_bench-essl.csv
		for ((T=1; T<=8; ++T))
		do
			OMP_NUM_THREADS=$T OMP_PROC_BIND=close OMP_PLACES=cores /usr/bin/numactl -m 0 -N 0 ./"$X"gemm_bench-essl_smp.exe 10 N N $N $N $N 1.0 1.0 >> "$X"gemm_bench-essl_smp_"$T".csv
		done
	done
done
