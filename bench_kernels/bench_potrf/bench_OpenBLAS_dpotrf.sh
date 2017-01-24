#!/bin/bash
export DATA_DIR=../../../ralna/TR-SpLLT_GPU/data/`hostname`
echo "DATA_DIR=$DATA_DIR"
#export OMP_PLACES="0,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30" # CPU0
#export OMP_PLACES="0,2,4,6,8,10,12,14" # CPU0,THR0
export OMP_PLACES="{16,18,20,22,24,26,28,30}" # CPU0,THR1
export OMP_PROC_BIND=spread
echo "OMP_PLACES=$OMP_PLACES"
echo "OMP_PROC_BIND=$OMP_PROC_BIND"
for ((I=$5; I<=$6; ++I))
do
	export OMP_NUM_THREADS=$I
	echo "OMP_NUM_THREADS=$OMP_NUM_THREADS"
	export OUT_CSV=OpenBLAS_dpotrf-$I-$1-$2-$3-$4-`date +%s`.csv
	echo "OUT_CSV=$OUT_CSV"
	/usr/bin/numactl -m 0 -N 0 ./bench_OpenBLAS_potrf.exe $1 $2 $3 $4 > $DATA_DIR/$OUT_CSV
	pushd $DATA_DIR
	git add -v $OUT_CSV
	popd
done
