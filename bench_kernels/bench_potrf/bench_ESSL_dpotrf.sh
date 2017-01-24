#!/bin/bash
export DATA_DIR=.
echo "DATA_DIR=$DATA_DIR"
export OMP_PLACES=cores
export OMP_PROC_BIND=CLOSE
echo "OMP_PLACES=$OMP_PLACES"
echo "OMP_PROC_BIND=$OMP_PROC_BIND"
for ((I=$5; I<=$6; ++I))
do
	export OMP_NUM_THREADS=$I
	echo "OMP_NUM_THREADS=$OMP_NUM_THREADS"
	export OUT_CSV=ESSL_dpotrf-$I-$1-$2-$3-$4-`date +%s`.csv
	echo "OUT_CSV=$OUT_CSV"
	./bench_ESSL_potrf.exe $1 $2 $3 $4 > $DATA_DIR/$OUT_CSV
done
