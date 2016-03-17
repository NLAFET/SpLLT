#!/bin/bash

build_dir=`pwd`
outdir=data

echo "[run_tests] build dir: $build_dir"
#matrices=(JGD_Trefethen/Trefethen_20000)

ncpu=4
nb=256
trace_dir=/home/flopez/traces
prof_file=prof_file_flopez_0
# for matrix in ${matrices[@]}
for matrix in `cat list.matrix`
do
    matname=`echo $matrix | sed -e 's/\//_/g'`
    echo "[run_tests] test matrix: $matname"
    # set up matrix
    ./prep.sh $matrix
    echo "[run_tests] run MA87"
    # ./run_ma87
    ./run_ma87 > $outdir/ma87/${matname}_NCPU-${ncpu}_NB-${nb}
    echo "[run_tests] run SPLLT_STARPU"
    rm -rf $trace_dir/$prof_file
    ./spllt_starpu_test --ncpu ${ncpu} --nb ${nb} > $outdir/spllt_starpu/${matname}_NCPU-${ncpu}_NB-${nb}
    
    if [ -f $trace_dir/$prof_file ];
    then
        starpu_fxt_tool -c -i $trace_dir/$prof_file
        mv $trace_dir/paje.trace data/spllt_starpu/traces/${matname}_NCPU-${ncpu}_NB-${nb}.trace
    fi
done
