#!/bin/bash

module purge
module load use.own
module load gcc/5.3.0
module load intel/mkl/11.3.1.150
module load hwloc/1.11.2
module load starpu/trunk-nogpu-nofxt
module load metis/4.0.3
module load hsl/latest
module load spral/trunk
module list

build_dir=`pwd`
id=`whoami`
outdir=data

echo "[run_tests] build dir: $build_dir"
#matrices=(JGD_Trefethen/Trefethen_20000)

ncpu_list=(14 27)
nb=256
trace_dir=/tmp
prof_file=prof_file_scarf462_0
# for matrix in ${matrices[@]}

mkdir -p $outdir
mkdir -p $outdir/spllt_starpu
mkdir -p $outdir/spllt_starpu/traces
mkdir -p $outdir/ma87

for matrix in `cat list.matrix`
do
    matname=`echo $matrix | sed -e 's/\//_/g'`
    echo "[run_tests] test matrix: $matname"
    # set up matrix
    echo "[run_tests] extract matrix: $matname"
    ./prep.sh $matrix

    for ncpu in ${ncpu_list[@]}
    do
        echo "[run_tests] ncpu: $ncpu"
        # ./run_ma87
        echo "[run_tests] run MA87"
        export OMP_NUM_THREADS=${ncpu} 
        ./run_ma87 --ncpu ${ncpu} --nb ${nb} > $outdir/ma87/${matname}_NCPU-${ncpu}_NB-${nb}
        echo "[run_tests] run SPLLT_STARPU"
        rm -rf $trace_dir/$prof_file
        # just to make sure
        export OMP_NUM_THREADS=1
        ./spllt_starpu_test --ncpu ${ncpu} --nb ${nb} > $outdir/spllt_starpu/${matname}_NCPU-${ncpu}_NB-${nb}
        
        if [ -f $trace_dir/$prof_file ];
        then
            mv $trace_dir/$prof_file $outdir/spllt_starpu/traces/${matname}_NCPU-${ncpu}_NB-${nb}.prof
            starpu_fxt_tool -c -i $outdir/spllt_starpu/traces/${matname}_NCPU-${ncpu}_NB-${nb}.prof -o $outdir/spllt_starpu/traces/${matname}_NCPU-${ncpu}_NB-${nb}.trace
        fi
    done
done
