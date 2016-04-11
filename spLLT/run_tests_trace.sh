#!/bin/bash

# list loaded modules
module list

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

        for nb in ${nb_list[@]}
        do
            echo "[run_tests] nb: $nb"

            for nemin in ${nemin_list[@]}
            do  
                echo "[run_tests] nemin: $nemin"

                rm -rf $trace_dir/$prof_file
                # just to make sure
                export OMP_NUM_THREADS=1
                # launch test
                ./spllt_starpu_test --ncpu ${ncpu} --nb ${nb} --nemin ${nemin}
                # generate trace
                if [ -f $trace_dir/$prof_file ];
                then
                    mv $trace_dir/$prof_file $outdir/spllt_starpu/traces/${matname}_NCPU-${ncpu}_NB-${nb}_NEMIN-${nemin}${outsuffix}.prof
                    starpu_fxt_tool -c -i $outdir/spllt_starpu/traces/${matname}_NCPU-${ncpu}_NB-${nb}_NEMIN-${nemin}${outsuffix}.prof -o $outdir/spllt_starpu/traces/${matname}_NCPU-${ncpu}_NB-${nb}_NEMIN-${nemin}${outsuffix}.trace
                fi
            done
        done
    done
done
