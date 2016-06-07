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

                # ./run_ma87
                echo "[run_tests] run SPLLT_STARPU"
                rm -rf $trace_dir/$prof_file
                # just to make sure
                export OMP_NUM_THREADS=1
                case $build in
                    ma87)
                        echo "[run_tests] run MA87"
                        export OMP_NUM_THREADS=${ncpu} 
                        ./run_ma87 --ncpu ${ncpu} --nb ${nb} --nemin ${nemin} > $outdir/ma87/${matname}_NCPU-${ncpu}_NB-${nb}_NEMIN-${nemin}
                        ;;
                    parsec)
                        ./spllt_parsec_test --ncpu ${ncpu} --nb ${nb} --nemin ${nemin} > $outdir/spllt_parsec/${matname}_NCPU-${ncpu}_NB-${nb}_NEMIN-${nemin}${outsuffix}
                        ;;
                    starpu)
                        ./spllt_starpu_test --ncpu ${ncpu} --nb ${nb} --nemin ${nemin} > $outdir/spllt_starpu/${matname}_NCPU-${ncpu}_NB-${nb}_NEMIN-${nemin}${outsuffix}
                        ;;
                    gnu_omp)
                        ./spllt_omp_test --ncpu ${ncpu} --nb ${nb} --nemin ${nemin} > $outdir/spllt_omp/gnu/${matname}_NCPU-${ncpu}_NB-${nb}_NEMIN-${nemin}${outsuffix}
                        ;;
                    intel_omp)
                        ./spllt_omp_test --ncpu ${ncpu} --nb ${nb} --nemin ${nemin} > $outdir/spllt_omp/intel/${matname}_NCPU-${ncpu}_NB-${nb}_NEMIN-${nemin}${outsuffix}
                        ;;
                esac

                if [ -f $trace_dir/$prof_file ];
                then
                    mv $trace_dir/$prof_file $outdir/spllt_starpu/traces/${matname}_NCPU-${ncpu}_NB-${nb}.prof
                    starpu_fxt_tool -c -i $outdir/spllt_starpu/traces/${matname}_NCPU-${ncpu}_NB-${nb}.prof -o $outdir/spllt_starpu/traces/${matname}_NCPU-${ncpu}_NB-${nb}.trace
                fi
            done
        done
    done
done
