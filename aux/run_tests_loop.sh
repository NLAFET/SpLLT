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
                        ../builds/ma87/spllt_test --ncpu ${ncpu} --nb ${nb} --nemin ${nemin} > $outdir/ma87/${matname}_NCPU-${ncpu}_NB-${nb}_NEMIN-${nemin}${outsuffix}
                        ;;
                    stf)
                        ../builds/stf/spllt_test --ncpu 1 --nb ${nb} --nemin ${nemin} > $outdir/stf/${matname}_NB-${nb}_NEMIN-${nemin}${outsuffix}
                        ;;
                    parsec)
                        ../builds/parsec/spllt_test --ncpu ${ncpu} --nb ${nb} --nemin ${nemin} > $outdir/parsec/${matname}_NCPU-${ncpu}_NB-${nb}_NEMIN-${nemin}${outsuffix}
                        ;;
                    starpu)
                        ../builds/starpu/spllt_test --ncpu ${ncpu} --nb ${nb} --nemin ${nemin} > $outdir/starpu/${matname}_NCPU-${ncpu}_NB-${nb}_NEMIN-${nemin}${outsuffix}
                        ;;
                    starpu_prune)
                        ../builds/starpu/spllt_test --prune-tree --ncpu ${ncpu} --nb ${nb} --nemin ${nemin} > $outdir/starpu/${matname}_NCPU-${ncpu}_NB-${nb}_NEMIN-${nemin}${outsuffix}
                        ;;
                    starpu_nested_stf)
                        ../builds/starpu-nested-stf/spllt_test --ncpu ${ncpu} --nb ${nb} --nemin ${nemin} > $outdir/starpu_nested_stf/${matname}_NCPU-${ncpu}_NB-${nb}_NEMIN-${nemin}${outsuffix}
                        ;;
                    starpu_trace)
                        ../builds/starpu-trace/spllt_test --ncpu ${ncpu} --nb ${nb} --nemin ${nemin} > $outdir/starpu/${matname}_NCPU-${ncpu}_NB-${nb}_NEMIN-${nemin}${outsuffix}
                        ;;
                    starpu_gpu)
                        ../builds/starpu-gpu/spllt_test --ncpu ${ncpu} --nb ${nb} --nemin ${nemin} > $outdir/starpu/${matname}_NCPU-${ncpu}_NB-${nb}_NEMIN-${nemin}${outsuffix}
                        ;;
                    gnu_omp)
                        ../builds/omp/gnu/spllt_test --ncpu ${ncpu} --nb ${nb} --nemin ${nemin} > $outdir/omp/gnu/${matname}_NCPU-${ncpu}_NB-${nb}_NEMIN-${nemin}${outsuffix}
                        ;;
                    gnu_omp_prune)
                        ../builds/omp/gnu/spllt_test --prune-tree --ncpu ${ncpu} --nb ${nb} --nemin ${nemin} > $outdir/omp/gnu/${matname}_NCPU-${ncpu}_NB-${nb}_NEMIN-${nemin}${outsuffix}
                        ;;
                    intel_omp)
                        ../builds/omp//spllt_test --ncpu ${ncpu} --nb ${nb} --nemin ${nemin} > $outdir/omp/intel/${matname}_NCPU-${ncpu}_NB-${nb}_NEMIN-${nemin}${outsuffix}
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
