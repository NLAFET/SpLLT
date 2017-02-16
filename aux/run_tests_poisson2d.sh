#!/bin/bash

# list loaded modules
module list

for msz in ${msz_list[@]}
do

    matname="poisson2d_${msz}x${msz}"
    matrix="${matname}.mtx"

    echo $matrix

    for ncpu in ${ncpu_list[@]}
    do    
        echo "[run_tests] ncpu: $ncpu"

        for nb in ${nb_list[@]}
        do
            echo "[run_tests] nb: $nb"
            
            for nemin in ${nemin_list[@]}
            do  
                echo "[run_tests] nemin: $nemin"
                
                export OMP_NUM_THREADS=1
                case $build in
                    ma87)
                        echo "[run_tests] run MA87"
                        export OMP_NUM_THREADS=${ncpu}
                        ../builds/ma87/spllt_test --mm ${matrix} --ncpu ${ncpu} --nb ${nb} --nemin ${nemin} > $outdir/ma87/${matname}_NCPU-${ncpu}_NB-${nb}_NEMIN-${nemin}${outsuffix}
                        ;;

                    # stf)
                    #     ./builds/stf/spllt_test --ncpu 1 --nb ${nb} --nemin ${nemin} > $outdir/spllt_stf/${matname}_NB-${nb}_NEMIN-${nemin}${outsuffix}
                    #     ;;
                    # parsec)
                    #     ./builds/parsec/spllt_test --ncpu ${ncpu} --nb ${nb} --nemin ${nemin} > $outdir/spllt_parsec/${matname}_NCPU-${ncpu}_NB-${nb}_NEMIN-${nemin}${outsuffix}
                    #     ;;

                    starpu)
                        ../builds/starpu/spllt_test --mm ${matrix} --ncpu ${ncpu} --nb ${nb} --nemin ${nemin} > $outdir/spllt_starpu/${matname}_NCPU-${ncpu}_NB-${nb}_NEMIN-${nemin}${outsuffix}
                        # ../builds/starpu/spllt_test --mm ${matrix} --ncpu ${ncpu} --nb ${nb} --nemin ${nemin} > $outdir/spllt_starpu/${matname}_NCPU-${ncpu}_NB-${nb}_NEMIN-${nemin}${outsuffix}
                        ;;

                    starpu_prune)
                        ../builds/starpu/spllt_test --prune-tree --mm ${matrix} --ncpu ${ncpu} --nb ${nb} --nemin ${nemin} > $outdir/spllt_starpu/${matname}_NCPU-${ncpu}_NB-${nb}_NEMIN-${nemin}${outsuffix}
                        ;;

                    # starpu_nested_stf)
                    #     ./builds/starpu-nested-stf/spllt_test --ncpu ${ncpu} --nb ${nb} --nemin ${nemin} > $outdir/spllt_starpu_nested_stf/${matname}_NCPU-${ncpu}_NB-${nb}_NEMIN-${nemin}${outsuffix}
                    #     ;;
                    # starpu_trace)
                    #     ./builds/starpu-trace/spllt_test --ncpu ${ncpu} --nb ${nb} --nemin ${nemin} > $outdir/spllt_starpu/${matname}_NCPU-${ncpu}_NB-${nb}_NEMIN-${nemin}${outsuffix}
                    #     ;;
                    # starpu_gpu)
                    #     ./builds/starpu-gpu/spllt_test --ncpu ${ncpu} --nb ${nb} --nemin ${nemin} > $outdir/spllt_starpu/${matname}_NCPU-${ncpu}_NB-${nb}_NEMIN-${nemin}${outsuffix}
                    #     ;;

                    gnu_omp)
                        ../builds/omp/gnu/spllt_test --mm ${matrix} --ncpu ${ncpu} --nb ${nb} --nemin ${nemin} > $outdir/spllt_omp/gnu/${matname}_NCPU-${ncpu}_NB-${nb}_NEMIN-${nemin}${outsuffix}
                        ;;

                    gnu_omp_prune)
                        ../builds/omp/gnu/spllt_test --prune-tree --mm ${matrix} --ncpu ${ncpu} --nb ${nb} --nemin ${nemin} > $outdir/spllt_omp/gnu/${matname}_NCPU-${ncpu}_NB-${nb}_NEMIN-${nemin}${outsuffix}
                        ;;
                    # intel_omp)
                    #     ./builds/omp//spllt_test --ncpu ${ncpu} --nb ${nb} --nemin ${nemin} > $outdir/spllt_omp/intel/${matname}_NCPU-${ncpu}_NB-${nb}_NEMIN-${nemin}${outsuffix}
                    #     ;;
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
