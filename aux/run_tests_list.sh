for mat in "${matrices[@]}"
do

    matrix=$(echo $mat | awk -F'matrix:' '{print $2}' | cut -f1 -d' ')
    ncpu=$(echo $mat | awk -F'ncpu:' '{print $2}' | cut -f1 -d' ')
    nb=$(echo $mat | awk -F'nb:' '{print $2}' | cut -f1 -d' ')
    nemin=$(echo $mat | awk -F'nemin:' '{print $2}' | cut -f1 -d' ')

    matname=`echo $matrix | sed -e 's/\//_/g'`
    echo "[run_tests] test matrix: $matname"

    echo "[run_tests] extract matrix: $matname"
    ./prep.sh $matrix

    case $build in
        ma87)
            echo "[run_tests] run MA87"
            export OMP_NUM_THREADS=${ncpu} 
            ./builds/ma87/spllt_test --ncpu ${ncpu} --nb ${nb} --nemin ${nemin} > $outdir/ma87/${matname}_NCPU-${ncpu}_NB-${nb}_NEMIN-${nemin}${outsuffix}
            ;;
        stf)
            ./builds/stf/spllt_test --ncpu ${ncpu} --nb ${nb} --nemin ${nemin} > $outdir/spllt_stf/${matname}_NB-${nb}_NEMIN-${nemin}${outsuffix}
            ;;
        parsec)
            ./builds/parsec/spllt_test --ncpu ${ncpu} --nb ${nb} --nemin ${nemin} > $outdir/spllt_parsec/${matname}_NCPU-${ncpu}_NB-${nb}_NEMIN-${nemin}${outsuffix}
            ;;
        starpu)
            ./builds/starpu/spllt_test --ncpu ${ncpu} --nb ${nb} --nemin ${nemin} > $outdir/spllt_starpu/${matname}_NCPU-${ncpu}_NB-${nb}_NEMIN-${nemin}${outsuffix}
            ;;
        starpu_prune)
            ./builds/starpu/spllt_test --ncpu ${ncpu} --nb ${nb} --nemin ${nemin} --prune-tree > $outdir/spllt_starpu/${matname}_NCPU-${ncpu}_NB-${nb}_NEMIN-${nemin}${outsuffix}
            ;;
        starpu_nested_stf)
            ./builds/starpu-nested-stf/spllt_test --ncpu ${ncpu} --nb ${nb} --nemin ${nemin} > $outdir/spllt_starpu_nested_stf/${matname}_NCPU-${ncpu}_NB-${nb}_NEMIN-${nemin}${outsuffix}
            ;;
        starpu_trace)
            ./builds/starpu-trace/spllt_test --ncpu ${ncpu} --nb ${nb} --nemin ${nemin} > $outdir/spllt_starpu/${matname}_NCPU-${ncpu}_NB-${nb}_NEMIN-${nemin}${outsuffix}
            ;;
        gnu_omp)
            ./builds/omp/spllt_test --ncpu ${ncpu} --nb ${nb} --nemin ${nemin} > $outdir/spllt_omp/gnu/${matname}_NCPU-${ncpu}_NB-${nb}_NEMIN-${nemin}${outsuffix}
            ;;
        intel_omp)
            ./builds/omp/spllt_test --ncpu ${ncpu} --nb ${nb} --nemin ${nemin} > $outdir/spllt_omp/intel/${matname}_NCPU-${ncpu}_NB-${nb}_NEMIN-${nemin}${outsuffix}
            ;;
    esac

    if [ -f $trace_dir/$prof_file ];
    then
        # mv $trace_dir/$prof_file $outdir/spllt_starpu/traces/${matname}_NCPU-${ncpu}_NB-${nb}.prof
        # starpu_fxt_tool -c -i $outdir/spllt_starpu/traces/${matname}_NCPU-${ncpu}_NB-${nb}.prof -o $outdir/spllt_starpu/traces/${matname}_NCPU-${ncpu}_NB-${nb}.trace
        # starpu_fxt_tool -c -i $outdir/spllt_starpu/traces/${matname}_NCPU-${ncpu}_NB-${nb}.prof
        starpu_fxt_tool -c -i $trace_dir/$prof_file
        mv paje.trace $outdir/spllt_starpu/traces/${matname}_NCPU-${ncpu}_NB-${nb}_NEMIN-${nemin}.trace
        mv tasks.rec $outdir/spllt_starpu/traces/${matname}_NCPU-${ncpu}_NB-${nb}_NEMIN-${nemin}_tasks.rec
        mv trace.rec $outdir/spllt_starpu/traces/${matname}_NCPU-${ncpu}_NB-${nb}_NEMIN-${nemin}.rec
    fi

done
