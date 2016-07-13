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
            ./run_ma87 --ncpu ${ncpu} --nb ${nb} --nemin ${nemin} > $outdir/ma87/${matname}_NCPU-${ncpu}_NB-${nb}_NEMIN-${nemin}${outsuffix}
            ;;
        stf)
            ./spllt_stf_test --ncpu ${ncpu} --nb ${nb} --nemin ${nemin} > $outdir/spllt_stf/${matname}_NB-${nb}_NEMIN-${nemin}${outsuffix}
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
        # mv $trace_dir/$prof_file $outdir/spllt_starpu/traces/${matname}_NCPU-${ncpu}_NB-${nb}.prof
        # starpu_fxt_tool -c -i $outdir/spllt_starpu/traces/${matname}_NCPU-${ncpu}_NB-${nb}.prof -o $outdir/spllt_starpu/traces/${matname}_NCPU-${ncpu}_NB-${nb}.trace
        # starpu_fxt_tool -c -i $outdir/spllt_starpu/traces/${matname}_NCPU-${ncpu}_NB-${nb}.prof
        starpu_fxt_tool -c -i $trace_dir/$prof_file
        mv paje.trace $outdir/spllt_starpu/traces/${matname}_NCPU-${ncpu}_NB-${nb}.trace
        mv tasks.rec $outdir/spllt_starpu/traces/${matname}_NCPU-${ncpu}_NB-${nb}_tasks.rec
        mv trace.rec $outdir/spllt_starpu/traces/${matname}_NCPU-${ncpu}_NB-${nb}.rec
    fi

done
