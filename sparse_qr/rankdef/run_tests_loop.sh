for matrix in `cat list.matrix`
do

    echo "[test] matrix: $matrix"

    for sigma in ${sigma_list[@]}
    do

        echo "[test] sigma: $sigma"

        ./run_rd ${matrix} ${sigma} >> $outdir/${matname}_SIGMA-1e${sigma}  
    done
done
