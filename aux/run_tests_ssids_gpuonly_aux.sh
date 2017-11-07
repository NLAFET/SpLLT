#!/bin/bash

for matrix in `cat list.matrix`
do    
    matname=`echo $matrix | sed -e 's/\//_/g'`

    ./prep.sh $matrix

    ~/builds/spral-gpuonly/spral_ssids --pos > $outdir/ssids-gpuonly/${matname}
done
