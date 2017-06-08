#!/bin/bash

for matrix in `cat list.matrix`
do

    ./prep.sh $matrix

    ~/builds/spral-gpuonly/spral_ssids --pos
done
