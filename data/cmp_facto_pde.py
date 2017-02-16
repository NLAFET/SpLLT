#!/usr/bin/env python

import sys
import os
import fileinput
import re
import csv
import subprocess
import StringIO

import vextractor
import latexprint as lp
from datafile import Datafile

# List of blocksizes
blocksizes = [256, 384, 512, 768, 1024]

# List of mesh sizes
msh_sizes = [32, 64, 96, 128, 160]

# Directory containing the outputs
outputdir = sys.argv[1]

# StarPU scheduler 
starpu_sched = 'ws'

# Parsing parameters
spllt_facto_time = 'Factor took' # '\[>\] \[factorize\] time:'
ma87_facto_time  = 'Factor took'
spllt_flops_str =  'Predict nflop =' # '\[>\] \[analysis\] num flops :'

for msh_size in msh_sizes:

    pbl = 'poisson3d_' + str(msh_size) + 'x' + str(msh_size) + 'x'  + str(msh_size)
    # print(pbl)

    # Flop count
    flops = []

    # StarPU
    starpu_t = []
    starpu_t_insert = []
    for blocksize in blocksizes:
        # print "blocksize: ",blocksize
        datafile = outputdir + '/' + 'spllt_starpu' + '/' + starpu_sched + '/' + pbl + '_NCPU-27' + '_NB-' + str(blocksize) + '_NEMIN-32'
        # print(datafile)
        # Create data structure
        df = Datafile(datafile)
        # Get factor time
        v = df.get_value(spllt_facto_time)
        starpu_t.append(float(v))
        
        # # Get task insert time
        # vi = df.get_value(spllt_task_insert_time)
        # spllt_t_insert.append(float(vi))

        # Get flops
        vf = df.get_value(spllt_flops_str)
        flops.append(float(vf))

    # StarPU with pruning
    starpu_prune_t = []
    starpu_prune_insert = []
    for blocksize in blocksizes:
        # print "blocksize: ",blocksize
        datafile = outputdir + '/' + 'spllt_starpu' + '/' + starpu_sched + '_prune/' + pbl + '_NCPU-27' + '_NB-' + str(blocksize) + '_NEMIN-32'
        # print(datafile)
        # Create data structure
        df = Datafile(datafile)
        # Get factor time
        v = df.get_value(spllt_facto_time)
        starpu_prune_t.append(float(v))
        
        # # Get task insert time
        # vi = df.get_value(spllt_task_insert_time)
        # spllt_t_insert.append(float(vi))

        # Get flops
        # vf = df.get_value(spllt_flops_str)
        # flops.append(float(vf))


    # OpenMP (gnu)
    gnu_omp_t = []
    for blocksize in blocksizes:
        datafile = outputdir + '/' + 'spllt_omp' + '/' + 'gnu' + '/' + pbl + '_NCPU-27' + '_NB-' + str(blocksize) + '_NEMIN-32'
        if os.path.exists(datafile):
            # print datafile
            f = open(datafile)
            v = vextractor.get_value(f, spllt_facto_time)
            # print v
            f.close()
            gnu_omp_t.append(float(v))

    # OpenMP (gnu) with pruning
    gnu_omp_prune_t = []
    for blocksize in blocksizes:
        datafile = outputdir + '/' + 'spllt_omp' + '/' + 'gnu_prune' + '/' + pbl + '_NCPU-27' + '_NB-' + str(blocksize) + '_NEMIN-32'
        if os.path.exists(datafile):
            # print datafile
            f = open(datafile)
            v = vextractor.get_value(f, spllt_facto_time)
            # print v
            f.close()
            gnu_omp_prune_t.append(float(v))

    # MA87
    ma87_t = []
    for blocksize in blocksizes:
        # print "blocksize: ",blocksize
        datafile = outputdir + '/' + 'ma87' + '/' + pbl + '_NCPU-28' + '_NB-' + str(blocksize) + '_NEMIN-32'
        if os.path.exists(datafile):

            # print datafile
            f = open(datafile)
            v = vextractor.get_value(f, ma87_facto_time)
            f.close()

            ma87_t.append(float(v))

    # parsec_t_idx = spllt_parsec_t_facto.index(min(spllt_parsec_t_facto))
    gnu_omp_idx = gnu_omp_t.index(min(gnu_omp_t))
    gnu_omp_prune_idx = gnu_omp_prune_t.index(min(gnu_omp_prune_t))
    starpu_idx = starpu_t.index(min(starpu_t))
    starpu_prune_idx = starpu_prune_t.index(min(starpu_prune_t))
    ma87_idx  = ma87_t.index(min(ma87_t))

    # OMP
    gnu_omp_nb = blocksizes[gnu_omp_idx]
    best_gnu_omp_t  = gnu_omp_t[gnu_omp_idx]

    # OMP with pruning
    gnu_omp_prune_nb = blocksizes[gnu_omp_prune_idx]
    best_gnu_omp_prune_t  = gnu_omp_prune_t[gnu_omp_prune_idx]

    # StarPU
    starpu_nb = blocksizes[starpu_idx]
    best_starpu_t  = starpu_t[starpu_idx]
    # best_spllt_t_insert = spllt_t_insert[best_spllt_t_facto_idx]
    best_flops    = flops[starpu_idx]
    # # GFlops 
    best_gflops    = best_flops/1e9

    # StarPU with pruning
    starpu_prune_nb = blocksizes[starpu_prune_idx]
    best_starpu_prune_t  = starpu_prune_t[starpu_prune_idx]
    # best_starpu_prune_insert = starpu_prune_insert[starpu_prune_idx]

    # MA87
    ma87_nb = blocksizes[ma87_idx]
    best_ma87_t  = ma87_t[ma87_idx]

    # # data print (GFlop/s)
    # print("%4s %10.3f %10.3f %10.3f" % (msh_size,
    #                                     (best_gflops/best_ma87_t), 
    #                                     (best_gflops/best_gnu_omp_t), 
    #                                     (best_gflops/best_starpu_t)))

    # # data print (GFlop/s) with prunning
    # print("%4s %10.3f %10.3f %10.3f %10.3f %10.3f" % (msh_size,
    #                                                   (best_gflops/best_ma87_t), 
    #                                                   (best_gflops/best_gnu_omp_t), 
    #                                                   (best_gflops/best_starpu_t),
    #                                                   (best_gflops/best_gnu_omp_prune_t),
    #                                                   (best_gflops/best_starpu_prune_t)))

    # data print (Timing) with prunning
    print("%4s %10.3f %4d %10.3f %4d %10.3f %4d %10.3f %4d %10.3f %4d" % (msh_size,
                                                                          best_ma87_t, ma87_nb,
                                                                          best_gnu_omp_t, gnu_omp_prune_nb,
                                                                          best_starpu_t, starpu_nb,
                                                                          best_gnu_omp_prune_t, gnu_omp_prune_nb,
                                                                          best_starpu_prune_t, starpu_prune_nb))

