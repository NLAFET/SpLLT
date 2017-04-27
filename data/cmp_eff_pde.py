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
# msh_sizes = [96]

# Directory containing the outputs
outputdir = sys.argv[1]

matcount = 1

# StarPU scheduler 
starpu_sched = 'ws'
ncpus = 27

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
        datafile = outputdir + '/' + 'starpu' + '/' + starpu_sched + '/' + pbl + '_NCPU-27' + '_NB-' + str(blocksize) + '_NEMIN-32'
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
        datafile = outputdir + '/' + 'starpu' + '/' + starpu_sched + '_prune/' + pbl + '_NCPU-27' + '_NB-' + str(blocksize) + '_NEMIN-32'
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

    # STF
    t_stf = []
    for blocksize in blocksizes:
        # print "blocksize: ",blocksize
        datafile = outputdir + '/' + 'stf' + '/' + pbl + '_NB-' + str(blocksize) + '_NEMIN-32'

        df = Datafile(datafile)
        
        # Get factor time
        v = df.get_value(spllt_facto_time)
        t_stf.append(float(v))


    starpu_idx = starpu_t.index(min(starpu_t))

    # StarPU
    best_starpu_nb  = blocksizes[starpu_idx]
    best_starpu_t   = starpu_t[starpu_idx]
    # best_spllt_t_insert = spllt_t_insert[best_spllt_t_facto_idx]
    # best_spllt_flops    = spllt_flops[best_spllt_t_facto_idx]
    # # GFlops 
    # best_spllt_flops    = best_spllt_flops/1e9
    # # print best_spllt_t_facto

    # STF
    # best_t_stf_idx = t_stf.index(min(t_stf))
    best_t_stf = min(t_stf)
    # print best_t_stf

    cvsfile = outputdir + '/' + 'starpu' + '/' + starpu_sched + '/traces/' + pbl + '_NCPU-27' + '_NB-' + str(best_starpu_nb) + '_NEMIN-32.csv'
    
    # If the CSV file dosn't exists, we use the StarPU script for
    # retreiving the information from the traces.
    if not os.path.exists(cvsfile):

        tracefile = outputdir + '/' + 'starpu' + '/' + starpu_sched + '/traces/' + pbl + '_NCPU-27' + '_NB-' + str(best_starpu_nb) + '_NEMIN-32.rec'
        cmd = "starpu_trace_state_stats.py -te -s %f %s" % (best_t_stf*1000, tracefile)
    
        output = subprocess.check_output(cmd.split()) 

        cvsfile = 'times.csv'

    tt = 0.0
    tr = 0.0
    ti = 0.0
    ts = 0.0

    with open(cvsfile, 'rb') as csvfile:
        reader = csv.reader(csvfile)

        for row in reader:
            if row[0] == 'Runtime':
                tr = float(row[1])
            elif row[0] == 'Task':
                tt = float(row[1])
            elif row[0] == 'Idle':
                ti = float(row[1])
            elif row[0] == 'Scheduling':
                ts = float(row[1])


    t_tot = tr + tt + ti + ts

    pr = tr/t_tot
    pt = tt/t_tot
    pi = ti/t_tot
    ps = ts/t_tot

    # compute actual cumulative times
    tr_real = best_starpu_t * ncpus * pr
    tt_real = best_starpu_t * ncpus * pt
    ti_real = best_starpu_t * ncpus * pi
    ts_real = best_starpu_t * ncpus * ps
    # compute time spent in runtime + scheduling
    to_real = tr_real + ts_real

    # print("tt_real: %f, best_t_stf: %f" % (tt_real, best_t_stf))

    et = best_t_stf / tt_real
    eo = tt_real / (tt_real+to_real)
    ep = (tt_real+to_real) / (tt_real+to_real+ti_real)
    # compute efficiency
    e = et * eo * ep

    # print("%f %f" % (best_t_stf, best_spllt_t_facto))
    print("%d %f %f %f %f" % (msh_size, e, et, eo, ep))

    matcount = matcount+1 

