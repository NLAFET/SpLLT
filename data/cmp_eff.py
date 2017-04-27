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

# spllt_facto_time = '\[>\] \[factorize\] time:'
spllt_facto_time = 'Factor took'
ma87_facto_time  = 'Factor took'
spllt_task_insert_time = '\[>\] \[spllt_stf_factorize\] task insert time:'
spllt_flops_str = '\[>\] \[analysis\] num flops :'

blocksizes = [256, 384, 512, 768, 1024]
# ncpu = 24

# print '% data directory: ', sys.argv[1]
outputdir = sys.argv[1]
# listmat  = outputdir + '/list.matrix'  
listmat = sys.argv[2]
flistmat = open(listmat)

matcount = 1

# LWS scheduler without tree pruning 
starpu_sched = 'lws'
ncpus = 27
# data print (factorization times)
# print("%4s %40s %10s %10s %10s" % ('#', 'Name', 'MA87', 'OMP', 'StarPU'))

# data print (factorization times and block sizes)
# print("%4s %40s %17s %17s %17s" % ('#', 'Name', 'MA87', 'OMP', 'StarPU'))
# print("%4s %40s %6s %10s %6s %10s %6s %10s" % ( ' ', ' ', 'nb', 'time', 'nb', 'time', 'nb', 'time'))

for mat in flistmat:

    mat = mat.rstrip()

    # print mat
    pbl = re.sub(r'/', '_', mat)
    pbl = pbl.rstrip()

    v = 0.0

    # StarPU
    spllt_t_facto = []
    spllt_t_insert = []
    spllt_flops = []
    for blocksize in blocksizes:
        # print "blocksize: ",blocksize
        datafile = outputdir + '/' + 'spllt_starpu' + '/' + starpu_sched + '/' + pbl + '_NCPU-27' + '_NB-' + str(blocksize) + '_NEMIN-32'
        # Create data structure
        df = Datafile(datafile)

        # Get factor time
        v = df.get_value(spllt_facto_time)
        spllt_t_facto.append(float(v))
        
        # Get task insert time
        vi = df.get_value(spllt_task_insert_time)
        spllt_t_insert.append(float(vi))

        # Get flops
        vf = df.get_value(spllt_flops_str)
        spllt_flops.append(float(vf))


    # STF
    t_stf = []
    for blocksize in blocksizes:
        # print "blocksize: ",blocksize
        datafile = outputdir + '/' + 'spllt_stf' + '/' + pbl + '_NB-' + str(blocksize) + '_NEMIN-32'

        df = Datafile(datafile)
        
        # Get factor time
        v = df.get_value(spllt_facto_time)
        t_stf.append(float(v))

    # print spllt_t_facto
    # print spllt_t_insert
    # print ma87_t_facto
        
    best_spllt_t_facto_idx = spllt_t_facto.index(min(spllt_t_facto))

    # StarPU
    best_spllt_nb       = blocksizes[best_spllt_t_facto_idx]
    best_spllt_t_facto  = spllt_t_facto[best_spllt_t_facto_idx]
    best_spllt_t_insert = spllt_t_insert[best_spllt_t_facto_idx]
    best_spllt_flops    = spllt_flops[best_spllt_t_facto_idx]
    # GFlops 
    best_spllt_flops    = best_spllt_flops/1e9
    # print best_spllt_t_facto

    # STF
    # best_t_stf_idx = t_stf.index(min(t_stf))
    best_t_stf = min(t_stf)
    # print best_t_stf

    tracefile = outputdir + '/' + 'spllt_starpu' + '/' + starpu_sched + '/traces/' + pbl + '_NCPU-27' + '_NB-' + str(best_spllt_nb) + '_NEMIN-32.rec'
    cmd = "starpu_trace_state_stats.py -te -s %f %s" % (best_t_stf*1000, tracefile)
    # print cmd 
    # subprocess.call(cmd.split())
    output = subprocess.check_output(cmd.split()) 
    # process = subprocess.Popen(cmd.split())
    # output = process.communicate()[0]
    # print(output)
    # outputIO = StringIO.StringIO(output)
    
    # data print (GFlop/s)
    # print("%4s %10.3f %10.3f %10.3f" % (matcount, 
    #                                     (best_spllt_flops/best_ma87_t_facto), 
    #                                     (best_spllt_flops/best_spllt_gnu_omp_t_facto), 
    #                                     (best_spllt_flops/best_spllt_t_facto)))

    # e = 0.0
    # et = 0.0
    # er = 0.0
    # es = 0.0
    # ep = 0.0
    # with open('efficiencies.csv', 'rb') as csvfile:
    #     reader = csv.reader(csvfile)

    #     for row in reader:
    #         # print row
    #         if row[0] == 'Parallel':
    #             e = float(row[1])
    #         elif row[0] == 'Task':
    #             et = float(row[1])
    #         elif row[0] == 'Runtime':
    #             er = float(row[1])
    #         elif row[0] == 'Scheduling':
    #             es = float(row[1])
    #         elif row[0] == 'Pipeline':
    #             ep = float(row[1])
            

    # eo = er*es

    # print("%d %f %f %f %f" % (matcount, e, et, eo, ep))

    tt = 0.0
    tr = 0.0
    ti = 0.0
    ts = 0.0

    with open('times.csv', 'rb') as csvfile:
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
    tr_real = best_spllt_t_facto * ncpus * pr
    tt_real = best_spllt_t_facto * ncpus * pt
    ti_real = best_spllt_t_facto * ncpus * pi
    ts_real = best_spllt_t_facto * ncpus * ps
    # compute time spent in runtime + scheduling
    to_real = tr_real + ts_real

    et = best_t_stf / tt_real
    eo = tt_real / (tt_real+to_real)
    ep = (tt_real+to_real) / (tt_real+to_real+ti_real)
    # compute efficiency
    e = et * eo * ep

    # print("%f %f" % (best_t_stf, best_spllt_t_facto))
    print("%d %f %f %f %f" % (matcount, e, et, eo, ep))    

    matcount = matcount+1 
