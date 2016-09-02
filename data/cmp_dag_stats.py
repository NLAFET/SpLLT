#!/usr/bin/env python

import sys
import os
import fileinput
import re
import subprocess
import csv
import StringIO

import vextractor
import latexprint as lp

spllt_facto_time = '\[>\] \[factorize\] time:'
ma87_facto_time  = 'Factor took'
spllt_task_insert_time = '\[>\] \[spllt_stf_factorize\] task insert time:'

blocksizes = [256, 384, 512, 768, 1024]

# print '% data directory: ', sys.argv[1]
outputdir = sys.argv[1]
# listmat  = outputdir + '/list.matrix'  
listmat = sys.argv[2]
flistmat = open(listmat)

# StarPU scheduler
starpu_sched = 'lws'

matcount = 1

for mat in flistmat:

    mat = mat.rstrip()

    # print mat
    pbl = re.sub(r'/', '_', mat)
    pbl = pbl.rstrip()

    # StarPU
    spllt_t_facto = []
    spllt_t_insert = []
    for blocksize in blocksizes:
        # print "blocksize: ",blocksize
        datafile = outputdir + '/' + 'spllt_starpu' + '/' + starpu_sched + '/' + pbl + '_NCPU-27' + '_NB-' + str(blocksize) + '_NEMIN-32'
        # print datafile
        f = open(datafile)
        v = vextractor.get_value(f, spllt_facto_time)
        f.seek(0)
        vi = vextractor.get_value(f, spllt_task_insert_time)
        # print vi
        f.close()

        spllt_t_facto.append(float(v))
        spllt_t_insert.append(float(vi))


    best_spllt_t_facto_idx = spllt_t_facto.index(min(spllt_t_facto))
    
    # StarPU
    best_spllt_nb       = blocksizes[best_spllt_t_facto_idx]
    best_spllt_t_facto  = spllt_t_facto[best_spllt_t_facto_idx]
    best_spllt_t_insert = spllt_t_insert[best_spllt_t_facto_idx]

    # print("%d" % (best_spllt_nb))

    tracefile = outputdir + '/' + 'spllt_starpu' + '/' + starpu_sched + '/traces/' + pbl + '_NCPU-27' + '_NB-' + str(best_spllt_nb) + '_NEMIN-32.rec'
    # print("tracesfile: %s" % (tracefile))
    cmd = "starpu_trace_state_stats.py -t %s" % (tracefile)
    # subprocess.call(cmd.split())
    output = subprocess.check_output(cmd.split()) 
    # process = subprocess.Popen(cmd.split())
    # output = process.communicate()[0]
    # print(output)
    outputIO = StringIO.StringIO(output)

    ntask = 0
    ttask = 0.0
    # timefile = open('dag_stats.csv', 'r')
    reader = csv.reader(outputIO)
    for row in reader:
        if row[2] == 'Task':
            ntask = ntask + int(row[1])
            ttask = ttask + float(row[3])

    # print("Number of tasks: %d, Time spent in tasks: %f, avg time per task: %f" % (ntask, ttask, ttask/ntask))
    print("%d %f" % (ntask, ttask/ntask))

    # try:
    # finally:
        # timefile.close()
