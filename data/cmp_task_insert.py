#!/usr/bin/env python

import sys
import os
import fileinput
import re

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
        
    # OpenMP (gnu)
    spllt_gnu_omp_t_facto = []
    spllt_gnu_omp_t_insert = []
    for blocksize in blocksizes:
        datafile = outputdir + '/' + 'spllt_omp' + '/' + 'gnu' + '/' + pbl + '_NCPU-28' + '_NB-' + str(blocksize) + '_NEMIN-32'
        if os.path.exists(datafile):
            # print datafile
            f = open(datafile)
            v = vextractor.get_value(f, spllt_facto_time)
            f.seek(0)
            vi = vextractor.get_value(f, spllt_task_insert_time)
            # print v
            f.close()
            spllt_gnu_omp_t_facto.append(float(v))
            spllt_gnu_omp_t_insert.append(float(vi))

    # MA87
    # ma87_t_facto = []
    # for blocksize in blocksizes:
    #     # print "blocksize: ",blocksize
    #     datafile = outputdir + '/' + 'ma87' + '/' + pbl + '_NCPU-28' + '_NB-' + str(blocksize) + '_NEMIN-32'
    #     # print datafile
    #     f = open(datafile)
    #     v = vextractor.get_value(f, ma87_facto_time)
    #     f.close()

    #     ma87_t_facto.append(float(v))
    
    # print spllt_t_facto
    # print spllt_t_insert
    # print ma87_t_facto
        
    best_spllt_gnu_omp_t_facto_idx = spllt_gnu_omp_t_facto.index(min(spllt_gnu_omp_t_facto))
    best_spllt_t_facto_idx = spllt_t_facto.index(min(spllt_t_facto))
    # best_ma87_t_facto_idx  = ma87_t_facto.index(min(ma87_t_facto))
    
    # StarPU
    best_spllt_nb       = blocksizes[best_spllt_t_facto_idx]
    best_spllt_t_facto  = spllt_t_facto[best_spllt_t_facto_idx]
    best_spllt_t_insert = spllt_t_insert[best_spllt_t_facto_idx]

    # OpenMP
    best_spllt_gnu_omp_nb       = blocksizes[best_spllt_gnu_omp_t_facto_idx]
    best_spllt_gnu_omp_t_facto  = spllt_gnu_omp_t_facto[best_spllt_gnu_omp_t_facto_idx]
    best_spllt_gnu_omp_t_insert = spllt_gnu_omp_t_insert[best_spllt_gnu_omp_t_facto_idx]

    # MA87
    # best_ma87_nb       = blocksizes[best_ma87_t_facto_idx]
    # best_ma87_t_facto  = ma87_t_facto[best_ma87_t_facto_idx]
    
    print("%4s & %40s & %10s & %10s \\\\" % (matcount, lp.escape(mat),
                                             lp.print_float(best_spllt_gnu_omp_t_insert),
                                             lp.print_float(best_spllt_t_insert)))

    # print("%40s & %10s & %10s & %10s \\\\" % (lp.escape(mat), 
    #                                           lp.print_float(best_spllt_t_insert),
    #                                           lp.print_float(best_spllt_t_facto, (best_spllt_t_facto<best_ma87_t_facto)), 
    #                                           lp.print_float(best_ma87_t_facto , (best_ma87_t_facto<best_spllt_t_facto))))
    
    # print("%40s & %10s \\\\" % (lp.escape(mat), lp.print_float(best_spllt_t_facto, (best_spllt_t_facto<best_ma87_t_facto)))

    matcount = matcount+1 
