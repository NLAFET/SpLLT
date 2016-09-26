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
spllt_flops_str = '\[>\] \[analysis\] num flops :'

blocksizes = [256, 384, 512, 768, 1024]
# ncpu = 24

# print '% data directory: ', sys.argv[1]
outputdir = sys.argv[1]
# listmat  = outputdir + '/list.matrix'  
listmat = sys.argv[2]
flistmat = open(listmat)

matcount = 1

starpu_dir = 'lws'
starpu_prune_dir = 'lws_prune'

for mat in flistmat:

    mat = mat.rstrip()

    # print mat
    pbl = re.sub(r'/', '_', mat)
    pbl = pbl.rstrip()

    v = 0.0

    # MA87
    ma87_t_facto = []
    for blocksize in blocksizes:
        # print "blocksize: ",blocksize
        datafile = outputdir + '/' + 'ma87' + '/' + pbl + '_NCPU-28' + '_NB-' + str(blocksize) + '_NEMIN-32'
        if os.path.exists(datafile):

            # print datafile
            f = open(datafile)
            v = vextractor.get_value(f, ma87_facto_time)
            f.close()

            ma87_t_facto.append(float(v))

    # StarPU
    spllt_t_facto = []
    spllt_t_insert = []
    spllt_flops = []
    for blocksize in blocksizes:
        # print "blocksize: ",blocksize
        datafile = outputdir + '/' + 'spllt_starpu' + '/' + starpu_dir + '/' + pbl + '_NCPU-27' + '_NB-' + str(blocksize) + '_NEMIN-32'
        if os.path.exists(datafile):
            # print datafile
            f = open(datafile)
            v = vextractor.get_value(f, spllt_facto_time)
            spllt_t_facto.append(float(v))
            # get tasks insert time
            f.seek(0)
            vi = vextractor.get_value(f, spllt_task_insert_time)
            spllt_t_insert.append(float(vi))
            # get flops
            f.seek(0)
            vf = vextractor.get_value(f, spllt_flops_str)
            # print vf
            spllt_flops.append(float(vf))
            
            f.close()

    # StarPU with tree pruning
    spllt_prune_t_facto = []
    spllt_prune_t_insert = []
    spllt_prune_flops = []
    for blocksize in blocksizes:
        # print "blocksize: ",blocksize
        datafile = outputdir + '/' + 'spllt_starpu' + '/' + starpu_prune_dir + '/' + pbl + '_NCPU-27' + '_NB-' + str(blocksize) + '_NEMIN-32'
        if os.path.exists(datafile):
            # print datafile
            f = open(datafile)
            v = vextractor.get_value(f, spllt_facto_time)
            spllt_prune_t_facto.append(float(v))
            # get tasks insert time
            f.seek(0)
            vi = vextractor.get_value(f, spllt_task_insert_time)
            spllt_prune_t_insert.append(float(vi))
            # get flops
            f.seek(0)
            vf = vextractor.get_value(f, spllt_flops_str)
            # print vf
            spllt_prune_flops.append(float(vf))
            
            f.close()

    # OpenMP (gnu)
    spllt_gnu_omp_t_facto = []
    for blocksize in blocksizes:
        datafile = outputdir + '/' + 'spllt_omp' + '/' + 'gnu' + '/' + pbl + '_NCPU-28' + '_NB-' + str(blocksize) + '_NEMIN-32'
        if os.path.exists(datafile):
            # print datafile
            f = open(datafile)
            v = vextractor.get_value(f, spllt_facto_time)
            # print v
            f.close()
            spllt_gnu_omp_t_facto.append(float(v))

    # OpenMP (gnu) with tree pruning
    spllt_gnu_omp_prune_t_facto = []
    for blocksize in blocksizes:
        datafile = outputdir + '/' + 'spllt_omp' + '/' + 'gnu_prune' + '/' + pbl + '_NCPU-28' + '_NB-' + str(blocksize) + '_NEMIN-32'
        if os.path.exists(datafile):
            # print datafile
            f = open(datafile)
            v = vextractor.get_value(f, spllt_facto_time)
            # print v
            f.close()
            spllt_gnu_omp_prune_t_facto.append(float(v))


    best_spllt_t_facto_idx = spllt_t_facto.index(min(spllt_t_facto))
    best_spllt_prune_t_facto_idx  = spllt_prune_t_facto.index(min(spllt_prune_t_facto))
    best_ma87_t_facto_idx  = ma87_t_facto.index(min(ma87_t_facto))
    best_spllt_gnu_omp_t_facto_idx = spllt_gnu_omp_t_facto.index(min(spllt_gnu_omp_t_facto))
    best_spllt_gnu_omp_prune_t_facto_idx = spllt_gnu_omp_prune_t_facto.index(min(spllt_gnu_omp_prune_t_facto))

    # MA87
    best_ma87_nb       = blocksizes[best_ma87_t_facto_idx]
    best_ma87_t_facto  = ma87_t_facto[best_ma87_t_facto_idx]

    # StarPU
    best_spllt_nb       = blocksizes[best_spllt_t_facto_idx]
    best_spllt_t_facto  = spllt_t_facto[best_spllt_t_facto_idx]
    best_spllt_t_insert = spllt_t_insert[best_spllt_t_facto_idx]
    best_spllt_flops    = spllt_flops[best_spllt_t_facto_idx]
    # GFlops 
    best_spllt_flops    = best_spllt_flops/(1024*1024*1024)

    # StarPU with tree pruning
    best_spllt_prune_nb       = blocksizes[best_spllt_prune_t_facto_idx]
    best_spllt_prune_t_facto  = spllt_prune_t_facto[best_spllt_prune_t_facto_idx]
    best_spllt_prune_t_insert = spllt_prune_t_insert[best_spllt_prune_t_facto_idx]
    best_spllt_prune_flops    = spllt_prune_flops[best_spllt_prune_t_facto_idx]

    # OMP (gnu)
    best_spllt_gnu_omp_nb       = blocksizes[best_spllt_gnu_omp_t_facto_idx]
    best_spllt_gnu_omp_t_facto  = spllt_gnu_omp_t_facto[best_spllt_gnu_omp_t_facto_idx]

    # OMP (gnu) with tree pruning
    best_spllt_gnu_omp_prune_nb       = blocksizes[best_spllt_gnu_omp_prune_t_facto_idx]
    best_spllt_gnu_omp_prune_t_facto  = spllt_gnu_omp_prune_t_facto[best_spllt_gnu_omp_prune_t_facto_idx]

    # print("%10.3f" % best_spllt_flops)

    # data print (GFlop/s) with Parsec
    print("%4s %10.3f %10.3f %10.3f %10.3f %10.3f" % (matcount,
                                                      (best_spllt_flops/best_ma87_t_facto), 
                                                      (best_spllt_flops/best_spllt_gnu_omp_t_facto),
                                                      (best_spllt_flops/best_spllt_gnu_omp_prune_t_facto),
                                                      (best_spllt_flops/best_spllt_t_facto),
                                                      (best_spllt_flops/best_spllt_prune_t_facto)))

    matcount = matcount+1 
