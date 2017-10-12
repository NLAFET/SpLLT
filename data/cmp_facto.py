#!/usr/bin/env python

import sys
import os
import fileinput
import re

import vextractor
import latexprint as lp
from datafile import Datafile

# spllt_facto_time = '\[>\] \[factorize\] time:'
t_facto_str = 'Factor took'
ma87_facto_time  = 'Factor took'
spllt_task_insert_time = '\[>\] \[spllt_stf_factorize\] task insert time:'
flops_str =  'Predict nflop =' # '\[>\] \[analysis\] num flops :'

# blocksizes = [256, 384, 512, 768, 1024]
blocksizes = [256, 384, 512, 768, 1024, 1536]
# ncpu = 24

# print '% data directory: ', sys.argv[1]
outputdir = sys.argv[1]
listmat = sys.argv[2]
flistmat = open(listmat)

starpu_sched = 'ws'

matcount = 1

for mat in flistmat:

    mat = mat.rstrip()

    # print mat
    pbl = re.sub(r'/', '_', mat)
    pbl = pbl.rstrip()

    # Flop count
    flops = []

    v = 0.0

    # SpLLT StarPU
    starpu_t_facto = []
    for blocksize in blocksizes:
        # print "blocksize: ",blocksize
        datafile = outputdir + '/' + 'starpu' + '/' + str(starpu_sched) + '/' + pbl + '_NCPU-27' + '_NB-' + str(blocksize) + '_NEMIN-32'
        # print "datafile: ", datafile

        df = Datafile(datafile)
        # Get factor time
        v = df.get_value(t_facto_str)
        starpu_t_facto.append(float(v))
        # Get flops
        vf = df.get_value(flops_str)
        flops.append(float(vf))

    # SpLLT Parsec
    parsec_t_facto = []
    for blocksize in blocksizes:
        # print "blocksize: ",blocksize
        datafile = outputdir + '/' + 'parsec' + '/' + pbl + '_NCPU-28' + '_NB-' + str(blocksize) + '_NEMIN-32'
        # print "datafile: ", datafile

        df = Datafile(datafile)
        # Get factor time
        v = df.get_value(t_facto_str)
        parsec_t_facto.append(float(v))
        # Get flops
        vf = df.get_value(flops_str)
        flops.append(float(vf))

    # SpLLT OpenMP (gnu)
    omp_t_facto = []
    for blocksize in blocksizes:
        # print "blocksize: ",blocksize
        datafile = outputdir + '/' + 'omp' + '/' + 'gnu' + '/' + pbl + '_NCPU-27' + '_NB-' + str(blocksize) + '_NEMIN-32'

        df = Datafile(datafile)
        # Get factor time
        v = df.get_value(t_facto_str)
        omp_t_facto.append(float(v))


    # MA87
    ma87_t_facto = []
    for blocksize in blocksizes:
        # print "blocksize: ",blocksize
        datafile = outputdir + '/' + 'ma87' + '/' + pbl + '_NCPU-28' + '_NB-' + str(blocksize) + '_NEMIN-32'
        
        df = Datafile(datafile)
        # Get factor time
        v = df.get_value(t_facto_str)
        ma87_t_facto.append(float(v))
    
    # print spllt_t_facto
    # print spllt_t_insert
    # print ma87_t_facto
        
    starpu_t_facto_idx = starpu_t_facto.index(min(starpu_t_facto))
    parsec_t_facto_idx = parsec_t_facto.index(min(parsec_t_facto))
    ma87_t_facto_idx = ma87_t_facto.index(min(ma87_t_facto))
    omp_t_idx = omp_t_facto.index(min(omp_t_facto))

    best_parsec_nb = blocksizes[parsec_t_facto_idx]
    best_parsec_t = parsec_t_facto[parsec_t_facto_idx]

    best_flops = flops[parsec_t_facto_idx]
    # # GFlops 
    best_gflops = best_flops/(1e9)

    best_ma87_nb = blocksizes[ma87_t_facto_idx]
    best_ma87_t = ma87_t_facto[ma87_t_facto_idx]

    best_omp_nb = blocksizes[omp_t_idx]
    best_omp_t = omp_t_facto[omp_t_idx]

    best_starpu_nb = blocksizes[starpu_t_facto_idx]
    best_starpu_t = starpu_t_facto[starpu_t_facto_idx]

    # Parsec and MA87, nb:time (txt)
    # print("%4s %10s %10s %10s %10s" % (matcount,
    #                                    best_parsec_nb,
    #                                    best_parsec_t,
    #                                    best_ma87_nb,
    #                                    best_ma87_t))

    # MA87 | OpenMP | StarPU | Parsec; nb:gflops (txt)
    print("%4s %10s %10.3f %10s %10.3f %10s %10.3f %10s %10.3f" % 
          (matcount,
           best_ma87_nb,
           best_gflops/best_ma87_t,
           best_omp_nb,
           best_gflops/best_omp_t,
           best_starpu_nb,
           best_gflops/best_starpu_t,
           best_parsec_nb,
           best_gflops/best_parsec_t,))

    # MA87, Parsec and OpenMP, nb:time (txt)
    # print("%4s %10s %10s %10s %10s %10s %10s" % (matcount,
    #                                              best_ma87_nb,
    #                                              best_ma87_t,
    #                                              best_parsec_nb,
    #                                              best_parsec_t,
    #                                              best_omp_nb,
    #                                              best_omp_t))

    # Parsec, OpenMP and MA87 nb:time (Latex)
    # print("%4s & %40s & %5s & %10.3f & %5s & %10.3f & %5s & %10.3f \\\\" % (matcount, lp.escape(mat),
    #                                                                         best_parsec_nb,
    #                                                                         best_parsec_t,
    #                                                                         best_omp_nb,
    #                                                                         best_omp_t,
    #                                                                         best_ma87_nb,
    #                                                                         best_ma87_t))
    
    # Parsec and MA87, nb:time (Latex)
    # print("%40s & %10s & %10s & %10s & %10s \\\\" % (lp.escape(mat),
    #                                                  best_spllt_nb,
    #                                                  lp.print_float(best_spllt_t_facto, 
    #                                                                 (best_spllt_t_facto<best_ma87_t_facto)),
    #                                                  best_ma87_nb,
    #                                                  lp.print_float(best_ma87_t_facto,
    #                                                                 (best_ma87_t_facto<best_spllt_t_facto))))

    # Parsec and MA87, nb:time (Latex)
    # print("%40s & %10s & %10s & %10s & %10s \\\\" % (lp.escape(mat),
    #                                                  best_spllt_nb,
    #                                                  lp.print_float(best_spllt_t_facto, 
    #                                                                 (best_spllt_t_facto<best_ma87_t_facto)),
    #                                                  best_ma87_nb,
    #                                                  lp.print_float(best_ma87_t_facto,
    #                                                                 (best_ma87_t_facto<best_spllt_t_facto))))

    matcount = matcount+1
