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

starpu_sched = 'lws'

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

    # spLLT PaRSEC
    spllt_parsec_t_facto = []
    spllt_parsec_t_insert = []
    for blocksize in blocksizes:
        # print "blocksize: ",blocksize
        datafile = outputdir + '/' + 'spllt_parsec' + '/' + pbl + '_NCPU-28' + '_NB-' + str(blocksize) + '_NEMIN-32'
        if os.path.exists(datafile):
            # print datafile
            f = open(datafile)
            v = vextractor.get_value(f, spllt_facto_time)
            f.seek(0)
            # print vi
            f.close()

            spllt_parsec_t_facto.append(float(v))

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

    # StarPU
    spllt_t_facto = []
    spllt_t_insert = []
    spllt_flops = []
    for blocksize in blocksizes:
        # print "blocksize: ",blocksize
        datafile = outputdir + '/' + 'spllt_starpu' + '/' + starpu_sched + '/' + pbl + '_NCPU-27' + '_NB-' + str(blocksize) + '_NEMIN-32'
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
    
    # print spllt_t_facto
    # print spllt_t_insert
    # print ma87_t_facto
        
    best_spllt_parsec_t_facto_idx = spllt_parsec_t_facto.index(min(spllt_parsec_t_facto))
    best_spllt_gnu_omp_t_facto_idx = spllt_gnu_omp_t_facto.index(min(spllt_gnu_omp_t_facto))
    best_spllt_t_facto_idx = spllt_t_facto.index(min(spllt_t_facto))
    best_ma87_t_facto_idx  = ma87_t_facto.index(min(ma87_t_facto))

    # Parsec
    best_spllt_parsec_nb       = blocksizes[best_spllt_parsec_t_facto_idx]
    best_spllt_parsec_t_facto  = spllt_parsec_t_facto[best_spllt_parsec_t_facto_idx]

    # omp
    best_spllt_gnu_omp_nb       = blocksizes[best_spllt_gnu_omp_t_facto_idx]
    best_spllt_gnu_omp_t_facto  = spllt_gnu_omp_t_facto[best_spllt_gnu_omp_t_facto_idx]

    # StarPU
    best_spllt_nb       = blocksizes[best_spllt_t_facto_idx]
    best_spllt_t_facto  = spllt_t_facto[best_spllt_t_facto_idx]
    best_spllt_t_insert = spllt_t_insert[best_spllt_t_facto_idx]
    best_spllt_flops    = spllt_flops[best_spllt_t_facto_idx]
    # GFlops 
    best_spllt_flops    = best_spllt_flops/(1024*1024*1024)

    # MA87
    best_ma87_nb       = blocksizes[best_ma87_t_facto_idx]
    best_ma87_t_facto  = ma87_t_facto[best_ma87_t_facto_idx]
    
    # print("%40s & %10s & %10s & %10s \\\\" % (lp.escape(mat), 
    #                                           lp.print_float(best_spllt_t_insert),
    #                                           lp.print_float(best_spllt_t_facto, (best_spllt_t_facto<best_ma87_t_facto)), 
    #                                           lp.print_float(best_ma87_t_facto , (best_ma87_t_facto<best_spllt_t_facto))))

    # Latex print

    # print("%40s & %10s & %10s & %10s & %10s & %10s & %10s \\\\" % (lp.escape(mat), 
    #                                                                best_spllt_gnu_omp_nb,
    #                                                                lp.print_float(best_spllt_gnu_omp_t_facto, 
    #                                                                               ((best_spllt_gnu_omp_t_facto<best_ma87_t_facto) and (best_spllt_gnu_omp_t_facto<best_spllt_t_facto))),
    #                                                                best_spllt_nb,
    #                                                                lp.print_float(best_spllt_t_facto, 
    #                                                                               ((best_spllt_t_facto<best_ma87_t_facto) and (best_spllt_t_facto<best_spllt_gnu_omp_t_facto))),
    #                                                                best_ma87_nb,
    #                                                                lp.print_float(best_ma87_t_facto,
    #                                                                               ((best_ma87_t_facto<best_spllt_t_facto) and (best_ma87_t_facto<best_spllt_gnu_omp_t_facto)))))

    # Latex print
    # print("%40s & %10s \\\\" % (lp.escape(mat), lp.print_float(best_spllt_t_facto, (best_spllt_t_facto<best_ma87_t_facto)))
    
    # data print
    # print("%2s %10s %10s %10s %10s" % (matcount, best_ma87_t_facto, best_spllt_gnu_omp_t_facto, best_spllt_t_facto, best_spllt_parsec_t_facto))

    # data print (factorization times)
    # print("%4s %10.3f %10.3f %10.3f" % (matcount, best_ma87_t_facto, best_spllt_gnu_omp_t_facto, best_spllt_t_facto))

    # data print (GFlop/s)
    print("%4s %10.3f %10.3f %10.3f" % (matcount, 
                                        (best_spllt_flops/best_ma87_t_facto), 
                                        (best_spllt_flops/best_spllt_gnu_omp_t_facto), 
                                        (best_spllt_flops/best_spllt_t_facto)))

    # data print (GFlop/s) with Parsec
    # print("%4s %10.3f %10.3f %10.3f %10.3f" % (matcount,
    #                                            (best_spllt_flops/best_ma87_t_facto), 
    #                                            (best_spllt_flops/best_spllt_gnu_omp_t_facto), 
    #                                            (best_spllt_flops/best_spllt_t_facto),
    #                                            (best_spllt_flops/best_spllt_parsec_t_facto)))

    # data print (factorization times and block sizes)
    # print("%4s %40s %6d %10.3f %6d %10.3f %6d %10.3f" % (matcount, lp.escape(mat), 
    #                                                      best_ma87_nb, best_ma87_t_facto, 
    #                                                      best_spllt_gnu_omp_nb, best_spllt_gnu_omp_t_facto, 
    #                                                      best_spllt_nb, best_spllt_t_facto))

    # Latex print (factorization times and block siezes)    
    # print("%4s & %40s & %6d & %10.3f & %6d & %10.3f & %6d & %10.3f \\\\" % (matcount, lp.escape(mat), 
    #                                                                         best_ma87_nb, best_ma87_t_facto, 
    #                                                                         best_spllt_gnu_omp_nb, best_spllt_gnu_omp_t_facto, 
    #                                                                         best_spllt_nb, best_spllt_t_facto))

    matcount = matcount+1 
