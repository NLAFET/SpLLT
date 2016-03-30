import fileinput
import re

import vextractor
import latexprint as lp

spllt_facto_time = '\[>\] \[factorize\] time:'
ma87_facto_time  = 'Factor took'
spllt_task_insert_time = '\[>\] \[spllt_stf_factorize\] task insert time:'

blocksizes = [64, 128, 256, 512]

outputdir = 'cn255'

for mat in fileinput.input():

    mat = mat.rstrip()

    # print mat
    pbl = re.sub(r'/', '_', mat)
    pbl = pbl.rstrip()

    # spLLT
    spllt_t_facto = []
    spllt_t_insert = []
    for blocksize in blocksizes:
        # print "blocksize: ",blocksize
        datafile = outputdir + '/' + 'spllt_starpu' + '/' + pbl + '_NCPU-27' + '_NB-' + str(blocksize)
        # print datafile
        f = open(datafile)
        v = vextractor.get_value(f, spllt_facto_time)
        f.seek(0)
        vi = vextractor.get_value(f, spllt_task_insert_time)
        # print vi
        f.close()

        spllt_t_facto.append(float(v))
        spllt_t_insert.append(float(vi))

    # MA87
    ma87_t_facto = []
    for blocksize in blocksizes:
        # print "blocksize: ",blocksize
        datafile = outputdir + '/' + 'ma87' + '/' + pbl + '_NCPU-27' + '_NB-' + str(blocksize)
        # print datafile
        f = open(datafile)
        v = vextractor.get_value(f, ma87_facto_time)
        f.close()

        ma87_t_facto.append(float(v))
    
    # print spllt_t_facto
    # print spllt_t_insert
    # print ma87_t_facto
        
    best_spllt_t_facto_idx = spllt_t_facto.index(min(spllt_t_facto))
    best_ma87_t_facto_idx  = ma87_t_facto.index(min(ma87_t_facto))

    best_spllt_nb       = blocksizes[best_spllt_t_facto_idx]
    best_spllt_t_facto  = spllt_t_facto[best_spllt_t_facto_idx]
    best_spllt_t_insert = spllt_t_insert[best_spllt_t_facto_idx]

    best_ma87_nb       = blocksizes[best_ma87_t_facto_idx]
    best_ma87_t_facto  = ma87_t_facto[best_ma87_t_facto_idx]
    
    print("%40s & %10s & %10s \\\\" % (lp.escape(mat), 
                                       lp.print_float(best_spllt_t_facto, (best_spllt_t_facto<best_ma87_t_facto)), 
                                       lp.print_float(best_ma87_t_facto , (best_ma87_t_facto<best_spllt_t_facto))))
    
    # print("%40s & %10s \\\\" % (lp.escape(mat), lp.print_float(best_spllt_t_facto, (best_spllt_t_facto<best_ma87_t_facto)))
