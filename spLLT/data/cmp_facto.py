import fileinput
import re

import vextractor
import latexprint as lp

spllt_facto_time = '\[>\] \[factorize\] time:'
ma87_facto_time  = 'Factor took'

blocksizes = [64, 128, 256, 512]

outputdir = 'cn255'

for mat in fileinput.input():

    print mat
    pbl = re.sub(r'/', '_', mat)
    pbl = pbl.rstrip()
    print pbl

    # spLLT
    spllt_t_facto = []
    for blocksize in blocksizes:
        print "blocksize: ",blocksize
        datafile = outputdir + '/' + 'spllt_starpu' + '/' + pbl + '_NCPU-27' + '_NB-' + str(blocksize)
        # print datafile
        f = open(datafile)
        v = vextractor.get_value(f, spllt_facto_time)
        f.close()

        spllt_t_facto.append(float(v))

    # MA87
    ma87_t_facto = []
    for blocksize in blocksizes:
        print "blocksize: ",blocksize
        datafile = outputdir + '/' + 'ma87' + '/' + pbl + '_NCPU-27' + '_NB-' + str(blocksize)
        # print datafile
        f = open(datafile)
        v = vextractor.get_value(f, ma87_facto_time)
        f.close()

        ma87_t_facto.append(float(v))
    
    print spllt_t_facto
    print ma87_t_facto
        
