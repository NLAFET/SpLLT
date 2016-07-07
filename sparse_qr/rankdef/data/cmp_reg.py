#!/usr/bin/env python

import sys
import os
import fileinput
import re

import vextractor as ve

sigmas = [-6, -12]

flop_str     = 'Estimated total flops at facto               :'
flop_reg_str = 'Estimated total flops at facto (regularized) :'

for mat in open('list.matrix'):

    mat = mat.rstrip()
    pbl = re.sub(r'.mtx', '', mat)
    print pbl
    
    for sigma in sigmas:

        # print sigma

        datafile = pbl + '_SIGMA-1e' + str(sigma)
        f = open(datafile)
        # fl = ve.get_value(f, flop_str)
        # f.seek(0)
        flr = ve.get_value(f, flop_reg_str)
        f.close()
                    
        print flr
        
        # flops = ve.SingleFloat(r'Estimated total flops at facto               :')
        # print flops
