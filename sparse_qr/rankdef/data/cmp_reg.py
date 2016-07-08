#!/usr/bin/env python

import sys
import os
import fileinput
import re

import vextractor as ve

sigmas = [-6, -12]

flop_str     = 'Estimated total flops at facto               :'
flop_reg_str = 'Estimated total flops at facto \(regularized\) :'
xnorm_str     = 'xnorm :'
resid_str     = 'resid :'

for mat in open('list.matrix'):

    mat = mat.rstrip()
    pbl = re.sub(r'.mtx', '', mat)
    # print pbl
    
    for sigma in sigmas:

        # print sigma

        datafile = pbl + '_SIGMA-1e' + str(sigma)
        f = open(datafile)
        fl = ve.get_int_value(f, flop_str)
        f.seek(0)
        flr = ve.get_int_value(f, flop_reg_str)
        f.seek(0)
        xnorm = ve.get_value(f, xnorm_str)
        f.seek(0)
        resid = ve.get_value(f, resid_str)
        f.close()
                    
        # print fl, flr
        # print xnorm, resid
        
        # flops = ve.SingleFloat(r'Estimated total flops at facto               :')
        # print flops

        print("1e%2s %20s %20s %10s %10s" % (sigma, fl, flr, xnorm, resid))
