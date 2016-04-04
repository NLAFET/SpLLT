def parsetrace(input_file):

    alltimes = {}
    allstatu = {}
    allbegin = {}
    allworke  = {}
     
    for line in input_file:
        sl = line.split()
        if len(sl)==0:
            continue
        if sl[0]=='7' :
            # print sl[2].find('t', 0, 0), sl[2]
            # if sl[2].find('t') == 0:
            if sl[3] == 'W':
                allworke[sl[2]] = sl[5]
            if sl[2].find('w') == 0:
            #if sl[2].find(argv[1]) == 0:
                if not alltimes.has_key(sl[2]):
                    alltimes[sl[2]] = {} #{'I':0, 'D':0, 'Fi':0, 'Po':0, 'C':0, 'B':0, 'E':0, 'Sc':0, 'Sl':0, 'P':0, 'U':0}
                    allstatu[sl[2]] = 'no'
                    allbegin[sl[2]] = -1
            if sl[2].find('t') == 0:
            #if sl[2].find(argv[1]) == 0:
                if not alltimes.has_key(sl[2]):
                    alltimes[sl[2]] = {} #{'I':0, 'D':0, 'Fi':0, 'Po':0, 'C':0, 'B':0, 'E':0, 'Sc':0, 'Sl':0, 'P':0, 'U':0}
                    allstatu[sl[2]] = 'no'
                    allbegin[sl[2]] = -1
        elif sl[0]=='10':
            if alltimes.has_key(sl[2]):
                # print 'Fucking error 1 '+sl[2]
            # else:
                thread = sl[2]
                t      = float(sl[1])
                if allstatu[thread] != 'no':
                    if alltimes[thread].has_key(allstatu[thread]):
                        alltimes[thread][allstatu[thread]] += t-float(allbegin[thread])
                    else:
                        alltimes[thread][allstatu[thread]]  = t-float(allbegin[thread])
                else:
                    alltimes[thread]['P']  = t-float(allbegin[thread])
     
                allbegin[thread] = t
                allstatu[thread] = sl[4]
                # print 'Fucking error 1 '+sl[4]                
        elif sl[0]=='9':
            if sl[4]=='stop_profiling':
                stopt = float(sl[1])
     
    # for th in alltimes:
        # if alltimes[thread].has_key('Sl'):
            # alltimes[th]['Sl'] += stopt-float(allbegin[th])
        # else:
            # alltimes[th]['Sl'] = stopt-float(allbegin[th])
            
    glob_tt_time = 0
    glob_ov_time = 0
    glob_sl_time = 0
    glob_ex_time = 0

    for th in alltimes:
        if  th.find('w') == 0:
            if allworke[th].find('CUDA') == 0:
                alltimes[th] = {} 
                # print th ,alltimes[th]

    # for th in alltimes:
    #     print th, alltimes[th]

    for th in alltimes:
        if alltimes[th]:
            # print th
            th_tt_time  = sum(alltimes[th].values())

            th_ex_time  = 0

            # if(argv[1]=='w'):
            if(alltimes[th].has_key('FACTO_BLK'     )): th_ex_time += alltimes[th]['FACTO_BLK'     ]
            if(alltimes[th].has_key('SOLVE_BLK'     )): th_ex_time += alltimes[th]['SOLVE_BLK'     ]
            if(alltimes[th].has_key('UPDATE_BLK'    )): th_ex_time += alltimes[th]['UPDATE_BLK'    ]
            if(alltimes[th].has_key('UPDATE_BETWEEN')): th_ex_time += alltimes[th]['UPDATE_BETWEEN']

            if (alltimes[th].has_key('E')):
                th_ex_time = alltimes[th]['E']
            
            if (alltimes[th].has_key('Sl')):
                th_sl_time  = alltimes[th]['Sl']

            #elif(argv[1]=='t'):
            #    th_ex_time = alltimes[th]['E']
            th_ov_time = th_tt_time-th_sl_time-th_ex_time
            glob_tt_time += th_tt_time
            glob_ex_time += th_ex_time
            glob_ov_time += th_ov_time
            glob_sl_time += th_sl_time
            # print th,' ',th_ex_time/th_tt_time,th_sl_time/th_tt_time,th_ov_time/th_tt_time

        # else:
        #     print th
    # print '\nGlobal: ',glob_ex_time/glob_tt_time,glob_sl_time/glob_tt_time,glob_ov_time/glob_tt_time

    r_ex_time = glob_ex_time/glob_tt_time
    r_sl_time = glob_sl_time/glob_tt_time
    r_ov_time = glob_ov_time/glob_tt_time

    return r_ex_time, r_sl_time, r_ov_time, glob_tt_time
