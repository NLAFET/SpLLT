#!/usr/bin/python
import sys, getopt
import traceutils as tu

def main(argv):

    input_file = open(argv[0], 'r')
    
    r_ex_time = 0.0
    r_sl_time = 0.0
    r_ov_time = 0.0    
    glob_tt_time = 0.0

    r_ex_time, r_sl_time, r_ov_time, glob_tt_time = tu.parsetrace(input_file)
    
    glob_ex_time = r_ex_time*glob_tt_time
    glob_sl_time = r_sl_time*glob_tt_time
    glob_ov_time = r_ov_time*glob_tt_time
        
    print("%-50s        (exec., idle, ohead): %10.0f (%6.3f),  %10.0f (%6.3f),   %10.0f (%6.3f)  = %10.0f" % (argv[0],
                                                                                                              glob_ex_time, r_ex_time,
                                                                                                              glob_sl_time, r_sl_time,
                                                                                                              glob_ov_time, r_ov_time,
                                                                                                              glob_tt_time))

if __name__ == "__main__":
   main(sys.argv[1:])    
