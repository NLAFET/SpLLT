import re

def print_float(v, bold=False):
    if (bold):
        return "\\bf {0:.3f}".format(float(v))
    else:
        return "{0:.3f}".format(float(v))

# scientific notation
def print_float_sn(v):
    print ("{0:.3e}".format(float(v)))

# scientific notation
def escape(s):
    return re.sub(r'_', '\_', s)


