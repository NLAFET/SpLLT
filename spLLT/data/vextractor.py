#!/usr/bin/python

import re

def get_value(f, pattern):
    for line in f:
        line = line.rstrip()
        if re.search(pattern, line):
            value = re.findall('[^0-9]*([0-9E.+-]*).*', line)
            # value = re.findall('-?[1-9]+[0-9]*.?[0-9]*E-?\+?[0-9]+', line)
            # value = re.findall('\s-?[1-9]+[0-9]*.?[0-9]*E-?\+?[0-9]+\s?', line)
            print line
            return value[0]
