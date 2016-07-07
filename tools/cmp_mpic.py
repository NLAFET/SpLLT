#!/usr/bin/python3

import sys

import ValueExtractor as ve

class Table(ve.SeriesParser):
    class meta(ve.FileParser):
        extension = '.diag'
        n = ve.SingleInt(r'n, ne =', set_context=True)
        ne = ve.SingleInt(r'n, ne = *[0-9]+ ', set_context=True)
    class diag(ve.FileParser):
        extension = '.diag'
        nfact = ve.SingleInt(r'Number of entries in factors =')
        precond_time = ve.SingleFloat(r'Finding preconditioner took')
        bwderr = ve.SingleFloat(r'bwderr')
        pcg_time = ve.SingleFloat(r'PCG took')
        nitr = ve.SingleInt(r'Number of m-v products')
        ops = ve.Expression(lambda ne, nfact, nitr: (ne+nfact)*nitr, '%10.2e')
    class maxplus(ve.FileParser):
        extension = '.maxplus'
        nfact = ve.SingleInt(r'Number of entries in factors =')
        precond_time = ve.SingleFloat(r'Finding preconditioner took')
        bwderr = ve.SingleFloat(r'bwderr')
        pcg_time = ve.SingleFloat(r'PCG took')
        nitr = ve.SingleInt(r'Number of m-v products')
        ops = ve.Expression(lambda ne, nfact, nitr: (ne+nfact)*nitr, '%10.2e')
    class ic0(ve.FileParser):
        extension = '.mi22.nlvl0'
        nfact = ve.SingleInt(r'Number of entries in factors =')
        precond_time = ve.SingleFloat(r'Finding preconditioner took')
        bwderr = ve.SingleFloat(r'bwderr')
        pcg_time = ve.SingleFloat(r'PCG took')
        nitr = ve.SingleInt(r'Number of m-v products')
        ops = ve.Expression(lambda ne, nfact, nitr: (ne+nfact)*nitr, '%10.2e')
    class ic1(ve.FileParser):
        extension = '.mi22.nlvl1'
        nfact = ve.SingleInt(r'Number of entries in factors =')
        precond_time = ve.SingleFloat(r'Finding preconditioner took')
        bwderr = ve.SingleFloat(r'bwderr')
        pcg_time = ve.SingleFloat(r'PCG took')
        nitr = ve.SingleInt(r'Number of m-v products')
        ops = ve.Expression(lambda ne, nfact, nitr: (ne+nfact)*nitr, '%10.2e')
    class mi28(ve.FileParser):
        extension = '.mi28'
        nfact = ve.SingleInt(r'Number of entries in factors =')
        precond_time = ve.SingleFloat(r'Finding preconditioner took')
        bwderr = ve.SingleFloat(r'bwderr')
        pcg_time = ve.SingleFloat(r'PCG took')
        nitr = ve.SingleInt(r'Number of m-v products')
        ops = ve.Expression(lambda ne, nfact, nitr: (ne+nfact)*nitr, '%10.2e')

# Get list of problems
if sys.stdin.isatty():
    filename = 'list.mpic'
    if len(sys.argv) > 1:
        filename = sys.argv[1]
    with open(filename) as f:
        problems = [x.strip() for x in f.readlines()]
else:
    problems = [x.strip() for x in sys.stdin.readlines()]

# Print table
table = Table('mpic_1', problems)
# Everything
#table.print(['n', 'ne', 'nfact', 'precond_time', 'bwderr', 'pcg_time', 'nitr', 'ops'])
# Just interesting stuff
table.print(['n', 'ne', 'bwderr', 'pcg_time', 'ops'])
