reset

# Colors
set style line 1 lt 1 ps 0.5 pt 1 lc rgb '#0072bd' # blue
set style line 2 lt 1 ps 0.5 pt 2 lc rgb '#a2142f' # red
set style line 3 lt 1 ps 0.5 pt 3 lc rgb '#77ac30' # green
set style line 4 lt 1 ps 0.5 pt 4 lc rgb '#edb120' # yellow
set style line 5 lt 1 lc rgb '#4dbeee' # light-blue
set style line 6 lt 1 lc rgb '#d95319' # orange
set style line 7 lt 1 lc rgb '#7e2f8e' # purple
set style line 8 lt 1 lc rgb '#000000' # black
set style line 10 lt 4 ps 0.5 pt 2 lc rgb '#0072bd' # blue
set style line 20 lt 4 ps 0.5 pt 2 lc rgb '#a2142f' # red
set style line 30 lt 4 ps 0.5 pt 3 lc rgb '#77ac30' # green
set style line 102 lc rgb '#999999' lt 0 lw 1


set term pdf dashed font "Courier,12" 

set key box left top w 1.1 font "Courier,12"
# set key box right bottom w 1.1 font "Courier,12"

set output "cmp_prune.pdf"
# set output "cmp_prune_stf-starpu.pdf"
# set output "cmp_prune_stf-openmp.pdf"
# set output "cmp_prune_stf-openmp_big.pdf"
# set output "cmp_prune_stf-openmp_small.pdf"

set title "Factorization GFlop/s - 28 cores"

# Set axis label
set xlabel "Matrix \#"
set ylabel "GFlop/s"

set xrange [1:37]
# set xrange [1:38]
# set xrange [1:25]
# set xrange [26:38]

# set yrange [*:700]

set boxwidth 1.0 absolute

set grid ytics lc rgbcolor "#000000" lt 0 lw 1

# StarPU (wo/ and w/ pruning) and OpenMP (wo/ and w/ pruning)
plot 'cn255/data_cmp_prune_perf.dat' using ($0+1):3 ls 1  w lp t 'SpLLT OpenMP           ', \
     ''                              using ($0+1):5 ls 10 w lp t 'SpLLT OpenMP no prunnig', \
     ''                              using ($0+1):7 ls 2  w lp t 'SpLLT StarPU           ', \
     ''                              using ($0+1):9 ls 20 w lp t 'SpLLT StarPU no prunnig', \


# StarPU (wo/ and w/ pruning) vs MA87
# plot 'cn255/data_cmp_prune_perf.dat' using ($0+1):2 ls 1 w lp t 'MA87', \
#      ''                              using ($0+1):5 ls 2 w lp t 'SpLLT-STF (StarPU)           ', \
#      ''                              using ($0+1):6 ls 3 w lp t 'SpLLT-STF (StarPU) \w pruning'

# # OMP (wo/ and w/ pruning) vs MA87
# plot 'cn255/data_cmp_prune_perf.dat' using ($0+1):2 ls 1 w lp t 'HSL_MA87', \
#      ''                              using ($0+1):3 ls 2 w lp t 'SpLLT-STF (OpenMP)             ', \
     # ''                              using ($0+1):4 ls 3 w lp t 'SpLLT-STF (OpenMP) with pruning'

# # OMP (wo/ and w/ pruning) vs MA87 (big)
# plot 'cn255/data_cmp_prune_perf_big.dat' using ($0+26):2 ls 1 w lp t 'MA87', \
#      ''                                  using ($0+26):3 ls 2 w lp t 'SpLLT-STF (OpenMP)           ', \
#      ''                                  using ($0+26):4 ls 3 w lp t 'SpLLT-STF (OpenMP) \w pruning'

# OMP (wo/ and w/ pruning) vs MA87 (small)
# plot 'cn255/data_cmp_prune_perf_small.dat' using ($0+1):2 ls 1 w lp t 'MA87', \
#      ''                                    using ($0+1):3 ls 2 w lp t 'SpLLT-STF (OpenMP)           ', \
#      ''                                    using ($0+1):4 ls 3 w lp t 'SpLLT-STF (OpenMP) \w pruning'
