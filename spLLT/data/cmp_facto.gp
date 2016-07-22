reset

set style line 1 lt 1 ps 0.5 pt 1 lc rgb '#0072bd' # blue
set style line 2 lt 2 ps 0.5 pt 2 lc rgb '#a2142f' # red
set style line 3 lt 3 ps 0.5 pt 3 lc rgb '#77ac30' # green
set style line 4 lt 4 lc rgb '#edb120' # yellow
set style line 5 lt 5 lc rgb '#4dbeee' # light-blue
set style line 6 lt 1 lc rgb '#d95319' # orange
set style line 7 lt 1 lc rgb '#7e2f8e' # purple
set style line 8 lt 1 lc rgb '#000000' # black
set style line 102 lc rgb '#999999' lt 0 lw 1

set term pdf font "Courier,12" 

# Style
# set style data histograms
# set style histogram clustered gap 1
# set style fill solid 0.5

set key box left top w 1.1 font "Courier,12"

# set output "cmp_facto_all.pdf"
# set output "cmp_facto_all0.pdf"
# set output "cmp_facto_stf.pdf"
set output "cmp_perf_stf.pdf"
# set output "cmp_facto_rel_stf.pdf"

# set title "Factorization times - 28 cores"
set title "Factorization GFlop/s - 28 cores"
# set title "Relative performance with MA87 - 28 cores"

set xlabel "Matrix \#"
set ylabel "GFlop/s"
# set ylabel "time (s)"
# set ylabel "perf SpLLT / perf MA87"

# set yrange [0:1.5]
# set logscale y 10

set xrange [1:38]
# set xtics 1
# set ytics 0.2

set boxwidth 1.0 absolute

# set arrow from 1,1 to 38,1 nohead  lt 0 lc rgb "#000000"  lw 2

# set grid ytics lc rgbcolor "#000000" lt 0 lw 1

# STF (times)
# plot 'cn255/data_cmp_facto.dat' using 2:xtic(1) ls 1 t 'MA87', \
#      ''                         using 3:xtic(1) ls 2 t 'SpLLT-STF (OpenMP)', \
#      ''                         using 4:xtic(1) ls 3 t 'SpLLT-STF (StarPU)'

# STF (perf GFlop/s)
plot 'cn255/data_cmp_perf.dat' using ($0+1):2 ls 1 w lp t 'MA87', \
     ''                        using ($0+1):3 ls 2 w lp t 'SpLLT-STF (OpenMP)', \
     ''                        using ($0+1):4 ls 3 w lp t 'SpLLT-STF (StarPU)'

# STF perf/times relative to MA87
# plot 'cn255/data_cmp_perf.dat' using ($0+1):($3/$2) ls 2 w lp t 'SpLLT-STF (OpenMP)', \
#      ''                        using ($0+1):($4/$2) ls 3 w lp t 'SpLLT-STF (StarPU)'

# all 0
# plot 'cn255/data_cmp_facto.dat' using 2:xtic(1) ls 1 t 'MA87', \
#      ''                         using 3:xtic(1) ls 2 t 'SpLLT-STF (OpenMP)', \
#      ''                         using 4:xtic(1) ls 3 t 'SpLLT-STF (StarPU)', \
#      ''                         using 6:xtic(1) ls 7 t ''

# all
# plot 'data_cmp_facto.dat' using 2:xtic(1) ls 1 t 'MA87', \
#      ''                   using 3:xtic(1) ls 2 t 'SpLLT-STF (OpenMP)', \
#      ''                   using 4:xtic(1) ls 3 t 'SpLLT-STF (StarPU)', \
#      ''                   using 5:xtic(1) ls 7 t 'SpLLT-PTG (PaRSEC)'

set output
