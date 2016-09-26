reset

set style line 1 lt 1 ps 0.5 pt 1 lc rgb '#0072bd' # blue
set style line 2 lt 2 ps 0.5 pt 2 lc rgb '#a2142f' # red
set style line 3 lt 3 ps 0.5 pt 3 lc rgb '#77ac30' # green
set style line 4 lt 4 ps 0.5 pt 4 lc rgb '#edb120' # yellow
set style line 5 lt 5 ps 0.5 pt 5 lc rgb '#4dbeee' # light-blue
set style line 6 lt 1 lc rgb '#d95319' # orange
set style line 7 lt 1 lc rgb '#7e2f8e' # purple
set style line 8 lt 1 lc rgb '#000000' # black
set style line 102 lc rgb '#999999' lt 0 lw 1

set term pdf font "Courier,12" 

set key box left top w 1.1 font "Courier,12"

set output "cmp_prune.pdf"

set title "Factorization GFlop/s - 28 cores"

# Set axis label
set xlabel "Matrix \#"
set ylabel "GFlop/s"

set xrange [1:38]

set boxwidth 1.0 absolute

set grid ytics lc rgbcolor "#000000" lt 0 lw 1

plot 'cn255/data_cmp_prune_perf.dat' using ($0+1):2 ls 1 w lp t 'MA87', \
     ''                              using ($0+1):3 ls 2 w lp t 'SpLLT-STF (OpenMP)', \
     ''                              using ($0+1):4 ls 3 w lp t 'SpLLT-STF (OpenMP with tree pruning)', \
     ''                              using ($0+1):5 ls 4 w lp t 'SpLLT-STF (StarPU)', \
     ''                              using ($0+1):6 ls 5 w lp t 'SpLLT-STF (StarPU with tree pruning)'     
