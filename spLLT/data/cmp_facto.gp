reset

set style line 1 lt 1 lc rgb '#0072bd' # blue
set style line 2 lt 1 lc rgb '#a2142f' # red
set style line 3 lt 1 lc rgb '#77ac30' # green
set style line 4 lt 1 lc rgb '#edb120' # yellow
set style line 5 lt 1 lc rgb '#4dbeee' # light-blue
set style line 6 lt 1 lc rgb '#d95319' # orange
set style line 7 lt 1 lc rgb '#7e2f8e' # purple
set style line 8 lt 1 lc rgb '#000000' # black
set style line 102 lc rgb '#999999' lt 0 lw 1

set term pdf font "Courier,12" 
set style data histograms
set style histogram clustered gap 1
set style fill solid 0.5

set key box left top w 1.1 font "Courier,12"

# set output "cmp_facto_all.pdf"
set output "cmp_facto_all0.pdf"
# set output "cmp_facto_stf.pdf"

set title "Factorization times"

set xlabel "Matrix \#"
set ylabel "time (s)"

set yrange [0:14]

set boxwidth 1.0 absolute

# STF
# plot 'data_cmp_facto.dat' using 2:xtic(1) ls 1 t 'MA87', \
#      ''                   using 3:xtic(1) ls 2 t 'SpLLT-STF (OpenMP)', \
#      ''                   using 4:xtic(1) ls 3 t 'SpLLT-STF (StarPU)'

# all 0
plot 'data_cmp_facto.dat' using 2:xtic(1) ls 1 t 'MA87', \
     ''                   using 3:xtic(1) ls 2 t 'SpLLT-STF (OpenMP)', \
     ''                   using 4:xtic(1) ls 3 t 'SpLLT-STF (StarPU)', \
     ''                   using 6:xtic(1) ls 7 t ''

# all
# plot 'data_cmp_facto.dat' using 2:xtic(1) ls 1 t 'MA87', \
#      ''                   using 3:xtic(1) ls 2 t 'SpLLT-STF (OpenMP)', \
#      ''                   using 4:xtic(1) ls 3 t 'SpLLT-STF (StarPU)', \
#      ''                   using 5:xtic(1) ls 7 t 'SpLLT-PTG (PaRSEC)'

set output
