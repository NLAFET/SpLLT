reset

set style line 1 lt 1 ps 0.5 pt 1 lc rgb '#0072bd' # blue
set style line 2 lt 2 ps 0.5 pt 2 lc rgb '#a2142f' # red
set style line 3 lt 3 ps 0.5 pt 3 lc rgb '#77ac30' # green
set style line 4 lt 4 ps 0.5 pt 4 lc rgb '#edb120' # yellow
set style line 5 lt 5 lc rgb '#4dbeee' # light-blue
set style line 6 lt 1 lc rgb '#d95319' # orange
set style line 7 lt 1 lc rgb '#7e2f8e' # purple
set style line 8 lt 1 lc rgb '#000000' # black
set style line 102 lc rgb '#999999' lt 0 lw 1

set term svg size 1400, 1000 font "Courier,16" 
# filename
# set output "cmp_eff.pdf"
set output "cmp_eff.svg"

set multiplot layout 2, 2

# Style
set style data histograms
set style histogram clustered gap 1
set style fill solid 0.5

# put label box
# set key box left top w 1.1 font "Courier,12"


# figure title
set title "Factorization GFlop/s - 28 cores"

# Axes
set xlabel "Matrix \#"

set yrange [0:1.1]
# set xrange [-1:10]

set arrow 1 from -1,1 to 12,1 nohead  lt 1 lc rgb "#000000"  lw 2

# set boxwidth 1.0 absolute

set grid ytics lc rgbcolor "#000000" lt 0 lw 1

set title "Parallel efficiency"
plot 'cn255/data_cmp_eff2.dat' using 2:xtic(1) ls 8 fill solid 0.5 t ''

set title "Task efficiency"
plot 'cn255/data_cmp_eff2.dat' using 3:xtic(1) ls 1 fill solid 0.5 t ''

set title "Runtime efficiency"
plot 'cn255/data_cmp_eff2.dat' using 4:xtic(1) ls 3 fill solid 0.5 t ''

set title "Pipeline efficiency"
plot 'cn255/data_cmp_eff2.dat' using 5:xtic(1) ls 4 fill solid 0.5 t ''

unset multiplot
