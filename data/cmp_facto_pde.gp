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
set style line 20 lt 2 ps 0.5 pt 2 lc rgb '#a2142f' # red
set style line 30 lt 2 ps 0.5 pt 3 lc rgb '#77ac30' # green
set style line 102 lc rgb '#999999' lt 0 lw 1

# Style
# set style data histograms
# set style histogram clustered gap 1
# set style fill solid 0.5

# font
set term pdf dashed font "Courier,12" 

# Legend
set key box bottom right w 1.1 font "Courier,12"

# Output
# set output "cmp_perf_helmholtz3d.pdf"
set output "cmp_perf_poisson3d.pdf"

# Title
set title "Poisson 3D - Factor. GFlop/s - 28 cores"

# Labels
set xlabel "Mesh size"
set ylabel "GFlop/s"

# Axis
set xrange [0:180]

# Set grid
set grid ytics lc rgbcolor "#000000" lt 0 lw 1

# # Plot
# plot 'cn255/data_cmp_perf_pde.dat' using 1:2 ls 1 w lp t 'MA87', \
#      ''                            using 1:3 ls 2 w lp t 'SpLLT-STF (OpenMP)', \
#      ''                            using 1:4 ls 3 w lp t 'SpLLT-STF (StarPU)'

# Plot
# With pruning
plot 'cn255/data_cmp_perf_poisson3d.dat' using 1:2 ls 1 w lp t 'MA87', \
     ''                                  using 1:4 ls 2 w lp t 'SpLLT-STF (OpenMP)', \
     ''                                  using 1:8 ls 20 w lp t 'SpLLT-STF w/ pruning (OpenMP)', \
     ''                                  using 1:6 ls 3 w lp t 'SpLLT-STF (StarPU)', \
     ''                                  using 1:10 ls 30 w lp t 'SpLLT-STF w/ pruning (StarPU)',

# plot Helmholtz 3D
# With pruning
# plot 'cn255/data_cmp_perf_helmholtz3d.dat' using 1:2 ls 1 w lp t 'MA87', \
#      ''                                    using 1:4 ls 2 w lp t 'SpLLT-STF (OpenMP)', \
#      ''                                    using 1:8 ls 20 w lp t 'SpLLT-STF w/ pruning (OpenMP)', \
#      ''                                    using 1:6 ls 3 w lp t 'SpLLT-STF (StarPU)', \
#      ''                                    using 1:10 ls 30 w lp t 'SpLLT-STF w/ pruning (StarPU)',
