set terminal pngcairo size 1350,900 enhanced
set xlabel 't'
set ylabel 'Rg (promedio)'
set grid
set style line 1 lc rgb '#228B22' lw 2
set output 'results/Radio_giro.png'
plot \
'results/Variables_00001.dat' using 1:4 with lines ls 1 title 'Rg (promedio)',\
'results/Variables_00002.dat' using 1:4 with lines ls 1 notitle,\
'results/Variables_00003.dat' using 1:4 with lines ls 1 notitle,\
'results/Variables_00004.dat' using 1:4 with lines ls 1 notitle,\
'results/Variables_00005.dat' using 1:4 with lines ls 1 notitle,\
'results/Variables_00006.dat' using 1:4 with lines ls 1 notitle,\
'results/Variables_00007.dat' using 1:4 with lines ls 1 notitle,\
'results/Variables_00008.dat' using 1:4 with lines ls 1 notitle,\
'results/Variables_00009.dat' using 1:4 with lines ls 1 notitle,\
'results/Variables_00010.dat' using 1:4 with lines ls 1 notitle
unset output
