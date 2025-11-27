set terminal pngcairo size 1350,900 enhanced
set xlabel 't'
set ylabel 'Longitud normalizada (promedio)'
set grid
set style line 1 lc rgb '#228B22' lw 2
stats 'results/Variables_00010.dat' using 4 nooutput
set output 'results/L_eff.png'
plot \
'results/Variables_00001.dat' using 1:6 with lines ls 1 title 'Longitud normalizada (promedio)',\
'results/Variables_00002.dat' using 1:6 with lines ls 1 notitle,\
'results/Variables_00003.dat' using 1:6 with lines ls 1 notitle,\
'results/Variables_00004.dat' using 1:6 with lines ls 1 notitle,\
'results/Variables_00005.dat' using 1:6 with lines ls 1 notitle,\
'results/Variables_00006.dat' using 1:6 with lines ls 1 notitle,\
'results/Variables_00007.dat' using 1:6 with lines ls 1 notitle,\
'results/Variables_00008.dat' using 1:6 with lines ls 1 notitle,\
'results/Variables_00009.dat' using 1:6 with lines ls 1 notitle,\
'results/Variables_00010.dat' using 1:6 with lines ls 1 notitle
unset output
