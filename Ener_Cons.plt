set terminal pngcairo size 1350,900 enhanced
set xlabel 't'
set ylabel 'Energía'
set grid
set yrange [0:100]
set style line 1 lc rgb '#dd181f' lw 2 lt 1
set style line 2 lc rgb '#0060ad' lw 2 lt 1
set style line 3 lc rgb '#800080' lw 2 lt 1
set output 'results/Energias_Medias.png'
plot \
'results/Variables_00001.dat' using 1:2 with lines ls 1 title 'E. cinética (promedio)',\
'results/Variables_00001.dat' using 1:3 with lines ls 2 title 'E. potencial (promedio)',\
'results/Variables_00001.dat' using 1:($2+$3) with lines ls 3 title 'E. total (promedio)',\
'results/Variables_00002.dat' using 1:2 with lines ls 1 notitle,\
'results/Variables_00002.dat' using 1:3 with lines ls 2 notitle,\
'results/Variables_00002.dat' using 1:($2+$3) with lines ls 3 notitle,\
'results/Variables_00003.dat' using 1:2 with lines ls 1 notitle,\
'results/Variables_00003.dat' using 1:3 with lines ls 2 notitle,\
'results/Variables_00003.dat' using 1:($2+$3) with lines ls 3 notitle,\
'results/Variables_00004.dat' using 1:2 with lines ls 1 notitle,\
'results/Variables_00004.dat' using 1:3 with lines ls 2 notitle,\
'results/Variables_00004.dat' using 1:($2+$3) with lines ls 3 notitle,\
'results/Variables_00005.dat' using 1:2 with lines ls 1 notitle,\
'results/Variables_00005.dat' using 1:3 with lines ls 2 notitle,\
'results/Variables_00005.dat' using 1:($2+$3) with lines ls 3 notitle,\
'results/Variables_00006.dat' using 1:2 with lines ls 1 notitle,\
'results/Variables_00006.dat' using 1:3 with lines ls 2 notitle,\
'results/Variables_00006.dat' using 1:($2+$3) with lines ls 3 notitle,\
'results/Variables_00007.dat' using 1:2 with lines ls 1 notitle,\
'results/Variables_00007.dat' using 1:3 with lines ls 2 notitle,\
'results/Variables_00007.dat' using 1:($2+$3) with lines ls 3 notitle,\
'results/Variables_00008.dat' using 1:2 with lines ls 1 notitle,\
'results/Variables_00008.dat' using 1:3 with lines ls 2 notitle,\
'results/Variables_00008.dat' using 1:($2+$3) with lines ls 3 notitle,\
'results/Variables_00009.dat' using 1:2 with lines ls 1 notitle,\
'results/Variables_00009.dat' using 1:3 with lines ls 2 notitle,\
'results/Variables_00009.dat' using 1:($2+$3) with lines ls 3 notitle,\
'results/Variables_00010.dat' using 1:2 with lines ls 1 notitle,\
'results/Variables_00010.dat' using 1:3 with lines ls 2 notitle,\
'results/Variables_00010.dat' using 1:($2+$3) with lines ls 3 notitle
unset output
