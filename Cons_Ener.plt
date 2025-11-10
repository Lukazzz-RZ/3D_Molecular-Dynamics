set terminal pngcairo size 1350,900 enhanced font 'Helvetica,16'
set output 'results/Energias_verlet_estocastico_dt0.0100_gamma1.0000.png'
set xlabel 't [s]'
set ylabel 'En [J]'
set autoscale x
set autoscale y
set style line 1 lc rgb '#dd181f' lw 3 lt 1
set style line 2 lc rgb '#0060ad' lw 3 lt 1
set style line 3 lc rgb '#800080' lw 3 lt 1
plot \
'results/verlet_estocastico_dt0.0100_gamma1.0000.txt' every ::1 using 1:4 with lines ls 1 title 'Energía cinética', \
'results/verlet_estocastico_dt0.0100_gamma1.0000.txt' every ::1 using 1:5 with lines ls 2 title 'Energía potencial', \
'results/verlet_estocastico_dt0.0100_gamma1.0000.txt' every ::1 using 1:6 with lines ls 3 title 'Energía total'
unset output
set output 'results/PosicionVel_verlet_estocastico_dt0.0100_gamma1.0000.png'
set xlabel 't [s]'
set ylabel 'Pos/Vel []'
set autoscale x
set autoscale y
set style line 1 lc rgb '#dd181f' lw 3 lt 1
set style line 3 lc rgb '#800080' lw 3 lt 1
plot \
'results/verlet_estocastico_dt0.0100_gamma1.0000.txt' every ::1 using 1:2 with lines ls 1 title 'Posición', \
'results/verlet_estocastico_dt0.0100_gamma1.0000.txt' every ::1 using 1:3 with lines ls 3 title 'Velocidad'
unset output
set output 'results/Energias_Medias_verlet_estocastico_dt0.0100_gamma1.0000.png'
set xlabel 't [s]'
set ylabel 'En [J]'
set autoscale x
set autoscale y
set style line 1 lc rgb '#dd181f' lw 3 lt 1
set style line 2 lc rgb '#0060ad' lw 3 lt 1
set style line 3 lc rgb '#800080' lw 3 lt 1
plot \
'results/verlet_estocastico_dt0.0100_gamma1.0000.txt' every ::1 using 1:7 with lines ls 1 title '<EC>', \
'results/verlet_estocastico_dt0.0100_gamma1.0000.txt' every ::1 using 1:8 with lines ls 2 title '<EP>', \
'results/verlet_estocastico_dt0.0100_gamma1.0000.txt' every ::1 using 1:9 with lines ls 3 title '<ET>'
unset output
