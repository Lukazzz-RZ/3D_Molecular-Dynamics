set boxwidth 0.9 relative
set style fill solid 0.5
set yrange [0:*]
set xlabel 'Bins'
set ylabel 'Frecuencia Normalizada'
set title 'Histograma de posicion'
plot 'histograma.dat' using 1:2 with boxes notitle
