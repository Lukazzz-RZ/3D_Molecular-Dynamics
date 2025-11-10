set ytics nomirror
plot 'results/rk2_estocastico_dt0.0010_gamma0.0_A1.0_KbT0.2000.txt' using 1:6 with lines notitle,\
     'results/rk2_estocastico_dt0.0010_gamma0.0_A1.0_KbT0.2000.txt' using 1:6 with lines notitle
