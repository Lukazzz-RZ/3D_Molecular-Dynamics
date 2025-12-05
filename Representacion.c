#include "head.h"

// FUNCIONES DE ESCRITURA //

FILE* crear_archivo_xyz(int bloque) {
    char nombre[64];
    sprintf(nombre, "results/Data_%05d.xyz", bloque);
    FILE* f = fopen(nombre, "w");
    if (!f) {
        perror("Error al crear archivo xyz");
        exit(1);
    }
    return f;
}

void guardar_bloque_xyz(Particula* P, FILE* f, int paso_global) {
    fprintf(f, "%d\nFrame %d\n", N_particulas, paso_global);
    for (int i = 0; i < N_particulas; i++) {
        
        char tipo;
        if (P[i].q == 0.0)
            tipo = 'R';   
        else
            tipo = 'B';    

        fprintf(f, "%c %.6f %.6f %.6f\n",
            tipo,
            P[i].pos.x,
            P[i].pos.y,
            P[i].pos.z);
    }
}

FILE* crear_archivo_variables(int bloque, const char* cabecera) {
    char nombre_fichero[64];
    sprintf(nombre_fichero, "results/Variables_%05d.dat", bloque);

    FILE* file = fopen(nombre_fichero, "w");
    if (!file) {
        perror("Error al crear el fichero de variables");
        exit(1);
    }
    fprintf(file, "%s\n", cabecera);

    return file;
}

void Escribir_buffer(FILE* f, double* buffer) {
    if (!f) {
        fprintf(stderr, "Error: fichero no abierto en escribir_buffer.\n");
        return;
    }

    fprintf(f, "%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\n",
        buffer[0], buffer[1], buffer[2],
        buffer[3], buffer[4], buffer[5]);
}

void Gnuplot_EnerCons(int bloques) {
    FILE* pfileout = fopen("Ener_Cons.plt", "w");
    if (!pfileout) {
        perror("Error al crear Ener_Cons.plt");
        return;
    }

    fprintf(pfileout,
        "set terminal pngcairo size 1350,900 enhanced\n"
        "set xlabel 't'\n"
        "set ylabel 'Energía'\n"
        "set grid\n"
        "set style line 1 lc rgb '#dd181f' lw 2 lt 1\n"
        "set style line 2 lc rgb '#0060ad' lw 2 lt 1\n"
        "set style line 3 lc rgb '#800080' lw 2 lt 1\n"
        "set output 'results/Energias_Medias.png'\n"
        "plot \\\n"
    );

    for (int b = 1; b <= bloques; b++) {
        char nombre_fich[64];
        snprintf(nombre_fich, sizeof(nombre_fich), "results/Variables_%05d.dat", b);

        fprintf(pfileout,
            "'%s' using 1:2 with lines ls 1 %s,\\\n" // Ecin
            "'%s' using 1:3 with lines ls 2 %s,\\\n" // Epot
            "'%s' using 1:($2+$3) with lines ls 3 %s",  // Etotal = Ecin + Epot
            nombre_fich, (b == 1 ? "title 'E. cinética (promedio)'" : "notitle"),
            nombre_fich, (b == 1 ? "title 'E. potencial (promedio)'" : "notitle"),
            nombre_fich, (b == 1 ? "title 'E. total (promedio)'" : "notitle")
        );

        if (b < bloques)
            fprintf(pfileout, ",\\\n");
        else
            fprintf(pfileout, "\n");
    }

    fprintf(pfileout, "unset output\n");
    fclose(pfileout);

    system("gnuplot -persist Ener_Cons.plt");
    printf("Grafica de energias generadas con exito.\n");
}

void Gnuplot_Rg(int bloques) {
    FILE* pfileout = fopen("Rg.plt", "w");
    if (!pfileout) {
        perror("Error al crear Rg.plt");
        return;
    }

    fprintf(pfileout,
        "set terminal pngcairo size 1350,900 enhanced\n"
        "set xlabel 't'\n"
        "set ylabel 'Rg (promedio)'\n"
        "set grid\n"
        "set style line 1 lc rgb '#228B22' lw 2\n"
    );

    // >>> Correcto: calcular rango ANTES del output y plot <<<
    fprintf(pfileout,
        "stats 'results/Variables_%05d.dat' using 4 nooutput\n"
        "set yrange [STATS_min*0.9:STATS_max*1.1]\n",
        bloques
    );

    // >>> Después ya viene el output y el plot <<<
    fprintf(pfileout,
        "set output 'results/Radio_giro.png'\n"
        "plot \\\n"
    );

    for (int b = 1; b <= bloques; b++) {
        char nombre_fich[64];
        snprintf(nombre_fich, sizeof(nombre_fich), "results/Variables_%05d.dat", b);

        fprintf(pfileout,
            "'%s' using 1:4 with lines ls 1 %s",
            nombre_fich,
            (b == 1 ? "title 'Rg (promedio)'" : "notitle")
        );

        if (b < bloques)
            fprintf(pfileout, ",\\\n");
        else
            fprintf(pfileout, "\n");
    }

    fprintf(pfileout, "unset output\n");
    fclose(pfileout);

    system("gnuplot -persist Rg.plt");
    printf("Grafica del radio de giro generada con exito.\n");
}

void Gnuplot_L_eff(int bloques) {
    FILE* pfileout = fopen("L_eff.plt", "w");
    if (!pfileout) {
        perror("Error al crear L_eff.plt");
        return;
    }

    fprintf(pfileout,
        "set terminal pngcairo size 1350,900 enhanced\n"
        "set xlabel 't'\n"
        "set ylabel 'Longitud normalizada (promedio)'\n"
        "set grid\n"
        "set style line 1 lc rgb '#228B22' lw 2\n"
    );

    // >>> Correcto: calcular rango ANTES del output y plot <<<
    fprintf(pfileout,
        "stats 'results/Variables_%05d.dat' using 4 nooutput\n",
        bloques
    );

    // >>> Después ya viene el output y el plot <<<
    fprintf(pfileout,
        "set output 'results/L_eff.png'\n"
        "plot \\\n"
    );

    for (int b = 1; b <= bloques; b++) {
        char nombre_fich[64];
        snprintf(nombre_fich, sizeof(nombre_fich), "results/Variables_%05d.dat", b);

        fprintf(pfileout,
            "'%s' using 1:6 with lines ls 1 %s",
            nombre_fich,
            (b == 1 ? "title 'Longitud normalizada (promedio)'" : "notitle")
        );

        if (b < bloques)
            fprintf(pfileout, ",\\\n");
        else
            fprintf(pfileout, "\n");
    }

    fprintf(pfileout, "unset output\n");
    fclose(pfileout);

    system("gnuplot -persist L_eff.plt");
    printf("Grafica de longitud efectiva generada con exito.\n");
}

void Ajuste_L_eff_en_Fx() {

    double beta = 1.0 / KbT; // 1/(k_B T)

    FILE* pfileout = fopen("L_eff_ajuste_Fx.plt", "w");
    if (!pfileout) {
        perror("Error al crear L_eff_ajuste_Fx.plt");
        return;
    }

    fprintf(pfileout,
        "set terminal pngcairo size 1350,900 enhanced font 'Arial,14'\n"
        "set output 'results/L_eff_vs_Fx.png'\n"
        "set title 'Longitud Efectiva vs Fx'\n"
        "set xlabel 'Fx'\n"
        "set ylabel 'L_{eff}'\n"
        "set grid\n"
        "set key bottom right\n"
        "set yrange [0:1.5]\n"
        "\n"
        "# Definir coth para Gnuplot\n"
        "coth(x) = cosh(x)/sinh(x)\n"
        "\n"
        "# Funcion de ajuste Langevin\n"
        "f_ajuste(Fx) = (coth(%f * Fx) - 1/(%f * Fx) + Fx/%f)\n"
        "\n"
        "plot f_ajuste(x) with lines lw 2 lc rgb '#00aa00' title 'Teoria Langevin', \\\n"
        "     'results/L_eff_vs_Fx.txt' using 1:2 with points pt 2 ps 3 lw 3 lc rgb '#dd181f' title 'N=4', \\\n"
        "     'results/L_eff_vs_Fx.txt' using 1:3 with points pt 1 ps 4 lw 3 lc rgb '#0060ad' title 'N=10', \\\n"
        "     'results/L_eff_vs_Fx.txt' using 1:4 with points pt 7 ps 1 lw 3 lc rgb '#FF8800' title 'N=30'\n",
        beta * b, beta * b, b*k
    );

    fclose(pfileout);
    system("gnuplot L_eff_ajuste_Fx.plt");

    printf("Grafica generada: results/L_eff_vs_Fx.png\n");
}


void Ajuste_Rg_en_N() {

    double beta = 1.0 / KbT;
    double l2 = (b * b * b * b + 6 * b * b / (beta * k) + 3.0 / (beta * k * beta * k)) /
        (b * b + 1.0 / (beta * k));

    FILE* pfileout = fopen("Rg_ajuste_N.plt", "w");
    if (!pfileout) {
        perror("Error al crear Rg_ajuste_N.plt");
        return;
    }

    fprintf(pfileout,
        "set terminal pngcairo size 1350,900 enhanced font 'Arial,14'\n"
        "set output 'results/Rg_vs_N.png'\n"
        "set title 'Radio de giro vs N'\n"
        "set xlabel 'N'\n"
        "set ylabel 'Rg^2'\n"
        "set grid\n"
        "set key top left\n"
        "\n"
        "f_ajuste(x) = %.10f*(x - 1.0/x)/6\n"
        "f_teorica(x) = %.10f * x / 6\n"
        "\n"
        "plot 'results/Rg_vs_N.txt' using 1:2 with points pt 2 ps 3 lw 3 lc rgb '#dd181f' title '<Rg^2>', \\\n"
        "     'results/Rg_vs_N.txt' using 1:3 with points pt 1 ps 3 lw 3 lc rgb '#9900FF' title '<Rg>^2', \\\n"
        "     f_ajuste(x) with lines lw 2 lc rgb '#0060ad' title 'Dependencia Teorica', \\\n"
        "     f_teorica(x) with lines lw 2 lc rgb '#00aa00' title 'Dependencia Aproximada'\n",
        l2, b * b
    );

    fclose(pfileout);
    system("gnuplot Rg_ajuste_N.plt");

    printf("Grafica generada: results/Rg_vs_N.png\n");
}

void Ajuste_Rg_en_k() {

    FILE* pfileout = fopen("Rg_ajuste_ke.plt", "w");
    if (!pfileout) {
        perror("Error al crear Rg_ajuste_ke.plt");
        return;
    }

    double beta = 1.0 / KbT;
    double b2 = b * b;
    double b4 = b2 * b2;

    fprintf(pfileout,
        "set terminal pngcairo size 1350,900 enhanced font 'Arial,14'\n"
        "set output 'results/Rg_vs_ke.png'\n"
        "set title 'Radio de giro vs k'\n"
        "set xlabel 'k'\n"
        "set ylabel 'Rg^2'\n"
        "set grid\n"
        "set key top left\n"
        "\n"
        "f_ajuste(x) = ((%f + 6*%f/(%f*x) + 3/( (%f*x)**2 )) / (%f + 1/(%f*x))) * ((%d - 1.0/%d)/6)\n"
        "f_teorica(x) = %f / 6 * (%d - 1.0/%d)\n"
        "\n"
        "plot 'results/Rg_vs_ke.txt' using 1:2 with points pt 2 ps 3 lw 3 lc rgb '#dd181f' title '<Rg^2>', \\\n"
        "     'results/Rg_vs_ke.txt' using 1:3 with points pt 1 ps 3 lw 3 lc rgb '#0060ad' title '<Rg>^2', \\\n"
        "     f_ajuste(x) with lines lw 2 lc rgb '#0060ad' title 'Dependencia Teorica', \\\n"
        "     f_teorica(x) with lines lw 2 lc rgb '#00aa00' title 'Dependencia Aproximada'\n",
        b4, b2, beta, beta, b2, beta,
        N_particulas, N_particulas,
        b2, N_particulas, N_particulas
    );

    fclose(pfileout);
    system("gnuplot Rg_ajuste_ke.plt");

    printf("Grafica generada: results/Rg_vs_ke.png\n");
}

void Ajuste_Rg_en_b() {

    FILE* pfileout = fopen("Rg_ajuste_b.plt", "w");
    if (!pfileout) {
        perror("Error al crear Rg_ajuste_b.plt");
        return;
    }

    double beta = 1.0 / KbT;

    fprintf(pfileout,
        "set terminal pngcairo size 1350,900 enhanced font 'Arial,14'\n"
        "set output 'results/Rg_vs_b.png'\n"
        "set title 'Radio de giro vs b'\n"
        "set xlabel 'b'\n"
        "set ylabel 'Rg^2'\n"
        "set grid\n"
        "set key top left\n"
        "\n"
        "f_ell2(x) = (x**4 + 6*x*x/(%f*%f) + 3/( (%f*%f)**2 )) / (x*x + 1/(%f*%f))\n"
        "f_ajuste(x) = f_ell2(x) * (%d - 1.0/%d) / 6\n"
        "f_teorica(x) = (x*x) * (%d - 1.0/%d) / 6\n"
        "\n"
        "plot 'results/Rg_vs_b.txt' using 1:2 with points pt 2 ps 3 lw 3 lc rgb '#dd181f' title '<Rg^2>', \\\n"
        "     'results/Rg_vs_b.txt' using 1:3 with points pt 1 ps 3 lw 3 lc rgb '#0060ad' title '<Rg>^2', \\\n"
        "     f_ajuste(x) with lines lw 2 lc rgb '#0060ad' title 'Dependencia Teorica', \\\n"
        "     f_teorica(x) with lines lw 2 lc rgb '#00aa00' title 'Dependencia Aproximada'\n",
        beta, k, beta, k, beta, k,
        N_particulas, N_particulas,
        N_particulas, N_particulas
    );

    fclose(pfileout);
    system("gnuplot Rg_ajuste_b.plt");

    printf("Grafica generada: results/Rg_vs_b.png\n");
}

void crear_script_vmd(int N_bloques) {
    FILE* f = fopen("ver_polimero.vmd", "w");
    if (!f) {
        perror("Error al crear ver_polimero.vmd");
        exit(1);
    }

    // Cargar archivos XYZ
    fprintf(f, "mol new results/Data_00001.xyz type xyz first 0 last 0 step 1 waitfor all\n");
    for (int i = 2; i <= N_bloques; i++)
        fprintf(f, "mol addfile results/Data_%05d.xyz type xyz first 0 last -1 step 1 waitfor all 0\n", i);

    // Crear enlaces consecutivos
    fprintf(f, "color Display Background white\n");
    fprintf(f, "topo clearbonds\n");
    fprintf(f, "for {set i 0} {$i < %d} {incr i} {\n", N_particulas - 1);
    fprintf(f, "    topo addbond $i [expr {$i+1}]\n");
    fprintf(f, "}\n");

    fprintf(f, "mol delrep 0 0\n");
    fprintf(f, "mol delrep 0 0\n");
    // Licorice para todos los enlaces
    fprintf(f, "mol representation Licorice 0.05 0.2 12.0\n");
    fprintf(f, "mol color ColorID 1\n");
    fprintf(f, "mol selection all\n");
    fprintf(f, "mol addrep 0\n");

    // CPK para partículas tipo R
    fprintf(f, "mol representation CPK 0.7 0.2\n");
    fprintf(f, "mol color ColorID 1\n");
    fprintf(f, "mol selection {type R}\n");
    fprintf(f, "mol addrep 0\n");

    // CPK para partículas tipo B
    fprintf(f, "mol representation CPK 0.7 0.2\n");
    fprintf(f, "mol color ColorID 0\n");
    fprintf(f, "mol selection {type B}\n");
    fprintf(f, "mol addrep 0\n");

    fprintf(f, "set sel [atomselect top all]\n");
    fprintf(f, "set com [measure center $sel]\n");
    fprintf(f, "$sel moveby [vecinvert $com]\n");
    fprintf(f, "$sel delete\n");

    // -------------------------------
    // Animación
    // -------------------------------
    fprintf(f, "animate goto 0\n");
    fprintf(f, "animate style Loop\n");
    fprintf(f, "animate speed 1.0\n");
    fprintf(f, "animate forward\n");

    fclose(f);
    printf("Script VMD 'ver_polimero.vmd' creado correctamente.\n");
}








