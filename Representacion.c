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
        fprintf(f, "P %.6f %.6f %.6f\n",
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

void crear_script_vmd(int N_bloques) {
    FILE* f = fopen("ver_polimero.vmd", "w");
    if (!f) {
        perror("Error al crear ver_polimero.vmd");
        exit(1);
    }

    // Cargar el primer archivo
    fprintf(f, "mol new results/Data_00001.xyz type xyz first 0 last 0 step 1 waitfor all\n");

    // Añadir los demás bloques
    for (int i = 2; i <= N_bloques; i++) {
        fprintf(f, "mol addfile results/Data_%05d.xyz type xyz first 0 last -1 step 1 waitfor all 0\n", i);
    }

    // Configuración visual
    fprintf(f, "\ndisplay resetview\n");
    fprintf(f, "display resize 1000 800\n");
    fprintf(f, "display projection Perspective\n");
    fprintf(f, "color Display Background white\n");
    fprintf(f, "axes location off\n");
    fprintf(f, "stage location off\n");

    fprintf(f, "\nmol delrep 0 0\n");
    fprintf(f, "mol representation Licorice 0.4 12 12\n");
    fprintf(f, "mol color ColorID 1\n");
    fprintf(f, "mol addrep 0\n");

    // Centrar molécula
    fprintf(f, "\nset sel [atomselect top all]\n");
    fprintf(f, "set com [measure center $sel]\n");
    fprintf(f, "$sel moveby [vecinvert $com]\n");
    fprintf(f, "$sel delete\n");

    // Animación
    fprintf(f, "\nanimate goto 0\n");
    fprintf(f, "animate style Loop\n");
    fprintf(f, "animate speed 1.0\n");
    fprintf(f, "animate forward\n");

    fclose(f);
    printf("Script VMD 'ver_polimero.vmd' creado correctamente.\n");
}

