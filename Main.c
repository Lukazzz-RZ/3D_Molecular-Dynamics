#include "head.h"

int main() {

#ifdef COMPLEXMODEL
    printf("Ejecutando en COMPLEXMODEL...\n");
#endif

    Inicializar();

    int Ndata = (int)(tmax / dt);
    int Nbloques = 20;
    int pasos_por_bloque = Ndata / Nbloques;

    double Ecin_acum = 0.0;
    double Epot_acum = 0.0;
    double Rg_acum = 0.0;
    double Rg2_acum = 0.0;
    double L_eff_acum = 0.0;
    long int nsteps = 0;

    double buffer[6];

    // BLOQUE DE ACTUALIZACION //
    for (int j = 0; j < Nbloques; j++) {

        FILE* f_xyz = crear_archivo_xyz(j + 1);
        FILE* f_variables = crear_archivo_variables(j + 1, "t\tEcin\tEpot\tRg\tRg2\tLeff");

        for (int local_steps = 0; local_steps < pasos_por_bloque; local_steps++) {
            int total_steps = j * pasos_por_bloque + local_steps;
            double t_actual = total_steps * dt;

            // INTEGRADO
            nsteps++;
            verlet_estocastico_3D(P);

            // CALCULO ENERGIAS
            Actualizar_Energias(P);
            double Ecin_step = 0.0, Epot_step = 0.0;
            for (int n = 0; n < N_particulas; n++) {
                Ecin_step += P[n].Ecin / (3 * N_particulas); //normalizado por grados de libertad
                Epot_step += P[n].Epot / (N_particulas - 1); //normalizado por numero de muelles
            }
            // RADIO DE GIRO Y LONGITUD EFECTIVA
            double Rg = Radio_giro(P);
            double Rg2 = Rg * Rg;
            double L_eff = Longitud_efectiva(P, 0);

            // ACUMULADOS
            Rg_acum += (Rg - Rg_acum) / nsteps;
            Rg2_acum += (Rg2 - Rg2_acum) / nsteps;
            L_eff_acum += (L_eff - L_eff_acum) / nsteps;
            Ecin_acum += (Ecin_step - Ecin_acum) / nsteps;
            Epot_acum += (Epot_step - Epot_acum) / nsteps;

            // LLENADO DEL BUFFER
            buffer[0] = t_actual;
            buffer[1] = Ecin_acum;
            buffer[2] = Epot_acum;
            buffer[3] = Rg_acum;
            buffer[4] = Rg2_acum;
            buffer[5] = L_eff_acum;

            // VOLCADO DE DATOS
            guardar_bloque_xyz(P, f_xyz, total_steps);
            if (total_steps % 100 == 0) {
                Escribir_buffer(f_variables, buffer);
                //printf("Paso %d: t=%.2f, Ecin=%.5f, Epot=%.5f, Rg=%.5f\n",
                    //total_steps, t_actual, Ecin_acum, Epot_acum, Rg_acum);
            }
        }
        fclose(f_xyz);
        fclose(f_variables);
        printf("Bloque %d guardado.\n", j + 1);
    }


    // ESTIMADORES DE TERMALIZADO
    printf("Energia cinetica promedio: %.5f\n", Ecin_acum);
    printf("Energia potencial promedio: %.5f\n", Epot_acum);
    printf("Longitud efectiva promedio: %.5f\n", L_eff_acum);
    printf("Radio de giro promedio al cuadrado: %.5f\n", Rg_acum * Rg_acum);
    printf("Radio de giro al cuadrado promedio: %.5f\n", Rg2_acum);



    // DIBUJADO
    //Ajuste_Rg_en_b();
    //Ajuste_Rg_en_N();
    //Ajuste_Rg_en_k();
    Ajuste_L_eff_en_Fx();
    //Gnuplot_EnerCons(Nbloques);
    //Gnuplot_Rg(Nbloques);
    //Gnuplot_L_eff(Nbloques);
    crear_script_vmd(Nbloques);
    system("vmd -e ver_polimero.vmd");

    printf("Simulacion completa.\n");
    return 0;
}