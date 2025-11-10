#include "head.h"

int main() {
    
    Inicializar();

    int Ndata = (int)(tmax / dt);
    int Nbloques = 10;
    int pasos_por_bloque = Ndata / Nbloques;

    double Ecin_acum = 0.0;
    double Epot_acum = 0.0;
	double Rg_acum = 0.0;
    int nsteps = 0;

    double buffer[6];

    // BLOQUE DE ACTUALIZACION //
    for (int j = 0; j < Nbloques; j++) {

        FILE* f_xyz = crear_archivo_xyz(j + 1);
        FILE* f_variables = crear_archivo_variables(j + 1, "t\tEcin\tEpot\tRg\tbuf4\tbuf5");

        for (int local_steps = 0; local_steps < pasos_por_bloque; local_steps++) {
            int total_steps = j * pasos_por_bloque + local_steps;
            double t_actual = total_steps * dt;

            // INTEGRADO
            nsteps++;
            verlet_estocastico_3D_extremo(&P[0], P[1], Fuerza_Extremo);
            for (int i = 1; i < N_particulas - 1; i++)
                verlet_estocastico_3D_intermedio(P[i - 1], &P[i], P[i + 1], Fuerza_Intermedio);
            verlet_estocastico_3D_extremo(&P[N_particulas - 1], P[N_particulas - 2], Fuerza_Extremo);
            
            
            // CALCULO ENERGIAS
            Actualizar_Energias(P);
            double Ecin_step = 0.0, Epot_step = 0.0;
            for (int n = 0; n < N_particulas; n++) {
                Ecin_step += P[n].Ecin;
                Epot_step += 0.5 * P[n].Epot; // 0.5 para evitar doble conteo de muelles
            }
            // RADIO DE GIRO
            double Rg = Radio_giro(P);

            
            // ACUMULADOS
			Rg_acum += (Rg - Rg_acum) / nsteps;
            Ecin_acum += (Ecin_step - Ecin_acum) / nsteps;
            Epot_acum += (Epot_step - Epot_acum) / nsteps;

			
			// LLENADO DEL BUFFER
            buffer[0] = t_actual;
            buffer[1] = Ecin_acum;
            buffer[2] = Epot_acum;
            buffer[3] = Rg_acum;
            buffer[4] = 0.0;
            buffer[5] = 0.0;

            // VOLCADO DE DATOS
            guardar_bloque_xyz(P, f_xyz, total_steps);
            Escribir_buffer(f_variables, buffer);
        }
        fclose(f_xyz);
        fclose(f_variables);
        printf("Bloque %d guardado.\n", j + 1);
    }

    // ESTIMADORES DE TERMALIZADO
    printf("Energia cinetica promedio: %.5f\n", Ecin_acum);
    printf("Energia potencial promedio: %.5f\n", Epot_acum);
	printf("Radio de giro promedio: %.5f\n", Rg_acum);

    // DIBUJADO
    //Gnuplot_EnerCons(Nbloques);
	//Gnuplot_Rg(Nbloques);
    /*crear_script_vmd(Nbloques);
    system("vmd -e ver_polimero.vmd");*/

    printf("Simulacion completa.\n");
    return 0;
}