#include "head.h"

// METODOS NUMERICOS PARA ECUACIONES DIFERENCIALES //

void verlet_estocastico_3D_extremo(Particula* P, Particula P2, FuncionFuerzaExtremo Fuerza) {

    double b_coef = 1.0 / (1.0 + gamma_DP * dt / (2.0 * m));
    double a_coef = (1.0 - gamma_DP * dt / (2.0 * m)) / (1.0 + gamma_DP * dt / (2.0 * m));
    Vector F0 = Fuerza(*P, P2);
    double prefactor = sqrt(2.0 * gamma_DP * KbT * dt / m);
    Vector ruido = {
        prefactor * N_Rand_Gauss(),
        prefactor * N_Rand_Gauss(),
        prefactor * N_Rand_Gauss()
    };

    P->pos.x += b_coef * dt * P->vel.x + b_coef * dt * dt / (2.0 * m) * F0.x + b_coef * dt / (2.0 * m) * ruido.x;
    P->pos.y += b_coef * dt * P->vel.y + b_coef * dt * dt / (2.0 * m) * F0.y + b_coef * dt / (2.0 * m) * ruido.y;
    P->pos.z += b_coef * dt * P->vel.z + b_coef * dt * dt / (2.0 * m) * F0.z + b_coef * dt / (2.0 * m) * ruido.z;

    Vector F1 = Fuerza(*P ,P2);

    P->vel.x = P->vel.x * a_coef + dt / (2.0 * m) * (a_coef * F0.x + F1.x) + b_coef / m * ruido.x;
    P->vel.y = P->vel.y * a_coef + dt / (2.0 * m) * (a_coef * F0.y + F1.y) + b_coef / m * ruido.y;
    P->vel.z = P->vel.z * a_coef + dt / (2.0 * m) * (a_coef * F0.z + F1.z) + b_coef / m * ruido.z;
}

void verlet_estocastico_3D_intermedio(Particula P2, Particula* P, Particula P3, FuncionFuerzaIntermedio Fuerza) {

    double b_coef = 1.0 / (1.0 + gamma_DP * dt / (2.0 * m));
    double a_coef = (1.0 - gamma_DP * dt / (2.0 * m)) / (1.0 + gamma_DP * dt / (2.0 * m));
    Vector F0 = Fuerza(P2, *P, P3);
    double prefactor = sqrt(2.0 * gamma_DP * KbT * dt / m);
    Vector ruido = {
        prefactor * N_Rand_Gauss(),
        prefactor * N_Rand_Gauss(),
        prefactor * N_Rand_Gauss()
    };

    P->pos.x += b_coef * dt * P->vel.x + b_coef * dt * dt / (2.0 * m) * F0.x + b_coef * dt / (2.0 * m) * ruido.x;
    P->pos.y += b_coef * dt * P->vel.y + b_coef * dt * dt / (2.0 * m) * F0.y + b_coef * dt / (2.0 * m) * ruido.y;
    P->pos.z += b_coef * dt * P->vel.z + b_coef * dt * dt / (2.0 * m) * F0.z + b_coef * dt / (2.0 * m) * ruido.z;

    Vector F1 = Fuerza(P2, *P, P3);

    P->vel.x = P->vel.x * a_coef + dt / (2.0 * m) * (a_coef * F0.x + F1.x) + b_coef / m * ruido.x;
    P->vel.y = P->vel.y * a_coef + dt / (2.0 * m) * (a_coef * F0.y + F1.y) + b_coef / m * ruido.y;
    P->vel.z = P->vel.z * a_coef + dt / (2.0 * m) * (a_coef * F0.z + F1.z) + b_coef / m * ruido.z;
}

void Polimer_Updater(Particula *P_new){
    Particula P_old[N_particulas];
    Copia_Polimero(P_new, P_old);

    //Polimeros idénticos. Basta con pasar en pos no actualizadas el _old y en las que se actualizen el _new
    
    //Actualizamos polímero
    verlet_estocastico_3D_extremo(&P_new[0], P_old[1], Fuerza_Extremo);
         for (int i = 1; i < N_particulas - 1; i++)
            verlet_estocastico_3D_intermedio(P_old[i - 1], &P_new[i], P_old[i + 1], Fuerza_Intermedio);
    verlet_estocastico_3D_extremo(&P_new[N_particulas - 1], P_old[N_particulas - 2], Fuerza_Extremo);
}