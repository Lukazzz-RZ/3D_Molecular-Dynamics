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