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

void verlet_estocastico_3D(Particula* P) {

    double b = 1.0 / (1.0 + gamma_DP * dt / (2.0 * m));
    double a = (1.0 - gamma_DP * dt / (2.0 * m)) / (1.0 + gamma_DP * dt / (2.0 * m));
    double prefactor = sqrt(2.0 * gamma_DP * KbT * dt / m);

    Vector* F0 = malloc(N_particulas * sizeof(Vector));
    Vector* F1 = malloc(N_particulas * sizeof(Vector));
    Vector* ruido = malloc(N_particulas * sizeof(Vector));

    Vector* delta_pos = malloc(N_particulas * sizeof(Vector));
    Vector* delta_vel = malloc(N_particulas * sizeof(Vector));

    // ---------- 1) Fuerzas iniciales ----------
    Calcular_fuerzas(F0, P);

    // ---------- 2) Ruido ----------
    for (int i = 0; i < N_particulas; i++) {
        ruido[i].x = prefactor * N_Rand_Gauss();
        ruido[i].y = prefactor * N_Rand_Gauss();
        ruido[i].z = prefactor * N_Rand_Gauss();
    }

    // ---------- 3) Calcular dr ----------
    for (int i = 0; i < N_particulas; i++) {
        delta_pos[i].x = b * dt * P[i].vel.x + b * dt * dt / (2 * m) * F0[i].x + b * dt / (2 * m) * ruido[i].x;
        delta_pos[i].y = b * dt * P[i].vel.y + b * dt * dt / (2 * m) * F0[i].y + b * dt / (2 * m) * ruido[i].y;
        delta_pos[i].z = b * dt * P[i].vel.z + b * dt * dt / (2 * m) * F0[i].z + b * dt / (2 * m) * ruido[i].z;
    }

    // ---------- 4) Actualizar posiciones ----------
    for (int i = 0; i < N_particulas; i++) {
        P[i].pos.x += delta_pos[i].x;
        P[i].pos.y += delta_pos[i].y;
        P[i].pos.z += delta_pos[i].z;
    }

    // ---------- 5) Fuerzas nuevas ----------
    Calcular_fuerzas(F1, P);

    // ---------- 6) Calcular dv ----------
    for (int i = 0; i < N_particulas; i++) {
        delta_vel[i].x = a * P[i].vel.x + dt / (2 * m) * (a * F0[i].x + F1[i].x) + b / m * ruido[i].x;
        delta_vel[i].y = a * P[i].vel.y + dt / (2 * m) * (a * F0[i].y + F1[i].y) + b / m * ruido[i].y;
        delta_vel[i].z = a * P[i].vel.z + dt / (2 * m) * (a * F0[i].z + F1[i].z) + b / m * ruido[i].z;
    }

    // ---------- 7) Actualizar velocidades ----------
    for (int i = 0; i < N_particulas; i++) {
        P[i].vel = delta_vel[i];
    }

    free(F0);
    free(F1);
    free(ruido);
    free(delta_pos);
    free(delta_vel);
}