#include "head.h"

// DECLARACIÓN DE VARIABLES FISICAS GLOBALES //
Particula P[N_particulas];
double k;
double b;
double m;

double gamma_DP;
double KbT;
Vector Fext;

double dt;
double tmax;

double sigma;
double eps;



// FUNCIONES DE ENERGIA //

double Potencial_Extremo(Particula P1, Particula P2) {

    Vector r_rel = {
        P1.pos.x - P2.pos.x,
        P1.pos.y - P2.pos.y,
        P1.pos.z - P2.pos.z
    };
    double modr_rel = modulo(r_rel);

    double V_elastico = 0.5 * k * (modr_rel - b) * (modr_rel - b);
    double V_externo = -(Fext.x * P1.pos.x + Fext.y * P1.pos.y + Fext.z * P1.pos.z);

    #ifdef COMPLEXMODEL

    double V_LJ= V_LennardJones(P1);
    V_externo +=V_LJ;

    #endif

    return V_elastico + V_externo;
}

double Potencial_Intermedio(Particula P_ant, Particula P, Particula P_sig) {
    Vector r1 = {
        P.pos.x - P_ant.pos.x,
        P.pos.y - P_ant.pos.y,
        P.pos.z - P_ant.pos.z
    };
    Vector r2 = {
        P.pos.x - P_sig.pos.x,
        P.pos.y - P_sig.pos.y,
        P.pos.z - P_sig.pos.z
    };

    double modr1 = modulo(r1);
    double modr2 = modulo(r2);

    double V_elastico = 0.5 * k * (modr1 - b) * (modr1 - b)
        + 0.5 * k * (modr2 - b) * (modr2 - b);

    double V_externo = -(Fext.x * P.pos.x + Fext.y * P.pos.y + Fext.z * P.pos.z);

    #ifdef COMPLEXMODEL

    double V_LJ= V_LennardJones(P);
    V_externo +=V_LJ;

    #endif

    return V_elastico + V_externo;
}

double V_LennardJones(Particula pi){
    //Potencial que observa la partícula pi por el resto de las pj
    double V_Aux = 0;
    double r;
    for (int j=0; j<N_particulas; j++){
        r=modulo(resta(pi.pos,P[j].pos));
            if(r !=0 ) //Es decir, no es la misma partícula
            V_Aux += 4*eps*(pow(sigma/r,12)-pow(sigma/r,6));
    }
    return V_Aux;
}

// DEFINICIÓN DEL POTENCIAL //

Vector Fuerza_Extremo(Particula P1, Particula P2) {
    Vector r_rel = {
        P1.pos.x - P2.pos.x,
        P1.pos.y - P2.pos.y,
        P1.pos.z - P2.pos.z
    };
    double modr_rel = modulo(r_rel);

    Vector F;

    if (modr_rel <= 1e-10) {
        F = Fext;
    }
    else {
        double factor = -k * (modr_rel - b) / modr_rel;
        F.x = factor * r_rel.x + Fext.x;
        F.y = factor * r_rel.y + Fext.y;
        F.z = factor * r_rel.z + Fext.z;
    }

    #ifdef COMPLEXMODEL

    F.x += dV_LennardJones(P1,1);
    F.y += dV_LennardJones(P1,2);
    F.z += dV_LennardJones(P1,3);

    #endif

    return F;
}

Vector Fuerza_Intermedio(Particula Pant, Particula P, Particula Psig) {

    Vector r1 = {
        P.pos.x - Pant.pos.x,
        P.pos.y - Pant.pos.y,
        P.pos.z - Pant.pos.z
    };
    double modr1 = modulo(r1);

    Vector r2 = {
        P.pos.x - Psig.pos.x,
        P.pos.y - Psig.pos.y,
        P.pos.z - Psig.pos.z
    };
    double modr2 = modulo(r2);

    Vector F;
    F.x = F.y = F.z = 0.0;

    if (modr1 > 1e-10) {
        double factor1 = -k * (modr1 - b) / modr1;
        F.x += factor1 * r1.x;
        F.y += factor1 * r1.y;
        F.z += factor1 * r1.z;
    }
    if (modr2 > 1e-10) {
        double factor2 = -k * (modr2 - b) / modr2;
        F.x += factor2 * r2.x;
        F.y += factor2 * r2.y;
        F.z += factor2 * r2.z;
    }

    F.x += Fext.x;
    F.y += Fext.y;
    F.z += Fext.z;

    #ifdef COMPLEXMODEL

    F.x += dV_LennardJones(P,1);
    F.y += dV_LennardJones(P,2);
    F.z += dV_LennardJones(P,3);

    #endif

    return F;
    }

    double dV_LennardJones(Particula pi, int Coor){
    //Potencial que observa la partícula pi por el resto de las pj
    //Si, he reutilizado el otro descaradamente
    double V_Aux = 0;
    double r;
    int part_Nominal;
    for (int j=0; j<N_particulas; j++){
        r=modulo(resta(pi.pos,P[j].pos));
            if(r !=0 ) //Es decir, no es la misma partícula
            V_Aux += 1/r/r/4*eps*(12*pow(sigma/r,12)-6*pow(sigma/r,6));
            else part_Nominal = j;
    }
    if (Coor == 1) return V_Aux*P[part_Nominal].pos.x;
    if (Coor == 2) return V_Aux*P[part_Nominal].pos.y;
    if (Coor == 3) return V_Aux*P[part_Nominal].pos.z;

}

// ACTUALIZACION DE ENERGIAS //

void Actualizar_Energias(Particula* P) {

    for (int j = 0; j < N_particulas; j++) {
        double v2 = modulo(P[j].vel) * modulo(P[j].vel);
        P[j].Ecin = 0.5 * m * v2;
    }
    for (int j = 0; j < N_particulas; j++) {
        if (j == 0) {
            P[j].Epot = Potencial_Extremo(P[j], P[j + 1]);
        }
        else if (j == N_particulas - 1) {
            P[j].Epot = Potencial_Extremo(P[j], P[j - 1]);
        }
        else {
            P[j].Epot = Potencial_Intermedio(P[j - 1], P[j], P[j + 1]);
        }
    }
}

// MAGNITUDES DEL POLIMERO

double Radio_giro(Particula* P) {

    Vector posiciones[N_particulas];
    for (int i = 0; i < N_particulas; i++) {
        posiciones[i] = P[i].pos;
	}

    Vector r_cdm = CDM_uniforme(posiciones, N_particulas);
    double suma = 0.0;
    for (int i = 0; i < N_particulas; i++) {
        Vector diff = {
            posiciones[i].x - r_cdm.x,
            posiciones[i].y - r_cdm.y,
            posiciones[i].z - r_cdm.z
        };
        double mod_diff = modulo(diff);
        suma += mod_diff * mod_diff;
    }
    double Rg_squared = suma / N_particulas;
    return sqrt(Rg_squared);
    
}