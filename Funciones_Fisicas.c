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

double rc;
double sigma;
double eps;


double theta_0;
double kb;

char CMode;

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
        #ifdef ALPHATEST
            V_externo += CoulombV(P1);
        #endif

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
    //El 2 ya que hay que tener en cuenta la vuelta
    double V_theta = 2*kb*(1-Pesc(resta(P.pos,P_ant.pos),resta(P_sig.pos,P.pos)));
    V_externo +=V_LJ + V_theta; // Compensar el factor 1/2
    
        #ifdef ALPHATEST
            V_externo += CoulombV(P);
        #endif

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

double CoulombV(Particula Pi){
    double V_Aux = 0;
    double r;
    if (Pi.q==0.) return 0;
    else {
        for (int j=0; j<N_particulas; j++){
        r=modulo(resta(Pi.pos,P[j].pos));

            if(r !=0 && P[j].q !=0) //Es decir, no es la misma partícula y la carga es distinta de 0
            V_Aux += P[j].q*Pi.q/r;
    }
    return V_Aux;


}
}
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

    Vector F_LJ = Fuerza_LennardJones(P1);
    F.x += F_LJ.x;
    F.y += F_LJ.y;
    F.z += F_LJ.z;

        #ifdef ALPHATEST
        Vector F_Cl = Fuerza_CoulombV(P1);
        F.x += F_Cl.x;
        F.y += F_Cl.y;
        F.z += F_Cl.z;
        #endif

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

	Vector F_LJ = Fuerza_LennardJones(P);
    Vector v1 = resta(Pant.pos, P.pos);
    Vector vBF = Fuerza_Bending(Pant, P, Psig);

    F.x += F_LJ.x + vBF.x;
	F.y += F_LJ.y + vBF.y;
	F.z += F_LJ.z + vBF.z;

        #ifdef ALPHATEST
        Vector F_Cl = Fuerza_CoulombV(P);
        F.x += F_Cl.x;
        F.y += F_Cl.y;
        F.z += F_Cl.z;
        #endif

    #endif

    return F;
    }

Vector Fuerza_LennardJones(Particula pi) {
    Vector F = { 0.0f, 0.0f, 0.0f };
    double sigma6 = pow(sigma, 6);
    double sigma12 = sigma6 * sigma6;

    for (int j = 0; j < N_particulas; j++) {
        if (&P[j] == &pi) continue;

        // vector separacion
        Vector rij = resta(pi.pos, P[j].pos);
        double r = modulo(rij);

        if (r == 0 || r > rc) continue;

        double inv_r2 = 1.0 / (r * r);
        double inv_r6 = inv_r2 * inv_r2 * inv_r2;
        double inv_r12 = inv_r6 * inv_r6;

        double Fscalar = 24 * eps * (2 * sigma12 * inv_r12 - sigma6 * inv_r6) * inv_r2;
        F.x += Fscalar * rij.x;
        F.y += Fscalar * rij.y;
        F.z += Fscalar * rij.z;
    }

    return F;
}

Vector Fuerza_CoulombV(Particula Pi){
    Vector F = { 0.0f, 0.0f, 0.0f };

    if (Pi.q==0) return F;
    double Vi = CoulombV(Pi);


    for (int j = 0; j < N_particulas; j++) {
        if (&P[j] == &Pi) continue;

        // vector separacion
        Vector rij = resta(Pi.pos, P[j].pos);
        double r = modulo(rij);

        if (r == 0 || r > b*8 || P[j].q ==0) continue;

        F.x += Vi*rij.x/r/r;
        F.y += Vi*rij.x/r/r;
        F.z += Vi*rij.x/r/r;

}
}

Vector Fuerza_Bending(Particula Pant, Particula P, Particula Psig){

    Vector F;
        F.x=0;
        F.y=0;
        F.z=0;

    Vector vP_Ant_1 = resta(P.pos,Pant.pos);
    Vector vSig_P_1 = resta (Psig.pos,P.pos);
    
;
    //Ten en mente luego copiar todos intercambiarlsos y ya que te ahorras pensar

    //VUELTA 1
    
    //Redefinir en vuelta 2
    double theta_1 = acos(Pesc(vP_Ant_1,vSig_P_1));
    double Rot_cte_1 = cos (theta_1 - theta_0)/cos(theta_1);
    double r_1 = modulo(vSig_P_1);
    Vector z_new_1 = Normalizador(vP_Ant_1);
    Vector a_1 = z_new_1;
        if (a_1.x !=0 ) a_1.x = -a_1.x;
        else a_1.y = -a_1.y;
    Vector x_new_1 = Normalizador (resta (a_1, Scalar_mult(z_new_1, Pesc(a_1,z_new_1))));
    Vector y_new_1 = Vprod(x_new_1,z_new_1);

    Vector Fx_new_1 = Scalar_mult(x_new_1, kb*Rot_cte_1*vSig_P_1.z/r_1/r_1/r_1*vSig_P_1.z);
    Vector Fy_new_1 = Scalar_mult(y_new_1, kb*Rot_cte_1*vSig_P_1.z/r_1/r_1/r_1*vSig_P_1.y);
    Vector Fz_new_1 = Scalar_mult(z_new_1, kb*Rot_cte_1*(vSig_P_1.y*vSig_P_1.y+vSig_P_1.x*vSig_P_1.x)/r_1/r_1/r_1);
    
    F.x += Fx_new_1.x + Fy_new_1.x + Fz_new_1.x;
    F.y += Fx_new_1.y + Fy_new_1.y + Fz_new_1.y;
    F.x += Fx_new_1.z + Fy_new_1.z + Fz_new_1.z;
    
    //Ya tenemos los vectores de fuerza en la base que toca
    
    //VUELTA 2
    //Mismo troncho que antes, solo que en sentido contrario
    
    Vector vP_Ant_2 =  resta (Psig.pos, P.pos);
    Vector vSig_P_2 = resta(P.pos, Pant.pos);

    vP_Ant_2 = Scalar_mult(vP_Ant_2, -1.);
    vSig_P_2 = Scalar_mult(vSig_P_2, -1.);

    double theta_2 = acos(Pesc(vP_Ant_2,vSig_P_2));

    double c = cos(theta_2);
    if (fabs(c) < 0.01) c = (c >= 0 ? 0.01 : -0.01);
    double Rot_cte_2 = cos(theta_2 - theta_0) / c;


    double r_2 = modulo(vSig_P_2);
    Vector z_new_2 = Normalizador(vP_Ant_2);

        Vector pa_1;
            pa_1.x = 0.;
            pa_1.y = 0.;
            pa_1.z = 1.;

        Vector pa_2;
            pa_2.x = 0.;
            pa_2.y = 1.;
            pa_2.z = 0.;

    Vector a_2 = (fabs(z_new_2.x) < 0.9 ? pa_1 : pa_2);//Podría haber 2 v paralelos?

    Vector x_new_2 = Normalizador (resta (a_2, Scalar_mult(z_new_2, Pesc(a_2,z_new_2))));
    Vector y_new_2 = Vprod(x_new_2,z_new_2);

    Vector Fx_new_2 = Scalar_mult(x_new_2, kb*Rot_cte_2*vSig_P_2.z/r_2/r_2/r_2*vSig_P_2.z);
    Vector Fy_new_2 = Scalar_mult(y_new_2, kb*Rot_cte_2*vSig_P_2.z/r_2/r_2/r_2*vSig_P_2.y);
    Vector Fz_new_2 = Scalar_mult(z_new_2, kb*Rot_cte_2*(vSig_P_2.y*vSig_P_2.y+vSig_P_2.x*vSig_P_2.x)/r_2/r_2/r_2);

    F.x += Fx_new_2.x + Fy_new_2.x + Fz_new_2.x;
    F.y += Fx_new_2.y + Fy_new_2.y + Fz_new_2.y;
    F.z += Fx_new_2.z + Fy_new_2.z + Fz_new_2.z;

   
    return F;
}
// ACTUALIZACION DE ENERGIAS //

void Actualizar_Energias(Particula* P) {

    for (int j = 0; j < N_particulas; j++) {
        double v2 = modulo(P[j].vel) * modulo(P[j].vel);
        P[j].Ecin = 0.5 * m * v2;
    }
    for (int j = 0; j < N_particulas; j++) { // 0.5 para evitar doble conteo de muelles
        if (j == 0) {
            P[j].Epot = 0.5*Potencial_Extremo(P[j], P[j + 1]);
        }
        else if (j == N_particulas - 1) {
            P[j].Epot = 0.5*Potencial_Extremo(P[j], P[j - 1]);
        }
        else {
            P[j].Epot = 0.5*Potencial_Intermedio(P[j - 1], P[j], P[j + 1]);
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