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

double kb;
double theta_0;
double Fold_Ct;
double Q_ct;

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

    double V_LJ= Potencial_LennardJones(P1);
    V_externo += 2 * V_LJ;

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

    double V_LJ= Potencial_LennardJones(P);
    V_externo +=2*V_LJ; // Compensar el factor 1/2

    #endif

    return V_elastico + V_externo;
}

double Potencial_LennardJones(Particula pi){
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
/*
double Vdih(Particula mm, Particula m, Particula i, Particula p, Particula pp){
    #ifdef DIHEDRICAL
        double C = 1.2;
        double D = 1.2;
        double E = 1.2;
        double Dih_1 = Phi_Dd(m, i, p, pp);
        double Dih_2 = Phi_Dd(p, i, m, mm);
    return  (C*(1-cos(Dih_1))+D*(1-cos(3*Dih_1))-E*(1+cos(Dih_1 +PI/4))) + //Vuelta 1
                (C*(1-cos(Dih_2))+D*(1-cos(3*Dih_2))-E*(1+cos(Dih_2 +PI/4))); //Vuelta 2
    #else

    return 0.;

    #endif}
*/

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


// FUNCIONES DE FUERZA //

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

    Vector F_Cl = Fuerza_CoulombV(P1);
    F.x += F_Cl.x;
    F.y += F_Cl.y;
    F.z += F_Cl.z;

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
    F.x += F_LJ.x;
	F.y += F_LJ.y;
	F.z += F_LJ.z;
    
    Vector F_Cl = Fuerza_CoulombV(P);
    F.x += F_Cl.x;
    F.y += F_Cl.y;
    F.z += F_Cl.z;
    #endif

    #ifdef ALPHATEST

    Vector F_Bd = Fuerza_Bending(Pant, P, Psig);
    F.x += F_Bd.x;
    F.y += F_Bd.y;
    F.z += F_Bd.z;
    #endif

    return F;
    }

// TERMINOS EXTRA DE FUERZAS

Vector Fuerza_LennardJones(Particula pi) {
    Vector F = { 0.0, 0.0, 0.0 };
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
    double epsilon = 10e-8;
    Vector F;
        F.x=0;
        F.y=0;
        F.z=0;

            Vector pa_1;
            pa_1.x = 0.;
            pa_1.y = 0.;
            pa_1.z = 1.;

        Vector pa_2;
            pa_2.x = 0.;
            pa_2.y = 1.;
            pa_2.z = 0.;

    Vector vP_Ant_1 = resta(P.pos,Pant.pos);
    Vector vSig_P_1 = resta (Psig.pos,P.pos);
    
;
    //Ten en mente luego copiar todos intercambiarlsos y ya que te ahorras pensar

    //VUELTA 1
    
    //Redefinir en vuelta 2
    double theta_1 = acos(Pesc(vP_Ant_1,vSig_P_1));

    double c_1 = cos(theta_1-theta_0);
    if (fabs(c_1) < epsilon) c_1 = (c_1 >= 0 ? epsilon : -epsilon);
    double Rot_cte_1 = cos(theta_1)/c_1;

    double r_1 = modulo(vSig_P_1);
    Vector z_new_1 = Normalizador(vP_Ant_1);
    Vector a_1 = (fabs(z_new_1.x) < 0.9 ? pa_1 : pa_2);

    Vector x_new_1 = Normalizador (resta (a_1, Scalar_mult(z_new_1, Pesc(a_1,z_new_1))));
    Vector y_new_1 = Vprod(x_new_1,z_new_1);

    Vector Fx_new_1 = Scalar_mult(x_new_1, kb*Rot_cte_1*Pesc_NoNorm(vSig_P_1,z_new_1)/r_1/r_1/r_1*Pesc_NoNorm(vSig_P_1, x_new_1));
    Vector Fy_new_1 = Scalar_mult(y_new_1, kb*Rot_cte_1*Pesc_NoNorm(vSig_P_1, z_new_1)/r_1/r_1/r_1*Pesc_NoNorm(vSig_P_1, y_new_1));
    Vector Fz_new_1 = Scalar_mult(z_new_1, kb*Rot_cte_1*(Pesc_NoNorm(vSig_P_1, y_new_1)*Pesc_NoNorm(vSig_P_1, y_new_1)+Pesc_NoNorm(vSig_P_1, x_new_1)*Pesc_NoNorm(vSig_P_1, x_new_1))/r_1/r_1/r_1);
    
    F.x += Fx_new_1.x + Fy_new_1.x + Fz_new_1.x;
    F.y += Fx_new_1.y + Fy_new_1.y + Fz_new_1.y;
    F.z += Fx_new_1.z + Fy_new_1.z + Fz_new_1.z;
    
    //Ya tenemos los vectores de fuerza en la base que toca
    
    //VUELTA 2
    //Mismo troncho que antes, solo que en sentido contrario
    
    Vector vP_Ant_2 =  resta (Psig.pos, P.pos);
    Vector vSig_P_2 = resta(P.pos, Pant.pos);

    vP_Ant_2 = Scalar_mult(vP_Ant_2, -1.);
    vSig_P_2 = Scalar_mult(vSig_P_2, -1.);

    double theta_2 = acos(Pesc(vP_Ant_2,vSig_P_2));

    double c_2 = cos(theta_2-theta_0);
    if (fabs(c_2) < epsilon) c_2 = (c_2 >= 0 ? epsilon : -epsilon);
    double Rot_cte_2 = cos(theta_2)/c_2;


    double r_2 = modulo(vSig_P_2);
    Vector z_new_2 = Normalizador(vP_Ant_2);

    Vector a_2 = (fabs(z_new_2.x) < 0.9 ? pa_1 : pa_2);//Podría haber 2 v paralelos?

    Vector x_new_2 = Normalizador (resta (a_2, Scalar_mult(z_new_2, Pesc(a_2,z_new_2))));
    Vector y_new_2 = Vprod(x_new_2,z_new_2);

    Vector Fx_new_2 = Scalar_mult(x_new_2, kb*Rot_cte_2*Pesc_NoNorm(vSig_P_2,z_new_2)/r_2/r_2/r_2*Pesc_NoNorm(vSig_P_2, x_new_2));
    Vector Fy_new_2 = Scalar_mult(y_new_2, kb*Rot_cte_2*Pesc_NoNorm(vSig_P_2, z_new_2)/r_2/r_2/r_2*Pesc_NoNorm(vSig_P_2, y_new_2));
    Vector Fz_new_2 = Scalar_mult(z_new_2, kb*Rot_cte_2*(Pesc_NoNorm(vSig_P_2, y_new_2)*Pesc_NoNorm(vSig_P_2, y_new_2)+Pesc_NoNorm(vSig_P_2, x_new_2)*Pesc_NoNorm(vSig_P_2, x_new_2))/r_2/r_2/r_2);

    F.x += Fx_new_2.x + Fy_new_2.x + Fz_new_2.x;
    F.y += Fx_new_2.y + Fy_new_2.y + Fz_new_2.y;
    F.z += Fx_new_2.z + Fy_new_2.z + Fz_new_2.z;

   
    return F;
}

Vector F_dih (Particula Pmm, Particula Pm, Particula Pi, Particula Pp, Particula Ppp, double LimIz, double LimDer){

    double eps = 1e-12;

    Vector F_Aux = {0., 0., 0.};

#ifdef DIHEDRICAL

    // vectores auxiliares para construir la base local
    Vector pa_1 = {1., 0., 0.};
    Vector pa_2 = {0., 1., 0.};

    //----------------------------------------------
    //                VUELTA 1
    //----------------------------------------------

    // Definimos b1, b2, b3 (primer grupo  i-1, i, i+1, i+2)
    Vector b1 = resta(Pi.pos, Pm.pos);      // r_i - r_{i-1}
    Vector b2 = resta(Pp.pos, Pi.pos);      // r_{i+1} - r_i
    Vector b3 = resta(Ppp.pos, Pp.pos);     // r_{i+2} - r_{i+1}

    // Normales correctas: n1 = b1 × b2,   n2 = b2 × b3
    Vector n1 = Normalizador((b1, b2));
    Vector n2 = Normalizador((b2, b3));

    double B = modulo(n1);
    double C = modulo(n2);

    if (B >= eps && C >= eps)  // si no, fuerza ≈ 0 (dihedral colapsado)
    {
        // Ejes locales
        Vector z_new = Normalizador(n1);

        // Elegir vector auxiliar "a" NO colineal con z_new
        Vector a = (fabs(z_new.x) < 0.9 ? pa_1 : pa_2);

        // Construcción de x_new ortogonal a z_new
        Vector tmp = resta(a, Scalar_mult(z_new, Pesc(a, z_new)));
        double tmp_norm = modulo(tmp);
        if (tmp_norm < eps) tmp_norm = eps;
        Vector x_new = Scalar_mult(tmp, 1.0/tmp_norm);

        Vector y_new = Vprod(x_new, z_new);

        double r = C;   // ||n2||
        if (r < eps) r = eps;

        // Componentes en la base local
        double c_z = Pesc(n2, z_new);
        double c_x = Pesc(n2, x_new);
        double c_y = Pesc(n2, y_new);

        // Construcción de la fuerza local
        double inv_r3 = 1.0 / (r*r*r);

        Vector Fx_new = Scalar_mult(x_new, c_z * c_x * inv_r3);
        Vector Fy_new = Scalar_mult(y_new, c_z * c_y * inv_r3);
        Vector Fz_new = Scalar_mult(z_new, (c_x*c_x + c_y*c_y) * inv_r3);

        // Transformar y sumar
        F_Aux.x += (Fx_new.x + Fy_new.x + Fz_new.x) * LimIz;
        F_Aux.y += (Fx_new.y + Fy_new.y + Fz_new.y) * LimIz;
        F_Aux.z += (Fx_new.z + Fy_new.z + Fz_new.z) * LimIz;
    }


    //----------------------------------------------
    //                VUELTA 2
    //----------------------------------------------
    //
    // La vuelta 2 debe ser *simétrica* a la primera,
    // usando el cuarteto (i-2, i-1, i, i+1)
    //----------------------------------------------

    b1 = resta(Pm.pos, Pmm.pos);     // r_{i-1} - r_{i-2}
    b2 = resta(Pi.pos, Pm.pos);      // r_i - r_{i-1}
    b3 = resta(Pp.pos, Pi.pos);      // r_{i+1} - r_i

    n1 = Normalizador(Vprod(b1, b2));
    n2 = Normalizador(Vprod(b2, b3));

    B = modulo(n1);
    C = modulo(n2);

    if (B >= eps && C >= eps)
    {
        Vector z_new = Normalizador(n1);

        Vector a = (fabs(z_new.x) < 0.9 ? pa_1 : pa_2);

        Vector tmp = resta(a, Scalar_mult(z_new, Pesc(a, z_new)));
        double tmp_norm = modulo(tmp);
        if (tmp_norm < eps) tmp_norm = eps;
        Vector x_new = Scalar_mult(tmp, 1.0/tmp_norm);

        Vector y_new = Vprod(x_new, z_new);

        double r = C;
        if (r < eps) r = eps;

        double c_z = Pesc(n2, z_new);
        double c_x = Pesc(n2, x_new);
        double c_y = Pesc(n2, y_new);

        double inv_r3 = 1.0 / (r*r*r);

        Vector Fx_new = Scalar_mult(x_new, c_z * c_x * inv_r3);
        Vector Fy_new = Scalar_mult(y_new, c_z * c_y * inv_r3);
        Vector Fz_new = Scalar_mult(z_new, (c_x*c_x + c_y*c_y) * inv_r3);

        F_Aux.x += (Fx_new.x + Fy_new.x + Fz_new.x) * LimDer;
        F_Aux.y += (Fx_new.y + Fy_new.y + Fz_new.y) * LimDer;
        F_Aux.z += (Fx_new.z + Fy_new.z + Fz_new.z) * LimDer;
    }

#endif

    return F_Aux;

}

// CALCULOS DINAMICOENERGETICOS DE POLIMERO //

void Calcular_fuerzas(Vector* F, Particula* P) {

    for (int i = 0; i < N_particulas; i++) {
        F[i].x = 0.0;
        F[i].y = 0.0;
        F[i].z = 0.0;
    }

    if (N_particulas > 1) {
        Vector Fe0 = Fuerza_Extremo(P[0], P[1]);
        F[0].x += Fe0.x;   F[0].y += Fe0.y;   F[0].z += Fe0.z;
    }
    if (N_particulas > 1) {
        Vector FeN = Fuerza_Extremo(P[N_particulas - 1], P[N_particulas - 2]);
        F[N_particulas - 1].x += FeN.x; F[N_particulas - 1].y += FeN.y; F[N_particulas - 1].z += FeN.z;
    }

    for (int i = 1; i < N_particulas - 1; i++) {

        Vector Fi = Fuerza_Intermedio(P[i - 1], P[i], P[i + 1]);
            #ifdef DIHEDRICAL

        Vector Dih_Force;
        if ( i == 1 ) Dih_Force = F_dih( P[i-1], P[i-1], P[i], P[i+1], P[i+2], 0., Fold_Ct);
        if ( i == N_particulas-2) Dih_Force = F_dih(P[i-2], P[i-1], P[i], P[i+1], P[i+1], Fold_Ct, 0.);
        else Dih_Force = F_dih(P[i-2], P[i-1], P[i], P[i+1], P[i+2], Fold_Ct, Fold_Ct);

        Fi.x += Dih_Force.x;
        Fi.y += Dih_Force.y;
        Fi.z += Dih_Force.z;

            #endif

        F[i].x += Fi.x;
        F[i].y += Fi.y;
        F[i].z += Fi.z;
    }
}

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

double Longitud_efectiva(Particula* P, int norm) {
	double L = (resta(P[N_particulas - 1].pos, P[0].pos)).x; // PREPARADO PARA FUERZAS EN X, PENDIENTE DE GENERALIZAR
    if (norm == 0) {
        return L/(N_particulas-1)/b;
	}
    else return L;
}