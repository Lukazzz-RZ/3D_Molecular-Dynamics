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

double Q_ct;
double permitivity;

double kb;
double theta_0;

double kb_dih;
double phi_0;
int mult;

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
	double V_Coul = Potencial_Coulomb(P1);
	V_externo += 2 * (V_LJ + V_Coul); // Compensar el factor 1/2

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
    double V_Coul = Potencial_Coulomb(P);
    V_externo += 2 * (V_LJ + V_Coul); // Compensar el factor 1/2

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

double Potencial_Coulomb(Particula Pi) {
    double V_Aux = 0;
    double r;
    if (Pi.q == 0.) return 0;
    else {
        for (int j = 0; j < N_particulas; j++) {
            r = modulo(resta(Pi.pos, P[j].pos));

            if (r != 0 && P[j].q != 0) //Es decir, no es la misma partícula y la carga es distinta de 0
                V_Aux += P[j].q * Pi.q / r;
        }
        return permitivity*V_Aux;

    }
}


////// FUNCIONES DE FUERZA ///////

//   FUERZAS ELÁSTICAS DEL ENLACE

Vector Fuerza_ElasticaExtremo(Particula P1, Particula P2) {

    Vector r_rel = {
        P1.pos.x - P2.pos.x,
        P1.pos.y - P2.pos.y,
        P1.pos.z - P2.pos.z
    };
    double modr_rel = modulo(r_rel);

    Vector F = { 0.0, 0.0, 0.0 };

    if (modr_rel > 1e-10) {
        double factor = -k * (modr_rel - b) / modr_rel;
        F.x = factor * r_rel.x;
        F.y = factor * r_rel.y;
        F.z = factor * r_rel.z;
    }

    return F;
}

Vector Fuerza_ElasticaIntermedio(Particula Pant, Particula P, Particula Psig) {

    Vector F = { 0.0, 0.0, 0.0 };

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

    return F;
}

// BENDING Y TORSION //

void Fuerza_Bending(Particula Pant, Particula P, Particula Psig, Vector F[3])
{
    // ---- Vectores desde el centro ----
    Vector a = resta(Pant.pos, P.pos);   // Pant - P
    Vector b = resta(Psig.pos, P.pos);   // Psig - P

    double a_len = modulo(a);
    double b_len = modulo(b);

    Vector a_hat = Normalize(a);
    Vector b_hat = Normalize(b);

    // ---- Ángulo ----
    double cos_t = Pesc_Norm(a, b);
    if (cos_t > 1.0) cos_t = 1.0;
    if (cos_t < -1.0) cos_t = -1.0;

    double t = acos(cos_t);
    double sin_t = sqrt(1.0 - cos_t * cos_t);

    // ---- derivada ----
    double dtheta = t - theta_0;
    double dVdt = kb * sin(dtheta);

    const double eps = 1e-12;
    Vector F1 = { 0,0,0 };
    Vector F3 = { 0,0,0 };

    if (sin_t > eps) {

        double c1 = -dVdt / (a_len * sin_t);
        double c3 = -dVdt / (b_len * sin_t);

        // Fuerza en Particula anterior (Pant)
        F1.x = c1 * (b_hat.x - cos_t * a_hat.x);
        F1.y = c1 * (b_hat.y - cos_t * a_hat.y);
        F1.z = c1 * (b_hat.z - cos_t * a_hat.z);

        // Fuerza en Particula siguiente (Psig)
        F3.x = c3 * (a_hat.x - cos_t * b_hat.x);
        F3.y = c3 * (a_hat.y - cos_t * b_hat.y);
        F3.z = c3 * (a_hat.z - cos_t * b_hat.z);
    }

    // Fuerza en la partícula central
    Vector F2;
    F2.x = -(F1.x + F3.x);
    F2.y = -(F1.y + F3.y);
    F2.z = -(F1.z + F3.z);

    // Guardar en array de salida
    F[0] = F1;  // fuerza en Pant
    F[1] = F2;  // fuerza en P
    F[2] = F3;  // fuerza en Psig
}

void Fuerza_Dihedral(Particula P1, Particula P2, Particula P3, Particula P4, Vector F[4])
{
    // Vectores de enlaces
    Vector b1 = resta(P2.pos, P1.pos);
    Vector b2 = resta(P3.pos, P2.pos);
    Vector b3 = resta(P4.pos, P3.pos);

    // Producto cruzado de planos
    Vector c1 = Pvect(b1, b2);  // plano P1-P2-P3
    Vector c2 = Pvect(b2, b3);  // plano P2-P3-P4

    double c1_norm2 = Pesc_NoNorm(c1, c1); // |c1|^2
    double c2_norm2 = Pesc_NoNorm(c2, c2); // |c2|^2
    const double eps = 1e-12;

    if (c1_norm2 < eps || c2_norm2 < eps) {
        // fuerzas = 0 si colineales
        for (int i = 0; i < 4; i++) F[i] = (Vector){ 0,0,0 };
        return;
    }

    double cos_phi = Pesc_Norm(c1, c2) / (sqrt(c1_norm2) * sqrt(c2_norm2));
    if (cos_phi > 1.0) cos_phi = 1.0;
    if (cos_phi < -1.0) cos_phi = -1.0;
    double phi = acos(cos_phi);

    // derivada del potencial
    double dVdphi = mult * kb_dih * sin(mult * phi - phi_0); // dV/dphi

    // escalado de fuerzas exacto
    double b2_len2 = Pesc_NoNorm(b2, b2);

    // F1 y F4
    Vector f1 = Scalar_mult(c1, -dVdphi * sqrt(Pesc_NoNorm(b2, b2)) / c1_norm2);
    Vector f4 = Scalar_mult(c2, dVdphi * sqrt(Pesc_NoNorm(b2, b2)) / c2_norm2);

    // F2 y F3: distribuimos la fuerza según la derivada completa
    // Vectores auxiliares
    double b1_dot_b2 = Pesc_NoNorm(b1, b2);
    double b3_dot_b2 = Pesc_NoNorm(b3, b2);

    Vector f2;
    f2.x = -f1.x * (1 - b1_dot_b2 / b2_len2) - f4.x * (b3_dot_b2 / b2_len2);
    f2.y = -f1.y * (1 - b1_dot_b2 / b2_len2) - f4.y * (b3_dot_b2 / b2_len2);
    f2.z = -f1.z * (1 - b1_dot_b2 / b2_len2) - f4.z * (b3_dot_b2 / b2_len2);

    Vector f3;
    f3.x = -(f1.x + f2.x + f4.x);
    f3.y = -(f1.y + f2.y + f4.y);
    f3.z = -(f1.z + f2.z + f4.z);

    // Guardar en array de salida
    F[0] = f1;  // P1
    F[1] = f2;  // P2
    F[2] = f3;  // P3
    F[3] = f4;  // P4
}

// TERMINOS EXTRA DE FUERZAS

Vector Fuerza_LennardJones(Particula pi) {

    Vector F = { 0,0,0 };
    double sigma6 = pow(sigma, 6);
    double sigma12 = sigma6 * sigma6;

    for (int j = 0; j < N_particulas; j++) {

        if (j == pi.idx) continue;   // <-- así evitas self-interaction

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

Vector Fuerza_Coulomb(Particula Pi) {

    Vector F = { 0.0, 0.0, 0.0 };
    if (Pi.q == 0.0) return F;

    for (int j = 0; j < N_particulas; j++) {

        if (j == Pi.idx) continue;  // <-- DE NUEVO, EVITA SELF

        if (P[j].q == 0.0) continue;

        Vector rij = resta(Pi.pos, P[j].pos);
        double r = modulo(rij);
        if (r < 1e-8) continue;

        Vector n = { rij.x / r, rij.y / r, rij.z / r };
        double coef = permitivity * (Pi.q * P[j].q) / (r * r);

        F.x += coef * n.x;
        F.y += coef * n.y;
        F.z += coef * n.z;
    }

    return F;
}

// CALCULOS DINAMICOENERGETICOS DE POLIMERO //

void Calcular_fuerzas(Vector* F, Particula* P) {

    for (int i = 0; i < N_particulas; i++) {
        F[i].x = 0.0;
        F[i].y = 0.0;
        F[i].z = 0.0;
    }

    // ELASTICAS
    if (N_particulas > 1) {
        Vector Fe0 = Fuerza_ElasticaExtremo(P[0], P[1]);
        F[0].x += Fe0.x; F[0].y += Fe0.y; F[0].z += Fe0.z;
    }

    if (N_particulas > 1) {
        Vector FeN = Fuerza_ElasticaExtremo(P[N_particulas - 1], P[N_particulas - 2]);
        F[N_particulas - 1].x += FeN.x;
        F[N_particulas - 1].y += FeN.y;
        F[N_particulas - 1].z += FeN.z;
    }

    for (int i = 1; i < N_particulas - 1; i++) {
        Vector Fi = Fuerza_ElasticaIntermedio(P[i - 1], P[i], P[i + 1]);
        F[i].x += Fi.x;
        F[i].y += Fi.y;
        F[i].z += Fi.z;

        #ifdef COMPLEXMODEL
        Vector Fb[3];
        Fuerza_Bending(P[i - 1], P[i], P[i + 1], Fb);

        F[i - 1].x += Fb[0].x;
        F[i - 1].y += Fb[0].y;
        F[i - 1].z += Fb[0].z;

        F[i].x += Fb[1].x;
        F[i].y += Fb[1].y;
        F[i].z += Fb[1].z;

        F[i + 1].x += Fb[2].x;
        F[i + 1].y += Fb[2].y;
        F[i + 1].z += Fb[2].z;
        #endif
    }

	// EXTERNAS (SOLO SE APLICA AL EXTREMO)
        F[N_particulas-1].x += Fext.x;
        F[N_particulas-1].y += Fext.y;
        F[N_particulas-1].z += Fext.z;

#ifdef COMPLEXMODEL
	
        // LENNARD-JONES
    for (int i = 0; i < N_particulas; i++) {
        Vector FLJ = Fuerza_LennardJones(P[i]);
        F[i].x += FLJ.x;
        F[i].y += FLJ.y;
        F[i].z += FLJ.z;
		
        Vector FCl = Fuerza_Coulomb(P[i]);
        F[i].x += FCl.x;
        F[i].y += FCl.y;
		F[i].z += FCl.z;
    }
    for (int i = 1; i < N_particulas - 2; i++) {
        Vector fdih[4];
        Fuerza_Dihedral(P[i - 1], P[i], P[i + 1], P[i + 2], fdih);

        for (int j = 0; j < 4; j++) {
            F[i - 1 + j].x += fdih[j].x;
            F[i - 1 + j].y += fdih[j].y;
            F[i - 1 + j].z += fdih[j].z;
        }
    }

#endif
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