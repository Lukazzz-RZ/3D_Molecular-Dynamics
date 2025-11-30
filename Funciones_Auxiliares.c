#include "head.h"

// DECLARACIÓN DE VARIABLES AUXILIARES GLOBALES //

unsigned int irr[256];
unsigned int ir1;
unsigned char ind_ran, ig1, ig2, ig3;


// FUNCIONES DE INICIALIZACION DE VARIABLES GLOBALES //

void Inicializar(void) {
    _mkdir("results");
    Ini_N_Rand(time(NULL));
    m = 1.0;
    k = 200.0;
    KbT = 1.0;
    b = 1.0;
    gamma_DP = 1.0;
    dt = 0.001;
    tmax = 25.0;
    Fext.x = 0.0;
    Fext.y = 0.0;
    Fext.z = 0.0;

    sigma = b / pow(2.0, 1.0 / 6.0);
    rc = 3 * sigma;
    eps = 1.0;

	permitivity = 0.1;
	Q_ct = 1.0; 

	kb = 100.0;
    theta_0 = PI + 30 * PI/180;

	kb_dih = 10.0;
    phi_0 = 10 * PI / 180;
    mult = 1;

    double radius = 10 * b;       
    double pitch = 5 * b;    
    double angle_step = 45.0;
    double q_ant = 1.;
    for (int j = 0; j < N_particulas; j++) {
        P[j].idx = j;
        double theta = PI/180 * (j * angle_step);
        P[j].pos.x = radius * cos(theta);
        P[j].pos.y = radius * sin(theta);
        P[j].pos.z = pitch * j / (2 * PI);
        P[j].vel.x = 0.0;
        P[j].vel.y = 0.0;
        P[j].vel.z = 0.0;

        P[j].Ecin = 0.5 * m * modulo(P[j].vel) * modulo(P[j].vel);
       
        if (j == 0 || j == N_particulas - 1) {
            P[j].Epot = (j == 0) ?
                Potencial_Extremo(P[j], P[j + 1]) :
                Potencial_Extremo(P[j], P[j - 1]);
            P[j].q = 0.;
        }
        else {
            P[j].Epot = Potencial_Intermedio(P[j - 1], P[j], P[j + 1]);

            int puente = j % 8;
            if (puente == 0) {
                P[j].q = +1.;      
            }
            else if (puente == 4) {
                P[j].q = -1.;      
            }
            else {
                P[j].q = 0.;     
            }
        }
        P[j].q *= Q_ct;
    }
}

// FUNCIONES DE CALCULO VECTORIAL //

double modulo(Vector r){
    return sqrt(r.x*r.x + r.y* r.y + r.z* r.z);
}

double Pesc_Norm(Vector r1, Vector r2){
    //VA NORM
    return (r1.x*r2.x+r1.y*r2.y+r1.z*r2.z)/modulo(r1)/modulo(r2);
}

double Pesc_NoNorm(Vector r1, Vector r2) {
    return (r1.x * r2.x + r1.y * r2.y + r1.z * r2.z);
}

Vector Pvect(Vector r1, Vector r2) {
    Vector result;
    result.x = r1.y * r2.z - r1.z * r2.y;
    result.y = r1.z * r2.x - r1.x * r2.z;
    result.z = r1.x * r2.y - r1.y * r2.x;

    return result;
}

Vector Normalize(Vector v1) {
    Vector result;
    double r = modulo(v1);
    result.x = v1.x / r;
    result.y = v1.y / r;
    result.z = v1.z / r;

    return result;

}

Vector Scalar_mult(Vector r, double lambda) {
    Vector result;
    result.x = r.x * lambda;
    result.y = r.y * lambda;
    result.z = r.z * lambda;
    return result;
}

Vector resta(Vector r1, Vector r2){
    Vector result;
    result.x = r1.x -r2.x;
    result.y = r1.y -r2.y;
    result.z = r1.z -r2.z;

    return result;
}

Vector CDM_uniforme(Vector* r, int N) {
    Vector r_cdm = { 0.0, 0.0, 0.0 };
    for (int i = 0; i < N; i++) {
        r_cdm.x += r[i].x;
        r_cdm.y += r[i].y;
        r_cdm.z += r[i].z;
    }
    r_cdm.x /= (double)N;
    r_cdm.y /= (double)N;
    r_cdm.z /= (double)N;
    return r_cdm;
}

double theta(Vector r1, Vector r2){
    //Devuelve angulo entre 2 vectores en rads entiendo
    return acos(Pesc_Norm(r1,r2));
}

// FUNCIONES DE NUMEROS ALEATORIOS //

void Ini_N_Rand(int Seed) {

    int INI, FACTOR, SUM, i;

    INI = Seed;
    FACTOR = 67397;
    SUM = 7364893;
    srand(Seed);

    for (i = 0; i < 256; i++) {
        INI = (INI * FACTOR + SUM);
        irr[i] = INI;
    }
    ind_ran = ig1 = ig2 = ig3 = 0;

    return;
}

double N_Rand(void) {

    double r;

    ig1 = ind_ran - 24;
    ig2 = ind_ran - 55;
    ig3 = ind_ran - 61;

    irr[ind_ran] = irr[ig1] + irr[ig2];
    ir1 = (irr[ind_ran] ^ irr[ig3]);
    ind_ran++;
    r = ir1 * NormRANu;

    return r;
}

double N_Rand_AB(double a, double b) {
    return a + (b - a) * N_Rand();
}

double N_Rand_Gauss() {
    return -sqrt(-2.0 * log(N_Rand())) * cos(2.0 * PI * N_Rand());
}


