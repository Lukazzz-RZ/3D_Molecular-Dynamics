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
    k = 100.0;
    KbT = 1.0;
    b = 1.0;
    gamma_DP = 1.0;
    dt = 0.01;
    tmax = 1000.0;
    Fext.x = 0.0;
    Fext.y = 0.0;
    Fext.z = 0.0;

    sigma = b / pow(2.0, 1.0 / 6.0);
    //sigma = 0;
    rc = 3 * sigma;
    eps = 1.2;
    kb = 1.;
    theta_0 = PI/2;

    double q_ant=1.;

    for (int j = 0; j < N_particulas; j++) {
        P[j].pos.x = (double)j;
        P[j].pos.y = (double)0;
        P[j].pos.z = (double)0;
        P[j].vel.x = (double) j*0.2;
        P[j].vel.y = (double)j * 0.2;
        P[j].vel.z = (double)j * 0.2;
        P[j].Ecin = 0.5 * m * modulo(P[j].vel) * modulo(P[j].vel);
        if (j == 0) {
            P[j].Epot = Potencial_Extremo(P[j], P[j+1]);
            P[j].q = -1.;
            q_ant = P[j].q;

        }
        else if (j == N_particulas - 1) {
            P[j].Epot = Potencial_Extremo(P[j], P[j-1]);
        }
        else {
            P[j].Epot = Potencial_Intermedio(P[j-1], P[j], P[j+1]);
            if((j+1)%3==0){
                P[j].q= -1.*q_ant;
                q_ant = P[j].q;

            } 
		}

        P[j].q = P[j].q*1.;
        
    }
}

// FUNCIONES DE CALCULO VECTORIAL //

double modulo(Vector r){
    return sqrt(fabs(r.x*r.x + r.y*r.y + r.z*r.z));
}

double Pesc(Vector r1, Vector r2){
    //NO NORM
    return (r1.x*r2.x+r1.y*r2.y+r1.z*r2.z)/modulo(r1)/modulo(r2);
}

Vector Vprod(Vector r1, Vector r2){
    Vector result;
    result.x = r1.y*r2.z - r1.z*r2.y;
    result.y = r1.z*r2.x - r1.x*r2.z;
    result.z = r1.x*r2.y - r1.y*r2.x;

    return result;
}

Vector resta(Vector r1, Vector r2){
    Vector result;
    result.x = r1.x -r2.x;
    result.y = r1.y -r2.y;
    result.z = r1.z -r2.z;

    return result;
}

Vector Normalizador (Vector v1){
    Vector result;
    double r = modulo(v1);
    result.x = v1.x/r;
    result.y = v1.y/r;
    result.z = v1.z/r;

    return result;

}

Vector Scalar_mult(Vector r, double lambda){
    Vector result;
        result.x = r.x*lambda;
        result.y = r.y*lambda;
        result.z = r.z*lambda;
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
    return acos(Pesc(r1,r2));
}

// FUNCIONES DE NUMEROS ALEATORIOS //

    // Inicializa el generador de numeros aleatorios
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

    // Generador de numeros aleatorios uniforme en [0,1)
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

    // Generador de numeros aleatorios uniforme en [a,b)
double N_Rand_AB(double a, double b) {
    return a + (b - a) * N_Rand();
}

    // Generador de numeros aleatorios con distribucion gaussiana N(0,1)
double N_Rand_Gauss() {
    return -sqrt(-2.0 * log(N_Rand())) * cos(2.0 * PI * N_Rand());
}


