
// LIBRERIAS
#define _CRT_SECURE_NO_WARNINGS // Para VS 2022 (si usais otro entorno no creo que os dé problemas)
#include "stdio.h"
#include <string.h>
#include <stdlib.h>
#include <direct.h> 
#include <math.h>
#include <time.h>
//#define SWEEPMODE
#define COMPLEXMODEL


// FORMATO DE TRABAJO
typedef struct {
	double x;
	double y;
	double z;
} Vector;

typedef struct {
	int idx;
	double q;
	Vector pos;
	Vector vel;
	double Ecin;
	double Epot;
} Particula;

typedef Vector(*FuncionFuerzaExtremo)(Particula, Particula);
typedef Vector(*FuncionFuerzaIntermedio)(Particula, Particula, Particula);
// NUMEROS ALEATORIOS

#define NormRANu (2.3283063671E-10F)
extern unsigned int irr[256];
extern unsigned int ir1;
extern unsigned char ind_ran, ig1, ig2, ig3;
double N_Rand(void);
double N_Rand_Gauss(); 
void Ini_N_Rand(int SEMILLA);

//CONSTANTES DEL PROGRAMA Y TIPOS

#define PI acos(-1.0)
#define N_bins 50
# define N_particulas 64

// VARIABLES GLOBALES

extern Particula P[N_particulas];
extern double k;
extern double b;
extern Vector Fext;
extern double gamma_DP;
extern double KbT;
extern double m;
extern double dt;
extern double x0;
extern double v0;
extern double tmax;

extern double rc;
extern double sigma;
extern double eps;

extern double Q_ct;
extern double permitivity;

extern double kb;
extern double theta_0;

extern double kb_dih;
extern double phi_0;
extern int mult;
// FUNCIONES

void Inicializar(void);
double modulo(Vector r);
double Pesc_Norm(Vector r1, Vector r2);
double Pesc_NoNorm(Vector r1, Vector r2);
Vector Scalar_mult(Vector r, double lambda);
Vector Normalize(Vector v1);
Vector Pvect(Vector r1, Vector r2);
void verlet_estocastico_3D_extremo(Particula* P, Particula P2, FuncionFuerzaExtremo Fuerza);
void verlet_estocastico_3D_intermedio(Particula P2, Particula* P, Particula P3, FuncionFuerzaIntermedio Fuerza);
double Potencial_Extremo(Particula P1, Particula P2);
double Potencial_Intermedio(Particula P_ant, Particula P, Particula Psig);
double Potencial_LennardJones(Particula pi);
Vector Fuerza_LennardJones(Particula pi);
void Actualizar_Energias(Particula* P);
Vector Fuerza_ElasticaIntermedio(Particula Pant, Particula P, Particula Psig);
Vector Fuerza_ElasticaExtremo(Particula P1, Particula P2);
Vector CDM_uniforme(Vector* r, int N);
Vector resta(Vector r1, Vector r2);
double Radio_giro(Particula* P);
FILE* crear_archivo_xyz(int bloque);
FILE* crear_archivo_variables(int bloque, const char* cabecera);
void guardar_bloque_xyz(Particula* P, FILE* f, int paso_global);
void Escribir_buffer(FILE* f, double* buffer);
void Gnuplot_EnerCons(int bloques);
void Gnuplot_Rg(int bloques);
void crear_script_vmd(int N_bloques);
void Ajuste_Rg_en_N();
void Ajuste_Rg_en_k();
void Ajuste_Rg_en_b();
void Calcular_fuerzas(Vector* F, Particula* P);
void verlet_estocastico_3D(Particula* P);
double Longitud_efectiva(Particula* P, int norm);
void Gnuplot_L_eff(int bloques);
void verlet_estocastico_3D_extremofijo(Particula* P);
void Ajuste_L_eff_en_Fx();
double Potencial_Coulomb(Particula Pi);
Vector Fuerza_Coulomb(Particula pi);
void Fuerza_Dihedral(Particula P1, Particula P2, Particula P3, Particula P4, Vector F[4]);