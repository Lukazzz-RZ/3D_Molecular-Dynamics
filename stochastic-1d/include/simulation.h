#ifndef SIMULATION_H
#define SIMULATION_H

#define _CRT_SECURE_NO_WARNINGS

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#ifdef _WIN32
    #include <direct.h>
    #define MAKE_DIR(path) _mkdir(path)
#else
    #include <sys/stat.h>
    #include <sys/types.h>
    #define MAKE_DIR(path) mkdir(path, 0777)
#endif

#define PI acos(-1.0)
#define N_BINS 50
#define RANDOM_NORMALIZATION 2.3283063671E-10F

typedef double (*ForceFunction)(double x);

/* Random number generator state */
extern unsigned int irr[256];
extern unsigned int ir1;
extern unsigned char ind_ran;
extern unsigned char ig1;
extern unsigned char ig2;
extern unsigned char ig3;

/* Physical parameters */
extern double spring_constant;
extern double well_depth;
extern double external_force;
extern double friction;
extern double thermal_energy;
extern double mass;
extern double time_step;
extern double total_time;
extern double initial_position;
extern double initial_velocity;

/* Initialization */
void initialize_simulation(void);

/* Random numbers */
void seed_random_generator(int seed);
double uniform_random(void);
double gaussian_random(void);
double uniform_random_range(double min, double max);

/* Forces */
double harmonic_force(double x);
double double_well_force(double x);
double forced_double_well_force(double x);

/* Energies */
double harmonic_potential(double x);
double double_well_potential(double x);
double forced_double_well_potential(double x);
double kinetic_energy(double mass, double velocity);

/* Numerical integrators */
void euler_step(double *x, double *v, ForceFunction force);
void stochastic_euler_step(double *x, double *v, ForceFunction force);
void verlet_step(double *x, double *v, ForceFunction force);
void stochastic_verlet_step(double *x, double *v, ForceFunction force);
void rk4_step(double *x, double *v, ForceFunction force);
void stochastic_rk2_step(double *x, double *v, ForceFunction force);

/* Output and analysis */
FILE *create_data_file(const char *filename, const char *header);
void write_simulation_row(FILE *file, double *buffer);

void plot_energy(int number_of_blocks);
void plot_position_velocity(int number_of_blocks);

void generate_distribution_histogram(
    char *column_name,
    ForceFunction potential,
    int number_of_blocks,
    int number_of_data_points
);

void generate_residence_time_histogram(int log_scale);
void process_residence_times(double *residence_times, int number_of_switches);

double calculate_right_well_occupation(const char *filename);
void plot_occupation_probability_force_sweep(void);

void test_random_number_generator(void);

#endif
