#include "../include/simulation.h"

/* Physical parameters */
double spring_constant;
double well_depth;
double external_force;
double friction;
double thermal_energy;
double mass;
double time_step;
double total_time;
double initial_position;
double initial_velocity;

/* Energies */
double harmonic_potential(double x) {
    return 0.5 * spring_constant * x * x;
}

double double_well_potential(double x) {
    return well_depth * (x * x - 1.0) * (x * x - 1.0);
}

double forced_double_well_potential(double x) {
    return double_well_potential(x) - external_force * x;
}

double kinetic_energy(double mass_value, double velocity) {
    return 0.5 * mass_value * velocity * velocity;
}

/* Forces. These functions return -dV/dx. */
double harmonic_force(double x) {
    return -spring_constant * x;
}

double double_well_force(double x) {
    return -4.0 * well_depth * x * (x * x - 1.0);
}

double forced_double_well_force(double x) {
    return double_well_force(x) + external_force;
}

/* Analysis */
double calculate_right_well_occupation(const char *filename) {
    FILE *file = fopen(filename, "r");

    if (!file) {
        perror(filename);
        return -1.0;
    }

    double position;
    double frequency;

    double total_probability = 0.0;
    double right_well_probability = 0.0;

    while (fscanf(file, "%lf %lf", &position, &frequency) == 2) {
        total_probability += frequency;

        if (position > 0.0) {
            right_well_probability += frequency;
        }
    }

    fclose(file);

    if (total_probability == 0.0) {
        fprintf(stderr, "Error: total histogram probability is zero.\n");
        return 0.0;
    }

    double probability = right_well_probability / total_probability;
    double standard_error = sqrt(
        probability * (1.0 - probability) / total_probability
    );

    printf("Right-well occupation standard error: %f\n", standard_error);

    return probability;
}
