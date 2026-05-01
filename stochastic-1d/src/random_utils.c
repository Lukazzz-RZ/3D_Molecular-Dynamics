#include "../include/simulation.h"

/* Random number generator state */
unsigned int irr[256];
unsigned int ir1;
unsigned char ind_ran;
unsigned char ig1;
unsigned char ig2;
unsigned char ig3;

void initialize_simulation(void) {
    MAKE_DIR("results");

    seed_random_generator((int)time(NULL));

    mass = 1.0;
    spring_constant = 1.0;
    thermal_energy = 0.2;
    well_depth = 0.5;
    friction = 5.0;

    time_step = 0.01;
    total_time = 100.0;

    initial_position = 1.0;
    initial_velocity = 1.0;

    external_force = 0.6;
}

void seed_random_generator(int seed) {
    int value = seed;
    const int factor = 67397;
    const int offset = 7364893;

    srand(seed);

    for (int i = 0; i < 256; i++) {
        value = value * factor + offset;
        irr[i] = value;
    }

    ind_ran = 0;
    ig1 = 0;
    ig2 = 0;
    ig3 = 0;
}

double uniform_random(void) {
    ig1 = ind_ran - 24;
    ig2 = ind_ran - 55;
    ig3 = ind_ran - 61;

    irr[ind_ran] = irr[ig1] + irr[ig2];
    ir1 = irr[ind_ran] ^ irr[ig3];

    ind_ran++;

    return ir1 * RANDOM_NORMALIZATION;
}

double uniform_random_range(double min, double max) {
    return min + (max - min) * uniform_random();
}

double gaussian_random(void) {
    double u1 = uniform_random();
    double u2 = uniform_random();

    if (u1 <= 0.0) {
        u1 = 1e-12;
    }

    return -sqrt(-2.0 * log(u1)) * cos(2.0 * PI * u2);
}

void test_random_number_generator(void) {
    const int number_of_points = 100000;

    FILE *data_file = fopen("results/random_test.txt", "w");
    if (!data_file) {
        perror("results/random_test.txt");
        return;
    }

    double previous = uniform_random();

    for (int i = 0; i < number_of_points; i++) {
        double current = uniform_random();
        fprintf(data_file, "%f\t%f\n", previous, current);
        previous = current;
    }

    fclose(data_file);

    FILE *plot_file = fopen("random_test.plt", "w");
    if (!plot_file) {
        perror("random_test.plt");
        return;
    }

    fprintf(plot_file, "set terminal pngcairo size 1000,800\n");
    fprintf(plot_file, "set output 'results/random_test.png'\n");
    fprintf(plot_file, "set xlabel 'u_n'\n");
    fprintf(plot_file, "set ylabel 'u_{n+1}'\n");
    fprintf(plot_file, "set title 'Random number generator test'\n");
    fprintf(plot_file, "plot 'results/random_test.txt' using 1:2 with points ps 0.1 notitle\n");

    fclose(plot_file);

    system("gnuplot random_test.plt");
}
