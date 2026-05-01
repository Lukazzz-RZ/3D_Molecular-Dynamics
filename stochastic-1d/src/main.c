#include "../include/simulation.h"

int main(void) {
    initialize_simulation();

    const char *file_header =
        "time\tposition\tvelocity\tkinetic_energy\tpotential_energy\ttotal_energy";

    const int number_of_data_points = (int)(total_time / time_step);
    const int number_of_blocks = 10;
    const int block_size = number_of_data_points / number_of_blocks;

    int current_block = 1;
    int number_of_switches = 0;

    double buffer[6] = {0.0};
    double last_switch_time = 0.0;
    double previous_position = initial_position;

    double *residence_times = malloc((number_of_data_points + 1) * sizeof(double));
    if (!residence_times) {
        perror("Error allocating residence time array");
        return 1;
    }

    char filename[128];

    snprintf(
        filename,
        sizeof(filename),
        "results/data_%03d.txt",
        current_block
    );

    FILE *data_file = create_data_file(filename, file_header);
    if (!data_file) {
        free(residence_times);
        return 1;
    }

    buffer[0] = 0.0;
    buffer[1] = initial_position;
    buffer[2] = initial_velocity;
    buffer[3] = kinetic_energy(mass, initial_velocity);
    buffer[4] = forced_double_well_potential(initial_position);
    buffer[5] = buffer[3] + buffer[4];

    for (int step = 0; step < number_of_data_points; step++) {
        buffer[0] += time_step;

        double normalization = (buffer[0] / time_step) + 1.0;

        stochastic_verlet_step(
            &buffer[1],
            &buffer[2],
            forced_double_well_force
        );

        buffer[3] += (
            kinetic_energy(mass, buffer[2]) - buffer[3]
        ) / normalization;

        buffer[4] += (
            forced_double_well_potential(buffer[1]) - buffer[4]
        ) / normalization;

        buffer[5] = buffer[3] + buffer[4];

        if (step > 0 && previous_position * buffer[1] < 0.0) {
            residence_times[number_of_switches] = buffer[0] - last_switch_time;
            number_of_switches++;
            last_switch_time = buffer[0];
        }

        previous_position = buffer[1];

        if (step % 10 == 0) {
            write_simulation_row(data_file, buffer);
        }

        if (
            (step + 1) % block_size == 0 &&
            current_block < number_of_blocks
        ) {
            fclose(data_file);

            current_block++;

            snprintf(
                filename,
                sizeof(filename),
                "results/data_%03d.txt",
                current_block
            );

            data_file = create_data_file(filename, file_header);
            if (!data_file) {
                free(residence_times);
                return 1;
            }
        }
    }

    if (number_of_switches > 0) {
        residence_times[number_of_switches] = buffer[0] - last_switch_time;
        number_of_switches++;
    }

    fclose(data_file);

    printf("Final kinetic energy: %.6f\n", buffer[3]);
    printf("Final potential energy: %.6f\n", buffer[4]);
    printf("Final total energy: %.6f\n", buffer[5]);
    printf("Number of well transitions: %d\n", number_of_switches);

    printf("\nGenerating analysis plots...\n");

    generate_distribution_histogram(
        "position",
        forced_double_well_potential,
        number_of_blocks,
        number_of_data_points
    );

    process_residence_times(residence_times, number_of_switches);

    double occupation = calculate_right_well_occupation(
        "results/histograms/position_histogram.txt"
    );

    printf("Right-well occupation probability: %.6f\n", occupation);

    generate_residence_time_histogram(1);
    plot_energy(number_of_blocks);
    plot_position_velocity(number_of_blocks);
    plot_occupation_probability_force_sweep();

    free(residence_times);

    printf("Simulation completed.\n");

    return 0;
}
