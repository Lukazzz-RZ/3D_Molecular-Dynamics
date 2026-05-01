#include "../include/simulation.h"

FILE *create_data_file(const char *filename, const char *header) {
    FILE *file = fopen(filename, "w");

    if (!file) {
        perror(filename);
        return NULL;
    }

    fprintf(file, "%s\n", header);

    return file;
}

void write_simulation_row(FILE *file, double *buffer) {
    fprintf(
        file,
        "%f\t%f\t%f\t%f\t%f\t%f\n",
        buffer[0],
        buffer[1],
        buffer[2],
        buffer[3],
        buffer[4],
        buffer[5]
    );
}

void plot_energy(int number_of_blocks) {
    FILE *plot_file = fopen("plot_energy.plt", "w");

    if (!plot_file) {
        perror("plot_energy.plt");
        return;
    }

    fprintf(plot_file, "set terminal pngcairo size 1000,800\n");
    fprintf(plot_file, "set output 'results/energy.png'\n");
    fprintf(plot_file, "set xlabel 'Time'\n");
    fprintf(plot_file, "set ylabel 'Energy'\n");
    fprintf(plot_file, "set title 'Energy vs Time'\n");
    fprintf(plot_file, "plot ");

    for (int i = 1; i <= number_of_blocks; i++) {
        fprintf(
            plot_file,
            "'results/data_%03d.txt' using 1:6 with lines title 'Block %d'%s",
            i,
            i,
            (i < number_of_blocks) ? ", " : "\n"
        );
    }

    fclose(plot_file);
    system("gnuplot plot_energy.plt");
}

void plot_position_velocity(int number_of_blocks) {
    FILE *plot_file = fopen("plot_position_velocity.plt", "w");

    if (!plot_file) {
        perror("plot_position_velocity.plt");
        return;
    }

    fprintf(plot_file, "set terminal pngcairo size 1000,800\n");
    fprintf(plot_file, "set output 'results/position_velocity.png'\n");
    fprintf(plot_file, "set xlabel 'Time'\n");
    fprintf(plot_file, "set ylabel 'Value'\n");
    fprintf(plot_file, "set title 'Position and Velocity vs Time'\n");
    fprintf(plot_file, "plot ");

    for (int i = 1; i <= number_of_blocks; i++) {
        fprintf(
            plot_file,
            "'results/data_%03d.txt' using 1:2 with lines title 'Position (block %d)', "
            "'results/data_%03d.txt' using 1:3 with lines title 'Velocity (block %d)'%s",
            i, i, i, i,
            (i < number_of_blocks) ? ", " : "\n"
        );
    }

    fclose(plot_file);
    system("gnuplot plot_position_velocity.plt");
}

void generate_distribution_histogram(
    char *column_name,
    ForceFunction potential,
    int number_of_blocks,
    int number_of_data_points
) {
    (void)column_name;
    (void)potential;
    (void)number_of_data_points;

    MAKE_DIR("results/histograms");

    char histogram_file[] = "results/histograms/position_histogram.txt";

    FILE *hist_file = fopen(histogram_file, "w");
    if (!hist_file) {
        perror(histogram_file);
        return;
    }

    double min = -2.0;
    double max = 2.0;
    double bin_width = (max - min) / N_BINS;

    int histogram[N_BINS];
    memset(histogram, 0, sizeof(histogram));

    char filename[128];

    for (int block = 1; block <= number_of_blocks; block++) {
        snprintf(filename, sizeof(filename), "results/data_%03d.txt", block);

        FILE *data = fopen(filename, "r");
        if (!data) {
            perror(filename);
            continue;
        }

        char header[512];
        fgets(header, sizeof(header), data);

        double time, position, velocity, ke, pe, te;

        while (fscanf(data, "%lf %lf %lf %lf %lf %lf",
                      &time, &position, &velocity, &ke, &pe, &te) == 6) {

            int bin = (int)((position - min) / bin_width);

            if (bin >= 0 && bin < N_BINS) {
                histogram[bin]++;
            }
        }

        fclose(data);
    }

    for (int i = 0; i < N_BINS; i++) {
        double center = min + (i + 0.5) * bin_width;
        fprintf(hist_file, "%f %d\n", center, histogram[i]);
    }

    fclose(hist_file);

    FILE *plot = fopen("plot_histogram.plt", "w");
    if (!plot) {
        perror("plot_histogram.plt");
        return;
    }

    fprintf(plot, "set terminal pngcairo size 1000,800\n");
    fprintf(plot, "set output 'results/histograms/position_histogram.png'\n");
    fprintf(plot, "set style fill solid\n");
    fprintf(plot, "set boxwidth %f\n", bin_width);
    fprintf(plot, "plot '%s' using 1:2 with boxes notitle\n", histogram_file);

    fclose(plot);
    system("gnuplot plot_histogram.plt");
}

void process_residence_times(double *residence_times, int number_of_switches) {
    MAKE_DIR("results/histograms");

    FILE *file = fopen("results/histograms/residence_times.txt", "w");

    if (!file) {
        perror("residence_times.txt");
        return;
    }

    for (int i = 0; i < number_of_switches; i++) {
        fprintf(file, "%f\n", residence_times[i]);
    }

    fclose(file);
}

void generate_residence_time_histogram(int log_scale) {
    FILE *plot = fopen("plot_residence.plt", "w");
    if (!plot) {
        perror("plot_residence.plt");
        return;
    }

    fprintf(plot, "set terminal pngcairo size 1000,800\n");
    fprintf(plot, "set output 'results/histograms/residence_times.png'\n");

    if (log_scale) {
        fprintf(plot, "set logscale y\n");
    }

    fprintf(plot, "bin_width = 0.1\n");
    fprintf(plot, "bin(x) = bin_width * floor(x/bin_width)\n");

    fprintf(plot,
        "plot 'results/histograms/residence_times.txt' using "
        "(bin($1)):(1.0) smooth freq with boxes notitle\n"
    );

    fclose(plot);
    system("gnuplot plot_residence.plt");
}

void plot_occupation_probability_force_sweep(void) {
    FILE *plot = fopen("plot_force_sweep.plt", "w");
    if (!plot) {
        perror("plot_force_sweep.plt");
        return;
    }

    fprintf(plot, "set terminal pngcairo size 1000,800\n");
    fprintf(plot, "set output 'results/force_sweep.png'\n");
    fprintf(plot, "set xlabel 'Force'\n");
    fprintf(plot, "set ylabel 'Occupation probability'\n");

    fprintf(plot,
        "plot 'results/force_sweep.txt' using 1:2 with linespoints title 'P(x>0)'\n"
    );

    fclose(plot);
    system("gnuplot plot_force_sweep.plt");
}
