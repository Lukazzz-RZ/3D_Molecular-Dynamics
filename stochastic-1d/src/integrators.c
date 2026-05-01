#include "../include/simulation.h"

void euler_step(double *x, double *v, ForceFunction force) {
    double acceleration = force(*x) / mass;

    *x += (*v) * time_step;
    *v += acceleration * time_step;
}

void stochastic_euler_step(double *x, double *v, ForceFunction force) {
    double noise_amplitude = sqrt(
        2.0 * friction * thermal_energy * time_step / mass
    );

    double acceleration = force(*x) / mass - friction * (*v);

    *x += (*v) * time_step;
    *v += acceleration * time_step + noise_amplitude * gaussian_random();
}

void verlet_step(double *x, double *v, ForceFunction force) {
    double acceleration = force(*x) / mass;

    *x += (*v) * time_step + 0.5 * acceleration * time_step * time_step;

    double new_acceleration = force(*x) / mass;

    *v += 0.5 * (acceleration + new_acceleration) * time_step;
}

void stochastic_verlet_step(double *x, double *v, ForceFunction force) {
    double noise_amplitude = sqrt(
        2.0 * friction * thermal_energy * time_step / mass
    );

    double acceleration = force(*x) / mass - friction * (*v);

    *x += (*v) * time_step + 0.5 * acceleration * time_step * time_step;

    double new_acceleration = force(*x) / mass - friction * (*v);

    *v += 0.5 * (acceleration + new_acceleration) * time_step;
    *v += noise_amplitude * gaussian_random();
}

void rk4_step(double *x, double *v, ForceFunction force) {
    double k1_x = *v;
    double k1_v = force(*x) / mass;

    double k2_x = *v + 0.5 * time_step * k1_v;
    double k2_v = force(*x + 0.5 * time_step * k1_x) / mass;

    double k3_x = *v + 0.5 * time_step * k2_v;
    double k3_v = force(*x + 0.5 * time_step * k2_x) / mass;

    double k4_x = *v + time_step * k3_v;
    double k4_v = force(*x + time_step * k3_x) / mass;

    *x += (time_step / 6.0) * (k1_x + 2.0 * k2_x + 2.0 * k3_x + k4_x);
    *v += (time_step / 6.0) * (k1_v + 2.0 * k2_v + 2.0 * k3_v + k4_v);
}

void stochastic_rk2_step(double *x, double *v, ForceFunction force) {
    double noise_amplitude = sqrt(
        2.0 * friction * thermal_energy * time_step / mass
    );

    double a1 = force(*x) / mass - friction * (*v);

    double x_mid = *x + 0.5 * (*v) * time_step;
    double v_mid = *v + 0.5 * a1 * time_step;

    double a2 = force(x_mid) / mass - friction * v_mid;

    *x += v_mid * time_step;
    *v += a2 * time_step + noise_amplitude * gaussian_random();
}
