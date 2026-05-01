# Molecular Simulation

This project explores stochastic dynamics and molecular simulation, progressing from simple 1D systems to polymer modeling.

It is divided into two independent modules:

- Stochastic 1D dynamics (Langevin systems)
- Polymer molecular dynamics

---

## Part I — Stochastic 1D Dynamics

This module studies the motion of a particle under Langevin dynamics, combining deterministic forces, friction, and thermal noise.

### Model

The system is governed by:

- External potential
- Friction proportional to velocity
- Gaussian stochastic force

This leads to stochastic differential equations describing the particle trajectory.

---

### Numerical Methods

Several integration schemes are implemented and compared:

- Stochastic Euler
- Stochastic Verlet
- Runge-Kutta (2nd order)

The comparison focuses on stability and convergence properties.

---

### Harmonic Oscillator

The harmonic potential is used as a validation test:

- Energy equipartition is recovered
- Kinetic and potential energy converge to expected values


::contentReference[oaicite:0]{index=0}


---

### Double-Well Potential

A bistable system is analyzed:

- Two stable states separated by an energy barrier
- Noise induces transitions between wells

Main observations:

- Position distributions depend on temperature and friction
- Residence times characterize switching dynamics
- Equilibrium occupation emerges at long times


::contentReference[oaicite:1]{index=1}


---

### External Force

Adding a constant force:

- Breaks symmetry between wells
- Produces biased occupation probabilities

The simulation agrees with theoretical predictions for small perturbations.

---

## Part II — Polymer Dynamics

This module simulates a polymer chain using molecular dynamics.

---

### Model

The polymer is represented as a chain of particles:

- Connected by elastic bonds
- Subject to thermal fluctuations
- Governed by Langevin dynamics

This corresponds to a bead-spring model inspired by the FJC framework.

---

### Thermodynamics

The simulation verifies:

- Equipartition of energy
- Stability of kinetic and potential contributions

---

### Radius of Gyration

A key observable is the radius of gyration:

- Measures spatial extension of the polymer
- Compared with theoretical predictions


::contentReference[oaicite:2]{index=2}


---

### Force-Extension Behavior

The response of the polymer to external force is studied:

- Compared with the Langevin function
- Distinction between rigid and elastic regimes


::contentReference[oaicite:3]{index=3}


---

### Extended Model

Additional interactions are included:

- Lennard-Jones (repulsion + attraction)
- Electrostatic forces
- Bending energy (chain stiffness)

These extensions allow more realistic molecular behavior.

---

## Project Structure
