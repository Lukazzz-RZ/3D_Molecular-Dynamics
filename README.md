# Molecular Simulation

This repository contains a computational physics project focused on stochastic dynamics and molecular simulation.

The project is divided into two independent modules:

- `stochastic-1d`: 1D Langevin dynamics for simple stochastic systems.
- `polymer-dynamics`: molecular dynamics simulation of polymer chains.

The first module studies simple one-dimensional systems such as the harmonic oscillator and the double-well potential. The second module extends the same physical ideas to a polymer chain model, including elastic bonds, thermal noise, structural observables, and additional molecular interactions.

---

## Table of Contents

- [Repository Structure](#repository-structure)
- [Part I: Stochastic 1D Dynamics](#part-i-stochastic-1d-dynamics)
  - [Physical Model](#physical-model)
  - [Numerical Integration Methods](#numerical-integration-methods)
  - [Random Number Generation](#random-number-generation)
  - [Harmonic Oscillator](#harmonic-oscillator)
  - [Double-Well Potential](#double-well-potential)
  - [Double-Well Position Distribution](#double-well-position-distribution)
  - [Residence Times and Well Occupation](#residence-times-and-well-occupation)
  - [External Force in the Double-Well Potential](#external-force-in-the-double-well-potential)
  - [Summary of the Stochastic 1D Module](#summary-of-the-stochastic-1d-module)
- [Part II: Polymer Dynamics](#part-ii-polymer-dynamics)
  - [Polymer Model](#polymer-model)
  - [Thermodynamic Validation](#thermodynamic-validation)
  - [Radius of Gyration](#radius-of-gyration)
  - [Analytical Scaling](#analytical-scaling)
  - [Radius of Gyration Simulation](#radius-of-gyration-simulation)
  - [Parameter Dependence](#parameter-dependence)
  - [Force-Extension Behavior](#force-extension-behavior)
- [Extended Polymer Model](#extended-polymer-model)
  - [Lennard-Jones Interactions](#lennard-jones-interactions)
  - [Electrostatic Interactions](#electrostatic-interactions)
  - [Bending Energy](#bending-energy)
  - [Combined Extended Model](#combined-extended-model)
- [Analytical Background](#analytical-background)
  - [Double-Well Occupation Correction](#double-well-occupation-correction)
  - [Analytical Radius of Gyration](#analytical-radius-of-gyration)
  - [Bending Force Derivation](#bending-force-derivation)
  - [Dihedral Force: Future Extension](#dihedral-force-future-extension)
- [Documentation](#documentation)
- [How to Run](#how-to-run)
- [Requirements](#requirements)
- [Main Concepts](#main-concepts)
- [Notes](#notes)
- [Figure List](#figure-list)

---

## Repository Structure

```text
molecular-simulation/
├── stochastic-1d/
│   ├── include/
│   ├── src/
│   ├── results/
│   └── Makefile
│
├── polymer-dynamics/
│   ├── include/
│   ├── src/
│   ├── results/
│   └── Makefile
│
├── figures/
│   ├── stochastic_methods.png
│   ├── random_gaussian_noise.png
│   ├── harmonic_oscillator_energy.png
│   ├── double_well_trajectory.png
│   ├── double_well_distribution.png
│   ├── double_well_external_force.png
│   ├── polymer_bead_spring_model.png
│   ├── polymer_radius_gyration_definition.png
│   ├── polymer_radius_gyration_simulation.png
│   ├── polymer_parameter_sweep.png
│   ├── polymer_force_extension.png
│   └── polymer_extended_model.png
│
└── README.md
```

---

# Part I: Stochastic 1D Dynamics

The `stochastic-1d` module studies the dynamics of a particle in one dimension under the effect of deterministic forces, viscous friction, and thermal stochastic noise.

The main goal is to simulate Langevin dynamics and analyze how numerical integration methods reproduce expected statistical and thermodynamic behavior.

---

## Physical Model

The motion is modeled using a Langevin-type equation. The particle is affected by:

- A deterministic force derived from a potential.
- A friction term proportional to velocity.
- A random Gaussian force representing thermal fluctuations.

In simplified form:

```text
m d2x/dt2 = F(x) - eta dx/dt + xi(t)
```

where:

- `m` is the mass.
- `F(x)` is the deterministic force.
- `eta` is the friction coefficient.
- `xi(t)` is the stochastic thermal force.

The thermal noise is controlled by the temperature through `kBT`.

---

## Numerical Integration Methods

Several stochastic integration methods are implemented and compared:

- Stochastic Euler
- Stochastic Verlet
- Second-order Runge-Kutta

These methods are used to integrate Langevin dynamics under different potentials and friction regimes.

![Stochastic integration methods](figures/stochastic_methods.png)

The comparison between methods is important because stochastic systems can be sensitive to timestep size, friction, and noise amplitude.

---

## Random Number Generation

Thermal noise is generated from Gaussian random variables. The project uses the Box-Muller transform to convert uniformly distributed random numbers into normally distributed variables.

![Gaussian noise generation](figures/random_gaussian_noise.png)

The stochastic force is essential for reproducing thermal equilibrium and allowing transitions between metastable states in systems such as the double-well potential.

---

## Harmonic Oscillator

The harmonic oscillator is used as the first validation system.

The potential is:

```text
V(x) = 1/2 k x^2
```

This system is useful because its equilibrium behavior is well known. The simulation checks whether the numerical methods recover the expected equipartition of energy.

Main checks:

- Kinetic energy approaches the expected thermal value.
- Potential energy approaches the expected thermal value.
- Total energy stabilizes statistically over long simulations.
- The position and velocity distributions are compatible with theoretical predictions.

![Harmonic oscillator energy](figures/harmonic_oscillator_energy.png)

This validation step is important before moving to more complex potentials.

---

## Double-Well Potential

The second system is a particle in a double-well potential.

The potential has two stable minima separated by an energy barrier:

```text
V(x) = A (x^2 - 1)^2
```

This system is useful for studying thermally activated transitions. Noise allows the particle to jump from one well to the other.

![Double-well trajectory](figures/double_well_trajectory.png)

Main quantities analyzed:

- Particle trajectory.
- Position distribution.
- Velocity distribution.
- Residence times in each well.
- Number of transitions between wells.
- Occupation probability of each well.

---

## Double-Well Position Distribution

The particle distribution depends strongly on temperature, friction, and barrier height.

At low temperature, the particle tends to remain trapped in one well for longer times. At higher temperature, transitions become more frequent.

![Double-well distribution](figures/double_well_distribution.png)

This allows the simulation to capture the relationship between thermal noise and barrier crossing.

---

## Residence Times and Well Occupation

For the symmetric double-well system, both wells should be occupied equally in the long-time limit.

The simulation measures:

```text
T_right
T_left
Number of transitions
```

These observables provide a statistical description of the switching process.

The residence time distribution is especially relevant because it characterizes how long the system remains in one metastable state before crossing the barrier.

---

## External Force in the Double-Well Potential

The model is extended by adding a constant external force:

```text
V(x) = A (x^2 - 1)^2 - F x
```

This breaks the symmetry between the two wells.

Consequences:

- One well becomes energetically favored.
- The occupation probability becomes biased.
- The probability of finding the particle at `x > 0` depends on the external force.

For small forces, the occupation probability can be approximated by a sigmoid-like expression.

![External force in double-well potential](figures/double_well_external_force.png)

This part connects stochastic simulation with equilibrium statistical mechanics.

---

## Summary of the Stochastic 1D Module

The `stochastic-1d` module demonstrates:

- Numerical integration of Langevin dynamics.
- Generation of thermal Gaussian noise.
- Validation through the harmonic oscillator.
- Barrier-crossing dynamics in a double-well potential.
- Statistical analysis of residence times and well occupation.
- Effect of a constant external force on equilibrium probabilities.

This module provides the foundation for the polymer molecular dynamics simulations.

---

# Part II: Polymer Dynamics

The `polymer-dynamics` module extends the simulation framework to a bead-spring polymer model.

The polymer is represented as a chain of particles connected by elastic bonds and evolving under molecular dynamics with thermal fluctuations.

---

## Polymer Model

The polymer is modeled as a chain of particles with neighboring particles connected by harmonic springs.

The bond potential is:

```text
V(r_i, r_i+1) = 1/2 k (|r_i - r_i+1| - b)^2
```

where:

- `k` is the bond stiffness.
- `b` is the equilibrium bond length.
- `r_i` is the position of particle `i`.

![Polymer bead-spring model](figures/polymer_bead_spring_model.png)

This model is inspired by the freely jointed chain framework, but uses elastic bonds instead of perfectly rigid segments.

---

## Thermodynamic Validation

The simulation checks whether the polymer model reaches a physically consistent thermal state.

The main validation criterion is energy equipartition.

The code analyzes:

- Mean kinetic energy.
- Mean potential energy.
- Stability of the total energy.
- Dependence on simulation parameters.

In the reference simulations, kinetic and potential contributions fluctuate around their expected equilibrium values.

This step ensures that the polymer model is numerically stable before analyzing structural observables.

---

## Radius of Gyration

The radius of gyration measures the spatial size of the polymer chain.

It is defined as the mean squared distance of the particles from the center of mass:

```text
Rg^2 = (1/N) sum_i |R_i - R_CM|^2
```

where:

- `N` is the number of particles.
- `R_i` is the position of particle `i`.
- `R_CM` is the center of mass.

![Radius of gyration definition](figures/polymer_radius_gyration_definition.png)

The radius of gyration is one of the most important observables in polymer physics because it quantifies the extension of the chain in space.

---

## Analytical Scaling

For an ideal polymer chain, the radius of gyration is expected to scale with the number of segments.

In the rigid-bond limit, the expected behavior is approximately:

```text
Rg^2 ~ b^2 (N - 1/N) / 6
```

The simulation compares numerical results with this theoretical prediction.

This makes the polymer module more than a visualization tool: it directly tests statistical physics predictions.

---

## Radius of Gyration Simulation

The simulation studies the behavior of `Rg^2` for fixed bond parameters.

Typical parameters include:

```text
kBT = 1.0
b = 1.0
k = 10.0 or 100.0
dt = 1e-2
```

![Radius of gyration simulation](figures/polymer_radius_gyration_simulation.png)

The results show the expected dependence of chain size on the number of particles and bond properties.

---

## Parameter Dependence

The project also studies how the polymer size changes when varying model parameters.

The main sweeps are:

- Bond stiffness `k`
- Equilibrium bond length `b`
- Number of particles `N`

![Polymer parameter sweep](figures/polymer_parameter_sweep.png)

These sweeps help verify how microscopic model parameters affect macroscopic chain observables.

---

## Force-Extension Behavior

The polymer is also studied under an external pulling force.

The goal is to compare the simulated end-to-end extension with the theoretical Langevin force-extension relation.

For a freely jointed chain, the extension is related to the Langevin function:

```text
L(F) / L_nat = coth(Fb / kBT) - kBT / (Fb)
```

For elastic bonds, an additional correction appears due to bond stretching.

![Polymer force-extension behavior](figures/polymer_force_extension.png)

This allows comparison between:

- Rigid polymer behavior.
- Elastic polymer behavior.
- Theoretical force-extension predictions.

---

# Extended Polymer Model

The polymer model is extended by adding more realistic molecular interactions.

These interactions include:

- Lennard-Jones interactions.
- Electrostatic interactions.
- Bending energy.

---

## Lennard-Jones Interactions

The Lennard-Jones potential models short-range repulsion and longer-range attraction:

```text
V_LJ(r) = 4 epsilon [ (sigma/r)^12 - (sigma/r)^6 ]
```

where:

- `sigma` controls the characteristic interaction distance.
- `epsilon` controls the energy scale.

This interaction is useful for representing excluded volume effects and attractive interactions between non-bonded particles.

A cutoff radius `rc` is also introduced to reduce computational cost and control the interaction range.

---

## Electrostatic Interactions

Electrostatic interactions are included by assigning charges to particles.

The model allows charged monomers to interact through Coulomb-like forces.

This is useful for representing simplified pH-dependent charge effects in polymer or protein-like chains.

---

## Bending Energy

The bending energy introduces stiffness into the polymer chain.

It depends on the angle between consecutive bonds:

```text
V_b(theta) = k_b [1 - cos(theta - theta_0)]
```

where:

- `k_b` controls bending stiffness.
- `theta_0` is the preferred angle.

This interaction involves three consecutive particles and controls the persistence length of the chain.

Higher bending stiffness produces more rigid polymer conformations.

---

## Combined Extended Model

The extended model combines:

- Elastic bonds.
- Lennard-Jones interactions.
- Electrostatic interactions.
- Bending forces.

![Extended polymer model](figures/polymer_extended_model.png)

This produces more realistic polymer behavior and allows comparison between the base bead-spring model and a more physically rich molecular model.

---

# Analytical Background

This project is supported by several analytical derivations used to validate and interpret the numerical simulations.

The main analytical results are related to:

- Occupation probability in a biased double-well potential.
- Radius of gyration of an ideal bead-spring polymer.
- Bending force derivation for three consecutive polymer particles.
- Dihedral force derivation as a possible future extension.

These results provide theoretical references for the stochastic and polymer simulations.

---

## Double-Well Occupation Correction

For the forced double-well system, the potential is:

```text
V(x) = A (x^2 - a^2)^2 - F x
```

When a small external force is applied, the two wells are no longer energetically equivalent.

For weak forces, the minimum positions can be approximated as:

```text
x_plus  =  a + F / (8 A a^2)
x_minus = -a + F / (8 A a^2)
```

The energy difference between both minima is approximately:

```text
Delta V ≈ -2 F a
```

Using a harmonic approximation around each minimum, the occupation probability of the right well can be corrected as:

```text
P_plus ≈ 1 / (1 + exp[-2 beta F a + 3F / (8 A a^3)])
```

This expression improves the simple two-state approximation by including the change in curvature of the potential wells caused by the applied force.

This analytical correction is used as a reference for the numerical occupation probability obtained from the stochastic simulation.

---

## Analytical Radius of Gyration

For the polymer module, the chain is modeled as an open bead-spring chain with harmonic bonds between neighboring particles:

```text
V(r_i, r_i+1) = 1/2 k (|r_i+1 - r_i| - b)^2
```

The radius of gyration is defined as:

```text
Rg^2 = (1/N) sum_i |r_i - R_CM|^2
```

An equivalent expression in terms of pairwise distances is:

```text
Rg^2 = (1 / 2N^2) sum_i sum_j < |r_i - r_j|^2 >
```

Using bond vectors:

```text
u_i = r_i+1 - r_i
```

and assuming independent, isotropically distributed bonds at equilibrium, one obtains:

```text
< |r_i - r_j|^2 > = |i - j| l^2
```

where:

```text
l^2 = < |u_i|^2 >
```

This leads to:

```text
Rg^2 = (l^2 / 6) (N - 1/N)
```

The effective squared bond length is obtained from the thermal radial distribution of the harmonic bond:

```text
l^2 = [b^4 + 6b^2/(beta k) + 3/(beta k)^2] / [b^2 + 1/(beta k)]
```

In the rigid-bond limit:

```text
beta k >> 1
```

this reduces to:

```text
l^2 ≈ b^2
```

and therefore:

```text
Rg ≈ (b / sqrt(6)) sqrt(N)
```

This analytical scaling is used to validate the simulated radius of gyration.

---

## Bending Force Derivation

The extended polymer model includes bending interactions involving three consecutive particles.

For three particles:

```text
P_previous, P, P_next
```

define the bond vectors:

```text
a = r_previous - r_P
b = r_next - r_P
```

The bending potential is:

```text
V_b(theta) = k_b [1 - cos(theta - theta_0)]
```

where `theta` is the angle between the two bond vectors.

The derivative of the potential is:

```text
dV_b / dtheta = k_b sin(theta - theta_0)
```

The resulting forces on the three particles can be written as:

```text
F_previous = (dV_b/dtheta) [b_hat - cos(theta) a_hat] / [|a| sin(theta)]

F_next = (dV_b/dtheta) [a_hat - cos(theta) b_hat] / [|b| sin(theta)]

F_center = -(F_previous + F_next)
```

This guarantees conservation of total force for the bending interaction.

The bending term introduces chain stiffness and allows the model to move beyond a purely freely jointed chain description.

---

## Dihedral Force: Future Extension

A possible future extension is the inclusion of dihedral interactions involving four consecutive particles.

The proposed dihedral potential has the form:

```text
V_dih(phi) = k_dih [1 - cos(n phi - phi_0)]
```

where:

- `phi` is the dihedral angle.
- `k_dih` is the dihedral stiffness.
- `n` is the multiplicity.
- `phi_0` is the preferred phase.

The derivation considers four particles:

```text
P1, P2, P3, P4
```

with bond vectors:

```text
b1 = r2 - r1
b2 = r3 - r2
b3 = r4 - r3
```

and plane normals:

```text
c1 = b1 x b2
c2 = b2 x b3
```

The angle can be defined from:

```text
cos(phi) = (c1 · c2) / (|c1| |c2|)
```

However, a robust implementation should use a signed definition of the dihedral angle, since the cosine alone does not distinguish chirality. For this reason, the dihedral force is left as future work rather than part of the current implemented model.

---

# Documentation

The repository includes two complementary levels of documentation:

- The source code, which implements the simulations.
- The analytical derivations, which justify and validate several of the numerical results.

The analytical material is especially useful for understanding:

- Why the biased double-well occupation follows a corrected sigmoid-like curve.
- Why the polymer radius of gyration scales as sqrt(N).
- How bending forces are obtained from an angular potential.
- Why a signed dihedral angle is needed for future torsional interactions.

---

# How to Run

Each module can be compiled independently.

## Stochastic 1D

```bash
cd stochastic-1d
make
./simulation
```

## Polymer Dynamics

```bash
cd polymer-dynamics
make
./simulation
```

Some plotting features may require `gnuplot`.

---

# Requirements

The project is written in C.

Recommended tools:

```text
gcc
make
gnuplot
```

Optional visualization tools:

```text
VMD
```

---

# Main Concepts

This repository covers:

- Langevin dynamics.
- Stochastic differential equations.
- Gaussian thermal noise.
- Numerical integration schemes.
- Harmonic oscillator validation.
- Double-well barrier crossing.
- Residence time statistics.
- Biased occupation in asymmetric double-well potentials.
- Harmonic approximation around potential minima.
- Polymer bead-spring models.
- Radius of gyration.
- Analytical scaling of polymer size.
- Force-extension behavior.
- Lennard-Jones interactions.
- Electrostatic interactions.
- Bending energy.
- Future extension to dihedral forces.

---

# Notes

This project is intended as a computational physics and molecular simulation study.

The code is designed to explore physical behavior rather than provide a production molecular dynamics engine.

The structure separates simple stochastic models from polymer simulations so that each part can be studied and extended independently.

The analytical derivations are included as theoretical support for the numerical results. They are not required to run the code, but they explain the physical formulas used to interpret the simulations.

---

# Figure List

The README expects the following files inside the `figures/` folder:

```text
figures/
├── stochastic_methods.png
├── random_gaussian_noise.png
├── harmonic_oscillator_energy.png
├── double_well_trajectory.png
├── double_well_distribution.png
├── double_well_external_force.png
├── polymer_bead_spring_model.png
├── polymer_radius_gyration_definition.png
├── polymer_radius_gyration_simulation.png
├── polymer_parameter_sweep.png
├── polymer_force_extension.png
└── polymer_extended_model.png
```
