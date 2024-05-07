# Finite Element Modeling and Simulation of Multiconductor Transmission-Line Problems

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://github.com/paulffm/Discrete-Time-Diffusion-Models-for-Discrete-Data/blob/main/LICENSE)

This repository analyzes multiconductor transmission line problems and provides a simulation tool to solve such problems through Finite Element Modeling. For further theoretical insight into multiconductor transmission lines, refer to the following paper:

- [Nonuniqueness of modal transformations for multiconductor transmission line problems](https://onlinelibrary.wiley.com/doi/10.1002/etep.2342)

<p align="center">
  <img src="ct_forwardrev_process2.png"  alt="1" width = 876px height = 621px >
</p>

## Overview

The repository is divided in seperate modules: Transmission Line Simulation, Multitransmission Line Simulation and Machine Slot Simulation.

### Transmission Line Simulation
#### Problem Setting
A model of a coaxial cable with axial length $l_z$ is considered. A current $I$ is homogeneously distributed along the cross section of a cylindrical copper wire with conductivity $\sigma_{Cu}$, radius $r_1$, aligned with the $z$-direction. The wire is enclosed within a cylindrical, non-conducting insulation shell with inner radius $r_1$, outer radius $r_2$ and with permeability μs and permittivity $\eps_s = \eps_0$, where $\mu_0$ and $\eps_0$ are the permeability and permittivity of vacuum. The outer surface of the insulation shell is considered as a perfect electric conductor.

In this module we solve the described magnetostatic problem in **Python**. We start by solving the problem analytical and plotting the respective fields.
[Pyrit](https://www.temf.tu-darmstadt.de/emft/forschung_emft/software_1/software.en.jsp)
### Multitransmission Line Simulation

### Machine Slot Simulation

