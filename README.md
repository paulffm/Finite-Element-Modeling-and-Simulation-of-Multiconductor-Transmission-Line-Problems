# Finite Element Modeling and Simulation of Multiconductor Transmission-Line Problems

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://github.com/paulffm/Finite-Element-Modeling-and-Simulation-of-Multiconductor-Transmission-Line-Problems/blob/main/LICENSE)

Welcome to the repository for Finite Element Modeling and Simulation of Multiconductor Transmission-Line Problems. This repository offers an in-depth analysis of multiconductor transmission line problems and provides a comprehensive simulation tool to effectively address such challenges through Finite Element Modeling. For a deeper theoretical understanding of multiconductor transmission lines, you can refer to the following paper:

- [Nonuniqueness of modal transformations for multiconductor transmission line problems](https://onlinelibrary.wiley.com/doi/10.1002/etep.2342)

<p align="center">
  <img src="mtlm_system.png"  alt="1" width = 595px height = 325px >
</p>

## Overview

The repository is structured into distinct modules: Transmission Line Simulation, Multitransmission Line Simulation, and Machine Slot Simulation.

### Transmission Line Simulation

#### Problem Setting

Consider a model of a coaxial cable with an axial length $l_z$. A current $I$ is uniformly distributed along the cross-section of a cylindrical copper wire with conductivity $\sigma_{Cu}$ and radius $r_1$, aligned with the $z$-direction. The wire is enclosed within a cylindrical, non-conducting insulation shell with inner radius $r_1$, outer radius $r_2$, permeability μs, and permittivity $\epsilon_s = \epsilon_0$, where $\mu_0$ and $\epsilon_0$ denote the permeability and permittivity of vacuum. The outer surface of the insulation shell is considered a perfect electric conductor.

<p align="center">
  <img src="coaxial_cable.png"  alt="1" width = 637px height = 198px >
</p>

In this module, we address the magnetostatic problem described above using **Python**. We begin by solving the problem analytically and visualizing the respective fields. Subsequently, we employ [Gmsh](https://gmsh.info), an open-source finite element mesh generator, to define the wire problem's geometry and specify regions of different materials. We then define boundary conditions and solve the magnetic field problem using finite element modeling on this meshed wire. Further, we conduct post-processing tasks such as comparing the numerical solution with the analytical one, evaluating convergence and discretization error. Finally, we utilize [Pyrit](https://www.temf.tu-darmstadt.de/emft/forschung_emft/software_1/software.en.jsp), a finite element solver for coupled electro- and magneto-quasistatic and thermal problems written in Python, to solve this problem and compare its solution with our previous results.

### Multitransmission Line Simulation

#### Problem Setting

<p align="center">
  <img src="power_cable.png"  alt="1" width = 736px height = 293px >
</p>

In this module, we develop our simulation tool for multiconductor transmission-line models (MTLMs). Similar to the previous problem, we mesh the geometry and specify the different materials using [Gmsh](https://gmsh.info). We then solve the electro- and magnetostatic, and the electro- and magnetoquasistatic problems with and without excitation in the frequency-domain. Specifically, we employ modal decomposition and the low-frequency approximation to numerically solve this problem.

### Machine Slot Simulation

#### Problem Setting

<p align="center">
  <img src="machine_slot.png"  alt="1" width = 584px height = 197px >
</p>

In this module, we implement our simulation tool for MTLMs. Similar to the previous problems, we mesh the geometry and specify the different materials using [Gmsh](https://gmsh.info). We then solve the electro- and magnetostatic problems on a machine slot with 36 single wires.

<p align="center">
  <img src="phidistrmitte.png"  alt="1" width = 640px height = 480px >
</p>

This README provides an overview of the repository's contents and highlights the problem settings and methodologies employed in each module. Feel free to explore the repository further for detailed implementations and simulations of multiconductor transmission-line problems.
