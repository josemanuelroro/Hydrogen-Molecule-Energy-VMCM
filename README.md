# Study of the Hydrogen Molecule using the MonteCarlo Variational Method

## Introduction

The calculation of energies for molecular systems is a problem of interest in chemistry and physics. The hydrogen molecule H2 is the simplest of the molecules, but it is of great importance. The first work on this system dates back to 1927 with Heitler and London [3] who were the first to apply quantum mechanics to molecules and whose results described the chemical bond. Since then, different studies have been carried out for this system such as those of James and Coolidge [4] in 1933 who used Hylleraas wave functions to study the H2 molecule or those carried out in 1960 by Kolos and Roothaan [5] whose work focused on wave functions describing the electronic, vibrational and rotational states of the fundamental state. Kolos and Wolniecwicz who achieved more accurate results than in their earlier work. More recent is the work of Cheng [1] 2018 which yielded the most accurate data to date and whose results represent a reference value for future relativistic and QED calculations of molecular energies.
The simplicity and the large amount of work done both theoretically and experimentally for the hydrogen molecule allows it to be used as a study tool for testing different methods, calculating physical constants and testing the theory of quantum electrodynamics.
In this work we will calculate the binding energy of the hydrogen molecule and the nuclear equilibrium distance for the ground state using the Variational Monte Carlo Method.
We will start by introducing the hydrogen molecule explaining its geometry and the approximations used to simplify the system and we will introduce the hamiltonian of the molecule. Then we will talk about the Variational Monte Carlo Method, what the Variational Method consists of and the reason for using a sampling algorithm such as the Metropolis Algorithm. We will also make a detailed study of the wave functions used since they play a crucial role in the implementation and accuracy of the method.


## Hydrogen Molecule

The neutral hydrogen molecule H2 is a diatomic molecule composed of two nuclei and two electrons. The nuclei are composed of two protons. It is one of the simplest systems we can use to test the validity of the Variational Monte Carlo Method.

![alt text](https://github.com/josemanuelroro/h2/blob/main/image.png?raw=true)

A modelization of the system for its study is the scheme shown in Figure 1. In it, each electron is bound to each of the nuclei A and B. The positions of electrons 1 and 2 are measured with respect to the origin of coordinates O, the position where nucleus A is located, and R represents the nuclear separation distance.
The Schrodïnger Equation for the hydrogen molecule, obviating relativistic and spin effects, has the form [2]:

![](https://latex.codecogs.com/gif.latex?%5BT_%7Bn%7D&plus;T_%7Be%7D&plus;V%5D%5Cpsi_%7BT%7D%3DE%5Cpsi_%7BT%7D)

Where Tn represents the kinetic energy operator of each nucleus A and B:

<div align="center">![](https://latex.codecogs.com/gif.latex?T_%7Bn%7D%3D-%5Cfrac%7B1%7D%7B2%7D%5Cnabla%5E%7B2%7D_%7BA%7D-%5Cfrac%7B1%7D%7B2%7D%5Cnabla%5E%7B2%7D_%7BB%7D)</div>

