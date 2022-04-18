# Study of the Hydrogen Molecule using the MonteCarlo Variational Method

## Introduction

The calculation of energies for molecular systems is a problem of interest in chemistry and physics. The hydrogen molecule H2 is the simplest of the molecules, but it is of great importance. The first work on this system dates back to 1927 with Heitler and London [3] who were the first to apply quantum mechanics to molecules and whose results described the chemical bond. Since then, different studies have been carried out for this system such as those of James and Coolidge [4] in 1933 who used Hylleraas wave functions to study the H2 molecule or those carried out in 1960 by Kolos and Roothaan [5] whose work focused on wave functions describing the electronic, vibrational and rotational states of the fundamental state. Kolos and Wolniecwicz who achieved more accurate results than in their earlier work. More recent is the work of Cheng [1] 2018 which yielded the most accurate data to date and whose results represent a reference value for future relativistic and QED calculations of molecular energies.
The simplicity and the large amount of work done both theoretically and experimentally for the hydrogen molecule allows it to be used as a study tool for testing different methods, calculating physical constants and testing the theory of quantum electrodynamics.
In this work we will calculate the binding energy of the hydrogen molecule and the nuclear equilibrium distance for the ground state using the Variational Monte Carlo Method.
We will start by introducing the hydrogen molecule explaining its geometry and the approximations used to simplify the system and we will introduce the hamiltonian of the molecule. Then we will talk about the Variational Monte Carlo Method, what the Variational Method consists of and the reason for using a sampling algorithm such as the Metropolis Algorithm. We will also make a detailed study of the wave functions used since they play a crucial role in the implementation and accuracy of the method.


## Hydrogen Molecule

The neutral hydrogen molecule H2 is a diatomic molecule composed of two nuclei and two electrons. The nuclei are composed of two protons. It is one of the simplest systems we can use to test the validity of the Variational Monte Carlo Method.
<p align="center">
<img width="460" height="300" src="https://github.com/josemanuelroro/h2/blob/main/image.png">
</p>
A modelization of the system for its study is the scheme shown in Figure 1. In it, each electron is bound to each of the nuclei A and B. The positions of electrons 1 and 2 are measured with respect to the origin of coordinates O, the position where nucleus A is located, and R represents the nuclear separation distance.
The Schrodïnger Equation for the hydrogen molecule, obviating relativistic and spin effects, has the form [2]:

<p align="center">
<img width="200" height="20" src="https://latex.codecogs.com/gif.latex?%5BT_%7Bn%7D&plus;T_%7Be%7D&plus;V%5D%5Cpsi_%7BT%7D%3DE%5Cpsi_%7BT%7D">
</p>

Where Tn represents the kinetic energy operator of each nucleus A and B:

<p align="center">
<img width="200" height="40" src="https://latex.codecogs.com/gif.latex?T_%7Bn%7D%3D-%5Cfrac%7B1%7D%7B2%7D%5Cnabla%5E%7B2%7D_%7BA%7D-%5Cfrac%7B1%7D%7B2%7D%5Cnabla%5E%7B2%7D_%7BB%7D">
</p>

Te represents the kinetic energy operator of electrons 1 and 2:
<p align="center">
<img width="200" height="40" src="https://latex.codecogs.com/gif.latex?T_%7Be%7D%3D-%5Cfrac%7B1%7D%7B2%7D%5Cnabla%5E%7B2%7D_%7B1%7D-%5Cfrac%7B1%7D%7B2%7D%5Cnabla%5E%7B2%7D_%7B2%7D">
</p>

V represents the columbian interactions between the different particles. Electron 1 with nucleus A, electron 2 with nucleus A, electron 1 with nucleus B, electron 2 with nucleus B, repulsion between electrons and repulsion between nuclei respectively:

<p align="center">
<img width="400" height="40" src="https://latex.codecogs.com/gif.latex?V%3D-%5Cfrac%7B1%7D%7B%7Cr_%7B1%7D%7C%7D-%5Cfrac%7B1%7D%7B%7Cr_%7B2%7D%7C%7D-%5Cfrac%7B1%7D%7B%7Cr_%7B1%7D-R%7C%7D-%5Cfrac%7B1%7D%7B%7Cr_%7B2%7D-R%7C%7D&plus;%5Cfrac%7B1%7D%7B%7Cr_%7B1%7D-r%7B2%7D%7C%7D&plus;%5Cfrac%7B1%7D%7B%7CR%7C%7D">
</p>

As a molecule, electronic, vibrational and rotational energy levels appear. To simplify the problem we use the Born-Oppenheimer approximation [2] which tells us that the wave function of nuclei and electrons within a molecule can be treated separately. This approximation is based on the high mass of the proton relative to that of the electron, about 1850 times larger. Therefore, we can consider the nuclei  fixed in space.
                                                                        𝑇𝑛=0
Moreover, most of the energy of the bound is contained in the electronic part of the wave function [2]. Then for the study of the hydrogen molecule, we can neglect the rotational and vibrational contributions.
Having considered all the above, we can write the Schrodïnger equation of the electronic contribution as:

<p align="center">
<img width="420" height="40" src="https://latex.codecogs.com/gif.latex?%28-%5Cfrac%7B1%7D%7B%7Cr_%7B1%7D%7C%7D-%5Cfrac%7B1%7D%7B%7Cr_%7B2%7D%7C%7D-%5Cfrac%7B1%7D%7B%7Cr_%7B1%7D-R%7C%7D-%5Cfrac%7B1%7D%7B%7Cr_%7B2%7D-R%7C%7D&plus;%5Cfrac%7B1%7D%7B%7Cr_%7B1%7D-r%7B2%7D%7C%7D&plus;%5Cfrac%7B1%7D%7B%7CR%7C%7D%29%5Cpsi%28%5Coverrightarrow%7Br%7D%29%3DE%5Cpsi%28%5Coverrightarrow%7Br%7D%29">
</p>

This is the hamiltonian we will work with. To calculate the energy of the fundamental level of the hydrogen molecule we must solve the eigenvalue problem, but this equation cannot be solved analytically and numerical integration methods are very computationally expensive. Therefore, to solve the problem we resort to the Variational Monte Carlo Method.
Atomic units will be used.
## MonteCarlo Variational Method
The Variational Method allows us to establish a bound on the energy of the fundamental state using test functions that depend on some variational parameters that we then minimize to make the energy a minimum. It tells us that the mean value of the Hamiltonian using a certain test wave function is greater than or equal to the ground state energy [6].

<p align="center">
<img width="100" height="50" src="https://latex.codecogs.com/gif.latex?E%5B%5CPsi%5D%3D%5Cfrac%7B%5Cint%20%5CPsi%5E%7B*%7D%5Cwidehat%7BH%7D%5CPsi%20dr%7D%7B%5Cint%20%5CPsi%5E%7B*%7D%5CPsi%20dr%7D%5Cgeq%20E_%7B0%7D">
</p>

If 𝛹(𝑟) depends on the variational parameter 𝛼 the energy will also depend on 𝛼 so we can vary this parameter to minimize the energy value.

Physically, the probabilities of finding a particle at a point in space are very small for most of space [7], therefore, using a distribution of uniformly generated random values to calculate the integral (2) is not a good method . This is where the Metropolis Algorithm comes into play which allows us to sample from a probability distribution and will be useful for calculating mean values [7][8].
We can rewrite equation (2) as:

<p align="center">
<img width="100" height="40" src="https://latex.codecogs.com/gif.latex?E_%7BL%7D%3D%5Cfrac%7B%5Cwidehat%7BH%7D%5CPsi%7D%7B%5CPsi%7D">
</p>

Because this method provides us with a discretization of the samples. This being the local energy, which only depends on the positions of the particles. This is the energy of each Monte Carlo simulation. The total energy of the system will be [8]:

<p align="center">
<img width="120" height="40" src="https://latex.codecogs.com/gif.latex?E%5B%5CPsi%5D%3D%5Clim_%7Bn%20%5Cto%20%5Cinfty%20%7D%5Cfrac%7B1%7D%7Bn%7D%5Csum_%7Bi%3D1%7D%5E%7Bn%7DE_%7BL%7D">
</p>

## Wave Functions

We will use wave functions that are combinations of hydrogenoids:
<p align="center">
<img width="120" height="20" src="https://latex.codecogs.com/gif.latex?%5Cpsi%28r_%7B1%7D%29%3Dexp%28-%5Calpha%7Cr_%7B1%7D%7C%29">
</p>
Taking into account the total antisymmetry that the wave function must present. We arrive at the following wave functions
<p align="center">
<img width="500" height="25" src="https://latex.codecogs.com/gif.latex?%5Cinline%20%5Cdpi%7B100%7D%20%5Cphi_%7B1%7D%3D%5Cfrac%7B1%7D%7B2%7D%5B%5Cpsi%28r_%7B1%7D%29&plus;%5Cpsi%28r_%7B1%7D-R%29%5D%5B%5Cpsi%28r_%7B2%7D%29&plus;%5Cpsi%28r_%7B2%7D-R%29%5D%7C00%3E">
</p>
<p align="center">
<img width="500" height="25" src="https://latex.codecogs.com/gif.latex?%5Cinline%20%5Cdpi%7B100%7D%20%5Cphi_%7B2%7D%3D%5Cfrac%7B1%7D%7B2%7D%5B%5Cpsi%28r_%7B1%7D%29-%5Cpsi%28r_%7B1%7D-R%29%5D%5B%5Cpsi%28r_%7B2%7D%29-%5Cpsi%28r_%7B2%7D-R%29%5D%7C00%3E">
</p>
<p align="center">
<img width="500" height="25" src="https://latex.codecogs.com/gif.latex?%5Cinline%20%5Cdpi%7B100%7D%20%5Cphi_%7B3%7D%3D%5Cfrac%7B1%7D%7B2%7D%5B%5Cpsi%28r_%7B1%7D%29%5Cpsi%28r_%7B2%7D%29-%5Cpsi%28r_%7B1%7D-R%29%5Cpsi%28r_%7B2%7D-R%29%5D%7C00%3E">
</p>
<p align="center">
<img width="500" height="25" src="https://latex.codecogs.com/gif.latex?%5Cinline%20%5Cdpi%7B100%7D%20%5Cphi_%7B4%7D%3D%5Cfrac%7B1%7D%7B2%7D%5B-%5Cpsi%28r_%7B1%7D%29%5Cpsi%28r_%7B2%7D-R%29-%5Cpsi%28r_%7B1%7D-R%29%5Cpsi%28r_%7B2%7D%29%5D%7C1MS%3E">
</p>

The wave functions seen above do not take into account the electronic repulsion of the hamiltonian. We have a 1/|𝑟1-𝑟2| term which is affected by the distance between the two electrons, electron density clouds are generated around the nuclei and unless the
distance of the nuclei is high, these densities can overlap causing the positions of the electrons to be very close causing the distance between them to tend to 0. If this distance tends to 0, a divergence occurs in the 1/|𝑟1-𝑟2| term of the Hamiltonian.
To avoid this divergence that occurs when the electron densities overlap  a correlation function called Padé-Jastrow function is introduced by [9][10]
<p align="center">
<img width="200" height="40" src="https://latex.codecogs.com/gif.latex?%5Cdpi%7B100%7D%20%5Cphi_%7Bj%7D%28r_%7B1%7D%2Cr_%7B2%7D%29%3Dexp%28-%5Cfrac%7BF%7D%7B2%281&plus;%5Cfrac%7Br_%7B12%7D%7D%7BF%7D%29%7D%29">
</p>
This Jastrow term behaves asymptotically as follows:
If the distance between electrons tends to 0, the term tends to some finite value. On the other hand, if the distance between the electrons becomes infinite the term tends to 1 and has no effect on the wave function of the system. Therefore, the wave function will have the following form where F will be a variational parameter.
<p align="center">
<img width="120" height="40" src="https://latex.codecogs.com/gif.latex?%5Cdpi%7B100%7D%20%5CPsi%28r_%7B1%7D%2Cr_%7B2%7D%29%3D%5Cphi_%7B1%7D%5Cphi_%7Bj%7D">
</p>
Another wave function can be proposed for the study of the hydrogen molecule by taking the functions and proposing a linear combination of them of the form:
<p align="center">
<img width="120" height="40" src="https://latex.codecogs.com/gif.latex?%5Cdpi%7B100%7D%20%5CPsi%28r_%7B1%7D%2Cr_%7B2%7D%29%3D%5Cphi_%7B1%7D&plus;%5Clambda%5Cphi_%7B2%7D">
With a variational parameter λ.
</p>
