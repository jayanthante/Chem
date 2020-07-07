# Quantum Dynamics
Various computational chemistry related codes

<u>Chebyshev.cpp</u> contains Chebyshev propagator I used for the propagation of wavepackets in multistate system.
It is written with the aim of including time-dependent terms in Hamiltonian though this is different from the propagator with time-dependent Hamiltnoian
This is best used when there aren't any non-adiabatic interactions OR diabatic interactions are present.

<u>main.cpp, time.cpp, numerov.cpp</u> - propagates excited states nuclear wavepackets using vibrational levels
Steps
1. For a given a potential energy surface, calculates the wavefunction Numerov method
2. To get the ith vibrational level, if the number of nodes is not equal to i, then goes to a higher/lower energy by bisection method
3. Repeats till all the vibrational levels are calculated
4. Calculates <chi_ij|H|chi_kl> for all the electronic state i and k and vibrational levels j and l
5. Calculates HC = i dC/dt
6. Solves the method using Gear Predictor-Corrector method (also verified using Runge-Kutta)
