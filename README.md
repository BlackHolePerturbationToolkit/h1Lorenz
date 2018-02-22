GSF_Lorenz_Schwarz_output: code to compute the first-order (in the mass-ratio) Lorenz-gauge metric perturbation
from a particle in a circular orbit about a Schwarzschild black hole

### Compile

Compile using `scons`

### Usage

In the notebook directory run:

1) ComputeRadialGrid.nb
2) h1ret_l0m0.nb

The run:

mpirun -n $numprocs ./GSF_circ_output r0 lmax