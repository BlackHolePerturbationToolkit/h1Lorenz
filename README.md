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


### ToDo

 - Add calculation of monopole
 - Add calculation of C_in/C_out and expansion near horizon
 - Add calculation of F^R and h^{R1}_{ab}
 
The above are current performed in Mathematica which is repeating some of this calculation of this code. Be nice to have it all done in a single, fast code run.