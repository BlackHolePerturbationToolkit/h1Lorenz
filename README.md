# h1Lorenz

This repository contains code to compute the first-order (in the mass-ratio) Lorenz-gauge metric perturbation
from a particle in a circular orbit about a Schwarzschild black hole.

This code is a modified version of the code developed for [arXiv:1308.5223](https://arxiv.org/abs/1308.5223)

### Compile

Compile using `scons`

### Usage

In the inputs directory you will find the grids we currently use in the project. If you want to generate your own grid, use Mathematica notebook in notebooks/ComputeRadialGrid.nb.

To run the code use:

mpirun -n $numprocs ./h1Lorenz r0 lmax

where $numprocs should be at least 2 as one core is used to distribute work to the other cores (note: I've not tested how well the code works with n >= 2 in a long time).

### Output

The data is saved into the /data in HDF5 format. The metric perturbation is outputted in the Barack-Lousto-Sago basis.

### ToDo

 - Add calculation of C_in/C_out and expansion near horizon
 - Add calculation of F^R and h^{R1}_{ab}
 
The above are current performed in a separate Mathematica code which is repeating some of this calculation of this code. Be nice to have it all done in a single, fast code run.

### Authors

Niels Warburton
Sarp Akcay

### Papers on related topics

Lorenz-gauge decomposition into tensor spherical harmonics: [arXiv:0510019](https://arxiv.org/abs/gr-qc/0510019)

Lorenz-gauge circular orbits in frequency domain [arXiv:1012.5860](https://arxiv.org/abs/1012.5860)  
Lorenz-gauge eccentric orbits in the frequency domain [arXiv:1308.5223](https://arxiv.org/abs/1308.5223), [arXiv:1409.4419](https://arxiv.org/abs/1409.4419)

Lorenz-gauge circular orbits in time-domain: [arXiv:gr-qc/0701069](https://arxiv.org/abs/gr-qc/0701069)
Lorenz-gauge eccentric orbits in time-domain: [arXiv:1002.2386](https://arxiv.org/abs/1002.2386)
