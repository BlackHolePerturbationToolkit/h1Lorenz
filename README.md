# h1Lorenz

This repository contains code to compute the first-order (in the mass-ratio) Lorenz-gauge metric perturbation
from a particle in a circular orbit about a Schwarzschild black hole.

The code computes the tensor spherical harmonic modes of the metric perturbation and its first radial derivative on a supplied grid of radial values.
The spherical harmonic basis used is the Barack-Lousto-Sago basis.
This code is a modified version of the code developed for [arXiv:1308.5223](https://arxiv.org/abs/1308.5223)

### Dependencies and compilation

This software needs the GNU Scientific Library, Scons, and OpenMPI to compile and run.

Compile using `scons`

### Usage

First you need to generate an input grid. For this you can use the Mathematica notebook in notebooks/ComputeRadialGrid.nb.

To run the code use:

`mpirun -n $numprocs ./h1Lorenz r0 lmax gridfile outdir`

where ``$numprocs`` should be at least 2 as one core is used to distribute work to the other cores (note: we've not tested how well the code works with n >= 2 in a long time), `r0` is the particle's radius (in Schwarzschild coordinates), `lmax` is the maximum l value to compute, gridfile is the location of the gridfile, and outdir is the output directory.

For example, after runnign the ComputeRadialGrid.nb notebook for r0=8.1 you could run

`./mpirun -n 2 15 input/radial_grid_r8.1.h5 data/fields_r8.1/`

### Output

The code output the all the components of the first-order metric perturbation in the Lorenz gauge. For some modes, some of the asymptotic amplitudes are also computed.
The data is saved using the HDF5 format. 
There is a Mathematica notebook (Loadh1Lorenz.nb) in the notebooks/ subfolder which shows how to read in the data.


In general the l=even, m=0 modes are not computed accurately (or at all) due to long integration regions for the homogeneous solutions. For now and so these should be replaced with another calculation.


### Authors

Niels Warburton
Sarp Akcay

### Papers on related topics

Lorenz-gauge decomposition into tensor spherical harmonics: [arXiv:0510019](https://arxiv.org/abs/gr-qc/0510019)

Lorenz-gauge circular orbits in frequency domain [arXiv:1012.5860](https://arxiv.org/abs/1012.5860)  
Lorenz-gauge eccentric orbits in the frequency domain [arXiv:1308.5223](https://arxiv.org/abs/1308.5223), [arXiv:1409.4419](https://arxiv.org/abs/1409.4419)

Lorenz-gauge circular orbits in time-domain: [arXiv:gr-qc/0701069](https://arxiv.org/abs/gr-qc/0701069)  
Lorenz-gauge eccentric orbits in time-domain: [arXiv:1002.2386](https://arxiv.org/abs/1002.2386)
