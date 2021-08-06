# Analysis Tools and Methods

## MDTraj
Much of the analysis code uses the classes provided in MDTraj.  MDTraj offers analysis functions
to compute the distances between atoms, as well as the radial distribution function (RDF) of a
system.  MDTraj can also be used as the backend to perform other analyses such as number densities
and charge densities across a pore.  See the SOFTWARE section for a further overview of MDTraj.

## Structure Factor
The static structure factor, S(q), is an analysis that can analyzed both by experiment and theory.
Functions are provided in [scattering](https://github.com/mattwthompson/scattering) to calculate
the structure factor by performing a fourier transform of the RDF.

## Van Hove Function
The Van Hove Function is defined as the probability distribution of particle pairs as a function
of both distance and time.  In 2017, our [experimental collaborators](https://advances.sciencemag.org/content/3/12/e1603079.abstract) were able to compute the Van
Hove function of water for the first time through inelastic x-ray scattering.
To help analyze these experimental efforts, molecular simulation can be performed alongside to
study the specific pairwise interactions in detail.

To compute the Van Hove function, we have developed a set of functions contained in the
[scattering](https://github.com/mattwthompson/scattering) GitHub repository.
This package takes advantage of the MDTraj `trajectory` object to keep track of atom informtion.

## Transport Properties

### Self-Diffusion Coefficients
Self-diffusion coefficients are often computed to characterize the dynamics of a fluid system.
This is done by computing the mean-squared displacement of particles. 

NOTE: MD Simulations often have periodic-boundary conditions applied.  Before computing the
mean-squared displacements, the trajectory coordinates must be "unwrapped", which reverses the
effects of the periodic boundary conditions.  In Gromacs, this is done with the `gmx trjconv -pbc
nojump` command.

Gromacs also contains a built-in
[tool](https://manual.gromacs.org/archive/5.0.7/programs/gmx-msd.html) to compute the mean-squared
displacement.

### Conductivity
Conductivity can be computed several ways using molecular dynamics.  Perhaps the easiest method is
through the use of the Nernst-Einstein equation, which provides an ideal ionic conductivity
through the use of self-diffusion coefficients.  Nernst-Einstein assumes all ions are
uncorrelated.

To achieve a more realistic conductivity value, the Einstein-Helfand or Green-Kubo approaches can
be used.  Rather than using self-diffusion coefficients, Einstein-Helfand computes conductivity
through the use of the mean-translational dipole moment.  Green-Kubo computes conductivity through
the current autocorrelation function and requires the use particle velocities rather than
coordinates.  Please see publications below for mathematic descriptions and simulation details for
these calculations.

- Nernst-Einstein: https://doi.org/10.1021/acs.jpcb.8b11527
- Green-Kubo: https://dx.doi.org/10.1021/acs.macromol.0c02001
- Einstein-Helfand: https://aip.scitation.org/doi/10.1063/1.2868752

