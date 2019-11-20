# Introduction to Molecular Dynamics

In the FIRST Center, we primarily run atomistic simulations of bulk electrolytes (ionic liquids) and systems
with liquid-solid interfaces (channel and pore simulations).  There are several simulation engines that are
primarily suited for running these types of simulations.

## A Note on Versions

You should strive to use the most recent versions of these packages if possible.  This is due to bug fixes
and performance increases that are included in the updates.  It is NOT acceptable to be using GROMACS 4 in
2019.

## GROMACS

[GROMACS](http://www.gromacs.org) is a MD simulation package developed primarily for biomolecules.  It has
support for various parallelization schemes, including MPI and CUDA, making it faster in comparison to many
other simulation engines.  On this project, we use GROMACS for bulk simulations of electrolytes, biased free
energy simulations via
[AWH](http://manual.gromacs.org/documentation/2019-beta1/reference-manual/special/awh.html), and some
polarizable water simulations with [BK3](https://github.com/Marcello-Sega/BK3-water-model).  On Rahman, you
can expect to get between 70 to 100 ns per day depending on the system size.

GROMACS is unique to other simulation engines in that it's run on the command line.  The simulation
parameters are contained in a `.mdp` file, the atom and coordinate information in a `.gro` file and the
topology information contained in a `.top` file.  The input files for a simulation file are compiled with
`gmx grompp`, and the simulation is run with `gmx mdrun`.  There are also analysis tools within GROMACS, but
we typically try to avoid using these due to their reliance on XMGrace and the fact they're not very
reproducible.

A popular set of tutorials by Justin Lemkul for Gromacs can be viewed [here](http://www.mdtutorials.com/gmx/).
These tutorials will give you a sense of how simulations are run in GROMACS, but keep in mind that our
workflow is slightly different due to our integration with MoSDeF.

### Quirks of GROMACS
- GROMACS does not allow for mixed scaling for Lennard-Jones and Coulomb interactions
- Be aware of `comb-rules` in the `.top` file.  OPLS uses geometric mixing rules, the Lopes force field uses
  lorentz-berthelot mixing rules
- GROMACS will automatically generate cross-interactions based on the combinining rules

## LAMMPS
[LAMMPS](https://lammps.sandia.gov) is a MD simulation package developed by Sandia National Lab.  Our main
use for LAMMPS currently is to run MXene simulations.  We've also used LAMMPS to run a few grand canonial
monte carlo (GCMC) simulations.  Unlike GROMACS, LAMMPS is run with input scripts created by the user.  The
process of creating input scripts is more involved in LAMMPS, where the user is responsible for setting the
majority of the simulation parameters.  However, this means that LAMMPS is more flexible than GROMACS.  The
main drawback of LAMMPS is that it is much less computationally efficient.

I haven't found a set of LAMMPS tutorials as useful as the ones for GROMACS, but there are a good number of
examlples on their [GitHub page](https://github.com/lammps/lammps/tree/master/examples)

### Quirks of LAMMPS
- LAMMPS does not automatically generate cross-interactions.  The user has to explictly specify
  cross-interactions.
- The lammpstrj format isn't very robust for post-processing.  It's better to write your trajectories to DCD
  or XTC formats.
