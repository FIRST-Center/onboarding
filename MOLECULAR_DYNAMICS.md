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

## OpenMM
[OpenMM](http://openmm.org/) is a MD simulation engine largely developed by the Pande
Group at Stanford University.  Is it much younger than other engines, about 10 years old, but is
fairly mature and offers a handfkl of interesting features that make it worth consideration.  It has
a Python wrapper, so you can do everything within a Python script; minimal command-line usage ano no
writing code in an application-specific language.  It also has good performance on GPUs (on par or
better than GROMACS and roughly 5-100 times faster than LAMMPS) OpenMM supports arbitrary non-bonded
functional forms in their
[CustomNonbondForce](http://docs.openmm.org/7.0.0/api-c++/generated/OpenMM.CustomNonbondedForce.html).
In principle, you can pass it an arbitrary algebraic expression and it will handle everything under
the hood for you. Something more complex than a 12-6 LJ potential can be difficult/impossible or not
performant to do in other engines

### Quirks of OpenMM
- OpenMM is almost exclusively used by biophysicists so things more relevant to materials & fluids
  simulations (i.e. surfaces with fixed position) may be tricky in practice.
- OpenMM only supports a few file formats for
  [input](http://docs.openmm.org/latest/api-python/app.html#loaders-and-setup)/
  [output](http://docs.openmm.org/latest/api-python/app.html#reporting-output) and it also takes a
  few lines of code to go from an input to actually running a simulation.
- There seems to be a very limited selection of [thermostats and
  barostats](http://docs.openmm.org/latest/api-python/library.html#core-objects) to choose from
  (although it seems possible to generate a new one, it does not ship with some common ones such as
  the Nose-Hoover thermostat or the Parrinello-Rahman barostat)

## Cassandra
[Cassandra](https://cassandra.nd.edu/) is a Monte Carlo (MC) engine developed by Ed Maginn since the
1990s and, more recently, his group at Notre Dame. It began sometime in the 1990s, although most of
its use and development has come in the past few years. It is designed primarily for thermodynamic
studies, i.e. fluid-phase equilibria, but can also be used to do conventional MD simulations of
liquids with translational moves. It includes a number of extended ensembles (at least GCMC and
GEMC) and a number of types of MC moves. One fairly unique feature it implements is configurational
bias Monte Carlo, which is much better at generation configurations for (and subsequently relaxing)
systems of long-chain molecules than blind packing (i.e. PACKMOL/`gmx solvate`) and MD simulations.
A postdoc (Ryan DeFever) is working on integrating it into MoSDeF and (as of writing this, in
December 2019) this is complete in principle but needs some more stress-testing and development to
work well. However, it can currently be used outside of MoSDeF, although you need to do some manual
edited of input files, as is typical of most other engines.

### Quirks of Cassandra
- No GPU support
- Fairly young engine and limited user base (positive spin: you can drive feature development!)
- Limited "glue," as in it only supports a few input/output formats (although this may change with
  MoSDeF)

## HOOMD-blue
[HOOMD-blue](http://glotzerlab.engin.umich.edu/hoomd-blue/) is a particule simulation toolkit that
is most often used for MD simulations. There is also a hard particle Monte Carlo (HPMC) module but
it is not currently clear how this could be used in FIRST. Its origin is around 2008 when Josh
Anderson wrote some CUDA code to run MD simulation on GPUs, which was extremely novel at the time.
It is currently developed in the Glotzer lab at the University of Michigan by him and several
others. Its primary draws are high performance on GPUs and a clean Python API. (OpenMM and HOOMD-blue
are similar in this respect, although they are typically applied in different domains.)

### Quirks of HOOMD-blue
- The Glotzer lab and most/all of the HOOMD-blue user base does no atomistic simulations. There are
  no fundamental limitations to using it to do atomistic MD, but there may be some practical issues
  trying to use it for such systems, particularly for interesting and exotic things.
- Limited "glue," it uses a fairly specific format (GSD, or general simulation data) for input and
  output. It can be read by a Python package of the same name but may require converting into other
  formats for particular analyses. As of December 2019, MDTraj and MDAnalysis have GSD support, but
  few or no other tools do.
- It is currently going through an overhaul to version 3.x; while 2.x may be supported for a while,
  there are some API changes coming in the near future.
