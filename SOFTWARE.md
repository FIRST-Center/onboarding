# Software

Besides MoSDeF, there are several main software packages that we use in our research.
These packages are all open-source, so I recommend checking out the source codes on GitHub.

## MDTraj

[MDTraj](http://mdtraj.org/1.9.3/) is an analysis library for molecular dynamics trajectories.
The main data structre in MDTraj is the `trajectory`, which is loaded in with the command: 
`mdtraj.load(trajectory, top=topology)`.  A trajectory can then be used for all types of analyses, including
mean squared displacements and radial distribution functions.  MDTraj trajectories are also the data
structure used for analyses in [pairing](https://github.com/mattwthompson/pairing) and 
[scattering](https://github.com/mattwthompson/scattering).  

Trajectories can also be sliced by specific atom selections.  When we compute diffusivities of ionic
liquids, we will often slice the trajectory separately by cations and anions to get separate diffusivities.
We can also slice trajectories by time, which is useful when we only want to analyze a chunk of the
trajectory.

MDTraj is especially good for analyzing GROMACS trajectories.  To do so, you can load in a gro and top file
as a trajectory.  If you need to perform an analysis that requires bonding information, you can load in a
MOL2 file.  MDTraj also handles LAMMPS trajectories, although its slightly more difficult.  Because of this,
it is recommended you also save out GRO and MOL2 files when you generate LAMMPS data files.  That way, you
can use these as structure files when loading in a file format such as `lammpstrj`.

One final note on MDTraj is that it's not as actively developed as it used to be and you may end up finding pieces
of code that have bugs in it.  I personally found several bugs in the MOL2 reader.  If you happen to find a
bug, go ahead and fix it with a pull request on GitHub.  There are still several developers on this project
that will go ahead and merge your pull request as long as the tests are passing.

## MDAnalysis

[MDAnalysis](https://www.mdanalysis.org/about/) is another analysis library for molecular dynamics
trajectories.  MDAnalysis is a [NumFOCUS](https://www.numfocus.org/affliated-projects.html) project that is
more actively developed compared to MDTraj.  Despite this, I don't have as much experience using this
package.

The main data structure of MDAnalysis is the `Universe`.  `Universe` objects can be instantiated
with the following command `MDAnalysis.Universe(topology, trajectory)`.

One great feature of MDAnalysis is that they support a large variety of topology file formats which can be
viewed [here](https://www.mdanalysis.org/docs/documentation_pages/coordinates/init.html#id2).  The few times
I've used MDAnalysis was because I was attempting to load in a topology that wasn't supported by MDTraj.

## ParmEd

## Signac

The [signac framework](https://signac.io) is a workflow management package developed by the Glotzer Group at
the University of Michigan.  We have been using signac since 2017, and was a huge component of our screening
study of ionic liquid and organic solvent mixtures.  
