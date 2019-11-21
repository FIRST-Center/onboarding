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
study of ionic liquid and organic solvent mixtures. Oftentimes as grad students, we are not organized 
with our file
management and we do things like inconsistently name files, move files around, etc. This is exactly what
signac can prevent.  Signac is designed to help the user store, generate, and analyze data for computational
studies.  In signac, a similarly structured data set is called a `project` and the
elements of the project's data space are called `jobs`.  In our screening project for example, our project
was called `il_solvent` because we were studying various solvated ionic liquids, and each job had a unique
solvent, and ionic liquid concentration.

The signac framework also includes a tool called
[signac-flow](https://signac-flow.readthedocs.io/en/v0.6.2/) which is useful for executing the various data
space operations in your project.  Signac-flow also makes it easy to submit operations to computer clusters,
which many of us are using to submit jobs on Rahman, NERSC, and previously Titan.  In signac-flow, the
operations of your project go into a `project.py` folder.  Some examples include initialization, sampling,
and calculating the mean squared displacement.  Some examples can be viewed
[here](https://github.com/csadorf/signac-examples/tree/master/projects).  To execute an operation in
`project.py` the basic syntax is like this: `python project.py [operation]`.  One of the more useful command
in signac-flow is `python project.py status`, which gives you an overview of the operations you've run, and
the operations you have left to run.  This is one main benefit of using signac, as you always have an idea
on the progress of your project.  To submit a job to the cluster, the command is `python project.py submit
-o operation`.  

This is just a brief overview of signac.  I pretty much use it for all of my projects now, and I recommend
learning it if you're a new grad student.
