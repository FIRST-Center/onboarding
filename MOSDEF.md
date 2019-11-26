# Molecular Simulation and Design Framework (MoSDeF)

Much of the FIRST-related work in our lab is well-suited for MoSDeF since we primarily run atomistic MD
simulations.  Within the last few years, we have heavily incorporated MoSDeF into our various research workflows.
For example, we used mBuild and foyer to paramterize and initialize 400+ systems for a ionic liquid - organic
solvent screening study.  Additionally we have developed the `Pore-builder` mBuild recipe to initialize
graphene pore systems and the `MXenes` package to initialize MXene systems.  Below I will outline a general
workflow that we use to initialize systems with MOSDEF

## Loading in molecules into mBuild

In mBuild, compounds can be created in three ways:

1. Develop compound classes with particles, ports, and connections
2. Specify a SMILES string
3. Load in a molecule from a structure file (PDB, MOL2)

For fluids, we primarily load in structure files to create mBuild compounds.  These structure files are created
with [Avogadro](https://avogadro.cc).  With Avogadro, you can draw the desired molecule and save it out to
desired file format.  We typically save these molecules to MOL2 files since they are easier to work with
compared to PDB files.

For solid materials, we develop compounds classes.  Examples of this again include `Pore-Builder` and
`MXenes`.  If you need to create a new solid material, you should use the mBuild Lattice builder to do so.

To create our system with our various mBuild compounds, we will use `mbuild.fill_box` to fill a simulation
box.

## Apply Force fields with foyer

Once we create our system as an mBuild compound, we can now use foyer to atomtype and apply the force field
paramters.  The force field information in foyer is stored in XML files.  Matt Thompson has created and
compiled a list of popular atomistic force fields for ionic liquids titled `il_forcefields`.  For alkali and
halide ions, we typically either use the force field by
[Dang](https://aip.scitation.org/doi/abs/10.1063/1.459714) or by [Joung and
Cheatham](https://doi.org/10.1021/jp8001614).  The XML files for these ion force fields can be found
[here](https://github.com/mattwthompson/aqueous-electrolytes).  There are a variety of models that can
be used to parameterize water, although we most often use SPC/E or TIP3P which are also contained in the
aqueous-electrolytes repository.  Please note that the paramters for Joung and Cheatham are slightly
different for TIP3P and SPCE, so make sure these force fields match.  Other solvents that we simulate can
often be paramterized with the OPLS-AA Force Field.

For ionic liquid force fields, we often have to scale the partial charges by a factor of 0.6-0.8 to get good
agreement of diffusivity with experimental results.  A paper on this method can be viewed
[here](https://pubs.rsc.org/en/content/articlelanding/2015/cp/c4cp05550k#!divAbstract).  You may have to
experiment with partial charge scaling factors to yield the correct ionic liquid diffusivities.

### Applying multiple force fields

In general, its not advisable to mix and match force fields.  You should take proper precaution when doing
this, and verify that the same combining rules, scaling factors, etc. are used.

The ionic liquid force fields contained in Matt's force field repository are compatible with OPLS, so its
okay in this case to combine these force fields.  Currently there isn't a great way to apply multiple force
fields to an mBuild compound.  To do this, we have to create separate mBuild compounds to apply the separate
force fields to.  We then loop through the children in the original mBuild compound and add these children
to the new, separate mBuild compounds.  Then, we can separately add the force fields to the individual
compounds.  When force fields are applied to compounds, a parameterized
[ParmEd](https://parmed.github.io/ParmEd/html/index.html) structure is returned.  ParmEd structure can be
added together, which is what we do with the individual ParmEd structure.

For more information on these steps, check out the mBuild and foyer documentation.

### Tutorials

You can learn the basics of MoSDeF by going through the
[tutorials](https://github.com/mosdef-hub/mosdef_tutorials).  In particular it would
be usefull to go through the `overview.ipynb` notebook.

### MoSDeF Workflows

To get a better idea of how MoSDeF is used for research, take a look at
[mosdef-workflows](https://github.com/mosdef-hub/mosdef-workflows).  This is a
collection of notebooks with sample workflows.  In particular, take a look at
`graphene-pore`, which goes through building and simulating a graphene slit pore
system with water, and calculating a number density profile of water within the slit
pores.
