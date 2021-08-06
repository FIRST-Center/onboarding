# Forcefields

Forcefields are the set of parameters that define the interactions between atoms in your
simulation.  Here, a list of common forcefields and sanity checks is provided.

## Useful Forcefields

### Ionic Liquids
The most widely used force field for ionic liquids is
[CL&P](https://link.springer.com/article/10.1007/s00214-012-1129-7).  This set of parameters is
based off of the OPLS-AA set of parameters and has been further tuned for ionic liquids.  These
parameters have been tuned to match densities and energies, and has been widely reported to
underestimte dynamics.  A common fix for this is to scale the partial charges described
[here](https://pubs.rsc.org/en/content/articlelanding/2012/CP/c2cp23329k). 

The [KPL](10.1002/cphc.200700552) force field by Koddermann, Paschek, and Ludwig is a further
refinement of CL&P for ionic liquids containing immidazolium cations and bistriflate anions.
Non-bonded interactions have been modified to achieve more accurate dynamics for this specific
ionic liquid.

### Solvents
For solvents and other small molecules, the OPLS-AA set of parameters is a good first force field
to try.
OPLS-AA is natively implemented in [foyer](https://github.com/mosdef-hub/foyer.git) as an XML
file.
The General Amber Force Field [(GAFF)](https://github.com/rsdefever/antefoyer) is another good
option to try.  A foyer implementation has been developed by Ryan DeFever.

### Water
The two most used water models are [SPC/E](https://pubs.acs.org/doi/10.1021/j100308a038) and [TIP\_/3P](https://aip.scitation.org/doi/10.1063/1.445869).
These are 3-site rigid water models with point charges and Lennard-Jones interactions.

[BK3](http://marcello-sega.github.io/BK3-water-model/) is a polarizable model that uses Drude
Oscillators to model dynamic charges.

### Ions
The parameters by [Joung-Cheatham](https://pubs.acs.org/doi/10.1021/jp8001614) is widely used for alkali and halide ions.  Specific parameters are provided for specific models.

## Checklist
When applying a force field, it is advised to check that the following parameters are correct:
- Combining rule
- 1-4 Scaling Factors for Electrostatic and non-bonded interactions
- Functional forms defined by force fields are consistent with the functional forms implemented by
  the simulation engine.
