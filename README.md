# BeamFE2
Beam Finite element suite under development with linear static and modal analysis.

The module provides the essentials for a 2D Beam FE calculation:

* 2D, 2 node 3 DOF cubic Bernoulli-Euler beam element with consistent mass matrix
* linear elastic materials
* nodal loads (force, moment) and masses
* beam internal loads
    * concentrated moment
    * concentrated force (axial or perpendicular)
    * force uniformly distributed along the length (axial or perpendicular)
* solvers for linear static and modal analysis
* about 60 tests, test coverage at cca. 90%
* no pre-processor
* some post-processing abilities: graphical display for displacements/shapes, internal actions, query of results

Requirements are:
* numpy (linear static solver), scipy (modal solver), matplotlib (display)

Further plans:
* enable hinges at element ends
* pinned ends
* temperature load
* support displacement as load

Currently the most welcome help would be to impolement a solver for the buckling analysis.
