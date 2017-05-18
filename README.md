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

Missing / todo:
* setup.py
* pypi
* query for the results as written in BeamFE2/__init__.py
* buckling analysis
* proper scaling of the displayed results
* moment loads need a graphical symbol
* global reaction forces should be displayed
* modal masses, modal participations
* jupyter notebooks?
* code optimization
* code for displaying needs a cleanup

Further plans:
* enable hinges at element ends

Currently the most welcome help would be to impolement a solver for the buckling analysis.
