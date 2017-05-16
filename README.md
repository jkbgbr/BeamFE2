# BeamFE2
Beam Finite element suite under development with linear static and modal analysis.

The module provides the essentials for a 2D Beam FE calculation:

* 2D cubic beam element (Bernoulli-Euler assumptions) with 3 DOF each node
* linear elastic materials
* nodal loads (force, moment) and masses
* beam internal loads
    * concentrated moment
    * concentrated force (axial or perpendicular)
    * force uniformly distributed along the length (axial or perpendicular)
* solvers
* graphical result display for displacements, internal actions, mode shapes etc.
* numerical results are easily accessible
* about 55 tests, test coverage at cca. 85%
* no pre-processor

Requirements are:
* numpy (linear static solver), scipy (modal solver), matplotlib (display)

Missing / todo:
* setup.py
* checking the signs of internal loads
* signs of reaction moments are not correct
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
