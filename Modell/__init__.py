import numpy as np
np.set_printoptions(precision=6, suppress=True, linewidth=250)

# todo: currently
# use logging?
# Tests:
#   linear and modal analysis, units should be checked. 100 % test coverage?
# Lin stat:
#   internal actions in internal points, also plot them
#   calculating stresses (Sections need stress points for this...)
# Modal:
#   check added mass etc. Chopra example with higher resolution
#   tests
#   modal masses, modal participation
#   EC-based stuff: modal antwort, summation...
# Buckling:
#   implement.

# todo: LATER
# Lin stat:
#   adding loads off-node (concentrated force and moment, bending)
#   http://12000.org/my_notes/stiffness_matrix/stiffness_matrix_report.htm for internal actions, at the end of the page
#   LATER: defining hinges?
#   LATER: 3D Beam element?
# Modal:
#   ability to convert loads to masses - through defining gravity and its direction
#   LATER: ability to choose which directions mass act?

#
# off-node loads:
# - solve the element to get nodal reactions -> apply nodal reactions as loads
# - handle cases like concentrated force or moment and distributed load perpendicular to the beam axis
# - distributed load: linearly changing, acting on arbitrary part of the beam as general case
#
#
