import numpy as np
np.set_printoptions(precision=6, suppress=True, linewidth=250)

# todo:
# Tests:
#   tests for internal displacements due to nodal loading
#   linear and modal analysis, units should be checked. use pint?
# Lin stat:
#   adding loads off-node (concentrated force and moment, bending)
#   http://12000.org/my_notes/stiffness_matrix/stiffness_matrix_report.htm for internal actions, at the end of the page
#   calculating stresses (Sections need stress points for this...)
#   LATER: defining hinges?
#   LATER: 3D Beam element?
# Modal:
#   change the solver to use elimination instead of added masses
#   tests, ideally opensess-based and with analytical formulae -> chopra book, example around 10.5
#   (http://opensees.berkeley.edu/wiki/index.php/Eigen_analysis_of_a_two-story_shear_frame)
#   lumped mass matrix?
#   modal masses, modal participation
#   EC-based stuff: modal antwort, summation...
#   ability to convert loads to masses - through defining gravity and its direction
#   LATER: ability to choose which directions mass act?
# Buckling:
#   implement. the geometrical stiffness matrix is needed for that.
