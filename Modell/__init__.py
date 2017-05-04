import numpy as np
np.set_printoptions(precision=6, suppress=True, linewidth=250)




# todo:
# Lin stat:
#   adding loads off-node (concentrated force and moment, bending)
#   http://12000.org/my_notes/stiffness_matrix/stiffness_matrix_report.htm for internal actions, at the end of the page
#   calculating stresses (Sections need stress points for this...)
#   defining hinges?
#   3D Beam element?
# Modal:
#   tests, ideally opensess-based and with analytical formulae
#   (http://opensees.berkeley.edu/wiki/index.php/Eigen_analysis_of_a_two-story_shear_frame)
#   lumped mass matrix?
#   plot mode shapes (as a matter of fact, plotting should be re-thought)
#   modal masses, modal participation
#   EC-based stuff: modal antwort, summation...
#   ability to add mass
#   ability to convert loads to masses
#   ability to choose which directions mass act?
# Buckling:
#   implement. the geometrical stiffness matrix is needed for that.
# Displacements:
#   displacements are to be stored in a way that multiple shapes can be handled (linstat: 1 list, all other: multiple lists)

# make results more egyseges for the analyses.
