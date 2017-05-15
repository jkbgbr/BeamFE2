import numpy as np
# np.set_printoptions(precision=6, suppress=True, linewidth=250)

"""
Units to be used to have correct results in the modal analysis. Examples only.
Length: [m]
Mass: [kg] -> Density: [kg/m3]
Young's module: [N/m2] 
Time: [s]

OR the engineering approach
Length: [mm]
Mass: [t] -> Density: [t/mm3]!
Young's module: [N/mm2, MPa]
Time: [s]
"""

# todo: currently
# todo: enable drawing the structure without having an analysis performed previously
# todo: BCs -> support objects
# todo: re-think plotting/drawing

# use logging?
# Tests:
#   linear and modal analysis, units should be checked. 100 % test coverage?
# Lin stat:
#   calculating stresses (Sections need stress points for this...)
# Modal:
#   check added mass etc. Chopra example with higher resolution
#   tests
#   modal masses, modal participation
#   EC-based stuff: modal antwort, summation...
# Buckling:
#   implement.
# draw:
#   beautify: scaling all drawn elements.
#   draw with subplots for lin stat
#   draw global reaction forces
#   draw moment loads


# todo: LATER
# Lin stat:
#   LATER: defining hinges?
#   LATER: 3D Beam element?
# Modal:
#   ability to convert loads to masses - through defining gravity and its direction
#   LATER: ability to choose which directions mass act?
