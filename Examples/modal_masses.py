# -*- coding: utf-8 -*-
from Modell import HermitianBeam_2D as HB
from Modell import BeamSections as sections
from Modell import Structure
from Modell import Node
from solver import solve

"""
A simple cantilever beam in vertical or horizontal position.
"""

NR_BEAMS = 3  # number of finite elements
LENGTH = 1400  # length of cantilever

# nodes
_nodes = [Node.Node.from_dict(adict={'ID': i+1, 'coords': (0, i*LENGTH/NR_BEAMS)}) for i in range(NR_BEAMS+1)]

# beams
section_column = sections.Circle(r=1)  # section
rho = 7.85e-9
EE = 2.1e5
_beams = [HB.HermitianBeam2D.from_dict(adict=
                                       {'ID': i, 'E': EE, 'I': section_column.I['x'], 'A': section_column.A,
                                        'rho': rho, 'i': _nodes[i], 'j': _nodes[i+1]})
          for i in range(NR_BEAMS)]

# supports
BCs = {1: ['ux', 'uy', 'rotz']}  # supports as dict

# this is the cantilever itself, composed of the beams, complete with supports
structure = Structure.Structure(beams=_beams, supports=BCs)

# adding loads, masses
_last = max([x.ID for x in _nodes])
structure.add_mass_to_node(nodeID=_last, mass=1./1000)
# # adding loads
# structure.add_single_dynam_to_node(nodeID=len(_nodes), dynam={'FX': F_HORIZONTAL, 'FY': F_VERTICAL}, clear=True)

# solving it
solve(structure, analysis='modal')

# posprocessing
structure.draw(analysistype='modal')

import numpy as np
alakok = 5
m = structure.M_with_masses
Sel = np.ones((1, len(m)))
unit = np.ones((1, len(m))).T

for i in range(len(m)):
    fi = structure.results['modal'].modeshape(i)
    print(fi)
    print((m*fi*(fi.T*m*unit)/(fi.T*m*fi))*Sel[i])

# m = np.matrix([[160299, 0], [0, 135214]])
# fi = np.matrix([[-0.537], [-0.843]])
# Sel = np.matrix([[1.491], [2.354]])
# unit = np.matrix([[1], [1]])
#
# # print(m)
# # print(fi)
# # print(Sel)
# # print(unit)
# #
# # print('')
# # print(fi.T)
# # print(fi.T*m*unit)
# # print(fi.T*m*fi)
# # print((fi.T*m*unit)/(fi.T*m*fi))
#
# print((m*fi*(fi.T*m*unit)/(fi.T*m*fi))*1.491)


#
# for b in structure.beams:
#     print('')
#     print(b)
#     disp = structure.results['linear static'].element_displacements(local=True, beam=b, asvector=True)
#     print(b.nodal_reactions(disps=disp))
