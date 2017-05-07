# -*- coding: utf-8 -*-

from Modell import HermitianBeam_2D as HB
from Modell import BeamSections as sections
from Modell import Structure
from Modell import Node
from solver import solve

"""
A simple cantilever beam in vertical or horizontal position.
"""

VERTICAL = False  # True/False for vertical/horizontal
NR_BEAMS = 3  # number of finite elements
LENGTH = 1400  # length of cantilever
F_HORIZONTAL = 1000000
F_VERTICAL = 1000000


# nodes
if VERTICAL:  # vertical beam
    _nodes = [Node.Node.from_dict(adict={'ID': i+1, 'coords': (0, i*LENGTH/NR_BEAMS)}) for i in range(NR_BEAMS+1)]
else:  # horizontal beam
    _nodes = [Node.Node.from_dict(adict={'ID': i+1, 'coords': (i*LENGTH/NR_BEAMS, 0)}) for i in range(NR_BEAMS+1)]

# beams
section_column = sections.Circle(r=55)  # section
rho = 7850000
EE = 2.1e11
_beams = [HB.HermitianBeam2D.from_dict(adict=
                                       {'ID': i, 'E': EE, 'I': section_column.I['x'], 'A': section_column.A,
                                        'rho': rho, 'i': _nodes[i], 'j': _nodes[i+1]})
          for i in range(NR_BEAMS)]

# supports
BCs = {1: ['ux', 'uy', 'rotz']}  # supports as dict

# this is the cantilever itself, composed of the beams, complete with supports
structure = Structure.Structure(beams=_beams, supports=BCs)

# adding loads
structure.add_single_dynam_to_node(nodeID=len(_nodes), dynam={'FX': F_HORIZONTAL, 'FY': F_VERTICAL}, clear=True)

# solving it
solve(structure, analysis='all')

# posprocessing
structure.draw(analysistype='linear static')
# for i in range(3):
#     structure.draw(analysistype='modal', mode=i)

for b in structure.beams:
    print('')
    print(b)
    disp = structure.results['linear static'].element_displacements(local=True, beam=b, asvector=True)
    print(b.nodal_reactions(disps=disp))