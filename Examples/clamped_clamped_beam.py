# -*- coding: utf-8 -*-

from Modell import HermitianBeam_2D as HB
from Modell import BeamSections as sections
from Modell import Material
from Modell import Structure
from Modell import Node
from solver import solve

"""
A simple cantilever beam in vertical or horizontal position.
"""

VERTICAL = False  # True/False for vertical/horizontal
NR_BEAMS = 2  # number of finite elements
LENGTH = 200  # length of cantilever
F_HORIZONTAL = 0
F_VERTICAL = -1000


# nodes
if VERTICAL:  # vertical beam
    _nodes = [Node.Node.from_dict(adict={'ID': i+1, 'coords': (0, i*LENGTH/NR_BEAMS)}) for i in range(NR_BEAMS+1)]
else:  # horizontal beam
    _nodes = [Node.Node.from_dict(adict={'ID': i+1, 'coords': (i*LENGTH/NR_BEAMS, 0)}) for i in range(NR_BEAMS+1)]

# beams
section_column = sections.Rectangle(height=3, width=10)  # section
mat = Material.Steel()
rho = mat.rho
EE = mat.E
_beams = [HB.HermitianBeam2D.from_dict(adict=
                                       {'ID': i, 'E': EE, 'I': section_column.I['x'], 'A': section_column.A,
                                        'rho': rho, 'i': _nodes[i], 'j': _nodes[i+1]})
          for i in range(NR_BEAMS)]




# supports
BCs = {1: ['ux', 'uy', 'rotz'], 3: ['ux', 'uy', 'rotz']}  # supports as dict

# this is the cantilever itself, composed of the beams, complete with supports
structure = Structure.Structure(beams=_beams, supports=BCs)

# adding loads
# directly defined nodal loads
# structure.add_single_dynam_to_node(nodeID=4, dynam={'FX': F_HORIZONTAL, 'FY': F_VERTICAL}, clear=True)
# beam internal loads
for b in structure.beams:
    structure.add_internal_loads(beam=b, loadtype='uniform perpendicular force', q=23.00)

print(structure.load_vector)

# solving it
solve(structure, analysis='all')


# posprocessing
structure.draw(analysistype='linear static')
# for i in range(3):
#     structure.draw(analysistype='modal', mode=i)

print(structure.results['linear static'].displacement_results)

#
# for b in structure.beams:
#     print('')
#     print(b)
#     disp = structure.results['linear static'].element_displacements(local=True, beam=b, asvector=True)
#     print(disp)
#     print(b.nodal_reactions(disps=disp))
#     # print('local reactions')
#     # print(b.nodal_reactions(disps=disp)-b.load_vector)