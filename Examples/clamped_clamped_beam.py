# -*- coding: utf-8 -*-

from Modell import HermitianBeam_2D as HB
from Modell import BeamSections as sections
from Modell import Material
from Modell import Structure
from Modell import Node
from solver import solve
from drawing import _plotting_available, plt

"""
A simple cantilever beam in vertical or horizontal position.
"""

VERTICAL = not True  # True/False for vertical/horizontal
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
_beams = [HB.HermitianBeam2D.from_dict(adict=
                                       {'ID': i, 'E': mat.E, 'I': section_column.I['x'], 'A': section_column.A,
                                        'rho': mat.rho, 'i': _nodes[i], 'j': _nodes[i+1]})
          for i in range(NR_BEAMS)]

# supports
BCs = {1: ['ux', 'uy', 'rotz'], NR_BEAMS+1: ['ux', 'uy', 'rotz']}  # supports as dict

# this is the cantilever itself, composed of the beams, complete with supports
structure = Structure.Structure(beams=_beams, supports=BCs)

# adding loads
# directly defined nodal loads
structure.add_nodal_load(nodeID=2, dynam={'FX': F_HORIZONTAL, 'FY': F_VERTICAL}, clear=True)
# structure.add_nodal_load(nodeID=3, dynam={'FX': -F_HORIZONTAL, 'FY': -F_VERTICAL})

# # beam internal loads
# for b in structure.beams:
#     pass
#     structure.add_internal_loads(beam=b, loadtype='uniform perpendicular force', value=-.5)
# structure.add_internal_loads(beam=structure.beams[2], loadtype='concentrated perpendicular force', value=-30.00, position=0.2)
# structure.add_internal_loads(beam=structure.beams[5], loadtype='concentrated perpendicular force', value=60.00, position=0.2)
# structure.add_internal_loads(beam=structure.beams[0], loadtype='concentrated moment', value=500, position=0.3)
# structure.add_internal_loads(beam=structure.beams[5], loadtype='concentrated moment', value=-500, position=0.8)

# solving it
solve(structure, analysis='linear static')

# posprocessing
structure.draw(analysistype='linear static')
structure.draw(analysistype='linear static', internal_action='moment')
structure.draw(analysistype='linear static', internal_action='shear')
# for i in range(3):
#     structure.draw(analysistype='modal', mode=i)

print(structure.results['linear static'].reaction_forces)
