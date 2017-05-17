# -*- coding: utf-8 -*-

from BeamFE2 import HermitianBeam_2D as HB
from BeamFE2 import BeamSections as sections
from BeamFE2 import Material
from BeamFE2 import Structure
from BeamFE2 import Node


"""
A simple cantilever beam in vertical or horizontal position.
"""

VERTICAL = False  # True/False for vertical/horizontal
NR_BEAMS = 1  # number of finite elements
LENGTH = 200  # length of cantilever
F_HORIZONTAL = 0
F_VERTICAL = 0


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
BCs = {1: ['ux', 'uy', 'rotz']}  # supports as dict

# this is the cantilever itself, composed of the beams, complete with supports
structure = Structure.Structure(beams=_beams, supports=BCs)

# adding loads
# directly defined nodal loads
# structure.add_nodal_load(nodeID=_nodes[-1].ID, dynam={'FX': -100, 'FY': -200, 'MZ': 0}, clear=True)
# beam internal loads
# structure.add_internal_loads(beam=structure.beams[0], loadtype='uniform axial force', value=1)
# structure.add_internal_loads(beam=structure.beams[0], loadtype='concentrated axial force', value=200, position=0.3)
# structure.add_internal_loads(beam=structure.beams[0], loadtype='concentrated moment', value=-1000000, position=0.5)
structure.add_internal_loads(beam=structure.beams[0], loadtype='uniform perpendicular force', value=-10)

# solving it
structure.solver['linear static'].solve()


# posprocessing
# structure.draw(analysistype='linear static')
structure.draw(analysistype='linear static', internal_action='axial')
# structure.draw(analysistype='linear static', internal_action='shear')
structure.draw(analysistype='linear static', internal_action='moment')

print(structure.results['linear static'].reaction_forces)
print(structure.results['linear static'].displacement_results)

for beam in structure.beams:
    disp = structure.results['linear static'].element_displacements(local=True, beam=beam, asvector=True)
    print(beam.internal_action(disp=disp, action='moment', pos=0.5))
    print(beam.internal_action(disp=disp, action='moment', pos=.99999))
    print(beam.internal_action_distribution(action='axial'))
    print(beam.internal_action_distribution(action='moment'))