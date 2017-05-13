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

VERTICAL = not False  # True/False for vertical/horizontal
NR_BEAMS = 2  # number of finite elements
LENGTH = 200  # length of cantilever
F_HORIZONTAL = 0
F_VERTICAL = 0


# nodes
if VERTICAL:  # vertical beam
    _nodes = [Node.Node.from_dict(adict={'ID': i+1, 'coords': (0, i*LENGTH/NR_BEAMS)}) for i in range(NR_BEAMS+1)]
else:  # horizontal beam
    _nodes = [Node.Node.from_dict(adict={'ID': i+1, 'coords': (i*LENGTH/NR_BEAMS, 0)}) for i in range(NR_BEAMS+1)]

# beams
section_column = sections.Rectangle(height=10, width=3)  # section
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
structure.add_nodal_load(nodeID=_nodes[-1].ID, dynam={'FX': F_HORIZONTAL, 'FY': F_VERTICAL, 'MZ': 0}, clear=True)
# beam internal loads
structure.add_internal_loads(beam=structure.beams[0], loadtype='uniform perpendicular force', q=2.3)
# structure.add_internal_loads(beam=structure.beams[0], loadtype='uniform perpendicular force', q=1.3)

# solving it
solve(structure, analysis='all')


# posprocessing
structure.draw(analysistype='linear static')

print(structure.results['linear static'].reaction_forces)

exit()


strudisp = structure.results['linear static'].displacement_results
print(strudisp)
print(structure.K_with_BC * strudisp)
#
for b in structure.beams:
    print('')
    print(b)
    disp = structure.results['linear static'].element_displacements(local=True, beam=b, asvector=True)
    print(b.nodal_reactions(disps=disp))
    print('local reactions')
    print(b.nodal_reactions(disps=disp)-structure.load_vector)