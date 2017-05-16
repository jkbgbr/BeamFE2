# -*- coding: utf-8 -*-

from BeamFE2 import HermitianBeam_2D as HB
from BeamFE2 import BeamSections as sections
from BeamFE2 import Structure
from BeamFE2 import Node

"""
A simple vertical cantilever beam with a horizontal endbeam.
"""

NR_BEAMS = 10  # number of finite elements in the vertical part
LENGTH = 1400  # height of cantilever, mm
LENGTH_ENDBEAM = 300  # legth of horizontal beam at the top of the cantilever, mm
F_HORIZONTAL = 100
F_VERTICAL = 100

# nodes
_nodes = [Node.Node.from_dict(adict={'ID': i+1, 'coords': (0, i*LENGTH/NR_BEAMS)}) for i in range(NR_BEAMS+1)]
_nodes.append(Node.Node.from_dict(adict={'ID': len(_nodes)+1, 'coords': (LENGTH_ENDBEAM, LENGTH)}))

# beams
section_column = sections.Circle(r=55)  # section
rho = 7.850e-8
EE = 2.1e6
_beams = [HB.HermitianBeam2D.from_dict(adict={'ID': i, 'E': EE, 'I': section_column.I['x'], 'A': section_column.A, 'rho': rho, 'i': _nodes[i], 'j': _nodes[i+1]}) for i in range(NR_BEAMS)]
_beams.append(HB.HermitianBeam2D.from_dict(adict={'ID': len(_beams)+1, 'E': EE, 'I': section_column.I['x'], 'A': section_column.A, 'rho': rho, 'i': _nodes[-2], 'j': _nodes[-1]}))

# supports
BCs = {1: ['ux', 'uy', 'rotz']}  # supports as dict

# this is the cantilever itself, composed of the beams, complete with supports
structure = Structure.Structure(beams=_beams, supports=BCs)

# adding loads
structure.add_single_dynam_to_node(nodeID=len(_nodes), dynam={'FX': F_HORIZONTAL, 'FY': F_VERTICAL}, clear=True)  # clears previous loads

# solving it
structure.solver['linear static'].solve()
structure.solver['modal'].solve()
# posprocessing

structure.draw(analysistype='linear static', internal_action='axial')
structure.draw(analysistype='linear static', internal_action='shear')
structure.draw(analysistype='linear static', internal_action='moment')
for i in range(10):
    structure.draw(analysistype='modal', mode=i)
