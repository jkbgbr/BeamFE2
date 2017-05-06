# -*- coding: utf-8 -*-

from Modell import HermitianBeam_2D as HB
from Modell import BeamSections as sections
from Modell import Structure
from Modell import Node
from solver import solve

VERTICAL = False

# nodes
_pieces = 3  # No. of element divisions
_length = 1400  # length of element [mm]

# nodes
if VERTICAL:  # vertical beam
    _nodes = [Node.Node.from_dict(adict={'ID': i+1, 'coords': (0, i*_length/_pieces)}) for i in range(_pieces+1)]
else:  # horizontal beam
    _nodes = [Node.Node.from_dict(adict={'ID': i+1, 'coords': (i*_length/_pieces, 0)}) for i in range(_pieces+1)]

# beams
section_column = sections.Circle(r=55)  # section
rho = 7850000
EE = 2.1e11
_beams = [HB.HermitianBeam2D.from_dict(adict=
                                       {'ID': i, 'E': EE, 'I': section_column.I['x'], 'A': section_column.A,
                                        'rho': rho, 'i': _nodes[i], 'j': _nodes[i+1]})
          for i in range(_pieces)]

# supports
BCs = {1: ['ux', 'uy', 'rotz']}  # supports as dict

# this is the cantilever itself, composed of the beams, complete with supports
structure = Structure.Structure(beams=_beams, supports=BCs)

# adding loads
structure.add_single_dynam_to_node(nodeID=len(_nodes), dynam={'FX': 1000000000000, 'FY': 0}, clear=True)  # clears previous loads

# solving it
solve(structure, analysis='all')
# posprocessing

structure.draw(analysistype='linear static')
for i in range(10):
    structure.draw(analysistype='modal', mode=i)
