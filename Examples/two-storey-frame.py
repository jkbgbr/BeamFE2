# -*- coding: utf-8 -*-

from Modell import HermitianBeam_2D as HB
from Modell import BeamSections as sections
from Modell import Structure
from Modell import Node
from solver import solve

h = 2000  # story heigth
L = 2 * h  # bay width

_nodes = [0]  # 0th element to have mathcing ID and list position
_nodes.append(Node.Node.from_dict(adict={'ID': 1, 'coords': (0, 0)}))
_nodes.append(Node.Node.from_dict(adict={'ID': 2, 'coords': (L, 0)}))
_nodes.append(Node.Node.from_dict(adict={'ID': 3, 'coords': (0, h)}))
_nodes.append(Node.Node.from_dict(adict={'ID': 4, 'coords': (L, h)}))
_nodes.append(Node.Node.from_dict(adict={'ID': 5, 'coords': (0, 2 * h)}))
_nodes.append(Node.Node.from_dict(adict={'ID': 6, 'coords': (L, 2 * h)}))

# beams
section_column_1 = sections.Recangle(height=100, width=100)  # section, 1st storey
section_column_2 = sections.Recangle(height=200, width=100)  # section, 2nd storey
section_beam = sections.Recangle(height=1e3, width=1e3)  # section, 2nd storey
rho = 7850000
EE = 2.1e11

_beams = []
_beams.append(HB.HermitianBeam2D.from_dict(adict={'ID': 1, 'E': EE, 'I': section_column_1.I['x'], 'A': section_column_1.A, 'rho': rho, 'i': _nodes[1], 'j': _nodes[3]}))
_beams.append(HB.HermitianBeam2D.from_dict(adict={'ID': 2, 'E': EE, 'I': section_column_2.I['x'], 'A': section_column_2.A, 'rho': rho, 'i': _nodes[3], 'j': _nodes[5]}))

_beams.append(HB.HermitianBeam2D.from_dict(adict={'ID': 3, 'E': EE, 'I': section_column_1.I['x'], 'A': section_column_1.A, 'rho': rho, 'i': _nodes[2], 'j': _nodes[4]}))
_beams.append(HB.HermitianBeam2D.from_dict(adict={'ID': 4, 'E': EE, 'I': section_column_2.I['x'], 'A': section_column_2.A, 'rho': rho, 'i': _nodes[4], 'j': _nodes[6]}))

_beams.append(HB.HermitianBeam2D.from_dict(adict={'ID': 5, 'E': EE, 'I': section_beam.I['x'], 'A': section_beam.A, 'rho': rho, 'i': _nodes[3], 'j': _nodes[4]}))
_beams.append(HB.HermitianBeam2D.from_dict(adict={'ID': 6, 'E': EE, 'I': section_beam.I['x'], 'A': section_beam.A, 'rho': rho, 'i': _nodes[5], 'j': _nodes[6]}))

# supports
BCs = {1: ['ux', 'uy', 'rotz'], 2: ['ux', 'uy', 'rotz']}  # supports as dict

# this is the cantilever itself, composed of the beams, complete with supports
structure = Structure.Structure(beams=_beams, supports=BCs)

# adding loads
structure.add_single_dynam_to_node(nodeID=5, dynam={'FX': 1000000, 'FY': 0}, clear=True)  # clears previous loads

# # adding loads
# mass_1 = 2000  # kg
# mass_2 = 2 * mass_1  # kg
# structure.add_mass_to_node(nodeID=3, mass=mass_1 / 2, clear=True)  # clears previous loads
# structure.add_mass_to_node(nodeID=4, mass=mass_1 / 2, clear=True)  # clears previous loads
# structure.add_mass_to_node(nodeID=5, mass=mass_2 / 2, clear=True)  # clears previous loads
# structure.add_mass_to_node(nodeID=6, mass=mass_2 / 2, clear=True)  # clears previous loads

# solving it
solve(structure, analysis='all')
# posprocessing

structure.draw(analysistype='linear static')
for i in range(10):
    structure.draw(analysistype='modal', mode=i)
