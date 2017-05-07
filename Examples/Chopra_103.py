# -*- coding: utf-8 -*-

from Modell import HermitianBeam_2D as HB
from Modell import Structure
from Modell import Node
from solver import solve
import math

"""
Anil K. Chopra, Example in 10.3
Kinda OK.
"""


L = 10
mass = 100000000 / 386.
rho = 7850  # steel / 1000
EE = 29000
I = 320
A = 63.41

_nodes = [0]  # 0th element to have mathcing ID and list position
_nodes.append(Node.Node.from_dict(adict={'ID': 1, 'coords': (0, 0)}))
_nodes.append(Node.Node.from_dict(adict={'ID': 2, 'coords': (0, L)}))
_nodes.append(Node.Node.from_dict(adict={'ID': 3, 'coords': (L, L)}))

# beams
_beams = []
_beams.append(HB.HermitianBeam2D.from_dict(adict={'ID': 1, 'E': EE, 'I': I, 'A': A, 'rho': rho, 'i': _nodes[1], 'j': _nodes[2]}))
_beams.append(HB.HermitianBeam2D.from_dict(adict={'ID': 2, 'E': EE, 'I': I, 'A': A, 'rho': rho, 'i': _nodes[2], 'j': _nodes[3]}))

# supports
BCs = {1: ['ux', 'uy', 'rotz']}  # supports as dict

# this is the cantilever itself, composed of the beams, complete with supports
structure = Structure.Structure(beams=_beams, supports=BCs)

# adding loads, masses
structure.add_mass_to_node(nodeID=1, mass=10000000 * mass / 2., clear=True)  # clears previous loads
structure.add_mass_to_node(nodeID=2, mass=2 * mass)
structure.add_mass_to_node(nodeID=3, mass=mass)

# solving it
solve(structure, analysis='modal')

# posprocessing
om_1 = (0.6987 * math.sqrt(EE * I / (mass * L ** 3))) / (2 * math.pi)
om_2 = (1.874 * math.sqrt(EE * I / (mass * L ** 3))) / (2 * math.pi)
print('expected / calculated values in [Hz]:')
print('Mode 1: %.8f / %.8f' % (om_1, structure.results['modal'].frequencies[0]))
print('Mode 2: %.8f / %.8f' % (om_2, structure.results['modal'].frequencies[1]))

for i in range(2):
    structure.draw(analysistype='modal', mode=i)
