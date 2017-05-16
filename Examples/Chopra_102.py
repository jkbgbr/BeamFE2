# -*- coding: utf-8 -*-

from Modell import HermitianBeam_2D as HB
from Modell import Structure
from Modell import Node
import math

"""
Anil K. Chopra, Example in 10.2
Kinda OK.
"""


L = 10
mass = 100000 / 386.
rho = 0.1  # consistent mass matrix but density=0 -> only the masses are considered.
EE = 29000
I = 320
A = 63.41

_nodes = [0]  # 0th element to have mathcing ID and list position
_nodes.append(Node.Node.from_dict(adict={'ID': 1, 'coords': (0, 0)}))
_nodes.append(Node.Node.from_dict(adict={'ID': 2, 'coords': (L/2., 0)}))
_nodes.append(Node.Node.from_dict(adict={'ID': 3, 'coords': (L, 0)}))

# beams
_beams = []
_beams.append(HB.HermitianBeam2D.from_dict(adict={'ID': 1, 'E': EE, 'I': I, 'A': A, 'rho': rho, 'i': _nodes[1], 'j': _nodes[2]}))
_beams.append(HB.HermitianBeam2D.from_dict(adict={'ID': 2, 'E': EE, 'I': I, 'A': A, 'rho': rho, 'i': _nodes[2], 'j': _nodes[3]}))

# supports
BCs = {1: ['ux', 'uy', 'rotz']}  # supports as dict

# this is the cantilever itself, composed of the beams, complete with supports
structure = Structure.Structure(beams=_beams, supports=BCs)

# adding loads, masses
structure.add_nodal_mass(nodeID=2, mass={'mx': mass * L / 2., 'my': mass * L / 2.})
structure.add_nodal_mass(nodeID=3, mass={'mx': mass * L / 4., 'my': mass * L / 4.})

# solving it
structure.solver['modal'].solve()

# posprocessing
om_1 = (3.15623 * math.sqrt(EE * I / (mass * L ** 4))) / (2 * math.pi)
om_2 = (16.258 * math.sqrt(EE * I / (mass * L ** 4))) / (2 * math.pi)
print('expected / calculated values in [Hz]:')
print('Mode 1: %.3f / %.3f' % (om_1, structure.results['modal'].frequencies[0]))
print('Mode 2 is axial, not shown here')
print('Mode 3: %.3f / %.3f' % (om_2, structure.results['modal'].frequencies[2]))

for i in range(3):
    structure.draw(analysistype='modal', mode=i)
