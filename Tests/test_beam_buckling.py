# import unittest
# from Modell import HermitianBeam_2D as HB
# from Modell import BeamSections as sections
# from Modell import Structure
# from Modell import Node
# from Modell import Material
# from solver import solve
# import math
#
#
# F_HORIZONTAL = 1000000
# F_VERTICAL = 1
#
# class Hermitian2D_Model(unittest.TestCase):
#     """
#     tests for the buckling analysis of the FE Modell
#     """
#
#     @classmethod
#     def setUpClass(cls):
#         steel = Material.Steel()
#         cls.rho = steel.rho
#         cls.E = steel.E
#
#     def test_clamped_rod(self):
#         """
#         A beam with clamped end under self weight.
#         Results are accepted if circular frequencies are within 1% of the analytical solution.
#         http://iitg.vlab.co.in/?sub=62&brch=175&sim=1080&cnt=1
#         """
#         _pieces = 2  # number of finite elements
#         _length = 450.  # mm
#         rho = 7.850e-8
#         section_column = sections.Rectangle(height=1., width=1.)
#         E = 2.1e5
#
#         # nodes
#         _nodes = []
#         for i in range(_pieces + 1):
#             _nodes.append(Node.Node.from_dict(adict={'ID': i + 1, 'coords': (0, i * _length / _pieces)}))
#
#         # beams
#         I = section_column.I['x']
#         A = section_column.A
#         I = 1
#         _beams = []
#         for i in range(_pieces):
#             _beams.append(HB.HermitianBeam2D.from_dict(adict={'ID': i, 'E': E, 'I': I,
#                                                               'A': A, 'rho': rho,
#                                                               'i': _nodes[i], 'j': _nodes[i + 1]}))
#
#         # supports
#         _last = max([x.ID for x in _nodes])
#         BCs = {1: ['ux', 'uy', 'rotz']}  # supports as dict
#
#         # this is the structure
#         structure = Structure.Structure(beams=_beams, supports=BCs)
#
#         # adding loads
#         structure.add_single_dynam_to_node(nodeID=_last, dynam={'FX': 0, 'FY': -1*F_VERTICAL}, clear=True)
#
#         # solver :-) whatever happens here is done by numpy.
#
#         solve(structure, analysis='linear static')
#         solve(structure, analysis='buckling')
#
