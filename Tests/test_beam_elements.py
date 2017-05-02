import unittest
from Beams import HermitianBeam as HB
from Beams import Sections as sections

from drawing import draw_beam
import numpy as np
import math
from Tests import ATOL


class Hermitian2D(unittest.TestCase):
    """
    tests for the hermitian 2D beam
    """

    @classmethod
    def setUpClass(cls):
        cls.EE = 1.
        cls.nodes = {1: HB.Node(ID=1, coords=(0, 0)), 2: HB.Node(ID=2, coords=(1, 0)), 3: HB.Node(ID=3, coords=(2, 0))}
        cls.beam1 = HB.HermitianBeam2D(E=cls.EE, I=1., A=1., i=cls.nodes[1], j=cls.nodes[2])
        cls.beam2 = HB.HermitianBeam2D(E=cls.EE, I=1., A=1., i=cls.nodes[2], j=cls.nodes[3])
        cls.beam_as_structure = HB.Structure(beams=[cls.beam1], supports={1: ['ux', 'uy', 'rotz']})

    def test_element_stiffness_matrix(self):
        k_asinteger = np.matrix([[1, 0, 0, -1, 0, 0],
                                 [0, 12, 6, 0, -12, 6],
                                 [0, 6, 4, 0, -6, 2],
                                 [-1, 0, 0, 1, 0, 0],
                                 [0, -12, -6, 0, 12, -6],
                                 [0, 6, 2, 0, -6, 4]])

        self.assertTrue(np.allclose(k_asinteger, self.beam2.Kg))
        self.assertTrue(np.allclose(self.beam1.Kg, self.beam2.Kg))

    def test_element_stiffness_matrix_symmetry(self):
        # tests if the matrix is symmetric
        diff = self.beam1.Kg.transpose() - self.beam1.Kg
        self.assertTrue(np.allclose(diff, np.zeros(6)))

    def test_element_stiffness_matrix_positive_definite(self):
        # tests if the matrix is positive definite. The element stiffness matrix without proper constraints is not.
        self.assertRaises(np.linalg.linalg.LinAlgError, np.linalg.cholesky, self.beam1.Kg)

    def test_element_stiffness_matrix_positive_definite_constrained(self):
        # tests if the matrix is positive definite. The element stiffness matrix without proper constraints is not.
        # therefore all DOFs iat node i are eliminated
        k = self.beam1.Kg
        for i in range(3):
            k = np.delete(k, 0, axis=0)
            k = np.delete(k, 0, axis=1)

        # the inverse of the constrained element
        inverted = np.matrix([[1., 0., 0.],
                              [0., 3.46410162, 0.],
                              [0., -1.73205081, 1.]])

        self.assertTrue(np.allclose(inverted, np.linalg.cholesky(k)))

    def test_FX_load(self):
        self.beam_as_structure.add_single_dynam_to_node(nodeID=2, dynam={'FX': 1}, clear=True)
        self.beam_as_structure.solve()
        disps = self.beam_as_structure.displacements
        _expected = np.matrix([[0.0],
                               [0.0],
                               [0.0],
                               [1.0],
                               [0.0],
                               [0.0]])
        self.assertTrue(np.allclose(disps, _expected, atol=ATOL))

    def test_FY_load(self):
        # tests a clamped beam for vertical load at the tip
        import random
        P = random.randrange(1, 200)
        L = self.beam1.l
        EI = self.beam1.EI
        self.beam_as_structure.add_single_dynam_to_node(nodeID=2, dynam={'FY': P}, clear=True)
        self.beam_as_structure.solve()
        disps = self.beam_as_structure.displacements
        _expected = np.matrix([[0.0],
                               [0.0],
                               [0.0],
                               [0.0],
                               [(P * L ** 3) / (3 * EI)],  # vertical displacement at the end
                               [(P * L ** 2) / (2 * EI)]])  # rotation at the end
        self.assertTrue(np.allclose(disps, _expected, atol=ATOL))

    def test_MZ_load(self):
        # tests a clamped beam for bending moment at the tip
        import random
        M = random.randrange(1, 200)
        L = self.beam1.l
        EI = self.beam1.EI
        self.beam_as_structure.add_single_dynam_to_node(nodeID=2, dynam={'MZ': M}, clear=True)
        self.beam_as_structure.solve()
        disps = self.beam_as_structure.displacements
        _expected = np.matrix([[0.0],
                               [0.0],
                               [0.0],
                               [0.0],
                               [(M * L ** 2) / (2 * EI)],  # vertical displacement at the end
                               [(M * L) / EI]])  # rotation at the end
        self.assertTrue(np.allclose(disps, _expected, atol=ATOL))


class Hermitian2D_Structure(unittest.TestCase):
    """
    tests for the hermitian 2D beam
    """

    @classmethod
    def setUpClass(cls):
        # cls.EE = 1.
        # cls.nodes = {1: HB.Node(ID=1, coords=(0, 0)), 2: HB.Node(ID=2, coords=(1, 0)), 3: HB.Node(ID=3, coords=(2, 0))}
        # cls.beam1 = HB.HermitianBeam2D(E=cls.EE, I=1., A=1., i=cls.nodes[1], j=cls.nodes[2])
        # cls.beam2 = HB.HermitianBeam2D(E=cls.EE, I=1., A=1., i=cls.nodes[2], j=cls.nodes[3])

        # structure #1: 3 elements in a row along the global X axis
        cls.n1 = HB.Node(ID=1, coords=(0, 0))
        cls.n2 = HB.Node(ID=2, coords=(100, 0))
        cls.n3 = HB.Node(ID=3, coords=(200, 0))
        cls.n4 = HB.Node(ID=4, coords=(300, 0))

        sect = sections.Recangle(a=10, b=10)
        b1 = HB.HermitianBeam2D(E=210000., ID=1, crosssection=sect, i=cls.n1, j=cls.n2)
        b2 = HB.HermitianBeam2D(E=210000., ID=2, crosssection=sect, i=cls.n2, j=cls.n3)
        b3 = HB.HermitianBeam2D(E=210000., ID=3, crosssection=sect, i=cls.n3, j=cls.n4)
        cls.structure_1 = HB.Structure(beams=[b1, b2, b3], supports={1: ['ux', 'uy', 'rotz']})

        # the same structure, with the last node not inline. to test rotation.
        cls.n5 = HB.Node(ID=4, coords=(300, 100))
        b4 = HB.HermitianBeam2D(E=210000., ID=3, I=833.33, A=100., i=cls.n3, j=cls.n5)
        cls.structure_2 = HB.Structure(beams=[b1, b2, b4], supports={1: ['ux', 'uy', 'rotz']})

    def test_nodal_displacements_1(self):
        # assertion #1: single Axial force on Node #4
        self.structure_1.add_single_dynam_to_node(nodeID=4, dynam={'FX': 1}, clear=True)
        self.structure_1.solve()
        disps = self.structure_1.displacements
        _expected = np.matrix([[0.0],
                               [0.0],
                               [0.0],
                               [4.76190476e-06],
                               [0.00000000e+00],
                               [0.00000000e+00],
                               [9.52380952e-06],
                               [0.00000000e+00],
                               [0.00000000e+00],
                               [1.42857143e-05],
                               [0.00000000e+00],
                               [0.00000000e+00]])
        self.assertTrue(np.allclose(disps, _expected, atol=ATOL))

    def test_nodal_displacements_11(self):
        # assertion #1: single Axial force on Node 2, the element is rotated, so are the loads.
        # but then at the end the resulting displacements are rotated back and we get the same results as in test _1

        # structure #1: 3 elements in a row along the global X axis
        _u = math.sqrt(2) / 2.  # unit; instead of 1 as unit we have cos(45) to account for the rotation
        n1 = HB.Node(ID=1, coords=(0, 0))
        n2 = HB.Node(ID=2, coords=(100 * _u, 100 * _u))
        n3 = HB.Node(ID=3, coords=(200 * _u, 200 * _u))
        n4 = HB.Node(ID=4, coords=(300 * _u, 300 * _u))

        # definition of the beams is unchanged
        b1 = HB.HermitianBeam2D(E=210000., ID=1, I=833.33, A=100., i=n1, j=n2)
        b2 = HB.HermitianBeam2D(E=210000., ID=2, I=833.33, A=100., i=n2, j=n3)
        b3 = HB.HermitianBeam2D(E=210000., ID=3, I=833.33, A=100., i=n3, j=n4)
        structure = HB.Structure(beams=[b1, b2, b3], supports={1: ['ux', 'uy', 'rotz']})

        # two loads instead of one as these are in the global system
        structure.add_single_dynam_to_node(nodeID=4, dynam={'FX': _u}, clear=True)
        structure.add_single_dynam_to_node(nodeID=4, dynam={'FY': _u})

        # since the elements were previously rotated by 45 degrees, now we rotate the results by -45 degrees back...
        T = HB.transfer_matrix(-45, asdegree=True, blocks=len(structure.nodes))
        structure.solve()
        disps = T * structure.displacements

        _expected = np.matrix([[0.0],
                               [0.0],
                               [0.0],
                               [4.76190476e-06],  # same as in _1
                               [0.00000000e+00],
                               [0.00000000e+00],
                               [9.52380952e-06],
                               [0.00000000e+00],
                               [0.00000000e+00],
                               [1.42857143e-05],
                               [0.00000000e+00],
                               [0.00000000e+00]])
        self.assertTrue(np.allclose(disps, _expected, atol=ATOL))

    def test_nodal_displacements_12(self):
        # same as _11, but the beams are rotated 60 degree clockwise
        # assertion #1: single Axial force on Node 2, the element is rotated, so are the loads.
        # but then at the end the resulting displacements are rotated back and we get the same results as in test _1

        # structure #1: 3 elements in a row along the global X axis
        _ux = 0.5  # unit in x direction; instead of 1 as unit we have cos(60) to account for the rotation
        _uy = -math.sqrt(3) / 2.  # unit in x direction; instead of 1 as unit we have cos(60) to account for the rotation
        n1 = HB.Node(ID=1, coords=(0, 0))
        n2 = HB.Node(ID=2, coords=(100 * _ux, 100 * _uy))
        n3 = HB.Node(ID=3, coords=(200 * _ux, 200 * _uy))
        n4 = HB.Node(ID=4, coords=(300 * _ux, 300 * _uy))

        # definition of the beams is unchanged
        b1 = HB.HermitianBeam2D(E=210000., ID=1, I=833.33, A=100., i=n1, j=n2)
        b2 = HB.HermitianBeam2D(E=210000., ID=2, I=833.33, A=100., i=n2, j=n3)
        b3 = HB.HermitianBeam2D(E=210000., ID=3, I=833.33, A=100., i=n3, j=n4)
        structure = HB.Structure(beams=[b1, b2, b3], supports={1: ['ux', 'uy', 'rotz']})

        # two loads instead of one as these are in the global system
        structure.add_single_dynam_to_node(nodeID=4, dynam={'FX': _ux}, clear=True)
        structure.add_single_dynam_to_node(nodeID=4, dynam={'FY': _uy})

        # since the elements were previously rotated by -60 degrees, now we rotate the results by 60 degrees back...
        T = HB.transfer_matrix(60, asdegree=True, blocks=len(structure.nodes))
        structure.solve()
        disps = T * structure.displacements
        _expected = np.matrix([[0.0],
                               [0.0],
                               [0.0],
                               [4.76190476e-06],  # same as in _1
                               [0.00000000e+00],
                               [0.00000000e+00],
                               [9.52380952e-06],
                               [0.00000000e+00],
                               [0.00000000e+00],
                               [1.42857143e-05],
                               [0.00000000e+00],
                               [0.00000000e+00]])
        self.assertTrue(np.allclose(disps, _expected, atol=ATOL))

    def test_nodal_displacements_2(self):
        # assertion #2: single shear load at Node 2
        self.structure_1.add_single_dynam_to_node(nodeID=3, dynam={'FY': -1}, clear=True)
        self.structure_1.solve()
        disps = self.structure_1.displacements
        _expected = np.matrix([[0.000000e+00],
                               [-1.000000e-21],
                               [-2.000000e-19],
                               [0.000000e+00],
                               [-4.761905e-03],
                               [-8.571429e-05],
                               [0.000000e+00],
                               [-1.523810e-02],
                               [-1.142857e-04],
                               [0.000000e+00],
                               [-2.666667e-02],
                               [-1.142857e-04]])
        self.assertTrue(np.allclose(disps, _expected, atol=ATOL))

    def test_nodal_displacements_3(self):
        # assertion #3: single concentrated moment at node 3
        self.structure_1.add_single_dynam_to_node(nodeID=3, dynam={'MZ': 1200000}, clear=True)
        disps = self.structure_1.solve()
        _expected = np.matrix([[0.0],
                               [0.0],
                               [0.0],
                               [4.76190476e-06],
                               [2.85667809e+01],
                               [5.71345143e-01],
                               [9.52380952e-06],
                               [1.14270933e+02],
                               [1.14274743e+00],
                               [1.42857143e-05],
                               [2.28545676e+02],
                               [1.14274743e+00]])
        # print(self.structure_1.draw())
        # print('_3')
        # print(disps)
        # print('')
        # print(disps-_expected)

        self.assertTrue(np.allclose(disps, _expected, atol=ATOL))

    def test_nodal_displacements_4(self):
        # assertion #4: lots of loads at Node 3
        self.structure_1.add_single_dynam_to_node(nodeID=3, dynam={'FX': 1000, 'FY': 1000, 'MZ': 1000000}, clear=True)
        self.structure_1.solve()
        disps = self.structure_1.displacements
        _expected = np.matrix([[1.000000e-48],
                               [1.000000e-48],
                               [1.200000e-45],
                               [4.761905e-03],
                               [3.333333e+01],
                               [6.571429e-01],
                               [9.523810e-03],
                               [1.295238e+02],
                               [1.257143e+00],
                               [9.523810e-03],
                               [2.552381e+02],
                               [1.257143e+00]])
        self.assertTrue(np.allclose(disps, _expected, atol=ATOL))

    def test_nodal_displacements_5(self):
        # assertion #4: lots of loads at Node 4
        self.structure_1.add_single_dynam_to_node(nodeID=4, dynam={'FX': 1000, 'FY': 1000, 'MZ': 1000000}, clear=True)
        self.structure_1.solve()
        disps = self.structure_1.displacements
        _expected = np.matrix([[1.000000e-48],
                               [1.000000e-48],
                               [1.300000e-45],
                               [4.761905e-03],
                               [3.619048e+01],
                               [7.142857e-01],
                               [9.523810e-03],
                               [1.409524e+02],
                               [1.371429e+00],
                               [1.428571e-02],
                               [3.085714e+02],
                               [1.971429e+00]])

        self.assertTrue(np.allclose(disps, _expected, atol=ATOL))

    # todo:
    # test added supports on differet structures


    # non-inline structures
    def test_nodal_displacements_6(self):
        # assertion #5: shear load on node 3
        self.structure_2.add_single_dynam_to_node(nodeID=3, dynam={'FY': -1000}, clear=True)
        self.structure_2.solve()
        disps = self.structure_2.displacements
        _expected = np.matrix([[0.],
                               [0.],
                               [0.],
                               [0.],
                               [-4.76192381e+00],
                               [-8.57146286e-02],
                               [0.],
                               [-1.52381562e+01],
                               [-1.14286171e-01],
                               [1.14286171e+01],
                               [-2.66667733e+01],
                               [-1.14286171e-01]])
        self.assertTrue(np.allclose(disps, _expected, atol=ATOL))

    def test_nodal_displacements_7(self):
        # assertion #5: shear load on node 2 and node 3
        self.structure_2.add_single_dynam_to_node(nodeID=3, dynam={'FX': 1000, 'FY': 1000, 'MZ': 1000000}, clear=True)
        self.structure_2.solve()
        disps = self.structure_2.displacements
        _expected = np.matrix([[1.000000e-48],
                               [1.000000e-48],
                               [1.200000e-45],
                               [4.761905e-03],
                               [3.333333e+01],
                               [6.571429e-01],
                               [9.523810e-03],
                               [1.295238e+02],
                               [1.257143e+00],
                               [-1.257048e+02],
                               [2.552381e+02],
                               [1.257143e+00]])
        self.structure_2.draw()
        self.assertTrue(np.allclose(disps, _expected, atol=ATOL))

    def test_nodal_displacements_8(self):
        # assertion #5: shear load on node 2 and node 3
        self.structure_2.add_single_dynam_to_node(nodeID=4, dynam={'FX': 0, 'FY': -1, 'MZ': 0}, clear=True)
        self.structure_2.solve()
        disps = self.structure_2.displacements
        # todo: this one fails, probably the transformation is not OK
        _expected = np.matrix([[0.000000e+00],
                               [-1.000000e-51],
                               [-3.000000e-49],
                               [0.000000e+00],
                               [-7.619048e-03],
                               [-1.428571e-04],
                               [0.000000e+00],
                               [-2.666667e-02],
                               [-2.285714e-04],
                               [2.554753e-02],
                               [-5.222093e-02],
                               [-2.689777e-04]])
        self.assertTrue(np.allclose(disps, _expected, atol=ATOL))
