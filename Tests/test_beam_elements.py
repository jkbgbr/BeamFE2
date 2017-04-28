import unittest
from Beams import HermitianBeam as HB
import numpy as np
import math


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

        # structure #1: 3 elements in a row along the global X axis
        cls.n1 = HB.Node(ID=1, coords=(0, 0))
        cls.n2 = HB.Node(ID=2, coords=(100, 0))
        cls.n3 = HB.Node(ID=3, coords=(200, 0))
        cls.n4 = HB.Node(ID=4, coords=(300, 0))
        b1 = HB.HermitianBeam2D(E=210000., ID=1, I=833.33, A=100., i=cls.n1, j=cls.n2)
        b2 = HB.HermitianBeam2D(E=210000., ID=2, I=833.33, A=100., i=cls.n2, j=cls.n3)
        b3 = HB.HermitianBeam2D(E=210000., ID=3, I=833.33, A=100., i=cls.n3, j=cls.n4)
        cls.structure_1 = HB.Structure(beams=[b1, b2, b3], BCs=None)

        # the same structure, with the last node not inline. to test rotation.
        cls.n5 = HB.Node(ID=4, coords=(300, 100))
        b4 = HB.HermitianBeam2D(E=210000., ID=3, I=833.33, A=100., i=cls.n3, j=cls.n5)
        cls.structure_2 = HB.Structure(beams=[b1, b2, b4], BCs=None)

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

    def test_nodal_displacements_1(self):
        # assertion #1: single Axial force on Node 2
        self.structure_1.add_single_dynam_to_node(nodeID=3, dynam={'FX': 1}, clear=True)
        disps = self.structure_1.solve()
        _expected = np.matrix([[4.76190476e-06],
                               [0.00000000e+00],
                               [0.00000000e+00],
                               [9.52380952e-06],
                               [0.00000000e+00],
                               [0.00000000e+00],
                               [1.42857143e-05],
                               [0.00000000e+00],
                               [0.00000000e+00]])
        self.assertTrue(np.allclose(disps, _expected))

    def test_nodal_displacements_11(self):
        # assertion #1: single Axial force on Node 2, the element is rotated

        # structure #1: 3 elements in a row along the global X axis
        _u = math.sqrt(2) / 2.  # unit
        n1 = HB.Node(ID=1, coords=(0, 0))
        n2 = HB.Node(ID=2, coords=(100 * _u, 100 * _u))
        n3 = HB.Node(ID=3, coords=(200 * _u, 200 * _u))
        n4 = HB.Node(ID=4, coords=(300 * _u, 300 * _u))
        b1 = HB.HermitianBeam2D(E=210000., ID=1, I=833.33, A=100., i=n1, j=n2)
        b2 = HB.HermitianBeam2D(E=210000., ID=2, I=833.33, A=100., i=n2, j=n3)
        b3 = HB.HermitianBeam2D(E=210000., ID=3, I=833.33, A=100., i=n3, j=n4)
        structure = HB.Structure(beams=[b1, b2, b3], BCs=None)
        structure.add_single_dynam_to_node(nodeID=3, dynam={'FX': _u}, clear=True)
        structure.add_single_dynam_to_node(nodeID=3, dynam={'FY': _u})
        disps = structure.solve()
        print(disps)
        _expected = np.matrix([[4.76190476e-06],
                               [0.00000000e+00],
                               [0.00000000e+00],
                               [9.52380952e-06],
                               [0.00000000e+00],
                               [0.00000000e+00],
                               [1.42857143e-05],
                               [0.00000000e+00],
                               [0.00000000e+00]])
        self.assertTrue(np.allclose(disps, _expected))


    def test_nodal_displacements_2(self):
        # assertion #2: single shear load at Node 2
        self.structure_1.add_single_dynam_to_node(nodeID=2, dynam={'FY': -1}, clear=True)
        disps = self.structure_1.solve()
        _expected = np.matrix([[4.76190476e-06],
                               [-4.76192381e-03],
                               [-8.57146286e-05],
                               [9.52380952e-06],
                               [-1.52381562e-02],
                               [-1.14286171e-04],
                               [1.42857143e-05],
                               [-2.66667733e-02],
                               [-1.14286171e-04]])
        self.assertTrue(np.allclose(disps, _expected))

    def test_nodal_displacements_3(self):
        # assertion #3: single concentrated moment at node 2
        self.structure_1.add_single_dynam_to_node(nodeID=2, dynam={'MZ': 1000000}, clear=True)
        disps = self.structure_1.solve()
        _expected = np.matrix([[4.76190476e-06],
                               [2.85667809e+01],
                               [5.71345143e-01],
                               [9.52380952e-06],
                               [1.14270933e+02],
                               [1.14274743e+00],
                               [1.42857143e-05],
                               [2.28545676e+02],
                               [1.14274743e+00]])
        self.assertTrue(np.allclose(disps, _expected))

    def test_nodal_displacements_4(self):
        # assertion #4: lots of loads at Node 2
        self.structure_1.add_single_dynam_to_node(nodeID=2, dynam={'FX': 1000, 'FY': 1000, 'MZ': 1000000}, clear=True)
        disps = self.structure_1.solve()
        _expected = np.matrix([[4.76666667e-03],
                               [3.33334667e+01],
                               [6.57145486e-01],
                               [9.53333333e-03],
                               [1.29524328e+02],
                               [1.25714789e+00],
                               [9.53809524e-03],
                               [2.55239116e+02],
                               [1.25714789e+00]])
        self.assertTrue(np.allclose(disps, _expected))

    # non-inline structures
    def test_nodal_displacements_5(self):
        # assertion #5: shear load on node 2
        self.structure_2.add_single_dynam_to_node(nodeID=2, dynam={'FX': 0, 'FY': 1000, 'MZ': 0}, clear=True)
        disps = self.structure_2.solve()
        _expected = np.matrix([[6.32938743e-16],
                               [4.76192381e+00],
                               [8.57146286e-02],
                               [1.26587749e-15],
                               [1.52381562e+01],
                               [1.14286171e-01],
                               [1.14286171e+01],
                               [2.66667733e+01],
                               [1.14286171e-01]])
        self.assertTrue(np.allclose(disps, _expected))

    def test_nodal_displacements_6(self):
        # assertion #5: shear load on node 2 and node 3
        self.structure_2.add_single_dynam_to_node(nodeID=2, dynam={'FX': 0, 'FY': 1000, 'MZ': 0}, clear=True)
        self.structure_2.add_single_dynam_to_node(nodeID=3, dynam={'FX': 0, 'FY': 1000, 'MZ': 0})
        disps = self.structure_2.solve()
        _expected = np.matrix([[1.26587749e-15],
                               [1.23810019e+01],
                               [2.28572343e-01],
                               [2.53175497e-15],
                               [4.19049295e+01],
                               [3.42858514e-01],
                               [3.69762351e+01],
                               [7.88878990e+01],
                               [3.83264778e-01]])
        self.assertTrue(np.allclose(disps, _expected))

    def test_nodal_displacements_7(self):
        # assertion #5: shear load on node 2 and node 3
        self.structure_2.add_single_dynam_to_node(nodeID=2, dynam={'FX': 0, 'FY': 1000, 'MZ': 0}, clear=True)
        self.structure_2.add_single_dynam_to_node(nodeID=3, dynam={'FX': 1000, 'FY': 1000, 'MZ': 0})
        disps = self.structure_2.solve()
        print(disps)
        # todo: this one fails, probably the transformation is not OK
        _expected = np.matrix([[1.26587749e-15],
                               [1.23810019e+01],
                               [2.28572343e-01],
                               [2.53175497e-15],
                               [4.19049295e+01],
                               [3.42858514e-01],
                               [3.69762351e+01],
                               [7.88878990e+01],
                               [3.83264778e-01]])
        self.assertTrue(np.allclose(disps, _expected))
