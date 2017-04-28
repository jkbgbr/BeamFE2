import unittest
from Beams import HermitianBeam as HB
import numpy as np


class Hermitian2D(unittest.TestCase):
    """
    tests for the hermitian 2D beam
    """

    @classmethod
    def setUpClass(cls):
        cls.EE = 1.
        cls.topol = {1: HB.Node(ID=1, coords=(0, 0)), 2: HB.Node(ID=2, coords=(1, 0)), 3: HB.Node(ID=3, coords=(2, 0))}
        cls.beam1 = HB.HermitianBeam2D(E=cls.EE, I=1., A=1., i=cls.topol[1], j=cls.topol[2])
        cls.beam2 = HB.HermitianBeam2D(E=cls.EE, I=1., A=1., i=cls.topol[2], j=cls.topol[3])
        cls.Nx = np.matrix([-1, 0, 0, 1, 0, 0]).T  # unit axial Force
        cls.Ny = np.matrix([0, -1, 0, 0, -1, 0]).T  # unit shear Force
        cls.Mz = np.matrix([0, 0, -1, 0, 0, 1]).T  # unit bending moment

        n1 = HB.Node(ID=1, coords=(0, 0))
        n2 = HB.Node(ID=2, coords=(100, 0))
        n3 = HB.Node(ID=3, coords=(200, 0))
        n4 = HB.Node(ID=4, coords=(300, 0))
        b1 = HB.HermitianBeam2D(ID=1, I=833.33, A=100., i=n1, j=n2)
        b2 = HB.HermitianBeam2D(ID=2, I=833.33, A=100., i=n2, j=n3)
        b3 = HB.HermitianBeam2D(ID=3, I=833.33, A=100., i=n3, j=n4)
        cls.structure = HB.Structure(beams=[b1, b2, b3], BCs=None)

    def test_element_stiffness_matrix(self):
        k_asinteger = np.matrix([[1, 0, 0, -1, 0, 0],
                                 [0, 12, 6, 0, -12, 6],
                                 [0, 6, 4, 0, -6, 2],
                                 [-1, 0, 0, 1, 0, 0],
                                 [0, -12, -6, 0, 12, -6],
                                 [0, 6, 2, 0, -6, 4]])
        self.assertTrue(np.allclose(k_asinteger, self.beam2.K))
        self.assertTrue(np.allclose(self.beam1.K, self.beam2.K))

    def test_element_stiffness_matrix_symmetry(self):
        # tests if the matrix is symmetric
        diff = self.beam1.K.transpose() - self.beam1.K
        self.assertTrue(np.allclose(diff, np.zeros(6)))

    def test_element_stiffness_matrix_positive_definite(self):
        # tests if the matrix is positive definite. The element stiffness matrix without proper constraints is not.
        self.assertRaises(np.linalg.linalg.LinAlgError, np.linalg.cholesky, self.beam1.K)

    def test_element_stiffness_matrix_positive_definite_constrained(self):
        # tests if the matrix is positive definite. The element stiffness matrix without proper constraints is not.
        # therefore all DOFs iat node i are eliminated
        k = self.beam1.K
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
        self.structure.add_single_dynam_to_node(nodeID=3, dynam={'FX': 1}, clear=True)
        disps = self.structure.solve()
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
        self.structure.add_single_dynam_to_node(nodeID=2, dynam={'FY': -1}, clear=True)
        disps = self.structure.solve()
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
        self.structure.add_single_dynam_to_node(nodeID=2, dynam={'MZ': 1000000}, clear=True)
        disps = self.structure.solve()
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
        self.structure.add_single_dynam_to_node(nodeID=2, dynam={'FX': 1000, 'FY': 1000, 'MZ': 1000000}, clear=True)
        disps = self.structure.solve()
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
