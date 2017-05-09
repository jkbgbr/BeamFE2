import unittest
from Modell import HermitianBeam_2D as HB
from Modell import BeamSections as sections
from Modell import Structure
from Modell import Node
from solver import solve

import numpy as np
import math
from Tests import ATOL


class Hermitian2D_Element(unittest.TestCase):
    """
    basic tests for the hermitian 2D beam
    """

    @classmethod
    def setUpClass(cls):
        cls.EE = 1.
        cls.nodeset_1 = {1: Node.Node(ID=1, coords=(0, 0)), 2: Node.Node(ID=2, coords=(1, 0)), 3: Node.Node(ID=3, coords=(2, 0))}
        cls.beam1 = HB.HermitianBeam2D(E=cls.EE, I=1., A=1., i=cls.nodeset_1[1], j=cls.nodeset_1[2])
        cls.beam2 = HB.HermitianBeam2D(E=cls.EE, I=1., A=1., i=cls.nodeset_1[2], j=cls.nodeset_1[3])
        cls.beam_as_structure = Structure.Structure(beams=[cls.beam1], supports={1: ['ux', 'uy', 'rotz']})

        # the same beam but inclined by atan(0,5) ~= 26.56505 degrees degrees CCW.
        # all tests are performed for these as well, results should reflect only the rotation
        cls.rotation_angle = math.atan(0.5)
        cls._r = 1 / math.sqrt(1.25)  # the lenght of the elements is modified to keep the beam length
        cls.nodeset_2 = {1: Node.Node(ID=1, coords=(0, 0)),
                         2: Node.Node(ID=2, coords=(1. * cls._r, 0.5 * cls._r)),
                         3: Node.Node(ID=3, coords=(2. * cls._r, 1. * cls._r))}
        cls.beam3 = HB.HermitianBeam2D(E=cls.EE, I=1., A=1., i=cls.nodeset_2[1], j=cls.nodeset_2[2])
        cls.beam4 = HB.HermitianBeam2D(E=cls.EE, I=1., A=1., i=cls.nodeset_2[2], j=cls.nodeset_2[3])
        cls.beam_as_structure_2 = Structure.Structure(beams=[cls.beam3], supports={1: ['ux', 'uy', 'rotz']})

    def test_element_stiffness_matrix(self):
        k_asinteger = np.matrix([[1, 0, 0, -1, 0, 0],
                                 [0, 12, 6, 0, -12, 6],
                                 [0, 6, 4, 0, -6, 2],
                                 [-1, 0, 0, 1, 0, 0],
                                 [0, -12, -6, 0, 12, -6],
                                 [0, 6, 2, 0, -6, 4]])

        self.assertTrue(np.allclose(k_asinteger, self.beam2.Ke))
        self.assertTrue(np.allclose(self.beam1.Ke, self.beam2.Ke))

    def test_element_stiffness_matrix_symmetry(self):
        # tests if the matrix is symmetric
        diff = self.beam1.Ke.transpose() - self.beam1.Ke
        self.assertTrue(np.allclose(diff, np.zeros(6)))

    def test_element_stiffness_matrix_positive_definite(self):
        # tests if the matrix is positive definite. The element stiffness matrix without proper constraints is not.
        self.assertRaises(np.linalg.linalg.LinAlgError, np.linalg.cholesky, self.beam1.Ke)

    def test_element_stiffness_matrix_positive_definite_constrained(self):
        # tests if the matrix is positive definite. The element stiffness matrix without proper constraints is not.
        # therefore all DOFs iat node i are eliminated
        k = self.beam1.Ke
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
        solve(self.beam_as_structure, analysis='linear static')
        disps = self.beam_as_structure.results['linear static'].element_displacements(local=True, beam=self.beam_as_structure.beams[0])
        _expected = {'ux': np.matrix([[0.], [1.]]),
                     'uy': np.matrix([[0.], [0.]]),
                     'rotz': np.matrix([[0.], [0.]])}
        self.assertTrue(np.allclose(disps['ux'], _expected['ux'], atol=ATOL))
        self.assertTrue(np.allclose(disps['uy'], _expected['uy'], atol=ATOL))
        self.assertTrue(np.allclose(disps['rotz'], _expected['rotz'], atol=ATOL))

    def test_FY_load(self):
        # tests a clamped beam for vertical load at the tip
        import random
        P = random.randrange(1, 200)
        L = self.beam1.l
        EI = self.beam1.EI
        self.beam_as_structure.add_single_dynam_to_node(nodeID=2, dynam={'FY': P}, clear=True)
        solve(self.beam_as_structure, analysis='linear static')
        disps = self.beam_as_structure.results['linear static'].element_displacements(local=True, beam=self.beam_as_structure.beams[0])
        _expected = {'ux': np.matrix([[0.], [0.]]),
                     'uy': np.matrix([[0.], [(P * L ** 3) / (3 * EI)]]),
                     'rotz': np.matrix([[0.], [(P * L ** 2) / (2 * EI)]])}

        self.assertTrue(np.allclose(disps['ux'], _expected['ux'], atol=ATOL))
        self.assertTrue(np.allclose(disps['uy'], _expected['uy'], atol=ATOL))
        self.assertTrue(np.allclose(disps['rotz'], _expected['rotz'], atol=ATOL))

    def test_MZ_load(self):
        # tests a clamped beam for bending moment at the tip
        import random
        M = random.randrange(1, 200)
        L = self.beam1.l
        EI = self.beam1.EI
        self.beam_as_structure.add_single_dynam_to_node(nodeID=2, dynam={'MZ': M}, clear=True)
        solve(self.beam_as_structure, analysis='linear static')
        disps = self.beam_as_structure.results['linear static'].element_displacements(local=True, beam=self.beam_as_structure.beams[0])
        _expected = {'ux': np.matrix([[0.], [0.]]),
                     'uy': np.matrix([[0.], [(M * L ** 2) / (2 * EI)]]),
                     'rotz': np.matrix([[0.], [(M * L) / EI]])}

        self.assertTrue(np.allclose(disps['ux'], _expected['ux'], atol=ATOL))
        self.assertTrue(np.allclose(disps['uy'], _expected['uy'], atol=ATOL))
        self.assertTrue(np.allclose(disps['rotz'], _expected['rotz'], atol=ATOL))

    def test_nodal_reactions(self):
        # FX
        self.beam_as_structure.add_single_dynam_to_node(nodeID=2, dynam={'FX': 1}, clear=True)
        solve(self.beam_as_structure, analysis='linear static')
        for beam in self.beam_as_structure.beams:
            disps = self.beam_as_structure.results['linear static'].element_displacements(local=True, beam=beam, asvector=True)
        #     print(beam.Ke * disps)
        # print('')
        # FY
        self.beam_as_structure.add_single_dynam_to_node(nodeID=2, dynam={'FY': 1}, clear=True)
        solve(self.beam_as_structure, analysis='linear static')
        for beam in self.beam_as_structure.beams:
            disps = self.beam_as_structure.results['linear static'].element_displacements(local=True, beam=beam, asvector=True)
        #     print(beam.Ke * disps)
        # print('')
        # MZ
        self.beam_as_structure.add_single_dynam_to_node(nodeID=2, dynam={'MZ': 1}, clear=True)
        solve(self.beam_as_structure, analysis='linear static')
        for beam in self.beam_as_structure.beams:
            disps = self.beam_as_structure.results['linear static'].element_displacements(local=True, beam=beam, asvector=True)
            solve(self.beam_as_structure, analysis='linear static')

    # tests on the rotated structure
    def test_element_stiffness_matrix_rotated(self):
        k_asinteger = np.matrix([[1, 0, 0, -1, 0, 0],
                                 [0, 12, 6, 0, -12, 6],
                                 [0, 6, 4, 0, -6, 2],
                                 [-1, 0, 0, 1, 0, 0],
                                 [0, -12, -6, 0, 12, -6],
                                 [0, 6, 2, 0, -6, 4]])

        _tm = HB.transfer_matrix(alpha=self.rotation_angle, asdegree=False, blocks=2, blocksize=3)
        k_asinteger = _tm * k_asinteger * _tm.T
        self.assertTrue(np.allclose(k_asinteger, self.beam3.matrix_in_global(mtrx=self.beam3.Ke)))
        self.assertTrue(np.allclose(self.beam3.matrix_in_global(mtrx=self.beam3.Ke), self.beam4.matrix_in_global(mtrx=self.beam4.Ke)))

    def test_element_stiffness_matrix_symmetry_rotated(self):
        # tests if the matrix is symmetric
        diff = self.beam3.Ke.transpose() - self.beam3.Ke
        self.assertTrue(np.allclose(diff, np.zeros(6)))

    def test_element_stiffness_matrix_positive_definite_rotated(self):
        # tests if the matrix is positive definite. The element stiffness matrix without proper constraints is not.
        self.assertRaises(np.linalg.linalg.LinAlgError, np.linalg.cholesky, self.beam3.Ke)

    def test_FX_load_rotated(self):
        self.beam_as_structure_2.add_single_dynam_to_node(nodeID=2, dynam={'FX': 1., 'FY': 0.5}, clear=True)
        solve(self.beam_as_structure_2, analysis='linear static')

        # checking the displacements in the global system - results of the structure
        disps = self.beam_as_structure_2.results['linear static'].element_displacements(local=False, beam=self.beam_as_structure_2.beams[0])

        _expected = {'ux': np.matrix([[0.], [1.]]),
                     'uy': np.matrix([[0.], [0.5]]),
                     'rotz': np.matrix([[0.], [0.]])}
        self.assertTrue(np.allclose(disps['ux'], _expected['ux'], atol=ATOL))
        self.assertTrue(np.allclose(disps['uy'], _expected['uy'], atol=ATOL))
        self.assertTrue(np.allclose(disps['rotz'], _expected['rotz'], atol=ATOL))

        # # checking the displacements in the local system - results of the first an only beam
        disps = self.beam_as_structure_2.results['linear static'].element_displacements(local=True, beam=self.beam_as_structure_2.beams[0])
        _expected = {'ux': np.matrix([[0.], [1. / self._r]]),
                     'uy': np.matrix([[0.], [0.]]),
                     'rotz': np.matrix([[0.], [0.]])}

        self.assertTrue(np.allclose(disps['ux'], _expected['ux'], atol=ATOL))
        self.assertTrue(np.allclose(disps['uy'], _expected['uy'], atol=ATOL))
        self.assertTrue(np.allclose(disps['rotz'], _expected['rotz'], atol=ATOL))

    #
    def test_FY_load_rotated(self):
        # tests a clamped beam for vertical load at the tip
        import random
        P = random.randrange(1, 200)
        L = self.beam3.l
        EI = self.beam3.EI
        _t = HB.transfer_matrix(alpha=self.rotation_angle, asdegree=False, blocks=1, blocksize=2)
        _load = _t * np.matrix([[0, P]]).T

        self.beam_as_structure_2.add_single_dynam_to_node(nodeID=2, dynam={'FX': _load[0], 'FY': _load[1]}, clear=True)
        solve(self.beam_as_structure_2, analysis='linear static')
        disps = self.beam_as_structure_2.results['linear static'].element_displacements(local=True, beam=self.beam_as_structure_2.beams[0], asvector=True)
        _expected = np.matrix([[0.0],
                               [0.0],
                               [0.0],
                               [0.0],
                               [(P * L ** 3) / (3 * EI)],  # vertical displacement at the end
                               [(P * L ** 2) / (2 * EI)]])  # rotation at the end
        self.assertTrue(np.allclose(disps, _expected, atol=ATOL))

    def test_MZ_load_rotated(self):
        # tests a clamped beam for bending moment at the tip
        import random
        M = random.randrange(1, 200)
        L = self.beam3.l
        EI = self.beam3.EI
        self.beam_as_structure_2.add_single_dynam_to_node(nodeID=2, dynam={'MZ': M}, clear=True)
        solve(self.beam_as_structure_2, analysis='linear static')
        disps = self.beam_as_structure_2.results['linear static'].element_displacements(local=True, beam=self.beam_as_structure_2.beams[0], asvector=True)
        _expected = np.matrix([[0.0],
                               [0.0],
                               [0.0],
                               [0.0],
                               [(M * L ** 2) / (2 * EI)],  # vertical displacement at the end
                               [(M * L) / EI]])  # rotation at the end
        self.assertTrue(np.allclose(disps, _expected, atol=ATOL))


class Hermitian2D_Structure(unittest.TestCase):
    """
    advanced tests for the hermitian 2D beam
    """

    @classmethod
    def setUpClass(cls):

        # structure #1: 3 elements in a row along the global X axis
        cls.n1 = Node.Node(ID=1, coords=(0, 0))
        cls.n2 = Node.Node(ID=2, coords=(100, 0))
        cls.n3 = Node.Node(ID=3, coords=(200, 0))
        cls.n4 = Node.Node(ID=4, coords=(300, 0))

        # section for the beams: 10 by 10 rectangle
        sect = sections.Rectangle(width=10, height=10)
        b1 = HB.HermitianBeam2D(E=210000., ID=1, crosssection=sect, i=cls.n1, j=cls.n2)
        b2 = HB.HermitianBeam2D(E=210000., ID=2, crosssection=sect, i=cls.n2, j=cls.n3)
        b3 = HB.HermitianBeam2D(E=210000., ID=3, crosssection=sect, i=cls.n3, j=cls.n4)
        cls.structure_1 = Structure.Structure(beams=[b1, b2, b3], supports={1: ['ux', 'uy', 'rotz']})

        # the same structure, with the last node not inline. to test rotation.
        cls.n5 = Node.Node(ID=4, coords=(300, 100))
        b4 = HB.HermitianBeam2D(E=210000., ID=3, I=833.33, A=100., i=cls.n3, j=cls.n5)
        cls.structure_2 = Structure.Structure(beams=[b1, b2, b4], supports={1: ['ux', 'uy', 'rotz']})

    def test_nodal_displacements_1(self):
        # assertion #1: single Axial force on Node #4
        self.structure_1.add_single_dynam_to_node(nodeID=4, dynam={'FX': 1}, clear=True)
        solve(self.structure_1, analysis='linear static')
        disps = self.structure_1.results['linear static'].global_displacements()
        _expected = {'ux': np.matrix([[0.], [4.76190476e-06], [9.52380952e-06], [1.42857143e-05]]),
                     'uy': np.matrix([[0.], [0.], [0.], [0.]]),
                     'rotz': np.matrix([[0.], [0.], [0.], [0.]])}
        self.assertTrue(np.allclose(disps['ux'], _expected['ux'], atol=ATOL))
        self.assertTrue(np.allclose(disps['uy'], _expected['uy'], atol=ATOL))
        self.assertTrue(np.allclose(disps['rotz'], _expected['rotz'], atol=ATOL))

    def test_nodal_displacements_11(self):
        # assertion #1: single Axial force on Node 2, the element is rotated, so are the loads.
        # but then at the end the resulting displacements are rotated back and we get the same results as in test _1

        # structure #1: 3 elements in a row along the global X axis
        _u = math.sqrt(2) / 2.  # unit; instead of 1 as unit we have cos(45) to account for the rotation
        n1 = Node.Node(ID=1, coords=(0, 0))
        n2 = Node.Node(ID=2, coords=(100 * _u, 100 * _u))
        n3 = Node.Node(ID=3, coords=(200 * _u, 200 * _u))
        n4 = Node.Node(ID=4, coords=(300 * _u, 300 * _u))

        # definition of the beams is unchanged
        b1 = HB.HermitianBeam2D(E=210000., ID=1, I=833.33, A=100., i=n1, j=n2)
        b2 = HB.HermitianBeam2D(E=210000., ID=2, I=833.33, A=100., i=n2, j=n3)
        b3 = HB.HermitianBeam2D(E=210000., ID=3, I=833.33, A=100., i=n3, j=n4)
        structure = Structure.Structure(beams=[b1, b2, b3], supports={1: ['ux', 'uy', 'rotz']})

        # two loads instead of one as these are in the global system
        structure.add_single_dynam_to_node(nodeID=4, dynam={'FX': _u}, clear=True)
        structure.add_single_dynam_to_node(nodeID=4, dynam={'FY': _u})

        # since the elements were previously rotated by 45 degrees, now we rotate the results by -45 degrees back...
        T = HB.transfer_matrix(-45, asdegree=True, blocks=len(structure.nodes))
        solve(structure, analysis='linear static')
        disps = structure.results['linear static'].global_displacements(asvector=True)
        disps = T * disps  # in linear static only one solution
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
        n1 = Node.Node(ID=1, coords=(0, 0))
        n2 = Node.Node(ID=2, coords=(100 * _ux, 100 * _uy))
        n3 = Node.Node(ID=3, coords=(200 * _ux, 200 * _uy))
        n4 = Node.Node(ID=4, coords=(300 * _ux, 300 * _uy))

        # definition of the beams is unchanged
        b1 = HB.HermitianBeam2D(E=210000., ID=1, I=833.33, A=100., i=n1, j=n2)
        b2 = HB.HermitianBeam2D(E=210000., ID=2, I=833.33, A=100., i=n2, j=n3)
        b3 = HB.HermitianBeam2D(E=210000., ID=3, I=833.33, A=100., i=n3, j=n4)
        structure = Structure.Structure(beams=[b1, b2, b3], supports={1: ['ux', 'uy', 'rotz']})

        # two loads instead of one as these are in the global system
        structure.add_single_dynam_to_node(nodeID=4, dynam={'FX': _ux}, clear=True)
        structure.add_single_dynam_to_node(nodeID=4, dynam={'FY': _uy})

        # since the elements were previously rotated by -60 degrees, now we rotate the results by 60 degrees back...
        T = HB.transfer_matrix(60, asdegree=True, blocks=len(structure.nodes))
        solve(structure, analysis='linear static')
        disps = structure.results['linear static'].global_displacements(asvector=True)
        disps = T * disps
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
        solve(self.structure_1, analysis='linear static')
        disps = self.structure_1.results['linear static'].global_displacements(asvector=True)
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
        self.structure_1.add_single_dynam_to_node(nodeID=3, dynam={'MZ': 1000000}, clear=True)
        solve(self.structure_1, analysis='linear static')
        disps = self.structure_1.results['linear static'].global_displacements(asvector=True)
        _expected = np.matrix([[0.000000e+00],
                               [-1.040834e-61],
                               [1.000000e-45],
                               [0.000000e+00],
                               [2.857143e+01],
                               [5.714286e-01],
                               [0.000000e+00],
                               [1.142857e+02],
                               [1.142857e+00],
                               [0.000000e+00],
                               [2.285714e+02],
                               [1.142857e+00]])
        self.assertTrue(np.allclose(disps, _expected, atol=ATOL))

    def test_nodal_displacements_4(self):
        # assertion #4: lots of loads at Node 3
        self.structure_1.add_single_dynam_to_node(nodeID=3, dynam={'FX': 1000, 'FY': 1000, 'MZ': 1000000}, clear=True)
        solve(self.structure_1, analysis='linear static')
        disps = self.structure_1.results['linear static'].global_displacements(asvector=True)
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
        solve(self.structure_1, analysis='linear static')
        disps = self.structure_1.results['linear static'].global_displacements(asvector=True)
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

    # non-inline structures
    def test_nodal_displacements_6(self):
        # assertion #5: shear load on node 3
        self.structure_2.add_single_dynam_to_node(nodeID=3, dynam={'FY': -1000}, clear=True)
        solve(self.structure_2, analysis='linear static')
        disps = self.structure_2.results['linear static'].global_displacements(asvector=True)
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
        solve(self.structure_2, analysis='linear static')
        disps = self.structure_2.results['linear static'].global_displacements(asvector=True)
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
        self.assertTrue(np.allclose(disps, _expected, atol=ATOL))

    def test_nodal_displacements_8(self):
        # assertion #5: shear load on node 2 and node 3
        self.structure_2.add_single_dynam_to_node(nodeID=4, dynam={'FX': 0, 'FY': -1, 'MZ': 0}, clear=True)
        solve(self.structure_2, analysis='linear static')
        disps = self.structure_2.results['linear static'].global_displacements(asvector=True)
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

#
# class Hermitian2D_Internal(unittest.TestCase):
#     """
#     advanced tests for the Beams: internal displacements
#     """
#
#     @classmethod
#     def setUpClass(cls):
#
#         # structure #1: 3 elements in a row along the global X axis
#         cls.n1 = Node.Node(ID=1, coords=(0, 0))
#         cls.n2 = Node.Node(ID=2, coords=(100, 0))
#         cls.n3 = Node.Node(ID=3, coords=(200, 0))
#
#         # section for the beams: 10 by 10 rectangle
#         sect = sections.Rectangle(height=30, width=30)
#         b1 = HB.HermitianBeam2D(E=210000., ID=1, crosssection=sect, i=cls.n1, j=cls.n2)
#         b2 = HB.HermitianBeam2D(E=210000., ID=2, crosssection=sect, i=cls.n2, j=cls.n3)
#         cls.structure_1 = Structure.Structure(beams=[b1, b2], supports={1: ['ux', 'uy', 'rotz'], 3: ['ux', 'uy', 'rotz']})
#
#     def test_internal_deflections(self):
#         # for beam in self.structure_1.beams:
#         #     beam.restore_DOFs()
#         self.structure_1.add_single_dynam_to_node(nodeID=2, dynam={'FY': 1000000}, clear=True)
#         solve(self.structure_1, analysis='linear static')
#         disps = self.structure_1.results['linear static'].global_displacements(asvector=True)
#         # self.structure_1.draw(analysistype='linear static')
#         # self.assertTrue(True)