import unittest
from BeamFE2 import HermitianBeam_2D as HB
from BeamFE2 import BeamSections as sections
from BeamFE2 import Structure
from BeamFE2 import Node
from BeamFE2 import Material
import math


class Hermitian2D_Model(unittest.TestCase):
    """
    tests for the modal analysis of the FE BeamFE2
    """

    @classmethod
    def setUpClass(cls):

        """
        A beam with clamped end under self weight.
        Results are accepted if circular frequencies are within 1% of the analytical solution.
        http://iitg.vlab.co.in/?sub=62&brch=175&sim=1080&cnt=1
        """

        # steel = Material.Steel()
        # cls.rho = steel.rho
        # cls.E = steel.E

        _pieces = 30  # number of finite elements
        for i in range(2):
            if i == 0:
                cls._length = 450.  # mm
                cls.rho = 7.850e-9
                section_column = sections.Rectangle(height=3., width=20.)
                cls.E = 2.1e5
            else:
                cls._length = 0.45  # m
                cls.rho = 7850
                section_column = sections.Rectangle(height=0.003, width=0.02)
                cls.E = 2.1e11

            # nodes
            _nodes = []
            for i in range(_pieces + 1):
                _nodes.append(Node.Node.from_dict(adict={'ID': i + 1, 'coords': (i * cls._length / _pieces, 0)}))

            # beams
            cls.I = section_column.I.x
            cls.A = section_column.A
            _beams = []
            for i in range(_pieces):
                _beams.append(HB.HermitianBeam2D.from_dict(adict={'ID': i, 'E': cls.E, 'I': section_column.I.x,
                                                                  'A': section_column.A, 'rho': cls.rho,
                                                                  'i': _nodes[i], 'j': _nodes[i + 1]}))

            # supports
            BCs = {1: ['ux', 'uy', 'rotz']}  # supports as dict

            # this is the structure
            cls.structure = Structure.Structure(beams=_beams, supports=BCs)

    def test_clamped_rod(self):
        # testing using the consistent mass matrix
        # note: would work with 8 elements as well

        self.structure.solver['modal'].solve()

        _base = math.sqrt((self.I * self.E) / ((self.A * self.rho) * self._length ** 4))
        omegas = [3.52, 22.03, 61.7, 121, 200]  # clamped at one end
        omegas = [x * _base for x in omegas]
        for pre, f in zip(omegas, self.structure.results['modal'].circular_frequencies[0:5]):
            self.assertAlmostEqual(pre, f, delta=pre*0.01)

    def test_clamped_rod_lumped(self):
        # testing using the lumped mass matrix
        # note: nr. of elements should be high, about 30

        self.structure.set_mass_matrix_type(matrixtype='lumped')
        self.structure.solver['modal'].solve()

        _base = math.sqrt((self.I * self.E) / ((self.A * self.rho) * self._length ** 4))
        omegas = [3.52, 22.03, 61.7, 121, 200]  # clamped at one end
        omegas = [x * _base for x in omegas]
        for pre, f in zip(omegas, self.structure.results['modal'].circular_frequencies[0:5]):
            self.assertAlmostEqual(pre, f, delta=pre*0.01)
