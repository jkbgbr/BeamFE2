import unittest
from Modell import HermitianBeam_2D as HB
from Modell import BeamSections as sections
from Modell import Structure
from Modell import Node
from Modell import Material
from solver import solve
import math


class Hermitian2D_Model(unittest.TestCase):
    """
    tests for the modal analysis of the FE Modell
    """

    @classmethod
    def setUpClass(cls):
        steel = Material.Steel()
        cls.rho = steel.rho
        cls.E = steel.E

    def test_clamped_rod(self):
        """
        A beam with clamped end under self weight.
        Results are accepted if circular frequencies are within 1% of the analytical solution.
        http://iitg.vlab.co.in/?sub=62&brch=175&sim=1080&cnt=1
        """
        _pieces = 8  # number of finite elements
        for i in range(2):
            if i == 0:
                _length = 450.  # mm
                rho = 7.850e-8
                section_column = sections.Rectangle(height=3., width=20.)
                E = 2.1e6
            else:
                _length = 0.45  # m
                rho = 7850
                section_column = sections.Rectangle(height=0.003, width=0.02)
                E = 2.1e11

            # nodes
            _nodes = []
            for i in range(_pieces + 1):
                _nodes.append(Node.Node.from_dict(adict={'ID': i + 1, 'coords': (i * _length / _pieces, 0)}))

            # beams
            I = section_column.I['x']
            A = section_column.A
            _beams = []
            for i in range(_pieces):
                _beams.append(HB.HermitianBeam2D.from_dict(adict={'ID': i, 'E': E, 'I': I,
                                                                  'A': A, 'rho': rho,
                                                                  'i': _nodes[i], 'j': _nodes[i + 1]}))

            # supports
            _last = max([x.ID for x in _nodes])
            BCs = {1: ['ux', 'uy', 'rotz']}  # supports as dict

            # this is the structure
            structure = Structure.Structure(beams=_beams, supports=BCs)

            # solver :-) whatever happens here is done by numpy.
            solve(structure, analysis='modal')

            _base = math.sqrt((I * E) / ((A * rho) * _length ** 4))
            omegas = [3.52, 22.03, 61.7, 121, 200]  # clamped at one end
            omegas = [x * _base for x in omegas]
            for pre, f in zip(omegas, structure.results['modal'].circular_frequencies[0:5]):
                self.assertAlmostEqual(pre, f, delta=pre*0.01)
