import unittest
from Modell import HermitianBeam_2D as HB
from Modell import BeamSections as sections
from Modell import Structure
from Modell import Node
from Modell import Material
from solver import solve
import math

from Tests import ATOL


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
        A beam with clamped end at Node #1, under self weight.
        L = 1400 mm,
        cross-section: r=55 mm rod
        steel: rho = 7.85 t/m3, E=210 GPa
        """

        _pieces = 9  # number of segments along the length
        _length = 1400  # mm

        # nodes
        _nodes = []
        for i in range(_pieces + 1):
            _nodes.append(Node.Node.from_dict(adict={'ID': i + 1, 'coords': (i * _length / _pieces, 0)}))

        # beams
        section_column = sections.Circle(r=55)
        _beams = []
        for i in range(_pieces):
            _beams.append(HB.HermitianBeam2D.from_dict(adict={'ID': i, 'E': self.E, 'I': section_column.I['x'],
                                                              'A': section_column.A, 'rho': self.rho,
                                                              'i': _nodes[i], 'j': _nodes[i + 1]}))

        # supports
        BCs = {1: ['ux', 'uy', 'rotz']}  # supports as dict

        # this is the structure
        structure = Structure.Structure(beams=_beams, supports=BCs)

        # adding loads
        # structure.add_single_dynam_to_node(nodeID=len(_nodes) - 1, dynam={'FY': -1000000}, clear=True)  # clears previous loads

        # solver :-) whatever happens here is done by numpy.
        solve(structure, analysis='modal')

        # pre-calculated values are from Axis, except for f3, that should be 708.3 Hz
        for pre, f in zip([40.56, 253.55, 711.3, 923.26], structure.results['modal'].frequencies[0:4]):
            self.assertAlmostEqual(pre, f, delta=2.)

    def test_doubly_clamped_rod(self):
        """
        A beam with clamped ends under self weight.
        """

        _pieces = 9  # number of segments along the length
        _length = 1  # mm
        rho = 7833

        # nodes
        _nodes = []
        for i in range(_pieces + 1):
            _nodes.append(Node.Node.from_dict(adict={'ID': i + 1, 'coords': (i * _length / _pieces, 0)}))

        # beams
        # section_column = sections.Circle(r=(55/1000.))
        I = 1/1e12
        A = 1/1e6
        E = 2.0684e11
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

        # adding loads
        # structure.add_single_dynam_to_node(nodeID=len(_nodes) - 1, dynam={'FY': -1000000}, clear=True)  # clears previous loads

        # solver :-) whatever happens here is done by numpy.
        solve(structure, analysis='modal')

        print(structure.results['modal'].circular_frequencies)
        print(structure.results['modal'].frequencies)
        print(I)
        print(E)
        print(rho)
        print(_length)
        q = rho * A * _length / _length
        mu = q / 9.81

        _grund = math.sqrt((I * E) / (mu * _length ** 4))
        # om = [22.4, 61.7, 121, 200, 298]  # doubly clamped
        om = [3.52, 22.4, 61.7, 121, 200]  # clamped at one end
        om = [x * _grund for x in om]
        print(om)
