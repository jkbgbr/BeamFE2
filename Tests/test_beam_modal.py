import unittest
from Beams import HermitianBeam as HB
from Beams import Sections as sections

from drawing import draw_beam
import numpy as np
import math
from Tests import ATOL


class Hermitian2D_Model(unittest.TestCase):
    """
    tests for the modal analysis of the FE Beams
    """

    @classmethod
    def setUpClass(cls):
        cls.rho = 7850000  # g/m3
        cls.E = 2.1e11  # N/m2

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
            _nodes.append(HB.Node.from_dict(adict={'ID': i + 1, 'coords': (i * _length / _pieces, 0)}))

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
        structure = HB.Structure(beams=_beams, supports=BCs)

        # adding loads
        # structure.add_single_dynam_to_node(nodeID=len(_nodes) - 1, dynam={'FY': -1000000}, clear=True)  # clears previous loads

        # solver :-) whatever happens here is done by numpy.
        structure.solve(analysis='modal')

        # pre-calculated values are from Axis, except for f3, that should be 708.3 Hz
        for pre, f in zip([40.56, 253.55, 711.3, 923.26], structure.frequencies[0:4]):
            self.assertAlmostEqual(pre, f, delta=2.)
