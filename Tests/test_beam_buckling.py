import unittest
from BeamFE2 import HermitianBeam_2D as HB
from BeamFE2 import BeamSections as sections
from BeamFE2 import Structure
from BeamFE2 import Node
from BeamFE2 import Material
import math


F_HORIZONTAL = 1
F_VERTICAL = 0


class Hermitian2D_Model(unittest.TestCase):
    """
    tests for the buckling analysis of the FE BeamFE2
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
        _pieces = 1  # number of finite elements
        _length = 850.  # mm
        mat = Material.Steel()
        section_column = sections.Rectangle(height=1., width=1.)

        # nodes
        _nodes = []
        for i in range(_pieces + 1):
            _nodes.append(Node.Node.from_dict(adict={'ID': i + 1, 'coords': (i * _length / _pieces, 0)}))

        # beams
        I = section_column.I['x']
        A = section_column.A
        I = 1
        _beams = []
        for i in range(_pieces):
            _beams.append(HB.HermitianBeam2D.from_dict(adict={'ID': i, 'E': mat.E, 'I': I,
                                                              'A': A, 'rho': mat.rho,
                                                              'i': _nodes[i], 'j': _nodes[i + 1]}))

        # supports
        _last = max([x.ID for x in _nodes])
        BCs = {1: ['ux', 'uy', 'rotz']}  # supports as dict

        # this is the structure
        structure = Structure.Structure(beams=_beams, supports=BCs)

        # adding loads
        structure.add_single_dynam_to_node(nodeID=_last, dynam={'FX': -1 * F_HORIZONTAL, 'FY': -1*F_VERTICAL}, clear=True)

        # solver :-) whatever happens here is done by numpy.

        # beam = structure.beams[0]
        # print(math.pi**2 * beam.EI / (2 * _length) ** 2)
        # L = beam.l
        # print(beam.E)
        # print(beam.A)
        # print(beam.l)
        # print(beam.EA/L)
        # print(12*beam.EI/L**3)

        structure.solver['linear static'].solve()
        structure.solver['modal'].solve()
        structure.solver['buckling'].solve()
