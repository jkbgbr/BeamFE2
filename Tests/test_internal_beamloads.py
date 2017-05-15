import unittest
from Modell import HermitianBeam_2D as HB
from Modell import BeamSections as sections
from Modell import Structure
from Modell import Node
from Modell import Material
from solver import solve
from drawing import _plotting_available, plt
import numpy as np
import math
from Tests import ATOL

VERTICAL = False  # True/False for vertical/horizontal
NR_BEAMS = 1  # number of finite elements
LENGTH = 100  # length of cantilever
F_HORIZONTAL = 0
F_VERTICAL = 1000


class Single_Beam_Internal_Loads(unittest.TestCase):
    """
    One-beam structure with various loads to test the calculated distributions of the internal actions
    """

    @classmethod
    def setUpClass(cls):

        # nodes
        _n1 = Node.Node(ID=1, coords=(0, 0))
        _n2 = Node.Node(ID=2, coords=(100, 0))

        # beams
        section_column = sections.Rectangle(height=3, width=10)  # section
        mat = Material.Steel()
        _b1 = HB.HermitianBeam2D.from_dict(adict={'ID': 1, 'E': mat.E, 'I': section_column.I['x'], 'A': section_column.A, 'rho': mat.rho, 'i': _n1, 'j': _n2})
        _b1.number_internal_points = 4  # setting this to have much less values for comparison

        # supports
        BCs = {1: ['ux', 'uy', 'rotz'], NR_BEAMS+1: ['ux', 'uy', 'rotz']}  # supports as dict

        # this is the cantilever itself, composed of the beams, complete with supports
        cls.structure = Structure.Structure(beams=[_b1], supports=BCs)

        # adding nodal loads
        # directly defined nodal loads
        # structure.add_nodal_load(nodeID=2, dynam={'FX': F_HORIZONTAL, 'FY': F_VERTICAL}, clear=True)
        # structure.add_nodal_load(nodeID=3, dynam={'FX': -F_HORIZONTAL, 'FY': -F_VERTICAL})

    def test_uniform_perpendicular_force_bending(self):
        # beam internal loads
        self.structure.clear_loads()
        for b in self.structure.beams:
            self.structure.add_internal_loads(beam=b, loadtype='uniform perpendicular force', value=-1.0)
        solve(self.structure, analysis='linear static')

        # moment distribution - values on the hinged-hinged beam to modify endnode-based values
        for b in self.structure.beams:
            _contour = tuple(b.internal_action_distribution(action='moment'))
            _expected = tuple([0.0, -937.5, -1250.0, -937.5, 0.0])
            self.assertEqual(_contour, _expected)

    def test_uniform_perpendicular_force_shear(self):
        # beam internal loads
        self.structure.clear_loads()
        for b in self.structure.beams:
            self.structure.add_internal_loads(beam=b, loadtype='uniform perpendicular force', value=-1.0)
        solve(self.structure, analysis='linear static')

        # shear force distribution - values on the hinged-hinged beam to modify endnode-based values
        for b in self.structure.beams:
            _contour = tuple(b.internal_action_distribution(action='shear'))
            _expected = tuple([0.0, 25.0, 50.0, 75.0, 100.0])
            self.assertEqual(_contour, _expected)

    def test_uniform_perpendicular_force_reactions(self):
        # beam internal loads
        self.structure.clear_loads()
        for b in self.structure.beams:
            self.structure.add_internal_loads(beam=b, loadtype='uniform perpendicular force', value=-1.0)
        solve(self.structure, analysis='linear static')

        # nodal reaction forces
        b = self.structure.beams[0]
        disp = self.structure.results['linear static'].element_displacements(local=True, beam=b, asvector=True)
        _soll = np.matrix([[0.], [50.], [833.33333333], [0.], [50.], [-833.33333333]])
        self.assertTrue(np.allclose(b.nodal_reactions_asvector(disps=disp), _soll))

    def test_concentrated_perpendicular_force_bending(self):
        # beam internal loads
        self.structure.clear_loads()
        self.structure.add_internal_loads(
            beam=self.structure.beams[0], loadtype='concentrated perpendicular force', value=-30.00, position=0.3)
        solve(self.structure, analysis='linear static')

        # moment distribution
        for b in self.structure.beams:
            _contour = tuple(b.internal_action_distribution(action='moment'))
            _expected = tuple([0.0, -630., -630., 0.0])
            self.assertTrue(np.allclose(_contour, _expected))

    def test_concentrated_perpendicular_force_shear(self):
        # beam internal loads
        self.structure.clear_loads()
        self.structure.add_internal_loads(
            beam=self.structure.beams[0], loadtype='concentrated perpendicular force', value=-30.00, position=0.3)
        solve(self.structure, analysis='linear static')

        # shear force distribution
        for b in self.structure.beams:
            _contour = tuple(b.internal_action_distribution(action='shear'))
            _expected = tuple([0, 0, 30.0, 30.0])
            self.assertTrue(np.allclose(_contour, _expected))

    def test_concentrated_perpendicular_force_reactions(self):
        # beam internal loads
        self.structure.clear_loads()
        self.structure.add_internal_loads(
            beam=self.structure.beams[0], loadtype='concentrated perpendicular force', value=-30.00, position=0.3)
        solve(self.structure, analysis='linear static')

        # nodal reaction forces
        b = self.structure.beams[0]
        disp = self.structure.results['linear static'].element_displacements(local=True, beam=b, asvector=True)
        _soll = np.matrix([[0.], [23.52], [441.], [0.], [6.48], [-189.]])
        self.assertTrue(np.allclose(b.nodal_reactions_asvector(disps=disp), _soll))

    def test_concentrated_moment_bending(self):
        # beam internal loads
        self.structure.clear_loads()
        self.structure.add_internal_loads(
            beam=self.structure.beams[0], loadtype='concentrated moment', value=500, position=0.3)
        solve(self.structure, analysis='linear static')

        # moment distribution
        for b in self.structure.beams:
            _contour = tuple(b.internal_action_distribution(action='moment'))
            _expected = tuple([0.0, -149.99999999995, 349.99999999995, 0.0])
            self.assertTrue(np.allclose(_contour, _expected))

    def test_concentrated_moment_shear(self):
        # beam internal loads
        self.structure.clear_loads()
        self.structure.add_internal_loads(
            beam=self.structure.beams[0], loadtype='concentrated moment', value=500, position=0.3)
        solve(self.structure, analysis='linear static')

        # shear force distribution
        for b in self.structure.beams:
            _contour = tuple(b.internal_action_distribution(action='shear'))
            _expected = tuple([0., 0., 0., 0.])
            self.assertTrue(np.allclose(_contour, _expected))

    def test_concentrated_moment_reactions(self):
        # beam internal loads
        self.structure.clear_loads()
        self.structure.add_internal_loads(
            beam=self.structure.beams[0], loadtype='concentrated moment', value=500, position=0.3)
        solve(self.structure, analysis='linear static')

        # nodal reaction forces
        b = self.structure.beams[0]
        disp = self.structure.results['linear static'].element_displacements(local=True, beam=b, asvector=True)
        _soll = np.matrix([[0.], [6.3], [-35.], [0.], [-6.3], [165.]])
        self.assertTrue(np.allclose(b.nodal_reactions_asvector(disps=disp), _soll))
