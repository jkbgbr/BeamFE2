import unittest
from BeamFE2 import HermitianBeam_2D as HB
from BeamFE2 import BeamSections as sections
from BeamFE2 import Structure
from BeamFE2 import Node
from BeamFE2 import Material
from Tests import ATOL
import numpy as np


class Test_All_Results_1(unittest.TestCase):

    """
    a continuous beam
    """

    @classmethod
    def setUpClass(cls):

        _nodes = [Node.Node(ID=x+1, coords=(500 * x, 0)) for x in range(7)]
        mat = Material.Steel()
        sect = sections.Rectangle(height=3, width=10)  # section
        _beams = [HB.HermitianBeam2D(ID=x+1, E=mat.E, rho=mat.rho, I=sect.I['x'], A=sect.A, i=_nodes[x], j=_nodes[x+1]) for x in range(6)]
        _supports = {1: ['uy'],
                     3: ['uy'],
                     5: ['uy'],
                     7: ['ux', 'uy', 'rotz'],
                     }
        cls.structure = Structure.Structure(beams=_beams, supports=_supports)

        # adding nodal loads
        cls.structure.clear_loads()
        cls.structure.add_nodal_load(nodeID=1, dynam={'FX': -1000}, clear=True)
        cls.structure.add_nodal_load(nodeID=2, dynam={'FY': -1000})
        cls.structure.add_nodal_load(nodeID=4, dynam={'MZ': 500000})

        # internal loads
        b = cls.structure.beams[4]
        cls.structure.add_internal_loads(beam=b, loadtype='uniform perpendicular force', value=-10.0)
        b = cls.structure.beams[5]
        cls.structure.add_internal_loads(beam=b, loadtype='uniform perpendicular force', value=-20.0)

    def test_reactions(self):
        self.structure.solver['linear static'].solve()

        _expected_reactions = np.matrix([[0.],
                                         [-481.971154],
                                         [0.],
                                         [0.],
                                         [-0.],
                                         [0.],
                                         [0.],
                                         [-608.173077],
                                         [-0.],
                                         [0.],
                                         [0.],
                                         [0.],
                                         [0.],
                                         [-4848.557692],
                                         [0.],
                                         [0.],
                                         [0.],
                                         [-0.],
                                         [-1000.],
                                         [-10061.298077],
                                         [1739182.692308]])

        self.assertTrue(np.allclose(self.structure.results['linear static'].reaction_forces, _expected_reactions))

    def test_deformations(self):
        self.structure.solver['linear static'].solve()

        _expected_deformations = [np.matrix([[-0.47619],
                                            [-0.],
                                            [-12.591575],
                                            [-0.396825],
                                            [-4170.694275],
                                            [0.158985],
                                            [-0.31746],
                                            [-0.],
                                            [11.955637],
                                            [-0.238095],
                                            [5898.326211],
                                            [12.432591],
                                            [-0.15873],
                                            [-0.],
                                            [-35.230973],
                                            [-0.079365],
                                            [-12671.067359],
                                            [6.052011],
                                            [-0.],
                                            [-0.],
                                            [0.]])]

        self.assertTrue(np.allclose(self.structure.results['linear static'].displacement_results, _expected_deformations))

    def test_frequencies(self):
        self.structure.solver['modal'].solve()

        _expected_frequencies = [7.613434040661384, 11.09323127793537, 15.326172330326486,
                                 32.693985190216615, 41.63737293315384, 54.555488726848836,
                                 84.23566717176273, 110.95431530916404, 138.32201256527927, 432.24810222320696,
                                 1326.485590168174, 2310.0531736149446, 3428.9430318147342, 4633.931252454926,
                                 5560.428871743147]

        self.assertTrue(np.allclose(self.structure.results['modal'].frequencies, _expected_frequencies, atol=1e-5))

    def test_frequencies_with_mass(self):
        self.structure.add_nodal_mass(nodeID=4, mass={'mx': 1, 'my': 1})
        self.structure.solver['modal'].solve()
        _expected_frequencies = [0.10469612329002663, 9.231045221012806, 10.311990761936627,
                                 13.586565471302825, 32.34979734510168]
        self.assertTrue(np.allclose(self.structure.results['modal'].frequencies[0:5], _expected_frequencies, atol=1e-5))

    # def test_plotall(self):
    #     self.structure.draw(analysistype='linear static')
    #     self.structure.draw(analysistype='linear static', internal_action='shear')
    #     self.structure.draw(analysistype='linear static', internal_action='moment')
    #     self.structure.draw(analysistype='modal', mode=0)


class Test_All_Results_2(unittest.TestCase):

    """
    a portal frame to test reaction forces, modal results
    """

    @classmethod
    def setUpClass(cls):

        _nodes = {1: Node.Node(ID=1, coords=(0, 0)),
                  2: Node.Node(ID=2, coords=(0, 2000)),
                  3: Node.Node(ID=3, coords=(2000, 4000)),
                  4: Node.Node(ID=4, coords=(5000, 4000)),
                  5: Node.Node(ID=5, coords=(5000, 3000)),
                  6: Node.Node(ID=6, coords=(5000, 2000)),
                  7: Node.Node(ID=7, coords=(5000, 1000)),
                  8: Node.Node(ID=8, coords=(5000, 0)),
                  9: Node.Node(ID=9, coords=(5000, -1000)),
                  }

        mat = Material.Steel()
        s1 = sections.Rectangle(height=100, width=100)  # section
        s2 = sections.Rectangle(height=200, width=200)  # section
        _sections = [s1, s1, s2, s1, s1, s1, s1, s1]
        _beams = [HB.HermitianBeam2D(ID=x+1, E=mat.E, rho=mat.rho, I=_sections[x].I['x'], A=_sections[x].A, i=_nodes[y], j=_nodes[y+1]) for x, y in zip(range(8), range(1, len(_nodes.keys())))]
        _supports = {1: ['ux', 'uy'],
                     9: ['ux', 'uy'],
                     }
        cls.structure = Structure.Structure(beams=_beams, supports=_supports)

        # adding nodal loads
        cls.structure.clear_loads()
        cls.structure.add_nodal_load(nodeID=3, dynam={'FX': 1000, 'FY': 500}, clear=True)
        cls.structure.add_nodal_load(nodeID=5, dynam={'FX': -1500})

        # internal loads
        b = cls.structure.beams[1]
        cls.structure.add_internal_loads(beam=b, loadtype='uniform perpendicular force', value=2.0)

    def test_reactions(self):
        self.structure.solver['linear static'].solve()

        _expected_reactions = np.matrix([[-2848.39394],
                                         [669.678788],
                                         [-1651.60606],
                                         [3830.321212]])
        refo = self.structure.results['linear static'].reaction_forces
        # _relevant_forces = np.concatenate(refo[0:2, 0], refo[-3:-1, 0])
        self.assertTrue(np.allclose(refo[0:2], _expected_reactions[0:2, 0]))
        self.assertTrue(np.allclose(refo[-3:-1], _expected_reactions[2:4, 0]))

    def test_deformations(self):
        self.structure.solver['linear static'].solve()

        _expected_deformations = [np.matrix([[0.],
                                             [0.],
                                             [0.011805],
                                             [-21.439837],
                                             [0.000638],
                                             [0.00855],
                                             [-27.697614],
                                             [6.255481],
                                             [-0.002034],
                                             [-27.69756],
                                             [0.00912],
                                             [-0.001925],
                                             [-27.706412],
                                             [0.007296],
                                             [0.001893],
                                             [-24.083021],
                                             [0.005472],
                                             [0.005196],
                                             [-17.628305],
                                             [0.003648],
                                             [0.007556],
                                             [-9.28604],
                                             [0.001824],
                                             [0.008971],
                                             [0.],
                                             [0.],
                                             [0.009443]])]

        self.assertTrue(np.allclose(self.structure.results['linear static'].displacement_results, _expected_deformations, atol=1e-5))

    def test_frequencies(self):
        self.structure.solver['modal'].solve()

        _expected_frequencies = np.matrix([1.8098134635564367, 11.886156983577484, 15.868886208391775, 46.8696497539431,
                                           56.91143213363919,82.10962203501121, 101.00378682424322, 146.20319107194706,
                                           172.68133458649052, 191.3400650417073, 278.6082063770441, 314.0478965189491,
                                           347.15683024183244, 429.28119919168194, 606.6507091264467, 626.4477915106148,
                                           705.2507665620494, 873.2857713706983, 922.718433097769, 1147.3740125153104,
                                           1231.6447012112453, 1861.2153655737682, 2520.6433978302216]
                                          )

        self.assertTrue(np.allclose(self.structure.results['modal'].frequencies, _expected_frequencies, atol=1e-5))

    def test_frequencies_with_mass(self):
        self.structure.add_nodal_mass(nodeID=4, mass={'mx': 15})
        self.structure.solver['modal'].solve()
        _expected_frequencies = [0.531783218671432, 11.802244251990954, 15.313592559544537,
                                 46.307812940599206, 56.707579685571524]
        self.assertTrue(np.allclose(self.structure.results['modal'].frequencies[0:5], _expected_frequencies, atol=1e-5))

    def test_plotall(self):
        self.structure.draw(analysistype='linear static')
        self.structure.draw(analysistype='linear static', internal_action='axial')
        self.structure.draw(analysistype='linear static', internal_action='shear')
        self.structure.draw(analysistype='linear static', internal_action='moment')
        self.structure.draw(analysistype='modal', mode=3)


class Test_All_Results_3(unittest.TestCase):

    """
    a clamped beam to test reaction forces. essentially the same structure as in #1 but no vertical supports
    """

    @classmethod
    def setUpClass(cls):

        _nodes = [Node.Node(ID=x+1, coords=(500 * x, 0)) for x in range(7)]
        mat = Material.Steel()
        sect = sections.Rectangle(height=3, width=10)  # section
        _beams = [HB.HermitianBeam2D(ID=x+1, E=mat.E, rho=mat.rho, I=sect.I['x'], A=sect.A, i=_nodes[x], j=_nodes[x+1]) for x in range(6)]
        _supports = {7: ['ux', 'uy', 'rotz']}
        cls.structure = Structure.Structure(beams=_beams, supports=_supports)

        # adding nodal loads
        cls.structure.clear_loads()
        cls.structure.add_nodal_load(nodeID=1, dynam={'FY': -1000}, clear=True)

    def test_reactions(self):
        self.structure.solver['linear static'].solve()
        _expected_reactions = np.matrix([[      -0.],
                                         [      -0.],
                                         [      -0.],
                                         [      -0.],
                                         [       0.],
                                         [      -0.],
                                         [      -0.],
                                         [      -0.],
                                         [       0.],
                                         [      -0.],
                                         [       0.],
                                         [      -0.],
                                         [      -0.],
                                         [       0.],
                                         [       0.],
                                         [      -0.],
                                         [      -0.],
                                         [       0.],
                                         [      -0.],
                                         [   -1000.],
                                         [ 3000000.]])

        # todo: the moment reactions are ALWAYS incxorrect
        self.assertTrue(np.allclose(self.structure.results['linear static'].reaction_forces, _expected_reactions, atol=1e-5))

    def test_plotall(self):
        self.structure.draw(analysistype='linear static')
        self.structure.draw(analysistype='linear static', internal_action='axial')
        self.structure.draw(analysistype='linear static', internal_action='shear')
        self.structure.draw(analysistype='linear static', internal_action='moment')
        self.structure.draw(analysistype='modal', mode=3)


class Test_All_Results_4(unittest.TestCase):

    """
    a clamped beam to test reaction forces. essentially the same structure as in #3 but with multiple loads
    """

    @classmethod
    def setUpClass(cls):

        _nodes = [Node.Node(ID=x+1, coords=(500 * x, 0)) for x in range(7)]
        mat = Material.Steel()
        sect = sections.Rectangle(height=3, width=10)  # section
        _beams = [HB.HermitianBeam2D(ID=x+1, E=mat.E, rho=mat.rho, I=sect.I['x'], A=sect.A, i=_nodes[x], j=_nodes[x+1]) for x in range(6)]
        _supports = {7: ['ux', 'uy', 'rotz']}
        cls.structure = Structure.Structure(beams=_beams, supports=_supports)

        # adding nodal loads
        cls.structure.clear_loads()
        cls.structure.add_nodal_load(nodeID=1, dynam={'FY': -1000, 'FX': -1000}, clear=True)
        cls.structure.add_internal_loads(beam=cls.structure.beams[0], loadtype='uniform axial force', value=1)
        cls.structure.add_internal_loads(beam=cls.structure.beams[2], loadtype='concentrated axial force', value=-1000, position=0.4)
        cls.structure.add_internal_loads(beam=cls.structure.beams[4], loadtype='concentrated axial force', value=1000, position=0.99)

    def test_reactions(self):
        self.structure.solver['linear static'].solve()
        _expected_reactions = np.matrix([[      -0.],
                                         [      -0.],
                                         [      -0.],
                                         [      -0.],
                                         [       0.],
                                         [      -0.],
                                         [      -0.],
                                         [      -0.],
                                         [       0.],
                                         [      -0.],
                                         [       0.],
                                         [      -0.],
                                         [      -0.],
                                         [       0.],
                                         [       0.],
                                         [      -0.],
                                         [      -0.],
                                         [       0.],
                                         [    -500.],
                                         [   -1000.],
                                         [ 3000000.]])

        # todo: the moment reactions are ALWAYS incxorrect
        self.assertTrue(np.allclose(self.structure.results['linear static'].reaction_forces, _expected_reactions, atol=1e-5))

    def test_internal_actions(self):
        # tests the query of the internal actions
        self.structure.solver['linear static'].solve()
        beam = self.structure.beams[0]
        disp = self.structure.results['linear static'].element_displacements(local=True, beam=beam, asvector=True)
        self.assertAlmostEqual(beam.internal_action(disp=disp, action='axial', pos=0.), 1000, delta=1e-8)
        self.assertAlmostEqual(beam.internal_action(disp=disp, action='axial', pos=0.5), 750, delta=1e-8)
        self.assertAlmostEqual(beam.internal_action(disp=disp, action='axial', pos=1.), 500, delta=1e-8)
        self.assertAlmostEqual(beam.internal_action(disp=disp, action='moment', pos=0.), 1.421431e-07, delta=1e-8)
        self.assertAlmostEqual(beam.internal_action(disp=disp, action='shear', pos=0.5), 1000, delta=1e-8)
        self.assertRaises(AssertionError, beam.internal_action, {'action': 'axial', 'pos': 1.2, 'disp': disp})
        self.assertRaises(AssertionError, beam.internal_action, {'action': 'woos', 'pos': 0.2, 'disp': disp})


    def test_plotall(self):
        self.structure.solver['linear static'].solve()
        self.structure.draw(analysistype='linear static')
        self.structure.draw(analysistype='linear static', internal_action='axial')
        self.structure.draw(analysistype='linear static', internal_action='shear')
        self.structure.draw(analysistype='linear static', internal_action='moment')
        self.structure.draw(analysistype='modal', mode=3)
