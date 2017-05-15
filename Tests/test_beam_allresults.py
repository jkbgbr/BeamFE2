import unittest
from Modell import HermitianBeam_2D as HB
from Modell import BeamSections as sections
from Modell import Structure
from Modell import Node
from Modell import Material
from solver import solve

import numpy as np


class Test_All_Results_1(unittest.TestCase):

    """
    a continuous beam to test reaction forces, modal results
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

    def test_do(self):
        solve(structure=self.structure, analysis='all')

        _expected_reactions = np.matrix([[0.],
                                         [481.971154],
                                         [0.],
                                         [0.],
                                         [-0.],
                                         [0.],
                                         [0.],
                                         [608.173077],
                                         [-0.],
                                         [0.],
                                         [0.],
                                         [0.],
                                         [0.],
                                         [4848.557692],
                                         [0.],
                                         [0.],
                                         [0.],
                                         [-0.],
                                         [1000.],
                                         [10061.298077],
                                         [-1739182.692308]])

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

        _expected_frequencies = [7.613434040661384, 11.09323127793537, 15.326172330326486,
                                 32.693985190216615, 41.63737293315384, 54.555488726848836,
                                 84.23566717176273, 110.95431530916404, 138.32201256527927, 432.24810222320696,
                                 1326.485590168174, 2310.0531736149446, 3428.9430318147342, 4633.931252454926,
                                 5560.428871743147]

        self.assertTrue(np.allclose(self.structure.results['modal'].frequencies, _expected_frequencies))
        self.assertTrue(np.allclose(self.structure.results['linear static'].displacement_results, _expected_deformations))
        self.assertTrue(np.allclose(self.structure.results['linear static'].reaction_forces, _expected_reactions))
