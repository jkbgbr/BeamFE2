# -*- coding: utf-8 -*-

import unittest
from BeamFE2.Material import *


class Material_Test(unittest.TestCase):
    """
    basic tests for the materials. ONLY the engineering units are tested.
    """

    def test_LinearElasticMaterial(self):
        # tests the assertions
        args = []
        kwargslist = []
        kwargslist.append({'name': 'hooke', 'E': 25, 'nu': 0.4, 'rho': -0.1})  # negative density
        kwargslist.append({'name': 'hooke', 'E': 0, 'nu': 0.4, 'rho': -0.1})  # E-module zero
        kwargslist.append({'name': 'hooke', 'E': 25, 'nu': 0.6, 'rho': 0.1})  # Poisson to large
        kwargslist.append({'name': 'hooke', 'E': 25, 'nu': 0, 'rho': 0.1})  # Poisson to small
        for kw in kwargslist:
            self.assertRaises(AssertionError, LinearElasticMaterial, *args, **kw)

    def test_Steel(self):
        steel = Steel()
        self.assertTrue(steel.E == 2.1e5)
        self.assertAlmostEqual(steel.G, 9.13e4, delta=10)
        self.assertTrue(steel.nu == 0.3)
        self.assertTrue(steel.rho == 7.85e-9)
        self.assertTrue(steel.name == 'steel')
        self.assertFalse(steel.E == 0)
        self.assertTrue(repr(steel) == 'Steel(E=210000.0, nu=0.3, rho=7.85e-09)')

    def test_set_rigid(self):
        steel = Steel()
        steel.set_rigid()
        self.assertFalse(steel.E == 2.1e5)
        self.assertTrue(steel.E == 1e40)
        self.assertTrue('rigid' in steel.name)
        steel.set_rigid()
        self.assertTrue(steel.name.count('rigid') == 1)

    def test_set_weightless(self):
        steel = Steel()
        steel.set_weightless()
        self.assertFalse(steel.rho == 7.85e-9)
        self.assertTrue(steel.rho == 0)
        self.assertTrue('weightless' in steel.name)
        steel.set_weightless()
        self.assertTrue(steel.name.count('weightless') == 1)

