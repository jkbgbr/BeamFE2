# -*- coding: utf-8 -*-

import unittest
from Modell.BeamSections import *


class Material_Test(unittest.TestCase):
    """
    basic tests for the materials. ONLY the engineering units are tested.
    """

    def test_Crosssection(self):
        # tests the assertions
        args = []
        kwargslist = []
        kwargslist.append({'shape': 'custom', 'A': 0, 'I': 1})  # A zero
        kwargslist.append({'shape': 'custom', 'A': 10, 'I': -1})  # I negative
        kwargslist.append({'shape': 'custom', 'I': -1})  # A is None
        kwargslist.append({'shape': 'custom', 'A': 1})  # I is None
        for kw in kwargslist:
            self.assertRaises(AssertionError, CustomCrosssection, *args, **kw)

        cs = CustomCrosssection(shape='someshape', A=10, I=35)
        self.assertEqual(cs._A, 10)
        self.assertEqual(cs._I, 35)