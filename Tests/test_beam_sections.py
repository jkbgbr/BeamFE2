# -*- coding: utf-8 -*-

import unittest
from BeamFE2.BeamSections import *


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
        self.assertEqual(cs.A, 10)
        self.assertEqual(cs.I, 35)

    def test_Rectange(self):
        cs = Rectangle(height=2, width=1.5)
        self.assertEqual(cs.A, 3)
        self.assertEqual(cs.I['x'], 1)
        self.assertEqual(cs.I['y'], 0.5625)

    def test_Circle(self):
        cs = Circle(r=1)
        self.assertEqual(cs.A, math.pi)
        self.assertEqual(cs.I['x'], math.pi/4.)
        self.assertEqual(cs.I['y'], math.pi/4.)

    def test_HollowCircle_1(self):
        # thinkness almost zero
        cs = HollowCircle(r=1, s=0.00001)
        self.assertAlmostEqual(cs.A, 0, delta=0.0001)
        self.assertAlmostEqual(cs.I['x'], 0., delta=0.0001)
        self.assertAlmostEqual(cs.I['y'], 0., delta=0.0001)

    def test_HollowCircle_2(self):
        # thickness almost r
        cs = HollowCircle(r=1, s=1-0.00001)
        print(cs.A)
        print(cs.I)
        self.assertAlmostEqual(cs.A, math.pi, delta=0.0001)
        self.assertAlmostEqual(cs.I['x'], math.pi/4., delta=0.0001)
        self.assertAlmostEqual(cs.I['y'], math.pi/4., delta=0.0001)

    def test_HollowCircle_3(self):
        # thickness almost r
        cs = HollowCircle(r=1, s=0.1)
        self.assertAlmostEqual(cs.A, 0.5969026041820605, delta=0.0001)
        self.assertAlmostEqual(cs.I['x'], 0.27009842839238246, delta=0.0001)
        self.assertAlmostEqual(cs.I['y'], 0.27009842839238246, delta=0.0001)

