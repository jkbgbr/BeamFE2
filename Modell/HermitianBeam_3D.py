# -*- coding: utf-8 -*-

import copy
import math
import numpy as np


class HermitianBeam3D(object):
    """
    A hermitian 3D FE Beam with 2 nodes and 6 DOFs each node: 3 translational, 3 rotational.
    A 2D version is also available with 3 DOFs each node: 2 translational, 1 rotational.
    Small deformations are assumed.
    Bernoulli-Euler beam theory (no shear deformations)
    Non-restrained torsion only (no warping)

    Nodes: i, j with direction i -> j
    Element coordinate system is right-handed, thumb: X, pointing finger Y, middle finger Z
    Displacements are positive when going away from the origin.
    Rotations are positive when clockwise when looking away from the origin.

    Displacement vector u:
    u1 : Node i, displacement in X direction
    u2 : Node i, displacement in Y direction
    u3 : Node i, displacement in Z direction
    u4 : Node i, rotation about X axis
    u5 : Node i, rotation about Y axis
    u6 : Node i, rotation about Z axis
    u7 - u12: same for Node j

    """

    dof = 6  # degrees of freedom
    dofnames = 'ux', 'uy', 'uz', 'rotx', 'roty', 'rotz'
    number_internal_points = 20

    def __init__(self, ID=None, i=None, j=None, E=None, nu=None, crosssection=None, Ix=None, Iy=None, A=None, rho=None):
        self.ID = ID
        self.E = E
        self.nu = nu
        self.i = i
        self.j = j
        if crosssection is not None:
            self.A = crosssection.A
            self.Ix = crosssection.I['x']
            self.Iy = crosssection.I['y']
        else:
            self.A = A
            self.Ix = Ix
            self.Iy = Iy
        self._end_DOFs = {'i': list(copy.deepcopy(self.dofnames)), 'j': list(copy.deepcopy(self.dofnames))}

        self.rho = rho
        self.displacements = None  # displacements in the GLOBAL system
    @property
    def l(self):
        return math.sqrt((self.i[0] - self.j[0]) ** 2 + (self.i[1] - self.j[1]) ** 2)

    @property
    def G(self):
        return self.E * self.nu

    @property
    def Ke(self):
        l = self.l
        E_x = self.E * self.A
        G_x = self.G * self.A
        EIz_1 = EIz_2 = self.E * self.Ix
        EIy_1 = EIy_2 = self.E * self.Iy

        _ret = np.array([
            [E_x / l, 0, 0, 0, 0, 0, -E_x / l, 0, 0, 0, 0, 0],
            [0, 6 * (EIz_1 + EIz_2) / (l ** 3), 0, 0, 0, 2 * (2 * EIz_1 + EIz_2) / (l ** 2), 0,
             -6 * (EIz_1 + EIz_2) / (l ** 3), 0, 0, 0, 2 * (EIz_1 + 2 * EIz_2) / (l ** 2)],
            [0, 0, 6 * (EIy_1 + EIy_2) / (l ** 3), 0, -2 * (2 * EIy_1 + EIy_2) / (l ** 2), 0, 0, 0,
             -6 * (EIy_1 + EIy_2) / (l ** 3), 0, -2 * (EIy_1 + 2 * EIy_2) / (l ** 2), 0],
            [0, 0, 0, G_x / l, 0, 0, 0, 0, 0, -G_x / l, 0, 0],
            [0, 0, -2 * (2 * EIy_1 + EIy_2) / (l ** 2), 0, (3 * EIy_1 + EIy_2) / l, 0, 0, 0,
             2 * (2 * EIy_1 + EIy_2) / (l ** 2), 0, (EIy_1 + EIy_2) / l, 0],
            [0, 2 * (2 * EIz_1 + EIz_2) / (l ** 2), 0, 0, 0, (3 * EIz_1 + EIz_2) / l, 0,
             -2 * (2 * EIz_1 + EIz_2) / (l ** 2), 0, 0, 0, (EIz_1 + EIz_2) / l],
            [-E_x / l, 0, 0, 0, 0, 0, E_x / l, 0, 0, 0, 0, 0],
            [0, -6 * (EIz_1 + EIz_2) / (l ** 3), 0, 0, 0, -2 * (2 * EIz_1 + EIz_2) / (l ** 2), 0,
             6 * (EIz_1 + EIz_2) / (l ** 3), 0, 0, 0, -2 * (EIz_1 + 2 * EIz_2) / (l ** 2)],
            [0, 0, -6 * (EIy_1 + EIy_2) / (l ** 3), 0, 2 * (2 * EIy_1 + EIy_2) / (l ** 2), 0, 0, 0,
             6 * (EIy_1 + EIy_2) / (l ** 3), 0, 2 * (EIy_1 + 2 * EIy_2) / (l ** 2), 0],
            [0, 0, 0, -G_x / l, 0, 0, 0, 0, 0, G_x / l, 0, 0],
            [0, 0, -2 * (EIy_1 + 2 * EIy_2) / (l ** 2), 0, (EIy_1 + EIy_2) / l, 0, 0, 0,
             2 * (EIy_1 + 2 * EIy_2) / (l ** 2), 0, (EIy_1 + 3 * EIy_2) / l, 0],
            [0, 2 * (EIz_1 + 2 * EIz_2) / (l ** 2), 0, 0, 0, (EIz_1 + EIz_2) / l, 0,
             -2 * (EIz_1 + 2 * EIz_2) / (l ** 2), 0, 0, 0, (EIz_1 + 3 * EIz_2) / l]
        ])
        return _ret.astype(int)
