from __future__ import division
import numpy as np
import math
import pprint as pp

EE = 210000.
GG = 81000.

topol = {1: (0, 0),
         2: (1, 0),
         3: (2, 0)}


class Node(object):
    """
    
    
    """

    def __init__(self, ID=None, coords=()):
        self.ID = ID
        self.coords = coords
        self.DOFs = {'ux': True, 'uy': True, 'uz': True}

    @property
    def x(self):
        return self.coords[0]

    @property
    def y(self):
        return self.coords[1]


class HermitianBeam2D(object):
    """
    A hermitian 3D FE Beam with 2 nodes and 3 DOFs each node: 2 translational, 1 rotational.
    Small deformations are assumed.
    Bernoulli-Euler beam theory is used (no shear deformations)
    Non-restrained torsion only (no warping)
    
    Nodes: i, j with direction i -> j
    Element coordinate system is right-handed, thumb: X, pointing finger Y, middle finger Z
    Displacements are positive when going away from the origin.
    Rotations are positive when clockwise when looking away from the origin. 
    
    Displacement vector u:
    u1 : Node i, displacement in X direction
    u2 : Node i, displacement in Y direction
    u3 : Node i, rotation about Z axis
    u4 - u6: same for Node j
    
    The stiffness matrix is composed as a numpy block matrix to make the compilation of the global stiffness
    matrix easier. However, due to this the load and displacement "vectors" also need to be numpy matrices.
    
    """
    dof = 3  # degrees of freedom

    def __init__(self, ID=None, i=None, j=None, E=EE, I=None, A=None):
        """
        :param ID: ID of the beam
        :param i: Node i, a Node instance
        :param j: Node j, a Node instance 
        :param E: Module of elastivity
        :param I: Second moment of inertia of the cross-section
        :param A: cross-sectional area
        """
        self.ID = ID
        self.E = E
        self.i = i
        self.j = j
        self.A = A
        self.I = I

    @property
    def l(self):
        return math.sqrt((self.i.x-self.j.x)**2 + (self.i.y-self.j.y)**2)

    @property
    def EA(self):
        return self.E * self.A

    @property
    def EI(self):
        return self.E * self.I

    @property
    def ul(self):
        # upper left block of the stiffness matrix
        return np.matrix([
            np.array([self.EA / self.l,     0,                                  0]),
            np.array([0,                    12 * self.EI / (self.l ** 3),       6 * self.EI / (self.l ** 2)]),
            np.array([0,                    6 * self.EI / (self.l ** 2),        4 * self.EI / self.l]),
        ])

    @property
    def ll(self):
        # lower left block of the stiffness matrix
        return np.matrix([
            np.array([-self.EA / self.l,    0,                                  0]),
            np.array([0,                    -12 * self.EI / (self.l ** 3),      -6 * self.EI / (self.l ** 2)]),
            np.array([0,                    6 * self.EI / (self.l ** 2),        2 * self.EI / self.l]),
        ])

    @property
    def lr(self):
        # lower right block of the stiffness matrix
        return np.matrix([
            np.array([self.EA / self.l,     0,                                  0]),
            np.array([0,                    12 * self.EI / (self.l ** 3),       -6 * self.EI / (self.l ** 2)]),
            np.array([0,                    -6 * self.EI / (self.l ** 2),       4 * self.EI / self.l]),
        ])

    @property
    def ur(self):
        # upper right block of the stiffness matrix
        return np.matrix([
            np.array([-self.EA / self.l,    0,                                  0]),
            np.array([0,                    -12 * self.EI / (self.l ** 3),      6 * self.EI / (self.l ** 2)]),
            np.array([0,                    -6 * self.EI / (self.l ** 2),        2 * self.EI / self.l]),
        ])

    @property
    def K(self):
        return np.bmat([[self.ul, self.ur], [self.ll, self.lr]])


    def remove_DOF(self, node=None, ux=False, uy=False, rotz=False):

        print(np.transpose(np.nonzero(k)))  # nonzero elements
        print(np.transpose(np.where(k == 0)))  # zero elements

        pass


if __name__ == '__main__':

    beam = HermitianBeam2D(I=1., A=1., i=1, j=2)
    k = beam.K
    for i in range(3):
        k = np.delete(k, 0, axis=0)
        k = np.delete(k, 0, axis=1)
    pp.pprint(k)

    print(beam.dof)

    try:
        np.linalg.cholesky(k)
    except np.linalg.linalg.LinAlgError:
        print('baszki a 2D nem pozitiv definit')
        print('szimmetrikus?', (k.transpose() == k).all())

    print(sorted(np.linalg.eigvals(k)))
    print(sorted(np.linalg.eigvalsh(k)))

    P = 1

    F = np.matrix([P, 0, 0]).T

    print(F)

    print(np.linalg.inv(k) * F)
    print(P * (beam.l) / (beam.E * beam.A))

    print('')
    F = np.matrix([0, P, 0]).T
    print(np.linalg.inv(k) * F)
    print('phi', P * (beam.l ** 2) / (2 * beam.E * beam.I))
    print('delta', P * (beam.l ** 3) / (3 * beam.E * beam.I))

    exit()



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

    def __init__(self, ID=None, i=(0, 0), j=(1., 0.), dof=6, E=EE, nu=0.3, Iy=None, Iz=None, A=None):
        self.dof = dof
        self.ID = ID
        self.E = E
        self.nu = nu
        self.i = i
        self.j = j
        self.A = A
        self.Iy = Iy
        self.Iz = Iz

    @property
    def l(self):
        return math.sqrt((self.i[0] - self.j[0]) ** 2 + (self.i[1] - self.j[1]) ** 2)

    @property
    def G(self):
        return GG

    @property
    def K(self):
        l = self.l
        E_x = self.E * self.A
        G_x = self.G * self.A
        EIz_1 = EIz_2 = self.E * self.Iz
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