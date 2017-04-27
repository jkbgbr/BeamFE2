from __future__ import division
import numpy as np
import math
import pprint as pp

EE = 1.
GG = 1.


class HermitianBeam2D(object):
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
    u3 : Node i, rotation about Z axis
    u4 - u6: same for Node j
    
    """

    def __init__(self, ID=None, i=(0, 0), j=(1., 0.), dof=3, E=EE, nu=0.3, Iz=None, A=None):
        self.dof = dof
        self.ID = ID
        self.E = E
        self.nu = nu
        self.i = i
        self.j = j
        self.A = A
        self.Iz = Iz

    @property
    def l(self):
        return math.sqrt((self.i[0]-self.j[0])**2 + (self.i[1]-self.j[1])**2)

    @property
    def K(self):
        l = self.l
        EA = self.E * self.A
        EI = self.E * self.Iz

        _ret = np.array([
            [EA/l,      0,              0,              -EA/l,      0,              0],
            [0,         12*EI/(l**3),   6*EI/(l**2),    0,          -12*EI/(l**3),  6*EI/(l**2)],
            [0,         6*EI/(l**2),    4*EI/l,         0,          -6*EI/(l**2),   2*EI/l],
            [-EA/l,     0,              0,              EA/l,       0,              0],
            [0,         -12*EI/(l**3),  -6*EI/(l**2),   0,          12*EI/(l**3),   -6*EI/(l**2)],
            [0,         6*EI/(l**2),    2*EI/l,         0,          -6*EI/(l**2),   4*EI/l],
        ])

        return _ret
        # return _ret.astype(int)


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
        return math.sqrt((self.i[0]-self.j[0])**2 + (self.i[1]-self.j[1])**2)

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
            [0, 6 * (EIz_1 + EIz_2) / (l ** 3), 0, 0, 0, 2 * (2 * EIz_1 + EIz_2) / (l ** 2), 0, -6 * (EIz_1 + EIz_2) / (l ** 3), 0, 0, 0, 2 * (EIz_1 + 2 * EIz_2) / (l ** 2)],
            [0, 0, 6 * (EIy_1 + EIy_2) / (l ** 3), 0, -2 * (2 * EIy_1 + EIy_2) / (l ** 2), 0, 0, 0, -6 * (EIy_1 + EIy_2) / (l ** 3), 0, -2 * (EIy_1 + 2 * EIy_2) / (l ** 2), 0],
            [0, 0, 0, G_x / l, 0, 0, 0, 0, 0, -G_x / l, 0, 0],
            [0, 0, -2 * (2 * EIy_1 + EIy_2) / (l ** 2), 0, (3 * EIy_1 + EIy_2) / l, 0, 0, 0, 2 * (2 * EIy_1 + EIy_2) / (l ** 2), 0, (EIy_1 + EIy_2) / l, 0],
            [0, 2 * (2 * EIz_1 + EIz_2) / (l ** 2), 0, 0, 0, (3 * EIz_1 + EIz_2) / l, 0, -2 * (2 * EIz_1 + EIz_2) / (l ** 2), 0, 0, 0, (EIz_1 + EIz_2) / l],
            [-E_x / l, 0, 0, 0, 0, 0, E_x / l, 0, 0, 0, 0, 0],
            [0, -6 * (EIz_1 + EIz_2) / (l ** 3), 0, 0, 0, -2 * (2 * EIz_1 + EIz_2) / (l ** 2), 0, 6 * (EIz_1 + EIz_2) / (l ** 3), 0, 0, 0, -2 * (EIz_1 + 2 * EIz_2) / (l ** 2)],
            [0, 0, -6 * (EIy_1 + EIy_2) / (l ** 3), 0, 2 * (2 * EIy_1 + EIy_2) / (l ** 2), 0, 0, 0, 6 * (EIy_1 + EIy_2) / (l ** 3), 0, 2 * (EIy_1 + 2 * EIy_2) / (l ** 2), 0],
            [0, 0, 0, -G_x / l, 0, 0, 0, 0, 0, G_x / l, 0, 0],
            [0, 0, -2 * (EIy_1 + 2 * EIy_2) / (l ** 2), 0, (EIy_1 + EIy_2) / l, 0, 0, 0, 2 * (EIy_1 + 2 * EIy_2) / (l ** 2), 0, (EIy_1 + 3 * EIy_2) / l, 0],
            [0, 2 * (EIz_1 + 2 * EIz_2) / (l ** 2), 0, 0, 0, (EIz_1 + EIz_2) / l, 0, -2 * (EIz_1 + 2 * EIz_2) / (l ** 2), 0, 0, 0, (EIz_1 + 3 * EIz_2) / l]
        ])
        return _ret.astype(int)

if __name__ == '__main__':

    beam2D = HermitianBeam2D(Iz=1., A=1.)
    Ke = beam2D.K
    for i in range(3):
        Ke = np.delete(Ke, 0, axis=0)
        Ke = np.delete(Ke, 0, axis=1)
    pp.pprint(Ke)

    try:
        np.linalg.cholesky(Ke)
    except np.linalg.linalg.LinAlgError:
        print('baszki a 2D nem pozitiv definit')
        print('szimmetrikus?', (Ke.transpose() == Ke).all())

    print(sorted(np.linalg.eigvals(Ke)))
    print(sorted(np.linalg.eigvalsh(Ke)))

    exit()

    beam3D = HermitianBeam3D(Iz=1., Iy=1., A=1.)
    Ke = beam3D.K

    pp.pprint(Ke)



    try:
        np.linalg.cholesky(Ke)
    except np.linalg.linalg.LinAlgError:
        print('baszki a 3D nem pozitiv definit')
        print('szimmetrikus?', (Ke.transpose() == Ke).all())

    print(sorted(np.linalg.eigvals(beam3D.K)))
    print(sorted(np.linalg.eigvalsh(beam3D.K)))
