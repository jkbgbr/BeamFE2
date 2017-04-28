from __future__ import division
import numpy as np
import math
import itertools
import copy


class Node(object):
    """
    Node objects to be used with the FE Model
    """

    def __init__(self, ID=None, coords=()):
        self.ID = ID
        self.coords = coords
        self.DOFs = {'ux': True, 'uy': True, 'uz': True}

    def __repr__(self):
        return 'Node(ID=%d, coords=(%.2f, %.2f)' % (self.ID, self.x, self.y)

    def set_coords(self, coords):
        self.coords = coords

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
    The stiffness matrix is the same as in Cook: finite Element Modeling for Stress Analysis (ISBN 0-471-10774-3), 2.3-9
    
    """
    dof = 3  # degrees of freedom

    def __init__(self, ID=None, i=None, j=None, E=None, I=None, A=None, rho=None):
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
        self.rho = rho

    def __repr__(self):
        return 'HermitianBeam2D(ID=%d, i=%r, j=%r, E=%d, I=%.2f, A=%.2f' \
               % (self.ID, self.i, self.j, self.E, self.I, self.A)

    @property
    def nodes(self):
        return self.i, self.j

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
    def M(self):
        l = self.l
        _ret = self.rho * self.A * self.l / 420 * np.matrix([
            [140, 0, 0, 70, 0, 0],
            [0, 156, 22*l, 0, 54, -13*l],
            [0, 22*l, 4*l**2, 0, 13*l, -3*l**2],
            [70, 0, 0, 140, 0, 0],
            [0, 54, 13*l, 0, 156, -22*l],
            [0, -13*l, -3*l**2, 0, -22*l, 4*l**2]])
        return _ret


    @property
    def ul(self):
        # upper left block of the stiffness matrix
        return np.matrix([
            np.array([self.EA / self.l,     0,                                  0]),
            np.array([0,                    12 * self.EI / (self.l ** 3),       6 * self.EI / (self.l ** 2)]),
            np.array([0,                    6 * self.EI / (self.l ** 2),        4 * self.EI / self.l]),
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
    def Ke(self):
        # full stiffness matrix of the beam element in the element coordinate system
        return np.bmat([[self.ul, self.ur], [self.ll, self.lr]])

    @property
    def direction(self):
        # the direction of the beam in the global coordinate system
        # this is a crucial point for the transfer from the local to the global systems

        _dy = self.j.y - self.i.y
        _dx = self.j.x - self.i.x

        # if _dy != 0:
        #     print(self.i)
        #     print(self.j)
        #     print(_dx)
        #     print(_dy)
        #     print(math.degrees(np.arctan2(_dx, _dy)))
        #     exit()

        _ret = np.arctan2(_dy, _dx)
        # if _ret < 0:
        #     _ret += 2 * math.pi

        return _ret

    @property
    def Kg(self):
        # full stiffness matrix of the Beam element in the global coordinate system
        return self.transfer_matrix * self.Ke * self.transfer_matrix.T

    @property
    def Kg_ul(self):
        # upper left block of the element stiffness matrix in the global coordinate system
        return self.Kg[:self.dof, :self.dof]

    @property
    def Kg_ur(self):
        # upper right block of the element stiffness matrix in the global coordinate system
        return self.Kg[:self.dof, self.dof:]

    @property
    def Kg_ll(self):
        # lower left block of the element stiffness matrix in the global coordinate system
        return self.Kg[self.dof:, :self.dof]

    @property
    def Kg_lr(self):
        # lower right block of the element stiffness matrix in the global coordinate system
        return self.Kg[self.dof:, self.dof:]

    @property
    def transfer_matrix(self):
        # matrix to rotate the stiffness matrix for compilation
        return transfer_matrix(alpha=self.direction, asdegree=False, dof=self.dof, blocks=2)

    # def remove_DOF(self, node=None, ux=False, uy=False, rotz=False):
    #     print(np.transpose(np.nonzero(k)))  # nonzero elements
    #     print(np.transpose(np.where(k == 0)))  # zero elements
    #     pass


class Structure(object):
    """
    The Structure, composed of the FE Beams.
    """
    def __init__(self, beams, BCs):
        self.BCs = BCs
        self.beams = beams
        self.dof = self.beams[0].dof
        self._load_vector = None

    @property
    def dofnames(self):
        assert self.dof in [3, 6]
        if self.dof == 3:
            return 'FX', 'FY', 'MZ'
        else:
            return 'FX', 'FY', 'FZ', 'MX', 'MY', 'MZ'

    @property
    def dofnumbers(self):
        assert self.dof in [3, 6]
        if self.dof == 3:
            return 0, 1, 2
        else:
            return 0, 1, 2, 3, 4, 5

    @property
    def nodes(self):
        # set of Nodes of the model
        return set(itertools.chain.from_iterable([x.nodes for x in self.beams]))

    @property
    def sumdof(self):
        # sum of DOFs, without eliminating for BCs
        return self.dof * len(self.nodes)

    @property
    def K(self):
        # the compiled stiffness matrix
        return self.compile_global_stiffness_matrix()

    @property
    def K_with_BC(self):
        # the stiffness matrix with the boundary conditions taken into account
        # CURRENTLY: NODE 1 is clamped
        glo = copy.deepcopy(self.K)
        for i in range(self.dof):
            glo = np.delete(glo, 0, axis=0)
            glo = np.delete(glo, 0, axis=1)
        return glo

    @property
    def q(self):
        # the load vector
        return self._load_vector.T

    @property
    def q_with_BC(self):
        # the load vector, changed to be usable with the BCs used in K_with_BC
        _q = copy.deepcopy(self.q)
        for i in range(self.dof):
            _q = np.delete(_q, 0, axis=0)
        return _q

    def add_single_dynam_to_node(self, nodeID=None, dynam=None, clear=False):
        """
        Adds a dynam (FX, FY, MZ) to the chosen Node of the model.
        Checks if the node exists.
        Checks if the name of the dynam is OK.
        clears previous loads on the node, if clear is True
        
        :param nodeID: 
        :param dynam: 
        :param clear: 
        :return: 
        """
        #

        assert nodeID in [x.ID for x in self.nodes]
        assert all([x in self.dofnames for x in dynam.keys()])

        if self._load_vector is None:
            self._load_vector = np.matrix(np.zeros(self.sumdof))

        for k, v in dynam.items():

            for name, number in zip(self.dofnames, self.dofnumbers):  # pairs e.g. 'FX' with 0, 'FY' with 1 etc.
                _sti = nodeID * self.dof + number  # starting index
                if k == name:
                    if clear:
                        self._load_vector[0, _sti: _sti + 1] = 0
                    self._load_vector[0, _sti: _sti + 1] += v

        # self._load_vector *= self.transfer_matrix

    @property
    def stiffness_matrix_is_symmetric(self):
        # checks if the global stiffness matrix (without BCs) is symmetric. MUST be.
        diff = self.K.transpose() - self.K
        if np.allclose(diff, np.zeros(self.sumdof)):
            return True
        else:
            return False

    @property
    def stiffness_matrix_is_nonsingular(self):
        # checks if the global stiffness matrix (with BCs) is positive definite. MUST be.
        try:
            np.linalg.cholesky(self.K_with_BC)
            return True
        except np.linalg.linalg.LinAlgError:
            return False

    def solve(self):
        """
        solves the system, returns the vector of displacements.
        :return: 
        """
        assert self.stiffness_matrix_is_nonsingular
        try:
            return np.linalg.inv(self.K_with_BC) * self.q_with_BC
        except np.linalg.linalg.LinAlgError:
            raise Exception('OK something went wrong')

    def compile_global_stiffness_matrix(self):
        """
        Compiles the global stiffness matrix from the element matrices.
        :return: 
        """
        _sumdof = self.sumdof
        _empty = np.zeros(_sumdof ** 2)
        _empty = np.matrix(_empty.reshape(_sumdof, _sumdof))
        for b in self.beams:
            _sti = (b.i.ID - 1) * self.dof  # starting element of the block for node i
            _eni = _sti + self.dof  # end element for the block if node i
            _stj = (b.j.ID - 1) * self.dof
            _enj = _stj + self.dof
            # upper left block
            _empty[_sti:_eni, _sti:_eni] += b.Kg_ul
            # left lower block
            _empty[_stj:_enj, _sti:_eni] += b.Kg_ll
            # upper right block
            _empty[_sti:_eni, _stj:_enj] += b.Kg_ur
            # lower right block
            _empty[_stj:_enj, _stj:_enj] += b.Kg_lr
        return _empty


def transfer_matrix(alpha, asdegree=False, blocks=2, dof=3):
    # matrix to rotate the stiffness matrix for compilation
    if asdegree:
        alpha = math.radians(alpha)
    cs = math.cos(alpha)
    ss = math.sin(alpha)
    _block = np.matrix([[cs,    -ss,    0],
                        [ss,    cs,     0],
                        [0,     0,      1]])
    _sumdof = blocks * dof
    _empty = np.zeros(_sumdof ** 2)
    _empty = np.matrix(_empty.reshape(_sumdof, _sumdof))

    for b in range(blocks):
        _sti = b * dof  # starting element of the block for node i
        _eni = _sti + dof  # end element for the block if node i
        _empty[_sti:_eni, _sti:_eni] += _block

    return _empty


if __name__ == '__main__':

    n1 = Node(ID=1, coords=(0, 0))
    n2 = Node(ID=2, coords=(100, 0))
    n3 = Node(ID=3, coords=(200, 0))
    n4 = Node(ID=4, coords=(300, 0))
    b1 = HermitianBeam2D(E=21000., ID=1, I=833.33, A=100., i=n1, j=n2, rho=7850)
    b2 = HermitianBeam2D(E=21000., ID=2, I=833.33, A=100., i=n2, j=n3, rho=7850)
    b3 = HermitianBeam2D(E=21000., ID=3, I=833.33, A=100., i=n3, j=n4, rho=7850)

    structure = Structure(beams=[b1, b2, b3], BCs=None)

    print(b1.M)

    exit()

    # structure.add_single_dynam_to_node(nodeID=3, dynam={'FX': 1})
    # structure.add_single_dynam_to_node(nodeID=2, dynam={'FY': -1})
    # structure.add_single_dynam_to_node(nodeID=2, dynam={'MZ': 1000000})
    structure.add_single_dynam_to_node(nodeID=2, dynam={'FX': 1000, 'FY': 1000, 'MZ': 1000000})
    # structure.add_single_dynam_to_node(nodeID=2, dynam={'FX': -1000, 'FY': -1000, 'MZ': -1000000})

    disps = structure.solve()
    print(disps)


# class HermitianBeam3D(object):
#     """
#     A hermitian 3D FE Beam with 2 nodes and 6 DOFs each node: 3 translational, 3 rotational.
#     A 2D version is also available with 3 DOFs each node: 2 translational, 1 rotational.
#     Small deformations are assumed.
#     Bernoulli-Euler beam theory (no shear deformations)
#     Non-restrained torsion only (no warping)
#
#     Nodes: i, j with direction i -> j
#     Element coordinate system is right-handed, thumb: X, pointing finger Y, middle finger Z
#     Displacements are positive when going away from the origin.
#     Rotations are positive when clockwise when looking away from the origin.
#
#     Displacement vector u:
#     u1 : Node i, displacement in X direction
#     u2 : Node i, displacement in Y direction
#     u3 : Node i, displacement in Z direction
#     u4 : Node i, rotation about X axis
#     u5 : Node i, rotation about Y axis
#     u6 : Node i, rotation about Z axis
#     u7 - u12: same for Node j
#
#     """
#
#     def __init__(self, ID=None, i=(0, 0), j=(1., 0.), dof=6, E=EE, nu=0.3, Iy=None, Iz=None, A=None):
#         self.dof = dof
#         self.ID = ID
#         self.E = E
#         self.nu = nu
#         self.i = i
#         self.j = j
#         self.A = A
#         self.Iy = Iy
#         self.Iz = Iz
#
#     @property
#     def l(self):
#         return math.sqrt((self.i[0] - self.j[0]) ** 2 + (self.i[1] - self.j[1]) ** 2)
#
#     @property
#     def G(self):
#         return GG
#
#     @property
#     def Ke(self):
#         l = self.l
#         E_x = self.E * self.A
#         G_x = self.G * self.A
#         EIz_1 = EIz_2 = self.E * self.Iz
#         EIy_1 = EIy_2 = self.E * self.Iy
#
#         _ret = np.array([
#             [E_x / l, 0, 0, 0, 0, 0, -E_x / l, 0, 0, 0, 0, 0],
#             [0, 6 * (EIz_1 + EIz_2) / (l ** 3), 0, 0, 0, 2 * (2 * EIz_1 + EIz_2) / (l ** 2), 0,
#              -6 * (EIz_1 + EIz_2) / (l ** 3), 0, 0, 0, 2 * (EIz_1 + 2 * EIz_2) / (l ** 2)],
#             [0, 0, 6 * (EIy_1 + EIy_2) / (l ** 3), 0, -2 * (2 * EIy_1 + EIy_2) / (l ** 2), 0, 0, 0,
#              -6 * (EIy_1 + EIy_2) / (l ** 3), 0, -2 * (EIy_1 + 2 * EIy_2) / (l ** 2), 0],
#             [0, 0, 0, G_x / l, 0, 0, 0, 0, 0, -G_x / l, 0, 0],
#             [0, 0, -2 * (2 * EIy_1 + EIy_2) / (l ** 2), 0, (3 * EIy_1 + EIy_2) / l, 0, 0, 0,
#              2 * (2 * EIy_1 + EIy_2) / (l ** 2), 0, (EIy_1 + EIy_2) / l, 0],
#             [0, 2 * (2 * EIz_1 + EIz_2) / (l ** 2), 0, 0, 0, (3 * EIz_1 + EIz_2) / l, 0,
#              -2 * (2 * EIz_1 + EIz_2) / (l ** 2), 0, 0, 0, (EIz_1 + EIz_2) / l],
#             [-E_x / l, 0, 0, 0, 0, 0, E_x / l, 0, 0, 0, 0, 0],
#             [0, -6 * (EIz_1 + EIz_2) / (l ** 3), 0, 0, 0, -2 * (2 * EIz_1 + EIz_2) / (l ** 2), 0,
#              6 * (EIz_1 + EIz_2) / (l ** 3), 0, 0, 0, -2 * (EIz_1 + 2 * EIz_2) / (l ** 2)],
#             [0, 0, -6 * (EIy_1 + EIy_2) / (l ** 3), 0, 2 * (2 * EIy_1 + EIy_2) / (l ** 2), 0, 0, 0,
#              6 * (EIy_1 + EIy_2) / (l ** 3), 0, 2 * (EIy_1 + 2 * EIy_2) / (l ** 2), 0],
#             [0, 0, 0, -G_x / l, 0, 0, 0, 0, 0, G_x / l, 0, 0],
#             [0, 0, -2 * (EIy_1 + 2 * EIy_2) / (l ** 2), 0, (EIy_1 + EIy_2) / l, 0, 0, 0,
#              2 * (EIy_1 + 2 * EIy_2) / (l ** 2), 0, (EIy_1 + 3 * EIy_2) / l, 0],
#             [0, 2 * (EIz_1 + 2 * EIz_2) / (l ** 2), 0, 0, 0, (EIz_1 + EIz_2) / l, 0,
#              -2 * (EIz_1 + 2 * EIz_2) / (l ** 2), 0, 0, 0, (EIz_1 + 3 * EIz_2) / l]
#         ])
#         return _ret.astype(int)
