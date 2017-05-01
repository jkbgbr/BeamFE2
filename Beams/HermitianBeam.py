from __future__ import division
import numpy as np
# print precision, no scientific notation, wide lines
np.set_printoptions(precision=6, suppress=False, linewidth=140)
import math
import itertools
import copy
from drawing import draw_beam


# http://12000.org/my_notes/stiffness_matrix/stiffness_matrix_report.htm


class Node(object):
    """
    Node objects to be used with the FE Model
    """

    def __init__(self, ID=None, coords=()):
        self.ID = ID
        self.coords = coords
        # self.DOFs = {'i': {'ux': True, 'uy': True, 'rotz': True},
        #              'j': {'ux': True, 'uy': True, 'rotz': True}}

    def __repr__(self):
        return 'Node(ID=%d, coords=(%.2f, %.2f)' % (self.ID, self.x, self.y)

    @classmethod
    def from_dict(cls, adict):
        try:
            node = Node(
                ID=adict['ID'],
                coords=adict['coords'],
            )
            assert all([hasattr(node, x) for x in adict.keys()])  # check if the dict is too long
        except KeyError as e:  # check if the dict is too short
            raise Exception('Missing key from dict when creating %s: %s' % (cls.__name__, e))
        return node

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
    Node numbering MUST start with 0 and MUST be sequential, without missing values. Renumber Ndes if necessary.
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
        self.displacements = None

    def __repr__(self):
        return 'HermitianBeam2D(ID=%d, i=%r, j=%r, E=%d, I=%.2f, A=%.2f' \
               % (self.ID, self.i, self.j, self.E, self.I, self.A)

    @classmethod
    def from_dict(cls, adict):
        try:
            beam = HermitianBeam2D(
                ID=adict['ID'],
                i=adict['i'],
                j=adict['j'],
                E=adict['E'],
                I=adict['I'],
                A=adict['A'],
                rho=adict['rho'],
            )
            assert all([hasattr(beam, x) for x in adict.keys()])  # check if the dict is too long
        except KeyError as e:  # check if the dict is too short
            raise Exception('Missing key from dict when creating %s: %s' % (cls.__name__, e))

        return beam

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

    def _Me(self):
        l = self.l
        _ret = self.rho * self.A * self.l / 420 * \
               np.matrix([
                   [140, 0, 0, 70, 0, 0],
                   [0, 156, 22*l, 0, 54, -13*l],
                   [0, 22*l, 4*l**2, 0, 13*l, -3*l**2],
                   [70, 0, 0, 140, 0, 0],
                   [0, 54, 13*l, 0, 156, -22*l],
                   [0, -13*l, -3*l**2, 0, -22*l, 4*l**2]])
        return _ret

    @property
    def Me(self):
        return self._Me()

    def _Ke(self):
        # full stiffness matrix of the beam element in the element coordinate system
        return np.matrix([
            np.array([self.EA / self.l, 0, 0, -self.EA / self.l, 0, 0]),
            np.array([0, 12 * self.EI / (self.l ** 3), 6 * self.EI / (self.l ** 2), 0, -12 * self.EI / (self.l ** 3), 6 * self.EI / (self.l ** 2)]),
            np.array([0, 6 * self.EI / (self.l ** 2), 4 * self.EI / self.l, 0, -6 * self.EI / (self.l ** 2), 2 * self.EI / self.l]),
            np.array([-self.EA / self.l, 0, 0, self.EA / self.l,  0, 0]),
            np.array([0, -12 * self.EI / (self.l ** 3), -6 * self.EI / (self.l ** 2), 0, 12 * self.EI / (self.l ** 3), -6 * self.EI / (self.l ** 2)]),
            np.array([0, 6 * self.EI / (self.l ** 2), 2 * self.EI / self.l, 0, -6 * self.EI / (self.l ** 2), 4 * self.EI / self.l]),
            ])

    @property
    def Ke(self):
        return self._Ke()

    @property
    def direction(self):
        # the direction of the beam in the global coordinate system
        # this is a crucial point for the transfer from the local to the global systems
        _dy = self.j.y - self.i.y
        _dx = self.j.x - self.i.x
        return np.arctan2(_dy, _dx)

    @property
    def Kg(self):
        # full stiffness matrix of the Beam element in the global coordinate system
        return self.transfer_matrix * self.Ke * self.transfer_matrix.T

    @property
    def Mg(self):
        # full stiffness matrix of the Beam element in the global coordinate system
        return self.transfer_matrix * self.Me * self.transfer_matrix.T

    @property
    def transfer_matrix(self):
        # matrix to rotate the stiffness matrix for compilation
        return transfer_matrix(alpha=self.direction, asdegree=False, dof=self.dof, blocks=2)


class Structure(object):
    """
    The Structure, composed of the FE Beams.
    """
    def __init__(self, beams=None, supports=None):
        self.beams = beams
        self._load_vector = None
        self.supports = supports
    #
    # @classmethod
    # def from_dict(cls, adict):
    #
    #     _nodes = []
    #     for added, values in adict.items():
    #         if added == 'node':
    #             _nodes.append(Node.from_dict(values))
    #         if added == 'beams':
    #             HermitianBeam2D.from_dict(values)
    #         if added == 'BC':
    #

    @property
    def dof(self):
        return self.beams[0].dof

    @property
    def loadnames(self):
        assert self.dof in [3, 6]
        if self.dof == 3:
            return 'FX', 'FY', 'MZ'
        else:
            return 'FX', 'FY', 'FZ', 'MX', 'MY', 'MZ'

    @property
    def dofnames(self):
        assert self.dof in [3, 6]
        if self.dof == 3:
            return 'ux', 'uy', 'rotz'
        else:
            return 'ux', 'uy', 'uz', 'rotx', 'roty', 'rotz'

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
        # copy of the the stiffness matrix with the boundary conditions taken into account
        _K = copy.deepcopy(self.K)
        for k, v in self.supports.items():
            for dof in v:
                _pos = self.position_in_matrix(nodeID=k, DOF=dof)
                _K[_pos, _pos] += 10e20
        return _K

    @property
    def q(self):
        # the load vector
        return self._load_vector.T

    @property
    def q_with_BC(self):
        return self.q

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
        assert all([x in self.loadnames for x in dynam.keys()])

        if self._load_vector is None:
            self._load_vector = np.matrix(np.zeros(self.sumdof))

        # clear all loads defined previously
        if clear:
            # print('loads cleared')
            self._load_vector[0, :-1] = 0

        for k, v in dynam.items():
            for name, number in zip(self.loadnames, range(self.dof)):  # pairs e.g. 'FX' with 0, 'FY' with 1 etc.
                # _sti = nodeID * self.dof + number  # starting index
                if k == name:
                    _sti = self.position_in_matrix(nodeID=nodeID, dynam=k)
                    self._load_vector[0, _sti] += v
                    print('added: Node %d, dynam %s = %.2f' % (nodeID, k, v))
        # print(self._load_vector)
        # print(self.sumdof)
        # self._load_vector *= self.transfer_matrix

    @property
    def stiffness_matrix_is_symmetric(self):
        # checks if the global stiffness matrix (without BCs) is symmetric. MUST be.
        diff = self.K.transpose() - self.K
        if np.allclose(diff, np.zeros(self.sumdof)):
            return True
        else:
            print('The stiffness matrix is not symmetric')
            return False

    @property
    def stiffness_matrix_is_nonsingular(self):
        # checks if the global stiffness matrix (with BCs) is positive definite. MUST be.
        try:
            np.linalg.cholesky(self.K_with_BC)
            return True
        except np.linalg.linalg.LinAlgError:
            return False

    @property
    def node_numbering_is_ok(self):
        # checks node numbering: they must be sequential, without missing values
        if set([x.ID for x in self.nodes]) == set(range(1, len(self.nodes)+1)):
            return True
        else:
            print('Node numbering is not ok')
            return False

    @property
    def stiffness_matrix_is_ok(self):
        if self.stiffness_matrix_is_nonsingular:
            return True
        else:
            print('The stiffness matrix is singular')
            return False

    def solve(self):
        """
        solves the system, returns the vector of displacements.
        :return: 
        """
        assert self.stiffness_matrix_is_ok
        assert self.node_numbering_is_ok
        self.displacements = np.linalg.inv(self.K_with_BC) * self.q_with_BC
        return True

    def compile_global_stiffness_matrix(self):
        """
        Compiles the global stiffness matrix from the element matrices.
        :return: 
        """
        return compile_global_matrix(self.beams, stiffness=True, mass=False)

    def position_in_matrix(self, nodeID=None, DOF=None, dynam=None):
        """
        Tells the index of the given nodeID, DOF in a global K or M matrix. 
        :param nodeID: 
        :param DOF: 
        :return: 
        """
        # assert nodeID in [x.nodeID for x in self.nodes]
        if DOF is not None:
            assert DOF in self.dofnames
            return (nodeID - 1) * self.dof + self.dofnames.index(DOF)

        if dynam is not None:
            assert dynam in self.loadnames
            return (nodeID - 1) * self.dof + self.loadnames.index(dynam)
        # print('Adding support to Node %d, DOF %s at index %d' % (nodeID, DOF, _ret))
        # return _ret


def compile_global_matrix(beams, stiffness=True, mass=False):
    """
    Compiles the global stiffness matrix from the element matrices.
    :return: 
    """
    assert stiffness is not mass  # either or

    nodes = set(itertools.chain.from_iterable([x.nodes for x in beams]))  # nodes of the beams
    dof = beams[0].dof  # nr. of dofs
    _sumdof = len(nodes) * dof
    _empty = np.zeros(_sumdof ** 2)
    _empty = np.matrix(_empty.reshape(_sumdof, _sumdof))
    for b in beams:
        if stiffness:
            mtrx = b.Kg
        else:
            mtrx = b.Mg
        _sti = (b.i.ID - 1) * dof  # starting element of the block for node i
        _eni = _sti + dof  # end element for the block if node i
        _stj = (b.j.ID - 1) * dof
        _enj = _stj + dof
        # upper left block
        _empty[_sti:_eni, _sti:_eni] += mtrx[:dof, :dof]
        # left lower block
        _empty[_stj:_enj, _sti:_eni] += mtrx[dof:, :dof]
        # upper right block
        _empty[_sti:_eni, _stj:_enj] += mtrx[:dof, dof:]
        # lower right block
        _empty[_stj:_enj, _stj:_enj] += mtrx[dof:, dof:]
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

#


if __name__ == '__main__':

    # structure = Structure()

    # nodes
    n1 = Node.from_dict(adict={'ID': 1, 'coords': (0, 0)})
    n2 = Node.from_dict(adict={'ID': 2, 'coords': (250, 0)})
    n3 = Node(ID=3, coords=(400, 0))

    # beams
    b1 = HermitianBeam2D.from_dict(adict={'ID': 1, 'E': 210000., 'I': 39760.78, 'A': 706.5, 'rho': 7.85e-5, 'i': n1, 'j': n2})
    b2 = HermitianBeam2D(E=21000., ID=2, I=39760.78, A=706.5, i=n2, j=n3, rho=7.85e-5)

    # supports
    BCs = {1: ['ux', 'uy', 'rotz']}  # supports as dict

    # this is the structure
    structure = Structure(beams=[b1, b2], supports=BCs)

    print(structure)

    print(structure.K_with_BC)

    # adding loads
    structure.add_single_dynam_to_node(nodeID=3, dynam={'FX': 10000000}, clear=True)

    print(structure.q_with_BC)

    # solver :-)
    structure.solve()
    disps = structure.displacements



    # sorry, no postprocessng, however you can print the displacements
    print(disps)

    draw_beam.draw_structure(structure)


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
