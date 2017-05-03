from __future__ import division
import numpy as np
from scipy.linalg import eigh
# print precision, no scientific notation, wide lines
np.set_printoptions(precision=6, suppress=True, linewidth=250)
import math
import itertools
import copy
from drawing import draw_beam
from Beams import Sections as sections
import pprint as pp
from drawing import _plotting_available, plt


# todo:
# Lin stat:
#   adding loads off-node (concentrated force and moment, bending)
#   http://12000.org/my_notes/stiffness_matrix/stiffness_matrix_report.htm for internal actions, at the end of the page
#   calculating stresses (Sections need stress points for this...)
#   defining hinges?
#   3D Beam element?
# Modal:
#   tests, ideally opensess-based and with analytical formulae
#   (http://opensees.berkeley.edu/wiki/index.php/Eigen_analysis_of_a_two-story_shear_frame)
#   lumped mass matrix?
#   plot mode shapes (as a matter of fact, plotting should be re-thought)
#   modal masses, modal participation
#   EC-based stuff: modal antwort, summation...
#   ability to add mass
#   ability to convert loads to masses
#   ability to choose which directions mass act?
# Buckling:
#   implement. the geometrical stiffness matrix is needed for that.
# Displacements:
#   displacements are to be stored in a way that multiple shapes can be handled (linstat: 1 set, all other: multiple sets)


class Node(object):
    """
    Node objects to be used with the FE Model
    """

    def __init__(self, ID=None, coords=()):
        self.ID = ID
        self.coords = coords

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
    dofnames = 'ux', 'uy', 'rotz'
    number_internal_points = 20

    def __init__(self, ID=None, i=None, j=None, E=None, crosssection=None, I=None, A=None, rho=None):
        """
        :type crosssection: object
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
        if crosssection is not None:
            self.A = crosssection.A
            self.I = crosssection.I['x']
        else:
            self.A = A
            self.I = I
        self._end_DOFs = {'i': list(copy.deepcopy(self.dofnames)), 'j': list(copy.deepcopy(self.dofnames))}

        self.rho = rho
        self.displacements = None  # displacements in the GLOBAL system

    def __repr__(self):
        return 'HermitianBeam2D(ID=%d, i=%r, j=%r, E=%d, I=%.2f, A=%.2f' \
               % (self.ID, self.i, self.j, self.E, self.I, self.A)

    # shape functions. the input argument x is the "running parameter", a value normalized with the beam length
    def N1(self, x, L=1.):
        # shape function for axial displacements, node 'i'
        x *= L
        return x / L

    def N2(self, x, L=1.):
        # shape function for shear displacements, node 'i'
        x *= L
        return 1 - (3 * (x ** 2) / (L ** 2)) + (2 * (x ** 3) / (L ** 3))

    def N3(self, x, L=1.):
        # shape function for bending, node 'i'
        x *= L
        return x - (2 * (x ** 2) / L) + (x ** 3) / (L ** 2)

    def N4(self, x, L=1.):
        # shape function for axial displacements, node 'j'
        x *= L
        return 1 - x / L
        # return -1 * self.N1(x, L=L) + 1

    def N5(self, x, L=1.):
        # shape function for shear displacements, node 'j'
        x *= L
        return (3 * (x ** 2) / (L ** 2)) - (2 * (x ** 3) / (L ** 3))
        # return -1 * self.N2(x, L=L) + 1

    def N6(self, x, L=1.):
        # shape function for bending, node 'j'
        x *= L
        return - (x ** 2) / L + (x ** 3) / (L ** 2)

    def N(self, x, L=1.):
        # the N matrix, assembled. 2 x 6 matrix, multiplied with self.local_displacements yields the deflections
        # for the internal points as a tuple with ux, uy values in the local system.
        return np.matrix([[self.N1(x=x, L=L), 0,                    0,
                           self.N4(x=x, L=L),  0,                  0],
                          [0,                 self.N2(x=x, L=L),    self.N3(x=x, L=L),
                           0,                  self.N5(x=x, L=L),  self.N6(x=x, L=L)]])

    def deflected_shape(self, local=True, scale=1.):
        """
        The deflected shape of the beam. based on the local_displacements, results are in the local coordinate system.
        (ux, uy) = N * q
        if local = False, the values are transferred in the global coordinate system
        :return: 
        """
        _ip = self.number_internal_points
        _deflected_shape = []

        _t = transfer_matrix(-self.direction, asdegree=False, blocks=1, dof=2)  # 2x2 transform matrix

        for i in range(_ip + 1):
            _val = np_matrix_tolist(self.N(x=i / _ip, L=self.l) * self.local_displacements)  # displacements in the local system
            _val = [x * scale for x in _val]
            if not local:
                _val = ((i / _ip) * self.l + _val[0], _val[1])  # the deflected shape in the local coordinate system
                _val *= _t  # the deflected shape rotated
                _val = np_matrix_tolist(_val)
                _val = (_val[0] + self.i.x, _val[1] + self.i.y)
            _deflected_shape.append(_val)

        return _deflected_shape

    def nodal_reactions(self):
        """ reactions in the local coordinate system of the beam """
        return self.Ke * self.local_displacements

    @property
    def local_displacements(self):
        """
        calculates the nodal displacements in the elements local system.
        :return: 
        """
        # de-rotating matrix to transfer the displacements from the global to the elements local system
        T = transfer_matrix(-self.direction, asdegree=False, blocks=len(self.nodes))
        return T * self.displacements

    def displacement_component(self, component=None, localsystem=False):
        """
        Returns a list with the displacements of the component defined
        These are the nodal displacements calculated 
        :param localsystem: boolean telling if the results to be provided are in the global or in the local system
        :param component: any value from the dofnames list
        :return: 
        """
        assert component in self.dofnames
        if localsystem:
            _ret = self.local_displacements[self.dofnames.index(component)::self.dof]
        else:
            _ret = self.displacements[self.dofnames.index(component)::self.dof]
        return np_matrix_tolist(_ret)
        # return list(itertools.chain.from_iterable(_ret.tolist()))

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

    # def _Me(self):
    #     L = self.l
    #     _ret = (self.rho / 10.) * self.A * L * \
    #            np.matrix([
    #                [1/2., 0, 0, 0, 0, 0],
    #                [0, -0*(L**2)/24, 0, 0, 0, 0],
    #                [0, 0, -0*(L**2)/24, 0, 0, 0],
    #                [0, 0, 0, 1/2., 0, 0],
    #                [0, 0, 0, 0, -0*(L**2)/24, 0],
    #                [0, 0, 0, 0, 0, -0*(L**2)/24]])
    #     return _ret

    def _Me(self):
        L = self.l
        _ret = self.rho * (self.A / 1e6) * (self.l / 1e3) / 420 * \
               np.matrix([
                   [140, 0, 0, 70, 0, 0],
                   [0, 156, 22*L, 0, 54, -13*L],
                   [0, 22*L, 4*(L**2), 0, 13*L, -3*(L**2)],
                   [70, 0, 0, 140, 0, 0],
                   [0, 54, 13*L, 0, 156, -22*L],
                   [0, -13*L, -3*(L**2), 0, -22*L, 4*(L**2)]])
        return _ret

    @property
    def Me(self):
        return self._Me()

    def _Ke_geom(self):
        # the geometrical stiffness matrix, from H-P. Gavin CEE 421L. Matrix Structural Anyalsis - Duke University
        _locdisp = self.local_displacements
        T = _locdisp[3] - _locdisp[0]  # change of the axial length
        L = self.l
        _ret = np.matrix([
            [0, 0, 0, 0, 0, 0],
            [0, 6./5., L/10., 0, -(6./5.), L/10.],
            [0, L/10., 2*(L**2)/15., 0, -(L/10.), -(L/30.)],
            [0, 0, 0, 0, 0, 0],
            [0, -(6./5.), -(L/10.), 0, (6./5.), -(L/10.)],
            [0, (L/10.), -(L**2)/30., 0, -(L/10.), -(2*(L**2))/15.]
        ])
        _ret = np.multiply((T / L), _ret)

        return _ret

    @property
    def Ke_geom(self):
        return self._Ke_geom()

    def _Ke(self):
        # full stiffness matrix of the beam element in the element coordinate system
        _ret = np.matrix([
            np.array([self.EA / self.l, 0, 0, -self.EA / self.l, 0, 0]),
            np.array([0, 12 * self.EI / (self.l ** 3), 6 * self.EI / (self.l ** 2), 0, -12 * self.EI / (self.l ** 3), 6 * self.EI / (self.l ** 2)]),
            np.array([0, 6 * self.EI / (self.l ** 2), 4 * self.EI / self.l, 0, -6 * self.EI / (self.l ** 2), 2 * self.EI / self.l]),
            np.array([-self.EA / self.l, 0, 0, self.EA / self.l,  0, 0]),
            np.array([0, -12 * self.EI / (self.l ** 3), -6 * self.EI / (self.l ** 2), 0, 12 * self.EI / (self.l ** 3), -6 * self.EI / (self.l ** 2)]),
            np.array([0, 6 * self.EI / (self.l ** 2), 2 * self.EI / self.l, 0, -6 * self.EI / (self.l ** 2), 4 * self.EI / self.l]),
            ])

        return _ret

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

    def matrix_in_global(self, mtrx=None):
        # matrix of the Beam element in the global coordinate system
        return self.transfer_matrix * mtrx * self.transfer_matrix.T

    # @property
    # def Kg(self):
    #     # full stiffness matrix of the Beam element in the global coordinate system
    #     return self.transfer_matrix * self.Ke * self.transfer_matrix.T
    #
    # @property
    # def Mg(self):
    #     # full stiffness matrix of the Beam element in the global coordinate system
    #     return self.transfer_matrix * self.Me * self.transfer_matrix.T

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
        self.displacements = None  # global dispacements from the linear elastic analysis
        self.criticals = None  # critical forces from the buckling analysis
        self.buckling_shapes = None  # buckling shapes
        self.frequencies = None  # eigenfrequencies
        self.nodal_shapes = None  # modal shapes

    def displacements_as_dict(self, local=False):

        uxs = self.displacement_component(component='ux')
        uys = self.displacement_component(component='uy')
        rotzs = self.displacement_component(component='rotz')

        _ret = {}

        for beam in self.beams:
            _ret[beam.ID] = {'Node i': {'ux': uxs[beam.i.ID-1], 'uy': uys[beam.i.ID-1], 'rotz': rotzs[beam.i.ID-1]},
                             'Node j': {'ux': uxs[beam.j.ID-1], 'uy': uys[beam.j.ID-1], 'rotz': rotzs[beam.j.ID-1]}
                             }

        return _ret

    def draw(self, show=True):
        draw_beam.draw_structure(self, show=show, displacements=True)

    def displacements_for_beams(self):
        """
        De-compiles the global displacement vector and assigns the results to the beams
        :return: True, if succesful
        """
        uxs = self.displacement_component(component='ux')
        uys = self.displacement_component(component='uy')
        rotzs = self.displacement_component(component='rotz')

        for beam in self.beams:
            # displacements in the global system
            beam.displacements = np.matrix([uxs[beam.i.ID-1], uys[beam.i.ID-1], rotzs[beam.i.ID-1],
                                            uxs[beam.j.ID-1], uys[beam.j.ID-1], rotzs[beam.j.ID-1]]).T

        return True

    def displacement_component(self, component=None):
        """
        Returns a list with the displacements of the component defined
        These are the nodal displacements from the structure
        :param component: any value from the dofnames list
        :return: 
        """
        assert component in self.dofnames
        _ret = self.displacements[self.dofnames.index(component)::self.dof]
        return np_matrix_tolist(_ret)
        # return list(itertools.chain.from_iterable(_ret.tolist()))

    @property
    def resulting_displacement(self):
        _ux = self.displacement_component(component='ux')
        _uy = self.displacement_component(component='uy')
        return [math.sqrt(x ** 2 + y ** 2) for x, y in zip(_ux, _uy)]

    def displacements_of_beam(self, beam):
        pass

    def node_by_ID(self, id=None):
        # return the Node object that has the ID
        _ret = [x for x in self.nodes if x.ID == id]
        assert len(_ret) == 1
        return _ret[0]

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

    def zero_BC(self, mtrx=None):
        """
        Eliminates the rows and columns of the BC
        """
        mtrx = copy.deepcopy(mtrx)
        print(mtrx.size)
        _to_eliminate = []  # list of rows and columns to eliminate
        for nodeID, DOFs in self.supports.items():
            print(nodeID)
            for DOF in DOFs:
                print(DOF)
                _to_eliminate.append(self.position_in_matrix(nodeID=nodeID, DOF=DOF))

        _to_eliminate.sort()

        for _ in _to_eliminate[::-1]:
            mtrx[_] = 0
            mtrx[:, _] = 0

        return mtrx

    def eliminate_BC(self, mtrx=None):
        """
        Eliminates the rows and columns of the BC
        """
        mtrx = copy.deepcopy(mtrx)
        print(mtrx.size)
        _to_eliminate = []  # list of rows and columns to eliminate
        for nodeID, DOFs in self.supports.items():
            print(nodeID)
            for DOF in DOFs:
                print(DOF)
                _to_eliminate.append(self.position_in_matrix(nodeID=nodeID, DOF=DOF))

        _to_eliminate.sort()

        for _ in _to_eliminate[::-1]:
            mtrx = np.delete(mtrx, _, axis=0)
            mtrx = np.delete(mtrx, _, axis=1)

        return mtrx

    @property
    def M(self):
        # the compiled mass matrix
        return compile_global_matrix(self.beams, mass=True)

    @property
    def K_geom(self):
        # the compiled geometrical stiffness matrix
        return compile_global_matrix(self.beams, geometrical=True)

    @property
    def K(self):
        # the compiled stiffness matrix, without supports
        return compile_global_matrix(self.beams, stiffness=True)

    @property
    def K_with_BC(self):
        # copy of the the stiffness matrix with the boundary conditions taken into account
        _K = copy.deepcopy(self.K)
        for k, v in self.supports.items():  # k is the nodeID that has support
            for dof in v:  # dof to be fixed: 'ux', 'rotz' etc.
                _pos = self.position_in_matrix(nodeID=k, DOF=dof)
                # check, if the DOF has been released previously
                _K[_pos, _pos] += 10e20
                # print('added support: nodeID %d, DOF %s' % (k, dof))
        return _K

    @property
    def q(self):
        # the load vector
        return self._load_vector.T

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
                    # print('added: Node %d, dynam %s = %.2f' % (nodeID, k, v))

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

    def solve(self, analysis=None):
        """
        solves the system, returns the vector of displacements.
        :return: 
        """
        assert analysis in ['linear_elastic', 'modal', 'all']
        assert self.stiffness_matrix_is_ok
        assert self.node_numbering_is_ok

        if analysis in ['linear_elastic', 'all']:
            # linear static analysis
            self.displacements = np.linalg.inv(self.K_with_BC) * self.q
            self.displacements_for_beams()

        elif analysis in ['buckling']:
            raise NotImplementedError
            # Linear buckling
            # not implemented yet!

        # modal analyse
        elif analysis in ['modal', 'all']:
            eigvals, eigvecs = eigh(self.K_with_BC, self.M)
            self.frequencies = [math.sqrt(x)/(2*math.pi) for x in eigvals if x > 0]
            self.nodal_shapes = eigvecs.T

            # print(eigvals[0:5])
            # print([math.sqrt(x) for x in eigvals if x > 0][0:5])
            # print([math.sqrt(x)/(2*3.1415) for x in eigvals if x > 0][0:5])
            # print('')
            # for shind, sh in enumerate(eigvecs[0:5]):
            #
            #     print('')
            #     _ux = sh[0::3]
            #     _uy = sh[1::3]
            #     _rotz = sh[2::3]
            #     print('ux:', _ux)
            #     print('uy:', _uy)
            #     print('rotz:', _rotz)
            #
            #     # Two subplots, the axes array is 1-d
            #     f, axarr = plt.subplots(2)
            #     axarr[0].plot(list(range(len(_uy))), _ux)
            #     axarr[0].set_title('#%d, f=%.2f Hz' % (shind+1, math.sqrt(eigvals[shind]) / 2*3.1415))
            #     axarr[1].plot(list(range(len(_uy))), _uy)
            #     plt.show()


        return True

    def compile_global_geometrical_stiffness_matrix(self):
        """
        Compiles the global stiffness matrix from the element matrices.
        :return: 
        """
        return compile_global_matrix(self.beams, geometrical=True)

    def compile_global_stiffness_matrix(self):
        """
        Compiles the global stiffness matrix from the element matrices.
        :return: 
        """
        return compile_global_matrix(self.beams, stiffness=True)

    def nodenr_dof_from_position(self, position=None):
        """
        Tells the node number, DOF from the position provided.
        :param position: 
        :return: 
        """
        _nodeID, DOFnr = divmod(position, self.dof)
        return _nodeID+1, self.loadnames[DOFnr]

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


def compile_global_matrix(beams, stiffness=False, mass=False, geometrical=False):
    """
    Compiles the global stiffness matrix from the element matrices.
    :return: 
    """
    assert stiffness or mass or geometrical  # either or

    nodes = set(itertools.chain.from_iterable([x.nodes for x in beams]))  # nodes of the beams
    dof = beams[0].dof  # nr. of dofs
    _sumdof = len(nodes) * dof
    _empty = np.zeros(_sumdof ** 2)
    _empty = np.matrix(_empty.reshape(_sumdof, _sumdof))
    for b in beams:
        if stiffness:
            mtrx = b.matrix_in_global(mtrx=b.Ke)
        elif mass:
            mtrx = b.matrix_in_global(mtrx=b.Me)
        elif geometrical:
            mtrx = b.matrix_in_global(mtrx=b.Ke_geom)
        else:
            raise Exception('matrix not specified')
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
    if dof == 3:
        _block = np.matrix([[cs,    -ss,    0],
                            [ss,    cs,     0],
                            [0,     0,      1]])
    elif dof == 2:
        _block = np.matrix([[cs,    -ss],
                            [ss,    cs]])

    else:
        raise Exception('not implementeds, dof should be 2 or 3')

    _sumdof = blocks * dof
    _empty = np.zeros(_sumdof ** 2)
    _empty = np.matrix(_empty.reshape(_sumdof, _sumdof))

    for b in range(blocks):
        _sti = b * dof  # starting element of the block for node i
        _eni = _sti + dof  # end element for the block if node i
        _empty[_sti:_eni, _sti:_eni] += _block

    return _empty


def np_matrix_tolist(mtrx):
    """
    casts a numpy matrix into a flat list
    :param mtrx: 
    :return: 
    """
    return list(itertools.chain.from_iterable(mtrx.tolist()))

#


if __name__ == '__main__':

    # # nodes
    # n1 = Node.from_dict(adict={'ID': 1, 'coords': (0, 0)})  # from dict
    # n2 = Node.from_dict(adict={'ID': 2, 'coords': (250, 0)})
    # n3 = Node(ID=3, coords=(400, 0))  # direct
    #
    # # beams
    # b1 = HermitianBeam2D.from_dict(adict={'ID': 1, 'E': 210000., 'I': 39760.78, 'A': 706.5, 'rho': 7.85e-5, 'i': n1, 'j': n2})  # from dict
    # b2 = HermitianBeam2D(E=21000., ID=2, I=39760.78, A=706.5, i=n2, j=n3, rho=7.85e-5)  # direct
    #
    # # supports
    # BCs = {1: ['ux', 'uy'], 3: ['uy']}  # supports as dict
    # # BCs = {1: ['ux', 'uy', 'rotz']}  # supports as dict
    #
    # # this is the structure
    # structure = Structure(beams=[b1, b2], supports=BCs)
    #
    # # adding loads
    # structure.add_single_dynam_to_node(nodeID=2, dynam={'FX': 10000}, clear=True)  # clears previous loads
    # structure.add_single_dynam_to_node(nodeID=2, dynam={'FY': 1000})  # no clearing, just adding
    #
    # # solver :-) whatever happens here is done by numpy.
    # structure.solve()
    #
    # # sorry, no postprocessng, however you can only plot the displacements
    # pp.pprint(structure.displacements_as_dict())
    # draw_beam.draw_structure(structure)

    import time
    sta = time.time()
    # nodes

    _pieces = 20  # stk.
    _length = 1400  # mm

    _nodes = []
    for i in range(_pieces+1):
        _nodes.append(Node.from_dict(adict={'ID': i+1, 'coords': (i*_length/_pieces, 0)}))

    # beams
    # section_column = sections.Recangle(height=3, width=20)
    section_column = sections.Circle(r=55)
    print(section_column.I)
    print(section_column.A)
    # rho = 0.007850  # kg/m3
    rho = 7850000  # g/m3
    # b1 = HermitianBeam2D.from_dict(adict={'ID': 1, 'E': 2.1e11, 'I': section_column.I['x'], 'A': section_column.A, 'rho': rho, 'i': n1, 'j': n2})  # from dict
    _beams = []
    for i in range(_pieces):
        _beams.append(HermitianBeam2D.from_dict(adict={'ID': i, 'E': 2.1e11, 'I': section_column.I['x'], 'A': section_column.A, 'rho': rho, 'i': _nodes[i], 'j': _nodes[i+1]}))

    for n in _nodes:
        print(n)
    for b in _beams:
        print(b)

    # supports
    BCs = {1: ['ux', 'uy', 'rotz']}  # supports as dict
    # BCs = {1: ['ux', 'uy', 'rotz']}  # supports as dict

    # this is the structure
    structure = Structure(beams=_beams, supports=BCs)
    # structure = Structure(beams=[b1], supports=BCs)

    # adding loads
    structure.add_single_dynam_to_node(nodeID=len(_nodes)-1, dynam={'FY': -1000000}, clear=True)  # clears previous loads
    # structure.add_single_dynam_to_node(nodeID=3, dynam={'FX': 20000})

    # solver :-) whatever happens here is done by numpy.
    structure.solve(analysis='all')

    _gew = 0
    for b in structure.beams:
        _gew += b.A / 1e6 * b.rho / 1e3 * b.l / 1e3
    print('Gewicht:', _gew)
    # print(structure.M)
    # print(structure.M.diagonal())
    # print(np.sum(structure.M.diagonal()))

    print(time.time()-sta)


    # structure.draw()


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
