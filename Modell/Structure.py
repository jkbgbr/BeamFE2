# -*- coding: utf-8 -*-
from drawing import draw_beam
import copy
from Modell.helpers import *


class Structure(object):
    """
    The Structure, composed of the FE Beams.
    """
    def __init__(self, beams=None, supports=None):
        self.beams = beams
        self._load_vector = None
        self._mass_vector = None
        self.supports = supports
        self.results = {'linear static': None, 'modal': None, 'buckling': None}

    def mass(self):
        """ Structural mass """
        return [x.mass for x in self.beams]

    def draw(self, show=True, analysistype=None, mode=0):
        draw_beam.draw_structure(self, show=show, analysistype=analysistype, mode=mode)

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

    # todo: string helyett valtozonev, konstanssal
    # type hinting-et hasznalni, de csak ott igazan fontos amiket kivulrol is hivhatnak
    # egyszerubb refraktor

    @property
    def dofnames(self):
        assert self.dof in [3, 6]
        if self.dof == 3:
            return 'ux', 'uy', 'rotz'
        else:
            return 'ux', 'uy', 'uz', 'rotx', 'roty', 'rotz'

    # @property
    # def dofnumbers(self):
    #     assert self.dof in [3, 6]
    #     if self.dof == 3:
    #         return 0, 1, 2
    #     else:
    #         return 0, 1, 2, 3, 4, 5

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
    def M_with_masses(self):
        # copy of the the stiffness matrix with the boundary conditions taken into account
        _M = copy.deepcopy(self.M)

        if self._mass_vector is None:
            self._mass_vector = np.matrix(np.zeros(self.sumdof))

        for mindex, m in enumerate(np_matrix_tolist(self.mass_vector)):
            _M[mindex, mindex] += m

        # for k, v in self.supports.items():  # k is the nodeID that has support
        #     for dof in v:  # dof to be fixed: 'ux', 'rotz' etc.
        #         _pos = self.position_in_matrix(nodeID=k, DOF=dof)
        #         # check, if the DOF has been released previously
        #         _M[_pos, _pos] += 10e20

        return _M

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
    def load_vector(self):
        # the load vector
        return self._load_vector.T

    @property
    def mass_vector(self):
        # a vector containing the masses
        # currently no distinction between directions, defined massesact in all directions
        return self._mass_vector

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

    def add_mass_to_node(self, nodeID=None, mass=None, clear=False):
        """
        Adds a dynam (FX, FY, MZ) to the chosen Node of the model.
        Checks if the node exists.
        Checks if the name of the dynam is OK.
        clears previous loads on the node, if clear is True
        
        :param nodeID: ID of the node the mass is added to
        :param mass: mass to be added in [kg]
        :param clear: if True, all previously defined masses will be deleted
        :return: 
        """
        #

        assert nodeID in [x.ID for x in self.nodes]

        if self._mass_vector is None:
            self._mass_vector = np.matrix(np.zeros(self.sumdof))

        # clear all loads defined previously
        if clear:
            # print('loads cleared')
            self._mass_vector[0, :-1] = 0

        _sti = self.position_in_matrix(nodeID=nodeID, DOF='ux')
        for p in range(3):
            self._mass_vector[0, _sti + p] += mass
            print('added: Node %d, mass %.2f' % (nodeID, mass))

        #
        # for k, v in mass.items():
        #
        #     for name, number in zip(self.loadnames, range(self.dof)):  # pairs e.g. 'FX' with 0, 'FY' with 1 etc.
        #         # _sti = nodeID * self.dof + number  # starting index
        #         if k == name:
        #             _sti = self.position_in_matrix(nodeID=nodeID, dynam=k)
        #             self._load_vector[0, _sti] += v

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
    def mass_matrix_is_ok(self):
        return True
        # todo: check for symmetry, positive definiteness

    @property
    def stiffness_matrix_is_ok(self):
        if self.stiffness_matrix_is_nonsingular:
            return True
        else:
            print('The stiffness matrix is singular')
            return False

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

