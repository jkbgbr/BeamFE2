# -*- coding: utf-8 -*-
from drawing import draw_beam
import copy
from Modell.helpers import *
from Modell.Loads import *
from Modell import Loads as BL
from Modell import Results
from Modell import Solver


class Structure(object):
    """
    The Structure, composed of the FE Beams.
    """
    def __init__(self, beams=None, supports=None):
        self.beams = beams
        self.nodal_loads = []
        self._load_vector = None
        self._mass_vector = None
        self.supports = supports
        self.results = {'linear static': Results.LinearStaticResult(),
                        'modal': Results.ModalResult(),
                        'buckling': Results.BucklingResult()
                        }
        self.solver = {'linear static': Solver.LinearStaticSolver(self),
                       'modal': Solver.ModalSolver(self),
                       'buckling': Solver.BucklingSolver(self)
                       }

    def set_mass_matrix_type(self, matrixtype='consistent', beam_IDs='all'):
        if beam_IDs == 'all':
            for beam in self.beams:
                beam.mass_matrix = matrixtype
        else:
            assert all([x in [y.ID for y in self.beams] for x in beam_IDs])
            for beam in self.beams:
                if beam.ID in beam_IDs:
                    beam.mass_matrix = matrixtype

    def add_nodal_load(self, nodeID=None, dynam=None, clear=False):
        """
        
        :param nodeID: ID of the node the last acts on
        :param dynam: values of the load components
        :param local_input: True: input is understood in the beams local system. False: in the global system
        :param clear: clear all previously defined loads before applying this
        :return: 
        """

        # finding the node
        node = [x for x in self.nodes if nodeID == x.ID]
        if len(node) != 1:
            raise Exception('There is no node with ID %d' % nodeID)
        else:
            node = node[0]
        # making sure the dynam is full, if some component is mssing we replace it with a zero
        for ln in self.loadnames:
            if ln not in dynam.keys():
                dynam[ln] = 0
        self.add_single_dynam_to_node(nodeID=nodeID, dynam=dynam, clear=clear)
        self.nodal_loads.append(BL.NodalLoad(node=node, dynam=dynam))


    def mass(self):
        """ Structural mass """
        return [x.mass for x in self.beams]

    def draw(self, show=True, analysistype=None, mode=0, internal_action=None):
        if self.results[analysistype].solved:
            draw_beam.draw_structure(self, show=show, analysistype=analysistype, mode=mode, intac=internal_action)
        else:
            print('no results available, no printing')

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

    # def zero_BC(self, mtrx=None):
    #     """
    #     Eliminates the rows and columns of the BC
    #     """
    #     mtrx = copy.deepcopy(mtrx)
    #     print(mtrx.size)
    #     _to_eliminate = []  # list of rows and columns to eliminate
    #     for nodeID, DOFs in self.supports.items():
    #         print(nodeID)
    #         for DOF in DOFs:
    #             print(DOF)
    #             _to_eliminate.append(self.position_in_matrix(nodeID=nodeID, DOF=DOF))
    #
    #     _to_eliminate.sort()
    #
    #     for _ in _to_eliminate[::-1]:
    #         mtrx[_] = 0
    #         mtrx[:, _] = 0
    #
    #     return mtrx

    @property
    def positions_to_eliminate(self):
        """
        Numbers of rows and columns to be eliminated when condensing the K and M matrices for the modal analysis, based
        on the boundary conditions defined.
        Possible zero-rows are not considered here.
        Returned is a sorted list of these.
        When deleting the rows, one should begin with the highest number, that is, the reversed list of positions
        When re-populating the displacement vectors, the list should not be reversed.
        :return: the list with the numbers of the rows, columns
        """

        _to_eliminate = []  # list of rows and columns to eliminate
        for nodeID, DOFs in self.supports.items():
            for DOF in DOFs:
                _to_eliminate.append(self.position_in_matrix(nodeID=nodeID, DOF=DOF))
        _to_eliminate.sort()
        return _to_eliminate

    def condense(self, mtrx=None):
        """
        Eliminates the rows and columns of the BC
        """
        for _ in self.positions_to_eliminate[::-1]:
            mtrx = np.delete(mtrx, _, axis=0)
            mtrx = np.delete(mtrx, _, axis=1)

        return mtrx

    # def repopulate(self, mtrx):
    #     """
    #     The opposite of condense
    #     """
    #     for _ in self.positions_to_eliminate:
    #         mtrx = np.insert(mtrx, _, axis=0)
    #         mtrx = np.insert(mtrx, _, axis=1)
    #
    #     return mtrx



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
                _K[_pos, _pos] += 10e40
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

    def add_internal_loads(self, beam=None, **kwargs):
        beam.add_internal_load(**kwargs)
        for node in beam.nodes:
            dynam_as_dict = beam.reduce_internal_load(load=beam.internal_loads[-1])
            self.add_single_dynam_to_node(nodeID=node.ID, dynam=dynam_as_dict[node])

    # def reduced_internal_loads(self, load):
    #     """
    #     Reduces the internal loads to the nodes.
    #     :return:
    #     """
    #     _ret = {load.beam.i: {k: 0 for k in load.beam.dynamnames}, load.beam.j: {k: 0 for k in load.beam.dynamnames}}
    #     for component in load.beam.dynamnames:
    #         for node in self.nodes:
    #             _ret[node][component] += sum([x.reactions[node][component] for x in load])
    #     return _ret
    #     # for b in self.beams:
    #     #     print(b.reduced_internal_loads)


    #
    # @property
    # def reduced_internal_loads(self):
    #     """
    #     Summing the nodal foreces from the internal loads
    #     :return:
    #     """
    #     _ret = {self.i: {k: 0 for k in self.dynamnames}, self.j: {k: 0 for k in self.dynamnames}}
    #     for component in self.dynamnames:
    #         for node in self.nodes:
    #             _ret[node][component] += sum([x.reactions[node][component] for x in self.internal_loads])
    #
    #     return _ret

    def clear_loads(self):
        # clear all loads defined previously
        # loads defined on beams
        for b in self.beams:
            b.internal_loads = []
        self.nodal_loads = []  # deleting all nodal loads
        self._load_vector = None  # zeroing the load vector

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

        # clear all loads defined previously
        if clear:
            self.clear_loads()

        if self._load_vector is None:
            self._load_vector = np.matrix(np.zeros(self.sumdof))

        for k, v in dynam.items():
            for name, number in zip(self.loadnames, range(self.dof)):  # pairs e.g. 'FX' with 0, 'FY' with 1 etc.
                # _sti = nodeID * self.dof + number  # starting index
                if k == name:
                    _sti = self.position_in_matrix(nodeID=nodeID, dynam=k)
                    self._load_vector[0, _sti] += v

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
            # print('added: Node %d, mass %.2f' % (nodeID, mass))


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

