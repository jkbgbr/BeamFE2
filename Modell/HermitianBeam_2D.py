# -*- coding: utf-8 -*-

from __future__ import division
from Modell.helpers import *
import copy
from Modell import BeamSections as sections
from Modell import Structure
from Modell import Node
from solver import solve


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

        self.rho = rho  # density in g/m3 for the nodal analysis

        # displacements are the result from the analysis. Its a dictionary, where the keys are
        # the names of the analyses, and the results themselves are in the LOCAL system, separated
        # according to DOF.
        self.displacements = None  # displacements in the GLOBAL system, provided by the solver

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

    @property
    def direction(self):
        # the direction of the beam in the global coordinate system
        # this is a crucial point for the transfer from the local to the global systems
        _dy = self.j.y - self.i.y
        _dx = self.j.x - self.i.x
        return np.arctan2(_dy, _dx)

    #
    # diverse matrices
    #

    @property
    def transfer_matrix(self):
        # matrix to rotate the stiffness matrix for compilation
        return transfer_matrix(alpha=self.direction, asdegree=False, dof=self.dof, blocks=2)

    def matrix_in_global(self, mtrx=None):
        # matrix of the Beam element in the global coordinate system
        return self.transfer_matrix * mtrx * self.transfer_matrix.T

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
        return np.matrix([[self.N1(x=x, L=L),  0,                  0,
                           self.N4(x=x, L=L),  0,                  0],
                          [0,                  self.N2(x=x, L=L),  self.N3(x=x, L=L),
                           0,                  self.N5(x=x, L=L),  self.N6(x=x, L=L)]])

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

    #
    # results
    #

    def deflected_shape(self, local=True, scale=1.):
        """
        The deflected shape of the beam. based on the local_displacements.
        (ux, uy) = N * q
        if local = False, the values are transferred in the global coordinate system e.g. for plotting
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
    # solve(structure, analysis='linear static')
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
        _nodes.append(Node.Node.from_dict(adict={'ID': i+1, 'coords': (i*_length/_pieces, 0)}))

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
    structure = Structure.Structure(beams=_beams, supports=BCs)
    # structure = Structure(beams=[b1], supports=BCs)

    # adding loads
    structure.add_single_dynam_to_node(nodeID=len(_nodes)-1, dynam={'FY': -1000000}, clear=True)  # clears previous loads
    # structure.add_single_dynam_to_node(nodeID=3, dynam={'FX': 20000})

    # solver :-) whatever happens here is done by numpy.
    solve(structure, analysis='all')

    _gew = 0
    for b in structure.beams:
        _gew += b.A / 1e6 * b.rho / 1e3 * b.l / 1e3
    print('Gewicht:', _gew)
    # print(structure.M)
    # print(structure.M.diagonal())
    # print(np.sum(structure.M.diagonal()))

    print(time.time()-sta)

    structure.draw()

