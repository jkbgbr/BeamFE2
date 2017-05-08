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
    number_internal_points = 10

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
        # self.displacements = None  # displacements in the GLOBAL system, provided by the solver
        self.results = {'linear static': None, 'modal': None, 'buckling': None}

    def __repr__(self):
        return 'HermitianBeam2D(ID=%d, i=%r, j=%r, E=%d, I=%.2f, A=%.2f' \
               % (self.ID, self.i, self.j, self.E, self.I, self.A)

    def __str__(self):
        return 'HB2D_%d' % self.ID

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
    def mass(self):
        return self.A / 1e6 * self.rho / 1e3 * self.l / 1e3

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
        return transfer_matrix(alpha=self.direction, asdegree=False, blocksize=self.dof, blocks=2)

    def matrix_in_global(self, mtrx=None):
        # matrix of the Beam element in the global coordinate system
        return self.transfer_matrix * mtrx * self.transfer_matrix.T

    # shape functions. the input argument x is the "running parameter", a value normalized with the beam length
    def N1(self, x, L=1.):
        """
        Shape function for axial displacements, node 'i'
        :param x: position of the internal point along the length, 0.0 <= x <= 1.0
        :param L: length of the beam
        :return: value of the shape function at x
        """
        x *= L
        return 1 - x / L

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
        return x / L
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
        return np.matrix([[self.N1(x=x, L=L),  0,                  0,                     self.N4(x=x, L=L),  0,                  0],
                          [0,                  self.N2(x=x, L=L),  self.N3(x=x, L=L),     0,                  self.N5(x=x, L=L),  self.N6(x=x, L=L)]])

    def _Me(self):
        L = self.l
        _ret = self.rho * self.A * self.l / 420 * \
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

        # http://opensees.berkeley.edu/WebSVN/filedetails.php?repname=OpenSees&path=%2Ftrunk%2FSRC%2Felement%2FelasticBeamColumn%2FElasticTimoshenkoBeam2d.cpp

        # // compute initial stiffness matrix in local system
        # kl.Zero();
        # kl(0,0) = kl(3,3) = E*A/L;
        # kl(0,3) = kl(3,0) = -kl(0,0);
        # double a1z = E*Iz/(1.0 + phi);
        # kl(1,1) = kl(4,4) = 12.0*a1z/(L*L*L);
        # kl(1,4) = kl(4,1) = -kl(1,1);
        # kl(2,2) = kl(5,5) = (4.0 + phi)*a1z/L;
        # kl(2,5) = kl(5,2) = (2.0 - phi)*a1z/L;
        # kl(1,2) = kl(2,1) = kl(1,5) = kl(5,1) = 6.0*a1z/(L*L);
        # kl(2,4) = kl(4,2) = kl(4,5) = kl(5,4) = -kl(1,2);
        #
        # // compute geometric stiffness matrix in local system
        # klgeo.Zero();
        # if (nlGeo == 1)  {
        #     double b1z = 1.0/(30.0*L*pow(1.0 + phi,2));
        #     klgeo(1,1) = klgeo(4,4) = b1z*(30.0*phi*phi + 60.0*phi + 36.0);
        #     klgeo(1,4) = klgeo(4,1) = -klgeo(1,1);
        #     klgeo(2,2) = klgeo(5,5) = b1z*L*L*(2.5*phi*phi + 5.0*phi + 4.0);
        #     klgeo(2,5) = klgeo(5,2) = -b1z*L*L*(2.5*phi*phi + 5.0*phi + 1.0);
        #     klgeo(1,2) = klgeo(2,1) = klgeo(1,5) = klgeo(5,1) = 3.0*L;
        #     klgeo(2,4) = klgeo(4,2) = klgeo(4,5) = klgeo(5,4) = -klgeo(1,2);
        # }
    # def _Ke_geom(self):
    #     # the geometrical stiffness matrix, from H-P. Gavin CEE 421L. Matrix Structural Anyalsis - Duke University
    #     _locdisp = self.local_displacements
    #     T = _locdisp[3] - _locdisp[0]  # change of the axial length
    #     L = self.l
    #     _ret = np.matrix([
    #         [0, 0, 0, 0, 0, 0],
    #         [0, 6./5., L/10., 0, -(6./5.), L/10.],
    #         [0, L/10., 2*(L**2)/15., 0, -(L/10.), -(L/30.)],
    #         [0, 0, 0, 0, 0, 0],
    #         [0, -(6./5.), -(L/10.), 0, (6./5.), -(L/10.)],
    #         [0, (L/10.), -(L**2)/30., 0, -(L/10.), -(2*(L**2))/15.]
    #     ])
    #     _ret = np.multiply((T / L), _ret)
    #
    #     return _ret
    #
    # @property
    # def Ke_geom(self):
    #     return self._Ke_geom()

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

    def deflected_shape(self, local=True, scale=1., disps=None):
        """
        The deflected shape of the beam, based on the local, non-partitioned displacements of the end nodes, 
        provided in disps.
        (ux, uy) = N * load_vector
        :param local: if local = False, the values are transferred in the global coordinate system e.g. for plotting
        :param scale: displacement magnification scale
        :param disps: nodal displacements in the local system
        :return: 
        """
        _ip = self.number_internal_points
        _deflected_shape = []

        if local:
            _t = transfer_matrix(alpha=0, asdegree=False, blocks=1, blocksize=2)  # 2x2 transform matrix
        else:
            _t = transfer_matrix(self.direction, asdegree=False, blocks=1, blocksize=2)  # 2x2 transform matrix

        for i in range(_ip + 1):  # e.g. 0, 1, 2, 3 for _ip == 3
            _val = self.N(x=i / _ip, L=self.l) * disps  # displacements in the local system
            _val = np.multiply(scale, _val)
            _val[0, 0] += (i / _ip) * self.l  # the deflected shape in the local coordinate system
            _val = _t * _val  # the deflected shape rotated (possibly by 0 degree)
            if not local:
                _val = np_matrix_tolist(_val)  # type casting
                _val = (_val[0] + self.i.x, _val[1] + self.i.y)  # translating the start point to node i
            _deflected_shape.append(_val)
        return _deflected_shape

    def nodal_reactions(self, disps):
        """ reactions in the local coordinate system of the beam """
        return self.Ke * disps


if __name__ == '__main__':
    pass
