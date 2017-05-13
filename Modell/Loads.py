
import numpy as np
import math
from matplotlib.patches import Circle, Wedge, Polygon
from matplotlib.collections import PatchCollection
from Modell.helpers import *

from drawing import _plotting_available, plt
import sys
import inspect

BEAM_LOAD_TYPES = 'uniform perpendicular force'
NODAL_LOAD_TYPES = ['force', 'moment']


class NodalLoad(object):
    def __init__(self, node=None, dynam=None, *args, **kwargs):
        self.node = node  # length of beam
        self.dynam = dynam

    @property
    def asvector(self):
        d = self.dynam
        return np.matrix([d['FX'], d['FY'], d['MZ']])

    def draw_load(self, scale=1.):
        mp = self.node.coords  # starting point of the arrow
        for component, load in self.dynam.items():
            if load:
                if component == 'FX':
                    _norm = (load * scale / abs(load), 0,)
                elif component == 'FY':
                    _norm = (0, load * scale / abs(load))
                else:
                    _norm = False
                # plotting, if there is a norm
                ax = plt.gca()
                ax.arrow(mp[0], mp[1], _norm[0], _norm[1], head_width=0.5 * scale, head_length=scale, fc='r', ec='r')


class BeamLoad(object):
    def __init__(self, loadtype=None, beam=None, *args, **kwargs):
        self.loadtype = loadtype  # type of load
        self.beam = beam  # length of beam

    def reduce(self):
        raise NotImplementedError


class UniformPerpendicularForce(BeamLoad):
    """
    Distributed, perpendicular load acting on a beam.
    Everything is calculated in the elements local system.
    The load starting and end positions a, c are understood as values normalized over the beam length L. 
    The distributed load acts on the full length of the beam with uniform intensity.
    """

    def __init__(self, loadtype='uniform perpendicular force', q=1, beam=None):
        super(UniformPerpendicularForce, self).__init__(loadtype, beam)
        self.q = q
        self.nr_points = beam.number_internal_points
        self.beam = beam

    def draw_load(self, scale=1.):
        p1 = [0, 0]
        p2 = [0, scale * -self.q]
        p3 = [self.beam.l, scale * -self.q]
        p4 = [self.beam.l, 0]
        pts = [p1, p2, p3, p4]
        _tr = transfer_matrix(alpha=-self.beam.direction, asdegree=False, blocks=1, blocksize=2)
        pts = [np_matrix_tolist(x * _tr + self.beam.i.coords) for x in pts]
        polygon = Polygon(pts, True)
        patches = [polygon]
        p = PatchCollection(patches, alpha=0.4)
        # p.set_array(np.array('b'))
        ax = plt.gca()
        ax.add_collection(p)

    @property
    def internal_points(self):
        """
        The points to calculate the internal actions - points of interest.
        :return: 
        """
        return [x / self.nr_points for x in range(self.nr_points+1)]  # points where the load is

    @property
    def reactions_asvector(self):
        # the nodal forces resulting from the beam internal load,
        # acting on the nodes of the beam, in the local system
        _ret = np.matrix([[0, self.q * self.beam.l / 2, self.q * (self.beam.l ** 2) / 12,
                           0, self.q * self.beam.l / 2, - 1 * self.q * (self.beam.l ** 2) / 12.]])
        return _ret

    @property
    def reactions(self):
        """
        reactions at node i for the DOFs, in the global system these will only be used in the assembly of the global load vector
        :return: 3-tuple for local FX, FY, MZ 
        """
        _av = self.reactions_asvector * transfer_matrix(alpha=-self.beam.direction)
        _ret = {
            self.beam.i: {'FX': _av[0, 0], 'FY': _av[0, 1], 'MZ': _av[0, 2]},
            self.beam.j: {'FX': _av[0, 3], 'FY': _av[0, 4], 'MZ': _av[0, 5]}
                }

        return _ret

    def deflection_at_position(self, xi):
        # Schneider Bautabellen 4.2, einfeldträger mit verteilten Last
        xi_bar = 1 - xi
        return ((1 + xi * xi_bar) / (24 * self.beam.EI)) * (xi * xi_bar) * self.q * self.beam.l

    def moment_at_position(self, xi):
        # Schneider Bautabellen 4.2, einfeldträger mit verteilten Last
        xi_bar = 1 - xi
        return ((xi * xi_bar) * self.q * self.beam.l ** 2) / 2.

    def shear_at_position(self, xi):
        return (self.q * self.beam.l) / 2. - xi * self.beam.l * self.q

    @property
    def deflections(self):
        return (self.deflection_at_position(x) for x in self.internal_points)

    @property
    def moments(self):
        return (self.moment_at_position(x) for x in self.internal_points)

    @property
    def shears(self):
        return (self.shear_at_position(x) for x in self.internal_points)


if __name__ == '__main__':


    exit()
#
#
#     _a = 0.0
#     L = 5.
#     _c = 0.1
#     q1 = 10000
#     q2 = 10000
#     a = _a * L
#     c = _c * L
#     b = L - a - c
#     A0 = q2 * b / 2 * (b / 3 + c) / L + q1 * b / 2 * (2 * b / 3 + c) / L
#     C0 = q2 * b / 2 * (a + 2 / 3 * b) / L + q1 * b / 2 * (a + 1 / 3 * b) / L
#
#     print(A0)
#     print(C0)
#     _form = []
#     for xi in range(101):
#
#         x0 = (xi / 100) * L
#
#         if x0 <= a:
#             dyl = (1 / 2 * A0 * L * a * a - 1 / 8 * L * q1 * b * b * b - 1 / 24 * L * q2 * b * b * b + 1 / 8 * q1 * a * b * b * b - 2 / 3 * C0 * L * L * L + 11 / 120 * q1 * b * b * b * b + 1 / 30 * q2 * b * b * b * b + 1 / 24 * q2 * a * b * b * b - 1 / 6 * A0 * x0 * x0 * L - 1 / 3 * A0 * a * a * a - 1 / 3 * A0 * b * b * b - 1 / 3 * C0 * a * a * a + L * L * C0 * c - A0 * b * a * a - A0 * b * b * a + 1 / 2 * L * A0 * b * b - C0 * a * a * b - C0 * a * b * b + C0 * L * a * a + C0 * L * b * b + L * A0 * b * a + 2 * C0 * L * a * b - 1 / 3 * C0 * b * b * b) *x0 / L
#
#         elif x0 <= (a + b):
#
#             dyl = 1 / 120 * (L * q2 - L * q1) / b / L * (x0 ** 5) +1 / 120 * (5 * L * q1 * a + 5 * L * q1 * b - 5 * L * q2 * a) / b / L * (x0 ** 4) +1 / 120 * (-20 * L * A0 * b - 10 * L * q1 * a * a - 20 * L * q1 * a * b + 10 * L * q2 * a * a) / b / L * (x0 ** 3) +1 / 120 * (10 * L * q1 * a * a * a - 10 * L * q2 * a * a * a + 30 * L * q1 * b * a * a) / b / L * x0 * x0 +1 / 120 * (
#                 -5 * L * q1 * a * a * a * a - 20 * L * q1 * a * a * a * b - 5 * L * q2 * b * b * b * b + 5 * L * q2 * a * a * a * a - 120 * A0 * b * b * a * a
#                 - 120 * A0 * b * b * b * a - 40 * A0 * b * b * b * b - 15 * L * q1 * b * b * b * b - 40 * A0 * a * a * a * b + 60 * L * A0 * b * a * a
#                 + 120 * L * A0 * b * b * a + 60 * L * A0 * b * b * b - 40 * b * b * b * b * C0 + 4 * q2 * b * b * b * b * b + 11 * q1 * b * b * b * b * b
#                 - 40 * a * a * a * b * C0 + 5 * q2 * a * b * b * b * b - 120 * a * a * b * b * C0 - 120 * a * b * b * b * C0 + 15 * q1 * a * b * b * b * b
#                 - 80 * C0 * b * L * L * L + 120 * C0 * b * b * b * L + 120 * L * L * C0 * c * b + 120 * C0 * b * L * a * a + 240 * C0 * b * b * L * a) / b / L * x0 +1 / 120 * (L * q1 * a * a * a * a * a - L * q2 * a * a * a * a * a + 5 * L * q1 * b * a * a * a * a) / b / L
#
#         else:
#
#             dyl = 1 / 120 * (L - x0) * (40 * L * L * C0 * x0 - 120 * L * C0 * a * b - 60 * L * C0 * a * a - 60 * L * C0 * b * b - 20 * L * C0 * x0 * x0 - 5 * q2 * a * b * b * b - 4 * q2 * b * b * b * b + 120 * a * a * b * C0 - 11 * q1 * b * b * b * b + 120 * A0 * b * a * a + 120 * a * b * b * C0 + 40 * C0 * a * a * a + 40 * b * b * b * C0 - 15 * q1 * a * b * b * b + 40 * A0 * b * b * b + 40 * a * a * a * A0 + 120 * A0 * b * b * a) / L
#
#         _form.append(dyl)
#
#     print(_form)
#     plt.plot(_form)
#     plt.show()
#
#     s = 0.5 * (q1 + q2) * b
#     R1 = (q2 * b / 2 * (b / 3 + c) + q1 * b / 2 * (2 * b / 3 + c)) / L
#     _R2 = (q2 * b / 2 * (a + 2 / 3 * b) + q1 * b / 2 * (a + b / 3)) / L
#     R2 = s - R1
#
#     print(_R2 == R2)
#
#     nx = 0
#
#     _txlist = []
#     _mxlist = []
#     for xi in range(101):
#
#         x = (xi / 100) * L
#
#         if x <= a:
#             tx = 0
#             mx = R1 * x
#
#         if a < x < (a + b):
#             tx = 0.5 * (q1 + (q1 + (q2 - q1) * (x - a) / b)) * (x - a)
#             mx = x * R1 - q2 / 6 / b * ((x - a) ** 3) - q1 / 2 * ((x - a) ** 2) + q1 / 6 / b * ((x - a) ** 3)
#
#         if x >= (a + b):
#             tx = s
#             mx = R2 * (L - x)
#
#         _txlist.append(tx)
#         _mxlist.append(mx)
#
#
#     plt.plot(_txlist)
#     plt.plot(_mxlist)
#     plt.show()
#
#
# # 3: begin // dist
# # load
# # perp
# # q1: = loadprops[ilc, nb, j, 1];
# # q2: = loadprops[ilc, nb, j, 2];
# # a: = loadprops[ilc, nb, j, 3] * L;
# # c: = (1 - loadprops[ilc, nb, j, 4]) * L;
# # b: = L - a - c;
# # s: = 0.5 * (q1 + q2) * b;
# # R1: = (q2 * b / 2 * (b / 3 + c) + q1 * b / 2 * (2 * b / 3 + c)) / L;
# # R2: = (q2 * b / 2 * (a + 2 / 3 * b) + q1 * b / 2 * (a + b / 3)) / L;
# # R2: = s - R1;
# # nx: = 0;
# # if x <= a then tx:=
# #     0;
# # if (x > a) and (x < (a + b)) then tx:=
# #     0.5 * (q1 + (q1 + (q2 - q1) * (x - a) / b)) * (x - a);
# # if x >= (a + b) then tx:=
# #     s;
# # if x <= a then mx:=
# #     R1 * x;
# # if (x > a) and (x < (a + b)) then mx:=
# #     x * R1 - q2 / 6 / b * PowInt(x - a, 3) - q1 / 2 * PowInt(x - a, 2) + q1 / 6 / b * PowInt(x - a, 3);
# # if x >= (a + b) then mx:=
# #     R2 * (L - x);
# # end;
#
#
#
#
#
#
#
#
# #
# # #
# #     if x0 <= a then
# #     dyl = dyl + (
# #                  1 / 2 * A0 * L * a * a - 1 / 8 * L * q1 * b * b * b - 1 / 24 * L * q2 * b * b * b + 1 / 8 * q1 * a * b * b * b
# #                  - 2 / 3 * C0 * L * L * L + 11 / 120 * q1 * b * b * b * b + 1 / 30 * q2 * b * b * b * b + 1 / 24 * q2 * a * b * b * b
# #                  - 1 / 6 * A0 * x0 * x0 * L - 1 / 3 * A0 * a * a * a - 1 / 3 * A0 * b * b * b - 1 / 3 * C0 * a * a * a + L * L * C0 * c
# #                  - A0 * b * a * a - A0 * b * b * a + 1 / 2 * L * A0 * b * b - C0 * a * a * b - C0 * a * b * b + C0 * L * a * a
# #                  + C0 * L * b * b + L * A0 * b * a + 2 * C0 * L * a * b - 1 / 3 * C0 * b * b * b) * x0 / L
# # if (x0 > a) and (x0 < (a + b)) then
# # dyl = dyl + 1 / 120 * (L * q2 - L * q1) / b / L * PowInt(x0, 5)
# # +1 / 120 * (5 * L * q1 * a + 5 * L * q1 * b - 5 * L * q2 * a) / b / L * PowInt(x0, 4)
# # +1 / 120 * (-20 * L * A0 * b - 10 * L * q1 * a * a - 20 * L * q1 * a * b + 10 * L * q2 * a * a) / b / L * PowInt(x0, 3)
# # +1 / 120 * (10 * L * q1 * a * a * a - 10 * L * q2 * a * a * a + 30 * L * q1 * b * a * a) / b / L * x0 * x0
# # +1 / 120 * (
# # -5 * L * q1 * a * a * a * a - 20 * L * q1 * a * a * a * b - 5 * L * q2 * b * b * b * b + 5 * L * q2 * a * a * a * a - 120 * A0 * b * b * a * a
# # - 120 * A0 * b * b * b * a - 40 * A0 * b * b * b * b - 15 * L * q1 * b * b * b * b - 40 * A0 * a * a * a * b + 60 * L * A0 * b * a * a
# # + 120 * L * A0 * b * b * a + 60 * L * A0 * b * b * b - 40 * b * b * b * b * C0 + 4 * q2 * b * b * b * b * b + 11 * q1 * b * b * b * b * b
# # - 40 * a * a * a * b * C0 + 5 * q2 * a * b * b * b * b - 120 * a * a * b * b * C0 - 120 * a * b * b * b * C0 + 15 * q1 * a * b * b * b * b
# # - 80 * C0 * b * L * L * L + 120 * C0 * b * b * b * L + 120 * L * L * C0 * c * b + 120 * C0 * b * L * a * a + 240 * C0 * b * b * L * a) / b / L * x0
# # +1 / 120 * (L * q1 * a * a * a * a * a - L * q2 * a * a * a * a * a + 5 * L * q1 * b * a * a * a * a) / b / L
# # if x0 >= (a + b) then
# # dyl = dyl + 1 / 120 * (L - x0) * (40 * L * L * C0 * x0 - 120 * L * C0 * a * b - 60 * L * C0 * a * a
# #                                    - 60 * L * C0 * b * b - 20 * L * C0 * x0 * x0 - 5 * q2 * a * b * b * b
# #                                    - 4 * q2 * b * b * b * b + 120 * a * a * b * C0 - 11 * q1 * b * b * b * b
# #                                    + 120 * A0 * b * a * a + 120 * a * b * b * C0 + 40 * C0 * a * a * a
# #                                    + 40 * b * b * b * C0 - 15 * q1 * a * b * b * b + 40 * A0 * b * b * b
# #                                    + 40 * a * a * a * A0 + 120 * A0 * b * b * a) / L
# # end