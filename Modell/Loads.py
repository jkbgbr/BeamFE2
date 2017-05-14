
import numpy as np
import math
from drawing import _plotting_available, plt, Polygon, PatchCollection
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
        p = PatchCollection(patches, alpha=0.4, facecolors=['lightblue'], edgecolors=['blue'])
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

    def axial_at_position(self, xi):
        # no axial action from perpendicular load
        return 0

    def moment_at_position(self, xi):
        # Schneider Bautabellen 4.2, einfeldträger mit verteilten Last
        xi_bar = 1 - xi
        return ((xi * xi_bar) * self.q * self.beam.l ** 2) / 2.

    def shear_at_position(self, xi):
        return 0

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
    pass
