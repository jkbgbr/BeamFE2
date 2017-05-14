
import numpy as np
import math
from drawing import _plotting_available, plt, Polygon, PatchCollection
from Modell.helpers import *

from drawing import _plotting_available, plt
import sys
import inspect

BEAM_LOAD_TYPES = 'uniform perpendicular force', 'concentrated perpendicular force', 'concentrated moment'
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
                ax.arrow(mp[0], mp[1], _norm[0], _norm[1], head_width=0.5 * scale, head_length=scale, fc='blue', ec='blue')


class BeamLoad(object):
    def __init__(self, loadtype=None, beam=None, value=1, *args, **kwargs):
        self.loadtype = loadtype  # type of load
        self.beam = beam  # length of beam
        self.value = value

    def reduce(self):
        raise NotImplementedError

    def draw_load(self, scale=1.):
        """
        Draws the load
        :param scale: 
        :return: 
        """
        raise NotImplementedError

    @property
    def internal_points(self):
        """
        The points to calculate the internal actions - points of interest.
        :return: 
        """
        raise NotImplementedError

    @property
    def reactions_asvector(self):
        """
        Nodal reactions as a vector in the local system. 
        :return: 
        """
        raise NotImplementedError

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

    """
    The following methods are used to calculate the appropriate values at the point xi along the length of the 
    beam element. xi is in [0, 1], starting from node i.
    """

    def deflection_at_position(self, xi):
        raise NotImplementedError

    def axial_at_position(self, xi):
        raise NotImplementedError

    def moment_at_position(self, xi):
        raise NotImplementedError

    def shear_at_position(self, xi):
        raise NotImplementedError

    @property
    def deflections(self):
        return (self.deflection_at_position(x) for x in self.internal_points)

    @property
    def moments(self):
        return (self.moment_at_position(x) for x in self.internal_points)

    @property
    def shears(self):
        return (self.shear_at_position(x) for x in self.internal_points)


class UniformPerpendicularForce(BeamLoad):
    """
    Distributed, perpendicular load acting on a beam.
    Everything is calculated in the elements local system.
    The load starting and end positions a, c are understood as values normalized over the beam length L. 
    The distributed load acts on the full length of the beam with uniform intensity.
    """

    def __init__(self, loadtype='uniform perpendicular force', value=None, beam=None):
        super(UniformPerpendicularForce, self).__init__(loadtype, beam, value)
        self.nr_points = beam.number_internal_points
        self.beam = beam

    def draw_load(self, scale=1.):
        p1 = [0, 0]
        p2 = [0, scale * -self.value]
        p3 = [self.beam.l, scale * -self.value]
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
        _ret = np.matrix([[0, self.value * self.beam.l / 2, self.value * (self.beam.l ** 2) / 12,
                           0, self.value * self.beam.l / 2, - 1 * self.value * (self.beam.l ** 2) / 12.]])
        return _ret

    def deflection_at_position(self, xi):
        # Schneider Bautabellen 4.2, einfeldträger mit verteilten Last
        xi_bar = 1 - xi
        return ((1 + xi * xi_bar) / (24 * self.beam.EI)) * (xi * xi_bar) * self.value * self.beam.l

    def axial_at_position(self, xi):
        # no axial action from perpendicular load
        return 0

    def moment_at_position(self, xi):
        # Schneider Bautabellen 4.2, einfeldträger mit verteilten Last
        # the value to be added to the value interpolated between the values at node i, node j
        xi_bar = 1 - xi
        return ((xi * xi_bar) * self.value * self.beam.l ** 2) / 2.

    def shear_at_position(self, xi):
        # the value to be added to the value at node i
        return xi * self.beam.l * -self.value


class ConcentratedPerpendicularForce(BeamLoad):
    """
    Single perpendicular load acting on a beam.
    Everything is calculated in the elements local system.
    The load position a, c are understood as values normalized over the beam length L. 
    The distributed load acts on the full length of the beam with uniform intensity.
    """

    def __init__(self, loadtype='concentrated perpendicular force', value=None, position=0, beam=None):
        super(ConcentratedPerpendicularForce, self).__init__(loadtype, beam, value)
        self.position = position
        self.nr_points = beam.number_internal_points
        self.beam = beam

    @property
    def alpha(self):
        # ratio of node i to position to 1.0
        return self.position

    @property
    def beta(self):
        # ratio of node j to position to 1.0
        return 1 - self.position

    def draw_load(self, scale=1.):
        mp = [self.beta * self.beam.i.coords[0] + self.alpha * self.beam.j.coords[0],
              self.beta * self.beam.i.coords[1] + self.alpha * self.beam.j.coords[1]]
        _norm = [0, 1 * self.value / abs(self.value)] * transfer_matrix(alpha=-self.beam.direction, blocks=1, blocksize=2)
        _norm = np_matrix_tolist(_norm)
        ax = plt.gca()
        ax.arrow(mp[0], mp[1], _norm[0], _norm[1], head_width=0.5 * scale, head_length=scale, fc='blue', ec='blue')

    @property
    def internal_points(self):
        """
        The points to calculate the internal actions - points of interest.
        :return: 
        """
        # both ends and the point before and after the position of the force
        return [0, self.position-0.001, self.position+0.001, 1]

    @property
    def reactions_asvector(self):
        # the nodal forces resulting from the beam internal load,
        # acting on the nodes of the beam, in the local system
        # Schneider 4.8
        _a = self.alpha
        _b = self.beta
        _ret = np.matrix([[0, (3 - 2 * _b) * (_b ** 2) * self.value, -_a * (_b ** 2) * self.value * self.beam.l,
                           0, (3 - 2 * _a) * (_a ** 2) * self.value, -_b * (_a ** 2) * self.value * self.beam.l]])
        return _ret

    def deflection_at_position(self, xi):
        return 0

    def axial_at_position(self, xi):
        # no axial action from perpendicular load
        return 0

    def moment_at_position(self, xi):
        # added to the "baseline"
        _ret = self.beta * self.value * xi
        if xi > self.position:
            _ret -= abs((xi - self.position)) * self.value
        _ret *= self.beam.l
        return _ret

    def shear_at_position(self, xi):
        # the value to be added to the value at node i
        if xi < self.position:
            _ret = 0
        else:
            _ret = -self.value
        return _ret


class ConcentratedMoment(BeamLoad):
    """
    Single perpendicular load acting on a beam.
    Everything is calculated in the elements local system.
    The load position a, c are understood as values normalized over the beam length L. 
    The distributed load acts on the full length of the beam with uniform intensity.
    """

    def __init__(self, loadtype='concentrated moment', value=None, position=0, beam=None):
        super(ConcentratedMoment, self).__init__(loadtype, beam, value)
        self.position = position
        self.nr_points = beam.number_internal_points
        self.beam = beam

    @property
    def alpha(self):
        # ratio of node i to position to 1.0
        return self.position

    @property
    def beta(self):
        # ratio of node j to position to 1.0
        return 1 - self.position

    def draw_load(self, scale=1.):
        mp = [self.beta * self.beam.i.coords[0] + self.alpha * self.beam.j.coords[0],
              self.beta * self.beam.i.coords[1] + self.alpha * self.beam.j.coords[1]]
        ax = plt.gca()
        ax.scatter(mp[0], mp[1], marker='o', color='y', s=30, zorder=3)

    @property
    def internal_points(self):
        """
        The points to calculate the internal actions - points of interest.
        :return: 
        """
        # both ends and the point before and after the position of the force
        return [0, self.position-0.001, self.position+0.001, 1]

    @property
    def reactions_asvector(self):
        # the nodal forces resulting from the beam internal load,
        # acting on the nodes of the beam, in the local system
        # Schneider 4.8
        _a = self.alpha
        _b = self.beta
        MR = 6 * _a * _b * self.value / self.beam.l
        _ret = np.matrix([[0, -MR, -1 * (3 * _b - 2) * _b * self.value,
                           0, MR, 1 * (3 * _a - 2) * _a * self.value]])
        return _ret

    def deflection_at_position(self, xi):
        return 0

    def axial_at_position(self, xi):
        # no axial action from perpendicular load
        return 0

    def moment_at_position(self, xi):
        # added to the "baseline"
        _ret = -(self.value / self.beam.l) * xi
        if xi > self.position:
            _ret += self.value
        return _ret

    def shear_at_position(self, xi):
        # the value to be added to the value at node i
        return 0


if __name__ == '__main__':
    pass
