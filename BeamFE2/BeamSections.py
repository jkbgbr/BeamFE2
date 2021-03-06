
import math

# second moment of inertia
from collections import namedtuple
I = namedtuple('I', ['x', 'y'])


class Crossection(object):

    def __init__(self, *args, **kwargs):
        self.shape = self.__class__.__name__

    @property
    def A(self):
        raise NotImplementedError

    @property
    def I(self):
        raise NotImplementedError


class CustomCrosssection(Crossection):
    def __init__(self, A=None, I=None):
        assert A is not None
        assert I is not None
        assert A > 0
        assert I > 0
        super(CustomCrosssection, self).__init__()
        self._A = A
        self._I = I

    @property
    def A(self):
        return self._A

    @property
    def I(self):
        return self._I


class Rectangle(Crossection):
    def __init__(self, height=None, width=None):
        super(Rectangle, self).__init__()
        self.height = height
        self.width = width

    def __repr__(self):
        return 'Rectangle(height=%.2f, width=%.2f)' % (self.height, self.width)

    def __str__(self):
        return 'Rectangle %.2f x %.2f, A=%.2f, Ix=%.2f, Iy=%.2f' % (self.height, self.width, self.A, self.I.x, self.I.y)

    @property
    def A(self):
        return self.height * self.width

    @property
    def I(self):
        return I(x=(self.height ** 3) * self.width / 12, y=(self.width ** 3) * self.height / 12)
        # _ret = {'x': (self.height ** 3) * self.width / 12, 'y': (self.width ** 3) * self.height / 12}
        # return _ret


class Circle(Crossection):
    def __init__(self, r=None):
        super(Circle, self).__init__()
        self.r = r

    def __repr__(self):
        return 'Circle(r=%.2f)' % self.r

    def __str__(self):
        return 'Circle r = %.2f, A=%.2f, Ix=%.2f, Iy=%.2f' % (self.r, self.A, self.I.x, self.I.y)

    @property
    def A(self):
        return self.r ** 2 * math.pi

    @property
    def I(self):
        _val = (self.r ** 4) * math.pi / 4
        return I(x=_val, y=_val)

    # @property
    # def I(self):
    #     _ret = {}
    #     _ret['x'] = _ret['y'] = (self.r ** 4) * math.pi / 4
    #     return _ret


class HollowCircle(Crossection):
    def __init__(self, r=None, s=None):
        super(HollowCircle, self).__init__()
        self.r_out = r
        self.s = s
        self.r_in = r - s

    def __repr__(self):
        return 'HollowCircle(r=%.2f, s=%.2f)' % (self.r_out, self.s)

    def __str__(self):
        return 'Hollow Circle r = %.2f, s = %.2f, A=%.2f, Ix=%.2f, Iy=%.2f' % (self.r_out, self.s, self.A, self.I.x, self.I.y)

    @property
    def A(self):
        return (self.r_out ** 2 - self.r_in ** 2) * math.pi

    @property
    def I(self):
        _val = (self.r_out ** 4 - self.r_in ** 4) * math.pi / 4
        return I(x=_val, y=_val)


def run():
    pass

if __name__ == '__main__':

    run()