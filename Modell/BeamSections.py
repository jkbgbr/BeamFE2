
import math


class Crossection(object):
    def __init__(self, shape=None, *args, **kwargs):
        self.shape = shape

    def A(self):
        raise NotImplementedError

    def I(self):
        raise NotImplementedError


class Rectangle(Crossection):
    def __init__(self, shape='rectangle', height=None, width=None):
        super(Rectangle, self).__init__(shape=shape)
        # assert a >= b
        self.height = height
        self.width = width

    def __repr__(self):
        return 'Rectangle(height=%.2f, width=%.2f)' % (self.height, self.width)

    def __str__(self):
        return 'Rectangle %.2f x %.2f, A=%.2f, Ix=%.2f, Iy=%.2f' % (self.height, self.width, self.A, self.I['x'], self.I['y'])

    @property
    def A(self):
        return self.height * self.width

    @property
    def I(self):
        _ret = {'x': (self.height ** 3) * self.width / 12, 'y': (self.width ** 3) * self.height / 12}
        return _ret


class Circle(Crossection):
    def __init__(self, shape='circle', r=None):
        super(Circle, self).__init__(shape=shape)
        self.r = r

    def __repr__(self):
        return 'Circle(r=%.2f)' % self.r

    @property
    def A(self):
        return self.r ** 2 * math.pi

    @property
    def I(self):
        _ret = {}
        _ret['x'] = _ret['y'] = (self.r ** 4) * math.pi / 4
        return _ret


class HollowCircle(Crossection):
    def __init__(self, shape='circle', r=None, s=None):
        super(HollowCircle, self).__init__(shape=shape)
        self.r_out = r
        self.s = s
        self.r_in = s

    def __repr__(self):
        return 'HollowCircle(r=%.2f, s=%.2f)' % (self.r_out, self.s)

    @property
    def A(self):
        return (self.r_out ** 2 - self.r_in ** 2) * math.pi

    @property
    def I(self):
        _ret = {}
        _ret['x'] = _ret['y'] = (self.r_out ** 4 - self.r_in ** 4) * math.pi / 4
        return _ret
