
import math

# general

class Crossection(object):
    def __init__(self, shape=None, *args, **kwargs):
        self.shape = shape

    def A(self):
        raise NotImplementedError

    def I(self):
        raise NotImplementedError

    def W(self):
        raise NotImplementedError


class Recangle(Crossection):
    def __init__(self, shape='rectangle', a=None, b=None):
        super(Recangle, self).__init__(shape=shape)
        assert a >= b
        self.a = a
        self.b = b

    def __repr__(self):
        return 'Rectangle(a=%.2f, b=%.2f)' % (self.a, self.b)

    def __str__(self):
        return 'Rectangle %.2f x %.2f, A=%.2f, Ix=%.2f, Iy=%.2f' % (self.a, self.b, self.A, self.I['x'], self.I['y'])

    @property
    def A(self):
        return self.a * self.b

    @property
    def I(self):
        _ret = {'x': (self.a ** 3) * self.b / 12, 'y': (self.b ** 3) * self.a / 12}
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
        return 'HollowCircle(r=%.2f, s=%.2f)' % (self.r, self.s)

    @property
    def A(self):
        return (self.r_out ** 2 - self.r_in ** 2) * math.pi

    @property
    def I(self):
        _ret = {}
        _ret['x'] = _ret['y'] = (self.r_out ** 4 - self.r_in ** 4) * math.pi / 4
        return _ret
