# -*- coding: utf-8 -*-

"""
Material models for the FE elements.
Consult __init__.py for units
"""

VERY_LARGE_NUMBER = 1e40

class Material(object):
    def __init__(self):
        self.name = self.__class__.__name__


class LinearElasticMaterial(Material):
    # abstract class
    def __init__(self, E=None, nu=None, rho=None):
        assert E is not None
        assert 0 < E
        assert nu is not None
        assert 0 < nu
        assert nu < 0.5
        assert rho >= 0
        super(LinearElasticMaterial, self).__init__()
        self._E = E  # the elastic modulus in N/m2
        self._nu = nu  # Poisson's ration, dimensionless
        self._rho = rho

    def __repr__(self):
        return "LinearElasticMaterial(E=%r, nu=%r, rho=%r)" % (self.E, self.nu, self.rho)

    @property
    def E(self):
        return self._E

    @property
    def nu(self):
        return self._nu

    @property
    def rho(self):
        return self._rho

    @property
    def G(self):
        return self.E/(2 * 1 + self.nu)

    def set_weightless(self):
        self._rho = 0
        if 'weightless' not in self.name:
            self.name = ' '.join(['weightless', self.name])

    def set_rigid(self):
        self._E = VERY_LARGE_NUMBER
        if 'rigid' not in self.name:
            self.name = ' '.join(['rigid', self.name])


class Steel(LinearElasticMaterial):
    # ideal elastic steel
    def __init__(self):
        super(Steel, self).__init__(E=2.1e5, nu=0.3, rho=7.85e-9)

    def __repr__(self):
        return "Steel(E=%r, nu=%r, rho=%r)" % (self.E, self.nu, self.rho)
