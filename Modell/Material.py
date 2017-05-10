# -*- coding: utf-8 -*-

"""
Material models for the FE elements.
Consult __init__.py for units
"""


class Material(object):
    def __init__(self, name=None):
        self.name = name


class CustomLinearElasticMaterial(Material):
    def __init__(self, name='custom', E=None, nu=None, rho=None):
        assert E is not None
        assert 0 < E
        assert nu is not None
        assert 0 < nu
        assert nu < 0.5
        assert rho > 0
        super(CustomLinearElasticMaterial, self).__init__(name=name)
        self.E = E  # the elastic modulus
        self.nu = nu  # Poisson's ration, dimensionless
        self.rho = rho  # Density

    def G(self):
        return self.E/(2 * 1 + self.nu)


class LinearElasticMaterial(Material):
    # abstract class
    def __init__(self, name=None):
        super(LinearElasticMaterial, self).__init__(name=name)
        self.E = None  # the elastic modulus in N/m2
        self.nu = None  # Poisson's ration, dimensionless
        self.rho = None

    def G(self):
        return self.E/(2 * 1 + self.nu)


class RigidLinearElasticMaterial(LinearElasticMaterial):
    def __init__(self, name='rigid'):
        super(RigidLinearElasticMaterial, self).__init__(name=name)
        self.E = 1e20  # the elastic modulus in N/m2
        self.nu = 0.5  # Poisson's ration, dimensionless

    def G(self):
        return 1e20


class Steel(LinearElasticMaterial):
    # normal ideal elastic steel
    def __init__(self, E=2.1e5, nu=0.3, name='steel'):
        super(Steel, self).__init__(name=name)
        self.rho = 7.850e-9  # to result correct modal results together with N, mm, second
        self.E = E  # N/m2
        self.nu = nu


class WeightlessSteel(Steel):
    # steel without weight
    def __init__(self, name='weightless steel'):
        super(WeightlessSteel, self).__init__(name=name)
        self.rho = 0  # to result correct modal results together with N, mm, second


class RigidSteel(Steel):
    # steel considered infinitely rigid
    def __init__(self, name='rigid steel'):
        super(RigidSteel, self).__init__(name=name)
        self.rho = 7.85e-9  # to result correct modal results together with N, mm, second
