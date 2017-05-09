# -*- coding: utf-8 -*-

"""
Material models for the FE elements.
Consult __init__.py for units
"""


class Material(object):
    def __init__(self):
        pass


class LinearElasticMaterial(Material):
    def __init__(self):
        super(LinearElasticMaterial, self).__init__()
        self.E = None  # the elastic modulus in N/m2
        self.nu = None  # Poisson's ration, dimensionless

    def G(self):
        return self.E/(2 * 1 + self.nu)


class Steel(LinearElasticMaterial):
    def __init__(self, E=2.1e5, nu=0.3):
        super(Steel, self).__init__()
        self.rho = 7.850e-9  # to result correct modal results together with N, mm, second
        self.E = E  # N/m2
        self.nu = nu
