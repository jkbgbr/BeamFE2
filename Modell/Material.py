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


class RigidMaterial(Material):
    def __init__(self):
        super(RigidMaterial, self).__init__()
        self.E = 1e20  # the elastic modulus in N/m2
        self.nu = 0.5  # Poisson's ration, dimensionless

    def G(self):
        return 1e20


class Steel(LinearElasticMaterial):
    # normal ideal elastic steel
    def __init__(self, E=2.1e5, nu=0.3):
        super(Steel, self).__init__()
        self.rho = 7.850e-9  # to result correct modal results together with N, mm, second
        self.E = E  # N/m2
        self.nu = nu


class WeightlessSteel(Steel):
    # steel without weight
    def __init__(self):
        super(WeightlessSteel, self).__init__()
        self.rho = 0  # to result correct modal results together with N, mm, second


class RigidSteel(Steel):
    # steel considered infinitely rigid
    def __init__(self):
        super(RigidSteel, self).__init__()
        self.rho = 7.85e-9  # to result correct modal results together with N, mm, second
