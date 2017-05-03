# -*- coding: utf-8 -*-


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
    def __init__(self, E=2.1e11, nu=0.3):
        super(Steel, self).__init__()
        self.rho = 7850000  # g/m3
        self.E = E  # N/m2
        self.nu = nu
