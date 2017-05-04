# -*- coding: utf-8 -*-

import math
from Modell.helpers import *


class AnalysisResult(object):
    def __init__(self, parent=None):
        self.parent = parent
        self.displacements = None

    @property
    def resulting_displacements(self):
        disps = self.displacements
        print(disps)
        for disp in disps:
            print(disp)

        exit()
            # print([math.sqrt(x ** 2 + y ** 2) for x, y in zip(disp['ux'], disp['uy'])])
        return [[math.sqrt(x ** 2 + y ** 2) for x, y in zip(disp['ux'], disp['uy'])] for disp in disps]

    #
    # @property
    # def displacement_vector(self):
    #     """ die Verschiebungen als vector damit es mit  """
    #     disps = self.displacements
    #     print(disps)
    #     print(disps[0]['ux'])
    #     _ret = []
    #     for disp in disps:
    #         _ret.append([])
    #         for i in range(len(disp['ux'])):
    #             for dof in ['ux', 'uy', 'rotz']:
    #                 _ret[-1].append(np_matrix_tolist(disp[dof][i]))
    #     print(_ret)
    #     return _ret

# def resulting_displacement(self):
#     _ux = self.displacement_component(component='ux')
#     _uy = self.displacement_component(component='uy')
#     return [math.sqrt(x ** 2 + y ** 2) for x, y in zip(_ux, _uy)]


class LinearStaticResult(AnalysisResult):
    def __init__(self, parent=None, displacements=None, displacement_vector=None):
        super(LinearStaticResult, self).__init__(parent=parent)
        self.parent = parent
        self.displacements = displacements
        self.displacements_as_vector = displacement_vector

    @property
    def displacement_vector(self):
        return self.displacements_as_vector


class ModalResult(AnalysisResult):
    def __init__(self, parent=None, circular_freq=None, modalshapes=None, modalshape_vector=None):
        super(ModalResult, self).__init__(parent=parent)
        self.parent = parent
        self.circular_frequencies = circular_freq
        self.displacements = modalshapes
        self.modalshape_vector = modalshape_vector

    @property
    def displacement_vector(self):
        return self.modalshape_vector

    @property
    def frequencies(self):
        return [x/(2*math.pi) for x in self.circular_frequencies]

    @property
    def periods(self):
        return [1./x for x in self.frequencies]


class BucklingResult(AnalysisResult):
    def __init__(self, parent=None, criticals=None, bucklingshapes=None, bucklingshape_vector=None):
        super(BucklingResult, self).__init__(parent=parent)
        self.parent = parent
        self.criticals = criticals
        self.displacements = bucklingshapes
        self.bucklingshape_vector = bucklingshape_vector

    @property
    def displacement_vector(self):
        return self.bucklingshape_vector
