# -*- coding: utf-8 -*-

import math


class AnalysisResult(object):
    def __init__(self, parent=None):
        self.parent = parent


class LinearStaticResult(AnalysisResult):
    def __init__(self, parent=None):
        super(LinearStaticResult, self).__init__(parent=parent)
        self.parent = parent
        self.displacements = None


class ModalResult(AnalysisResult):
    def __init__(self, parent=None, circular_freq=None, modalshapes=None):
        super(ModalResult, self).__init__(parent=parent)
        self.parent = parent
        self.circular_frequencies = circular_freq
        self.displacements = modalshapes

    @property
    def frequencies(self):
        return [x/(2*math.pi) for x in self.circular_frequencies]

    @property
    def periods(self):
        return [1./x for x in self.frequencies]


class BucklingResult(AnalysisResult):
    def __init__(self, parent=None, criticals=None, bucklingshapes=None):
        super(BucklingResult, self).__init__(parent=parent)
        self.parent = parent
        self.criticals = criticals
        self.displacements = bucklingshapes
