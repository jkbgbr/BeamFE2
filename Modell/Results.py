# -*- coding: utf-8 -*-


class AnalysisResult(object):
    def __init__(self, parent=None):
        self.parent = parent
        self.linear_static = None
        self.modal = None
        self.buckling = None
