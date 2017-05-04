# -*- coding: utf-8 -*-

import math
from Modell.helpers import *


class AnalysisResult(object):
    def __init__(self, structure=None):
        self.structure = structure
        self.displacement_results = []  # will be a list, even if only with one element

    def global_displacements(self, mode=0, asvector=False):
        # the displacements as a vector or partitioned in a dict by the dofs
        disps = self.global_displacement_vector(mode)
        if asvector:
            return disps  # the full displacement vector as a numpy matrix
        else:
            globaldisps = {}
            for dindex, dofname in enumerate(self.structure.dofnames):
                globaldisps[dofname] = disps[dindex::self.structure.dof]
            return globaldisps  # displacements partitioned by DOF, each a matrix

    def element_displacements(self, mode=None, beam=None, asvector=False):
        # displacements an element for a mode, partitioned in a dict. For each elem an NxDOF numpy matrix
        assert beam in self.structure.beams
        struct = self.structure
        disp = self.global_displacement_vector(mode)
        _T = transfer_matrix(-beam.direction, asdegree=False, blocks=1, blocksize=3)
        _sta = struct.position_in_matrix(nodeID=beam.i.ID, DOF='ux')  # starting position in the global displacement vector for node i
        _end = struct.position_in_matrix(nodeID=beam.j.ID, DOF='ux')  # starting position in the global displacement vector for node j

        disp = np.concatenate([_T * np.matrix(disp[_sta:_sta+beam.dof]), _T * np.matrix(disp[_end:_end+beam.dof])], axis=0)

        if asvector:  # the displacement vector in the local coordinate system as a numpy matrix
            return disp
        else:
            _ret = {}
            for dindex, dofname in enumerate(beam.dofnames):
                _ret[dofname] = disp[dindex::beam.dof]
            return _ret  # displacements partitioned by DOF, each a matrix

    def global_displacement_vector(self, mode=None):
        # the displacement vector of the structure in the global system, as calculated by the solver
        # returned are the dsplacements for a mode (#index), in form of the displacement vector
        return self.displacement_results[mode]




class LinearStaticResult(AnalysisResult):
    def __init__(self, structure=None, displacements=None):
        super(LinearStaticResult, self).__init__(structure=structure)
        self.displacement_results = displacements  # the result as a one-element vector

    # def global_displacements(self, mode=0, asvector=False):
    #     # the displacements as a vector or partitioned in a dict by the dofs
    #     disps = self.global_displacement_vector(mode)
    #     if asvector:
    #         return disps  # the full displacement vector as a numpy matrix
    #     else:
    #         globaldisps = {}
    #         for dindex, dofname in enumerate(self.structure.dofnames):
    #             globaldisps[dofname] = disps[dindex::self.structure.dof]
    #         return globaldisps  # displacements partitioned by DOF, each a matrix
    #
    # def element_displacements(self, mode=None, beam=None, asvector=False):
    #     # displacements an element for a mode, partitioned in a dict. For each elem an NxDOF numpy matrix
    #     assert beam in self.structure.beams
    #     struct = self.structure
    #     disp = self.global_displacement_vector(mode)
    #     _T = transfer_matrix(-beam.direction, asdegree=False, blocks=1, blocksize=3)
    #     _sta = struct.position_in_matrix(nodeID=beam.i.ID, DOF='ux')  # starting position in the global displacement vector for node i
    #     _end = struct.position_in_matrix(nodeID=beam.j.ID, DOF='ux')  # starting position in the global displacement vector for node j
    #     disp = np.concatenate([_T * np.matrix(disp[_sta:_sta+beam.dof]), _T * np.matrix(disp[_end:_end+beam.dof])], axis=0)
    #
    #     if asvector:  # the displacement vector in the local coordinate system as a numpy matrix
    #         return disp
    #     else:
    #         _ret = {}
    #         for dindex, dofname in enumerate(beam.dofnames):
    #             _ret[dofname] = disp[dindex::beam.dof]
    #         return _ret  # displacements partitioned by DOF, each a matrix
    #
    # def global_displacement_vector(self, mode=None):
    #     # the displacement vector of the structure in the global system, as calculated by the solver
    #     # returned are the dsplacements for a mode (#index), in form of the displacement vector
    #     mode = 0  # we only have one, so index is there only to provide a common interface for all analyses
    #     return self.displacement_results[mode]


class ModalResult(AnalysisResult):
    def __init__(self, structure=None, circular_freq=None, modalshapes=None):
        super(ModalResult, self).__init__(structure=structure)
        self.circular_frequencies = circular_freq
        self.displacement_results = modalshapes  # the result as a matrix




    #
    # def displacements_detailed(self, index=0):
    #     # the displacements partitioned for the dofs, for the mode #index
    #     globaldisps = {}
    #     for dindex, dofname in enumerate(self.structure.dofnames):
    #         globaldisps[dofname] = self.displacement_results[index][dindex::self.structure.dof]
    #     yield globaldisps
    #
    # @property
    # def displacement_vector(self, index=None):
    #     # the displacement vector of the structure in the global system, as calculated by the solver, for the mode in #index
    #     yield self.displacement_results[index]

    @property
    def frequencies(self):
        return [x/(2*math.pi) for x in self.circular_frequencies]

    @property
    def periods(self):
        return [1./x for x in self.frequencies]


class BucklingResult(AnalysisResult):
    def __init__(self, structure=None, criticals=None, bucklingshapes=None):
        super(BucklingResult, self).__init__(structure=structure)
        self.criticals = criticals
        self.displacements = bucklingshapes

    @property
    def displacement_vector(self):
        raise NotImplementedError
