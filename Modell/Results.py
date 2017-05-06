# -*- coding: utf-8 -*-

import math
from Modell.helpers import *


class AnalysisResult(object):
    """
    This function creates the somehow ordered results of the analyses.
    The resulting three element tuple contains:
    - globaldisps: displacements of the structure in the global coordinate system
    - beamdisps: displacements of the beams (per beam) in the beams local coordinate system
    - beam_dispvectors: the displacement vectors of the beams in the local system
    In all cases the result is a list, with as many elements as many results the analysis provides. For static this is 
    just one, for buckling, modal more.
    Each element of the list contains a dict with the keys the DOFs of the element. The number of elements in a value is
    the num,ber of nodes the element has.
    
    beamdisps[5][beam obj.]['ux'] is a matrix with both ux results of beam obj. in mode 5
    beam_dispvectors[3][beam obj.] is a matrix with all displacements of beam obj. in mode 3

    :param struct: 
    :param disps: 
    :return: 
    """

    def __init__(self, structure=None):
        self.structure = structure
        self.displacement_results = []  # will be a list, even if only with one element

    def global_displacements(self, mode=0, asvector=False):
        # the displacements as a vector or partitioned in a dict by the dofs
        disps = self.global_displacement_vector(mode)
        if asvector:
            # disps is a vector of all nodal displacements in
            return disps  # the full displacement vector as a numpy matrix
        else:
            # globaldisps['ux'] is a matrix with ux results of the whole strucutre in the mode requested
            globaldisps = {}
            for dindex, dofname in enumerate(self.structure.dofnames):
                globaldisps[dofname] = disps[dindex::self.structure.dof]
            return globaldisps  # displacements partitioned by DOF, each a matrix

    def element_displacements(self, local=True, mode=0, beam=None, asvector=False):
        """
        displacements an element for a mode, partitioned in a dict. For each elem an NxDOF numpy matrix
        :param local: boolean. For True, the results are in the element local system, for False in the global.
        :param mode: No. of mode
        :param beam: beam to get the displacements for
        :param asvector: for True, a load_vector vector is provided (NxDOF numpy matrix), for False the components are 
        partitioned in a dict with keys corresponfing the DOFs
        :return: the displacements in the desired format
        """

        assert beam in self.structure.beams
        struct = self.structure
        disp = self.global_displacement_vector(mode)

        if local:
            _T = transfer_matrix(-beam.direction, asdegree=False, blocks=1, blocksize=3)
        else:
            _T = transfer_matrix(alpha=0, asdegree=False, blocks=1, blocksize=3)

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


class ModalResult(AnalysisResult):
    def __init__(self, structure=None, circular_freq=None, modalshapes=None):
        super(ModalResult, self).__init__(structure=structure)
        self.circular_frequencies = circular_freq
        self.displacement_results = modalshapes  # the result as a matrix

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
