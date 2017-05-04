# -*- coding: utf-8 -*-
from scipy.linalg import eigh, inv
import numpy as np
import pprint as pp
import math
import Modell.Results as Results
from Modell.helpers import *


def displacement_postprocessor(struct, disps):
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
    globaldisps[1]['ux'][4] is a matrix with both ux results of beam ID 4 in mode 1
    beamdisps[5][beam obj.]['ux'] is a matrix with both ux results of beam obj. in mode 5
    beam_dispvectors[3][beam obj.] is a matrix with all displacements of beam obj. in mode 3
    
    :param struct: 
    :param disps: 
    :return: 
    """


    # here the result is a dict with keys: DOFs, values: an array, containing the shape's deflections for the given DOF.
    # however, we needed a list of dicts where each dict contains all three components.
    globaldisps = []
    beamdisps = []
    beam_dispvectors = []
    for disp in disps:

        # storing the displacements for the structure. Values are in the global coordinate system.
        _akt = {}
        for dindex, dofname in enumerate(struct.dofnames):
            _akt[dofname] = disp[dindex::struct.dof]
        globaldisps.append(_akt)

        # storing the displacements for the beams, for each in its own local system
        _akt_1 = {}
        _akt_2 = {}
        for beam in struct.beams:
            _akt_1[beam] = {}
            T = transfer_matrix(-beam.direction, asdegree=False, blocks=1, blocksize=2)
            for dindex, dofname in enumerate(beam.dofnames):
                _aktvalues = disp[dindex::beam.dof]
                _akt_1[beam][dofname] = T * np.matrix([_aktvalues[beam.i.ID - 1], _aktvalues[beam.j.ID - 1]]).T

            # displacement vectors for the beam elements. For each beam a 2xDOF numpy matrix
            _akt_2[beam] = []
            _T = transfer_matrix(-beam.direction, asdegree=False, blocks=1, blocksize=3)
            _sta = struct.position_in_matrix(nodeID=beam.i.ID, DOF='ux')  # starting position in the global displacement vector for node i
            _end = struct.position_in_matrix(nodeID=beam.j.ID, DOF='ux')  # starting position in the global displacement vector for node j
            _akt_2[beam] = np.concatenate((_T * np.matrix(disp[_sta:_sta+beam.dof]).T, _T * np.matrix(disp[_end:_end+beam.dof]).T), axis=0)

        beamdisps.append(_akt_1)
        beam_dispvectors.append(_akt_2)

    return globaldisps, beamdisps, beam_dispvectors


def solve(structure, analysis=None):
    """
    solves the system, returns the vector of displacements.
    :return: 
    """
    assert analysis in ['linear static', 'modal', 'all']
    assert structure.stiffness_matrix_is_ok
    assert structure.node_numbering_is_ok

    if analysis in ['linear static', 'all']:
        # linear static analysis
        disps = inv(structure.K_with_BC) * structure.q
        globaldisps, beamdisps, beam_dispvectors = displacement_postprocessor(structure, disps)

        structure.results['linear static'] = Results.LinearStaticResult(parent=structure, displacements=[globaldisps], displacement_vector=disps)
        for beam in structure.beams:
            beam.results['linear static'] = Results.LinearStaticResult(parent=beam, displacements=[beamdisps[beam]], displacement_vector=beam_dispvectors[beam])

    if analysis in ['buckling']:
        raise NotImplementedError
        # Linear buckling
        # not implemented yet!

    # modal analyse
    if analysis in ['modal', 'all']:
        eigvals, eigvecs = eigh(structure.K_with_BC, structure.M)
        circfreq = [math.sqrt(x) for x in eigvals]
        globaldisps, beamdisps, beam_dispvectors = displacement_postprocessor(struct=structure, disps=eigvecs.T)

        # globaldisps[1]['ux'][4] is a matrix with both ux results of beam ID 4 in mode 1
        print(globaldisps[0]['uy'][0])
        print(globaldisps[0]['uy'][1])
        print(globaldisps[0]['uy'][2])
        print(globaldisps[0]['uy'][3])
        print(globaldisps[0]['uy'][4])
        print(globaldisps[0]['uy'][5])
        print(globaldisps[0]['uy'][6])
        exit()

        structure.results['modal'] = Results.ModalResult(parent=structure, circular_freq=circfreq, modalshapes=globaldisps)
        for beam in structure.beams:

            beam.results['modal'] = Results.ModalResult(parent=beam, circular_freq=circfreq, modalshapes=beamdisps[beam], modalshape_vector=beam_dispvectors[beam])


        # structure.frequencies = [math.sqrt(x ) /( 2 *math.pi) for x in eigvals if x > 0]
        # structure.nodal_shapes = eigvecs.T

        # print(eigvals[0:5])
        # print([math.sqrt(x) for x in eigvals if x > 0][0:5])
        # print([math.sqrt(x)/(2*3.1415) for x in eigvals if x > 0][0:5])
        # print('')
        # for shind, sh in enumerate(eigvecs[0:5]):
        #
        #     print('')
        #     _ux = sh[0::3]
        #     _uy = sh[1::3]
        #     _rotz = sh[2::3]
        #     print('ux:', _ux)
        #     print('uy:', _uy)
        #     print('rotz:', _rotz)
        #
        #     # Two subplots, the axes array is 1-d
        #     f, axarr = plt.subplots(2)
        #     axarr[0].plot(list(range(len(_uy))), _ux)
        #     axarr[0].set_title('#%d, f=%.2f Hz' % (shind+1, math.sqrt(eigvals[shind]) / 2*3.1415))
        #     axarr[1].plot(list(range(len(_uy))), _uy)
        #     plt.show()


    return True