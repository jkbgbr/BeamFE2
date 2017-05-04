# -*- coding: utf-8 -*-
from scipy.linalg import eigh, inv
import numpy as np
import pprint as pp
import math
import Modell.Results as Results
from Modell.helpers import *


def displacement_postprocessor(struct, disps):

    # storing the displacements for the structure. Values are in the global coordinate system.
    globaldisps = {}
    for dindex, dofname in enumerate(struct.dofnames):
        globaldisps[dofname] = disps[dindex::struct.dof]

    # storing the displacements for the beams, for each in its own local system
    beamdisps = {}
    for beam in struct.beams:
        beamdisps[beam] = {}
        T = transfer_matrix(-beam.direction, asdegree=False, blocks=1, blocksize=2)
        for dindex, dofname in enumerate(beam.dofnames):
            _akt = np_matrix_tolist(disps[dindex::struct.dof])
            beamdisps[beam][dofname] = T * np.matrix([_akt[beam.i.ID - 1], _akt[beam.j.ID - 1]]).T

    # displacement vectors for the beam elements. For each beam a 2xDOF numpy matrix
    beam_dispvectors = {}
    for beam in struct.beams:
        beam_dispvectors[beam] = []
        T = transfer_matrix(-beam.direction, asdegree=False, blocks=1, blocksize=3)
        _sta = struct.position_in_matrix(nodeID=beam.i.ID, DOF='ux')  # starting position in the global displacement vector for node i
        _end = struct.position_in_matrix(nodeID=beam.j.ID, DOF='ux')  # starting position in the global displacement vector for node j
        beam_dispvectors[beam] = np.concatenate((T * disps[_sta:_sta+beam.dof], T * disps[_end:_end+beam.dof]), axis=0)

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