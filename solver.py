# -*- coding: utf-8 -*-
from scipy.linalg import eigh, inv
import math


def solve(structure, analysis=None):
    """
    solves the system, returns the vector of displacements.
    :return: 
    """
    assert analysis in ['linear_elastic', 'modal', 'all']
    assert structure.stiffness_matrix_is_ok
    assert structure.node_numbering_is_ok

    if analysis in ['linear_elastic', 'all']:
        # linear static analysis
        structure.displacements = inv(structure.K_with_BC) * structure.q
        structure.displacements_for_beams()

    elif analysis in ['buckling']:
        raise NotImplementedError
        # Linear buckling
        # not implemented yet!

    # modal analyse
    elif analysis in ['modal', 'all']:
        eigvals, eigvecs = eigh(structure.K_with_BC, structure.M)
        structure.frequencies = [math.sqrt(x ) /( 2 *math.pi) for x in eigvals if x > 0]
        structure.nodal_shapes = eigvecs.T

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