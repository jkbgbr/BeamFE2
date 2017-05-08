# -*- coding: utf-8 -*-
from scipy.linalg import eigh, inv
import numpy as np
import pprint as pp
import math
import Modell.Results as Results
from Modell.helpers import *


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
        disps = inv(structure.K_with_BC) * structure.load_vector
        structure.results['linear static'] = Results.LinearStaticResult(structure=structure, displacements=[disps])

    # modal analyse
    if analysis in ['modal', 'all']:
        eigvals, eigvecs = eigh(structure.K_with_BC, structure.M_with_masses)
        try:
            circfreq = [math.sqrt(x) for x in eigvals if x > 0]
            # circfreq = [math.sqrt(x) for x in eigvals]
        except ValueError:
            print(eigvals)
            print('negative eigenvalues found')
            for indi, i in enumerate([x for x in eigvals if x < 0]):
                print(indi+1, i)
            print('')
            import time
            time.sleep(0.5)
            raise

        shapes = [np.matrix([eigvecs[:, x]]).T for x in range(len(eigvecs))]  # casting to list of coulmn matrices

        structure.results['modal'] = Results.ModalResult(structure=structure, circular_freq=circfreq, modalshapes=shapes)

    if analysis in ['buckling']:
        raise NotImplementedError
        # Linear buckling
        # not implemented yet!



    return True
