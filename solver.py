# -*- coding: utf-8 -*-
import scipy
# from scipy.linalg import eigh, inv, eig
import scipy.linalg as sp
from scipy.linalg import eig as namivan
import Modell.Results as Results
from Modell.helpers import *

EPS = 1e-10

def round_results(mtrx, EPS=1e-10):
    """
    Rounds the displacement values smaller than threshold to zero.
    It is very much advised to use this ONLY on results like displacements, otherwise who know what happens.
    :param mtrx: numpy array like
    :return: 
    """
    mtrx = np.around(a=mtrx, decimals=10)  # rounding
    print(mtrx[abs(mtrx) < EPS])
    mtrx[abs(mtrx) < EPS] = 0  # values smaller than EPS will be zeroed
    return mtrx



def solve(structure, analysis=None):
    """
    solves the system, fills the structure.results attribute with values
    :return: False, if no appropriate loads are present. True, if the analysis gets done.
    """
    assert analysis in ['nonlinear static', 'linear static', 'modal', 'buckling', 'all']
    assert structure.stiffness_matrix_is_ok
    assert structure.node_numbering_is_ok

    # check if solving makes sense. if not, return False so we know
    if analysis in ['linear static', 'all']:
        if np.count_nonzero(structure._load_vector) == 0:
            return False
    elif analysis in ['modal', 'all']:
        if np.count_nonzero(structure.M_with_masses) == 0:
            return False
    elif analysis in ['buckling', 'all']:
        pass
    else:
        raise NotImplementedError

    if analysis in ['linear static', 'all']:
        # linear static analysis
        disps = sp.inv(structure.K_with_BC) * structure.load_vector

        structure.results['linear static'] = Results.LinearStaticResult(structure=structure, displacements=[disps])
        structure.results['linear static'].solved = True

    # modal analyse
    if analysis in ['modal', 'all']:
        # for details on unit choice see Modell/__init__.py
        K = structure.condense(mtrx=structure.K_with_BC)
        M = structure.condense(mtrx=structure.M_with_masses)
        eigvals, eigvecs = sp.eigh(K, M)
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

        shapes = [np.matrix([eigvecs[:, x]]).T for x in range(len(eigvecs))]  # casting to list of column matrices

        # shapes is to be updated (re-populated) to account for the rows and columns removed during condensing
        _ret = None
        positions_to_eliminate = structure.positions_to_eliminate
        for sh in shapes:
            for position in positions_to_eliminate:
                sh = np.insert(sh, position, [0], axis=0)
            if _ret is None:
                _ret = [sh]
            else:
                _ret.append(sh)

        structure.results['modal'] = Results.ModalResult(structure=structure, circular_freq=circfreq, modalshapes=_ret)
        structure.results['modal'].solved = True

    if analysis in ['buckling']:
        solve(structure, analysis='linear static')
        beam = structure.beams[0]
        K = structure.condense(mtrx=structure.K_with_BC)
        KG = structure.condense(mtrx=structure.K_geom)

        print('itt')

        # eigvals = namivan(a=K-KG)[0]
        # eigvecs = namivan(a=K-KG)[1]

        eigvals = namivan(a=K, b=KG)[0]
        eigvecs = namivan(a=K, b=KG)[1]

        # print(structure.K)
        # print('')
        # print(structure.K_with_BC)
        # print('')
        # print(K)
        #
        # print(KG)
        # print(K)

        # _eigs = sorted([x for x in eigvals.real if x>0])
        # beam = structure.beams[0]
        # ncr = math.pi ** 2 * beam.EI / (2 * beam.l) ** 2
        # print([x/ncr for x in _eigs])


        exit()

        # print((math.pi**2)*beam.EI/(900**2))
        # print([x for x in sorted(eigvals.real) if x>0])
        # print('')

        shapes = [np.matrix([eigvecs[:, x]]).T for x in range(len(eigvecs))]  # casting to list of column matrices

        # shapes is to be updated (re-populated) to account for the rows and columns removed during condensing
        _ret = None
        positions_to_eliminate = structure.positions_to_eliminate
        for sh in shapes:
            for position in positions_to_eliminate:
                sh = np.insert(sh, position, [0], axis=0)
            if _ret is None:
                _ret = [sh]
            else:
                _ret.append(sh)

        # print(_ret[0][0::3])
        # print(_ret[0][1::3])
        # print(_ret[0][2::3])

        structure.results['buckling'] = Results.BucklingResult(structure=structure, criticals=eigvals, bucklingshapes=_ret)

        # structure.draw(analysistype='linear static')
        # structure.draw(analysistype='buckling', mode=0)
        # structure.draw(analysistype='buckling', mode=1)
        # structure.draw(analysistype='buckling', mode=2)
        # print(sorted(scipy.linalg.eigvalsh(K, b=KG)))

        # eigvals, eigvecs = eigh(K, KG)

        # solve(structure, analysis='buckling')
        # print(structure.results['linear static'].displacement_results)

        # Linear buckling
        # not implemented yet!

    return True
