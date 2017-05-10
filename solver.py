# -*- coding: utf-8 -*-
import scipy
# from scipy.linalg import eigh, inv, eig
import scipy.linalg as sp
from scipy.linalg import eig as namivan
import Modell.Results as Results
from Modell.helpers import *


def solve(structure, analysis=None):
    """
    solves the system, returns the vector of displacements.
    :return: 
    """
    assert analysis in ['nonlinear static', 'linear static', 'modal', 'buckling', 'all']
    assert structure.stiffness_matrix_is_ok
    assert structure.node_numbering_is_ok

    if analysis in ['linear static', 'all']:
        # linear static analysis
        disps = sp.inv(structure.K_with_BC) * structure.load_vector
        structure.results['linear static'] = Results.LinearStaticResult(structure=structure, displacements=[disps])

    # if analysis in ['nonlinear static', 'all']:
    #     # https: // www.ethz.ch / content / dam / ethz / special - interest / baug / ibk / structural - mechanics - dam / education / femI / Lecture_2b.pdf
    #     # do a linear
    #     # calculate the geometrical stiffness matrix
    #     # add them, do a linear
    #     # linear static analysis
    #     disps = inv(structure.K_with_BC) * structure.load_vector
    #     structure.results['linear static'] = Results.LinearStaticResult(structure=structure, displacements=[disps])

    # modal analyse
    if analysis in ['modal', 'all']:
        """
        units: 
        Length: [m]
        Density: [kg/m3]
        Young's module: [N/m2] 
        Time: [s]
        OR
        Length: [mm]
        Mass: [t] -> Density: [t/mm3]!
        Young's module: [N/mm2] 
        Time: [s]
        """
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

    if analysis in ['buckling']:
        solve(structure, analysis='linear static')
        beam = structure.beams[0]
        K = structure.condense(mtrx=structure.K_with_BC)
        KG = structure.condense(mtrx=structure.K_geom)
        eigvals = namivan(a=K, b=KG)[0]
        eigvecs = namivan(a=K, b=KG)[1]

        print([x for x in sorted(eigvals.real)])
        print([x for x in sorted(eigvals.real) if x>0])
        print('')

        shapes = [np.matrix([eigvecs[:, x]]).T for x in range(len(eigvecs))]  # casting to list of column matrices
        print(shapes)

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
        structure.draw(analysistype='buckling', mode=0)
        structure.draw(analysistype='buckling', mode=1)
        structure.draw(analysistype='buckling', mode=2)
        # print(sorted(scipy.linalg.eigvalsh(K, b=KG)))

        # eigvals, eigvecs = eigh(K, KG)

        # solve(structure, analysis='buckling')
        print(structure.results['linear static'].displacement_results)
        print((math.pi**2)*beam.EI/(900**2))
        print((math.pi**2)*beam.EI/((0.7*450)**2))

        # Linear buckling
        # not implemented yet!



    return True
