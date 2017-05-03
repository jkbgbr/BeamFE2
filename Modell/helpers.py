# -*- coding: utf-8 -*-
import itertools
import numpy as np
import math


def compile_global_matrix(beams, stiffness=False, mass=False, geometrical=False):
    """
    Compiles the global stiffness matrix from the element matrices.
    :return: 
    """
    assert stiffness or mass or geometrical  # either or

    nodes = set(itertools.chain.from_iterable([x.nodes for x in beams]))  # nodes of the beams
    dof = beams[0].dof  # nr. of dofs
    _sumdof = len(nodes) * dof
    _empty = np.zeros(_sumdof ** 2)
    _empty = np.matrix(_empty.reshape(_sumdof, _sumdof))
    for b in beams:
        if stiffness:
            mtrx = b.matrix_in_global(mtrx=b.Ke)
        elif mass:
            mtrx = b.matrix_in_global(mtrx=b.Me)
        elif geometrical:
            mtrx = b.matrix_in_global(mtrx=b.Ke_geom)
        else:
            raise Exception('matrix not specified')
        _sti = (b.i.ID - 1) * dof  # starting element of the block for node i
        _eni = _sti + dof  # end element for the block if node i
        _stj = (b.j.ID - 1) * dof
        _enj = _stj + dof
        # upper left block
        _empty[_sti:_eni, _sti:_eni] += mtrx[:dof, :dof]
        # left lower block
        _empty[_stj:_enj, _sti:_eni] += mtrx[dof:, :dof]
        # upper right block
        _empty[_sti:_eni, _stj:_enj] += mtrx[:dof, dof:]
        # lower right block
        _empty[_stj:_enj, _stj:_enj] += mtrx[dof:, dof:]
    return _empty


def transfer_matrix(alpha, asdegree=False, blocks=2, dof=3):
    # matrix to rotate the stiffness matrix for compilation
    if asdegree:
        alpha = math.radians(alpha)
    cs = math.cos(alpha)
    ss = math.sin(alpha)
    if dof == 3:
        _block = np.matrix([[cs,    -ss,    0],
                            [ss,    cs,     0],
                            [0,     0,      1]])
    elif dof == 2:
        _block = np.matrix([[cs,    -ss],
                            [ss,    cs]])

    else:
        raise Exception('not implementeds, dof should be 2 or 3')

    _sumdof = blocks * dof
    _empty = np.zeros(_sumdof ** 2)
    _empty = np.matrix(_empty.reshape(_sumdof, _sumdof))

    for b in range(blocks):
        _sti = b * dof  # starting element of the block for node i
        _eni = _sti + dof  # end element for the block if node i
        _empty[_sti:_eni, _sti:_eni] += _block

    return _empty


def np_matrix_tolist(mtrx):
    """
    casts a numpy matrix into a flat list
    :param mtrx: 
    :return: 
    """
    return list(itertools.chain.from_iterable(mtrx.tolist()))
