#!/bin/env python3
#
# ba-proj
#
# This project is the implementation of the small
# CRT-RSA attack of Takayusa, Lu and Peng utilizing
# the coppersmith and L3 algorithm.
#
# (c) 2017 - Daniel Jankowski


from fpylll import LLL
from matrix import Matrix


def test_lll_reduction():
    """
    Test the fpylll library and its L3-algorithm
    """
    print('==> L3-Test')
    # define a matrix
    z = [[1, -1, 3], [1, 0, 5], [1, 2, 6]]

    # generate the matrix object for fpylll from our matrix
    matrix = Matrix.from_matrix(z)

    # print the matrix and if its reduced
    print('  -> Start matrix')
    print(matrix)
    print('  -> is reduced: {}'.format(LLL.is_reduced(matrix)))
    print()

    # reduce the matrix with fpylll
    reduction = LLL.reduction(matrix)

    # print the reduced matrix and if its reduced
    print('  -> Reduced matrix')
    print(reduction)
    print('  -> is reduced: {}'.format(LLL.is_reduced(reduction)))
    print('')
