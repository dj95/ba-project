#!/bin/env python3
#
# ba-proj
#
# This project is the implementation of the small
# CRT-RSA attack of Takayusa, Lu and Peng utilizing
# the coppersmith and L3 algorithm.
#
# (c) 2017 - Daniel Jankowski


import coppersmith
import numpy


from fpylll import LLL
from matrix import Matrix


def main():
    print('#####################################')
    print('# Attack an small CRT-RSA exponents #')
    print('#####################################')
    print('')

    print('==> L3-Test')
    z = [[1, -1, 3], [1, 0, 5], [1, 2, 6]]
    matrix = Matrix.from_matrix(z)
    print('  -> Start matrix')
    print(matrix)
    print('  -> is reduced: {}'.format(LLL.is_reduced(matrix)))
    print()

    reduction = LLL.reduction(matrix)

    print('  -> Reduced matrix')
    print(reduction)
    print('  -> is reduced: {}'.format(LLL.is_reduced(reduction)))
    pass


if __name__ == '__main__':
    main()

