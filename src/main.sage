#!/bin/env sage
#
# ba-proj
#
# This project is the implementation of the small
# CRT-RSA attack of Takayusa, Lu and Peng utilizing
# the coppersmith and L3 algorithm.
#
# (c) 2017 - Daniel Jankowski


def main():
    print('#####################################')
    print('# Attack an small CRT-RSA exponents #')
    print('#####################################')
    print('')

    # import sage modules
    load('./coppersmith.sage')
    load('./test.sage')
    load('./keygen.sage')

    # test the lll algorithm
    test_lll_reduction()

    keys = generate_keys(1024)

    # generate the lattice for our parameters
    matrix = generate_lattice(
            keys['N'],
            keys['e'],
            8
            )

    print('\n==> Got an {}x{} matrix'.format(matrix.nrows(), matrix.ncols()))
    print('==> matrix is reduced: {}'.format(matrix.is_LLL_reduced()))
    print('==> Reducing matrix with LLL-algorithm')

    # reduce it
    reduced_matrix = matrix.LLL()

    # use groebner
    pass


if __name__ == '__main__':
    main()