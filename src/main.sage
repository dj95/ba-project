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

    # test the lll algorithm
    test_lll_reduction()

    # generate the lattice for our parameters
    matrix = generate_lattice(
            630548215070129547156718332495889632234434145411971275888376987603260225252787926135276738944105689100036295535868141424386536403649578707699128189491432138631900590774729214990015369102760964884776344849717811484309528915040117952098061886881,
            65535,
            8
            )

    # reduce it
    #reduced_matrix = matrix.LLL()

    # use groebner
    pass


if __name__ == '__main__':
    main()
