#!/bin/env python3
#
# ba-proj
#
# This project is the implementation of the small
# CRT-RSA attack of Takayusa, Lu and Peng utilizing
# the coppersmith and L3 algorithm.
#
# (c) 2017 - Daniel Jankowski


import numpy


from coppersmith import generate_lattice
from fpylll import LLL
from matrix import Matrix
from test import test_lll_reduction


def main():
    print('#####################################')
    print('# Attack an small CRT-RSA exponents #')
    print('#####################################')
    print('')

    #NOTE: test l3 reduction in order to check if its working
    test_lll_reduction()

    generate_lattice(
            630548215070129547156718332495889632234434145411971275888376987603260225252787926135276738944105689100036295535868141424386536403649578707699128189491432138631900590774729214990015369102760964884776344849717811484309528915040117952098061886881,
            65535,
            8
            )
    pass


if __name__ == '__main__':
    main()

