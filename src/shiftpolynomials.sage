#!/usr/bin/env sage
#
# ba-proj
#
# This project is the implementation of the small
# CRT-RSA attack of Takayusa, Lu and Peng utilizing
# the coppersmith and L3 algorithm.
#
# (c) 2018 - Daniel Jankowski


def g(i1, i2, j1, j2, u, N, e, m):
    load('./polynomials.sage')
    load('./substitute.sage')

    # define the polynomial ring
    R.<xp1, xp2, xq1, xq2, yp, yq > = PolynomialRing(ZZ, order='deglex')

    fp1_val = fp1_sub_xp1(N, e, m)
    fp2_val = fp2_sub_xp2(N, e, m)
    h_val = h(N, e, m)

    polynomial = 1

    polynomial *= xp1**j1
    polynomial *= xp2**j2
    polynomial *= yq**(int(floor((i1 + i2) / 2)))
    polynomial *= fp1_val**i1
    polynomial *= fp2_val**i2
    polynomial *= h_val**u
    polynomial *= e**(m - (i1 + i2 + u))

    if polynomial == 0:
        return 0

    polynomial = substitute_y(polynomial, N, e, m)
    polynomial = substitute_x(polynomial, e, m)

    return polynomial


def gp(i1, i2, j1, N, e, m):
    load('./polynomials.sage')

    # define the polynomial ring
    R.<xp1, xp2, xq1, xq2, yp, yq > = PolynomialRing(ZZ, order='deglex')

    fp1_val = fp1_sub_xp1(N, e, m)
    fp2_val = fp2_sub_xp2(N, e, m)

    polynomial = 1

    polynomial *= yq**(int(floor((i1 + i2) / 2)) - j1)
    polynomial *= (fp1_val)**i1
    polynomial *= (fp2_val)**i2
    polynomial *= e**(m - (i1 + i2))

    polynomial = substitute_y(polynomial, N, e, m)
    polynomial = substitute_x(polynomial, e, m)

    return polynomial


def gq(i1, i2, j2, N, e, m):
    load('./polynomials.sage')

    # define the polynomial ring
    R.<xp1, xp2, xq1, xq2, yp, yq > = PolynomialRing(ZZ, order='deglex')

    fp1_val = fp1_sub_xp1(N, e, m)
    fp2_val = fp2_sub_xp2(N, e, m)

    polynomial = 1

    polynomial *= yq**(int(floor((i1 + i2) / 2)) + j2)
    polynomial *= fp1_val**i1
    polynomial *= fp2_val**i2
    polynomial *= e**(m - (i1 + i2))

    polynomial = substitute_y(polynomial, N, e, m)
    polynomial = substitute_x(polynomial, e, m)

    return polynomial
